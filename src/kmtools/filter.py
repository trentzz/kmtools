from __future__ import annotations

from pathlib import Path
from typing import Optional

from kmtools.filter_types import (
    FilterCondition,
    FilterResult,
    KmVariant,
    TargetSequenceLocation,
)
from kmtools.utils import Utils
import pandas as pd


class Filter:
    def __init__(
        self,
        reference: str,
        km_output: str,
        output: str,
        output_type: str,
        count_threshold: int,
        use_alt: bool = False,
        verbose: bool = False,
    ) -> None:
        self.reference = reference
        self.km_output = km_output
        self.output = output
        self.output_type = output_type
        self.count_threshold = count_threshold
        self.use_alt = use_alt
        self.verbose = verbose

        self.reference_df: Optional[pd.DataFrame] = None
        self.km_output_df: Optional[pd.DataFrame] = None

        self.sample_name = Path(self.km_output).stem
        self.output_df: list[dict] = []

    def run(self) -> None:
        mode = "use-alt" if self.use_alt else "reference"
        Utils.log(
            f"Filtering {self.km_output} using {self.reference} (mode: {mode})",
            self.verbose,
        )

        if self.use_alt:
            self.verify_alt_reference()
        else:
            self.verify_reference()

        self.verify_km_output()

        if self.use_alt:
            self.run_alt_filtering()
        else:
            self.run_filtering()

        self.write_output()

        print(f"Filtered results written to {self.output}")

    # === Shared helpers ===

    def _read_reference_file(self) -> pd.DataFrame:
        """Read a reference file (csv or tsv) and return as DataFrame."""
        if self.reference.endswith(".csv"):
            return pd.read_csv(self.reference)
        elif self.reference.endswith(".tsv"):
            return pd.read_csv(self.reference, sep="\t")
        else:
            raise ValueError(
                f"Unsupported file format for {self.reference}: Must be .csv or .tsv"
            )

    def get_ref_alt_pos_from_variant(self, variant_name) -> Optional[KmVariant]:
        if not isinstance(variant_name, str) or variant_name.strip() == "":
            return None

        ref, alt = variant_name.split("/")
        ref_position, ref_allele = ref.split(":")
        alt_allele, alt_position = alt.split(":")

        return KmVariant(
            int(ref_position), ref_allele.upper(), int(alt_position), alt_allele.upper()
        )

    def split_query(self, query: str) -> TargetSequenceLocation:
        chrom, start, end = query.split("_")
        return TargetSequenceLocation(chromosome=chrom, start=int(start), end=int(end))

    def get_km_alt(self, km_row) -> str:
        return km_row["Sequence"].upper()

    # === Reference mode (original) ===

    def get_calculated_reference_alt(self, reference_row, km_row) -> str:
        target_seq_loc = self.split_query(km_row["Query"])
        ref_seq = str(km_row["Reference_sequence"]).upper()

        pos = int(reference_row["POS"])

        ref_allele = str(reference_row["REF"]).upper()
        alt_allele = str(reference_row["ALT"]).upper()

        if not (target_seq_loc.start <= pos <= target_seq_loc.end):
            return ""

        offset = pos - target_seq_loc.start

        calculated = ref_seq[:offset] + alt_allele + ref_seq[offset + len(ref_allele):]

        return calculated

    def filter_line(
        self, km_row, reference_row, km_alt, calculated_reference_alt
    ) -> FilterResult:
        filtering_conditions = [
            FilterCondition(
                "TYPE",
                km_row["Type"] == reference_row["TYPE"],
                f"Type mismatch: KM Type {km_row['Type']} != Reference Type {reference_row['TYPE']}",
            ),
            FilterCondition(
                "COUNT",
                int(km_row["Min_coverage"]) >= self.count_threshold,
                f"Insufficient KM count: {km_row['Min_coverage']} < {self.count_threshold}",
            ),
            FilterCondition(
                "ALT",
                km_alt == calculated_reference_alt,
                f"ALT mismatch: KM ALT {km_alt} != Calculated Reference ALT {calculated_reference_alt}",
            ),
            FilterCondition(
                "TYPE_NOT_REFERENCE",
                reference_row["TYPE"].upper() != "REFERENCE",
                "Reference TYPE is 'REFERENCE', skipping",
            ),
            FilterCondition(
                "INFO_IS_VS_REF",
                km_row["Info"] == "vs_ref",
                f"KM Info {km_row['Info']} is not 'vs_ref'",
            ),
        ]

        failed_count = [cond for cond in filtering_conditions if not cond.condition]
        failed_count_str = ""

        if len(failed_count) == 1 and failed_count[0].name == "COUNT":
            failed_count_str = failed_count[0].message

        return FilterResult(
            passed=all(cond.condition for cond in filtering_conditions),
            failed_count=failed_count_str,
        )

    def write_filtered_line(self, reference_row, km_row, filter_result) -> None:
        if not filter_result:
            return

        output_row = {
            "SAMPLE": self.sample_name,
            "CHROM": reference_row["CHROM"],
            "POS": reference_row["POS"],
            "REF": reference_row["REF"],
            "ALT": reference_row["ALT"],
            "TYPE": reference_row["TYPE"],
            "FOUND": "TRUE" if filter_result.passed else "FALSE",
            "FILTER_NOTES": filter_result.failed_count,
        }

        km_fields = {
            "KMER_VAF": "rVAF",
            "KMER_MIN_COVERAGE": "Min_coverage",
            "KMER_EXPRESSION": "Expression",
            "REF_SEQUENCE": "Reference_sequence",
            "VARIANT_SEQUENCE": "Sequence",
        }

        for out_field, km_field in km_fields.items():
            output_row[out_field] = km_row[km_field] if filter_result.passed else ""

        self.output_df.append(output_row)

    def run_filtering(self) -> None:
        for _, reference_row in self.reference_df.iterrows():
            matching_found = FilterResult(passed=False, failed_count="")
            matching_km_row = None
            for _, km_row in self.km_output_df.iterrows():
                km_alt = self.get_km_alt(km_row)
                calculated_reference_alt = self.get_calculated_reference_alt(
                    reference_row, km_row
                )
                filter_result = self.filter_line(
                    km_row, reference_row, km_alt, calculated_reference_alt
                )

                if filter_result.failed_count or filter_result.passed:
                    matching_found = filter_result
                    matching_km_row = km_row

                if filter_result.passed:
                    break

            self.write_filtered_line(reference_row, matching_km_row, matching_found)

    def verify_reference(self) -> None:
        required = {"CHROM", "POS", "REF", "ALT", "TYPE"}

        self.reference_df = self._read_reference_file()

        if not required.issubset(self.reference_df.columns):
            missing = required - set(self.reference_df.columns)
            raise ValueError(
                f"Reference file is missing required columns: {', '.join(sorted(missing))}"
            )

        Utils.log(
            f"Reference file {self.reference} contains the columns {list(self.reference_df.columns)}",
            self.verbose,
        )

    # === Use-alt mode ===

    def verify_alt_reference(self) -> None:
        """Verify the alt reference file has required columns: CHROM, ALT_SEQUENCE, TYPE."""
        required = {"CHROM", "ALT_SEQUENCE", "TYPE"}

        self.reference_df = self._read_reference_file()

        if not required.issubset(self.reference_df.columns):
            missing = required - set(self.reference_df.columns)
            raise ValueError(
                f"Alt reference file is missing required columns: {', '.join(sorted(missing))}. "
                f"--use-alt mode requires: CHROM, ALT_SEQUENCE, TYPE"
            )

        Utils.log(
            f"Alt reference file {self.reference} contains the columns "
            f"{list(self.reference_df.columns)}",
            self.verbose,
        )

    def filter_alt_line(
        self, km_row, alt_row, km_sequence: str, expected_alt_sequence: str
    ) -> FilterResult:
        """Filter a km row against a use-alt reference row."""
        filtering_conditions = [
            FilterCondition(
                "TYPE",
                km_row["Type"] == alt_row["TYPE"],
                f"Type mismatch: KM Type {km_row['Type']} != Alt Type {alt_row['TYPE']}",
            ),
            FilterCondition(
                "COUNT",
                int(km_row["Min_coverage"]) >= self.count_threshold,
                f"Insufficient KM count: {km_row['Min_coverage']} < {self.count_threshold}",
            ),
            FilterCondition(
                "ALT_SEQUENCE",
                km_sequence == expected_alt_sequence,
                f"Sequence mismatch: KM Sequence != expected ALT_SEQUENCE",
            ),
            FilterCondition(
                "INFO_IS_VS_REF",
                km_row["Info"] == "vs_ref",
                f"KM Info {km_row['Info']} is not 'vs_ref'",
            ),
        ]

        failed_count = [cond for cond in filtering_conditions if not cond.condition]
        failed_count_str = ""

        if len(failed_count) == 1 and failed_count[0].name == "COUNT":
            failed_count_str = failed_count[0].message

        return FilterResult(
            passed=all(cond.condition for cond in filtering_conditions),
            failed_count=failed_count_str,
        )

    def write_alt_filtered_line(self, alt_row, km_row, filter_result) -> None:
        """Write a filtered output line for use-alt mode."""
        if not filter_result:
            return

        output_row = {
            "SAMPLE": self.sample_name,
            "CHROM": alt_row["CHROM"],
            "ALT_SEQUENCE": alt_row["ALT_SEQUENCE"],
            "TYPE": alt_row["TYPE"],
            "FOUND": "TRUE" if filter_result.passed else "FALSE",
            "FILTER_NOTES": filter_result.failed_count,
        }

        km_fields = {
            "KMER_VAF": "rVAF",
            "KMER_MIN_COVERAGE": "Min_coverage",
            "KMER_EXPRESSION": "Expression",
            "REF_SEQUENCE": "Reference_sequence",
            "VARIANT_SEQUENCE": "Sequence",
        }

        for out_field, km_field in km_fields.items():
            output_row[out_field] = km_row[km_field] if filter_result.passed else ""

        self.output_df.append(output_row)

    def run_alt_filtering(self) -> None:
        """Run filtering in use-alt mode: match km Sequence against provided ALT_SEQUENCE."""
        for _, alt_row in self.reference_df.iterrows():
            expected_alt = str(alt_row["ALT_SEQUENCE"]).upper()
            chrom = str(alt_row["CHROM"])

            matching_found = FilterResult(passed=False, failed_count="")
            matching_km_row = None

            for _, km_row in self.km_output_df.iterrows():
                # Match chromosome from km Query
                km_query_loc = self.split_query(km_row["Query"])
                if km_query_loc.chromosome != chrom:
                    continue

                km_sequence = self.get_km_alt(km_row)
                filter_result = self.filter_alt_line(
                    km_row, alt_row, km_sequence, expected_alt
                )

                if filter_result.failed_count or filter_result.passed:
                    matching_found = filter_result
                    matching_km_row = km_row

                if filter_result.passed:
                    break

            self.write_alt_filtered_line(alt_row, matching_km_row, matching_found)

    # === Shared ===

    def write_output(self) -> None:
        output_df = pd.DataFrame(self.output_df)

        match self.output_type:
            case "tsv":
                output_df.to_csv(self.output, sep="\t", index=False)
            case "csv":
                output_df.to_csv(self.output, index=False)
            case "xlsx":
                output_df.to_excel(self.output, index=False)

    def verify_km_output(self) -> None:
        required = {
            "Database",
            "Query",
            "Type",
            "Variant_name",
            "rVAF",
            "Expression",
            "Min_coverage",
            "Start_offset",
            "Sequence",
            "Reference_expression",
            "Reference_sequence",
            "Info",
        }

        self.km_output_df = pd.read_csv(
            self.km_output, sep="\t", dtype=str, comment="#"
        )

        if not required.issubset(self.km_output_df.columns):
            missing = required - set(self.km_output_df.columns)
            raise ValueError(
                f"KM output file is missing required columns: {', '.join(sorted(missing))}"
            )

        Utils.log(
            f"KM output file {self.km_output} contains columns: {list(self.km_output_df.columns)}",
            self.verbose,
        )
