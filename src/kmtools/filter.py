from pathlib import Path

from kmtools.exceptions import FileValidationError, FilterError
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
        self, reference, km_output, output, output_type, count_threshold, verbose=False
    ):
        self.reference = reference
        self.km_output = km_output
        self.output = output
        self.output_type = output_type
        self.count_threshold = count_threshold
        self.verbose = verbose

        self.reference_df = None
        self.km_output_df = None

        # Extract sample name from km_output filename using Path for robustness
        self.sample_name = Path(self.km_output).stem
        self.output_df = []

    def run(self):
        Utils.log(
            f"Filtering {self.km_output} using reference {self.reference}", self.verbose
        )
        self._validate_input_files()
        self.verify_reference()
        self.verify_km_output()

        self.run_filtering()

        self.write_output()

        print(f"Filtered results written to {self.output}")

    def _validate_input_files(self):
        """Check that input files exist before attempting to read them."""
        if not Path(self.reference).exists():
            raise FileValidationError(
                f"Reference file not found: {self.reference}"
            )
        if not Path(self.km_output).exists():
            raise FileValidationError(
                f"KM output file not found: {self.km_output}"
            )

    def get_ref_alt_pos_from_variant(self, variant_name) -> KmVariant:
        if not isinstance(variant_name, str) or variant_name.strip() == "":
            return None, None, None, None

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

    def get_calculated_reference_alt(self, reference_row, km_row) -> str:
        target_seq_loc = self.split_query(km_row["Query"])
        ref_seq = str(km_row["Reference_sequence"]).upper()
        
        pos = int(reference_row["POS"])

        ref_allele = str(reference_row["REF"]).upper()
        alt_allele = str(reference_row["ALT"]).upper()

        if not (target_seq_loc.start <= pos <= target_seq_loc.end):
            return ""

        # compute 0-based offset of the variant within the target sequence
        offset = pos - target_seq_loc.start
        
        # Replace the reference allele at the computed offset with the ALT allele
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
        # This should not happen
        if not filter_result:
            return
        
        # Define common fields that are always present
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

        # Add KM-specific fields based on whether filter passed
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

    def write_output(self):
        output_df = pd.DataFrame(self.output_df)

        match self.output_type:
            case "tsv":
                output_df.to_csv(self.output, sep="\t", index=False)
            case "csv":
                output_df.to_csv(self.output, index=False)
            case "xlsx":
                output_df.to_excel(self.output, index=False)

    def run_filtering(self):
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

    def verify_reference(self):
        required = {"CHROM", "POS", "REF", "ALT", "TYPE"}

        if self.reference.endswith(".csv"):
            self.reference_df = pd.read_csv(self.reference)
        elif self.reference.endswith(".tsv"):
            self.reference_df = pd.read_csv(self.reference, sep="\t")
        else:
            raise FileValidationError(
                f"Unsupported reference file format: {self.reference}. Must be .csv or .tsv"
            )

        if not required.issubset(self.reference_df.columns):
            missing = required - set(self.reference_df.columns)
            raise FilterError(
                f"Reference file is missing required columns: {', '.join(sorted(missing))}"
            )

        Utils.log(
            f"Reference file {self.reference} contains the columns {list(self.reference_df.columns)}",
            self.verbose,
        )

    def verify_km_output(self):
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
            raise FilterError(
                f"KM output file is missing required columns: {', '.join(sorted(missing))}"
            )

        Utils.log(
            f"KM output file {self.km_output} contains columns: {list(self.km_output_df.columns)}",
            self.verbose,
        )
