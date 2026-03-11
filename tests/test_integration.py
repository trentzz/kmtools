"""
Integration tests for kmtools.

These tests exercise the full pipeline end-to-end, ensuring that
the filter, merge, plot, and use-alt features work together correctly
with realistic data flowing through multiple components.
"""

from __future__ import annotations

import subprocess
import sys

import pandas as pd
import pytest
from pathlib import Path

from kmtools.filter import Filter
from kmtools.merge import Merge
from kmtools.plot import Plot
from kmtools.exceptions import MergeError


# ---------------------------------------------------------------------------
# Shared test data helpers
# ---------------------------------------------------------------------------

KM_OUTPUT_HEADER = (
    "Database\tQuery\tType\tVariant_name\trVAF\tExpression\t"
    "Min_coverage\tStart_offset\tSequence\tReference_expression\t"
    "Reference_sequence\tInfo"
)


def make_km_row(
    query: str = "chr1_114716091_114716161",
    variant_type: str = "Substitution",
    variant_name: str = "36:c/T:37",
    rvaf: str = "0.500",
    expression: str = "100.0",
    min_coverage: str = "10",
    sequence: str = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
    ref_expression: str = "5033.0",
    ref_sequence: str = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
    info: str = "vs_ref",
) -> str:
    return (
        f"file.jf\t{query}\t{variant_type}\t{variant_name}\t{rvaf}\t{expression}\t"
        f"{min_coverage}\t0\t{sequence}\t{ref_expression}\t{ref_sequence}\t{info}"
    )


def write_km_output(path: Path, rows: list[str]) -> None:
    path.write_text(KM_OUTPUT_HEADER + "\n" + "\n".join(rows) + "\n")


def write_reference_tsv(path: Path, rows: list[dict]) -> None:
    header = "CHROM\tPOS\tREF\tALT\tTYPE"
    lines = [header]
    for row in rows:
        lines.append(f"{row['CHROM']}\t{row['POS']}\t{row['REF']}\t{row['ALT']}\t{row['TYPE']}")
    path.write_text("\n".join(lines) + "\n")


def write_alt_reference_tsv(path: Path, rows: list[dict]) -> None:
    header = "CHROM\tALT_SEQUENCE\tTYPE"
    lines = [header]
    for row in rows:
        lines.append(f"{row['CHROM']}\t{row['ALT_SEQUENCE']}\t{row['TYPE']}")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Integration: Filter (reference mode) end-to-end
# ---------------------------------------------------------------------------


class TestFilterIntegration:
    """Full end-to-end tests for the filter workflow in reference mode."""

    def test_single_variant_found(self, tmp_path):
        """A single variant that passes all filter conditions."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(
                query="chr1_114716091_114716161",
                variant_type="Substitution",
                min_coverage="10",
                sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
                ref_sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
            ),
        ])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "TRUE"
        assert result.iloc[0]["CHROM"] == "chr1"
        assert result.iloc[0]["TYPE"] == "Substitution"
        assert float(result.iloc[0]["KMER_VAF"]) == pytest.approx(0.5)
        assert int(result.iloc[0]["KMER_MIN_COVERAGE"]) == 10

    def test_multiple_variants_mixed_results(self, tmp_path):
        """Multiple variants where some pass and some fail."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        # Variant 1: should pass (matching substitution)
        # Variant 2: should fail (type mismatch - reference says Deletion but km says Substitution)
        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
            {"CHROM": "chr2", "POS": "50000010", "REF": "AT", "ALT": "A", "TYPE": "Deletion"},
        ])

        write_km_output(km_file, [
            make_km_row(
                query="chr1_114716091_114716161",
                variant_type="Substitution",
                min_coverage="10",
                sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
                ref_sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
            ),
        ])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 2
        assert result.iloc[0]["FOUND"] == "TRUE"
        assert result.iloc[1]["FOUND"] == "FALSE"

    def test_low_coverage_fails_with_note(self, tmp_path):
        """Variant with only count issue gets a FILTER_NOTES entry."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        # Coverage of 1 is below threshold of 5
        write_km_output(km_file, [
            make_km_row(
                query="chr1_114716091_114716161",
                variant_type="Substitution",
                min_coverage="1",
                sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
                ref_sequence="TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
            ),
        ])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=5,
            verbose=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "FALSE"
        assert "Insufficient KM count" in str(result.iloc[0]["FILTER_NOTES"])

    def test_output_csv_format(self, tmp_path):
        """Verify CSV output format works correctly."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.csv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(min_coverage="10"),
        ])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="csv",
            count_threshold=2,
        ).run()

        result = pd.read_csv(out_file)
        assert len(result) == 1
        assert "FOUND" in result.columns
        assert "SAMPLE" in result.columns

    def test_sample_name_from_km_filename(self, tmp_path):
        """SAMPLE column should be derived from the km output filename stem."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "patient_001_km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [make_km_row(min_coverage="10")])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert result.iloc[0]["SAMPLE"] == "patient_001_km_output"

    def test_reference_type_skipped(self, tmp_path):
        """Variants with TYPE='Reference' should be marked as not found."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Reference"},
        ])

        write_km_output(km_file, [
            make_km_row(variant_type="Reference", min_coverage="100"),
        ])

        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "FALSE"


# ---------------------------------------------------------------------------
# Integration: Filter (use-alt mode) end-to-end
# ---------------------------------------------------------------------------


class TestFilterUseAltIntegration:
    """Full end-to-end tests for the --use-alt filter workflow."""

    def test_alt_mode_single_match(self, tmp_path):
        """Single alt sequence that matches a km row."""
        alt_ref = tmp_path / "alt_reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered_alt.tsv"

        seq = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"

        write_alt_reference_tsv(alt_ref, [
            {"CHROM": "chr1", "ALT_SEQUENCE": seq, "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(
                query="chr1_114716091_114716161",
                variant_type="Substitution",
                min_coverage="10",
                sequence=seq,
            ),
        ])

        Filter(
            reference=str(alt_ref),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
            verbose=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "TRUE"
        assert result.iloc[0]["ALT_SEQUENCE"] == seq
        assert float(result.iloc[0]["KMER_VAF"]) == pytest.approx(0.5)

    def test_alt_mode_wrong_chromosome_no_match(self, tmp_path):
        """Alt sequence on wrong chromosome should not match."""
        alt_ref = tmp_path / "alt_reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered_alt.tsv"

        seq = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"

        write_alt_reference_tsv(alt_ref, [
            {"CHROM": "chr99", "ALT_SEQUENCE": seq, "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(query="chr1_114716091_114716161", sequence=seq),
        ])

        Filter(
            reference=str(alt_ref),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "FALSE"

    def test_alt_mode_multiple_variants(self, tmp_path):
        """Multiple alt variants with some matching and some not."""
        alt_ref = tmp_path / "alt_reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered_alt.tsv"

        seq_match = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"
        seq_nomatch = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

        write_alt_reference_tsv(alt_ref, [
            {"CHROM": "chr1", "ALT_SEQUENCE": seq_match, "TYPE": "Substitution"},
            {"CHROM": "chr1", "ALT_SEQUENCE": seq_nomatch, "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(query="chr1_114716091_114716161", sequence=seq_match),
        ])

        Filter(
            reference=str(alt_ref),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 2
        assert result.iloc[0]["FOUND"] == "TRUE"
        assert result.iloc[1]["FOUND"] == "FALSE"

    def test_alt_mode_output_columns(self, tmp_path):
        """Verify use-alt mode produces the correct output columns."""
        alt_ref = tmp_path / "alt_reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered_alt.tsv"

        seq = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"

        write_alt_reference_tsv(alt_ref, [
            {"CHROM": "chr1", "ALT_SEQUENCE": seq, "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(query="chr1_114716091_114716161", sequence=seq),
        ])

        Filter(
            reference=str(alt_ref),
            km_output=str(km_file),
            output=str(out_file),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        ).run()

        result = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        expected_cols = {
            "SAMPLE", "CHROM", "ALT_SEQUENCE", "TYPE", "FOUND",
            "FILTER_NOTES", "KMER_VAF", "KMER_MIN_COVERAGE",
            "KMER_EXPRESSION", "REF_SEQUENCE", "VARIANT_SEQUENCE",
        }
        assert set(result.columns) == expected_cols


# ---------------------------------------------------------------------------
# Integration: Merge end-to-end
# ---------------------------------------------------------------------------


class TestMergeIntegration:
    """Full end-to-end tests for the merge workflow."""

    def _write_km_tsv(self, path: Path, rows: list[str]) -> None:
        write_km_output(path, rows)

    def test_merge_multiple_files(self, tmp_path, monkeypatch):
        """Merge several km output files into one."""
        monkeypatch.chdir(tmp_path)

        for i in range(4):
            f = tmp_path / f"chunk_{i}.txt"
            self._write_km_tsv(f, [
                make_km_row(min_coverage=str(i + 1), rvaf=f"0.{i}00"),
            ])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            verbose=True,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 4
        # All files should have been kept
        for i in range(4):
            assert (tmp_path / f"chunk_{i}.txt").exists()

    def test_merge_with_sort(self, tmp_path, monkeypatch):
        """Merge with --sort-by column."""
        monkeypatch.chdir(tmp_path)

        f1 = tmp_path / "chunk_0.txt"
        f2 = tmp_path / "chunk_1.txt"
        self._write_km_tsv(f1, [make_km_row(min_coverage="50")])
        self._write_km_tsv(f2, [make_km_row(min_coverage="10")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            sort_by="Min_coverage",
            verbose=True,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 2
        # Should be sorted by Min_coverage (ascending)
        coverages = list(result["Min_coverage"])
        assert coverages == sorted(coverages)

    def test_merge_with_drop_duplicates(self, tmp_path, monkeypatch):
        """Merge with --drop-duplicates removes exact duplicate rows."""
        monkeypatch.chdir(tmp_path)

        # Write identical rows in two files
        row = make_km_row(min_coverage="10")
        f1 = tmp_path / "chunk_0.txt"
        f2 = tmp_path / "chunk_1.txt"
        self._write_km_tsv(f1, [row])
        self._write_km_tsv(f2, [row])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            drop_duplicates=True,
            verbose=True,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 1

    def test_merge_deletes_inputs_by_default(self, tmp_path, monkeypatch):
        """Without --keep, merged input files are deleted."""
        monkeypatch.chdir(tmp_path)

        f1 = tmp_path / "chunk_0.txt"
        f2 = tmp_path / "chunk_1.txt"
        self._write_km_tsv(f1, [make_km_row()])
        self._write_km_tsv(f2, [make_km_row(min_coverage="5")])

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=False,
            verbose=True,
        ).run()

        assert output.exists()
        assert not f1.exists()
        assert not f2.exists()

    def test_merge_no_files_raises(self, tmp_path):
        """Merging with no matching files raises MergeError."""
        output = tmp_path / "merged.txt"
        with pytest.raises(MergeError, match="No input files found"):
            Merge(
                inputs=[str(tmp_path / "nonexistent_*.txt")],
                output=str(output),
                keep=True,
            ).run()

    def test_merge_large_dataset(self, tmp_path, monkeypatch):
        """Merge a larger dataset (100 files x 10 rows each = 1000 rows)."""
        monkeypatch.chdir(tmp_path)

        for i in range(100):
            f = tmp_path / f"chunk_{i:03d}.txt"
            rows = [
                make_km_row(min_coverage=str(j + 1), rvaf=f"0.{j:03d}")
                for j in range(10)
            ]
            self._write_km_tsv(f, rows)

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["chunk_*.txt"],
            output=str(output),
            keep=True,
            verbose=True,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 1000


# ---------------------------------------------------------------------------
# Integration: Plot end-to-end
# ---------------------------------------------------------------------------


class TestPlotIntegration:
    """Full end-to-end tests for the plot workflow."""

    def _write_filtered_tsv(self, path: Path, rows: list[dict]) -> None:
        if not rows:
            path.write_text(
                "SAMPLE\tCHROM\tPOS\tREF\tALT\tTYPE\tFOUND\tFILTER_NOTES\t"
                "KMER_VAF\tKMER_MIN_COVERAGE\tKMER_EXPRESSION\tREF_SEQUENCE\tVARIANT_SEQUENCE\n"
            )
            return

        df = pd.DataFrame(rows)
        df.to_csv(path, sep="\t", index=False)

    def test_plot_all_charts(self, tmp_path):
        """Generate all chart types from a filtered results file."""
        filtered = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        rows = [
            {
                "SAMPLE": "patient_001", "CHROM": "chr1", "POS": 100,
                "REF": "C", "ALT": "T", "TYPE": "Substitution",
                "FOUND": "TRUE", "FILTER_NOTES": "",
                "KMER_VAF": 0.45, "KMER_MIN_COVERAGE": 10,
                "KMER_EXPRESSION": 100.0, "REF_SEQUENCE": "ATCG",
                "VARIANT_SEQUENCE": "ATCG",
            },
            {
                "SAMPLE": "patient_001", "CHROM": "chr2", "POS": 200,
                "REF": "A", "ALT": "G", "TYPE": "Substitution",
                "FOUND": "FALSE", "FILTER_NOTES": "count too low",
                "KMER_VAF": "", "KMER_MIN_COVERAGE": "",
                "KMER_EXPRESSION": "", "REF_SEQUENCE": "",
                "VARIANT_SEQUENCE": "",
            },
            {
                "SAMPLE": "patient_002", "CHROM": "chr1", "POS": 300,
                "REF": "GC", "ALT": "G", "TYPE": "Deletion",
                "FOUND": "TRUE", "FILTER_NOTES": "",
                "KMER_VAF": 0.30, "KMER_MIN_COVERAGE": 5,
                "KMER_EXPRESSION": 50.0, "REF_SEQUENCE": "ATCG",
                "VARIANT_SEQUENCE": "ATCG",
            },
        ]

        self._write_filtered_tsv(filtered, rows)

        Plot(
            file=str(filtered),
            output_dir=str(plot_dir),
            charts="all",
            verbose=True,
        ).run()

        # Check all chart files were created
        assert (plot_dir / "vaf_distribution.png").exists()
        assert (plot_dir / "type_distribution.png").exists()
        assert (plot_dir / "sample_summary.png").exists()
        assert (plot_dir / "overall_summary.png").exists()

    def test_plot_single_chart(self, tmp_path):
        """Generate only one specific chart type."""
        filtered = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        rows = [
            {
                "SAMPLE": "s1", "CHROM": "chr1", "POS": 100,
                "REF": "C", "ALT": "T", "TYPE": "Substitution",
                "FOUND": "TRUE", "FILTER_NOTES": "",
                "KMER_VAF": 0.5, "KMER_MIN_COVERAGE": 10,
                "KMER_EXPRESSION": 100.0, "REF_SEQUENCE": "A",
                "VARIANT_SEQUENCE": "B",
            },
        ]

        self._write_filtered_tsv(filtered, rows)

        Plot(
            file=str(filtered),
            output_dir=str(plot_dir),
            charts="vaf",
            verbose=True,
        ).run()

        assert (plot_dir / "vaf_distribution.png").exists()
        assert not (plot_dir / "type_distribution.png").exists()

    def test_plot_empty_data(self, tmp_path):
        """Plot with empty data doesn't crash."""
        filtered = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        self._write_filtered_tsv(filtered, [])

        Plot(
            file=str(filtered),
            output_dir=str(plot_dir),
            charts="all",
            verbose=True,
        ).run()

        # No charts should be generated for empty data
        assert not (plot_dir / "vaf_distribution.png").exists()


# ---------------------------------------------------------------------------
# Integration: Filter -> Plot pipeline
# ---------------------------------------------------------------------------


class TestFilterToPlotPipeline:
    """Tests that chain filter output directly into plot input."""

    def test_filter_output_feeds_plot(self, tmp_path):
        """Run filter then plot on the output - full pipeline."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        filtered_file = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(min_coverage="10", rvaf="0.450"),
        ])

        # Step 1: Filter
        Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(filtered_file),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        assert filtered_file.exists()
        filter_result = pd.read_csv(filtered_file, sep="\t", dtype={"FOUND": str})
        assert len(filter_result) == 1
        assert filter_result.iloc[0]["FOUND"] == "TRUE"

        # Step 2: Plot from filter output
        Plot(
            file=str(filtered_file),
            output_dir=str(plot_dir),
            charts="all",
            verbose=True,
        ).run()

        assert (plot_dir / "vaf_distribution.png").exists()
        assert (plot_dir / "overall_summary.png").exists()


# ---------------------------------------------------------------------------
# Integration: Merge -> Filter pipeline
# ---------------------------------------------------------------------------


class TestMergeToFilterPipeline:
    """Tests that chain merge output into filter input."""

    def test_merge_then_filter(self, tmp_path, monkeypatch):
        """Merge multiple chunk files, then filter the merged result."""
        monkeypatch.chdir(tmp_path)

        ref_file = tmp_path / "reference.tsv"
        merged_file = tmp_path / "merged.txt"
        filtered_file = tmp_path / "filtered.tsv"

        # Create reference with 2 variants
        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        # Create 3 chunk files, each with the same variant detected
        for i in range(3):
            f = tmp_path / f"chunk_{i}.txt"
            write_km_output(f, [
                make_km_row(min_coverage=str(10 + i)),
            ])

        # Step 1: Merge
        Merge(
            inputs=["chunk_*.txt"],
            output=str(merged_file),
            keep=True,
            verbose=True,
        ).run()

        assert merged_file.exists()
        merged_df = pd.read_csv(merged_file, sep="\t")
        assert len(merged_df) == 3

        # Step 2: Filter
        Filter(
            reference=str(ref_file),
            km_output=str(merged_file),
            output=str(filtered_file),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(filtered_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 1
        assert result.iloc[0]["FOUND"] == "TRUE"


# ---------------------------------------------------------------------------
# Integration: Merge -> Filter -> Plot (full pipeline without chunk)
# ---------------------------------------------------------------------------


class TestFullPipeline:
    """Tests the full pipeline: merge -> filter -> plot."""

    def test_full_pipeline_reference_mode(self, tmp_path, monkeypatch):
        """Complete pipeline in reference mode: merge -> filter -> plot."""
        monkeypatch.chdir(tmp_path)

        ref_file = tmp_path / "reference.tsv"
        merged_file = tmp_path / "merged.txt"
        filtered_file = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        # Reference with 3 variants
        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
            {"CHROM": "chr2", "POS": "50000010", "REF": "AT", "ALT": "A", "TYPE": "Deletion"},
        ])

        # Create chunk files
        for i in range(4):
            f = tmp_path / f"chunk_{i}.txt"
            write_km_output(f, [
                make_km_row(min_coverage=str(10 + i), rvaf=f"0.{i + 1}00"),
            ])

        # Step 1: Merge
        Merge(
            inputs=["chunk_*.txt"],
            output=str(merged_file),
            keep=False,
            verbose=True,
        ).run()

        assert merged_file.exists()
        # Chunk files should be deleted
        for i in range(4):
            assert not (tmp_path / f"chunk_{i}.txt").exists()

        # Step 2: Filter
        Filter(
            reference=str(ref_file),
            km_output=str(merged_file),
            output=str(filtered_file),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        filter_result = pd.read_csv(filtered_file, sep="\t", dtype={"FOUND": str})
        assert len(filter_result) == 3
        found_count = (filter_result["FOUND"] == "TRUE").sum()
        not_found_count = (filter_result["FOUND"] == "FALSE").sum()
        assert found_count >= 1
        assert not_found_count >= 0

        # Step 3: Plot
        Plot(
            file=str(filtered_file),
            output_dir=str(plot_dir),
            charts="all",
            verbose=True,
        ).run()

        assert (plot_dir / "overall_summary.png").exists()
        assert (plot_dir / "type_distribution.png").exists()

    def test_full_pipeline_use_alt_mode(self, tmp_path, monkeypatch):
        """Complete pipeline in use-alt mode: merge -> filter(use-alt) -> plot."""
        monkeypatch.chdir(tmp_path)

        seq = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"
        alt_ref_file = tmp_path / "alt_reference.tsv"
        merged_file = tmp_path / "merged.txt"
        filtered_file = tmp_path / "filtered.tsv"
        plot_dir = tmp_path / "plots"

        write_alt_reference_tsv(alt_ref_file, [
            {"CHROM": "chr1", "ALT_SEQUENCE": seq, "TYPE": "Substitution"},
            {"CHROM": "chr1", "ALT_SEQUENCE": "AAAAAAA", "TYPE": "Substitution"},
        ])

        for i in range(2):
            f = tmp_path / f"chunk_{i}.txt"
            write_km_output(f, [
                make_km_row(
                    query="chr1_114716091_114716161",
                    sequence=seq,
                    min_coverage=str(5 + i),
                    rvaf=f"0.{3 + i}00",
                ),
            ])

        # Step 1: Merge
        Merge(
            inputs=["chunk_*.txt"],
            output=str(merged_file),
            keep=True,
            verbose=True,
        ).run()

        # Step 2: Filter (use-alt mode)
        Filter(
            reference=str(alt_ref_file),
            km_output=str(merged_file),
            output=str(filtered_file),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
            verbose=True,
        ).run()

        result = pd.read_csv(filtered_file, sep="\t", dtype={"FOUND": str})
        assert len(result) == 2
        assert result.iloc[0]["FOUND"] == "TRUE"
        assert result.iloc[1]["FOUND"] == "FALSE"
        assert "ALT_SEQUENCE" in result.columns

        # Step 3: Plot (overall summary should work even with alt-mode columns)
        Plot(
            file=str(filtered_file),
            output_dir=str(plot_dir),
            charts="overall,type",
            verbose=True,
        ).run()

        assert (plot_dir / "overall_summary.png").exists()
        assert (plot_dir / "type_distribution.png").exists()


# ---------------------------------------------------------------------------
# Integration: CLI invocation
# ---------------------------------------------------------------------------


class TestCLIIntegration:
    """Tests that the CLI entry point works correctly."""

    def test_version_flag(self):
        """--version flag should work."""
        result = subprocess.run(
            [sys.executable, "-m", "kmtools.kmtools", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "kmtools" in result.stdout

    def test_help_flag(self):
        """--help should list all subcommands."""
        result = subprocess.run(
            [sys.executable, "-m", "kmtools.kmtools", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "chunk" in result.stdout
        assert "merge" in result.stdout
        assert "filter" in result.stdout
        assert "plot" in result.stdout

    def test_filter_via_cli(self, tmp_path):
        """Run filter through the CLI entry point."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "114716126", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [make_km_row(min_coverage="10")])

        result = subprocess.run(
            [
                sys.executable, "-m", "kmtools.kmtools",
                "--verbose",
                "filter",
                "--reference", str(ref_file),
                "--km-output", str(km_file),
                "--output", str(out_file),
                "--output-type", "tsv",
                "--count-threshold", "2",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert out_file.exists()

        df = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(df) == 1
        assert df.iloc[0]["FOUND"] == "TRUE"

    def test_merge_via_cli(self, tmp_path, monkeypatch):
        """Run merge through the CLI entry point."""
        monkeypatch.chdir(tmp_path)

        for i in range(3):
            f = tmp_path / f"chunk_{i}.txt"
            write_km_output(f, [make_km_row(min_coverage=str(i + 1))])

        out_file = tmp_path / "merged.txt"

        result = subprocess.run(
            [
                sys.executable, "-m", "kmtools.kmtools",
                "--verbose",
                "merge",
                str(tmp_path / "chunk_0.txt"),
                str(tmp_path / "chunk_1.txt"),
                str(tmp_path / "chunk_2.txt"),
                "--output", str(out_file),
                "--keep",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert out_file.exists()

        df = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(df) == 3

    def test_filter_use_alt_via_cli(self, tmp_path):
        """Run filter with --use-alt through the CLI."""
        alt_ref = tmp_path / "alt_reference.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "filtered_alt.tsv"

        seq = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"

        write_alt_reference_tsv(alt_ref, [
            {"CHROM": "chr1", "ALT_SEQUENCE": seq, "TYPE": "Substitution"},
        ])

        write_km_output(km_file, [
            make_km_row(query="chr1_114716091_114716161", sequence=seq, min_coverage="10"),
        ])

        result = subprocess.run(
            [
                sys.executable, "-m", "kmtools.kmtools",
                "--verbose",
                "filter",
                "--reference", str(alt_ref),
                "--km-output", str(km_file),
                "--output", str(out_file),
                "--output-type", "tsv",
                "--use-alt",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert out_file.exists()

        df = pd.read_csv(out_file, sep="\t", dtype={"FOUND": str})
        assert len(df) == 1
        assert df.iloc[0]["FOUND"] == "TRUE"

    def test_filter_nonexistent_file_error(self, tmp_path):
        """CLI should exit with error when input file doesn't exist."""
        result = subprocess.run(
            [
                sys.executable, "-m", "kmtools.kmtools",
                "filter",
                "--reference", str(tmp_path / "missing.tsv"),
                "--km-output", str(tmp_path / "missing.txt"),
                "--output", str(tmp_path / "out.tsv"),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "not found" in result.stderr


# ---------------------------------------------------------------------------
# Integration: Error handling across boundaries
# ---------------------------------------------------------------------------


class TestErrorHandlingIntegration:
    """Tests for error handling across module boundaries."""

    def test_filter_with_malformed_reference(self, tmp_path):
        """Filter should raise ValueError for missing reference columns."""
        ref_file = tmp_path / "bad_ref.tsv"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "out.tsv"

        ref_file.write_text("COL_A\tCOL_B\nfoo\tbar\n")
        write_km_output(km_file, [make_km_row()])

        with pytest.raises(ValueError, match="missing required columns"):
            Filter(
                reference=str(ref_file),
                km_output=str(km_file),
                output=str(out_file),
                output_type="tsv",
                count_threshold=2,
            ).run()

    def test_filter_with_malformed_km_output(self, tmp_path):
        """Filter should raise ValueError for bad km output files."""
        ref_file = tmp_path / "reference.tsv"
        km_file = tmp_path / "bad_km.txt"
        out_file = tmp_path / "out.tsv"

        write_reference_tsv(ref_file, [
            {"CHROM": "chr1", "POS": "100", "REF": "C", "ALT": "T", "TYPE": "Substitution"},
        ])
        km_file.write_text("Bad\tColumns\nfoo\tbar\n")

        with pytest.raises(ValueError, match="missing required columns"):
            Filter(
                reference=str(ref_file),
                km_output=str(km_file),
                output=str(out_file),
                output_type="tsv",
                count_threshold=2,
            ).run()

    def test_filter_unsupported_reference_format(self, tmp_path):
        """Filter should raise ValueError for unsupported reference format."""
        ref_file = tmp_path / "reference.xlsx"
        km_file = tmp_path / "km_output.txt"
        out_file = tmp_path / "out.tsv"

        ref_file.write_text("dummy")
        write_km_output(km_file, [make_km_row()])

        with pytest.raises(ValueError, match="Unsupported file format"):
            Filter(
                reference=str(ref_file),
                km_output=str(km_file),
                output=str(out_file),
                output_type="tsv",
                count_threshold=2,
            ).run()

    def test_merge_empty_files(self, tmp_path, monkeypatch):
        """Merge files with only headers (no data rows)."""
        monkeypatch.chdir(tmp_path)

        for i in range(2):
            f = tmp_path / f"empty_{i}.txt"
            f.write_text(KM_OUTPUT_HEADER + "\n")

        output = tmp_path / "merged.txt"
        Merge(
            inputs=["empty_*.txt"],
            output=str(output),
            keep=True,
            verbose=True,
        ).run()

        result = pd.read_csv(output, sep="\t")
        assert len(result) == 0

    def test_plot_missing_file_raises(self, tmp_path):
        """Plot should raise FileNotFoundError for missing input."""
        with pytest.raises(FileNotFoundError):
            Plot(
                file=str(tmp_path / "nonexistent.tsv"),
                output_dir=str(tmp_path / "plots"),
                verbose=True,
            ).run()
