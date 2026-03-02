"""
Integration tests using realistic test data files in tests/data/.

These tests exercise the filter, merge, and plot pipeline against
biologically plausible data files with internally consistent sequences.
"""

from __future__ import annotations

from pathlib import Path
import shutil

import pandas as pd
import pytest

from kmtools.filter import Filter
from kmtools.merge import Merge
from kmtools.plot import Plot

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Substitution data files
# ---------------------------------------------------------------------------


class TestSubstitutionDataFiles:
    """Filter with substitution reference and km output data files."""

    def test_filter_substitution_default_threshold(self, tmp_path):
        """Filter substitution data with default count threshold (2)."""
        out = tmp_path / "filtered.tsv"

        Filter(
            reference=str(DATA_DIR / "reference_substitution.tsv"),
            km_output=str(DATA_DIR / "km_output_substitution.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        assert len(result) == 3  # 3 reference variants
        # chr1 sub should be found (row 1 has coverage 12)
        chr1_row = result[result["CHROM"] == "chr1"].iloc[0]
        assert chr1_row["FOUND"] == "TRUE"
        assert float(chr1_row["KMER_VAF"]) == pytest.approx(0.45)
        # chr2:25457242 should be found (coverage 8)
        chr2_rows = result[result["CHROM"] == "chr2"]
        assert len(chr2_rows) == 2
        found_count = (chr2_rows["FOUND"] == "TRUE").sum()
        assert found_count == 2

    def test_filter_substitution_high_threshold(self, tmp_path):
        """Higher count threshold should reject lower-coverage variants."""
        out = tmp_path / "filtered.tsv"

        Filter(
            reference=str(DATA_DIR / "reference_substitution.tsv"),
            km_output=str(DATA_DIR / "km_output_substitution.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=10,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        # chr1 has coverage 12, should still pass
        chr1_row = result[result["CHROM"] == "chr1"].iloc[0]
        assert chr1_row["FOUND"] == "TRUE"
        # chr2:25457242 has coverage 8, should fail with threshold 10
        chr2a = result[(result["CHROM"] == "chr2") & (result["POS"].astype(str) == "25457242")]
        if len(chr2a) > 0:
            row = chr2a.iloc[0]
            assert row["FOUND"] == "FALSE"
            assert "Insufficient KM count" in str(row["FILTER_NOTES"])


# ---------------------------------------------------------------------------
# Deletion data files
# ---------------------------------------------------------------------------


class TestDeletionDataFiles:
    """Filter with deletion reference and km output data files."""

    def test_filter_deletion_default_threshold(self, tmp_path):
        """Filter deletion data with default count threshold."""
        out = tmp_path / "filtered.tsv"

        Filter(
            reference=str(DATA_DIR / "reference_deletion.tsv"),
            km_output=str(DATA_DIR / "km_output_deletion.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        assert len(result) == 2  # 2 deletion reference variants
        # chr3 del should be found (coverage 7)
        chr3_row = result[result["CHROM"] == "chr3"].iloc[0]
        assert chr3_row["FOUND"] == "TRUE"
        assert chr3_row["TYPE"] == "Deletion"
        # chr4 del should be found (coverage 5)
        chr4_row = result[result["CHROM"] == "chr4"].iloc[0]
        assert chr4_row["FOUND"] == "TRUE"

    def test_filter_deletion_high_threshold(self, tmp_path):
        """Higher threshold causes low-coverage deletions to fail."""
        out = tmp_path / "filtered.tsv"

        Filter(
            reference=str(DATA_DIR / "reference_deletion.tsv"),
            km_output=str(DATA_DIR / "km_output_deletion.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=6,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        # chr3 should pass (coverage 7 >= 6)
        chr3_row = result[result["CHROM"] == "chr3"].iloc[0]
        assert chr3_row["FOUND"] == "TRUE"
        # chr4 should fail (coverage 5 < 6)
        chr4_row = result[result["CHROM"] == "chr4"].iloc[0]
        assert chr4_row["FOUND"] == "FALSE"


# ---------------------------------------------------------------------------
# Mixed variant type data files
# ---------------------------------------------------------------------------


class TestMixedDataFiles:
    """Filter with mixed variant types (substitution, deletion, insertion)."""

    def test_filter_mixed_variants(self, tmp_path):
        """Filter mixed variant types against mixed km output."""
        out = tmp_path / "filtered.tsv"

        Filter(
            reference=str(DATA_DIR / "reference_mixed.tsv"),
            km_output=str(DATA_DIR / "km_output_mixed.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        assert len(result) == 4  # 4 reference variants

        # All 4 should be found with threshold 2
        found = result[result["FOUND"] == "TRUE"]
        assert len(found) == 4

        # Verify variant types
        types = set(result["TYPE"].values)
        assert "Substitution" in types
        assert "Deletion" in types
        assert "Insertion" in types

    def test_filter_mixed_csv_output(self, tmp_path):
        """Verify CSV output format with mixed data."""
        out = tmp_path / "filtered.csv"

        Filter(
            reference=str(DATA_DIR / "reference_mixed.tsv"),
            km_output=str(DATA_DIR / "km_output_mixed.txt"),
            output=str(out),
            output_type="csv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out)
        assert len(result) == 4
        assert "SAMPLE" in result.columns
        assert "FOUND" in result.columns


# ---------------------------------------------------------------------------
# Use-alt mode data files
# ---------------------------------------------------------------------------


class TestUseAltDataFiles:
    """Filter in use-alt mode with data files."""

    def test_filter_use_alt_mode(self, tmp_path):
        """Use-alt mode matches km sequences against ALT_SEQUENCE."""
        out = tmp_path / "filtered_alt.tsv"

        Filter(
            reference=str(DATA_DIR / "alt_reference.tsv"),
            km_output=str(DATA_DIR / "km_output_mixed.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        assert len(result) == 3  # 3 alt reference entries
        assert "ALT_SEQUENCE" in result.columns
        # chr1 should match (substitution sequence matches)
        chr1_row = result[result["CHROM"] == "chr1"].iloc[0]
        assert chr1_row["FOUND"] == "TRUE"

    def test_filter_use_alt_output_columns(self, tmp_path):
        """Verify use-alt mode produces correct output columns."""
        out = tmp_path / "filtered_alt.tsv"

        Filter(
            reference=str(DATA_DIR / "alt_reference.tsv"),
            km_output=str(DATA_DIR / "km_output_mixed.txt"),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        expected = {
            "SAMPLE", "CHROM", "ALT_SEQUENCE", "TYPE", "FOUND",
            "FILTER_NOTES", "KMER_VAF", "KMER_MIN_COVERAGE",
            "KMER_EXPRESSION", "REF_SEQUENCE", "VARIANT_SEQUENCE",
        }
        assert set(result.columns) == expected


# ---------------------------------------------------------------------------
# Merge then filter pipeline
# ---------------------------------------------------------------------------


class TestMergeThenFilterDataFiles:
    """Merge multiple km outputs, then filter the merged result."""

    def test_merge_substitution_then_filter(self, tmp_path):
        """Merge copies of substitution data then filter."""
        # Copy km_output_substitution.txt to multiple files
        for i in range(3):
            shutil.copy(
                DATA_DIR / "km_output_substitution.txt",
                tmp_path / f"chunk_{i}.txt",
            )

        merged = tmp_path / "merged.txt"
        Merge(
            inputs=[str(tmp_path / f"chunk_{i}.txt") for i in range(3)],
            output=str(merged),
            keep=True,
            verbose=True,
        ).run()

        merged_df = pd.read_csv(merged, sep="\t")
        assert len(merged_df) == 18  # 6 rows × 3 files

        # Filter the merged output
        out = tmp_path / "filtered.tsv"
        Filter(
            reference=str(DATA_DIR / "reference_substitution.tsv"),
            km_output=str(merged),
            output=str(out),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(out, sep="\t", dtype={"FOUND": str})
        assert len(result) == 3  # 3 reference variants


# ---------------------------------------------------------------------------
# Full pipeline: merge → filter → plot
# ---------------------------------------------------------------------------


class TestFullPipelineDataFiles:
    """End-to-end pipeline using test data files."""

    def test_merge_filter_plot_pipeline(self, tmp_path):
        """Full pipeline: merge → filter → plot with mixed data."""
        # Copy mixed km output to simulate chunk files
        for i in range(2):
            shutil.copy(
                DATA_DIR / "km_output_mixed.txt",
                tmp_path / f"chunk_{i}.txt",
            )

        # Step 1: Merge
        merged = tmp_path / "merged.txt"
        Merge(
            inputs=[str(tmp_path / f"chunk_{i}.txt") for i in range(2)],
            output=str(merged),
            keep=True,
            verbose=True,
        ).run()

        assert merged.exists()
        merged_df = pd.read_csv(merged, sep="\t")
        assert len(merged_df) == 20  # 10 rows × 2 files

        # Step 2: Filter
        filtered = tmp_path / "filtered.tsv"
        Filter(
            reference=str(DATA_DIR / "reference_mixed.tsv"),
            km_output=str(merged),
            output=str(filtered),
            output_type="tsv",
            count_threshold=2,
            verbose=True,
        ).run()

        result = pd.read_csv(filtered, sep="\t", dtype={"FOUND": str})
        assert len(result) == 4
        found_count = (result["FOUND"] == "TRUE").sum()
        assert found_count >= 1

        # Step 3: Plot
        plot_dir = tmp_path / "plots"
        Plot(
            file=str(filtered),
            output_dir=str(plot_dir),
            charts="all",
            verbose=True,
        ).run()

        assert (plot_dir / "overall_summary.png").exists()
        assert (plot_dir / "type_distribution.png").exists()
