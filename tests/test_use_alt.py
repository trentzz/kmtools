"""Tests for --use-alt filter mode."""

import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch
import sys

from kmtools.filter import Filter
from kmtools.filter_types import FilterResult


# === Fixtures ===


@pytest.fixture
def alt_reference_content():
    return "CHROM\tALT_SEQUENCE\tTYPE\nchr1\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\tSubstitution\n"


@pytest.fixture
def km_output_content():
    return (
        "Database\tQuery\tType\tVariant_name\trVAF\tExpression\t"
        "Min_coverage\tStart_offset\tSequence\tReference_expression\t"
        "Reference_sequence\tInfo\n"
        "file.jf\tchr1_114716092_114716162\tSubstitution\t36:c/T:37\t0.000\t0.0\t"
        "2\t0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\t"
        "5033.0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\tvs_ref\n"
    )


@pytest.fixture
def alt_filter_instance(tmp_path, alt_reference_content, km_output_content):
    ref_file = tmp_path / "alt_ref.tsv"
    ref_file.write_text(alt_reference_content)
    km_file = tmp_path / "km_output.txt"
    km_file.write_text(km_output_content)

    return Filter(
        reference=str(ref_file),
        km_output=str(km_file),
        output=str(tmp_path / "output.tsv"),
        output_type="tsv",
        count_threshold=2,
        use_alt=True,
        verbose=False,
    )


# === Tests for verify_alt_reference ===


class TestVerifyAltReference:
    def test_valid_alt_reference(self, tmp_path, alt_reference_content):
        ref_file = tmp_path / "alt_ref.tsv"
        ref_file.write_text(alt_reference_content)
        f = Filter(
            reference=str(ref_file),
            km_output="dummy.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        f.verify_alt_reference()
        assert f.reference_df is not None
        assert "ALT_SEQUENCE" in f.reference_df.columns
        assert len(f.reference_df) == 1

    def test_valid_alt_reference_csv(self, tmp_path):
        ref_file = tmp_path / "alt_ref.csv"
        ref_file.write_text("CHROM,ALT_SEQUENCE,TYPE\nchr1,ATCG,Substitution\n")
        f = Filter(
            reference=str(ref_file),
            km_output="dummy.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        f.verify_alt_reference()
        assert f.reference_df is not None

    def test_missing_alt_sequence_column(self, tmp_path):
        ref_file = tmp_path / "bad_ref.tsv"
        ref_file.write_text("CHROM\tTYPE\nchr1\tSubstitution\n")
        f = Filter(
            reference=str(ref_file),
            km_output="dummy.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        with pytest.raises(ValueError, match="ALT_SEQUENCE"):
            f.verify_alt_reference()

    def test_wrong_format_for_use_alt(self, tmp_path):
        """Using a standard reference file (with POS/REF/ALT) with --use-alt should fail."""
        ref_file = tmp_path / "standard_ref.tsv"
        ref_file.write_text("CHROM\tPOS\tREF\tALT\tTYPE\nchr1\t100\tA\tT\tSNP\n")
        f = Filter(
            reference=str(ref_file),
            km_output="dummy.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        with pytest.raises(ValueError, match="ALT_SEQUENCE"):
            f.verify_alt_reference()


# === Tests for filter_alt_line ===


class TestFilterAltLine:
    def test_all_conditions_pass(self, alt_filter_instance):
        km_row = {
            "Type": "Substitution",
            "Min_coverage": "5",
            "Info": "vs_ref",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "ATCG"
        )
        assert result.passed is True

    def test_type_mismatch(self, alt_filter_instance):
        km_row = {
            "Type": "Deletion",
            "Min_coverage": "5",
            "Info": "vs_ref",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "ATCG"
        )
        assert result.passed is False

    def test_sequence_mismatch(self, alt_filter_instance):
        km_row = {
            "Type": "Substitution",
            "Min_coverage": "5",
            "Info": "vs_ref",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "GCTA"
        )
        assert result.passed is False

    def test_count_below_threshold(self, alt_filter_instance):
        km_row = {
            "Type": "Substitution",
            "Min_coverage": "1",
            "Info": "vs_ref",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "ATCG"
        )
        assert result.passed is False
        assert "Insufficient" in result.failed_count

    def test_info_not_vs_ref(self, alt_filter_instance):
        km_row = {
            "Type": "Substitution",
            "Min_coverage": "5",
            "Info": "something_else",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "ATCG"
        )
        assert result.passed is False

    @pytest.mark.parametrize(
        "min_cov,expected_pass",
        [("1", False), ("2", True), ("10", True)],
    )
    def test_count_threshold_boundary(self, alt_filter_instance, min_cov, expected_pass):
        km_row = {
            "Type": "Substitution",
            "Min_coverage": min_cov,
            "Info": "vs_ref",
        }
        alt_row = {"TYPE": "Substitution"}
        result = alt_filter_instance.filter_alt_line(
            km_row, alt_row, "ATCG", "ATCG"
        )
        if expected_pass:
            assert result.passed is True
        else:
            assert result.passed is False


# === Tests for write_alt_filtered_line ===


class TestWriteAltFilteredLine:
    def test_passed_result(self, alt_filter_instance):
        alt_row = {
            "CHROM": "chr1",
            "ALT_SEQUENCE": "ATCG",
            "TYPE": "Substitution",
        }
        km_row = {
            "rVAF": "0.5",
            "Min_coverage": "10",
            "Expression": "100.0",
            "Reference_sequence": "TTCG",
            "Sequence": "ATCG",
        }
        result = FilterResult(passed=True, failed_count="")
        alt_filter_instance.write_alt_filtered_line(alt_row, km_row, result)

        assert len(alt_filter_instance.output_df) == 1
        row = alt_filter_instance.output_df[0]
        assert row["FOUND"] == "TRUE"
        assert row["ALT_SEQUENCE"] == "ATCG"
        assert row["KMER_VAF"] == "0.5"
        assert "POS" not in row
        assert "REF" not in row

    def test_failed_result(self, alt_filter_instance):
        alt_row = {
            "CHROM": "chr1",
            "ALT_SEQUENCE": "ATCG",
            "TYPE": "Substitution",
        }
        km_row = {
            "rVAF": "0.5",
            "Min_coverage": "1",
            "Expression": "100.0",
            "Reference_sequence": "TTCG",
            "Sequence": "ATCG",
        }
        result = FilterResult(passed=False, failed_count="Count too low")
        alt_filter_instance.write_alt_filtered_line(alt_row, km_row, result)

        row = alt_filter_instance.output_df[0]
        assert row["FOUND"] == "FALSE"
        assert row["KMER_VAF"] == ""
        assert row["FILTER_NOTES"] == "Count too low"

    def test_none_result_skipped(self, alt_filter_instance):
        alt_filter_instance.write_alt_filtered_line({}, {}, None)
        assert len(alt_filter_instance.output_df) == 0


# === Tests for full run_alt_filtering ===


class TestRunAltFiltering:
    def test_end_to_end_match(self, alt_filter_instance):
        """Test that the full alt filtering pipeline finds a matching sequence."""
        alt_filter_instance.verify_alt_reference()
        alt_filter_instance.verify_km_output()
        alt_filter_instance.run_alt_filtering()

        assert len(alt_filter_instance.output_df) == 1
        row = alt_filter_instance.output_df[0]
        assert row["FOUND"] == "TRUE"
        assert row["CHROM"] == "chr1"

    def test_no_match_wrong_sequence(self, tmp_path, km_output_content):
        ref_file = tmp_path / "alt_ref.tsv"
        ref_file.write_text("CHROM\tALT_SEQUENCE\tTYPE\nchr1\tNOMATCHSEQUENCE\tSubstitution\n")
        km_file = tmp_path / "km_output.txt"
        km_file.write_text(km_output_content)

        f = Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(tmp_path / "output.tsv"),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        f.verify_alt_reference()
        f.verify_km_output()
        f.run_alt_filtering()

        assert len(f.output_df) == 1
        assert f.output_df[0]["FOUND"] == "FALSE"

    def test_no_match_wrong_chromosome(self, tmp_path, km_output_content):
        ref_file = tmp_path / "alt_ref.tsv"
        ref_file.write_text(
            "CHROM\tALT_SEQUENCE\tTYPE\n"
            "chr99\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\tSubstitution\n"
        )
        km_file = tmp_path / "km_output.txt"
        km_file.write_text(km_output_content)

        f = Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(tmp_path / "output.tsv"),
            output_type="tsv",
            count_threshold=2,
            use_alt=True,
        )
        f.verify_alt_reference()
        f.verify_km_output()
        f.run_alt_filtering()

        assert len(f.output_df) == 1
        assert f.output_df[0]["FOUND"] == "FALSE"

    def test_write_output_tsv(self, tmp_path, alt_filter_instance):
        alt_filter_instance.verify_alt_reference()
        alt_filter_instance.verify_km_output()
        alt_filter_instance.run_alt_filtering()
        alt_filter_instance.write_output()

        output_path = Path(alt_filter_instance.output)
        assert output_path.exists()
        result = pd.read_csv(output_path, sep="\t")
        assert len(result) == 1
        assert "ALT_SEQUENCE" in result.columns
        assert "FOUND" in result.columns


# === Tests for run() dispatch ===


class TestRunDispatch:
    def test_run_dispatches_to_alt_mode(self, alt_filter_instance):
        """Verify that run() calls alt filtering when use_alt=True."""
        alt_filter_instance.run()

        output_path = Path(alt_filter_instance.output)
        assert output_path.exists()
        result = pd.read_csv(output_path, sep="\t")
        assert "ALT_SEQUENCE" in result.columns

    def test_run_dispatches_to_reference_mode(self, tmp_path):
        """Verify that run() calls reference filtering when use_alt=False."""
        ref_file = tmp_path / "ref.tsv"
        ref_file.write_text("CHROM\tPOS\tREF\tALT\tTYPE\nchr1\t114716126\tC\tT\tSubstitution\n")
        km_file = tmp_path / "km_output.txt"
        km_file.write_text(
            "Database\tQuery\tType\tVariant_name\trVAF\tExpression\t"
            "Min_coverage\tStart_offset\tSequence\tReference_expression\t"
            "Reference_sequence\tInfo\n"
            "file.jf\tchr1_114716091_114716161\tSubstitution\t36:c/T:37\t0.000\t0.0\t"
            "2\t0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\t"
            "5033.0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\tvs_ref\n"
        )

        f = Filter(
            reference=str(ref_file),
            km_output=str(km_file),
            output=str(tmp_path / "output.tsv"),
            output_type="tsv",
            count_threshold=2,
            use_alt=False,
        )
        f.run()

        result = pd.read_csv(tmp_path / "output.tsv", sep="\t")
        assert "POS" in result.columns
        assert "ALT_SEQUENCE" not in result.columns


# === CLI test ===


class TestUseAltCLI:
    def test_filter_parser_accepts_use_alt(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "filter",
            "--reference",
            "ref.tsv",
            "--km-output",
            "km.txt",
            "--output",
            "filtered.tsv",
            "--use-alt",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_filter") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.use_alt is True

    def test_filter_parser_default_no_use_alt(self):
        from kmtools.kmtools import main

        test_args = [
            "kmtools",
            "filter",
            "--reference",
            "ref.tsv",
            "--km-output",
            "km.txt",
            "--output",
            "filtered.tsv",
        ]
        with patch.object(sys, "argv", test_args):
            with patch("kmtools.kmtools.run_filter") as mock_run:
                main()
                args = mock_run.call_args[0][0]
                assert args.use_alt is False
