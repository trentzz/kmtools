"""Comprehensive tests for kmtools.filter."""

import pytest
import pandas as pd
from pathlib import Path

from kmtools.filter import Filter
from kmtools.filter_types import KmVariant, TargetSequenceLocation


class TestSplitQuery:
    def test_standard_query(self, filter_instance):
        result = filter_instance.split_query("chr1_114716091_114716161")
        assert result.chromosome == "chr1"
        assert result.start == 114716091
        assert result.end == 114716161

    def test_different_chromosome(self, filter_instance):
        result = filter_instance.split_query("chr4_54658038_54658108")
        assert result.chromosome == "chr4"
        assert result.start == 54658038
        assert result.end == 54658108


class TestGetRefAltPosFromVariant:
    def test_substitution(self, filter_instance):
        result = filter_instance.get_ref_alt_pos_from_variant("36:c/T:37")
        assert result.ref_pos == 36
        assert result.ref_allele == "C"
        assert result.alt_pos == 37
        assert result.alt_allele == "T"

    def test_empty_string(self, filter_instance):
        result = filter_instance.get_ref_alt_pos_from_variant("")
        assert result is None

    def test_none_input(self, filter_instance):
        result = filter_instance.get_ref_alt_pos_from_variant(None)
        assert result is None


class TestGetKmAlt:
    def test_returns_uppercase_sequence(self, filter_instance, sample_km_row):
        result = filter_instance.get_km_alt(sample_km_row)
        assert result == result.upper()
        assert result == sample_km_row["Sequence"].upper()

    def test_lowercase_input(self, filter_instance):
        km_row = {"Sequence": "atcg"}
        assert filter_instance.get_km_alt(km_row) == "ATCG"


class TestGetCalculatedReferenceAlt:
    def test_substitution(self, filter_instance, sample_reference_row, sample_km_row):
        result = filter_instance.get_calculated_reference_alt(
            sample_reference_row, sample_km_row
        )
        expected = "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"
        assert result == expected

    def test_deletion(
        self, filter_instance, sample_deletion_reference_row, sample_deletion_km_row
    ):
        result = filter_instance.get_calculated_reference_alt(
            sample_deletion_reference_row, sample_deletion_km_row
        )
        expected = "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC"
        assert result == expected

    def test_pos_outside_range(self, filter_instance):
        reference_row = {
            "CHROM": "chr1",
            "POS": "999999999",
            "REF": "C",
            "ALT": "T",
            "TYPE": "Substitution",
        }
        km_row = {
            "Query": "chr1_100_200",
            "Reference_sequence": "ATCGATCG",
        }
        result = filter_instance.get_calculated_reference_alt(reference_row, km_row)
        assert result == ""


class TestFilterLine:
    def test_all_conditions_pass(
        self, filter_instance, sample_km_row, sample_reference_row
    ):
        # Make reference TYPE match km Type
        sample_reference_row["TYPE"] = "Substitution"
        km_alt = "ATCG"
        calculated_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, calculated_alt
        )
        assert result.passed is True

    def test_type_mismatch(self, filter_instance, sample_km_row, sample_reference_row):
        sample_reference_row["TYPE"] = "Deletion"
        km_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, km_alt
        )
        assert result.passed is False

    def test_count_below_threshold(
        self, filter_instance, sample_km_row, sample_reference_row
    ):
        sample_reference_row["TYPE"] = "Substitution"
        sample_km_row["Min_coverage"] = "1"
        km_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, km_alt
        )
        assert result.passed is False
        assert "Insufficient" in result.failed_count

    def test_alt_mismatch(self, filter_instance, sample_km_row, sample_reference_row):
        sample_reference_row["TYPE"] = "Substitution"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, "ATCG", "GCTA"
        )
        assert result.passed is False

    def test_reference_type_skipped(
        self, filter_instance, sample_km_row, sample_reference_row
    ):
        sample_reference_row["TYPE"] = "Reference"
        km_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, km_alt
        )
        assert result.passed is False

    def test_info_not_vs_ref(self, filter_instance, sample_km_row, sample_reference_row):
        sample_reference_row["TYPE"] = "Substitution"
        sample_km_row["Info"] = "other_info"
        km_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, km_alt
        )
        assert result.passed is False

    @pytest.mark.parametrize(
        "min_cov,expected_pass",
        [("1", False), ("2", True), ("3", True), ("100", True)],
    )
    def test_count_threshold_boundary(
        self, filter_instance, sample_km_row, sample_reference_row, min_cov, expected_pass
    ):
        sample_reference_row["TYPE"] = "Substitution"
        sample_km_row["Min_coverage"] = min_cov
        km_alt = "ATCG"
        result = filter_instance.filter_line(
            sample_km_row, sample_reference_row, km_alt, km_alt
        )
        # When only COUNT fails, passed is False but other conditions pass
        if not expected_pass:
            assert result.passed is False
        else:
            assert result.passed is True


class TestVerifyReference:
    def test_valid_tsv(self, tmp_path):
        ref_file = tmp_path / "ref.tsv"
        ref_file.write_text("CHROM\tPOS\tREF\tALT\tTYPE\nchr1\t100\tA\tT\tSNP\n")
        f = Filter(
            reference=str(ref_file),
            km_output="tests/km_output_1.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        f.verify_reference()
        assert f.reference_df is not None
        assert len(f.reference_df) == 1

    def test_valid_csv(self, tmp_path):
        ref_file = tmp_path / "ref.csv"
        ref_file.write_text("CHROM,POS,REF,ALT,TYPE\nchr1,100,A,T,SNP\n")
        f = Filter(
            reference=str(ref_file),
            km_output="tests/km_output_1.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        f.verify_reference()
        assert f.reference_df is not None

    def test_missing_columns(self, tmp_path):
        ref_file = tmp_path / "ref.tsv"
        ref_file.write_text("CHROM\tPOS\nchr1\t100\n")
        f = Filter(
            reference=str(ref_file),
            km_output="tests/km_output_1.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        with pytest.raises(ValueError, match="missing required columns"):
            f.verify_reference()

    def test_unsupported_format(self, tmp_path):
        ref_file = tmp_path / "ref.json"
        ref_file.write_text("{}")
        f = Filter(
            reference=str(ref_file),
            km_output="tests/km_output_1.txt",
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        with pytest.raises(ValueError, match="Unsupported file format"):
            f.verify_reference()


class TestVerifyKmOutput:
    def test_valid_km_output(self, tmp_path, km_output_tsv_content):
        km_file = tmp_path / "km.txt"
        km_file.write_text(km_output_tsv_content)
        f = Filter(
            reference="tests/reference_1.tsv",
            km_output=str(km_file),
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        f.verify_km_output()
        assert f.km_output_df is not None

    def test_missing_columns(self, tmp_path):
        km_file = tmp_path / "km.txt"
        km_file.write_text("Database\tQuery\nfile.jf\tchr1_100_200\n")
        f = Filter(
            reference="tests/reference_1.tsv",
            km_output=str(km_file),
            output="out.tsv",
            output_type="tsv",
            count_threshold=2,
        )
        with pytest.raises(ValueError, match="missing required columns"):
            f.verify_km_output()


class TestWriteOutput:
    def test_tsv_output(self, tmp_path, filter_instance):
        filter_instance.output = str(tmp_path / "out.tsv")
        filter_instance.output_type = "tsv"
        filter_instance.output_df = [
            {"SAMPLE": "test", "CHROM": "chr1", "FOUND": "TRUE"}
        ]
        filter_instance.write_output()
        result = pd.read_csv(tmp_path / "out.tsv", sep="\t")
        assert len(result) == 1
        assert result.iloc[0]["SAMPLE"] == "test"

    def test_csv_output(self, tmp_path, filter_instance):
        filter_instance.output = str(tmp_path / "out.csv")
        filter_instance.output_type = "csv"
        filter_instance.output_df = [
            {"SAMPLE": "test", "CHROM": "chr1", "FOUND": "TRUE"}
        ]
        filter_instance.write_output()
        result = pd.read_csv(tmp_path / "out.csv")
        assert len(result) == 1


class TestWriteFilteredLine:
    def test_passed_filter(self, filter_instance, sample_reference_row, sample_km_row):
        from kmtools.filter_types import FilterResult

        result = FilterResult(passed=True, failed_count="")
        filter_instance.write_filtered_line(sample_reference_row, sample_km_row, result)
        assert len(filter_instance.output_df) == 1
        row = filter_instance.output_df[0]
        assert row["FOUND"] == "TRUE"
        assert row["KMER_VAF"] == sample_km_row["rVAF"]

    def test_failed_filter(self, filter_instance, sample_reference_row, sample_km_row):
        from kmtools.filter_types import FilterResult

        result = FilterResult(passed=False, failed_count="Count too low")
        filter_instance.write_filtered_line(sample_reference_row, sample_km_row, result)
        assert len(filter_instance.output_df) == 1
        row = filter_instance.output_df[0]
        assert row["FOUND"] == "FALSE"
        assert row["KMER_VAF"] == ""
        assert row["FILTER_NOTES"] == "Count too low"

    def test_none_result_skipped(
        self, filter_instance, sample_reference_row, sample_km_row
    ):
        filter_instance.write_filtered_line(
            sample_reference_row, sample_km_row, None
        )
        assert len(filter_instance.output_df) == 0
