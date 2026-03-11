"""Shared fixtures for kmtools tests."""

import pytest
import pandas as pd


@pytest.fixture
def sample_reference_row():
    return {
        "CHROM": "chr1",
        "POS": "114716126",
        "REF": "C",
        "ALT": "T",
        "TYPE": "Substitution",
    }


@pytest.fixture
def sample_km_row():
    return {
        "Database": "file.jf",
        "Query": "chr1_114716091_114716161",
        "Type": "Substitution",
        "Variant_name": "36:c/T:37",
        "rVAF": "0.000",
        "Expression": "0.0",
        "Min_coverage": "2",
        "Start_offset": "0",
        "Sequence": "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
        "Reference_expression": "5033.0",
        "Reference_sequence": "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
        "Info": "vs_ref",
    }


@pytest.fixture
def sample_deletion_reference_row():
    return {
        "CHROM": "chr4",
        "POS": "54658073",
        "REF": "TC",
        "ALT": "T",
        "TYPE": "Deletion",
    }


@pytest.fixture
def sample_deletion_km_row():
    return {
        "Database": "file.jf",
        "Query": "chr4_54658038_54658108",
        "Type": "Deletion",
        "Variant_name": "35:tc/T:36",
        "rVAF": "0.000",
        "Expression": "0.0",
        "Min_coverage": "5",
        "Start_offset": "0",
        "Sequence": "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC",
        "Reference_expression": "3000.0",
        "Reference_sequence": "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC",
        "Info": "vs_ref",
    }


@pytest.fixture
def filter_instance():
    """Create a Filter instance with test data files."""
    from kmtools.filter import Filter

    return Filter(
        reference="tests/reference_1.tsv",
        km_output="tests/km_output_1.txt",
        output="output.csv",
        output_type="csv",
        count_threshold=2,
        verbose=False,
    )


@pytest.fixture
def km_output_tsv_content():
    """Standard km output TSV content for creating test files."""
    return (
        "Database\tQuery\tType\tVariant_name\trVAF\tExpression\t"
        "Min_coverage\tStart_offset\tSequence\tReference_expression\t"
        "Reference_sequence\tInfo\n"
        "file.jf\tchr1_114716091_114716161\tSubstitution\t36:c/T:37\t0.000\t0.0\t"
        "2\t0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\t"
        "5033.0\tTAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT\tvs_ref\n"
    )


@pytest.fixture
def reference_tsv_content():
    """Standard reference TSV content for creating test files."""
    return "CHROM\tPOS\tREF\tALT\tTYPE\nchr1\t114716126\tC\tT\tSubstitution\n"
