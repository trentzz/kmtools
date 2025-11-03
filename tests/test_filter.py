import pytest

from kmtools.filter import Filter
from kmtools.filter_types import KmVariant, TargetSequenceLocation


def test_get_calculated_reference_alt():
    filter_instance = Filter(
        reference="tests/reference_1.tsv",
        km_output="tests/km_output_1.txt",
        output="output.csv",
        output_type="csv",
        count_threshold=2,
        verbose=False,
    )

    reference_row = {
        "CHROM": "chr1",
        "POS": "114716126",
        "REF": "C",
        "ALT": "T",
        "TYPE": "Substitution",
    }

    km_row = {
        "Query": "chr1_114716091_114716161",
        "Sequence": "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
        "Reference_sequence": "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACCTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT",
    }

    calculated_alt = filter_instance.get_calculated_reference_alt(reference_row, km_row)
    assert (
        calculated_alt
        == "TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT"
    )


def test_get_calculated_reference_alt_2():
    filter_instance = Filter(
        reference="tests/reference_1.tsv",
        km_output="tests/km_output_1.txt",
        output="output.csv",
        output_type="csv",
        count_threshold=2,
        verbose=False,
    )

    reference_row = {
        "CHROM": "chr4",
        "POS": "54658073",
        "REF": "TC",
        "ALT": "T",
        "TYPE": "Deletion",
    }

    km_row = {
        "Query": "chr4_54658038_54658108",
        "Sequence": "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC",
        "Reference_sequence": "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC",
    }

    calculated_alt = filter_instance.get_calculated_reference_alt(reference_row, km_row)
    assert (
        calculated_alt
        == "GGATTTTCTCTGCGTTCTGCTCCTACTGCTTCGCGTCAGACAGGTGGGACACCGCGGCTGGCACCCCGAC"
    )
