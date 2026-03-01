"""Tests for kmtools.filter_types."""

from kmtools.filter_types import (
    FilterCondition,
    FilterResult,
    KmVariant,
    TargetSequenceLocation,
)


class TestTargetSequenceLocation:
    def test_construction(self):
        loc = TargetSequenceLocation(chromosome="chr1", start=100, end=200)
        assert loc.chromosome == "chr1"
        assert loc.start == 100
        assert loc.end == 200


class TestKmVariant:
    def test_construction(self):
        variant = KmVariant(ref_pos=36, ref_allele="C", alt_pos=37, alt_allele="T")
        assert variant.ref_pos == 36
        assert variant.ref_allele == "C"
        assert variant.alt_pos == 37
        assert variant.alt_allele == "T"


class TestFilterCondition:
    def test_passing_condition(self):
        cond = FilterCondition(name="TYPE", condition=True, message="")
        assert cond.condition is True

    def test_failing_condition(self):
        cond = FilterCondition(name="TYPE", condition=False, message="Type mismatch")
        assert cond.condition is False
        assert cond.message == "Type mismatch"


class TestFilterResult:
    def test_passed(self):
        result = FilterResult(passed=True, failed_count="")
        assert result.passed is True
        assert result.failed_count == ""

    def test_failed(self):
        result = FilterResult(passed=False, failed_count="Insufficient count")
        assert result.passed is False
        assert result.failed_count == "Insufficient count"
