
from dataclasses import dataclass


@dataclass
class TargetSequenceLocation:
    chromosome: str
    start: int
    end: int
    
@dataclass
class KmVariant:
    ref_pos: int
    ref_allele: str
    alt_pos: int
    alt_allele: str
    
@dataclass
class FilterCondition:
    name: str
    condition: bool
    message: str
    
@dataclass
class FilterResult:
    passed: bool
    failed_count: str