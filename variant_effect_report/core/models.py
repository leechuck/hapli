"""
Data models for variant effect analysis.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Union
from enum import Enum


class VariantType(Enum):
    """Enumeration of variant types."""
    SNP = "SNP"
    SNV = "SNV"
    DEL = "DEL"
    INS = "INS"
    DUP = "DUP"
    INV = "INV"
    UNKNOWN = "UNKNOWN"


@dataclass
class Segment:
    """Represents a segment in a GFA file."""
    id: str
    sequence: str
    length: int
    tags: Dict[str, Tuple[str, str]] = field(default_factory=dict)
    variant_id: Optional[str] = None
    variant_type: Optional[VariantType] = None
    length_change: int = 0
    original_segment: Optional[str] = None


@dataclass
class Link:
    """Represents a link between segments in a GFA file."""
    from_id: str
    from_dir: str
    to_id: str
    to_dir: str
    overlap: str


@dataclass
class Path:
    """Represents a path in a GFA file."""
    name: str
    segments: List[Tuple[str, str]]  # (segment_id, orientation)
    tags: Dict[str, Tuple[str, str]] = field(default_factory=dict)
    sample: Optional[str] = None
    haplotype: Optional[str] = None
    variant_ids: List[str] = field(default_factory=list)


@dataclass
class Variant:
    """Represents a genomic variant."""
    id: str
    type: VariantType = VariantType.UNKNOWN
    segments: List[str] = field(default_factory=list)
    length_change: int = 0
    pos: Optional[int] = None
    end: Optional[int] = None
    paths: List[str] = field(default_factory=list)
    samples: List[str] = field(default_factory=list)
    repeat_info: Optional[Dict] = None


@dataclass
class FeatureAttribute:
    """Represents an attribute of a genomic feature."""
    key: str
    value: str


@dataclass
class Feature:
    """Represents a genomic feature from a GFF3 file."""
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str] = field(default_factory=dict)
    children: List['Feature'] = field(default_factory=list)
    is_first_cds: bool = False
    line_num: Optional[int] = None


@dataclass
class SegmentOffset:
    """Represents the offset of a segment in a path."""
    seg_id: str
    orientation: str
    start: int
    end: int
    length: int
    variant_id: Optional[str] = None


@dataclass
class AlignmentResult:
    """Represents the result of a sequence alignment."""
    alignment_score: float
    ref_aligned: str
    alt_aligned: str
    insertions: List[Dict]
    deletions: List[Dict]
    substitutions: List[Dict]
    inversions: List[Dict]
    complex_regions: List[Dict]
    length_change: int


@dataclass
class FeatureEffect:
    """Represents the effect of a variant on a feature."""
    feature: Feature
    feature_type: str
    ref_feature_seq: str = ""
    alt_feature_seq: str = ""
    effects: List[str] = field(default_factory=list)
    details: Dict = field(default_factory=dict)
    variants: List[Variant] = field(default_factory=list)
    zygosity: Optional[str] = None
    haplotype_effects: Optional[Dict[str, 'FeatureEffect']] = None
    combined_effects: Optional[List[str]] = None


@dataclass
class Sample:
    """Represents a sample with haplotypes."""
    name: str
    paths: List[str] = field(default_factory=list)
    haplotypes: Dict[str, str] = field(default_factory=dict)
    effects: Dict = field(default_factory=dict)
