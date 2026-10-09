"""ParSe v2: Prediction of protein liquid-liquid phase separation potential.

Based on Whitten et al., Protein Science 32, e4756 (2023).
"""

from parse_v2.core import parse_sequence, clean_sequence
from parse_v2.parallel import parse_fasta
from parse_v2.models import ParSeResult, ParSeSummary, ResidueResult, DomainRegion
from parse_v2.io import (
    stream_fasta,
    export_summary,
    export_residues,
    export_domains,
    export_filtered_fasta,
)

__version__ = "0.1.0"

__all__ = [
    "parse_sequence",
    "clean_sequence",
    "parse_fasta",
    "ParSeResult",
    "ParSeSummary",
    "ResidueResult",
    "DomainRegion",
    "stream_fasta",
    "export_summary",
    "export_residues",
    "export_domains",
    "export_filtered_fasta",
]
