"""Data models and result containers for ParSe v2 predictions."""

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Dict, Any
import pandas as pd


@dataclass(frozen=True)
class ResidueResult:
    """Prediction values for an individual amino acid residue."""
    index: int  # 1-based residue number
    amino_acid: str
    label: Literal["F", "D", "P"]
    dist: float
    label_dh: Literal["F", "D", "P"]
    dist_dh: float
    label_csat: Literal["F", "D", "P"]
    dist_csat: float


@dataclass(frozen=True)
class DomainRegion:
    """Continuous domain region parsed from window classifications."""
    region_id: int  # 1-based domain index
    domain_type: Literal["F", "D", "P"]
    first_residue: int  # 1-based start residue
    last_residue: int  # 1-based end residue (inclusive)
    length: int


@dataclass(frozen=True)
class ParSeSummary:
    """Summary metrics for a protein sequence."""
    id: str
    length: int
    sum_p_dist: float
    sum_p_dist_dh: float
    sum_p_dist_csat: float
    num_ps_regions: int
    num_id_regions: int
    num_folded_regions: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "length": self.length,
            "sum_p_dist": self.sum_p_dist,
            "sum_p_dist_dh": self.sum_p_dist_dh,
            "sum_p_dist_csat": self.sum_p_dist_csat,
            "num_ps_regions": self.num_ps_regions,
            "num_id_regions": self.num_id_regions,
            "num_folded_regions": self.num_folded_regions,
        }


@dataclass
class ParSeResult:
    """Complete prediction output for a protein sequence."""
    id: str
    sequence: str
    summary: ParSeSummary
    residues: List[ResidueResult] = field(default_factory=list)
    domains: List[DomainRegion] = field(default_factory=list)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert per-residue prediction array to a pandas DataFrame."""
        return pd.DataFrame([
            {
                "protein_id": self.id,
                "residue_num": r.index,
                "amino_acid": r.amino_acid,
                "label": r.label,
                "dist": r.dist,
                "label_dh": r.label_dh,
                "dist_dh": r.dist_dh,
                "label_csat": r.label_csat,
                "dist_csat": r.dist_csat,
            }
            for r in self.residues
        ])

    def domains_to_dataframe(self) -> pd.DataFrame:
        """Convert segmented domain regions to a pandas DataFrame."""
        return pd.DataFrame([
            {
                "protein_id": self.id,
                "region_id": d.region_id,
                "domain_type": d.domain_type,
                "first_residue": d.first_residue,
                "last_residue": d.last_residue,
                "length": d.length,
            }
            for d in self.domains
        ])
