"""Vectorized calculation engine for ParSe v2."""

import math
from typing import List, Literal, Optional, Tuple
import numpy as np

from parse_v2.constants import (
    STANDARD_AMINO_ACIDS,
    HYDROPHOBICITY,
    PPII_PROPENSITY,
    HELIX_PROPENSITY,
    CHARGE,
    WINDOW_SIZE,
    HALF_WINDOW,
    FOLDED_HYDROPHOBICITY_CUTOFF,
    FOLDED_NORMALIZATION_FACTOR,
    BOUNDARY_SLOPE,
    BOUNDARY_INTERCEPT,
    PERP_SLOPE,
    ID_REF_DIST,
    PS_REF_DIST,
    DEFAULT_MIN_REGION_LEN,
    DEFAULT_PURITY_CUTOFF,
)
from parse_v2.models import ResidueResult, ParSeSummary, ParSeResult
from parse_v2.segmentation import segment_domains


# Precompute upper-triangular distance matrix for SCD: sqrt(l - k) for 0 <= k < l < 25
_SCD_WEIGHTS = np.zeros((WINDOW_SIZE, WINDOW_SIZE), dtype=np.float64)
for _k in range(WINDOW_SIZE):
    for _l in range(_k + 1, WINDOW_SIZE):
        _SCD_WEIGHTS[_k, _l] = math.sqrt(_l - _k)

_SUBSTITUTIONS = {
    "U": "C",  # Selenocysteine -> Cysteine
    "B": "D",  # Asx -> Asp
    "Z": "E",  # Glx -> Glu
    "J": "L",  # Xle -> Leu
}


def _pair_potential(n1: int, n2: int) -> float:
    """Calculate aromatic/cation-pi pair interaction term."""
    if n1 != n2:
        return float(n1 * n2) / abs(float(n1 - n2))
    return float(n1 * n2)


def clean_sequence(
    sequence: str,
    on_invalid: Literal["error", "skip", "substitute"] = "error",
) -> Optional[str]:
    """Clean and validate an amino acid sequence.

    Args:
        sequence: Raw sequence string.
        on_invalid: Policy for handling non-standard characters ('error', 'skip', 'substitute').

    Returns:
        Cleaned uppercase sequence, or None if skipped.
    """
    cleaned = "".join(sequence.split()).upper()

    if on_invalid == "substitute":
        cleaned = "".join(_SUBSTITUTIONS.get(ch, ch) for ch in cleaned)

    invalid_chars = [ch for ch in cleaned if ch not in STANDARD_AMINO_ACIDS]
    if invalid_chars:
        if on_invalid == "error":
            raise ValueError(
                f"Sequence contains invalid non-standard characters: {set(invalid_chars)}"
            )
        elif on_invalid == "skip":
            return None

    return cleaned


def parse_sequence(
    sequence: str,
    id: str = "seq",
    on_invalid: Literal["error", "skip", "substitute"] = "error",
    min_domain_len: int = DEFAULT_MIN_REGION_LEN,
    domain_purity: float = DEFAULT_PURITY_CUTOFF,
) -> Optional[ParSeResult]:
    """Execute ParSe v2 prediction for a single protein sequence.

    Args:
        sequence: Primary amino acid sequence string.
        id: Identifier for the sequence.
        on_invalid: Policy for non-standard amino acid characters.
        min_domain_len: Minimum consecutive length for domain segmentation (default: 20).
        domain_purity: Label purity threshold for domain segmentation (default: 0.90).

    Returns:
        ParSeResult containing summary metrics, per-residue profiles, and domains.
        Returns None if on_invalid='skip' and invalid characters were present.
    """
    clean_seq = clean_sequence(sequence, on_invalid=on_invalid)
    if clean_seq is None:
        return None

    seq_len = len(clean_seq)
    if seq_len < WINDOW_SIZE:
        raise ValueError(
            f"Sequence length ({seq_len}) is shorter than the minimum window size ({WINDOW_SIZE})."
        )

    # Convert sequence to numerical property arrays
    hydr_profile = np.array([HYDROPHOBICITY[aa] for aa in clean_seq], dtype=np.float64)
    ppii_profile = np.array([PPII_PROPENSITY[aa] for aa in clean_seq], dtype=np.float64)
    helix_profile = np.array([HELIX_PROPENSITY[aa] for aa in clean_seq], dtype=np.float64)
    charge_profile = np.array([CHARGE[aa] for aa in clean_seq], dtype=np.float64)

    # Output arrays of length seq_len
    labels: List[Literal["F", "D", "P"]] = ["D"] * seq_len
    dist_norm = np.zeros(seq_len, dtype=np.float64)
    labels_dh: List[Literal["F", "D", "P"]] = ["D"] * seq_len
    dist_dh = np.zeros(seq_len, dtype=np.float64)
    labels_csat: List[Literal["F", "D", "P"]] = ["D"] * seq_len
    dist_csat = np.zeros(seq_len, dtype=np.float64)

    sum_p_dist = 0.0
    sum_p_dist_dh = 0.0
    sum_p_dist_csat = 0.0

    num_windows = seq_len - WINDOW_SIZE + 1

    # Main sliding window scan
    for j in range(num_windows):
        c = j + HALF_WINDOW  # Center residue index (0-based)
        win_seq = clean_seq[j : j + WINDOW_SIZE]

        # 1. Hydrophobicity
        hydr = hydr_profile[j : j + WINDOW_SIZE].mean()

        if hydr >= FOLDED_HYDROPHOBICITY_CUTOFF:
            # Folded window
            d = (hydr - FOLDED_HYDROPHOBICITY_CUTOFF) / FOLDED_NORMALIZATION_FACTOR
            labels[c] = "F"
            dist_norm[c] = d
            labels_dh[c] = "F"
            dist_dh[c] = d
            labels_csat[c] = "F"
            dist_csat[c] = d
            continue

        # 2. Aromatic & Cation-pi interaction potential (U_pi)
        n_R = win_seq.count("R")
        n_Y = win_seq.count("Y")
        n_F = win_seq.count("F")
        n_K = win_seq.count("K")

        RY = _pair_potential(n_R, n_Y)
        RF = _pair_potential(n_R, n_F)
        KY = _pair_potential(n_K, n_Y)
        KF = _pair_potential(n_K, n_F)
        FY = _pair_potential(n_F, n_Y)

        u_pi_raw = 3.0 * RY + 2.0 * KY + 2.0 * RF + 1.0 * KF + 1.0 * FY
        u_pi = 0.137 * u_pi_raw
        u_pi_csat = 0.28 * u_pi_raw

        # 3. Electrostatics (U_q): SCD and NCPR
        win_charge = charge_profile[j : j + WINDOW_SIZE]
        scd_w = float(win_charge @ _SCD_WEIGHTS @ win_charge) / float(WINDOW_SIZE)

        n_D = win_seq.count("D")
        n_E = win_seq.count("E")
        net_charge_w = abs(n_D + n_E - n_K - n_R)
        ncpr_w = float(net_charge_w) / float(WINDOW_SIZE)

        u_q = 8.4 * scd_w + 5.6 * ncpr_w
        u_q_csat = -16.0 * scd_w + 33.0 * ncpr_w

        # 4. PPII fraction & hydrodynamic scaling exponent (nu_model)
        fppii = ppii_profile[j : j + WINDOW_SIZE].mean()
        if fppii == 1.0:
            fppii = 0.98

        v_exp = 0.503 - 0.11 * math.log(1.0 - fppii)
        rh_w = 2.16 * (100.0 ** v_exp) + 0.26 * (4.0 * net_charge_w) - 0.29 * (100.0 ** 0.5)
        nu_model = math.log(rh_w / 2.16) / math.log(100.0)

        # 5. Helix propensity
        helix = helix_profile[j : j + WINDOW_SIZE].mean()

        # 6. Distance to P/D boundary
        b = nu_model - PERP_SLOPE * helix
        x_proj = (b - BOUNDARY_INTERCEPT) / (BOUNDARY_SLOPE - PERP_SLOPE)
        y_proj = PERP_SLOPE * x_proj + b
        dist_to_boundary = math.sqrt((helix - x_proj) ** 2 + (nu_model - y_proj) ** 2)

        # Check whether in D sector or P sector
        if ((nu_model - BOUNDARY_INTERCEPT) / BOUNDARY_SLOPE) <= helix:
            # D sector
            d = dist_to_boundary / ID_REF_DIST
            labels[c] = "D"
            dist_norm[c] = d
            labels_dh[c] = "D"
            dist_dh[c] = d
            labels_csat[c] = "D"
            dist_csat[c] = d

            # Extension for Delta h°
            u_ext_dh = u_pi + u_q
            if d < u_ext_dh:
                labels_dh[c] = "P"
                d_dh = u_ext_dh - d
                dist_dh[c] = d_dh
                sum_p_dist_dh += d_dh

            # Extension for csat
            u_ext_csat = u_pi_csat + u_q_csat
            if d < u_ext_csat:
                labels_csat[c] = "P"
                d_csat = u_ext_csat - d
                dist_csat[c] = d_csat
                sum_p_dist_csat += d_csat

        else:
            # P sector
            d = dist_to_boundary / PS_REF_DIST
            labels[c] = "P"
            dist_norm[c] = d

            labels_dh[c] = "P"
            d_dh = d + (u_pi + u_q)
            dist_dh[c] = d_dh

            labels_csat[c] = "P"
            d_csat = d + (u_pi_csat + u_q_csat)
            dist_csat[c] = d_csat

            sum_p_dist += d
            sum_p_dist_dh += d_dh
            sum_p_dist_csat += d_csat

    # Edge padding: first 12 residues inherit from residue index 12 (window 0)
    center_first = HALF_WINDOW
    for idx in range(center_first):
        labels[idx] = labels[center_first]
        dist_norm[idx] = dist_norm[center_first]
        labels_dh[idx] = labels_dh[center_first]
        dist_dh[idx] = dist_dh[center_first]
        labels_csat[idx] = labels_csat[center_first]
        dist_csat[idx] = dist_csat[center_first]

    # Edge padding: last 12 residues inherit from residue index seq_len - HALF_WINDOW - 1
    center_last = seq_len - HALF_WINDOW - 1
    for idx in range(center_last + 1, seq_len):
        labels[idx] = labels[center_last]
        dist_norm[idx] = dist_norm[center_last]
        labels_dh[idx] = labels_dh[center_last]
        dist_dh[idx] = dist_dh[center_last]
        labels_csat[idx] = labels_csat[center_last]
        dist_csat[idx] = dist_csat[center_last]

    # Build ResidueResult items
    residues: List[ResidueResult] = []
    for idx in range(seq_len):
        residues.append(
            ResidueResult(
                index=idx + 1,  # 1-based index
                amino_acid=clean_seq[idx],
                label=labels[idx],
                dist=float(dist_norm[idx]),
                label_dh=labels_dh[idx],
                dist_dh=float(dist_dh[idx]),
                label_csat=labels_csat[idx],
                dist_csat=float(dist_csat[idx]),
            )
        )

    # Segment domains from base labels
    domains = segment_domains(labels, min_len=min_domain_len, purity_cutoff=domain_purity)

    num_ps = sum(1 for d in domains if d.domain_type == "P")
    num_id = sum(1 for d in domains if d.domain_type == "D")
    num_f = sum(1 for d in domains if d.domain_type == "F")

    summary = ParSeSummary(
        id=id,
        length=seq_len,
        sum_p_dist=float(sum_p_dist),
        sum_p_dist_dh=float(sum_p_dist_dh),
        sum_p_dist_csat=float(sum_p_dist_csat),
        num_ps_regions=num_ps,
        num_id_regions=num_id,
        num_folded_regions=num_f,
    )

    return ParSeResult(
        id=id,
        sequence=clean_seq,
        summary=summary,
        residues=residues,
        domains=domains,
    )
