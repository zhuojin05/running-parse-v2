"""Physical property scales and mathematical constants for ParSe v2.

References:
- Bastolla et al., Proteins 58, 22-30 (2005) (Contact-based interactivity/hydrophobicity)
- Hilser group, Prot Sci 22, 405-417 (2013) (PPII propensities in peptides)
- Tanaka & Scheraga, Macromolecules 10, 9-20 (1977); Kidera normalized (Alpha-helix)
- Whitten et al., Protein Science 32, e4756 (2023) (ParSe 2.0)
"""

import math

STANDARD_AMINO_ACIDS = (
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
)

# Bastolla contact-matrix-based hydrophobicity scale
HYDROPHOBICITY = {
    "A": 0.0728, "C": 0.3557, "D": -0.0552, "E": -0.0295, "F": 0.4201,
    "G": -0.0589, "H": 0.0874, "I": 0.3805, "K": -0.0053, "L": 0.3819,
    "M": 0.1613, "N": -0.0390, "P": -0.0492, "Q": 0.0126, "R": 0.0394,
    "S": -0.0282, "T": 0.0239, "V": 0.2947, "W": 0.4114, "Y": 0.3113,
}

# Polyproline II (PPII) propensity measured in peptides by Hilser lab
PPII_PROPENSITY = {
    "A": 0.37, "C": 0.25, "D": 0.30, "E": 0.42, "F": 0.17,
    "G": 0.13, "H": 0.20, "I": 0.39, "K": 0.56, "L": 0.24,
    "M": 0.36, "N": 0.27, "P": 1.00, "Q": 0.53, "R": 0.38,
    "S": 0.24, "T": 0.32, "V": 0.39, "W": 0.25, "Y": 0.25,
}

# Normalized frequency of alpha-helix (Tanaka-Scheraga / Kidera)
HELIX_PROPENSITY = {
    "A": 1.42, "C": 0.73, "D": 1.01, "E": 1.63, "F": 1.16,
    "G": 0.50, "H": 1.20, "I": 1.12, "K": 1.24, "L": 1.29,
    "M": 1.21, "N": 0.71, "P": 0.65, "Q": 1.02, "R": 1.06,
    "S": 0.71, "T": 0.78, "V": 0.99, "W": 1.05, "Y": 0.67,
}

# Charge lambda for Sequence Charge Decoration (SCD)
CHARGE = {
    "A": 0.0, "C": 0.0, "D": -1.0, "E": -1.0, "F": 0.0,
    "G": 0.0, "H": 0.0, "I": 0.0, "K": 1.0, "L": 0.0,
    "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 1.0,
    "S": 0.0, "T": 0.0, "V": 0.0, "W": 0.0, "Y": 0.0,
}

# Window parameters
WINDOW_SIZE = 25
HALF_WINDOW = WINDOW_SIZE // 2  # 12

# Folded sector boundary
FOLDED_HYDROPHOBICITY_CUTOFF = 0.08280152
FOLDED_NORMALIZATION_FACTOR = 0.01679414 * 2.0  # 0.03358828

# Boundary line between P and D sectors: y = SLOPE * x + INTERCEPT
# where x = helix propensity, y = nu_model
BOUNDARY_SLOPE = -0.244078945
BOUNDARY_INTERCEPT = 0.7885823
PERP_SLOPE = -1.0 / BOUNDARY_SLOPE

# Reference centroids
ID_MEAN_HELIX = 1.022552
ID_MEAN_NU = 0.5582901

PS_MEAN_HELIX = 0.9327272
PS_MEAN_NU = 0.5416

def _compute_reference_distance(mean_x: float, mean_y: float) -> float:
    b = mean_y - PERP_SLOPE * mean_x
    x_proj = (b - BOUNDARY_INTERCEPT) / (BOUNDARY_SLOPE - PERP_SLOPE)
    y_proj = PERP_SLOPE * x_proj + b
    return math.sqrt((mean_x - x_proj) ** 2 + (mean_y - y_proj) ** 2)

ID_REF_DIST = _compute_reference_distance(ID_MEAN_HELIX, ID_MEAN_NU)
PS_REF_DIST = _compute_reference_distance(PS_MEAN_HELIX, PS_MEAN_NU)

# Segmentation defaults
DEFAULT_MIN_REGION_LEN = 20
DEFAULT_PURITY_CUTOFF = 0.90
