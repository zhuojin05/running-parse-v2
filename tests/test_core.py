"""Tests for edge cases and core calculation logic."""

import pytest
from parse_v2.core import parse_sequence, clean_sequence


def test_short_sequence():
    """Sequences under window size (25 residues) must raise ValueError."""
    with pytest.raises(ValueError, match="shorter than the minimum window size"):
        parse_sequence("ACDEFGHIKLMNPQRSTVWY", id="short")


def test_polyproline_clamp():
    """Polyproline (fppii = 1.0) must not cause math domain error log(1 - 1)."""
    seq = "P" * 30
    res = parse_sequence(seq, id="polyP")
    assert res is not None
    assert res.summary.length == 30
    # All residues must be classified cleanly
    assert all(r.label in ("P", "D", "F") for r in res.residues)


def test_invalid_characters_handling():
    """Test 'error', 'skip', and 'substitute' policies for non-standard amino acids."""
    seq_with_x = "MASNDYTQQATQSYGAYPTQPGQGX"  # Length 25 with X

    # Error policy
    with pytest.raises(ValueError, match="invalid non-standard characters"):
        parse_sequence(seq_with_x, on_invalid="error")

    # Skip policy
    assert parse_sequence(seq_with_x, on_invalid="skip") is None

    # Substitute policy with Selenocysteine (U -> C)
    seq_with_u = "MASNDYTQQATQSYGAYPTQPGQGU"
    res = parse_sequence(seq_with_u, on_invalid="substitute")
    assert res is not None
    assert res.sequence[-1] == "C"


def test_sequence_cleaning():
    """Whitespace and lower-case must be cleaned automatically."""
    seq = "  masn dytq\nq atqs ygay ptqp gqgy \n" + "A" * 15
    cleaned = clean_sequence(seq)
    assert cleaned is not None
    assert " " not in cleaned
    assert "\n" not in cleaned
    assert cleaned.isupper()
