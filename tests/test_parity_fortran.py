"""Verification test suite asserting parity between Python and Fortran 77 implementation."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import pytest
import numpy as np
import pandas as pd

from parse_v2.core import parse_sequence

REFERENCE_FORTRAN = Path(__file__).parent / "reference" / "Parse_v2.f"

# Benchmark sequences
BENCHMARKS = {
    "FUS_LCD": (
        "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYG"
        "QSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQ"
        "SSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNY"
        "GQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG"
    ),
    "UBIQUITIN": (
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYN"
        "IQKESTLHLVLRLRGG"
    ),
    "hnRNPA1": (
        "MSKSESPKEPEQLRKLFIGGLSFETTDESLRSHFEQWGTLTDCVVMRDPNTKRSRGFGFV"
        "TYATVEEVDAAMNARPHKVDGRVVEPKRAVSREDSQRPGAHLTVKKIFVGGIKEDTEEEH"
        "LRDYFEQYGKIEVIEIMTDRGSGKKRGFAFVTFDDHDSVDKIVIQKYHTVNGHNCEVRKA"
        "LSKQEMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGS"
        "GDGYNGFGNDGGYGGGGPGYSGGSRGYGSGGQGYGNQGSGYGGSGSYDSYNNGGGGGFGG"
        "GSGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF"
    ),
}


def _run_fortran(sequence: str) -> pd.DataFrame:
    """Compile and run the original Fortran program in an isolated directory.

    Returns the parsed residue_label_classifier_distance.csv as a DataFrame.
    """
    gfortran_bin = shutil.which("gfortran")
    if not gfortran_bin:
        pytest.skip("gfortran compiler not available for Fortran parity comparison.")

    with tempfile.TemporaryDirectory() as tmpdir:
        src_f = Path(tmpdir) / "Parse_v2.f"
        bin_out = Path(tmpdir) / "parse_v2"
        shutil.copy(REFERENCE_FORTRAN, src_f)

        compile_proc = subprocess.run(
            [gfortran_bin, "-O2", str(src_f), "-o", str(bin_out)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )
        assert compile_proc.returncode == 0, f"Fortran compilation failed: {compile_proc.stderr}"

        run_proc = subprocess.run(
            [str(bin_out), sequence],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )
        assert run_proc.returncode == 0, f"Fortran execution failed: {run_proc.stderr}"

        csv_path = Path(tmpdir) / "residue_label_classifier_distance.csv"
        assert csv_path.exists(), "Fortran did not produce residue_label_classifier_distance.csv"

        # CSV format: residue_num, aa, label, dist, label_dh, dist_dh, label_csat, dist_csat
        df = pd.read_csv(
            csv_path,
            header=None,
            names=["residue_num", "aa", "label", "dist", "label_dh", "dist_dh", "label_csat", "dist_csat"],
            skipinitialspace=True,
        )
        return df


@pytest.mark.parametrize("name,sequence", list(BENCHMARKS.items()))
def test_fortran_residue_parity(name: str, sequence: str):
    """Assert residue-by-residue label parity and floating-point distance agreement."""
    df_fortran = _run_fortran(sequence)
    py_result = parse_sequence(sequence, id=name)

    assert py_result is not None
    assert len(py_result.residues) == len(df_fortran)

    df_py = py_result.to_dataframe()

    # 1. Base labels must match 100% exactly
    assert (df_py["label"].values == df_fortran["label"].values).all(), f"Base labels mismatched for {name}"

    # 2. Delta h° labels must match 100% exactly
    assert (df_py["label_dh"].values == df_fortran["label_dh"].values).all(), f"Label dh mismatched for {name}"

    # 3. csat labels must match 100% exactly
    assert (df_py["label_csat"].values == df_fortran["label_csat"].values).all(), f"Label csat mismatched for {name}"

    # 4. Floating point distances match within single-precision tolerance (< 0.005 absolute)
    np.testing.assert_allclose(df_py["dist"].values, df_fortran["dist"].values, atol=5e-3)
    np.testing.assert_allclose(df_py["dist_dh"].values, df_fortran["dist_dh"].values, atol=5e-3)
    np.testing.assert_allclose(df_py["dist_csat"].values, df_fortran["dist_csat"].values, atol=5e-3)
