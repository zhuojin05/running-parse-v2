"""Tests for IO functions and format exporters."""

import gzip
import tempfile
from pathlib import Path
import pandas as pd
import pytest

from parse_v2.core import parse_sequence
from parse_v2.io import (
    stream_fasta,
    export_summary,
    export_residues,
    export_domains,
    export_filtered_fasta,
)


@pytest.fixture
def sample_results():
    fus_seq = (
        "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYG"
        "QSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQ"
        "SSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNY"
        "GQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG"
    )
    ubi_seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    r1 = parse_sequence(fus_seq, id="FUS")
    r2 = parse_sequence(ubi_seq, id="UBIQUITIN")
    return [r1, r2]


def test_stream_fasta_plain_and_gzip():
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_txt = Path(tmpdir) / "test.fasta"
        fasta_gz = Path(tmpdir) / "test.fasta.gz"

        content = ">seq1 First\nACDEFGHIKLMNPQRSTVWYACDEF\n>seq2 Second\nMQIFVKTLTGKTITLEVEPSDTIEN\n"
        fasta_txt.write_text(content)

        with gzip.open(fasta_gz, "wt") as f:
            f.write(content)

        # Test plain text
        records_txt = list(stream_fasta(fasta_txt))
        assert len(records_txt) == 2
        assert records_txt[0][0] == "seq1 First"
        assert records_txt[0][1] == "ACDEFGHIKLMNPQRSTVWYACDEF"

        # Test gzipped
        records_gz = list(stream_fasta(fasta_gz))
        assert len(records_gz) == 2
        assert records_gz[1][0] == "seq2 Second"
        assert records_gz[1][1] == "MQIFVKTLTGKTITLEVEPSDTIEN"


def test_export_summary(sample_results):
    with tempfile.TemporaryDirectory() as tmpdir:
        tsv_out = Path(tmpdir) / "summary.tsv"
        export_summary(sample_results, tsv_out)
        assert tsv_out.exists()

        df = pd.read_csv(tsv_out, sep="\t")
        assert len(df) == 2
        assert "sum_p_dist" in df.columns
        assert df.iloc[0]["id"] == "FUS"


def test_export_residues_parquet(sample_results):
    with tempfile.TemporaryDirectory() as tmpdir:
        pq_out = Path(tmpdir) / "residues.parquet"
        export_residues(sample_results, pq_out)
        assert pq_out.exists()

        df = pd.read_parquet(pq_out)
        # Total residues: 214 + 76 = 290
        assert len(df) == 290
        assert "label" in df.columns
        assert "dist_dh" in df.columns


def test_export_filtered_fasta(sample_results):
    with tempfile.TemporaryDirectory() as tmpdir:
        fa_out = Path(tmpdir) / "filtered.fasta"
        # FUS has a 214-residue P domain, Ubiquitin has none
        count = export_filtered_fasta(sample_results, fa_out, min_p_length=50)
        assert count == 1
        assert fa_out.exists()
        text = fa_out.read_text()
        assert ">FUS_P_1-214" in text
