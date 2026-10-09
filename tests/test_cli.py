"""Tests for the parse-v2 command line interface."""

import tempfile
from pathlib import Path
from click.testing import CliRunner
import pandas as pd

from parse_v2.cli import main


def test_cli_sequence_mode():
    runner = CliRunner()
    seq = "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG"
    result = runner.invoke(main, ["--sequence", seq, "--id", "FUS_LCD"])
    assert result.exit_code == 0
    assert "ParSe v2 Results: FUS_LCD" in result.output
    assert "835.946" in result.output


def test_cli_fasta_mode_with_all_outputs():
    runner = CliRunner()
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "test.fasta"
        fasta_path.write_text(
            ">P1\nMASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG\n"
            ">P2\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"
        )

        sum_out = Path(tmpdir) / "summary.tsv"
        res_out = Path(tmpdir) / "residues.parquet"
        dom_out = Path(tmpdir) / "domains.tsv"
        filt_out = Path(tmpdir) / "filtered.fasta"

        args = [
            str(fasta_path),
            "--output-summary", str(sum_out),
            "--output-residues", str(res_out),
            "--output-domains", str(dom_out),
            "--output-filtered-fasta", str(filt_out),
            "--threads", "2",
        ]

        result = runner.invoke(main, args)
        assert result.exit_code == 0
        assert sum_out.exists()
        assert res_out.exists()
        assert dom_out.exists()
        assert filt_out.exists()

        # Check summary
        df_sum = pd.read_csv(sum_out, sep="\t")
        assert len(df_sum) == 2

        # Check residues parquet
        df_res = pd.read_parquet(res_out)
        assert len(df_res) == 214 + 76

        # Check filtered fasta
        assert ">P1_P_1-214" in filt_out.read_text()
