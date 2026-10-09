"""Input/Output utilities for FASTA files and tabular export."""

import gzip
import os
from pathlib import Path
from typing import Generator, Iterable, List, Optional, Tuple, Union
import pandas as pd

from parse_v2.models import ParSeResult


def stream_fasta(
    source: Union[str, Path]
) -> Generator[Tuple[str, str], None, None]:
    """Stream (header, sequence) records from a FASTA file (plain text or gzipped).

    Args:
        source: Path to FASTA file (.fasta, .fa, .faa, .fasta.gz, etc.)

    Yields:
        (header, sequence) tuples.
    """
    path = Path(source)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    opener = gzip.open if path.suffix == ".gz" or str(path).endswith(".gz") else open
    mode = "rt"

    current_header: Optional[str] = None
    seq_parts: List[str] = []

    with opener(path, mode=mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(seq_parts)
                    seq_parts.clear()
                current_header = line[1:].strip()
            else:
                seq_parts.append(line)

        if current_header is not None:
            yield current_header, "".join(seq_parts)


def export_summary(
    results: Iterable[ParSeResult],
    output_path: Union[str, Path],
    format: Optional[str] = None,
) -> None:
    """Export summary metrics for multiple proteins.

    Args:
        results: Collection of ParSeResult objects.
        output_path: Destination file path (.tsv, .csv, .parquet).
        format: Optional explicit format override ('tsv', 'csv', 'parquet').
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = format.lower() if format else path.suffix.lstrip(".").lower()

    rows = [r.summary.to_dict() for r in results]
    df = pd.DataFrame(rows)

    if fmt in ("tsv", "tab"):
        df.to_csv(path, sep="\t", index=False)
    elif fmt == "csv":
        df.to_csv(path, sep=",", index=False)
    elif fmt in ("parquet", "pq"):
        df.to_parquet(path, index=False)
    else:
        # Default to TSV
        df.to_csv(path, sep="\t", index=False)


def export_residues(
    results: Iterable[ParSeResult],
    output_path: Union[str, Path],
    format: Optional[str] = None,
) -> None:
    """Export per-residue prediction profiles for multiple proteins.

    Args:
        results: Collection of ParSeResult objects.
        output_path: Destination file path (.tsv, .csv, .parquet).
        format: Optional explicit format override ('tsv', 'csv', 'parquet').
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = format.lower() if format else path.suffix.lstrip(".").lower()

    frames = [r.to_dataframe() for r in results]
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame()

    if fmt in ("parquet", "pq"):
        df.to_parquet(path, index=False)
    elif fmt in ("tsv", "tab"):
        df.to_csv(path, sep="\t", index=False)
    elif fmt == "csv":
        df.to_csv(path, sep=",", index=False)
    else:
        df.to_parquet(path, index=False)


def export_domains(
    results: Iterable[ParSeResult],
    output_path: Union[str, Path],
    format: Optional[str] = None,
) -> None:
    """Export parsed domain regions for multiple proteins.

    Args:
        results: Collection of ParSeResult objects.
        output_path: Destination file path (.tsv, .csv).
        format: Optional explicit format override.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = format.lower() if format else path.suffix.lstrip(".").lower()

    frames = [r.domains_to_dataframe() for r in results]
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame()

    sep = "\t" if fmt in ("tsv", "tab") else ","
    df.to_csv(path, sep=sep, index=False)


def export_filtered_fasta(
    results: Iterable[ParSeResult],
    output_path: Union[str, Path],
    min_p_length: int = 50,
) -> int:
    """Export phase-separating (P) domains of at least min_p_length to a FASTA file.

    Reproduces the filtering behavior of Steve Whitten's Parse_v2_FASTA web application.

    Returns:
        Number of filtered P domain segments written.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(path, "w", encoding="utf-8") as f:
        for r in results:
            for d in r.domains:
                if d.domain_type == "P" and d.length >= min_p_length:
                    # 1-based start and end indices: 0-based slice is [start - 1 : end]
                    subseq = r.sequence[d.first_residue - 1 : d.last_residue]
                    header = f"{r.id}_P_{d.first_residue}-{d.last_residue} residues={d.first_residue}-{d.last_residue} length={d.length}"
                    f.write(f">{header}\n{subseq}\n")
                    count += 1

    return count
