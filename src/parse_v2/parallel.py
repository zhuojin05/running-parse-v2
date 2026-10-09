"""Parallel batch processing for multi-sequence FASTA files."""

import functools
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from pathlib import Path
from typing import Generator, Iterable, List, Literal, Optional, Tuple, Union

from parse_v2.core import parse_sequence
from parse_v2.io import stream_fasta
from parse_v2.models import ParSeResult
from parse_v2.constants import DEFAULT_MIN_REGION_LEN, DEFAULT_PURITY_CUTOFF


def _worker_wrapper(
    args: Tuple[str, str, str, int, float]
) -> Optional[ParSeResult]:
    """Helper worker for multiprocessing pool."""
    header, sequence, on_invalid, min_domain_len, domain_purity = args
    # Extract clean sequence ID (up to first whitespace)
    seq_id = header.split()[0] if header else "seq"
    try:
        return parse_sequence(
            sequence=sequence,
            id=seq_id,
            on_invalid=on_invalid,  # type: ignore
            min_domain_len=min_domain_len,
            domain_purity=domain_purity,
        )
    except Exception as exc:
        if on_invalid == "error":
            raise
        return None


def parse_fasta(
    source: Union[str, Path],
    n_jobs: Optional[int] = None,
    on_invalid: Literal["error", "skip", "substitute"] = "error",
    min_domain_len: int = DEFAULT_MIN_REGION_LEN,
    domain_purity: float = DEFAULT_PURITY_CUTOFF,
    chunksize: int = 50,
) -> Generator[ParSeResult, None, None]:
    """Process a FASTA file in parallel or sequential mode.

    Args:
        source: Path to FASTA file (plain text or gzipped).
        n_jobs: Number of parallel worker processes. Default: os.cpu_count().
        on_invalid: Policy for handling invalid characters ('error', 'skip', 'substitute').
        min_domain_len: Minimum consecutive length for domain segmentation.
        domain_purity: Label purity threshold for domain segmentation.
        chunksize: Number of sequences to dispatch per worker chunk.

    Yields:
        ParSeResult for each valid protein.
    """
    if n_jobs is None:
        n_jobs = max(1, (os.cpu_count() or 1))

    # Sequential execution if n_jobs == 1
    if n_jobs == 1:
        for header, seq in stream_fasta(source):
            seq_id = header.split()[0] if header else "seq"
            try:
                res = parse_sequence(
                    sequence=seq,
                    id=seq_id,
                    on_invalid=on_invalid,
                    min_domain_len=min_domain_len,
                    domain_purity=domain_purity,
                )
                if res is not None:
                    yield res
            except Exception:
                if on_invalid == "error":
                    raise
                continue
        return

    # Parallel execution with ProcessPoolExecutor
    def _arg_generator():
        for header, seq in stream_fasta(source):
            yield (header, seq, on_invalid, min_domain_len, domain_purity)

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        for result in executor.map(_worker_wrapper, _arg_generator(), chunksize=chunksize):
            if result is not None:
                yield result
