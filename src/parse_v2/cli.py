"""Command-line interface for ParSe v2."""

import sys
import time
from pathlib import Path
from typing import List, Optional
import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from parse_v2.core import parse_sequence
from parse_v2.parallel import parse_fasta
from parse_v2.io import (
    export_summary,
    export_residues,
    export_domains,
    export_filtered_fasta,
)
from parse_v2.constants import DEFAULT_MIN_REGION_LEN, DEFAULT_PURITY_CUTOFF

console = Console()
error_console = Console(stderr=True)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("input_file", type=click.Path(exists=True, dir_okay=False), required=False)
@click.option(
    "--sequence",
    "-seq",
    type=str,
    default=None,
    help="Single protein sequence to predict directly from CLI.",
)
@click.option(
    "--id",
    "seq_id",
    type=str,
    default="seq",
    help="Identifier when predicting a single --sequence.",
)
@click.option(
    "--output-summary",
    "-s",
    type=click.Path(dir_okay=False),
    default=None,
    help="Output file for protein-level summary metrics (.tsv, .csv, .parquet).",
)
@click.option(
    "--output-residues",
    "-r",
    type=click.Path(dir_okay=False),
    default=None,
    help="Output file for per-residue prediction profiles (.parquet, .tsv, .csv).",
)
@click.option(
    "--output-domains",
    "-d",
    type=click.Path(dir_okay=False),
    default=None,
    help="Output file for parsed domain architecture (.tsv, .csv).",
)
@click.option(
    "--output-filtered-fasta",
    "-f",
    type=click.Path(dir_okay=False),
    default=None,
    help="Export filtered phase-separating (P) domains to FASTA (web tool parity).",
)
@click.option(
    "--filter-min-len",
    type=int,
    default=50,
    show_default=True,
    help="Minimum residue length for --output-filtered-fasta.",
)
@click.option(
    "--on-invalid",
    type=click.Choice(["error", "skip", "substitute"], case_sensitive=False),
    default="error",
    show_default=True,
    help="Handling policy for non-standard amino acid characters.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=None,
    help="Number of parallel processes for FASTA processing (defaults to CPU count).",
)
@click.option(
    "--min-domain-len",
    type=int,
    default=DEFAULT_MIN_REGION_LEN,
    show_default=True,
    help="Minimum length for domain segmentation.",
)
@click.option(
    "--domain-purity",
    type=float,
    default=DEFAULT_PURITY_CUTOFF,
    show_default=True,
    help="Label purity cutoff for domain segmentation.",
)
@click.option(
    "--quiet",
    "-q",
    is_flag=True,
    default=False,
    help="Suppress progress bars and status messages.",
)
def main(
    input_file: Optional[str],
    sequence: Optional[str],
    seq_id: str,
    output_summary: Optional[str],
    output_residues: Optional[str],
    output_domains: Optional[str],
    output_filtered_fasta: Optional[str],
    filter_min_len: int,
    on_invalid: str,
    threads: Optional[int],
    min_domain_len: int,
    domain_purity: float,
    quiet: bool,
) -> None:
    """ParSe v2: Predict protein liquid-liquid phase separation potential.

    Provide either an INPUT_FILE (FASTA format, plain text or .gz) or a --sequence.
    """
    if not input_file and not sequence:
        error_console.print("[bold red]Error:[/bold red] Please provide either an INPUT_FILE or --sequence.")
        sys.exit(1)

    start_time = time.time()

    # --- Mode 1: Single Sequence ---
    if sequence:
        try:
            result = parse_sequence(
                sequence=sequence,
                id=seq_id,
                on_invalid=on_invalid.lower(),  # type: ignore
                min_domain_len=min_domain_len,
                domain_purity=domain_purity,
            )
        except Exception as e:
            error_console.print(f"[bold red]Prediction Error:[/bold red] {e}")
            sys.exit(1)

        if result is None:
            error_console.print("[yellow]Sequence skipped due to non-standard residues.[/yellow]")
            sys.exit(0)

        # Print summary to console
        if not quiet:
            table = Table(title=f"ParSe v2 Results: {result.id}", show_header=True, header_style="bold cyan")
            table.add_column("Metric", style="bold")
            table.add_column("Value", justify="right")

            table.add_row("Sequence Length", str(result.summary.length))
            table.add_row("∑ Distance (P windows)", f"{result.summary.sum_p_dist:.3f}")
            table.add_row("∑ Distance (+Uπ +Uq, ∆h° trained)", f"{result.summary.sum_p_dist_dh:.3f}")
            table.add_row("∑ Distance (+Uπ +Uq, csat trained)", f"{result.summary.sum_p_dist_csat:.3f}")
            table.add_row("PS Regions (Blue, P)", str(result.summary.num_ps_regions))
            table.add_row("ID Regions (Red, D)", str(result.summary.num_id_regions))
            table.add_row("Folded Regions (Black, F)", str(result.summary.num_folded_regions))
            console.print(table)

            if result.domains:
                dom_table = Table(title="Parsed Domains", show_header=True, header_style="bold magenta")
                dom_table.add_column("Domain #")
                dom_table.add_column("Type")
                dom_table.add_column("Start Residue", justify="right")
                dom_table.add_column("End Residue", justify="right")
                dom_table.add_column("Length", justify="right")

                for d in result.domains:
                    dom_table.add_row(
                        str(d.region_id),
                        d.domain_type,
                        str(d.first_residue),
                        str(d.last_residue),
                        str(d.length),
                    )
                console.print(dom_table)

        results = [result]

    # --- Mode 2: FASTA File ---
    else:
        results = []
        if not quiet:
            console.print(f"[bold green]Processing FASTA:[/bold green] {input_file}")

        try:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                console=console,
                disable=quiet,
            ) as progress:
                task = progress.add_task("[cyan]Calculating ParSe v2 predictions...", total=None)
                count = 0
                for r in parse_fasta(
                    source=input_file,  # type: ignore
                    n_jobs=threads,
                    on_invalid=on_invalid.lower(),  # type: ignore
                    min_domain_len=min_domain_len,
                    domain_purity=domain_purity,
                ):
                    results.append(r)
                    count += 1
                    progress.update(task, description=f"[cyan]Parsed {count} proteins...")
        except Exception as e:
            error_console.print(f"[bold red]Batch Processing Error:[/bold red] {e}")
            sys.exit(1)

        if not quiet:
            elapsed = time.time() - start_time
            console.print(f"[bold green]✓ Completed[/bold green] {len(results)} proteins in {elapsed:.2f}s")

    # --- Export Files ---
    if output_summary:
        export_summary(results, output_summary)
        if not quiet:
            console.print(f"  [dim]Saved summary to:[/dim] {output_summary}")

    if output_residues:
        export_residues(results, output_residues)
        if not quiet:
            console.print(f"  [dim]Saved residue profiles to:[/dim] {output_residues}")

    if output_domains:
        export_domains(results, output_domains)
        if not quiet:
            console.print(f"  [dim]Saved domain regions to:[/dim] {output_domains}")

    if output_filtered_fasta:
        n_filtered = export_filtered_fasta(results, output_filtered_fasta, min_p_length=filter_min_len)
        if not quiet:
            console.print(f"  [dim]Saved {n_filtered} filtered P-regions (≥{filter_min_len} aa) to:[/dim] {output_filtered_fasta}")


if __name__ == "__main__":
    main()
