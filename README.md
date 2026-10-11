# ParSe v2 (Python / `uv`)

A high-performance, batch-capable Python package and command-line runner for **ParSe v2** (*Partitioning Sequence to predict protein liquid–liquid phase separation*), managed with [`uv`](https://github.com/astral-sh/uv).

Developed based on the algorithms of Steve Whitten's group (Texas State University):
> **ParSe 2.0: A web tool to identify drivers of protein phase separation at the proteome level**  
> Whitten, S. T. et al., *Protein Science* (2023) **32**(10): e4756.  
> DOI: [10.1002/pro.4756](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4756)

---

## Highlights

- **Pure NumPy Vectorized Engine**: Zero external Fortran compiler required. Runs seamlessly on macOS (Apple Silicon / Intel), Linux HPC, and Windows.
- **100% Mathematical Parity**: Verified against Steve Whitten's original Fortran 77 source (`Parse_v2.f`) and web tool (`parse.js`).
- **High-Throughput Batch FASTA Processing**: Multi-process parallel streaming across multi-core CPUs.
- **Modern Data Science Integration**: Export per-residue prediction profiles directly to **Parquet** (recommended for multi-million residue proteomes), TSV, or CSV, with native Pandas DataFrame conversion.
- **No Browser Limits**: Avoid browser memory limits, tab crashes, or Fortran file-collision bugs during automated screens.

---

## Installation & Quickstart

This repository uses [`uv`](https://docs.astral.sh/uv/) for instant environment management.

### Clone and Sync
```bash
git clone https://github.com/zhuojin05/running-parse-v2.git
cd running-parse-v2

# Install dependencies and sync virtual environment
uv sync
```

---

## Command-Line Usage (`parse-v2`)

### 1. Single Sequence Quick Assessment
```bash
uv run parse-v2 --sequence "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG" --id "FUS_LCD"
```

Output:
```text
           ParSe v2 Results: FUS_LCD            
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Metric                             ┃   Value ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ Sequence Length                    │     214 │
│ ∑ Distance (P windows)             │ 835.946 │
│ ∑ Distance (+Uπ +Uq, ∆h° trained)  │ 864.417 │
│ ∑ Distance (+Uπ +Uq, csat trained) │ 915.712 │
│ PS Regions (Blue, P)               │       1 │
│ ID Regions (Red, D)                │       0 │
│ Folded Regions (Black, F)          │       0 │
└────────────────────────────────────┴─────────┘
```

### 2. Batch FASTA Processing
Process a FASTA file (plain text or `.gz`) and generate summary metrics, per-residue profiles, domain architectures, and filtered phase-separating segments:

```bash
uv run parse-v2 proteome.fasta \
    --output-summary summary.tsv \
    --output-residues residues.parquet \
    --output-domains domains.tsv \
    --output-filtered-fasta ps_regions_ge50.fasta \
    --threads 8
```

### CLI Options

| Option | Flag | Description | Default |
| :--- | :--- | :--- | :--- |
| `--sequence` | `-seq` | Predict on a single raw sequence string | `None` |
| `--id` | | Identifier for single sequence | `seq` |
| `--output-summary` | `-s` | Output path for protein-level summary (`.tsv`, `.csv`, `.parquet`) | `None` |
| `--output-residues`| `-r` | Output path for residue-level predictions (`.parquet`, `.tsv`, `.csv`) | `None` |
| `--output-domains` | `-d` | Output path for segmented domains (`.tsv`, `.csv`) | `None` |
| `--output-filtered-fasta` | `-f` | Export filtered $P$ regions to FASTA (web tool parity) | `None` |
| `--filter-min-len` | | Minimum length of $P$ regions for filtered FASTA export | `50` |
| `--threads` | `-t` | Number of worker processes (defaults to CPU count) | All CPUs |
| `--on-invalid` | | Non-standard AA policy: `error`, `skip`, or `substitute` | `error` |
| `--min-domain-len` | | Minimum domain length for segmentation | `20` |
| `--domain-purity` | | Minimum label purity for domain segmentation | `0.90` |
| `--quiet` | `-q` | Suppress progress bars and table prints | `False` |

---

## Python API

You can easily integrate `parse_v2` into Python workflows, Jupyter notebooks, or Nextflow/Snakemake scripts:

### Single Sequence
```python
from parse_v2 import parse_sequence

fus_lcd = "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQ..."
result = parse_sequence(fus_lcd, id="FUS_LCD")

# Summary metrics
print(result.summary)
# ParSeSummary(id='FUS_LCD', length=214, sum_p_dist=835.946, sum_p_dist_dh=864.417, sum_p_dist_csat=915.712, ...)

# Convert per-residue predictions to Pandas DataFrame
df_residues = result.to_dataframe()
# Columns: ['protein_id', 'residue_num', 'amino_acid', 'label', 'dist', 'label_dh', 'dist_dh', 'label_csat', 'dist_csat']

# Inspect parsed domains
for domain in result.domains:
    print(domain)
    # DomainRegion(region_id=1, domain_type='P', first_residue=1, last_residue=214, length=214)
```

### Multi-Sequence Batch (FASTA Streaming)
```python
from parse_v2 import parse_fasta, export_summary, export_residues

# Process FASTA with multi-processing (returns a streaming generator)
results = list(parse_fasta("uniprot_proteome.fasta.gz", n_jobs=8))

# Export results
export_summary(results, "proteome_summary.tsv")
export_residues(results, "proteome_residues.parquet")
```

---

## Verification & Parity Testing

The package includes automated tests comparing Python outputs against the original Fortran 77 implementation (`Parse_v2.f`) compiled via `gfortran`:

```bash
uv run pytest -v
```

Tests verify:
1. **Residue-by-residue label parity** (`label`, `label_dh`, `label_csat`) on benchmark proteins (FUS LCD, Ubiquitin, hnRNPA1).
2. **Classifier distances** agreement within single-to-double floating-point precision ($< 0.005$ tolerance).
3. **Polyproline singularity clamp** ($f_\text{PPII} = 1.0 \to 0.98$).
4. **FASTA streaming & decompression** (plain text and `.gz`).
5. **Multi-process CLI workflows**.

---

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0), matching the original [ParSe_v2](https://github.com/stevewhitten/ParSe_v2) repository license.
