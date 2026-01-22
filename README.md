# mkado

A modern Python implementation of the McDonald-Kreitman test toolkit for detecting selection in molecular evolution.

## Features

- **Standard MK test**: Classic 2x2 contingency table with Fisher's exact test
- **Polarized MK test**: Uses a third outgroup to assign mutations to lineages
- **Asymptotic MK test**: Frequency-bin α estimates with exponential extrapolation (Messer & Petrov 2013)
- **Batch processing**: Process multiple genes with parallel execution
- **Multiple output formats**: Pretty-print, TSV, and JSON

## Installation

```bash
# Clone the repository
git clone https://github.com/andrewkern/mkado.git
cd mkado

# Install with uv
uv sync

# Or install with pip
pip install .
```

## Quick Start

```bash
# Standard MK test (combined alignment file)
mkado test alignment.fa -i "dmel" -o "dsim"

# Asymptotic MK test
mkado test alignment.fa -i "dmel" -o "dsim" -a

# Polarized MK test
mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"

# Batch process a directory
mkado batch alignments/ -i "dmel" -o "dsim"

# Batch with asymptotic test and 8 parallel workers
mkado batch alignments/ -i "dmel" -o "dsim" -a -w 8

# Get file info
mkado info sequences.fa
```

## Usage Modes

mkado supports two modes for specifying ingroup/outgroup sequences:

### Combined File Mode (Recommended)

Use `-i` and `-o` to filter sequences by name pattern from a single alignment file:

```bash
mkado test alignment.fa -i "speciesA" -o "speciesB"
mkado batch alignments/ -i "speciesA" -o "speciesB"
```

### Separate Files Mode

Provide separate FASTA files for ingroup and outgroup:

```bash
mkado test ingroup.fa outgroup.fa
mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
```

## Commands

### `mkado test`

Run MK test on a single alignment.

```bash
mkado test FASTA [OUTGROUP_FILE] [OPTIONS]
```

**Key Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--ingroup-match` | `-i` | Ingroup sequence name pattern (combined mode) |
| `--outgroup-match` | `-o` | Outgroup sequence name pattern (combined mode) |
| `--asymptotic` | `-a` | Use asymptotic MK test |
| `--polarize` | `-p` | Second outgroup file (separate files mode) |
| `--polarize-match` | | Second outgroup pattern (combined mode) |
| `--bins` | `-b` | Frequency bins for asymptotic test (default: 10) |
| `--format` | `-f` | Output format: pretty, tsv, json |
| `--reading-frame` | `-r` | Reading frame 1-3 (default: 1) |

**Examples:**

```bash
# Combined file mode
mkado test alignment.fa -i "dmel" -o "dsim"
mkado test alignment.fa -i "dmel" -o "dsim" -a -b 20
mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"

# Separate files mode
mkado test ingroup.fa outgroup.fa
mkado test ingroup.fa outgroup.fa -a
mkado test ingroup.fa outgroup.fa -p outgroup2.fa
```

### `mkado batch`

Run MK test on multiple alignment files.

```bash
mkado batch DIRECTORY [OPTIONS]
```

**Key Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--ingroup-match` | `-i` | Ingroup pattern (enables combined file mode) |
| `--outgroup-match` | `-o` | Outgroup pattern (required with -i) |
| `--asymptotic` | `-a` | Use asymptotic MK test |
| `--aggregate/--per-gene` | | Aggregate results or per-gene (asymptotic) |
| `--pattern` | | File glob pattern (default: auto-detect *.fa, *.fasta, *.fna) |
| `--workers` | `-w` | Parallel workers (0=auto, 1=sequential) |
| `--bins` | `-b` | Frequency bins for asymptotic test |
| `--format` | `-f` | Output format: pretty, tsv, json |

**Examples:**

```bash
# Combined file mode (recommended)
mkado batch alignments/ -i "dmel" -o "dsim"
mkado batch alignments/ -i "dmel" -o "dsim" -a
mkado batch alignments/ -i "dmel" -o "dsim" -a --per-gene
mkado batch alignments/ -i "dmel" -o "dsim" -w 8

# Separate files mode
mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
```

### `mkado info`

Display information about a FASTA file.

```bash
mkado info FASTA [-r READING_FRAME]
```

## Example Output

```
$ mkado test alignment.fa -i "kreitman" -o "mauritiana"

Found 11 ingroup, 1 outgroup sequences
MK Test Results:
  Divergence:    Dn=6, Ds=8
  Polymorphism:  Pn=1, Ps=8
  Fisher's exact p-value: 0.176
  Neutrality Index (NI):  0.1667
  Alpha (α):              0.8333
```

## Python API

```python
from mkado import mk_test, asymptotic_mk_test, SequenceSet

# Run MK test
result = mk_test("ingroup.fa", "outgroup.fa")
print(f"Alpha: {result.alpha}")
print(f"P-value: {result.p_value}")

# Run asymptotic MK test
result = asymptotic_mk_test("ingroup.fa", "outgroup.fa")
print(f"Asymptotic Alpha: {result.alpha_asymptotic}")
print(f"95% CI: {result.ci_low} - {result.ci_high}")

# Combined file mode - filter by sequence name
all_seqs = SequenceSet.from_fasta("combined.fa")
ingroup = all_seqs.filter_by_name("dmel")
outgroup = all_seqs.filter_by_name("dsim")
result = mk_test(ingroup, outgroup)
```

## Interpretation

### Neutrality Index (NI)

- **NI = 1**: Neutral evolution
- **NI > 1**: Excess polymorphism (negative/purifying selection)
- **NI < 1**: Excess divergence (positive selection)

### Alpha (α)

- **α = 0**: No adaptive substitutions
- **α > 0**: Proportion of substitutions driven by positive selection
- **α < 0**: Excess polymorphism relative to divergence

## Development

```bash
# Install dev dependencies
uv sync

# Run tests
uv run pytest

# Run linter
uv run ruff check src/

# Run formatter
uv run ruff format src/
```

## References

- McDonald, J. H., & Kreitman, M. (1991). Adaptive protein evolution at the Adh locus in Drosophila. Nature, 351(6328), 652-654.
- Messer, P. W., & Petrov, D. A. (2013). Frequent adaptation and the McDonald–Kreitman test. PNAS, 110(21), 8615-8620.

## License

MIT License
