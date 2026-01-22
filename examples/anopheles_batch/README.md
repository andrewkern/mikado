# Anopheles Batch Example

This example demonstrates batch processing of MK tests using aligned coding sequences from three Anopheles mosquito species.

## Data

Eleven genes with aligned codon sequences (27 sequences each: 12 gamb, 12 afun, 3 amin):
- `AGAP000010.fa`
- `AGAP000011.fa`
- `AGAP000012.fa`
- `AGAP000014.fa`
- `AGAP000016.fa`
- `AGAP000021.fa`
- `AGAP000022.fa`
- `AGAP000023.fa`
- `AGAP000038.fa`
- `AGAP000041.fa`
- `AGAP029936.fa`

Species:
- **gamb**: Anopheles gambiae (ingroup)
- **afun**: Anopheles funestus (outgroup)
- **amin**: Anopheles minimus (second outgroup for polarized tests)

## Example Commands

### Standard MK Test

```bash
# Batch MK test: gamb (ingroup) vs afun (outgroup)
mkado batch examples/anopheles_batch/ -i gamb -o afun

# TSV output for downstream analysis
mkado batch examples/anopheles_batch/ -i gamb -o afun -f tsv

# JSON output
mkado batch examples/anopheles_batch/ -i gamb -o afun -f json
```

### Asymptotic MK Test

```bash
# Aggregate asymptotic MK test (pools polymorphism across genes)
mkado batch examples/anopheles_batch/ -i gamb -o afun -a

# Per-gene asymptotic analysis
mkado batch examples/anopheles_batch/ -i gamb -o afun -a --per-gene

# Custom frequency bins
mkado batch examples/anopheles_batch/ -i gamb -o afun -a -b 5
```

### Polarized MK Test

```bash
# Use amin as second outgroup to polarize mutations
mkado batch examples/anopheles_batch/ -i gamb -o afun --polarize-match amin
```

### Single Gene Analysis

```bash
# Standard MK test on a single gene
mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun

# Asymptotic test on a single gene
mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a
```

## Expected Output

Running `mkado batch examples/anopheles_batch/ -i gamb -o afun` produces TSV output:

```
gene        Dn   Ds   Pn  Ps  p_value    p_value_adjusted  NI        alpha
AGAP000010  79   176  3   6   1          1                 1.113924  -0.113924
AGAP000011  14   140  2   7   0.21639    0.595072          2.857143  -1.857143
AGAP000012  40   159  5   9   0.179949   0.595072          2.208333  -1.208333
AGAP000014  59   129  5   14  0.796994   1                 0.780872  0.219128
AGAP000016  158  216  8   4   0.136763   0.595072          2.734177  -1.734177
AGAP000021  16   45   0   1   1          1                 0.000000  1.000000
AGAP000022  73   218  3   17  0.423897   0.932574          0.526994  0.473006
AGAP000023  24   128  2   13  1          1                 0.820513  0.179487
AGAP000038  11   145  0   14  0.602424   0.946666          0.000000  1.000000
AGAP000041  60   161  2   10  0.523748   0.946666          0.536667  0.463333
AGAP029936  229  221  40  20  0.0270896  0.297986          1.930131  -0.930131
```

### Visualization Options

```bash
# Volcano plot (p-value vs alpha)
mkado batch examples/anopheles_batch/ -i gamb -o afun --volcano volcano.png

# Asymptotic alpha plot (aggregate)
mkado batch examples/anopheles_batch/ -i gamb -o afun -a --plot-asymptotic asymptotic.png
```

See the [batch workflow documentation](../../docs/batch-workflow.rst) for example plots.

## Notes

- Sequences are filtered by name pattern: `-i gamb` matches all sequences containing "gamb"
- The alignments are codon-aligned (in-frame); use `-r 1` (default) for reading frame
- For polarized tests, use a more distant outgroup (amin) to assign mutations to lineages
