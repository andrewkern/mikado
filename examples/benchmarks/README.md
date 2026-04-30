# mkado benchmarks

Reproducible runtime benchmarks for the mkado test suite.

## Aggregated asymptotic runtime

`aggregated_asymptotic_runtime.py` characterizes wall time and alpha
estimates for `asymptotic_mk_test_aggregated` in the regime where the
test is statistically reliable: aggregate ~1,000 random genes per
replicate, repeat over many replicates (per Murga-Moreno et al. 2022,
[doi](https://doi.org/10.1093/g3journal/jkac206)).

The script sweeps the full 2x2 grid of:

- `sfs_mode` in `{at, above}` — per-bin SFS (Messer & Petrov 2013) vs.
  inclusive cumulative SFS (Uricchio et al. 2019).
- `ci_method` in `{monte-carlo, bootstrap}` — parametric covariance
  sampling vs. case-resampling refit.

### What is measured

Per-replicate wall time covers the **`asymptotic_mk_test_aggregated`
call only** — aggregation, curve fit, and CI computation,
end-to-end for that step. Polymorphism extraction (FASTA parsing,
codon-aligning the ingroup vs. outgroup, classifying nonsynonymous /
synonymous polymorphism, Nei-Gojobori site counting) runs **once**
across all genes in the input directory, is cached to disk, and is
**amortized** across all replicates. The one-time extraction wall time
is recorded in the CSV header (`# extraction_seconds=...`) and printed
on the figure caption, so the cost composition is unambiguous.

The `--scaling` flag adds a **worker-count sweep** to the run: at each
worker count (default `8,16,32,64`), the script runs the full pipeline
(fresh extraction + 400 timed asymptotic runs) `--scaling-reps` times
(default 3). This produces a third figure panel showing how extraction
parallelism and bootstrap-CI parallelism each scale with workers.

### Inputs

- `--alignments DIR`: a directory of per-gene combined-file FASTAs,
  each with ingroup and outgroup haplotypes named so they can be
  filtered by substring. Default points at the lab dataset of
  ~12,437 human protein-coding genes
  (`Homo_sapiens` ingroup, `Pan_troglodytes` outgroup).
- `--ingroup-match`, `--outgroup-match`: substring matchers passed to
  `SequenceSet.filter_by_name`. Defaults match the lab dataset.

### Outputs

- `results/aggregated_asymptotic_runtime.csv` —
  `replicate_idx, sfs_mode, ci_method, n_genes, wall_time_seconds,
  alpha_asymptotic, ci_low, ci_high, num_bins, ci_replicates`. Header
  comment line records the one-time extraction wall time.
- `results/aggregated_asymptotic_runtime_scaling.csv` — only when
  `--scaling` is passed: end-to-end pipeline timings at each worker
  count, broken down into extraction and asymptotic-runs phases.
- `figures/aggregated_asymptotic_runtime.png` — two-panel boxplot of
  wall time (log-scale y) and asymptotic alpha across the four
  conditions. Adds a scaling panel on the left when `--scaling` data
  is present.
- `cache/human_mk_poly_data.pkl` — pickled extraction output. Re-runs
  reuse it; delete it to force a re-extract.

The `results/`, `figures/`, and `cache/` directories are gitignored.

### Running

Full benchmark, lab machine, 8 workers:

```bash
uv run python examples/benchmarks/aggregated_asymptotic_runtime.py --workers 8
```

Defaults: 100 replicates per condition x 4 conditions = 400 timed runs.
Each replicate samples 1,000 genes without replacement; the same
per-replicate seeds are reused across the grid so across-condition
variance is method-driven, not sample-driven.

Full benchmark with worker-count scaling sweep (~50 min on the lab
machine: 4 worker counts × 3 fresh extractions × ~4 min each):

```bash
uv run python examples/benchmarks/aggregated_asymptotic_runtime.py --scaling
```

Smoke run (5-gene checked-in subset, 3 genes/replicate, 2 replicates;
verifies plumbing without committing to the full run):

```bash
uv run python examples/benchmarks/aggregated_asymptotic_runtime.py --smoke
```

### Reproducibility

The base seed is `42` (override with `--seed`). Per-replicate seeds are
`base_seed + replicate_idx`. Same-seed runs of the same dataset are
deterministic, modulo PYTHONHASHSEED-dependent per-gene polymorphism
ordering during extraction (set `PYTHONHASHSEED=0` if you need byte-exact
reproducibility across machines).
