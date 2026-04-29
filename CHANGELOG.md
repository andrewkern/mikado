# Changelog

## [Unreleased]

### Added
- Nei-Gojobori site totals on every result dataclass (`MKResult`,
  `AsymptoticMKResult`, `PolarizedMKResult`, `ImputedMKResult`,
  `AlphaTGResult`) as `ln` and `ls` fields. Aggregated and per-gene
  paths populate them automatically from analyzed codons.
- `omega` (dN/dS, with Nei-Gojobori site weighting) on every result
  dataclass: `omega = (Dn/Ds) * (Ls/Ln)`. Returns `None` when `Ds = 0`
  or site counts are unavailable.
- `omega_a` and `omega_na` (adaptive / non-adaptive substitution rates)
  on result types whose alpha estimator is appropriate for the
  decomposition: `AsymptoticMKResult` (per-gene and aggregated),
  `ImputedMKResult`, and `AlphaTGResult`. Following Gossmann, Keightley
  & Eyre-Walker (2012, *Genome Biol Evol* 4(5):658-667),
  `omega_a = alpha * omega` and `omega_na = (1 - alpha) * omega`.
  `MKResult` and `PolarizedMKResult` deliberately omit these — the
  per-gene Smith-Eyre-Walker / polarized alpha is too noisy for a
  meaningful rate decomposition. See `docs/omega.rst` for the
  rationale.
- 95% confidence intervals on the omega decomposition where a sampling
  procedure is already available. `AsymptoticMKResult` exposes
  `omega_a_ci_low`/`omega_a_ci_high` and `omega_na_ci_low`/`omega_na_ci_high`
  by scaling the existing alpha CI (Dn/Ds/Ln/Ls are constants under the
  asymptotic Monte Carlo, so omega itself has no sampling distribution).
  `AlphaTGResult` recomputes omega/omega_a/omega_na per gene-resampled
  bootstrap replicate and reports `omega_ci_low/high`,
  `omega_a_ci_low/high`, `omega_na_ci_low/high` from the percentile
  distribution. `ImputedMKResult` does not (yet) bootstrap and reports
  point estimates only.
- `omega_decomposition()` helper in `mkado.analysis.statistics`.
- `AlignedPair.count_total_sites()` aggregating Nei-Gojobori `(Ln, Ls)`
  totals over analyzed codons.
- `Ln`, `Ls`, `omega`, `omega_a`, and `omega_na` columns in pretty,
  TSV, and JSON output for all result types.
- `--ci-method {monte-carlo,bootstrap}` CLI flag on `test`, `batch`, and
  `vcf` commands. The new `bootstrap` mode uses case-resampling of the
  pooled polymorphism list with replacement, refitting the curve per
  replicate (closes #10). Default remains `monte-carlo` to preserve
  existing behavior. Each result type now records the CI method via a
  `ci_method` field on `AsymptoticMKResult` (default `"monte-carlo"`),
  `AlphaTGResult` (`"bootstrap"`), and `ImputedMKResult`
  (`"bootstrap"` when set, `None` when no CI was computed).
- Bootstrap CI on imputed-MK alpha. `imputed_mk_test` and
  `imputed_mk_test_multi` now accept `n_bootstrap: int = 0` (preserves
  legacy behavior when 0) and `seed: int = 42`. CLI: enabled
  automatically when `--bootstrap > 0` (the existing flag drives the
  replicate count). New fields on `ImputedMKResult`:
  `alpha_ci_low/high`, `omega_a_ci_low/high`, `omega_na_ci_low/high`,
  `ci_method`. Following the asymptotic CI pattern, omega-decomposition
  CIs are derived analytically from the alpha CI scaled by omega.
- `_compute_ci_bootstrap()` helper in `mkado.analysis.asymptotic`
  alongside the existing `_compute_ci_monte_carlo()`.
- `_bootstrap_imputed_alpha()` helper in `mkado.analysis.imputed`.

### Performance
- `_bootstrap_imputed_alpha()` is now vectorized: pre-extracts numpy
  arrays of frequencies and N/S type once, then per replicate uses
  `rng.integers` indexing and `np.count_nonzero` instead of calling
  `_compute_imputed` (which does four Python-level `sum()` passes per
  call). ~63× speedup on the imputed-bootstrap path measured on a
  ~1000-polymorphism gene with 1000 replicates (1.4 s → 22 ms). CIs are
  byte-identical to the prior implementation under the same seed
  (closes #16).
- Per-gene `asymptotic_mk_test` bootstrap loop replaced by a call to
  the shared `_compute_ci_bootstrap()` helper (closes #17). Removes
  ~60 lines of duplicate code, vectorizes the per-replicate binning,
  and aligns failure-handling between per-gene and aggregated paths.

### Fixed
- `GeneticCode.translate` no longer uses `@lru_cache` on a bound method
  (closes #14). The decorator's cache key included `self`, pinning every
  `GeneticCode` instance for the lifetime of the class object — a real
  leak in tests, notebooks, and any process that constructs many
  ephemeral instances. The translation table is now eagerly precomputed
  in `__init__` as a per-instance dict (~640 bytes), mirroring the
  proven `_compute_codon_paths` pattern. Same fix applied to
  `count_synonymous_sites`, retiring its lazy-population dict for an
  eager precompute over the bounded 64-codon keyspace. Behavior is
  byte-identical for callers; a new regression test in
  `tests/test_codons.py::TestGeneticCodeMemory` uses `weakref` to verify
  instances are garbage-collected after `del` + `gc.collect()`.
- Standardized missing-value rendering on the literal `"NA"` (no slash)
  across pretty and TSV output (closes #13). The dataclass `__str__`
  methods previously rendered `"N/A"` while the TSV formatters used
  `"NA"` — a user piping pretty output through R/awk that parses `"NA"`
  would have silently missed `"N/A"` rows. JSON renders `null`,
  unchanged. New regression suite in `tests/test_output_format.py`
  guards both the convention and the absence of `"N/A"` in any output
  format.

### Changed
- Per-gene `asymptotic_mk_test` CI failure handling: replicates whose
  curve fit fails are now **dropped** rather than imputed with the last
  bin's alpha (the prior behavior biased the CI when fits failed
  systematically). When fewer than half of replicates succeed, the CI
  degenerates to the point estimate. Matches the aggregated path's
  behavior added in #10.

## [0.4.0] - 2026-03-17

### Added
- `mkado vcf` command: run MK tests directly from VCF + reference FASTA + GFF3 annotation
- Chunked parallel VCF processing for improved performance
- GFF3 parser with GTF format detection, malformed line handling, and hierarchy validation warnings
- `--verbose` flag for VCF mode to surface htslib, GFF parsing, and volcano plot diagnostics
- Example VCF dataset in `examples/example_vcf/`
- VCF input documentation (`docs/vcf-input.rst`)

### Fixed
- `cds_length` property called as method in GFF parser logging

## [0.3.0] - 2026-03-10

### Added
- Alternate genetic code support (`--code-table`) for 24 NCBI genetic code tables
- Name-based code selection (e.g. `--code-table vertebrate-mito`)
- `mkado codes` command to list available genetic code tables

## [0.2.0] - 2026-03-01

### Added
- Imputed MK test (Murga-Moreno et al. 2022) for correcting slightly deleterious mutations by imputation rather than discarding low-frequency polymorphisms
- `--imputed` flag for `test` and `batch` CLI commands
- Aggregated and per-gene batch modes for imputed MK test
- DFE decomposition (d, b, f fractions) when site counts are provided
- Documentation for the imputed MK test

## [0.1.2] - 2025-05-14

### Fixed
- Use delta method for `min_frequency` to correctly exclude singletons
- Remove title from volcano plots

### Added
- PyPI publish workflow (triggered by version tags)

## [0.1.1] - 2025-05-13

### Added
- Direction of Selection (DoS) statistic
- Tarone-Greenland alpha (α_TG) weighted multi-gene estimator
- `--no-singletons` option for automatic singleton exclusion
- Polarized MK test polymorphism polarization
- Volcano plot visualization for batch results
- Asymptotic alpha(x) plot visualization
- Benjamini-Hochberg adjusted p-values in batch output
- Sphinx documentation

## [0.1.0] - 2025-01-22

- Initial release

[0.4.0]: https://github.com/kr-colab/mkado/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kr-colab/mkado/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kr-colab/mkado/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/kr-colab/mkado/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/kr-colab/mkado/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/kr-colab/mkado/releases/tag/v0.1.0
