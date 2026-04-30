"""Aggregated-genes asymptotic-MK runtime benchmark.

Mirrors the regime in which the asymptotic test is statistically reliable:
aggregate ~1,000 random genes per replicate, repeat over many replicates
(Murga-Moreno et al. 2022, doi:10.1093/g3journal/jkac206). Sweeps the 2x2
grid of ``sfs_mode`` x ``ci_method`` so the four implementation choices
that matter for runtime are characterized in one figure.

Polymorphism extraction over the full ``--alignments`` directory is the
dominant cost; it runs once and is cached to disk. Per-replicate timing
covers the ``asymptotic_mk_test_aggregated`` call only -- aggregation +
curve fit + CI computation, end-to-end. The one-time extraction wall
time is recorded in the CSV header and the figure caption so the
"what's measured" question is unambiguous.

Run end-to-end on the lab machine:

    uv run python examples/benchmarks/aggregated_asymptotic_runtime.py \\
        --workers 8

Smoke-run on the 5-gene checked-in subset:

    uv run python examples/benchmarks/aggregated_asymptotic_runtime.py --smoke
"""

from __future__ import annotations

import argparse
import csv
import logging
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
from mkado.analysis.asymptotic import (
    PolymorphismData,
    asymptotic_mk_test_aggregated,
    extract_polymorphism_data,
)
from mkado.cli import find_alignment_files
from mkado.core.sequences import SequenceSet

logger = logging.getLogger("mkado.benchmark")

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ALIGNMENTS = Path("/sietch_colab/adkern/human_mk_data/alignments")
SMOKE_ALIGNMENTS = REPO_ROOT / "data" / "mk_human_chimp"
DEFAULT_BENCH_DIR = REPO_ROOT / "examples" / "benchmarks"

CSV_FIELDNAMES = (
    "replicate_idx",
    "sfs_mode",
    "ci_method",
    "n_genes",
    "wall_time_seconds",
    "alpha_asymptotic",
    "ci_low",
    "ci_high",
    "num_bins",
    "ci_replicates",
)

ALPHA_TG_CSV_FIELDNAMES = (
    "replicate_idx",
    "n_genes",
    "wall_time_seconds",
    "alpha_tg",
    "ci_low",
    "ci_high",
    "bootstrap_replicates",
)

SCALING_CSV_FIELDNAMES = (
    "n_workers",
    "rep_idx",
    "extraction_seconds",
    "asymptotic_seconds",
    "total_seconds",
)


@dataclass
class _ExtractionTask:
    file_path: Path
    ingroup_match: str
    outgroup_match: str


def _extract_one(task: _ExtractionTask) -> PolymorphismData | None:
    """Load a combined-file FASTA and return PolymorphismData for that gene.

    Mirrors the combined-file branch of ``mkado.batch_workers.process_gene``
    so the benchmark hits the same code path the CLI does.
    """
    try:
        all_seqs = SequenceSet.from_fasta(task.file_path)
        ingroup = all_seqs.filter_by_name(task.ingroup_match)
        outgroup = all_seqs.filter_by_name(task.outgroup_match)
        if len(ingroup) == 0 or len(outgroup) == 0:
            return None
        return extract_polymorphism_data(
            ingroup=ingroup,
            outgroup=outgroup,
            gene_id=task.file_path.stem,
        )
    except Exception:
        logger.exception("extraction failed for %s", task.file_path.name)
        return None


def load_or_extract(
    alignment_dir: Path,
    ingroup_match: str,
    outgroup_match: str,
    cache_path: Path,
    workers: int,
) -> tuple[list[PolymorphismData], float, int]:
    """Return (polymorphism data, extraction seconds, workers used).

    On cache hit the wall time and worker count are the values recorded
    on the original extraction run, not the current invocation -- so the
    figure caption stays internally consistent.
    """
    if cache_path.exists():
        logger.info("loading polymorphism cache from %s", cache_path)
        with cache_path.open("rb") as fh:
            payload = pickle.load(fh)
        return payload["data"], payload["extraction_seconds"], payload.get("workers", workers)

    files = find_alignment_files(alignment_dir)
    if not files:
        sys.exit(f"no FASTA files found under {alignment_dir}")

    logger.info("extracting polymorphism data from %d genes (workers=%d)", len(files), workers)
    tasks = [
        _ExtractionTask(f, ingroup_match=ingroup_match, outgroup_match=outgroup_match)
        for f in files
    ]
    start = time.perf_counter()
    data: list[PolymorphismData] = []
    if workers <= 1:
        for t in tasks:
            poly = _extract_one(t)
            if poly is not None:
                data.append(poly)
    else:
        # ~4 chunks per worker balances IPC overhead against load imbalance from
        # variable-size FASTAs; matches the multiprocessing.Pool default heuristic.
        chunksize = max(len(tasks) // (workers * 4), 1)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for poly in pool.map(_extract_one, tasks, chunksize=chunksize):
                if poly is not None:
                    data.append(poly)
    extraction_seconds = time.perf_counter() - start

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("wb") as fh:
        pickle.dump(
            {"data": data, "extraction_seconds": extraction_seconds, "workers": workers},
            fh,
        )
    logger.info(
        "extracted %d genes in %.1fs; cached to %s",
        len(data),
        extraction_seconds,
        cache_path,
    )
    return data, extraction_seconds, workers


def run_replicates(
    poly_data: list[PolymorphismData],
    *,
    sfs_mode: str,
    ci_method: str,
    n_genes: int,
    n_replicates: int,
    base_seed: int,
    num_bins: int,
    workers: int = 1,
) -> list[dict]:
    """Sample n_genes without replacement and time the aggregated test.

    The same per-replicate seed is used across the (sfs_mode, ci_method)
    grid so each cell sees the same gene subsets -- across-condition
    variance is method-driven, not sample-driven.
    """
    n_total = len(poly_data)
    if n_total < n_genes:
        sys.exit(f"only {n_total} genes available; cannot sample {n_genes} without replacement")

    # Realistic CLI defaults: MC samples params from a covariance matrix
    # (cheap, ~10000x); bootstrap refits per replicate (~100x).
    ci_replicates = 10000 if ci_method == "monte-carlo" else 100

    rows: list[dict] = []
    for replicate_idx in range(n_replicates):
        seed = base_seed + replicate_idx
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_total, size=n_genes, replace=False)
        sample = [poly_data[i] for i in idx]

        start = time.perf_counter()
        result = asymptotic_mk_test_aggregated(
            gene_data=sample,
            num_bins=num_bins,
            ci_replicates=ci_replicates,
            ci_method=ci_method,
            seed=seed,
            sfs_mode=sfs_mode,
            workers=workers,
        )
        wall = time.perf_counter() - start

        rows.append(
            {
                "replicate_idx": replicate_idx,
                "sfs_mode": sfs_mode,
                "ci_method": ci_method,
                "n_genes": n_genes,
                "wall_time_seconds": wall,
                "alpha_asymptotic": result.alpha_asymptotic,
                "ci_low": result.ci_low,
                "ci_high": result.ci_high,
                "num_bins": num_bins,
                "ci_replicates": ci_replicates,
            }
        )
        if (replicate_idx + 1) % 10 == 0 or replicate_idx == n_replicates - 1:
            logger.info(
                "[%s/%s] replicate %d/%d  alpha=%.4f  wall=%.3fs",
                sfs_mode,
                ci_method,
                replicate_idx + 1,
                n_replicates,
                result.alpha_asymptotic,
                wall,
            )
    return rows


def run_replicates_alpha_tg(
    poly_data: list[PolymorphismData],
    *,
    n_genes: int,
    n_replicates: int,
    base_seed: int,
    bootstrap_replicates: int,
) -> list[dict]:
    """Per-replicate Tarone-Greenland alpha on the same 1,000-gene samples.

    Uses identical per-replicate seeds as ``run_replicates`` so each row in
    this CSV pairs 1:1 with the four asymptotic rows for the same sample
    -- the manuscript can show all five alpha estimators on a common axis.
    """
    n_total = len(poly_data)
    if n_total < n_genes:
        sys.exit(f"only {n_total} genes available; cannot sample {n_genes} without replacement")

    rows: list[dict] = []
    for replicate_idx in range(n_replicates):
        seed = base_seed + replicate_idx
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_total, size=n_genes, replace=False)
        sample = [poly_data[i] for i in idx]

        start = time.perf_counter()
        result = alpha_tg_from_gene_data(
            gene_data=sample,
            bootstrap_replicates=bootstrap_replicates,
            seed=seed,
        )
        wall = time.perf_counter() - start

        rows.append(
            {
                "replicate_idx": replicate_idx,
                "n_genes": n_genes,
                "wall_time_seconds": wall,
                "alpha_tg": result.alpha_tg,
                "ci_low": result.ci_low,
                "ci_high": result.ci_high,
                "bootstrap_replicates": bootstrap_replicates,
            }
        )
        if (replicate_idx + 1) % 10 == 0 or replicate_idx == n_replicates - 1:
            logger.info(
                "[alpha_tg] replicate %d/%d  alpha=%.4f  wall=%.3fs",
                replicate_idx + 1,
                n_replicates,
                result.alpha_tg,
                wall,
            )
    return rows


def run_scaling_sweep(
    *,
    alignment_dir: Path,
    ingroup_match: str,
    outgroup_match: str,
    cache_path: Path,
    worker_counts: list[int],
    n_reps: int,
    n_genes: int,
    n_replicates: int,
    base_seed: int,
    num_bins: int,
) -> list[dict]:
    """End-to-end pipeline timing across worker counts.

    For each (n_workers, rep), deletes the cache, runs a fresh extraction,
    and runs the full 4-cell asymptotic grid. Records extraction wall time,
    asymptotic-runs wall time, and total -- so the manuscript can show how
    the parallelizable upstream stage scales while the asymptotic stage
    stays sequential.
    """
    rows: list[dict] = []
    for n_workers in worker_counts:
        for rep_idx in range(n_reps):
            cache_path.unlink(missing_ok=True)
            logger.info("[scaling] workers=%d rep=%d/%d", n_workers, rep_idx + 1, n_reps)
            poly_data, extraction_seconds, _ = load_or_extract(
                alignment_dir=alignment_dir,
                ingroup_match=ingroup_match,
                outgroup_match=outgroup_match,
                cache_path=cache_path,
                workers=n_workers,
            )
            asym_start = time.perf_counter()
            for sfs_mode in ("at", "above"):
                for ci_method in ("monte-carlo", "bootstrap"):
                    run_replicates(
                        poly_data,
                        sfs_mode=sfs_mode,
                        ci_method=ci_method,
                        n_genes=n_genes,
                        n_replicates=n_replicates,
                        base_seed=base_seed,
                        num_bins=num_bins,
                        workers=n_workers,
                    )
            asymptotic_seconds = time.perf_counter() - asym_start
            total_seconds = extraction_seconds + asymptotic_seconds
            rows.append(
                {
                    "n_workers": n_workers,
                    "rep_idx": rep_idx,
                    "extraction_seconds": extraction_seconds,
                    "asymptotic_seconds": asymptotic_seconds,
                    "total_seconds": total_seconds,
                }
            )
            logger.info(
                "[scaling] workers=%d rep=%d  extract=%.1fs  asym=%.1fs  total=%.1fs",
                n_workers,
                rep_idx + 1,
                extraction_seconds,
                asymptotic_seconds,
                total_seconds,
            )
    return rows


def write_csv(rows: list[dict], path: Path, extraction_seconds: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        fh.write(f"# extraction_seconds={extraction_seconds:.3f}\n")
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    logger.info("wrote %d rows to %s", len(rows), path)


def write_simple_csv(rows: list[dict], path: Path, fieldnames: tuple[str, ...]) -> None:
    """Write a CSV without the leading extraction_seconds comment row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    logger.info("wrote %d rows to %s", len(rows), path)


def make_figure(
    rows: list[dict],
    path: Path,
    extraction_seconds: float,
    n_genes: int,
    n_total_genes: int,
    workers: int,
    alpha_tg_rows: list[dict] | None = None,
    scaling_rows: list[dict] | None = None,
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="darkgrid")
    palette = ["#2c3e50", "#e74c3c", "#3498db", "#e67e22"]
    alpha_tg_color = "#27ae60"

    times = np.array([r["wall_time_seconds"] for r in rows])
    alphas = np.array([r["alpha_asymptotic"] for r in rows])
    sfs_modes = np.array([r["sfs_mode"] for r in rows])
    ci_methods = np.array([r["ci_method"] for r in rows])
    conditions = np.array([f"{s}\n{c}" for s, c in zip(sfs_modes, ci_methods)])
    order = sorted(set(conditions.tolist()))

    def _panel(ax, y, *, ylabel: str, title: str, log_y: bool = False) -> None:
        sns.boxplot(
            x=conditions,
            y=y,
            hue=conditions,
            order=order,
            hue_order=order,
            ax=ax,
            palette=palette,
            legend=False,
        )
        sns.stripplot(
            x=conditions,
            y=y,
            order=order,
            ax=ax,
            color="white",
            edgecolor="#2c3e50",
            linewidth=0.5,
            alpha=0.6,
            size=3,
        )
        if log_y:
            ax.set_yscale("log")
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlabel("sfs_mode / ci_method", fontsize=12)
        ax.set_title(title, fontsize=13, fontweight="bold")

    n_panels = 2 + (1 if scaling_rows else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    axes_iter = iter(axes if n_panels > 1 else [axes])

    if scaling_rows:
        ax_scaling = next(axes_iter)
        sc_workers = np.array([r["n_workers"] for r in scaling_rows], dtype=int)
        sc_total = np.array([r["total_seconds"] for r in scaling_rows])
        sc_extract = np.array([r["extraction_seconds"] for r in scaling_rows])
        sc_asym = np.array([r["asymptotic_seconds"] for r in scaling_rows])
        sc_order = sorted(set(sc_workers.tolist()))
        pos_map = {w: i for i, w in enumerate(sc_order)}
        positions = list(range(len(sc_order)))

        ax_scaling.boxplot(
            [sc_total[sc_workers == w] for w in sc_order],
            positions=positions,
            widths=0.5,
            patch_artist=True,
            boxprops={"facecolor": palette[0], "alpha": 0.4, "edgecolor": palette[0]},
            medianprops={"color": "white", "linewidth": 1.5},
            whiskerprops={"color": palette[0]},
            capprops={"color": palette[0]},
            showfliers=False,
        )
        sc_offsets = np.array([pos_map[w] for w in sc_workers], dtype=float)
        jitter = np.random.default_rng(0).uniform(-0.08, 0.08, size=len(sc_offsets))
        ax_scaling.scatter(
            sc_offsets + jitter,
            sc_extract,
            color=palette[1],
            edgecolor="white",
            linewidth=0.4,
            s=40,
            label="extraction",
            zorder=3,
        )
        ax_scaling.scatter(
            sc_offsets + jitter,
            sc_asym,
            color=palette[2],
            edgecolor="white",
            linewidth=0.4,
            s=40,
            label="asymptotic runs",
            zorder=3,
        )
        ax_scaling.set_xticks(positions)
        ax_scaling.set_xticklabels([str(w) for w in sc_order])
        ax_scaling.set_yscale("log")
        ax_scaling.set_ylabel("End-to-end wall time (s, log scale)", fontsize=12)
        ax_scaling.set_xlabel("Workers (extraction parallelism)", fontsize=12)
        ax_scaling.set_title("Pipeline runtime vs worker count", fontsize=13, fontweight="bold")
        ax_scaling.legend(loc="upper right", fontsize=9)

    ax_time = next(axes_iter)
    _panel(
        ax_time,
        times,
        ylabel="Wall time per replicate (s, log scale)",
        title=f"Asymptotic MK aggregated runtime, {n_genes} random genes / replicate",
        log_y=True,
    )

    ax_alpha = next(axes_iter)
    _panel(
        ax_alpha,
        alphas,
        ylabel=r"Asymptotic $\alpha$",
        title="Estimated asymptotic alpha",
    )
    if alpha_tg_rows:
        # Overlay alpha_TG point estimates at a fifth, slightly-offset x position.
        tg_alphas = np.array([r["alpha_tg"] for r in alpha_tg_rows])
        tg_x = np.full(len(tg_alphas), len(order))
        ax_alpha.scatter(
            tg_x + np.random.default_rng(0).uniform(-0.15, 0.15, size=len(tg_x)),
            tg_alphas,
            color=alpha_tg_color,
            edgecolor="white",
            linewidth=0.4,
            s=30,
            alpha=0.7,
            zorder=3,
            label=r"$\alpha_{TG}$",
        )
        # Boxplot for alpha_TG at the same x position.
        ax_alpha.boxplot(
            tg_alphas,
            positions=[len(order)],
            widths=0.5,
            patch_artist=True,
            boxprops={"facecolor": alpha_tg_color, "alpha": 0.4, "edgecolor": alpha_tg_color},
            medianprops={"color": "white", "linewidth": 1.5},
            whiskerprops={"color": alpha_tg_color},
            capprops={"color": alpha_tg_color},
            showfliers=False,
        )
        labels = list(order) + [r"$\alpha_{TG}$"]
        ax_alpha.set_xticks(range(len(labels)))
        ax_alpha.set_xticklabels(labels)
        ax_alpha.set_xlabel("estimator / condition", fontsize=12)
        ax_alpha.legend(loc="upper right", fontsize=9)

    fig.suptitle(
        (
            f"Extraction (parsing + site counting + SFS construction) of "
            f"{n_total_genes} genes amortized across replicates: "
            f"{extraction_seconds:.1f}s on {workers} workers"
        ),
        fontsize=10,
        y=-0.02,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote figure to %s", path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--alignments",
        type=Path,
        default=DEFAULT_ALIGNMENTS,
        help="Directory of per-gene multi-species FASTA alignments.",
    )
    parser.add_argument("--ingroup-match", default="Homo_sapiens")
    parser.add_argument("--outgroup-match", default="Pan_troglodytes")
    parser.add_argument("--n-genes", type=int, default=1000)
    parser.add_argument("--n-replicates", type=int, default=100)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--num-bins", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--cache",
        type=Path,
        default=DEFAULT_BENCH_DIR / "cache" / "human_mk_poly_data.pkl",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_BENCH_DIR / "results" / "aggregated_asymptotic_runtime.csv",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=DEFAULT_BENCH_DIR / "figures" / "aggregated_asymptotic_runtime.png",
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--no-alpha-tg",
        action="store_true",
        help="Skip the per-replicate alpha_TG comparison run.",
    )
    parser.add_argument(
        "--alpha-tg-bootstrap",
        type=int,
        default=200,
        help=(
            "Bootstrap replicates for alpha_TG CI per sample (default 200). "
            "alpha_tg_from_gene_data's library default is 1000; 200 is plenty "
            "for the manuscript figure and keeps the comparison cheap."
        ),
    )
    parser.add_argument(
        "--scaling",
        action="store_true",
        help=(
            "Run the worker-count scaling sweep (end-to-end pipeline timing "
            "at each --scaling-workers level, --scaling-reps reps each). "
            "Adds an extra panel to the figure on the left."
        ),
    )
    parser.add_argument(
        "--scaling-workers",
        default="8,16,32,64",
        help="Comma-separated worker counts for the scaling sweep.",
    )
    parser.add_argument("--scaling-reps", type=int, default=3)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help=(
            "Run against the 5-gene checked-in subset with n_genes=3 and "
            "n_replicates=2 to verify the script wires up correctly."
        ),
    )
    args = parser.parse_args(argv)

    if args.smoke:
        args.alignments = SMOKE_ALIGNMENTS
        args.n_genes = 3
        args.n_replicates = 2
        args.workers = 2
        args.alpha_tg_bootstrap = 50
        args.scaling_workers = "1,2"
        args.scaling_reps = 2
        args.cache = DEFAULT_BENCH_DIR / "cache" / "smoke_poly_data.pkl"
        args.csv = DEFAULT_BENCH_DIR / "results" / "smoke_runtime.csv"
        args.plot = DEFAULT_BENCH_DIR / "figures" / "smoke_runtime.png"

    args.scaling_workers_list = [int(w) for w in args.scaling_workers.split(",") if w.strip()]
    return args


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    args = parse_args(argv)

    if not args.alignments.is_dir():
        sys.exit(f"alignments dir not found: {args.alignments}")

    poly_data, extraction_seconds, extraction_workers = load_or_extract(
        alignment_dir=args.alignments,
        ingroup_match=args.ingroup_match,
        outgroup_match=args.outgroup_match,
        cache_path=args.cache,
        workers=args.workers,
    )

    logger.info(
        "running %d replicates x 4 conditions on %d-gene samples (cache size: %d genes)",
        args.n_replicates,
        args.n_genes,
        len(poly_data),
    )

    rows: list[dict] = []
    for sfs_mode in ("at", "above"):
        for ci_method in ("monte-carlo", "bootstrap"):
            rows.extend(
                run_replicates(
                    poly_data,
                    sfs_mode=sfs_mode,
                    ci_method=ci_method,
                    n_genes=args.n_genes,
                    n_replicates=args.n_replicates,
                    base_seed=args.seed,
                    num_bins=args.num_bins,
                    workers=args.workers,
                )
            )

    alpha_tg_rows: list[dict] = []
    if not args.no_alpha_tg:
        alpha_tg_rows = run_replicates_alpha_tg(
            poly_data,
            n_genes=args.n_genes,
            n_replicates=args.n_replicates,
            base_seed=args.seed,
            bootstrap_replicates=args.alpha_tg_bootstrap,
        )

    scaling_rows: list[dict] = []
    if args.scaling:
        # The sweep clobbers the cache for each (workers, rep) run -- save
        # the current cache aside so the post-sweep CSV/figure work above
        # has its data, then restore it after the sweep.
        if args.cache.exists():
            saved_cache = args.cache.with_suffix(args.cache.suffix + ".saved")
            args.cache.rename(saved_cache)
        else:
            saved_cache = None
        try:
            scaling_rows = run_scaling_sweep(
                alignment_dir=args.alignments,
                ingroup_match=args.ingroup_match,
                outgroup_match=args.outgroup_match,
                cache_path=args.cache,
                worker_counts=args.scaling_workers_list,
                n_reps=args.scaling_reps,
                n_genes=args.n_genes,
                n_replicates=args.n_replicates,
                base_seed=args.seed,
                num_bins=args.num_bins,
            )
        finally:
            args.cache.unlink(missing_ok=True)
            if saved_cache is not None:
                saved_cache.rename(args.cache)

    write_csv(rows, args.csv, extraction_seconds)
    if alpha_tg_rows:
        alpha_tg_path = args.csv.with_name(args.csv.stem + "_alpha_tg.csv")
        write_simple_csv(alpha_tg_rows, alpha_tg_path, ALPHA_TG_CSV_FIELDNAMES)
    if scaling_rows:
        scaling_path = args.csv.with_name(args.csv.stem + "_scaling.csv")
        write_simple_csv(scaling_rows, scaling_path, SCALING_CSV_FIELDNAMES)

    if not args.no_plot:
        make_figure(
            rows,
            path=args.plot,
            extraction_seconds=extraction_seconds,
            n_genes=args.n_genes,
            n_total_genes=len(poly_data),
            workers=extraction_workers,
            alpha_tg_rows=alpha_tg_rows or None,
            scaling_rows=scaling_rows or None,
        )


if __name__ == "__main__":
    main()
