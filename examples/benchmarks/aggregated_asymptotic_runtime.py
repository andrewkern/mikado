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


def write_csv(rows: list[dict], path: Path, extraction_seconds: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        fh.write(f"# extraction_seconds={extraction_seconds:.3f}\n")
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDNAMES)
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
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="darkgrid")
    palette = ["#2c3e50", "#e74c3c", "#3498db", "#e67e22"]

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

    fig, (ax_time, ax_alpha) = plt.subplots(1, 2, figsize=(12, 5))
    _panel(
        ax_time,
        times,
        ylabel="Wall time per replicate (s, log scale)",
        title=f"Asymptotic MK aggregated runtime, {n_genes} random genes / replicate",
        log_y=True,
    )
    _panel(
        ax_alpha,
        alphas,
        ylabel=r"Asymptotic $\alpha$",
        title="Estimated asymptotic alpha",
    )

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
        args.cache = DEFAULT_BENCH_DIR / "cache" / "smoke_poly_data.pkl"
        args.csv = DEFAULT_BENCH_DIR / "results" / "smoke_runtime.csv"
        args.plot = DEFAULT_BENCH_DIR / "figures" / "smoke_runtime.png"
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
                )
            )

    write_csv(rows, args.csv, extraction_seconds)
    if not args.no_plot:
        make_figure(
            rows,
            path=args.plot,
            extraction_seconds=extraction_seconds,
            n_genes=args.n_genes,
            n_total_genes=len(poly_data),
            workers=extraction_workers,
        )


if __name__ == "__main__":
    main()
