#!/usr/bin/env python3
"""Validate mkado DFE implementation against GRAPES reference implementation.

This script compares mkado's DFE-based alpha estimation against GRAPES (Galtier 2016),
the reference C++ implementation. It creates test datasets by aggregating Anopheles
gene data, runs both implementations, and verifies they produce equivalent results.

External Dependencies
---------------------
This script requires the GRAPES binary, which is NOT included in the repository.
You must build it from source before running this validation.

Building GRAPES
---------------
1. Clone GRAPES into the vendor directory:

   cd vendor
   git clone https://github.com/BioPP/grapes.git

2. Build GRAPES using CMake:

   cd grapes
   cmake -B build
   cmake --build build

3. The binary will be at: vendor/grapes/build/grapes/grapes

GRAPES Reference
----------------
- Paper: Galtier N (2016) Adaptive Protein Evolution in Animals and the Effective
  Population Size Hypothesis. PLoS Genet 12(1): e1005774
- Repository: https://github.com/BioPP/grapes

Usage
-----
Run from the repository root:

    uv run python scripts/validate_grapes_comparison.py

Expected Output
---------------
The script tests GammaZero and GammaExpo models on three dataset sizes (50, 150,
400 genes). All alpha values should match within 0.01 (1% tolerance). Typical
results show mkado is 15-35x faster than GRAPES while producing equivalent results.
"""

import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mkado.analysis.asymptotic import extract_polymorphism_data
from mkado.analysis.dfe import polymorphism_data_to_dfe_input
from mkado.analysis.dfe_models import DFEInput, PrecomputedData, get_model
from mkado.analysis.dfe_validation import export_to_dofe
from mkado.core.sequences import SequenceSet

# Path to GRAPES binary (not committed to repo)
GRAPES_BINARY = Path(__file__).parent.parent / "vendor" / "grapes" / "build" / "grapes" / "grapes"

# Path to Anopheles test data
ANOPHELES_DATA_DIR = Path(__file__).parent.parent / "examples" / "anopheles_batch"


def get_anopheles_files() -> list[Path]:
    """Get list of Anopheles FASTA files from examples directory."""
    if not ANOPHELES_DATA_DIR.exists():
        return []
    return sorted(ANOPHELES_DATA_DIR.glob("*.fa"))


def extract_ingroup_outgroup(
    fasta_path: Path,
    n_samples: int = 20,
) -> tuple[SequenceSet | None, SequenceSet | None]:
    """Extract ingroup (gamb) and outgroup (afun) sequences from Anopheles file.

    Args:
        fasta_path: Path to FASTA file containing aligned sequences
        n_samples: Number of ingroup samples to use

    Returns:
        Tuple of (ingroup, outgroup) SequenceSets, or (None, None) if insufficient data
    """
    try:
        all_seqs = SequenceSet.from_fasta(fasta_path)

        # Split by species: A. gambiae (ingroup) vs A. funestus (outgroup)
        gamb_seqs = [s for s in all_seqs if "_gamb_" in s.name]
        afun_seqs = [s for s in all_seqs if "_afun_" in s.name]

        if len(gamb_seqs) < n_samples or len(afun_seqs) < 1:
            return None, None

        ingroup = SequenceSet(gamb_seqs[:n_samples])
        outgroup = SequenceSet(afun_seqs[:1])

        return ingroup, outgroup
    except Exception:
        return None, None


def aggregate_sfs_data(
    fasta_files: list[Path],
    n_samples: int = 20,
) -> tuple[np.ndarray, np.ndarray, int, int, float, float, int]:
    """Aggregate SFS data from multiple gene files.

    Args:
        fasta_files: List of FASTA file paths to process
        n_samples: Number of chromosomes to sample per gene

    Returns:
        Tuple of (sfs_neutral, sfs_selected, dn, ds, n_sites_selected, n_sites_neutral, n_genes)
    """
    n_bins = n_samples // 2  # Folded SFS bins
    sfs_neutral_total = np.zeros(n_bins)
    sfs_selected_total = np.zeros(n_bins)
    dn_total = 0
    ds_total = 0
    n_sites_selected_total = 0.0
    n_sites_neutral_total = 0.0
    n_genes = 0

    for fasta_path in fasta_files:
        ingroup, outgroup = extract_ingroup_outgroup(fasta_path, n_samples)
        if ingroup is None:
            continue

        try:
            poly_data = extract_polymorphism_data(
                ingroup=ingroup,
                outgroup=outgroup,
                reading_frame=1,
            )

            dfe_input = polymorphism_data_to_dfe_input(poly_data, n_samples)

            sfs_neutral_total += dfe_input.sfs_neutral
            sfs_selected_total += dfe_input.sfs_selected
            dn_total += dfe_input.divergence_selected
            ds_total += dfe_input.divergence_neutral

            # Estimate site counts (approximate)
            n_sites_selected_total += np.sum(dfe_input.sfs_selected) * 3
            n_sites_neutral_total += np.sum(dfe_input.sfs_neutral) * 3

            n_genes += 1

        except Exception:
            continue

    return (
        sfs_neutral_total,
        sfs_selected_total,
        dn_total,
        ds_total,
        n_sites_selected_total,
        n_sites_neutral_total,
        n_genes,
    )


def parse_grapes_csv(csv_path: Path, model: str) -> dict:
    """Parse GRAPES CSV output and return results for a specific model.

    GRAPES outputs multiple models per run; this extracts the row for the requested model.

    Args:
        csv_path: Path to GRAPES CSV output file
        model: Model name to extract (e.g., "GammaZero", "GammaExpo")

    Returns:
        Dictionary of column name -> value for the requested model
    """
    if not csv_path.exists():
        return {}

    with open(csv_path) as f:
        lines = f.readlines()

    if len(lines) < 2:
        return {}

    header = lines[0].strip().split(",")
    model_idx = header.index("model") if "model" in header else -1
    if model_idx < 0:
        return {}

    for line in lines[1:]:
        values = line.strip().split(",")
        if len(values) > model_idx and values[model_idx] == model:
            results = {}
            for h, v in zip(header, values):
                try:
                    results[h.strip()] = float(v.strip())
                except ValueError:
                    results[h.strip()] = v.strip()
            return results

    return {}


def run_grapes(dofe_path: Path, model: str) -> tuple[dict, float]:
    """Run GRAPES on a .dofe file and return results with timing.

    Args:
        dofe_path: Path to input .dofe file
        model: DFE model to fit

    Returns:
        Tuple of (results_dict, elapsed_time_seconds)
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        output_path = Path(f.name)

    try:
        start = time.perf_counter()
        result = subprocess.run(
            [
                str(GRAPES_BINARY),
                "-in", str(dofe_path),
                "-out", str(output_path),
                "-model", model,
                "-fold",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        elapsed = time.perf_counter() - start

        if result.returncode != 0:
            print(f"  GRAPES error: {result.stderr.strip()}")
            return {}, elapsed

        return parse_grapes_csv(output_path, model), elapsed

    finally:
        if output_path.exists():
            output_path.unlink()


def run_mkado(dfe_input: DFEInput, model: str) -> tuple[dict, float]:
    """Run mkado DFE fitting and return results with timing.

    Args:
        dfe_input: DFE input data
        model: DFE model to fit

    Returns:
        Tuple of (results_dict, elapsed_time_seconds)
    """
    precalc = PrecomputedData(dfe_input.n_samples)
    dfe_model = get_model(model)

    start = time.perf_counter()
    result = dfe_model.fit(dfe_input, precalc)
    elapsed = time.perf_counter() - start

    parsed = {
        "alpha": result.alpha,
        "omegaA": result.omega_a,
        "lnL": result.log_likelihood,
    }

    if model == "GammaZero":
        parsed["GammaZero:negGshape"] = result.dfe_params.get("shape")
        parsed["GammaZero:negGmean"] = result.dfe_params.get("mean")
    elif model == "GammaExpo":
        parsed["GammaExpo:negGshape"] = result.dfe_params.get("shape_del")
        parsed["GammaExpo:negGmean"] = result.dfe_params.get("mean_del")
        parsed["GammaExpo:pos_prop"] = result.dfe_params.get("prop_ben")
        parsed["GammaExpo:posGmean"] = result.dfe_params.get("mean_ben")

    return parsed, elapsed


def main() -> int:
    """Run the validation comparison."""
    # Check for GRAPES binary
    if not GRAPES_BINARY.exists():
        print("ERROR: GRAPES binary not found")
        print(f"Expected location: {GRAPES_BINARY}")
        print()
        print("To build GRAPES:")
        print("  cd vendor")
        print("  git clone https://github.com/BioPP/grapes.git")
        print("  cd grapes")
        print("  cmake -B build")
        print("  cmake --build build")
        return 1

    # Check for test data
    all_files = get_anopheles_files()
    if not all_files:
        print("ERROR: Anopheles test data not found")
        print(f"Expected location: {ANOPHELES_DATA_DIR}")
        return 1

    print(f"Found {len(all_files)} Anopheles FASTA files")
    print(f"GRAPES binary: {GRAPES_BINARY}")
    print()

    # Test datasets of increasing size
    datasets = [
        ("small", all_files[:50]),
        ("medium", all_files[:150]),
        ("large", all_files[:400]),
    ]

    n_samples = 20
    models = ["GammaZero", "GammaExpo"]
    results_summary = []

    for dataset_name, files in datasets:
        print("=" * 70)
        print(f"Dataset: {dataset_name} ({len(files)} files)")
        print("=" * 70)

        # Aggregate SFS data
        (
            sfs_neutral, sfs_selected, dn, ds,
            n_sites_selected, n_sites_neutral, n_genes
        ) = aggregate_sfs_data(files, n_samples)

        if n_genes < 5:
            print(f"  Skipping: only {n_genes} genes processed successfully")
            continue

        print(f"  Genes processed: {n_genes}")
        print(f"  Polymorphisms: Pn={sfs_selected.sum():.0f}, Ps={sfs_neutral.sum():.0f}")
        print(f"  Divergence: Dn={dn}, Ds={ds}")

        # Create DFE input
        dfe_input = DFEInput(
            sfs_neutral=sfs_neutral,
            sfs_selected=sfs_selected,
            divergence_neutral=ds,
            divergence_selected=dn,
            n_samples=n_samples,
            n_sites_neutral=n_sites_neutral,
            n_sites_selected=n_sites_selected,
        )

        # Export to .dofe for GRAPES
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dofe", delete=False) as f:
            dofe_path = Path(f.name)

        export_to_dofe(
            dfe_input,
            dofe_path,
            dataset_name=f"anopheles_{dataset_name}",
            header=f"Anopheles {dataset_name} dataset ({n_genes} genes)",
            n_sites_selected=n_sites_selected,
            n_sites_neutral=n_sites_neutral,
        )

        for model in models:
            print(f"\n  {model}:")

            grapes_result, grapes_time = run_grapes(dofe_path, model)
            mkado_result, mkado_time = run_mkado(dfe_input, model)

            grapes_alpha = grapes_result.get("alpha", float("nan"))
            mkado_alpha = mkado_result.get("alpha", float("nan"))

            if np.isnan(grapes_alpha) or np.isnan(mkado_alpha):
                alpha_diff = float("nan")
            else:
                alpha_diff = abs(grapes_alpha - mkado_alpha)

            speedup = grapes_time / mkado_time if mkado_time > 0 else 0

            print(f"    GRAPES: alpha={grapes_alpha:.4f} ({grapes_time:.3f}s)")
            print(f"    mkado:  alpha={mkado_alpha:.4f} ({mkado_time:.3f}s)")
            print(f"    Diff: {alpha_diff:.6f}, Speedup: {speedup:.1f}x")

            is_pass = alpha_diff < 0.01
            print(f"    Status: {'PASS' if is_pass else 'FAIL'}")

            results_summary.append({
                "dataset": dataset_name,
                "n_genes": n_genes,
                "model": model,
                "grapes_alpha": grapes_alpha,
                "mkado_alpha": mkado_alpha,
                "alpha_diff": alpha_diff,
                "grapes_time": grapes_time,
                "mkado_time": mkado_time,
                "speedup": speedup,
                "pass": is_pass,
            })

        dofe_path.unlink()
        print()

    # Summary table
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    header = f"{'Dataset':<10} {'Genes':>6} {'Model':<12} {'GRAPES':>10} {'mkado':>10} {'Diff':>10} {'Speedup':>10} {'Status':>8}"
    print(header)
    print("-" * 100)

    all_pass = True
    for r in results_summary:
        status = "PASS" if r["pass"] else "FAIL"
        if not r["pass"]:
            all_pass = False
        print(
            f"{r['dataset']:<10} {r['n_genes']:>6} {r['model']:<12} "
            f"{r['grapes_alpha']:>10.4f} {r['mkado_alpha']:>10.4f} "
            f"{r['alpha_diff']:>10.6f} {r['speedup']:>9.1f}x {status:>8}"
        )

    print("-" * 100)

    if all_pass:
        print("\nAll tests PASSED - mkado and GRAPES produce equivalent results")
        return 0
    else:
        print("\nSome tests FAILED - results differ between mkado and GRAPES")
        return 1


if __name__ == "__main__":
    sys.exit(main())
