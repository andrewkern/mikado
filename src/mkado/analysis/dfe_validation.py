"""Utilities for validating DFE implementation against GRAPES.

This module provides:
1. Export functions to create .dofe files compatible with GRAPES
2. Import functions to read GRAPES output
3. Comparison utilities to validate mkado results against GRAPES
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.analysis.asymptotic import PolymorphismData
from mkado.analysis.dfe_models import DFEInput, DFEResult

if TYPE_CHECKING:
    pass


def export_to_dofe(
    dfe_input: DFEInput,
    output_path: Path | str,
    dataset_name: str = "mkado_export",
    header: str = "Exported from mkado",
    folded: bool = True,
    n_sites_selected: float | None = None,
    n_sites_neutral: float | None = None,
) -> None:
    """Export DFE input data to GRAPES .dofe format.

    The .dofe format is:
    - Line 1: Header/comment
    - Line 2 (optional): #unfolded (if using unfolded SFS)
    - Data line: name n_samples n_sites_N SFS_N... n_sites_S SFS_S... n_sites_div_N Dn n_sites_div_S Ds

    Args:
        dfe_input: DFE input data from mkado
        output_path: Path to output .dofe file
        dataset_name: Name for the dataset in the output
        header: Header comment line
        folded: Whether the SFS is folded (True) or unfolded (False)
        n_sites_selected: Number of nonsynonymous sites (estimated if None)
        n_sites_neutral: Number of synonymous sites (estimated if None)
    """
    output_path = Path(output_path)

    lines = [header]

    if not folded:
        lines.append("#unfolded")

    n = dfe_input.n_samples

    # Calculate approximate number of sites if not provided
    if n_sites_selected is None:
        n_sites_selected = float(np.sum(dfe_input.sfs_selected)) * 10
    if n_sites_neutral is None:
        n_sites_neutral = float(np.sum(dfe_input.sfs_neutral)) * 10

    # Build SFS strings (tab-separated)
    sfs_n_parts = [f"{x:.6f}" for x in dfe_input.sfs_selected]
    sfs_s_parts = [f"{x:.6f}" for x in dfe_input.sfs_neutral]

    # Build data line with all fields tab-separated
    # Format: name n n_sites_N SFS_N[0] SFS_N[1]... n_sites_S SFS_S[0] SFS_S[1]... n_sites_div_N Dn n_sites_div_S Ds
    parts = [
        dataset_name,
        str(n),
        f"{n_sites_selected:.1f}",
    ]
    parts.extend(sfs_n_parts)
    parts.append(f"{n_sites_neutral:.1f}")
    parts.extend(sfs_s_parts)
    parts.extend([
        f"{n_sites_selected:.1f}",  # n_sites_div_N (same as poly sites for simplicity)
        str(dfe_input.divergence_selected),
        f"{n_sites_neutral:.1f}",   # n_sites_div_S
        str(dfe_input.divergence_neutral),
    ])

    data_line = "\t".join(parts)
    lines.append(data_line)

    output_path.write_text("\n".join(lines) + "\n")


def export_polymorphism_data_to_dofe(
    poly_data: PolymorphismData,
    n_samples: int,
    output_path: Path | str,
    dataset_name: str = "mkado_export",
    header: str = "Exported from mkado",
) -> None:
    """Export PolymorphismData to GRAPES .dofe format.

    Args:
        poly_data: Polymorphism data from mkado
        n_samples: Number of chromosomes sampled
        output_path: Path to output .dofe file
        dataset_name: Name for the dataset
        header: Header comment line
    """
    from mkado.analysis.dfe import polymorphism_data_to_dfe_input

    dfe_input = polymorphism_data_to_dfe_input(poly_data, n_samples)
    export_to_dofe(dfe_input, output_path, dataset_name, header)


def parse_grapes_output(csv_path: Path | str) -> dict:
    """Parse GRAPES CSV output file.

    Args:
        csv_path: Path to GRAPES output CSV

    Returns:
        Dictionary with parsed results
    """
    csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(f"GRAPES output not found: {csv_path}")

    results = {}

    with open(csv_path) as f:
        lines = f.readlines()

    # GRAPES CSV has header line followed by data
    if len(lines) < 2:
        return results

    header = lines[0].strip().split(",")
    values = lines[1].strip().split(",")

    for h, v in zip(header, values):
        try:
            results[h.strip()] = float(v.strip())
        except ValueError:
            results[h.strip()] = v.strip()

    return results


def compare_results(
    mkado_result: DFEResult,
    grapes_result: dict,
    tolerance: float = 1e-6,
) -> dict:
    """Compare mkado DFE results against GRAPES output.

    Args:
        mkado_result: Result from mkado dfe_alpha()
        grapes_result: Parsed GRAPES output from parse_grapes_output()
        tolerance: Numerical tolerance for comparison

    Returns:
        Dictionary with comparison results:
        - 'matches': bool - True if all values match within tolerance
        - 'differences': dict - Differences for each compared value
    """
    differences = {}

    # Map mkado fields to GRAPES fields
    field_mapping = {
        "alpha": "alpha",
        "omega_a": "omegaA",
        "omega_na": "omegaNA",
    }

    # Compare DFE parameters based on model
    if mkado_result.model == "GammaZero":
        field_mapping["dfe_params.shape"] = "negGshape"
        field_mapping["dfe_params.mean"] = "negGmean"
    elif mkado_result.model == "GammaExpo":
        field_mapping["dfe_params.shape_del"] = "negGshape"
        field_mapping["dfe_params.mean_del"] = "negGmean"
        field_mapping["dfe_params.prop_ben"] = "pos_prop"
        field_mapping["dfe_params.mean_ben"] = "posGmean"

    all_match = True

    for mkado_field, grapes_field in field_mapping.items():
        if grapes_field not in grapes_result:
            continue

        # Get mkado value
        if "." in mkado_field:
            parts = mkado_field.split(".")
            mkado_value = mkado_result.dfe_params.get(parts[1])
        else:
            mkado_value = getattr(mkado_result, mkado_field, None)

        if mkado_value is None:
            continue

        grapes_value = grapes_result[grapes_field]

        if isinstance(grapes_value, (int, float)):
            diff = abs(mkado_value - grapes_value)
            rel_diff = diff / max(abs(grapes_value), 1e-10)

            if diff > tolerance and rel_diff > tolerance:
                all_match = False
                differences[mkado_field] = {
                    "mkado": mkado_value,
                    "grapes": grapes_value,
                    "abs_diff": diff,
                    "rel_diff": rel_diff,
                }

    return {
        "matches": all_match,
        "differences": differences,
    }


# Reference values from GRAPES for validation
# These are computed from the Microtus example data using GRAPES built from source
# Command: grapes -in Microtus.folded.dofe -out output.csv -model <Model> -fold
GRAPES_REFERENCE_VALUES = {
    "Microtus.folded": {
        "GammaZero": {
            "negGshape": 0.316874,
            "negGmean": 19895.262355,
            "alpha": 0.6495,
            "omegaA": 0.0637,
            "logL": -62.8716,
        },
        "GammaExpo": {
            "negGshape": 0.377815,
            "negGmean": 10176.942064,
            "pos_prop": 0.008380,
            "posGmean": 9999.999756,
            "alpha": 0.7416,
            "omegaA": 0.0727,
            "logL": -62.8710,
        },
    },
}
