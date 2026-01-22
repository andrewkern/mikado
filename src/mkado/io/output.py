"""Output formatters for MK test results."""

from __future__ import annotations

import json
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult


class OutputFormat(Enum):
    """Supported output formats."""

    PRETTY = "pretty"
    TSV = "tsv"
    JSON = "json"


def format_result(
    result: MKResult | PolarizedMKResult | AsymptoticMKResult,
    format: OutputFormat = OutputFormat.PRETTY,
) -> str:
    """Format MK test results for output.

    Args:
        result: MK test result object
        format: Output format

    Returns:
        Formatted string representation
    """
    if format == OutputFormat.PRETTY:
        return str(result)

    elif format == OutputFormat.JSON:
        return json.dumps(result.to_dict(), indent=2)

    elif format == OutputFormat.TSV:
        return _format_tsv(result)

    else:
        raise ValueError(f"Unknown format: {format}")


def _format_tsv(result: MKResult | PolarizedMKResult | AsymptoticMKResult) -> str:
    """Format results as tab-separated values."""
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    if isinstance(result, MKResult):
        header = "Dn\tDs\tPn\tPs\tp_value\tNI\talpha"
        ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
        alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
        values = f"{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t{result.p_value:.6g}\t{ni_str}\t{alpha_str}"
        return f"{header}\n{values}"

    elif isinstance(result, PolarizedMKResult):
        header = "lineage\tDn\tDs\tPn\tPs\tp_value\tNI\talpha"
        lines = [header]

        ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
        alpha_str = f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
        lines.append(
            f"ingroup\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
            f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
            f"{result.p_value_ingroup:.6g}\t{ni_str}\t{alpha_str}"
        )
        lines.append(f"outgroup\t{result.dn_outgroup}\t{result.ds_outgroup}\tNA\tNA\tNA\tNA\tNA")
        lines.append(
            f"unpolarized\t{result.dn_unpolarized}\t{result.ds_unpolarized}\tNA\tNA\tNA\tNA\tNA"
        )
        return "\n".join(lines)

    elif isinstance(result, AsymptoticMKResult):
        if result.num_genes > 0:
            # Aggregated result
            header = "Dn\tDs\tPn\tPs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tnum_genes"
            values = (
                f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
                f"{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
                f"{result.model_type}\t{result.num_genes}"
            )
        else:
            header = "Dn\tDs\talpha_asymptotic\tCI_low\tCI_high"
            values = (
                f"{result.dn}\t{result.ds}\t{result.alpha_asymptotic:.6f}\t"
                f"{result.ci_low:.6f}\t{result.ci_high:.6f}"
            )
        return f"{header}\n{values}"

    else:
        raise TypeError(f"Unknown result type: {type(result)}")


def format_batch_results(
    results: list[tuple[str, MKResult | PolarizedMKResult | AsymptoticMKResult]],
    format: OutputFormat = OutputFormat.PRETTY,
) -> str:
    """Format multiple MK test results for batch output.

    Args:
        results: List of (name, result) tuples
        format: Output format

    Returns:
        Formatted string representation
    """
    if format == OutputFormat.PRETTY:
        lines = []
        for name, result in results:
            lines.append(f"=== {name} ===")
            lines.append(str(result))
            lines.append("")
        return "\n".join(lines)

    elif format == OutputFormat.JSON:
        data = {name: result.to_dict() for name, result in results}
        return json.dumps(data, indent=2)

    elif format == OutputFormat.TSV:
        from mkado.analysis.asymptotic import AsymptoticMKResult
        from mkado.analysis.mk_test import MKResult
        from mkado.analysis.polarized import PolarizedMKResult

        if not results:
            return ""

        # Check type of first result
        _, first_result = results[0]
        if isinstance(first_result, MKResult):
            header = "gene\tDn\tDs\tPn\tPs\tp_value\tNI\talpha"
            lines = [header]
            for name, result in results:
                if isinstance(result, MKResult):
                    ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
                    alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
                    lines.append(
                        f"{name}\t{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t"
                        f"{result.p_value:.6g}\t{ni_str}\t{alpha_str}"
                    )
            return "\n".join(lines)

        elif isinstance(first_result, AsymptoticMKResult):
            header = "gene\tDn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel"
            lines = [header]
            for name, result in results:
                if isinstance(result, AsymptoticMKResult):
                    lines.append(
                        f"{name}\t{result.dn}\t{result.ds}\t"
                        f"{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t"
                        f"{result.ci_high:.6f}\t{result.model_type}"
                    )
            return "\n".join(lines)

        elif isinstance(first_result, PolarizedMKResult):
            header = "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value\tNI\talpha"
            lines = [header]
            for name, result in results:
                if isinstance(result, PolarizedMKResult):
                    ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
                    alpha_str = f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
                    lines.append(
                        f"{name}\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
                        f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
                        f"{result.dn_outgroup}\t{result.ds_outgroup}\t"
                        f"{result.p_value_ingroup:.6g}\t{ni_str}\t{alpha_str}"
                    )
            return "\n".join(lines)

        else:
            # Fall back to JSON for unknown types
            data = {name: result.to_dict() for name, result in results}
            return json.dumps(data, indent=2)

    else:
        raise ValueError(f"Unknown format: {format}")
