"""Output formatters for MK test results."""

from __future__ import annotations

import json
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult


class OutputFormat(Enum):
    """Supported output formats."""

    PRETTY = "pretty"
    TSV = "tsv"
    JSON = "json"


def format_result(
    result: MKResult | PolarizedMKResult | AsymptoticMKResult | AlphaTGResult | ImputedMKResult,
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


def _format_tsv(result: MKResult | PolarizedMKResult | AsymptoticMKResult | AlphaTGResult | ImputedMKResult) -> str:
    """Format results as tab-separated values."""
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    if isinstance(result, ImputedMKResult):
        alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
        header = "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff"
        values = (
            f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
            f"{result.pwd:.2f}\t{result.pn_neutral:.2f}\t{alpha_str}\t"
            f"{result.p_value:.6g}\t{result.cutoff}"
        )
        return f"{header}\n{values}"

    elif isinstance(result, MKResult):
        header = "Dn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS"
        ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
        alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
        dos_str = f"{result.dos:.6f}" if result.dos is not None else "NA"
        values = f"{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t{result.p_value:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
        return f"{header}\n{values}"

    elif isinstance(result, PolarizedMKResult):
        header = "lineage\tDn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS"
        lines = [header]

        ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
        alpha_str = f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
        dos_str = f"{result.dos_ingroup:.6f}" if result.dos_ingroup is not None else "NA"
        lines.append(
            f"ingroup\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
            f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
            f"{result.p_value_ingroup:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
        )
        lines.append(f"outgroup\t{result.dn_outgroup}\t{result.ds_outgroup}\tNA\tNA\tNA\tNA\tNA\tNA")
        lines.append(
            f"unpolarized\t{result.dn_unpolarized}\t{result.ds_unpolarized}\tNA\tNA\tNA\tNA\tNA\tNA"
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

    elif isinstance(result, AlphaTGResult):
        header = "Dn\tDs\tPn\tPs\talpha_TG\tNI_TG\tCI_low\tCI_high\tnum_genes"
        values = (
            f"{result.dn_total}\t{result.ds_total}\t{result.pn_total}\t{result.ps_total}\t"
            f"{result.alpha_tg:.6f}\t{result.ni_tg:.6f}\t{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
            f"{result.num_genes}"
        )
        return f"{header}\n{values}"

    else:
        raise TypeError(f"Unknown result type: {type(result)}")


def format_batch_results(
    results: list[tuple[str, MKResult | PolarizedMKResult | AsymptoticMKResult]],
    format: OutputFormat = OutputFormat.PRETTY,
    adjusted_pvalues: list[float] | None = None,
) -> str:
    """Format multiple MK test results for batch output.

    Args:
        results: List of (name, result) tuples
        format: Output format
        adjusted_pvalues: Optional list of Benjamini-Hochberg adjusted p-values
            (must be same length as results if provided)

    Returns:
        Formatted string representation
    """
    if format == OutputFormat.PRETTY:
        lines = []
        for i, (name, result) in enumerate(results):
            lines.append(f"=== {name} ===")
            lines.append(str(result))
            if adjusted_pvalues is not None:
                lines.append(f"  p-value (BH adj):     {adjusted_pvalues[i]:.6g}")
            lines.append("")
        return "\n".join(lines)

    elif format == OutputFormat.JSON:
        data = {}
        for i, (name, result) in enumerate(results):
            result_dict = result.to_dict()
            if adjusted_pvalues is not None:
                result_dict["p_value_adjusted"] = adjusted_pvalues[i]
            data[name] = result_dict
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
            if adjusted_pvalues is not None:
                header = "gene\tDn\tDs\tPn\tPs\tp_value\tp_value_adjusted\tNI\talpha\tDoS"
            else:
                header = "gene\tDn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS"
            lines = [header]
            for i, (name, result) in enumerate(results):
                if isinstance(result, MKResult):
                    ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
                    alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
                    dos_str = f"{result.dos:.6f}" if result.dos is not None else "NA"
                    if adjusted_pvalues is not None:
                        lines.append(
                            f"{name}\t{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t"
                            f"{result.p_value:.6g}\t{adjusted_pvalues[i]:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
                        )
                    else:
                        lines.append(
                            f"{name}\t{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t"
                            f"{result.p_value:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
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
            if adjusted_pvalues is not None:
                header = "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value\tp_value_adjusted\tNI\talpha\tDoS"
            else:
                header = "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value\tNI\talpha\tDoS"
            lines = [header]
            for i, (name, result) in enumerate(results):
                if isinstance(result, PolarizedMKResult):
                    ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
                    alpha_str = f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
                    dos_str = f"{result.dos_ingroup:.6f}" if result.dos_ingroup is not None else "NA"
                    if adjusted_pvalues is not None:
                        lines.append(
                            f"{name}\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
                            f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
                            f"{result.dn_outgroup}\t{result.ds_outgroup}\t"
                            f"{result.p_value_ingroup:.6g}\t{adjusted_pvalues[i]:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
                        )
                    else:
                        lines.append(
                            f"{name}\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
                            f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
                            f"{result.dn_outgroup}\t{result.ds_outgroup}\t"
                            f"{result.p_value_ingroup:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}"
                        )
            return "\n".join(lines)

        else:
            # Fall back to JSON for unknown types
            data = {name: result.to_dict() for name, result in results}
            return json.dumps(data, indent=2)

    else:
        raise ValueError(f"Unknown format: {format}")
