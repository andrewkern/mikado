"""Imputed McDonald-Kreitman test implementation.

Based on Murga-Moreno et al. (2022) G3: Genes|Genomes|Genetics.
https://doi.org/10.1093/g3journal/jkac206

The impMKT corrects for slightly deleterious mutations by imputing the number
of weakly deleterious nonsynonymous polymorphisms segregating at low frequency,
rather than discarding all low-frequency variants. This retains more data and
increases power to detect positive selection at the gene level.
"""

from __future__ import annotations

from dataclasses import dataclass

from mkado.analysis.asymptotic import PolymorphismData
from mkado.analysis.statistics import fishers_exact


@dataclass
class ImputedMKResult:
    """Results from an imputed McDonald-Kreitman test."""

    alpha: float | None
    p_value: float
    pn_neutral: float
    pwd: float
    dn: int
    ds: int
    pn_total: int
    ps_total: int
    d: float | None = None
    b: float | None = None
    f: float | None = None
    cutoff: float = 0.15

    def __str__(self) -> str:
        alpha_str = f"{self.alpha:.4f}" if self.alpha is not None else "N/A"
        lines = [
            "Imputed MK Test Results:",
            f"  Divergence:    Dn={self.dn}, Ds={self.ds}",
            f"  Polymorphism:  Pn={self.pn_total}, Ps={self.ps_total}",
            f"  DAF cutoff:    {self.cutoff}",
            f"  Imputed Pwd:   {self.pwd:.2f}",
            f"  Pn (neutral):  {self.pn_neutral:.2f}",
            f"  Alpha (α):     {alpha_str}",
            f"  p-value:       {self.p_value:.4g}",
        ]
        if self.d is not None:
            lines.append(f"  DFE fractions: d={self.d:.4f}, b={self.b:.4f}, f={self.f:.4f}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        result = {
            "alpha": self.alpha,
            "p_value": self.p_value,
            "pn_neutral": self.pn_neutral,
            "pwd": self.pwd,
            "dn": self.dn,
            "ds": self.ds,
            "pn_total": self.pn_total,
            "ps_total": self.ps_total,
            "cutoff": self.cutoff,
        }
        if self.d is not None:
            result["d"] = self.d
            result["b"] = self.b
            result["f"] = self.f
        return result


def _compute_imputed(
    polymorphisms: list[tuple[float, str]],
    dn: int,
    ds: int,
    cutoff: float,
    num_synonymous_sites: float | None,
    num_nonsynonymous_sites: float | None,
) -> ImputedMKResult:
    """Core imputed MK algorithm on pooled data."""
    # Split polymorphisms at cutoff
    pn_low = sum(1 for freq, t in polymorphisms if t == "N" and freq <= cutoff)
    pn_high = sum(1 for freq, t in polymorphisms if t == "N" and freq > cutoff)
    ps_low = sum(1 for freq, t in polymorphisms if t == "S" and freq <= cutoff)
    ps_high = sum(1 for freq, t in polymorphisms if t == "S" and freq > cutoff)

    pn_total = pn_low + pn_high
    ps_total = ps_low + ps_high

    # Compute ratio of low to high frequency synonymous polymorphisms
    if ps_high == 0:
        ratio_ps = 0.0
    else:
        ratio_ps = ps_low / ps_high

    # Impute weakly deleterious nonsynonymous polymorphisms
    pwd = pn_low - pn_high * ratio_ps
    pwd = max(pwd, 0.0)

    # Neutral nonsynonymous polymorphisms
    pn_neutral = pn_total - pwd

    # Corrected alpha
    if dn == 0 or ps_total == 0:
        alpha_val = None
    else:
        alpha_val = 1.0 - (pn_neutral / ps_total) * (ds / dn)

    # Fisher's exact test on corrected table
    pn_neutral_int = round(pn_neutral)
    p_value = fishers_exact(dn, ds, pn_neutral_int, ps_total)

    # DFE fractions
    d_frac = None
    b_frac = None
    f_frac = None

    if num_synonymous_sites is not None and num_nonsynonymous_sites is not None:
        m0 = num_synonymous_sites
        mi = num_nonsynonymous_sites
        if mi > 0 and ps_total > 0:
            b_frac = (pwd / ps_total) * (m0 / mi)
            f_frac = (m0 * pn_neutral) / (mi * ps_total)
            d_frac = 1.0 - f_frac - b_frac

    return ImputedMKResult(
        alpha=alpha_val,
        p_value=p_value,
        pn_neutral=pn_neutral,
        pwd=pwd,
        dn=dn,
        ds=ds,
        pn_total=pn_total,
        ps_total=ps_total,
        d=d_frac,
        b=b_frac,
        f=f_frac,
        cutoff=cutoff,
    )


def imputed_mk_test(
    gene_data: PolymorphismData,
    cutoff: float = 0.15,
    num_synonymous_sites: float | None = None,
    num_nonsynonymous_sites: float | None = None,
) -> ImputedMKResult:
    """Perform the imputed McDonald-Kreitman test on a single gene.

    The impMKT corrects for slightly deleterious mutations by imputing the
    count of weakly deleterious nonsynonymous polymorphisms from the ratio
    of low- to high-frequency synonymous polymorphisms.

    Based on Murga-Moreno et al. (2022) G3.
    https://doi.org/10.1093/g3journal/jkac206

    Args:
        gene_data: Polymorphism and divergence data from extract_polymorphism_data()
        cutoff: DAF cutoff separating low and high frequency (default 0.15)
        num_synonymous_sites: Number of synonymous sites (m0) for DFE fractions
        num_nonsynonymous_sites: Number of nonsynonymous sites (mi) for DFE fractions

    Returns:
        ImputedMKResult with corrected alpha and optionally DFE fractions
    """
    return _compute_imputed(
        gene_data.polymorphisms,
        gene_data.dn,
        gene_data.ds,
        cutoff,
        num_synonymous_sites,
        num_nonsynonymous_sites,
    )


def imputed_mk_test_multi(
    gene_data: list[PolymorphismData],
    cutoff: float = 0.15,
    num_synonymous_sites: float | None = None,
    num_nonsynonymous_sites: float | None = None,
) -> ImputedMKResult:
    """Perform the imputed MK test on pooled data from multiple genes.

    Pools all polymorphisms and divergence counts across genes, then
    runs the imputed MK algorithm on the combined data.

    Args:
        gene_data: List of PolymorphismData from individual genes
        cutoff: DAF cutoff separating low and high frequency (default 0.15)
        num_synonymous_sites: Number of synonymous sites (m0) for DFE fractions
        num_nonsynonymous_sites: Number of nonsynonymous sites (mi) for DFE fractions

    Returns:
        ImputedMKResult with corrected alpha from pooled data
    """
    all_polymorphisms: list[tuple[float, str]] = []
    dn_total = 0
    ds_total = 0

    for g in gene_data:
        all_polymorphisms.extend(g.polymorphisms)
        dn_total += g.dn
        ds_total += g.ds

    return _compute_imputed(
        all_polymorphisms,
        dn_total,
        ds_total,
        cutoff,
        num_synonymous_sites,
        num_nonsynonymous_sites,
    )
