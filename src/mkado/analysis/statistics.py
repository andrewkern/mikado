"""Statistical functions for MK test analysis."""

from __future__ import annotations

import math

from scipy import stats


def fishers_exact(
    dn: int, ds: int, pn: int, ps: int, alternative: str = "two-sided"
) -> float:
    """Perform Fisher's exact test on a 2x2 contingency table.

    The table is:
                  | Non-synonymous | Synonymous
        -----------------------------------------------
        Divergence |      dn       |     ds
        Polymorphism |    pn       |     ps

    Args:
        dn: Non-synonymous divergence (fixed differences)
        ds: Synonymous divergence (fixed differences)
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms
        alternative: 'two-sided', 'less', or 'greater'

    Returns:
        p-value from Fisher's exact test
    """
    table = [[dn, ds], [pn, ps]]
    _, p_value = stats.fisher_exact(table, alternative=alternative)
    return float(p_value)


def neutrality_index(dn: int, ds: int, pn: int, ps: int) -> float | None:
    """Calculate the Neutrality Index (NI).

    NI = (Pn/Ps) / (Dn/Ds)

    Under neutrality, NI = 1.
    NI > 1 suggests negative selection (excess polymorphism).
    NI < 1 suggests positive selection (excess divergence).

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        Neutrality Index, or None if cannot be calculated (division by zero)
    """
    if ds == 0 or ps == 0 or dn == 0:
        return None

    return (pn / ps) / (dn / ds)


def alpha(dn: int, ds: int, pn: int, ps: int) -> float | None:
    """Calculate alpha, the proportion of adaptive substitutions.

    α = 1 - (Ds × Pn) / (Dn × Ps)

    Under neutrality, α = 0.
    α > 0 suggests positive selection.
    α < 0 suggests negative selection.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        Alpha value, or None if cannot be calculated
    """
    if dn == 0 or ps == 0:
        return None

    return 1.0 - (ds * pn) / (dn * ps)


def g_test(dn: int, ds: int, pn: int, ps: int) -> float:
    """Perform G-test (log-likelihood ratio test) on MK table.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        p-value from G-test
    """
    observed = [[dn, ds], [pn, ps]]
    total = dn + ds + pn + ps
    if total == 0:
        return 1.0

    row_totals = [dn + ds, pn + ps]
    col_totals = [dn + pn, ds + ps]

    g_stat = 0.0
    for i in range(2):
        for j in range(2):
            obs = observed[i][j]
            if obs > 0:
                expected = row_totals[i] * col_totals[j] / total
                if expected > 0:
                    g_stat += 2 * obs * math.log(obs / expected)

    # G-test statistic follows chi-square with df=1
    p_value = 1 - stats.chi2.cdf(g_stat, df=1)
    return float(p_value)


def confidence_interval_alpha(
    dn: int, ds: int, pn: int, ps: int, confidence: float = 0.95
) -> tuple[float, float] | None:
    """Calculate confidence interval for alpha using bootstrap-like approximation.

    Uses the delta method to approximate the variance of alpha.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms
        confidence: Confidence level (default 0.95)

    Returns:
        Tuple of (lower_bound, upper_bound), or None if cannot be calculated
    """
    a = alpha(dn, ds, pn, ps)
    if a is None:
        return None

    # Use Wilson score interval approximation for small counts
    # For simplicity, use normal approximation with pseudo-counts
    dn_p = max(dn, 0.5)
    ds_p = max(ds, 0.5)
    pn_p = max(pn, 0.5)
    ps_p = max(ps, 0.5)

    # Variance approximation using delta method
    # alpha = 1 - (ds*pn)/(dn*ps)
    # Let r = (ds*pn)/(dn*ps)
    r = (ds_p * pn_p) / (dn_p * ps_p)

    # Var(r) ≈ r^2 * (1/dn + 1/ds + 1/pn + 1/ps) using delta method
    var_r = r**2 * (1 / dn_p + 1 / ds_p + 1 / pn_p + 1 / ps_p)

    # Var(alpha) = Var(r) since alpha = 1 - r
    se_alpha = math.sqrt(var_r)

    # z-score for confidence level
    z = stats.norm.ppf((1 + confidence) / 2)

    lower = a - z * se_alpha
    upper = a + z * se_alpha

    return (lower, upper)
