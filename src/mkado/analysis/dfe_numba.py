"""Numba-compiled kernels for DFE-based alpha estimation.

Based on GRAPES implementation (Galtier 2016; Al-Saffar & Hahn 2022).
Reference: https://github.com/BioPP/grapes
"""

from __future__ import annotations

import math

import numba
import numpy as np


# =============================================================================
# Pre-computation grid setup
# =============================================================================


def precalculate_s_grid(n_points: int = 4000) -> np.ndarray:
    """Create S grid matching GRAPES structure.

    GRAPES uses a specific grid structure with different spacing in different regions:
    - Large negative: -500 to -0.5, 1000 points
    - Medium: -0.5 to +0.5, 1000 points
    - Large positive: +0.5 to +500, 1000 points
    - Very large positive: +500 to +10000, 1000 points

    Args:
        n_points: Number of grid points (default 4000, ignored - uses GRAPES structure)

    Returns:
        Array of S values (4*Ne*s) for numerical integration
    """
    # GRAPES constants
    LARGE_S_VALUE = 500
    VERY_LARGE_S_VALUE = 10000
    GAMMA_INTEGRAL_BREAK = 0.5
    k1 = 1000  # INTEGRAL_NB_RECTANGLES_EXTREME
    k2 = 1000  # INTEGRAL_NB_RECTANGLES_MEDIAN

    s_points = []

    # Large negative: from -LARGE_S_VALUE to -GAMMA_INTEGRAL_BREAK
    a1 = (LARGE_S_VALUE - GAMMA_INTEGRAL_BREAK) / k1
    for i in range(k1):
        s_points.append(-LARGE_S_VALUE + (i + 0.5) * a1)

    # Medium values: from -GAMMA_INTEGRAL_BREAK to GAMMA_INTEGRAL_BREAK
    a2 = 2.0 * GAMMA_INTEGRAL_BREAK / k2
    for i in range(k2):
        s_points.append(-GAMMA_INTEGRAL_BREAK + (i + 0.5) * a2)

    # Large positive: from GAMMA_INTEGRAL_BREAK to LARGE_S_VALUE
    for i in range(k1):
        s_points.append(GAMMA_INTEGRAL_BREAK + (i + 0.5) * a1)

    # Very large positive: from LARGE_S_VALUE to VERY_LARGE_S_VALUE
    a3 = (VERY_LARGE_S_VALUE - LARGE_S_VALUE) / k1
    for i in range(k1):
        s_points.append(LARGE_S_VALUE + (i + 0.5) * a3)

    return np.array(s_points)


def precalculate_s_widths(s_grid: np.ndarray) -> np.ndarray:
    """Compute integration widths for each S grid point.

    Following GRAPES, computes the width of each rectangle for numerical integration.
    For interior points, uses the average of distances to neighbors.
    For endpoints, uses the distance to the single neighbor.

    Args:
        s_grid: Array of S values from precalculate_s_grid()

    Returns:
        Array of widths for integration
    """
    n = len(s_grid)
    widths = np.zeros(n)

    for i in range(n):
        if i == 0:
            # First point: width is distance to next point
            widths[i] = s_grid[i + 1] - s_grid[i]
        elif i == n - 1:
            # Last point: width is distance to previous point
            widths[i] = s_grid[i] - s_grid[i - 1]
        else:
            # Interior point: average of distances to neighbors
            widths[i] = (s_grid[i + 1] - s_grid[i - 1]) / 2.0

    return widths


# =============================================================================
# Sojourn time and expected allele count computation
# =============================================================================


@numba.njit(cache=True, fastmath=True)
def sojourn_density(x: float, S: float) -> float:
    """Compute sojourn time density at frequency x under selection S.

    The sojourn time density gives the expected time a mutation spends
    at frequency x before fixation or loss.

    For S = 4*Ne*s:
        density(x) = 2*(1 - exp(-S*(1-x))) / (x*(1-x)*(1-exp(-S)))

    For S=0 (neutral):
        density(x) = 2 / (x*(1-x))

    Args:
        x: Allele frequency (0 < x < 1)
        S: Scaled selection coefficient (4*Ne*s)

    Returns:
        Sojourn time density at frequency x
    """
    if x <= 0.0 or x >= 1.0:
        return 0.0

    base = 2.0 / (x * (1.0 - x))

    if abs(S) < 1e-10:
        # Neutral case
        return base

    # For strong negative selection (S << 0), density -> 0
    # Use asymptotic form: density ≈ 2 * exp(-|S|*x) / (x*(1-x))
    if S < -700:
        # exp(-|S|*x) underflows to 0 for large |S|
        return 0.0

    # For strong positive selection (S >> 0)
    # density ≈ 2*S / (x*(1-x)*(exp(S)-1)) ≈ 2*S*exp(-S) / (x*(1-x))
    if S > 700:
        # When S is very large, numerator ≈ 1, denom ≈ 1
        # Actually for large positive S, (1-exp(-S)) ≈ 1
        # and (1-exp(-S*(1-x))) ≈ 1 for most x
        return base

    # Standard formula with numerical checks
    exp_neg_S = math.exp(-S)
    exp_neg_S_1mx = math.exp(-S * (1.0 - x))

    numer = 1.0 - exp_neg_S_1mx
    denom_factor = 1.0 - exp_neg_S

    if abs(denom_factor) < 1e-300:
        return 0.0

    return base * numer / denom_factor


@numba.njit(cache=True, fastmath=True)
def binomial_coeff(n: int, k: int) -> float:
    """Compute binomial coefficient C(n,k) using log-gamma for stability."""
    if k < 0 or k > n:
        return 0.0
    if k == 0 or k == n:
        return 1.0

    # Use log for numerical stability
    result = 0.0
    for i in range(k):
        result += math.log(n - i) - math.log(i + 1)

    return math.exp(result)


@numba.njit(cache=True, fastmath=True)
def binomial_prob(n: int, k: int, x: float) -> float:
    """Compute binomial probability P(X=k | n, x).

    Args:
        n: Number of trials (sample size)
        k: Number of successes (allele count)
        x: Probability (allele frequency)

    Returns:
        Binomial probability
    """
    if x <= 0.0:
        return 1.0 if k == 0 else 0.0
    if x >= 1.0:
        return 1.0 if k == n else 0.0

    coeff = binomial_coeff(n, k)
    return coeff * (x ** k) * ((1.0 - x) ** (n - k))


@numba.njit(cache=True, fastmath=True)
def _integrand(x: float, S: float, j: int, n: int) -> float:
    """Integrand for sojourn integral: sojourn_density(x, S) * binomial_prob(n, j, x)."""
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return sojourn_density(x, S) * binomial_prob(n, j, x)


@numba.njit(cache=True, fastmath=True)
def _sojourn_integral_romberg(S: float, j: int, n: int, max_iter: int = 15) -> float:
    """Compute integral of sojourn density * binomial sampling.

    Integrates: int_0^1 sojourn_density(x, S) * binom(n, j, x) dx

    For strong negative selection (S < -50), uses a log transformation to
    properly sample the sharply peaked region near x=0. For other cases,
    uses standard Romberg integration.

    Args:
        S: Scaled selection coefficient (4*Ne*s)
        j: Allele count in sample
        n: Sample size
        max_iter: Maximum Romberg iterations

    Returns:
        Integral value (expected contribution to SFS bin j)
    """
    # For strong negative selection, the integrand is sharply peaked near x=0.
    # Use a log transformation: t = log(x), so x = exp(t), dx = exp(t) dt
    # This samples more densely near x=0.
    if S < -50:
        return _sojourn_integral_log_transform(S, j, n, max_iter)

    # Standard Romberg integration for moderate selection
    return _sojourn_integral_romberg_standard(S, j, n, max_iter)


@numba.njit(cache=True, fastmath=True)
def _sojourn_integral_log_transform(S: float, j: int, n: int, max_iter: int = 15) -> float:
    """Compute sojourn integral using log transformation for strong negative S.

    Uses substitution t = log(x), so x = exp(t), dx = exp(t) dt.
    Integral becomes: int_{log(x_min)}^{log(x_max)} f(exp(t)) * exp(t) dt

    This transformation naturally puts more sample points near x=0.
    """
    # Integration bounds in x-space (following GRAPES QuasiZero)
    x_min = 0.0001
    x_max = 0.9999

    # Transform to t-space: t = log(x)
    t_min = math.log(x_min)  # ~ -9.2
    t_max = math.log(x_max)  # ~ -0.0001

    # Initialize Romberg tableau
    R = np.zeros((max_iter, max_iter))

    # Evaluate transformed integrand: f(exp(t)) * exp(t)
    def eval_at_t(t: float) -> float:
        x = math.exp(t)
        return _integrand(x, S, j, n) * x

    # First trapezoidal estimate
    h = t_max - t_min
    R[0, 0] = 0.5 * h * (eval_at_t(t_min) + eval_at_t(t_max))

    # Refine with more trapezoids
    for i in range(1, max_iter):
        h = h / 2.0
        n_new = 2 ** (i - 1)

        # Sum function values at new midpoints
        total = 0.0
        for k in range(n_new):
            t = t_min + h * (2 * k + 1)
            total += eval_at_t(t)

        R[i, 0] = 0.5 * R[i - 1, 0] + h * total

        # Richardson extrapolation
        for k in range(1, i + 1):
            factor = 4.0 ** k
            R[i, k] = (factor * R[i, k - 1] - R[i - 1, k - 1]) / (factor - 1.0)

        # Check convergence
        if i >= 2:
            if abs(R[i, i] - R[i - 1, i - 1]) < 1e-10 * abs(R[i, i]) + 1e-15:
                return R[i, i]

    return R[max_iter - 1, max_iter - 1]


@numba.njit(cache=True, fastmath=True)
def _sojourn_integral_romberg_standard(S: float, j: int, n: int, max_iter: int = 15) -> float:
    """Standard Romberg integration over [x_min, x_max].

    Following GRAPES, uses QuasiZero=0.0001 as bounds.
    """
    # Integration bounds (following GRAPES QuasiZero)
    x_min = 0.0001
    x_max = 0.9999

    # Initialize Romberg tableau
    R = np.zeros((max_iter, max_iter))

    # First trapezoidal estimate
    h = x_max - x_min
    R[0, 0] = 0.5 * h * (_integrand(x_min, S, j, n) + _integrand(x_max, S, j, n))

    # Refine with more trapezoids
    for i in range(1, max_iter):
        h = h / 2.0
        n_new = 2 ** (i - 1)

        # Sum function values at new midpoints
        total = 0.0
        for k in range(n_new):
            x = x_min + h * (2 * k + 1)
            total += _integrand(x, S, j, n)

        R[i, 0] = 0.5 * R[i - 1, 0] + h * total

        # Richardson extrapolation
        for k in range(1, i + 1):
            factor = 4.0 ** k
            R[i, k] = (factor * R[i, k - 1] - R[i - 1, k - 1]) / (factor - 1.0)

        # Check convergence
        if i >= 2:
            if abs(R[i, i] - R[i - 1, i - 1]) < 1e-10 * abs(R[i, i]) + 1e-15:
                return R[i, i]

    return R[max_iter - 1, max_iter - 1]


@numba.njit(cache=True, fastmath=True)
def expected_count_unfolded(S: float, j: int, n: int) -> float:
    """Compute expected unfolded SFS count at frequency j/n under selection S.

    Args:
        S: Scaled selection coefficient (4*Ne*s)
        j: Allele count (1 to n-1)
        n: Sample size

    Returns:
        Expected number of segregating sites at frequency j/n
    """
    if j <= 0 or j >= n:
        return 0.0

    return _sojourn_integral_romberg(S, j, n)


@numba.njit(cache=True, parallel=True)
def precalculate_expected_counts(
    s_grid: np.ndarray, n: int
) -> np.ndarray:
    """Pre-compute expected allele counts E[j | S] for all S values.

    Pre-computes the integral of sojourn time * binomial sampling for
    each combination of S value and allele count j.

    Args:
        s_grid: Array of S values from precalculate_s_grid()
        n: Sample size (number of chromosomes)

    Returns:
        Array of shape (len(s_grid), n-1) where [i, j-1] = E[count_j | S_i]
        Note: j is 1-indexed in biology, 0-indexed in array
    """
    n_s = len(s_grid)
    # Stores expected counts for j = 1, 2, ..., n-1
    counts = np.zeros((n_s, n - 1))

    for i in numba.prange(n_s):
        S = s_grid[i]
        for j in range(1, n):
            counts[i, j - 1] = expected_count_unfolded(S, j, n)

    return counts


def fold_sfs(unfolded: np.ndarray) -> np.ndarray:
    """Fold an SFS by combining complementary frequency classes.

    Args:
        unfolded: Unfolded SFS of length n-1 (counts for j=1 to n-1)

    Returns:
        Folded SFS of length floor(n/2)
    """
    n = len(unfolded) + 1
    folded_len = n // 2

    folded = np.zeros(folded_len)

    for j in range(1, folded_len + 1):
        complement = n - j
        if j == complement:
            # At midpoint (n even), count once
            folded[j - 1] = unfolded[j - 1]
        else:
            # Combine j and n-j
            folded[j - 1] = unfolded[j - 1] + unfolded[complement - 1]

    return folded


# =============================================================================
# Fixation probability computation
# =============================================================================


@numba.njit(cache=True, fastmath=True)
def fixation_probability(S: float) -> float:
    """Compute fixation probability for a mutation with selection S.

    For a new mutation with selection coefficient s:
        P_fix = (1 - exp(-S/N)) / (1 - exp(-S))

    where S = 4*Ne*s. For large N this simplifies to:
        P_fix = S / (1 - exp(-S))

    For S=0 (neutral): P_fix = 1/N (approximated as very small)

    Args:
        S: Scaled selection coefficient (4*Ne*s)

    Returns:
        Fixation probability (relative to neutral)
    """
    if abs(S) < 1e-10:
        # Neutral case: P_fix = 1 (relative)
        return 1.0

    # For numerical stability with large negative S
    if S < -500:
        return 0.0

    if S > 500:
        return S

    exp_neg_S = math.exp(-S)
    return S / (1.0 - exp_neg_S)


@numba.njit(cache=True, fastmath=True)
def fixation_correction(S: float, n: int) -> float:
    """Compute correction for apparent divergence from unsampled polymorphisms.

    When only one allele class is sampled from a polymorphic site,
    it appears as a fixed difference. This correction accounts for that.

    Args:
        S: Scaled selection coefficient
        n: Sample size

    Returns:
        Correction factor
    """
    if abs(S) < 1e-10:
        # Neutral case: 1/n
        return 1.0 / n

    # Numerical integration
    result = 0.0
    n_steps = 100
    dx = 1.0 / n_steps

    for i in range(1, n_steps):
        x = i * dx
        # Probability of sampling only derived alleles (x^n)
        # weighted by sojourn density
        density = sojourn_density(x, S)
        result += density * (x ** n) * dx

    return result


# =============================================================================
# DFE density functions
# =============================================================================


@numba.njit(cache=True, fastmath=True)
def gamma_density(s: float, shape: float, mean: float) -> float:
    """Gamma distribution density for deleterious mutations.

    The DFE is modeled as a gamma distribution for deleterious effects.
    Convention: s < 0 for deleterious, uses |s| for computation.

    Args:
        s: Selection coefficient (should be negative for deleterious)
        shape: Shape parameter (alpha) of gamma distribution
        mean: Mean of gamma distribution (positive value, mean |s|)

    Returns:
        Probability density at s
    """
    if s >= 0:
        return 0.0

    # Use absolute value
    abs_s = -s

    if abs_s <= 0 or mean <= 0 or shape <= 0:
        return 0.0

    # Gamma density: x^(α-1) * exp(-x/θ) / (θ^α * Γ(α))
    # where θ = mean/shape (scale parameter)
    scale = mean / shape

    # Use log for numerical stability
    log_density = ((shape - 1.0) * math.log(abs_s)
                   - abs_s / scale
                   - shape * math.log(scale)
                   - math.lgamma(shape))

    return math.exp(log_density)


@numba.njit(cache=True, fastmath=True)
def exponential_density(s: float, mean: float) -> float:
    """Exponential distribution density for beneficial mutations.

    Args:
        s: Selection coefficient (should be positive for beneficial)
        mean: Mean of exponential distribution

    Returns:
        Probability density at s
    """
    if s <= 0 or mean <= 0:
        return 0.0

    return math.exp(-s / mean) / mean


@numba.njit(cache=True, fastmath=True)
def gamma_expo_density(
    s: float,
    shape_del: float,
    mean_del: float,
    prop_ben: float,
    mean_ben: float,
) -> float:
    """Gamma (deleterious) + Exponential (beneficial) DFE.

    Models the DFE as:
    - (1 - prop_ben) * Gamma(|s|; shape_del, mean_del) for s < 0
    - prop_ben * Exponential(s; mean_ben) for s > 0

    This is the best-performing model according to Al-Saffar & Hahn (2022).

    Args:
        s: Selection coefficient
        shape_del: Shape parameter for deleterious gamma
        mean_del: Mean for deleterious gamma
        prop_ben: Proportion of beneficial mutations
        mean_ben: Mean for beneficial exponential

    Returns:
        DFE density at s
    """
    if s < 0:
        return (1.0 - prop_ben) * gamma_density(s, shape_del, mean_del)
    elif s > 0:
        return prop_ben * exponential_density(s, mean_ben)
    else:
        return 0.0


@numba.njit(cache=True, fastmath=True)
def gamma_gamma_density(
    s: float,
    shape_del: float,
    mean_del: float,
    prop_ben: float,
    shape_ben: float,
    mean_ben: float,
) -> float:
    """Gamma (deleterious) + Gamma (beneficial) DFE.

    Models the DFE as:
    - (1 - prop_ben) * Gamma(|s|; shape_del, mean_del) for s < 0
    - prop_ben * Gamma(s; shape_ben, mean_ben) for s > 0

    Args:
        s: Selection coefficient
        shape_del: Shape parameter for deleterious gamma
        mean_del: Mean for deleterious gamma
        prop_ben: Proportion of beneficial mutations
        shape_ben: Shape parameter for beneficial gamma
        mean_ben: Mean for beneficial gamma

    Returns:
        DFE density at s
    """
    if s < 0:
        return (1.0 - prop_ben) * gamma_density(s, shape_del, mean_del)
    elif s > 0:
        # Use gamma density for positive s (flip sign convention)
        abs_s = s
        if abs_s <= 0 or mean_ben <= 0 or shape_ben <= 0:
            return 0.0

        scale = mean_ben / shape_ben
        log_density = ((shape_ben - 1.0) * math.log(abs_s)
                       - abs_s / scale
                       - shape_ben * math.log(scale)
                       - math.lgamma(shape_ben))
        return prop_ben * math.exp(log_density)
    else:
        return 0.0


@numba.njit(cache=True, fastmath=True)
def displaced_gamma_density(
    s: float,
    shape: float,
    mean: float,
    displacement: float,
) -> float:
    """Displaced Gamma DFE.

    Shifts the gamma distribution by a displacement factor to allow
    some beneficial mutations.

    The density is Gamma(s - displacement) where displacement is typically
    negative, allowing s > 0 to have non-zero probability.

    Args:
        s: Selection coefficient
        shape: Shape parameter
        mean: Mean of the unshifted gamma
        displacement: Shift parameter (typically negative, in [-100, 0])

    Returns:
        DFE density at s
    """
    # Shifted coordinate
    shifted_s = s - displacement

    if shifted_s >= 0:
        return 0.0

    # Use gamma density on shifted coordinate
    abs_shifted = -shifted_s

    if abs_shifted <= 0 or mean <= 0 or shape <= 0:
        return 0.0

    scale = mean / shape
    log_density = ((shape - 1.0) * math.log(abs_shifted)
                   - abs_shifted / scale
                   - shape * math.log(scale)
                   - math.lgamma(shape))

    return math.exp(log_density)


# =============================================================================
# Vectorized DFE grid evaluation (for performance)
# =============================================================================


@numba.njit(cache=True, fastmath=True)
def gamma_density_grid(s_grid: np.ndarray, shape: float, mean: float) -> np.ndarray:
    """Evaluate gamma DFE density on entire S grid (vectorized).

    Much faster than calling gamma_density() in a loop.

    Args:
        s_grid: Array of S values
        shape: Shape parameter
        mean: Mean parameter

    Returns:
        Array of density values
    """
    n = len(s_grid)
    result = np.zeros(n)

    if mean <= 0 or shape <= 0:
        return result

    scale = mean / shape
    log_scale = math.log(scale)
    lgamma_shape = math.lgamma(shape)

    for i in range(n):
        s = s_grid[i]
        if s >= 0:
            result[i] = 0.0
        else:
            abs_s = -s
            log_density = ((shape - 1.0) * math.log(abs_s)
                           - abs_s / scale
                           - shape * log_scale
                           - lgamma_shape)
            result[i] = math.exp(log_density)

    return result


@numba.njit(cache=True, fastmath=True)
def gamma_expo_density_grid(
    s_grid: np.ndarray,
    shape_del: float,
    mean_del: float,
    prop_ben: float,
    mean_ben: float,
) -> np.ndarray:
    """Evaluate GammaExpo DFE density on entire S grid (vectorized)."""
    n = len(s_grid)
    result = np.zeros(n)

    # Pre-compute constants for gamma part
    if mean_del > 0 and shape_del > 0:
        scale_del = mean_del / shape_del
        log_scale_del = math.log(scale_del)
        lgamma_shape_del = math.lgamma(shape_del)
    else:
        scale_del = 0.0
        log_scale_del = 0.0
        lgamma_shape_del = 0.0

    for i in range(n):
        s = s_grid[i]
        if s < 0:
            # Deleterious (gamma)
            if mean_del > 0 and shape_del > 0:
                abs_s = -s
                log_density = ((shape_del - 1.0) * math.log(abs_s)
                               - abs_s / scale_del
                               - shape_del * log_scale_del
                               - lgamma_shape_del)
                result[i] = (1.0 - prop_ben) * math.exp(log_density)
        elif s > 0:
            # Beneficial (exponential)
            if mean_ben > 0:
                result[i] = prop_ben * math.exp(-s / mean_ben) / mean_ben

    return result


@numba.njit(cache=True, fastmath=True)
def gamma_gamma_density_grid(
    s_grid: np.ndarray,
    shape_del: float,
    mean_del: float,
    prop_ben: float,
    shape_ben: float,
    mean_ben: float,
) -> np.ndarray:
    """Evaluate GammaGamma DFE density on entire S grid (vectorized)."""
    n = len(s_grid)
    result = np.zeros(n)

    # Pre-compute constants for deleterious gamma
    if mean_del > 0 and shape_del > 0:
        scale_del = mean_del / shape_del
        log_scale_del = math.log(scale_del)
        lgamma_shape_del = math.lgamma(shape_del)
    else:
        scale_del = 0.0
        log_scale_del = 0.0
        lgamma_shape_del = 0.0

    # Pre-compute constants for beneficial gamma
    if mean_ben > 0 and shape_ben > 0:
        scale_ben = mean_ben / shape_ben
        log_scale_ben = math.log(scale_ben)
        lgamma_shape_ben = math.lgamma(shape_ben)
    else:
        scale_ben = 0.0
        log_scale_ben = 0.0
        lgamma_shape_ben = 0.0

    for i in range(n):
        s = s_grid[i]
        if s < 0:
            # Deleterious (gamma)
            if mean_del > 0 and shape_del > 0:
                abs_s = -s
                log_density = ((shape_del - 1.0) * math.log(abs_s)
                               - abs_s / scale_del
                               - shape_del * log_scale_del
                               - lgamma_shape_del)
                result[i] = (1.0 - prop_ben) * math.exp(log_density)
        elif s > 0:
            # Beneficial (gamma)
            if mean_ben > 0 and shape_ben > 0:
                log_density = ((shape_ben - 1.0) * math.log(s)
                               - s / scale_ben
                               - shape_ben * log_scale_ben
                               - lgamma_shape_ben)
                result[i] = prop_ben * math.exp(log_density)

    return result


@numba.njit(cache=True, fastmath=True)
def displaced_gamma_density_grid(
    s_grid: np.ndarray,
    shape: float,
    mean: float,
    displacement: float,
) -> np.ndarray:
    """Evaluate DisplacedGamma DFE density on entire S grid (vectorized)."""
    n = len(s_grid)
    result = np.zeros(n)

    if mean <= 0 or shape <= 0:
        return result

    scale = mean / shape
    log_scale = math.log(scale)
    lgamma_shape = math.lgamma(shape)

    for i in range(n):
        shifted_s = s_grid[i] - displacement
        if shifted_s >= 0:
            result[i] = 0.0
        else:
            abs_shifted = -shifted_s
            log_density = ((shape - 1.0) * math.log(abs_shifted)
                           - abs_shifted / scale
                           - shape * log_scale
                           - lgamma_shape)
            result[i] = math.exp(log_density)

    return result


# =============================================================================
# Expected SFS under DFE
# =============================================================================


@numba.njit(cache=True)
def _integrate_over_dfe_sfs(
    dfe_values: np.ndarray,
    s_grid: np.ndarray,
    precalc_counts: np.ndarray,
    s_widths: np.ndarray | None = None,
) -> np.ndarray:
    """Integrate expected SFS over DFE using pre-computed values.

    Uses rectangle/midpoint rule matching GRAPES implementation.

    Args:
        dfe_values: DFE density evaluated at each S in s_grid
        s_grid: Grid of S values (midpoints of rectangles)
        precalc_counts: Pre-computed expected counts, shape (n_s, n-1)
        s_widths: Pre-computed widths for each S point (optional, computed if None)

    Returns:
        Expected SFS (unfolded), length n-1
    """
    n_s = len(s_grid)
    n_bins = precalc_counts.shape[1]

    expected_sfs = np.zeros(n_bins)

    # Rectangle/midpoint integration over S grid
    for j in range(n_bins):
        integral = 0.0
        for i in range(n_s):
            # Compute width if not provided
            if s_widths is not None:
                width = s_widths[i]
            elif i == 0:
                width = s_grid[i + 1] - s_grid[i]
            elif i == n_s - 1:
                width = s_grid[i] - s_grid[i - 1]
            else:
                width = (s_grid[i + 1] - s_grid[i - 1]) / 2.0

            integral += dfe_values[i] * precalc_counts[i, j] * width

        expected_sfs[j] = integral

    return expected_sfs


@numba.njit(cache=True)
def _integrate_over_dfe_fixation(
    dfe_values: np.ndarray,
    s_grid: np.ndarray,
    s_widths: np.ndarray | None = None,
) -> float:
    """Integrate fixation probability over DFE.

    Uses rectangle/midpoint rule matching GRAPES implementation.

    Args:
        dfe_values: DFE density evaluated at each S in s_grid
        s_grid: Grid of S values
        s_widths: Pre-computed widths for each S point (optional)

    Returns:
        Expected fixation rate (relative to neutral)
    """
    n_s = len(s_grid)

    integral = 0.0
    for i in range(n_s):
        # Compute width if not provided
        if s_widths is not None:
            width = s_widths[i]
        elif i == 0:
            width = s_grid[i + 1] - s_grid[i]
        elif i == n_s - 1:
            width = s_grid[i] - s_grid[i - 1]
        else:
            width = (s_grid[i + 1] - s_grid[i - 1]) / 2.0

        integral += dfe_values[i] * fixation_probability(s_grid[i]) * width

    return integral


@numba.njit(cache=True)
def _integrate_over_dfe_adaptive_fixation(
    dfe_values: np.ndarray,
    s_grid: np.ndarray,
    threshold: float = 5.0,
    s_widths: np.ndarray | None = None,
) -> float:
    """Integrate fixation probability for adaptive mutations only.

    Only counts fixations where S > threshold (beneficial enough to be
    considered adaptive). Uses rectangle/midpoint rule.

    Args:
        dfe_values: DFE density evaluated at each S in s_grid
        s_grid: Grid of S values
        threshold: Minimum S to count as adaptive (default 5.0)
        s_widths: Pre-computed widths for each S point (optional)

    Returns:
        Expected adaptive fixation rate
    """
    n_s = len(s_grid)

    integral = 0.0
    for i in range(n_s):
        # Only integrate where S > threshold
        if s_grid[i] <= threshold:
            continue

        # Compute width if not provided
        if s_widths is not None:
            width = s_widths[i]
        elif i == 0:
            width = s_grid[i + 1] - s_grid[i]
        elif i == n_s - 1:
            width = s_grid[i] - s_grid[i - 1]
        else:
            width = (s_grid[i + 1] - s_grid[i - 1]) / 2.0

        integral += dfe_values[i] * fixation_probability(s_grid[i]) * width

    return integral


@numba.njit(cache=True)
def _integrate_over_dfe_non_adaptive_fixation(
    dfe_values: np.ndarray,
    s_grid: np.ndarray,
    threshold: float = 5.0,
    s_widths: np.ndarray | None = None,
) -> float:
    """Integrate fixation probability for non-adaptive mutations only.

    Only counts fixations where S <= threshold (weakly selected, not
    strongly beneficial). Uses rectangle/midpoint rule.

    This is used in the FWW-style alpha calculation:
    omega_na = integral(DFE(S) × P_fix(S)) for S <= threshold

    Args:
        dfe_values: DFE density evaluated at each S in s_grid
        s_grid: Grid of S values
        threshold: Maximum S to count as non-adaptive (default 5.0)
        s_widths: Pre-computed widths for each S point (optional)

    Returns:
        Expected non-adaptive fixation rate (omega_na)
    """
    n_s = len(s_grid)

    integral = 0.0
    for i in range(n_s):
        # Only integrate where S <= threshold
        if s_grid[i] > threshold:
            continue

        # Compute width if not provided
        if s_widths is not None:
            width = s_widths[i]
        elif i == 0:
            width = s_grid[i + 1] - s_grid[i]
        elif i == n_s - 1:
            width = s_grid[i] - s_grid[i - 1]
        else:
            width = (s_grid[i + 1] - s_grid[i - 1]) / 2.0

        integral += dfe_values[i] * fixation_probability(s_grid[i]) * width

    return integral


# =============================================================================
# Likelihood computation
# =============================================================================


@numba.njit(cache=True, fastmath=True)
def poisson_log_likelihood(observed: np.ndarray, expected: np.ndarray) -> float:
    """Compute Poisson log-likelihood for SFS counts.

    Args:
        observed: Observed SFS counts
        expected: Expected SFS counts under model

    Returns:
        Log-likelihood value
    """
    ll = 0.0
    for i in range(len(observed)):
        obs = observed[i]
        exp = expected[i]

        if exp <= 0:
            if obs > 0:
                return -1e10  # Very negative likelihood
            continue

        # Poisson log-likelihood: obs*log(exp) - exp - log(obs!)
        if obs > 0:
            ll += obs * math.log(exp) - exp - math.lgamma(obs + 1)
        else:
            ll += -exp

    return ll


@numba.njit(cache=True, fastmath=True)
def multinomial_log_likelihood(observed: np.ndarray, expected: np.ndarray) -> float:
    """Compute multinomial log-likelihood for SFS proportions.

    Args:
        observed: Observed SFS counts
        expected: Expected SFS counts (will be normalized)

    Returns:
        Log-likelihood value
    """
    total_obs = observed.sum()
    total_exp = expected.sum()

    if total_obs == 0 or total_exp <= 0:
        return 0.0

    # Normalize expected to proportions
    ll = 0.0
    for i in range(len(observed)):
        obs = observed[i]
        prop = expected[i] / total_exp

        if prop <= 0:
            if obs > 0:
                return -1e10
            continue

        if obs > 0:
            ll += obs * math.log(prop)

    return ll
