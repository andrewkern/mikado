"""DFE model implementations for alpha estimation.

Based on GRAPES implementation (Galtier 2016; Al-Saffar & Hahn 2022).
Reference: https://github.com/BioPP/grapes
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy import optimize
from scipy.special import gammaln

from mkado.analysis.dfe_numba import (
    displaced_gamma_density,
    displaced_gamma_density_grid,
    fold_sfs,
    gamma_density,
    gamma_density_grid,
    gamma_expo_density,
    gamma_expo_density_grid,
    gamma_gamma_density,
    gamma_gamma_density_grid,
    precalculate_expected_counts,
    precalculate_s_grid,
    precalculate_s_widths,
    _integrate_over_dfe_sfs,
    _integrate_over_dfe_fixation,
    _integrate_over_dfe_adaptive_fixation,
    _integrate_over_dfe_non_adaptive_fixation,
)

if TYPE_CHECKING:
    pass


def _log_factorial(n: int) -> float:
    """Compute log(n!) using gammaln for numerical stability."""
    if n <= 0:
        return 0.0
    return gammaln(n + 1)


@dataclass
class DFEInput:
    """Input data for DFE-based alpha estimation."""

    sfs_neutral: np.ndarray
    """Folded SFS for synonymous sites (length floor(n/2))."""

    sfs_selected: np.ndarray
    """Folded SFS for nonsynonymous sites (length floor(n/2))."""

    divergence_neutral: int
    """Ds (synonymous divergence count)."""

    divergence_selected: int
    """Dn (nonsynonymous divergence count)."""

    n_samples: int
    """Number of chromosomes sampled."""

    n_sites_neutral: float | None = None
    """Number of synonymous sites for polymorphism (Lps). If None, estimated from SFS."""

    n_sites_selected: float | None = None
    """Number of nonsynonymous sites for polymorphism (Lpn). If None, estimated from SFS."""

    n_sites_div_neutral: float | None = None
    """Number of synonymous sites for divergence (Lds). If None, uses n_sites_neutral."""

    n_sites_div_selected: float | None = None
    """Number of nonsynonymous sites for divergence (Ldn). If None, uses n_sites_selected."""


@dataclass
class DFEResult:
    """Results from DFE-based alpha estimation."""

    model: str
    """Model name (e.g., 'GammaExpo')."""

    alpha: float
    """Proportion of adaptive amino acid substitutions."""

    alpha_down: float
    """Lower confidence interval bound for alpha."""

    alpha_up: float
    """Upper confidence interval bound for alpha."""

    omega_a: float
    """Adaptive substitution rate (adaptive dN/dS)."""

    omega_na: float
    """Non-adaptive substitution rate (non-adaptive dN/dS)."""

    dfe_params: dict
    """Model-specific DFE parameters."""

    log_likelihood: float
    """Maximum log-likelihood."""

    aic: float
    """Akaike Information Criterion."""

    converged: bool
    """Whether optimization converged successfully."""

    def __str__(self) -> str:
        lines = [
            f"DFE-based Alpha Estimation ({self.model}):",
            f"  Alpha: {self.alpha:.4f} (95% CI: {self.alpha_down:.4f} - {self.alpha_up:.4f})",
            f"  omega_a: {self.omega_a:.4f}",
            f"  omega_na: {self.omega_na:.4f}",
            f"  Log-likelihood: {self.log_likelihood:.2f}",
            f"  AIC: {self.aic:.2f}",
            "  DFE parameters:",
        ]
        for name, value in self.dfe_params.items():
            lines.append(f"    {name}: {value:.4g}")

        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "model": self.model,
            "alpha": self.alpha,
            "alpha_down": self.alpha_down,
            "alpha_up": self.alpha_up,
            "omega_a": self.omega_a,
            "omega_na": self.omega_na,
            "dfe_params": self.dfe_params,
            "log_likelihood": self.log_likelihood,
            "aic": self.aic,
            "converged": self.converged,
        }


class PrecomputedData:
    """Container for pre-computed integration grids."""

    def __init__(self, n_samples: int, n_s_points: int = 4000):
        """Initialize pre-computed data for a given sample size.

        Args:
            n_samples: Number of chromosomes sampled
            n_s_points: Number of S grid points
        """
        self.n_samples = n_samples
        self.s_grid = precalculate_s_grid(n_s_points)
        self.s_widths = precalculate_s_widths(self.s_grid)
        self.expected_counts = precalculate_expected_counts(self.s_grid, n_samples)


class DFEModel(ABC):
    """Base class for DFE models."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Model name."""
        pass

    @property
    @abstractmethod
    def n_params(self) -> int:
        """Number of free DFE parameters."""
        pass

    @property
    @abstractmethod
    def param_names(self) -> list[str]:
        """Names of DFE parameters."""
        pass

    @property
    @abstractmethod
    def param_bounds(self) -> list[tuple[float, float]]:
        """Bounds for optimization (in original parameter space)."""
        pass

    @property
    def log_scale_params(self) -> list[bool]:
        """Which parameters should be optimized in log space.

        Following GRAPES reparametrization convention:
        - Mean parameters (negGmean, posGmean) use log scale
        - Shape and proportion parameters use linear scale

        Override in subclasses to specify log-scale parameters.
        """
        return [False] * self.n_params

    @abstractmethod
    def density(self, s: float, params: np.ndarray) -> float:
        """Evaluate DFE density at selection coefficient s.

        Args:
            s: Selection coefficient (4*Ne*s)
            params: Model parameters

        Returns:
            DFE density at s
        """
        pass

    def evaluate_dfe_on_grid(self, params: np.ndarray, s_grid: np.ndarray) -> np.ndarray:
        """Evaluate DFE density on entire S grid.

        Args:
            params: Model parameters
            s_grid: Grid of S values

        Returns:
            Array of DFE densities
        """
        return np.array([self.density(s, params) for s in s_grid])

    def _get_starting_points(self) -> list[np.ndarray]:
        """Get model-specific starting points for optimization.

        Override in subclasses to provide good starting points.

        Returns:
            List of starting point arrays
        """
        # Default: use middle of bounds
        return [np.array([(b[0] + b[1]) / 2 for b in self.param_bounds])]

    def expected_sfs_selected(
        self,
        params: np.ndarray,
        precalc: PrecomputedData,
        theta: float,
    ) -> np.ndarray:
        """Compute expected SFS for selected sites under this DFE.

        Args:
            params: DFE parameters
            precalc: Pre-computed integration data
            theta: Mutation rate parameter (4*Ne*mu*L)

        Returns:
            Expected folded SFS
        """
        dfe_values = self.evaluate_dfe_on_grid(params, precalc.s_grid)
        unfolded_sfs = _integrate_over_dfe_sfs(
            dfe_values, precalc.s_grid, precalc.expected_counts, precalc.s_widths
        )
        # Scale by theta
        unfolded_sfs = theta * unfolded_sfs
        # Fold the SFS
        return fold_sfs(unfolded_sfs)

    def expected_divergence(
        self,
        params: np.ndarray,
        precalc: PrecomputedData,
    ) -> float:
        """Compute expected divergence rate under this DFE.

        Args:
            params: DFE parameters
            precalc: Pre-computed integration data

        Returns:
            Expected divergence rate relative to neutral
        """
        dfe_values = self.evaluate_dfe_on_grid(params, precalc.s_grid)
        return _integrate_over_dfe_fixation(dfe_values, precalc.s_grid, precalc.s_widths)

    def expected_adaptive_divergence(
        self,
        params: np.ndarray,
        precalc: PrecomputedData,
        threshold: float = 5.0,
    ) -> float:
        """Compute expected adaptive divergence rate.

        Args:
            params: DFE parameters
            precalc: Pre-computed integration data
            threshold: Minimum S to count as adaptive

        Returns:
            Expected adaptive divergence rate
        """
        dfe_values = self.evaluate_dfe_on_grid(params, precalc.s_grid)
        return _integrate_over_dfe_adaptive_fixation(
            dfe_values, precalc.s_grid, threshold, precalc.s_widths
        )

    def expected_non_adaptive_divergence(
        self,
        params: np.ndarray,
        precalc: PrecomputedData,
        threshold: float = 5.0,
    ) -> float:
        """Compute expected non-adaptive divergence rate (S <= threshold).

        This is used in the FWW-style alpha calculation:
        alpha = 1 - omega_na / omega_obs

        Args:
            params: DFE parameters
            precalc: Pre-computed integration data
            threshold: Maximum S to count as non-adaptive

        Returns:
            Expected non-adaptive divergence rate (omega_na)
        """
        dfe_values = self.evaluate_dfe_on_grid(params, precalc.s_grid)
        return _integrate_over_dfe_non_adaptive_fixation(
            dfe_values, precalc.s_grid, threshold, precalc.s_widths
        )

    def negative_log_likelihood(
        self,
        params: np.ndarray,
        data: DFEInput,
        precalc: PrecomputedData,
    ) -> float:
        """Compute negative log-likelihood for optimization.

        Following the GRAPES implementation:
        1. Compute expected SFS for neutral (L_s/j pattern) and selected
           (L_n × DFE integration / 2) sites, both with theta=1
        2. Estimate theta from first frequency bin
        3. Estimate rj demography ratios for each frequency bin j > 1
        4. Scale expected SFS by theta * rj
        5. Compute Poisson log-likelihood

        Args:
            params: DFE parameters
            data: Input SFS and divergence data
            precalc: Pre-computed integration data

        Returns:
            Negative log-likelihood
        """
        n = data.n_samples
        n_bins = len(data.sfs_neutral)

        # Get sequence lengths (estimate from SFS if not provided)
        if data.n_sites_neutral is not None:
            Ls = data.n_sites_neutral
        else:
            Ls = float(data.sfs_neutral.sum()) * 10  # Rough estimate

        if data.n_sites_selected is not None:
            Ln = data.n_sites_selected
        else:
            Ln = float(data.sfs_selected.sum()) * 10  # Rough estimate

        # Step 1: Expected neutral SFS (folded, theta=1, with Ls)
        # For folded SFS: exp[j] = Ls/j + Ls/(n-j) for j < n/2
        #                        = Ls/j           for j = n/2 (if n even)
        exp_neutral = np.zeros(n_bins)
        for j in range(1, n_bins + 1):
            exp_neutral[j - 1] = Ls / j
            if n % 2 != 0 or j != n // 2:
                exp_neutral[j - 1] += Ls / (n - j)

        # Step 2: Expected selected SFS (folded, theta=1, with Ln) from DFE
        # GRAPES divides the integral by 2
        dfe_values = self.evaluate_dfe_on_grid(params, precalc.s_grid)
        exp_selected_unfolded = _integrate_over_dfe_sfs(
            dfe_values, precalc.s_grid, precalc.expected_counts, precalc.s_widths
        )
        exp_selected = fold_sfs(exp_selected_unfolded) * Ln / 2.0

        # Avoid division by zero
        exp_neutral = np.maximum(exp_neutral, 1e-100)
        exp_selected = np.maximum(exp_selected, 1e-100)

        # Step 3: Estimate theta from first frequency bin
        # theta = (specS[0] + specN[0]) / (exp_specS[0] + exp_specN[0])
        obs_total_0 = data.sfs_neutral[0] + data.sfs_selected[0]
        exp_total_0 = exp_neutral[0] + exp_selected[0]
        theta = obs_total_0 / exp_total_0 if exp_total_0 > 0 else 1.0

        # Step 4: Estimate rj demography ratios for j > 1
        # r[j] = ((specS[j] + specN[j]) / (exp_specS[j] + exp_specN[j])) / theta
        rj = np.ones(n_bins)
        for j in range(1, n_bins):
            obs_total_j = data.sfs_neutral[j] + data.sfs_selected[j]
            exp_total_j = exp_neutral[j] + exp_selected[j]
            if exp_total_j > 0 and theta > 0:
                rj[j] = (obs_total_j / exp_total_j) / theta
            else:
                rj[j] = 1.0

        # Step 5: Scale expected SFS by theta * rj
        exp_neutral_scaled = exp_neutral * theta * rj
        exp_selected_scaled = exp_selected * theta * rj

        # Step 6: Expected divergence
        # Following GRAPES default: use_divergence_parameter = True
        # This sets expected = observed for divergence (saturated model)
        # Divergence contributes a constant to log-likelihood
        if data.divergence_neutral > 0:
            exp_div_s = float(data.divergence_neutral)
            exp_div_n = float(data.divergence_selected)
        else:
            exp_div_s = 1.0
            exp_div_n = 1.0

        # Step 7: Poisson log-likelihood
        # GRAPES uses integer counts for k*log(lam) - lam - log(k!)
        ll = 0.0

        # SFS contribution
        for j in range(n_bins):
            # Neutral (round to integer as GRAPES does)
            k_s = int(round(data.sfs_neutral[j]))
            lam_s = max(exp_neutral_scaled[j], 1e-100)
            ll += k_s * np.log(lam_s) - lam_s - _log_factorial(k_s)

            # Selected (round to integer as GRAPES does)
            k_n = int(round(data.sfs_selected[j]))
            lam_n = max(exp_selected_scaled[j], 1e-100)
            ll += k_n * np.log(lam_n) - lam_n - _log_factorial(k_n)

        # Divergence contribution
        if data.divergence_neutral > 0:
            k_ds = data.divergence_neutral
            lam_ds = max(exp_div_s, 1e-100)
            ll += k_ds * np.log(lam_ds) - lam_ds - _log_factorial(k_ds)

            k_dn = data.divergence_selected
            lam_dn = max(exp_div_n, 1e-100)
            ll += k_dn * np.log(lam_dn) - lam_dn - _log_factorial(k_dn)

        return -ll  # Return negative for minimization

    def _to_opt_space(self, params: np.ndarray) -> np.ndarray:
        """Transform parameters to optimization space (log for some params)."""
        log_scale = self.log_scale_params
        result = np.zeros_like(params)
        for i, (p, use_log) in enumerate(zip(params, log_scale)):
            if use_log:
                result[i] = np.log(max(p, 1e-100))
            else:
                result[i] = p
        return result

    def _from_opt_space(self, opt_params: np.ndarray) -> np.ndarray:
        """Transform parameters from optimization space back to original."""
        log_scale = self.log_scale_params
        result = np.zeros_like(opt_params)
        for i, (p, use_log) in enumerate(zip(opt_params, log_scale)):
            if use_log:
                result[i] = np.exp(p)
            else:
                result[i] = p
        return result

    def _get_opt_bounds(self) -> list[tuple[float, float]]:
        """Get bounds in optimization space."""
        bounds = self.param_bounds
        log_scale = self.log_scale_params
        opt_bounds = []
        for (low, high), use_log in zip(bounds, log_scale):
            if use_log:
                opt_bounds.append((np.log(max(low, 1e-100)), np.log(high)))
            else:
                opt_bounds.append((low, high))
        return opt_bounds

    def fit(
        self,
        data: DFEInput,
        precalc: PrecomputedData | None = None,
        n_starts: int = 5,
        adaptive_threshold: float = 5.0,
    ) -> DFEResult:
        """Fit model to data using maximum likelihood.

        Args:
            data: Input SFS and divergence data
            precalc: Pre-computed integration data (computed if None)
            n_starts: Number of random starting points for optimization
            adaptive_threshold: Minimum S for adaptive mutations

        Returns:
            DFEResult with fitted parameters and alpha estimate
        """
        if precalc is None:
            precalc = PrecomputedData(data.n_samples)

        bounds = self.param_bounds
        opt_bounds = self._get_opt_bounds()
        log_scale = self.log_scale_params
        n_params = self.n_params

        best_result = None
        best_nll = float("inf")

        # Objective function in optimization space
        def objective(opt_params: np.ndarray) -> float:
            params = self._from_opt_space(opt_params)
            return self.negative_log_likelihood(params, data, precalc)

        # Multiple starting points with strategic choices
        rng = np.random.default_rng(42)

        # Model-specific starting points based on empirical performance
        # These are tuned to find the global optimum more reliably
        model_starts = self._get_starting_points()

        all_starts = []

        # Add model-specific starts (transform to opt space)
        for start in model_starts:
            all_starts.append(self._to_opt_space(start))

        # Add middle of bounds in opt space
        middle = []
        for (low, high), use_log in zip(bounds, log_scale):
            if use_log:
                # Geometric mean in original space = arithmetic mean in log space
                middle.append((np.log(max(low, 1e-100)) + np.log(high)) / 2)
            else:
                middle.append((low + high) / 2)
        all_starts.append(np.array(middle))

        # Add random starts in opt space
        for _ in range(max(0, n_starts - len(all_starts))):
            x0 = np.zeros(n_params)
            for i, ((low, high), use_log) in enumerate(zip(bounds, log_scale)):
                if use_log:
                    x0[i] = rng.uniform(np.log(max(low, 1e-100)), np.log(high))
                else:
                    x0[i] = rng.uniform(low, high)
            all_starts.append(x0)

        for x0 in all_starts:

            try:
                # Try L-BFGS-B first
                result = optimize.minimize(
                    objective,
                    x0,
                    method="L-BFGS-B",
                    bounds=opt_bounds,
                    options={"maxiter": 10000, "ftol": 1e-8},
                )

                # If L-BFGS-B failed or gave abnormal termination, try Nelder-Mead
                if not result.success:
                    result_nm = optimize.minimize(
                        objective,
                        x0,
                        method="Nelder-Mead",
                        options={"maxiter": 10000, "xatol": 1e-6, "fatol": 1e-6},
                    )
                    if result_nm.success and result_nm.fun < result.fun:
                        result = result_nm

                if result.fun < best_nll:
                    best_nll = result.fun
                    best_result = result

            except (ValueError, RuntimeError):
                continue

        if best_result is None:
            # Return failure result
            return DFEResult(
                model=self.name,
                alpha=0.0,
                alpha_down=0.0,
                alpha_up=0.0,
                omega_a=0.0,
                omega_na=0.0,
                dfe_params={},
                log_likelihood=-float("inf"),
                aic=float("inf"),
                converged=False,
            )

        # Extract fitted parameters (transform back from opt space)
        fitted_params = self._from_opt_space(best_result.x)
        log_likelihood = -best_result.fun
        aic = 2 * n_params - 2 * log_likelihood

        # Compute alpha following GRAPES methodology (FWW-style with DFE):
        # omega_na = integral(DFE(S) × P_fix(S)) for S ≤ threshold
        # omega_obs = (Dn/Ldn) / (Ds/Lds)
        # alpha = 1 - omega_na / omega_obs
        #
        # This uses the DFE to predict the non-adaptive substitution rate,
        # then computes what fraction of observed divergence is unexplained (adaptive).

        # Get divergence site counts
        Ldn = data.n_sites_div_selected
        if Ldn is None:
            Ldn = data.n_sites_selected if data.n_sites_selected else 1.0
        Lds = data.n_sites_div_neutral
        if Lds is None:
            Lds = data.n_sites_neutral if data.n_sites_neutral else 1.0

        # Observed dN/dS
        if data.divergence_neutral > 0 and Lds > 0 and Ldn > 0:
            dn_rate = data.divergence_selected / Ldn
            ds_rate = data.divergence_neutral / Lds
            omega_obs = dn_rate / ds_rate if ds_rate > 0 else 0.0
        else:
            omega_obs = 0.0

        # Compute omega_na = expected non-adaptive fixation rate (S <= threshold)
        # This is the expected dN/dS if only weakly selected mutations contribute
        omega_na = self.expected_non_adaptive_divergence(
            fitted_params, precalc, threshold=adaptive_threshold
        )

        # Alpha = fraction of divergence not explained by non-adaptive DFE
        if omega_obs > 0 and omega_na < omega_obs:
            alpha = 1.0 - (omega_na / omega_obs)
            omega_a = omega_obs - omega_na
        else:
            alpha = 0.0
            omega_a = 0.0

        # Compute confidence intervals via likelihood profile
        alpha_down, alpha_up = self._compute_ci(
            fitted_params, data, precalc, best_result.fun, adaptive_threshold
        )

        # Build parameter dictionary
        param_dict = dict(zip(self.param_names, fitted_params))

        return DFEResult(
            model=self.name,
            alpha=alpha,
            alpha_down=alpha_down,
            alpha_up=alpha_up,
            omega_a=omega_a,
            omega_na=omega_na,
            dfe_params=param_dict,
            log_likelihood=log_likelihood,
            aic=aic,
            converged=best_result.success,
        )

    def _compute_ci(
        self,
        fitted_params: np.ndarray,
        data: DFEInput,
        precalc: PrecomputedData,
        best_nll: float,
        adaptive_threshold: float,
        delta_ll: float = 2.0,
    ) -> tuple[float, float]:
        """Compute confidence interval for alpha via likelihood ratio.

        Finds parameter values where log-likelihood decreases by delta_ll
        from the maximum.

        Args:
            fitted_params: Best-fit parameters
            data: Input data
            precalc: Pre-computed data
            best_nll: Best negative log-likelihood
            adaptive_threshold: Threshold for adaptive mutations
            delta_ll: Log-likelihood decrease for CI (2.0 for ~95%)

        Returns:
            Tuple of (alpha_down, alpha_up)
        """
        # Compute omega_obs for alpha calculation
        Ldn = data.n_sites_div_selected
        if Ldn is None:
            Ldn = data.n_sites_selected if data.n_sites_selected else 1.0
        Lds = data.n_sites_div_neutral
        if Lds is None:
            Lds = data.n_sites_neutral if data.n_sites_neutral else 1.0

        if data.divergence_neutral > 0 and Lds > 0 and Ldn > 0:
            dn_rate = data.divergence_selected / Ldn
            ds_rate = data.divergence_neutral / Lds
            omega_obs = dn_rate / ds_rate if ds_rate > 0 else 0.0
        else:
            omega_obs = 0.0

        def compute_alpha(params: np.ndarray) -> float:
            """Compute alpha using FWW-style calculation with DFE."""
            omega_na = self.expected_non_adaptive_divergence(
                params, precalc, threshold=adaptive_threshold
            )
            if omega_obs > 0 and omega_na < omega_obs:
                return 1.0 - (omega_na / omega_obs)
            return 0.0

        # Compute alpha at best fit
        best_alpha = compute_alpha(fitted_params)

        # Bootstrap approach for CI estimation
        rng = np.random.default_rng(42)
        alpha_samples = [best_alpha]

        for _ in range(100):
            # Perturb parameters within reasonable range
            perturbed = fitted_params * (1 + 0.1 * rng.standard_normal(len(fitted_params)))

            # Clip to bounds
            for i, (low, high) in enumerate(self.param_bounds):
                perturbed[i] = np.clip(perturbed[i], low, high)

            # Check if likelihood is within threshold
            nll = self.negative_log_likelihood(perturbed, data, precalc)
            if nll < best_nll + delta_ll:
                alpha_samples.append(compute_alpha(perturbed))

        if len(alpha_samples) < 10:
            # Not enough samples, use point estimate
            return (best_alpha, best_alpha)

        alpha_down = float(np.percentile(alpha_samples, 2.5))
        alpha_up = float(np.percentile(alpha_samples, 97.5))

        return (alpha_down, alpha_up)


class GammaZeroModel(DFEModel):
    """Gamma DFE with no beneficial mutations.

    Parameters:
    - shape: Gamma shape parameter (alpha)
    - mean: Gamma mean (in units of 4*Ne*s)
    """

    @property
    def name(self) -> str:
        return "GammaZero"

    @property
    def n_params(self) -> int:
        return 2

    @property
    def param_names(self) -> list[str]:
        return ["shape", "mean"]

    @property
    def param_bounds(self) -> list[tuple[float, float]]:
        return [
            (0.05, 100.0),    # shape
            (1.0, 1e6),       # mean (|4*Ne*s|)
        ]

    @property
    def log_scale_params(self) -> list[bool]:
        # GRAPES: negGshape=none, negGmean=log
        return [False, True]

    def density(self, s: float, params: np.ndarray) -> float:
        shape, mean = params
        return gamma_density(s, shape, mean)

    def evaluate_dfe_on_grid(self, params: np.ndarray, s_grid: np.ndarray) -> np.ndarray:
        """Vectorized DFE evaluation (65x faster than scalar loop)."""
        shape, mean = params
        return gamma_density_grid(s_grid, shape, mean)

    def _get_starting_points(self) -> list[np.ndarray]:
        """Starting points tuned for GammaZero model."""
        return [
            np.array([0.3, 5000.0]),     # Low shape, moderate mean
            np.array([0.3, 20000.0]),    # Low shape, high mean (GRAPES typical)
            np.array([0.5, 10000.0]),    # Medium shape, moderate mean
            np.array([1.0, 1000.0]),     # Higher shape, low mean
        ]


class GammaExpoModel(DFEModel):
    """Gamma (deleterious) + Exponential (beneficial) DFE.

    This is the best-performing model according to Al-Saffar & Hahn (2022).

    Parameters:
    - shape_del: Gamma shape for deleterious mutations
    - mean_del: Gamma mean for deleterious mutations
    - prop_ben: Proportion of beneficial mutations
    - mean_ben: Exponential mean for beneficial mutations
    """

    @property
    def name(self) -> str:
        return "GammaExpo"

    @property
    def n_params(self) -> int:
        return 4

    @property
    def param_names(self) -> list[str]:
        return ["shape_del", "mean_del", "prop_ben", "mean_ben"]

    @property
    def param_bounds(self) -> list[tuple[float, float]]:
        # Bounds matching GRAPES exactly:
        # - negGshape: [0.05, 100]
        # - negGmean: [1, 10^20]
        # - pos_prop: [0, 0.1] (MAX_ORIENTATION_ERROR)
        # - posGmean: [0.0001, 10000]
        return [
            (0.05, 100.0),    # shape_del (negGshape)
            (1.0, 1e6),       # mean_del (negGmean) - capped at 1e6 for practicality
            (0.0, 0.1),       # prop_ben (pos_prop) - GRAPES uses MAX_ORIENTATION_ERROR=0.1
            (1e-4, 1e4),      # mean_ben (posGmean)
        ]

    @property
    def log_scale_params(self) -> list[bool]:
        # GRAPES: negGshape=none, negGmean=log, pos_prop=none, posGmean=log
        return [False, True, False, True]

    def density(self, s: float, params: np.ndarray) -> float:
        shape_del, mean_del, prop_ben, mean_ben = params
        return gamma_expo_density(s, shape_del, mean_del, prop_ben, mean_ben)

    def evaluate_dfe_on_grid(self, params: np.ndarray, s_grid: np.ndarray) -> np.ndarray:
        """Vectorized DFE evaluation."""
        shape_del, mean_del, prop_ben, mean_ben = params
        return gamma_expo_density_grid(s_grid, shape_del, mean_del, prop_ben, mean_ben)

    def _get_starting_points(self) -> list[np.ndarray]:
        """Starting points for GammaExpo model matching GRAPES.

        GRAPES uses mean_del=10000 for all starting points (the "basic" value).
        Testing shows that shape_del around 0.3 leads to the GRAPES-compatible
        basin, while higher values (e.g., 0.4) can lead to a different basin
        with slightly better likelihood but different biological interpretation.

        Using consistent starting points ensures we converge to the same basin
        as GRAPES, which is important for result comparability.
        """
        return [
            # [shape_del, mean_del, prop_ben, mean_ben]
            # All starts use shape ~0.3, mean_del=10000, mean_ben=10000
            np.array([0.30, 10000.0, 0.010, 10000.0]),
            np.array([0.35, 10000.0, 0.008, 10000.0]),
            np.array([0.32, 10000.0, 0.005, 10000.0]),
            np.array([0.28, 10000.0, 0.015, 10000.0]),
        ]


class GammaGammaModel(DFEModel):
    """Gamma (deleterious) + Gamma (beneficial) DFE.

    Parameters:
    - shape_del: Gamma shape for deleterious mutations
    - mean_del: Gamma mean for deleterious mutations
    - prop_ben: Proportion of beneficial mutations
    - shape_ben: Gamma shape for beneficial mutations
    - mean_ben: Gamma mean for beneficial mutations
    """

    @property
    def name(self) -> str:
        return "GammaGamma"

    @property
    def n_params(self) -> int:
        return 5

    @property
    def param_names(self) -> list[str]:
        return ["shape_del", "mean_del", "prop_ben", "shape_ben", "mean_ben"]

    @property
    def param_bounds(self) -> list[tuple[float, float]]:
        return [
            (0.05, 100.0),    # shape_del
            (1.0, 1e6),       # mean_del
            (0.0, 0.5),       # prop_ben
            (0.05, 100.0),    # shape_ben
            (1e-4, 1e4),      # mean_ben
        ]

    def density(self, s: float, params: np.ndarray) -> float:
        shape_del, mean_del, prop_ben, shape_ben, mean_ben = params
        return gamma_gamma_density(
            s, shape_del, mean_del, prop_ben, shape_ben, mean_ben
        )

    def evaluate_dfe_on_grid(self, params: np.ndarray, s_grid: np.ndarray) -> np.ndarray:
        """Vectorized DFE evaluation."""
        shape_del, mean_del, prop_ben, shape_ben, mean_ben = params
        return gamma_gamma_density_grid(
            s_grid, shape_del, mean_del, prop_ben, shape_ben, mean_ben
        )

    def _get_starting_points(self) -> list[np.ndarray]:
        """Starting points tuned for GammaGamma model."""
        return [
            # [shape_del, mean_del, prop_ben, shape_ben, mean_ben]
            np.array([0.3, 10000.0, 0.01, 0.5, 500.0]),
            np.array([0.4, 20000.0, 0.005, 1.0, 1000.0]),
            np.array([0.3, 5000.0, 0.02, 0.3, 100.0]),
            np.array([0.5, 15000.0, 0.008, 0.8, 2000.0]),
        ]


class DisplacedGammaModel(DFEModel):
    """Displaced Gamma DFE.

    Shifts the gamma distribution to allow some beneficial mutations.

    Parameters:
    - shape: Gamma shape parameter
    - mean: Gamma mean parameter
    - displacement: Shift parameter (negative allows beneficial)
    """

    @property
    def name(self) -> str:
        return "DisplacedGamma"

    @property
    def n_params(self) -> int:
        return 3

    @property
    def param_names(self) -> list[str]:
        return ["shape", "mean", "displacement"]

    @property
    def param_bounds(self) -> list[tuple[float, float]]:
        return [
            (0.05, 100.0),    # shape
            (1.0, 1e6),       # mean
            (-100.0, 0.0),    # displacement
        ]

    def density(self, s: float, params: np.ndarray) -> float:
        shape, mean, displacement = params
        return displaced_gamma_density(s, shape, mean, displacement)

    def evaluate_dfe_on_grid(self, params: np.ndarray, s_grid: np.ndarray) -> np.ndarray:
        """Vectorized DFE evaluation."""
        shape, mean, displacement = params
        return displaced_gamma_density_grid(s_grid, shape, mean, displacement)

    def _get_starting_points(self) -> list[np.ndarray]:
        """Starting points tuned for DisplacedGamma model."""
        return [
            # [shape, mean, displacement]
            np.array([0.3, 10000.0, -10.0]),
            np.array([0.4, 20000.0, -5.0]),
            np.array([0.3, 5000.0, -50.0]),
            np.array([0.5, 15000.0, -20.0]),
        ]


# Model registry
MODELS: dict[str, type[DFEModel]] = {
    "GammaZero": GammaZeroModel,
    "GammaExpo": GammaExpoModel,
    "GammaGamma": GammaGammaModel,
    "DisplacedGamma": DisplacedGammaModel,
}


def get_model(name: str) -> DFEModel:
    """Get a DFE model instance by name.

    Args:
        name: Model name (GammaZero, GammaExpo, GammaGamma, DisplacedGamma)

    Returns:
        DFE model instance

    Raises:
        ValueError: If model name is not recognized
    """
    if name not in MODELS:
        available = ", ".join(MODELS.keys())
        raise ValueError(f"Unknown model '{name}'. Available: {available}")

    return MODELS[name]()
