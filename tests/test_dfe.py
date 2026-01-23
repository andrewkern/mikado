"""Tests for DFE-based alpha estimation."""

from pathlib import Path

import numpy as np
import pytest

from mkado.analysis.dfe import (
    DFEResult,
    dfe_alpha,
    dfe_alpha_aggregated,
    dfe_alpha_from_sfs,
    compare_models,
    polymorphism_data_to_dfe_input,
    AVAILABLE_MODELS,
)
from mkado.analysis.dfe_models import (
    DFEInput,
    DFEModel,
    GammaZeroModel,
    GammaExpoModel,
    GammaGammaModel,
    DisplacedGammaModel,
    PrecomputedData,
    get_model,
    MODELS,
)
from mkado.analysis.dfe_numba import (
    binomial_coeff,
    binomial_prob,
    sojourn_density,
    fixation_probability,
    gamma_density,
    exponential_density,
    gamma_expo_density,
    gamma_gamma_density,
    displaced_gamma_density,
    fold_sfs,
    precalculate_s_grid,
    poisson_log_likelihood,
)
from mkado.analysis.asymptotic import PolymorphismData


class TestDFENumbaKernels:
    """Tests for numba-compiled kernel functions."""

    def test_binomial_coeff_basic(self) -> None:
        """Test binomial coefficient computation."""
        assert binomial_coeff(5, 0) == 1
        assert binomial_coeff(5, 5) == 1
        assert abs(binomial_coeff(5, 2) - 10) < 1e-6
        assert abs(binomial_coeff(10, 3) - 120) < 1e-6
        # Edge cases
        assert binomial_coeff(5, 6) == 0
        assert binomial_coeff(5, -1) == 0

    def test_binomial_prob_basic(self) -> None:
        """Test binomial probability computation."""
        # P(X=2 | n=5, p=0.5) = C(5,2) * 0.5^2 * 0.5^3 = 10 * 0.03125 = 0.3125
        prob = binomial_prob(5, 2, 0.5)
        assert abs(prob - 0.3125) < 1e-6

        # Edge cases
        assert binomial_prob(5, 0, 0.0) == 1.0
        assert binomial_prob(5, 5, 1.0) == 1.0
        assert binomial_prob(5, 3, 0.0) == 0.0

    def test_sojourn_density_neutral(self) -> None:
        """Test sojourn density for neutral case (S=0)."""
        # For S=0: density(x) = 2 / (x * (1-x))
        x = 0.5
        expected = 2.0 / (0.5 * 0.5)  # = 8
        density = sojourn_density(x, 0.0)
        assert abs(density - expected) < 1e-6

    def test_sojourn_density_selected(self) -> None:
        """Test sojourn density for selected case."""
        # For S != 0, density should be positive for 0 < x < 1
        x = 0.5
        for S in [-10, -1, 1, 10]:
            density = sojourn_density(x, S)
            assert density > 0

    def test_sojourn_density_boundaries(self) -> None:
        """Test sojourn density at boundaries."""
        assert sojourn_density(0.0, 1.0) == 0.0
        assert sojourn_density(1.0, 1.0) == 0.0

    def test_fixation_probability_neutral(self) -> None:
        """Test fixation probability for neutral case."""
        # For S=0, relative fixation probability = 1
        assert fixation_probability(0.0) == 1.0

    def test_fixation_probability_beneficial(self) -> None:
        """Test fixation probability for beneficial mutations."""
        # Beneficial mutations (S > 0) have higher fixation probability
        fix_neutral = fixation_probability(0.0)
        fix_ben = fixation_probability(10.0)
        assert fix_ben > fix_neutral

    def test_fixation_probability_deleterious(self) -> None:
        """Test fixation probability for deleterious mutations."""
        # Deleterious mutations (S < 0) have lower fixation probability
        fix_neutral = fixation_probability(0.0)
        fix_del = fixation_probability(-10.0)
        assert fix_del < fix_neutral

    def test_gamma_density_basic(self) -> None:
        """Test gamma density function."""
        # Should be 0 for positive s (deleterious only)
        assert gamma_density(1.0, shape=2.0, mean=10.0) == 0.0

        # Should be positive for negative s
        density = gamma_density(-5.0, shape=2.0, mean=10.0)
        assert density > 0

    def test_exponential_density_basic(self) -> None:
        """Test exponential density function."""
        # Should be 0 for negative s (beneficial only)
        assert exponential_density(-1.0, mean=5.0) == 0.0

        # Should be positive for positive s
        density = exponential_density(5.0, mean=10.0)
        assert density > 0

        # Mean of exponential = 1/rate, so at x=mean, density = 1/mean * exp(-1)
        density_at_mean = exponential_density(10.0, mean=10.0)
        expected = (1.0 / 10.0) * np.exp(-1.0)
        assert abs(density_at_mean - expected) < 1e-6

    def test_gamma_expo_density(self) -> None:
        """Test combined gamma-exponential density."""
        # Deleterious part (s < 0)
        del_density = gamma_expo_density(
            -5.0, shape_del=2.0, mean_del=10.0, prop_ben=0.1, mean_ben=5.0
        )
        assert del_density > 0

        # Beneficial part (s > 0)
        ben_density = gamma_expo_density(
            5.0, shape_del=2.0, mean_del=10.0, prop_ben=0.1, mean_ben=5.0
        )
        assert ben_density > 0

        # At s=0, density should be 0
        zero_density = gamma_expo_density(
            0.0, shape_del=2.0, mean_del=10.0, prop_ben=0.1, mean_ben=5.0
        )
        assert zero_density == 0.0

    def test_gamma_gamma_density(self) -> None:
        """Test gamma-gamma density."""
        # Should have positive density for both negative and positive s
        del_density = gamma_gamma_density(
            -5.0, shape_del=2.0, mean_del=10.0,
            prop_ben=0.1, shape_ben=1.0, mean_ben=5.0
        )
        assert del_density > 0

        ben_density = gamma_gamma_density(
            5.0, shape_del=2.0, mean_del=10.0,
            prop_ben=0.1, shape_ben=1.0, mean_ben=5.0
        )
        assert ben_density > 0

    def test_displaced_gamma_density(self) -> None:
        """Test displaced gamma density."""
        # With displacement=-10, s=-5 becomes shifted_s=-5-(-10)=-5+10=5
        # Since shifted_s > 0, this should return 0 (gamma is for negative shifted values)
        density1 = displaced_gamma_density(-5.0, shape=2.0, mean=10.0, displacement=-10.0)
        assert density1 == 0.0

        # With displacement=-20, s=-5 becomes shifted_s=-5-(-20)=15
        # Still positive, so 0
        density2 = displaced_gamma_density(-5.0, shape=2.0, mean=10.0, displacement=-20.0)
        assert density2 == 0.0

        # With displacement=0, s=-5 becomes shifted_s=-5
        # Negative shifted_s, so should have positive density
        density3 = displaced_gamma_density(-5.0, shape=2.0, mean=10.0, displacement=0.0)
        assert density3 > 0

    def test_fold_sfs(self) -> None:
        """Test SFS folding."""
        # For n=6, unfolded SFS has 5 bins (j=1,2,3,4,5)
        # Folded SFS has 3 bins
        # j=1 combines with j=5
        # j=2 combines with j=4
        # j=3 stays alone
        unfolded = np.array([10, 20, 15, 8, 5])
        folded = fold_sfs(unfolded)

        assert len(folded) == 3
        assert folded[0] == 10 + 5  # j=1 + j=5
        assert folded[1] == 20 + 8  # j=2 + j=4
        assert folded[2] == 15       # j=3

    def test_precalculate_s_grid(self) -> None:
        """Test S grid generation matches GRAPES structure."""
        s_grid = precalculate_s_grid(n_points=100)

        # Should have negative and positive values
        assert np.any(s_grid < 0)
        assert np.any(s_grid > 0)

        # Should have values close to 0 (within the medium region [-0.5, 0.5])
        assert np.any((s_grid > -0.5) & (s_grid < 0.5))

        # Should be sorted
        assert np.all(np.diff(s_grid) >= 0)

        # GRAPES uses 4000 points total (ignores n_points parameter now)
        assert len(s_grid) == 4000

        # Check bounds match GRAPES
        assert s_grid.min() > -500.1
        assert s_grid.max() < 10000.1

    def test_poisson_log_likelihood(self) -> None:
        """Test Poisson log-likelihood computation."""
        observed = np.array([5.0, 10.0, 3.0])
        expected = np.array([5.0, 10.0, 3.0])

        # Perfect match should give best likelihood
        ll_perfect = poisson_log_likelihood(observed, expected)

        # Mismatch should give worse likelihood
        expected_bad = np.array([10.0, 5.0, 6.0])
        ll_bad = poisson_log_likelihood(observed, expected_bad)

        assert ll_perfect > ll_bad


class TestDFEModels:
    """Tests for DFE model classes."""

    def test_model_registry(self) -> None:
        """Test that all models are registered."""
        assert "GammaZero" in MODELS
        assert "GammaExpo" in MODELS
        assert "GammaGamma" in MODELS
        assert "DisplacedGamma" in MODELS

    def test_get_model_valid(self) -> None:
        """Test getting valid models."""
        for name in MODELS:
            model = get_model(name)
            assert isinstance(model, DFEModel)
            assert model.name == name

    def test_get_model_invalid(self) -> None:
        """Test getting invalid model raises error."""
        with pytest.raises(ValueError, match="Unknown model"):
            get_model("InvalidModel")

    def test_gamma_zero_model(self) -> None:
        """Test GammaZero model properties."""
        model = GammaZeroModel()
        assert model.name == "GammaZero"
        assert model.n_params == 2
        assert len(model.param_names) == 2
        assert len(model.param_bounds) == 2

        # Test density
        params = np.array([2.0, 10.0])  # shape=2, mean=10
        assert model.density(-5.0, params) > 0
        assert model.density(5.0, params) == 0  # No beneficial

    def test_gamma_expo_model(self) -> None:
        """Test GammaExpo model properties."""
        model = GammaExpoModel()
        assert model.name == "GammaExpo"
        assert model.n_params == 4
        assert "shape_del" in model.param_names
        assert "prop_ben" in model.param_names

        # Test density
        params = np.array([2.0, 10.0, 0.1, 5.0])
        assert model.density(-5.0, params) > 0
        assert model.density(5.0, params) > 0

    def test_gamma_gamma_model(self) -> None:
        """Test GammaGamma model properties."""
        model = GammaGammaModel()
        assert model.name == "GammaGamma"
        assert model.n_params == 5
        assert "shape_ben" in model.param_names

    def test_displaced_gamma_model(self) -> None:
        """Test DisplacedGamma model properties."""
        model = DisplacedGammaModel()
        assert model.name == "DisplacedGamma"
        assert model.n_params == 3
        assert "displacement" in model.param_names


class TestDFEInput:
    """Tests for DFEInput dataclass."""

    def test_creation(self) -> None:
        """Test creating DFEInput."""
        dfe_input = DFEInput(
            sfs_neutral=np.array([10, 8, 5, 3]),
            sfs_selected=np.array([8, 6, 4, 2]),
            divergence_neutral=100,
            divergence_selected=80,
            n_samples=10,
        )

        assert len(dfe_input.sfs_neutral) == 4
        assert dfe_input.divergence_neutral == 100
        assert dfe_input.n_samples == 10


class TestDFEResult:
    """Tests for DFEResult dataclass."""

    def test_creation(self) -> None:
        """Test creating DFEResult."""
        result = DFEResult(
            model="GammaExpo",
            alpha=0.3,
            alpha_down=0.2,
            alpha_up=0.4,
            omega_a=0.1,
            omega_na=0.2,
            dfe_params={"shape_del": 2.0, "mean_del": 100.0},
            log_likelihood=-150.0,
            aic=310.0,
            converged=True,
        )

        assert result.alpha == 0.3
        assert result.converged is True

    def test_str_representation(self) -> None:
        """Test string representation."""
        result = DFEResult(
            model="GammaExpo",
            alpha=0.3,
            alpha_down=0.2,
            alpha_up=0.4,
            omega_a=0.1,
            omega_na=0.2,
            dfe_params={"shape_del": 2.0},
            log_likelihood=-150.0,
            aic=310.0,
            converged=True,
        )

        s = str(result)
        assert "GammaExpo" in s
        assert "0.3" in s or "0.30" in s
        assert "Alpha" in s

    def test_to_dict(self) -> None:
        """Test dictionary conversion."""
        result = DFEResult(
            model="GammaExpo",
            alpha=0.3,
            alpha_down=0.2,
            alpha_up=0.4,
            omega_a=0.1,
            omega_na=0.2,
            dfe_params={"shape_del": 2.0},
            log_likelihood=-150.0,
            aic=310.0,
            converged=True,
        )

        d = result.to_dict()
        assert d["model"] == "GammaExpo"
        assert d["alpha"] == 0.3
        assert d["converged"] is True


class TestPrecomputedData:
    """Tests for PrecomputedData class."""

    def test_precomputation(self) -> None:
        """Test pre-computation for a sample size."""
        precalc = PrecomputedData(n_samples=10, n_s_points=100)

        assert precalc.n_samples == 10
        assert len(precalc.s_grid) > 0
        assert precalc.expected_counts.shape[0] == len(precalc.s_grid)
        assert precalc.expected_counts.shape[1] == 9  # n-1 frequency classes


class TestPolymorphismDataConversion:
    """Tests for converting PolymorphismData to DFE input."""

    def test_conversion_basic(self) -> None:
        """Test basic conversion."""
        poly_data = PolymorphismData(
            polymorphisms=[(0.2, "N"), (0.4, "S"), (0.6, "N"), (0.8, "S")],
            dn=10,
            ds=15,
        )

        dfe_input = polymorphism_data_to_dfe_input(poly_data, n_samples=10)

        assert dfe_input.divergence_selected == 10
        assert dfe_input.divergence_neutral == 15
        assert dfe_input.n_samples == 10
        assert len(dfe_input.sfs_neutral) == 5  # floor(10/2)
        assert len(dfe_input.sfs_selected) == 5


class TestDFEAlpha:
    """Integration tests for dfe_alpha function."""

    def test_dfe_alpha_basic(self, tmp_path: Path) -> None:
        """Test basic DFE alpha estimation."""
        # Create test data with some variation
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
>seq6
ATGCTGATGCTGATGATG
>seq7
ATGATGATGATGATGCTG
>seq8
ATGATGCTGATGATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = dfe_alpha(ingroup_fa, outgroup_fa, model="GammaZero")

        assert isinstance(result, DFEResult)
        assert result.model == "GammaZero"
        assert 0 <= result.alpha <= 1 or result.alpha < 0  # Alpha can be negative

    def test_dfe_alpha_all_models(self, tmp_path: Path) -> None:
        """Test that all models can run without error."""
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
>seq6
ATGCTGATGCTGATGATG
>seq7
ATGATGATGATGATGCTG
>seq8
ATGATGCTGATGATGATG
>seq9
ATGATGATGATGATGATG
>seq10
ATGCTGATGATGCTGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        for model_name in AVAILABLE_MODELS:
            result = dfe_alpha(ingroup_fa, outgroup_fa, model=model_name)
            assert isinstance(result, DFEResult)
            assert result.model == model_name


class TestDFEAlphaFromSFS:
    """Tests for dfe_alpha_from_sfs function."""

    def test_from_sfs_basic(self) -> None:
        """Test alpha estimation from pre-computed SFS."""
        # Create synthetic SFS data
        sfs_neutral = np.array([20, 15, 10, 8, 5])
        sfs_selected = np.array([15, 12, 8, 6, 4])

        result = dfe_alpha_from_sfs(
            sfs_neutral=sfs_neutral,
            sfs_selected=sfs_selected,
            dn=50,
            ds=80,
            n_samples=12,
            model="GammaZero",
        )

        assert isinstance(result, DFEResult)
        assert result.model == "GammaZero"


class TestCompareModels:
    """Tests for compare_models function."""

    def test_compare_basic(self, tmp_path: Path) -> None:
        """Test model comparison."""
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
>seq6
ATGCTGATGCTGATGATG
>seq7
ATGATGATGATGATGCTG
>seq8
ATGATGCTGATGATGATG
>seq9
ATGATGATGATGATGATG
>seq10
ATGCTGATGATGCTGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        results = compare_models(
            ingroup_fa, outgroup_fa,
            models=["GammaZero", "GammaExpo"]
        )

        assert len(results) == 2
        # Results should be sorted by AIC
        assert results[0].aic <= results[1].aic


class TestDFEAlphaAggregated:
    """Tests for dfe_alpha_aggregated function."""

    def test_aggregated_basic(self) -> None:
        """Test aggregated DFE analysis."""
        # Create synthetic gene data
        gene_data = []
        for i in range(5):
            poly = []
            for freq in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
                poly.extend([(freq, "N")] * 2)
                poly.extend([(freq, "S")] * 3)
            gene_data.append(
                PolymorphismData(
                    polymorphisms=poly,
                    dn=10,
                    ds=15,
                    gene_id=f"gene{i}",
                )
            )

        result = dfe_alpha_aggregated(
            gene_data,
            n_samples=20,
            model="GammaZero"
        )

        assert isinstance(result, DFEResult)
        assert result.model == "GammaZero"
