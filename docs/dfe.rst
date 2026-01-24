DFE-Based Alpha Estimation
==========================

MKado implements Distribution of Fitness Effects (DFE) models for estimating alpha, following the GRAPES methodology (`Galtier 2016`_). This approach explicitly models the DFE of new mutations and uses it to predict the expected rate of non-adaptive substitutions.

Background
----------

The standard MK test and even the asymptotic MK test make simplifying assumptions about the distribution of fitness effects. DFE-based methods take a more principled approach:

1. **Model the DFE explicitly**: Assume new mutations follow a specific fitness distribution (e.g., gamma for deleterious, exponential for beneficial)
2. **Predict polymorphism**: Use population genetics theory to predict the expected site frequency spectrum (SFS) under this DFE
3. **Fit to data**: Find DFE parameters that best explain the observed SFS
4. **Compute alpha**: Use the fitted DFE to predict non-adaptive divergence and calculate what fraction of observed divergence is adaptive

This approach has been shown to outperform both standard and asymptotic MK methods, particularly when the true DFE is complex (`Tataru et al. 2017`_).

The GRAPES Framework
--------------------

GRAPES (Galtier 2016) introduced a unified framework for DFE-based inference that mkado implements. The key components are:

**Site Frequency Spectrum (SFS)**

The SFS describes how many polymorphisms are observed at each frequency class. Under neutrality, the expected count at frequency *i/n* is proportional to 1/i. Selection distorts this pattern:

- Deleterious mutations are enriched at low frequencies (selection removes them before they spread)
- Beneficial mutations are depleted at intermediate frequencies (they sweep rapidly to fixation)

**Fixation Probability**

The probability that a new mutation reaches fixation depends on its selection coefficient *S* = 4Ns (scaled by effective population size):

.. math::

   P_{fix}(S) = \frac{1 - e^{-2S/N}}{1 - e^{-2S}}

For neutral mutations (S=0), P\ :sub:`fix` = 1/2N. Beneficial mutations (S>0) have higher fixation probability; deleterious mutations (S<0) have lower.

**Likelihood Function**

GRAPES uses a Poisson likelihood for the SFS data:

.. math::

   \mathcal{L} = \prod_i \text{Poisson}(O_i | E_i(\theta))

Where O\ :sub:`i` is the observed count in frequency bin *i*, E\ :sub:`i`\ (θ) is the expected count given DFE parameters θ, and the product is over all frequency bins for both synonymous (neutral) and nonsynonymous (selected) sites.

Available Models
----------------

MKado implements four DFE models from GRAPES:

GammaZero
^^^^^^^^^

The simplest model assuming all nonsynonymous mutations are deleterious or neutral:

- **Deleterious DFE**: Gamma distribution with shape β and mean |S|
- **No beneficial mutations**

Parameters:

- ``shape``: Gamma shape parameter (β)
- ``mean``: Mean |S| of the deleterious gamma distribution

This model has 2 free parameters. Use it when positive selection is expected to be negligible.

GammaExpo
^^^^^^^^^

Extends GammaZero by adding a class of beneficial mutations:

- **Deleterious DFE**: Gamma distribution (shape β\ :sub:`d`, mean |S|\ :sub:`d`)
- **Beneficial DFE**: Exponential distribution (mean S\ :sub:`b`) with proportion p\ :sub:`b`

Parameters:

- ``shape_del``: Gamma shape for deleterious mutations
- ``mean_del``: Mean |S| for deleterious mutations
- ``prop_ben``: Proportion of new mutations that are beneficial
- ``mean_ben``: Mean S for beneficial mutations

This model has 4 free parameters. It's the recommended default for most analyses as it can capture both purifying and positive selection.

GammaGamma
^^^^^^^^^^

Uses gamma distributions for both deleterious and beneficial effects:

- **Deleterious DFE**: Gamma(β\ :sub:`d`, |S|\ :sub:`d`)
- **Beneficial DFE**: Gamma(β\ :sub:`b`, S\ :sub:`b`) with proportion p\ :sub:`b`

This model has 5 free parameters. Use it when you expect the beneficial DFE to have a shape different from exponential.

DisplacedGamma
^^^^^^^^^^^^^^

A single displaced gamma distribution spanning both deleterious and beneficial effects:

- **DFE**: Gamma(β, mean) shifted by displacement S\ :sub:`0`

Parameters:

- ``shape``: Gamma shape
- ``mean``: Gamma mean
- ``displacement``: Shift parameter (positive values shift toward beneficial)

This model has 3 free parameters. It's useful when the DFE is expected to be continuous across the neutral boundary.

Computing Alpha
---------------

Once the DFE is fitted, alpha is computed as:

.. math::

   \alpha = 1 - \frac{\omega_{na}}{\omega_{obs}}

Where:

- ω\ :sub:`obs` = D\ :sub:`n`/D\ :sub:`s` (observed dN/dS ratio)
- ω\ :sub:`na` = expected dN/dS from non-adaptive substitutions (integrated over the DFE for S ≤ threshold)

The threshold for "non-adaptive" is typically S < some small positive value (e.g., S < 1), meaning mutations that fix primarily by drift rather than positive selection.

Usage
-----

Python API
^^^^^^^^^^

.. code-block:: python

   from mkado.analysis.dfe import dfe_alpha, compare_models

   # Basic usage with GammaExpo model (recommended)
   result = dfe_alpha("ingroup.fa", "outgroup.fa", model="GammaExpo")

   print(f"Alpha: {result.alpha:.4f}")
   print(f"95% CI: [{result.alpha_down:.4f}, {result.alpha_up:.4f}]")
   print(f"omega_a: {result.omega_a:.4f}")
   print(f"Log-likelihood: {result.log_likelihood:.2f}")
   print(f"AIC: {result.aic:.2f}")

   # Access fitted DFE parameters
   for param, value in result.dfe_params.items():
       print(f"  {param}: {value:.4f}")

   # Compare multiple models
   results = compare_models("ingroup.fa", "outgroup.fa")
   for r in results:
       print(f"{r.model}: alpha={r.alpha:.4f}, AIC={r.aic:.2f}")

Working with Pre-computed SFS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have SFS data from another source:

.. code-block:: python

   import numpy as np
   from mkado.analysis.dfe import dfe_alpha_from_sfs

   # Folded SFS data (counts at each frequency class)
   sfs_neutral = np.array([500., 200., 100., 80., 60., 50.])  # synonymous
   sfs_selected = np.array([100., 30., 15., 10., 8., 5.])     # nonsynonymous

   result = dfe_alpha_from_sfs(
       sfs_neutral=sfs_neutral,
       sfs_selected=sfs_selected,
       dn=200,           # nonsynonymous divergence
       ds=500,           # synonymous divergence
       n_samples=12,     # number of chromosomes sampled
       model="GammaExpo"
   )

Aggregated Multi-Gene Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For genome-wide analyses, aggregate SFS across genes:

.. code-block:: python

   from mkado.analysis.dfe import dfe_alpha_aggregated
   from mkado.analysis.asymptotic import extract_polymorphism_data
   from mkado.core.sequences import SequenceSet
   from pathlib import Path

   # Collect polymorphism data from multiple genes
   gene_data = []
   for fasta in Path("alignments/").glob("*.fa"):
       seqs = SequenceSet.from_fasta(fasta)
       ingroup = seqs.filter_by_name("species1")
       outgroup = seqs.filter_by_name("species2")

       if len(ingroup) >= 4 and len(outgroup) >= 1:
           poly_data = extract_polymorphism_data(ingroup, outgroup)
           gene_data.append(poly_data)

   # Fit DFE to aggregated data
   result = dfe_alpha_aggregated(
       gene_data,
       n_samples=20,      # must be consistent across genes
       model="GammaExpo"
   )

Output
------

The ``DFEResult`` object contains:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Description
   * - ``alpha``
     - Proportion of adaptive amino acid substitutions
   * - ``alpha_down``
     - Lower bound of 95% confidence interval
   * - ``alpha_up``
     - Upper bound of 95% confidence interval
   * - ``omega_a``
     - Adaptive substitution rate (ω\ :sub:`a`)
   * - ``omega_na``
     - Non-adaptive substitution rate (ω\ :sub:`na`)
   * - ``dfe_params``
     - Dictionary of fitted DFE parameters
   * - ``log_likelihood``
     - Log-likelihood of the fitted model
   * - ``aic``
     - Akaike Information Criterion (for model comparison)
   * - ``model``
     - Name of the fitted model

Model Selection
---------------

Use AIC to compare models:

.. code-block:: python

   results = compare_models("ingroup.fa", "outgroup.fa")

   # Results are sorted by AIC (best first)
   best = results[0]
   print(f"Best model: {best.model} (AIC={best.aic:.2f})")

General guidelines:

- **GammaZero**: Use when positive selection is not expected (e.g., highly constrained genes)
- **GammaExpo**: Recommended default; captures both purifying and positive selection
- **GammaGamma**: Use when beneficial mutations may have a non-exponential distribution
- **DisplacedGamma**: Use when the DFE is expected to be unimodal and continuous

Comparison with Other Methods
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Method
     - Corrects for
     - Best used when
   * - Standard α
     - Nothing
     - Quick assessment; single genes
   * - Asymptotic α
     - Weakly deleterious mutations
     - When frequency data is available
   * - **DFE-based α**
     - Full DFE shape
     - Genome-wide analysis (most accurate)

**Key advantages of DFE-based methods:**

- Explicit modeling of the full fitness effect distribution
- Can estimate the DFE itself, not just alpha
- More statistically principled than curve-fitting approaches
- Validated against simulations (`Tataru et al. 2017`_)

**Limitations:**

- Requires sufficient polymorphism data for reliable SFS estimation
- Assumes a specific parametric form for the DFE
- Computationally more intensive than simpler methods

Validation Against GRAPES
-------------------------

MKado's DFE implementation has been validated against the reference GRAPES C++ implementation. On identical input data, mkado produces alpha estimates within 0.001 of GRAPES while being 15-35x faster.

A validation script is provided at ``scripts/validate_grapes_comparison.py``. See the script docstring for instructions on building GRAPES and running the comparison.

References
----------

.. _Galtier 2016: https://doi.org/10.1371/journal.pgen.1005774
.. _Tataru et al. 2017: https://doi.org/10.1534/genetics.117.300323

Galtier N (2016) Adaptive Protein Evolution in Animals and the Effective Population Size Hypothesis. *PLoS Genetics* 12(1):e1005774. https://doi.org/10.1371/journal.pgen.1005774

Tataru P, Mollion M, Glémin S, Bataillon T (2017) Inference of Distribution of Fitness Effects and Proportion of Adaptive Substitutions from Polymorphism Data. *Genetics* 207(3):1103-1119. https://doi.org/10.1534/genetics.117.300323
