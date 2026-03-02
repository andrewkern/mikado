Imputed MK Test
===============

MKado implements the imputed MK test (impMKT) from `Murga-Moreno et al. (2022)`_, which corrects for slightly deleterious mutations by imputing their count from the synonymous frequency spectrum rather than discarding low-frequency variants.

Background
----------

The standard MK test is biased by **slightly deleterious mutations** (SDMs) that inflate nonsynonymous polymorphism (Pn) without contributing to divergence (Dn). Several approaches exist to address this:

- **Asymptotic MK test**: Fits a curve across the full frequency spectrum and extrapolates to x=1
- **fwwMKT** (frequency-weighted): Discards all low-frequency polymorphisms below a cutoff

The imputed MK test takes a different approach: instead of discarding low-frequency data, it **estimates how many** low-frequency nonsynonymous polymorphisms are slightly deleterious, using the synonymous frequency spectrum as a neutral reference.

The Key Insight
---------------

Under neutrality, nonsynonymous and synonymous polymorphisms should have the same ratio of low-frequency to high-frequency variants. An excess of low-frequency nonsynonymous polymorphisms (relative to synonymous) indicates segregating slightly deleterious mutations.

By imputing this excess, the test:

- Retains more data than methods that discard all low-frequency variants
- Increases statistical power at the gene level
- Decomposes the distribution of fitness effects (DFE) into interpretable fractions

The Imputation Formula
----------------------

Polymorphisms are split at a derived allele frequency (DAF) cutoff (default 15%):

1. Count low-frequency (DAF <= cutoff) and high-frequency (DAF > cutoff) variants separately for nonsynonymous (P\ :sub:`n`) and synonymous (P\ :sub:`s`) classes

2. Compute the neutral ratio from synonymous polymorphisms:

   .. math::

      r = \frac{P_{s,low}}{P_{s,high}}

3. Impute the number of weakly deleterious nonsynonymous polymorphisms:

   .. math::

      P_{wd} = P_{n,low} - P_{n,high} \times r

   This is clamped to >= 0.

4. Compute neutral nonsynonymous polymorphisms:

   .. math::

      P_{n,neutral} = P_n - P_{wd}

5. Calculate corrected alpha:

   .. math::

      \alpha = 1 - \frac{P_{n,neutral}}{P_s} \times \frac{D_s}{D_n}

6. Significance is assessed with Fisher's exact test on the corrected 2x2 table.

DFE Fractions
-------------

When the number of synonymous (m\ :sub:`0`) and nonsynonymous (m\ :sub:`i`) sites are provided, the test decomposes the DFE into four fractions:

- **alpha (a)**: Fraction of adaptive substitutions
- **f**: Fraction of effectively neutral nonsynonymous mutations
- **b**: Fraction of weakly deleterious mutations
- **d**: Fraction of strongly deleterious mutations (d = 1 - f - b)

These fractions sum to 1 and describe the full distribution of fitness effects for new nonsynonymous mutations.

Usage
-----

Single Gene Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Basic imputed MK test (default 15% DAF cutoff)
   mkado test alignment.fa -i ingroup -o outgroup --imputed

   # With custom DAF cutoff (10%)
   mkado test alignment.fa -i ingroup -o outgroup --imputed --min-freq 0.10

   # Separate ingroup/outgroup files
   mkado test ingroup.fa outgroup.fa --imputed

.. note::

   The ``--min-freq`` option is reused as the DAF cutoff when ``--imputed`` is set. If ``--min-freq`` is not specified, the default cutoff of 0.15 (15%) is used.

Batch Analysis (Aggregated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For multi-gene analyses, pooling data across genes increases power:

.. code-block:: bash

   # Pool polymorphisms and divergence across all genes
   mkado batch alignments/ -i ingroup -o outgroup --imputed

Batch Analysis (Per-Gene)
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Run imputed test separately for each gene
   mkado batch alignments/ -i ingroup -o outgroup --imputed --per-gene

Output
------

The imputed test reports:

- **alpha**: Corrected proportion of adaptive substitutions
- **p_value**: Fisher's exact test on the corrected contingency table
- **Pwd**: Imputed count of weakly deleterious nonsynonymous polymorphisms
- **Pn_neutral**: Nonsynonymous polymorphisms after removing imputed SDMs
- **Dn, Ds, Pn, Ps**: Raw counts
- **cutoff**: The DAF cutoff used

Example output (pretty format):

.. code-block:: text

   Imputed MK Test Results:
     Divergence:    Dn=6, Ds=8
     Polymorphism:  Pn=11, Ps=17
     DAF cutoff:    0.15
     Imputed Pwd:   4.82
     Pn (neutral):  6.18
     Alpha:         0.5154
     p-value:       0.0891

Comparison with Other Methods
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Method
     - Approach
     - Best used when
   * - Standard MK
     - No correction
     - Quick assessment; comparing specific genes
   * - Asymptotic MK
     - Curve fitting across frequency spectrum
     - Genome-wide analyses with many polymorphisms
   * - fwwMKT
     - Discards low-frequency variants
     - Simple correction, but loses data
   * - **Imputed MK**
     - Imputes SDM count from synonymous spectrum
     - Gene-level analyses; maximizing statistical power

The imputed test is particularly useful when:

- You want gene-level significance (p-values) rather than only genome-wide estimates
- You want to retain as much data as possible
- You want to decompose the DFE into interpretable fractions

When to Use the Imputed MK Test
-------------------------------

**Use imputed MK when:**

- Analyzing individual genes or small gene sets where power matters
- You want per-gene corrected alpha with p-values
- You want DFE decomposition (with site count information)
- The asymptotic test lacks sufficient data for curve fitting

**Consider alternatives when:**

- You have genome-wide data with thousands of polymorphisms (asymptotic MK may be more robust)
- You want frequency-bin visualization of alpha(x) (use asymptotic MK with ``--plot-asymptotic``)
- You need unbiased multi-gene weighting without frequency modeling (use alpha_TG)

Interpreting Results
--------------------

- **alpha > 0**: Evidence for positive selection (proportion of adaptive substitutions)
- **alpha ~ 0**: Consistent with neutral evolution
- **alpha < 0**: Unusual; may indicate model misspecification
- **Pwd > 0**: Low-frequency nonsynonymous excess detected and corrected
- **Pwd = 0**: No evidence of segregating slightly deleterious mutations (or all polymorphisms are high-frequency)

Choosing the DAF Cutoff
-----------------------

The default cutoff of 15% follows the recommendation of `Murga-Moreno et al. (2022)`_. The cutoff defines the boundary between "low-frequency" and "high-frequency" polymorphisms:

- **Lower cutoff** (e.g., 5%): Only the rarest variants are considered low-frequency. More conservative imputation.
- **Higher cutoff** (e.g., 25%): More variants classified as low-frequency. More aggressive imputation.

The optimal cutoff depends on the effective population size and strength of selection in your system.

Reference
---------

.. _Murga-Moreno et al. (2022): https://doi.org/10.1093/g3journal/jkac206

Murga-Moreno J, Coronado-Zamora M, Casillas S, Barbadilla A (2022) impMKT: the imputed McDonald and Kreitman test, a straightforward correction that significantly increases the evidence of positive selection of the McDonald and Kreitman test at the gene level. *G3: Genes, Genomes, Genetics* 12(10):jkac206. https://doi.org/10.1093/g3journal/jkac206
