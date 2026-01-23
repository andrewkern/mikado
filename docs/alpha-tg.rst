Tarone-Greenland Alpha (α_TG)
=============================

MKado implements the weighted α_TG estimator from `Stoletzki & Eyre-Walker (2011)`_, which provides an unbiased estimate of the proportion of adaptive substitutions when analyzing multiple genes.

Background
----------

When analyzing many genes, a common approach is to calculate alpha (α = 1 - NI) for each gene and take the mean. However, this simple average is **heavily biased** by genes with small sample sizes, where alpha can take extreme values (e.g., -30 or +5) due to sampling noise.

`Stoletzki & Eyre-Walker (2011)`_ showed that averaging across genes produces biased estimates even with large sample sizes, and introduced a weighted estimator that corrects this problem.

The NI_TG Formula
-----------------

The weighted neutrality index is calculated as:

.. math::

   NI_{TG} = \frac{\sum_i (D_{si} \times P_{ni}) / (P_{si} + D_{si})}{\sum_i (D_{ni} \times P_{si}) / (P_{si} + D_{si})}

Where for each gene *i*:

- D\ :sub:`ni` = nonsynonymous divergence (fixed differences)
- D\ :sub:`si` = synonymous divergence
- P\ :sub:`ni` = nonsynonymous polymorphism
- P\ :sub:`si` = synonymous polymorphism

The weighting by 1/(P\ :sub:`si` + D\ :sub:`si`) downweights genes with small denominators, where estimates are unreliable.

Alpha is then: **α_TG = 1 - NI_TG**

Usage
-----

Use the ``--alpha-tg`` flag with batch processing:

.. code-block:: bash

   # Basic usage
   mkado batch alignments/ -i ingroup -o outgroup --alpha-tg

   # With more bootstrap replicates for tighter CIs
   mkado batch alignments/ -i ingroup -o outgroup --alpha-tg --bootstrap 1000

Example with the included Anopheles data:

.. code-block:: bash

   mkado batch examples/anopheles_batch/ -i gamb -o afun --alpha-tg

Output
------

The output includes:

- **alpha_TG**: Proportion of adaptive substitutions (1 - NI_TG)
- **NI_TG**: The weighted neutrality index
- **CI_low, CI_high**: 95% bootstrap confidence interval
- **num_genes**: Number of genes analyzed
- **Dn, Ds, Pn, Ps**: Total counts across all genes

Example output (TSV format):

.. code-block:: text

   Dn      Ds      Pn    Ps      alpha_TG  NI_TG     CI_low    CI_high   num_genes
   18828   49857   7843  25083   0.022781  0.977219  -0.053529 0.088672  400

Comparison with Other Methods
-----------------------------

Different methods for estimating alpha correct for different biases:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Method
     - Corrects for
     - Best used when
   * - Simple mean α
     - Nothing
     - Never recommended for multi-gene analyses
   * - **α_TG**
     - Sample size heterogeneity
     - Comparing species with little slightly deleterious load
   * - Asymptotic α
     - Slightly deleterious mutations
     - Most genome-wide analyses

**Example comparison** (Anopheles gambiae vs. A. funestus, 400 genes):

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Method
     - Alpha estimate
     - 95% CI
   * - Simple mean
     - -1.19
     - —
   * - α_TG (weighted)
     - +0.02
     - -0.05 to +0.09
   * - Asymptotic α
     - +0.57
     - +0.49 to +0.66

The large gap between α_TG and asymptotic α suggests substantial slightly deleterious polymorphism — a common finding. The asymptotic method extrapolates to high frequencies where deleterious variants have been purged, revealing adaptive substitutions masked by segregating deleterious mutations.

When to Use α_TG
----------------

Use α_TG when:

- You want an unbiased multi-gene estimate without frequency spectrum modeling
- Your species pair has minimal slightly deleterious load
- You want to compare with published NI_TG values

Use asymptotic α (``-a``) when:

- Slightly deleterious mutations are a concern (most cases)
- You have sufficient polymorphism data for frequency binning
- You want the most accurate estimate of adaptive substitution rate

Reference
---------

.. _Stoletzki & Eyre-Walker (2011): https://doi.org/10.1093/molbev/msq249

Stoletzki N, Eyre-Walker A (2011) Estimation of the Neutrality Index. *Molecular Biology and Evolution* 28(1):63-70. https://doi.org/10.1093/molbev/msq249
