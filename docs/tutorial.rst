Getting Started Tutorial
========================

This tutorial walks through using MKado for McDonald-Kreitman analysis, from running your first test to interpreting results.

Prerequisites
-------------

1. Install MKado (see :doc:`installation`)
2. Have aligned coding sequences in FASTA format

The McDonald-Kreitman Test
--------------------------

The MK test (`McDonald & Kreitman 1991`_) compares the ratio of non-synonymous to synonymous changes within species (polymorphism) versus between species (divergence). Under neutral evolution, these ratios should be equal.

The test produces a 2x2 contingency table:

============  ===============  =============
              Non-synonymous   Synonymous
============  ===============  =============
Divergence    Dn               Ds
Polymorphism  Pn               Ps
============  ===============  =============

Your First MK Test
------------------

Let's run a basic MK test using the example data.

**Step 1: Examine your data**

.. code-block:: bash

   # Get information about an alignment file
   mkado info examples/anopheles_batch/AGAP000010.fa

This shows the number of sequences, sequence lengths, and sequence names.

**Step 2: Run the MK test**

.. code-block:: bash

   # Standard MK test: gamb (ingroup) vs afun (outgroup)
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun

**Step 3: Interpret the output**

The output shows:

- **Dn, Ds**: Fixed non-synonymous and synonymous differences between species
- **Pn, Ps**: Polymorphic non-synonymous and synonymous sites within species
- **Fisher's p-value**: Statistical significance of the departure from neutrality
- **Neutrality Index (NI)**: (Pn/Ps) / (Dn/Ds) - ratio of ratios
- **Alpha**: Proportion of substitutions driven by positive selection

Interpreting Results
--------------------

Neutrality Index (NI)
^^^^^^^^^^^^^^^^^^^^^

The Neutrality Index (`Rand & Kann 1996`_) measures the direction and degree of departure from neutral evolution:

- **NI = 1**: Consistent with neutral evolution
- **NI > 1**: Excess polymorphism (segregating weakly deleterious variants)
- **NI < 1**: Excess divergence (positive selection)

Alpha
^^^^^

Alpha, the proportion of adaptive substitutions (`Smith & Eyre-Walker 2002`_):

- **Alpha = 0**: No adaptive substitutions
- **Alpha > 0**: Proportion of fixed differences due to positive selection
- **Alpha < 0**: Excess polymorphism relative to divergence

Asymptotic MK Test
------------------

The standard MK test can be biased by slightly deleterious mutations that segregate as polymorphisms but rarely reach fixation. These inflate Pn relative to Dn, causing alpha to be underestimated. The asymptotic MK test (`Messer & Petrov 2013`_) addresses this by examining how alpha varies with **derived allele frequency** in the ingroup.

The key insight is that weakly deleterious mutations are more common at low frequencies (where they haven't yet been purged by selection) and rare at high frequencies. By extrapolating alpha to a derived frequency of 1.0, we estimate what alpha would be if all deleterious polymorphisms were removed.

Determining Derived Allele Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To calculate derived allele frequency, we need to identify which allele is **ancestral** (the original state) versus **derived** (the new mutation). MKado uses the outgroup to polarize mutations:

1. At each polymorphic codon position in the ingroup, examine which codons are present
2. Compare with the outgroup codon at the same position
3. The allele shared between ingroup and outgroup is assumed to be **ancestral**
4. **Derived frequency** = 1.0 − (frequency of ancestral allele in ingroup)

For example, if a codon position has allele A at 80% and allele B at 20% in the ingroup, and the outgroup has allele A, then A is ancestral and the derived allele frequency is 0.20.

.. note::

   Sites where no ingroup allele matches the outgroup cannot be polarized and are excluded from the frequency spectrum analysis.

The Asymptotic Method
^^^^^^^^^^^^^^^^^^^^^

Once derived frequencies are computed, the test proceeds as follows:

1. **Bin polymorphisms by derived frequency** — Group Pn and Ps into frequency bins (e.g., 0.0–0.05, 0.05–0.10, etc.)
2. **Calculate alpha at each bin** — α(x) = 1 − (Ds/Dn) × (Pn(x)/Ps(x))
3. **Fit a curve** — Fit a model to the per-bin alpha estimates (see Model Selection below)
4. **Extrapolate to x = 1** — The asymptotic alpha is the fitted value at derived frequency = 1.0

Model Selection
^^^^^^^^^^^^^^^

MKado fits two candidate models to the α(x) curve and automatically selects the best one:

- **Exponential model**: α(x) = a + b·exp(−c·x) — 3 parameters
- **Linear model**: α(x) = a + b·x — 2 parameters

The exponential model captures the expected shape when weakly deleterious mutations cause alpha to increase with frequency. However, when data are noisy or sparse, the exponential fit can be unstable.

Model selection follows the `asymptoticMK <https://github.com/MesserLab/asymptoticMK>`_ R package convention (`Haller & Messer 2017`_):

1. If the exponential fit produces a confidence interval width > 100, it is considered unstable and the **linear model** is used
2. Otherwise, both models are compared using **AIC** (Akaike Information Criterion), which balances fit quality against model complexity
3. The model with the lower AIC is selected
4. If both fits fail, the highest-frequency bin's alpha value is reported as a fallback

The output reports which model was selected (exponential or linear) along with the fitted parameters.

.. code-block:: bash

   # Run asymptotic MK test
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a

   # Customize number of frequency bins
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a -b 20

The output includes:

- **Alpha asymptotic**: Extrapolated alpha value at derived frequency = 1
- **95% CI**: Confidence interval for alpha asymptotic
- **Model type**: Whether exponential or linear model was selected
- **Per-bin alpha values**: Alpha estimates at each frequency class

Polarized MK Test
-----------------

The polarized MK test uses a second outgroup to determine the direction of mutations:

.. code-block:: bash

   # Use amin as second outgroup for polarization
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun --polarize-match amin

This allows you to distinguish:

- Mutations that occurred on the ingroup lineage
- Mutations that occurred on the outgroup lineage

Frequency Filtering Options
---------------------------

MKado provides two different frequency-related options that serve distinct purposes:

``--min-freq`` (Standard and Polarized MK tests)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option **excludes** rare polymorphisms below a frequency threshold from the Pn/Ps counts.

- **Applies to**: Standard MK test, Polarized MK test
- **Purpose**: Filter out singletons or very rare variants that may be sequencing errors or very recent mutations
- **How it works**: Polymorphisms with derived allele frequency < ``min_freq`` are not counted

.. code-block:: bash

   # Exclude polymorphisms below 5% frequency
   mkado test alignment.fa -i dmel -o dsim --min-freq 0.05

``--no-singletons`` (Standard and Polarized MK tests)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A convenience option that automatically sets ``--min-freq`` to exclude singletons (variants appearing only once).

- **Applies to**: Standard MK test, Polarized MK test
- **How it works**: Calculates ``1/n`` where n is the sample size, and uses that as the minimum frequency threshold
- **Sample size**: Uses ingroup count, or ingroup + outgroup if ``--pool-polymorphisms`` is enabled

.. code-block:: bash

   # Exclude singletons automatically
   mkado test alignment.fa -i dmel -o dsim --no-singletons

   # For batch processing (threshold calculated per gene)
   mkado batch alignments/ -i dmel -o dsim --no-singletons

.. note::

   ``--no-singletons`` cannot be used with ``--min-freq`` (they are mutually exclusive), ``--asymptotic``, or ``--alpha-tg``.

``--freq-cutoffs`` (Asymptotic MK test)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option defines the frequency range used for **curve fitting** in the asymptotic test — it does not exclude data.

- **Applies to**: Aggregated asymptotic MK test (``mkado batch -a``)
- **Default**: ``0.1,0.9``
- **Purpose**: Avoid fitting to extreme frequency bins where data may be sparse or noisy
- **How it works**:

  1. All polymorphisms are counted and binned by derived frequency
  2. Only bins within the ``[low, high]`` range are used to fit the exponential/linear model
  3. The curve is still extrapolated to frequency = 1.0

.. code-block:: bash

   # Fit model using only bins between 15% and 85% frequency
   mkado batch alignments/ -i dmel -o dsim -a --freq-cutoffs 0.15,0.85

Key Differences
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Aspect
     - ``--min-freq``
     - ``--freq-cutoffs``
   * - Effect
     - **Excludes** polymorphisms from counts
     - **Restricts** which bins are used for fitting
   * - Scope
     - Per-polymorphism filtering
     - Curve-fitting range only
   * - Data
     - Polymorphisms below threshold not counted
     - All polymorphisms counted; fitting uses subset
   * - Analysis types
     - Standard MK, Polarized MK
     - Asymptotic MK (aggregated batch mode)

Output Formats
--------------

MKado supports multiple output formats:

.. code-block:: bash

   # Pretty-printed output (default)
   mkado test alignment.fa -i species1 -o species2

   # Tab-separated values
   mkado test alignment.fa -i species1 -o species2 -f tsv

   # JSON format
   mkado test alignment.fa -i species1 -o species2 -f json

Next Steps
----------

- Learn about :doc:`batch-workflow` for processing multiple genes
- Review :doc:`file-formats` for input requirements
- Explore the :doc:`api` for programmatic access

References
----------

.. _Haller & Messer 2017: https://doi.org/10.1534/g3.117.039693
.. _McDonald & Kreitman 1991: https://doi.org/10.1038/351652a0
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Rand & Kann 1996: https://doi.org/10.1093/oxfordjournals.molbev.a025634
.. _Smith & Eyre-Walker 2002: https://doi.org/10.1038/4151022a
