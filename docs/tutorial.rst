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

The standard MK test can be biased by slightly deleterious mutations. The asymptotic MK test (`Messer & Petrov 2013`_) addresses this by:

1. Binning polymorphisms by frequency
2. Estimating alpha at each frequency bin
3. Fitting an exponential curve
4. Extrapolating to asymptotic alpha at frequency = 1

.. code-block:: bash

   # Run asymptotic MK test
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a

   # Customize number of frequency bins
   mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a -b 20

The output includes:

- **Alpha asymptotic**: Extrapolated alpha value
- **95% CI**: Confidence interval for alpha asymptotic
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

.. _McDonald & Kreitman 1991: https://doi.org/10.1038/351652a0
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Rand & Kann 1996: https://doi.org/10.1093/oxfordjournals.molbev.a025634
.. _Smith & Eyre-Walker 2002: https://doi.org/10.1038/4151022a
