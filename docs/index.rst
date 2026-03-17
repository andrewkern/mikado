MKado Documentation
===================

**MKado** is a modern Python implementation of the McDonald-Kreitman test toolkit for detecting selection in molecular evolution.

.. note::

   MKado requires Python 3.12 or later.

Features
--------

- **Standard MK test**: Classic 2x2 contingency table with Fisher's exact test
- **Polarized MK test**: Uses a third outgroup to assign mutations to lineages
- **Asymptotic MK test**: Frequency-bin alpha estimates with exponential extrapolation (`Messer & Petrov 2013`_)
- **Imputed MK test**: Corrects for slightly deleterious mutations by imputation (`Murga-Moreno et al. 2022`_)
- **Tarone-Greenland α_TG**: Weighted multi-gene estimator (`Stoletzki & Eyre-Walker 2011`_)
- **Alternate genetic codes**: Support for 24 NCBI genetic code tables (mitochondrial, plastid, etc.)
- **VCF input**: Go directly from VCF + reference + GFF3 to MK test results
- **Batch processing**: Process multiple genes with parallel execution
- **Volcano plots**: Visualize batch results with publication-ready volcano plots
- **Multiple output formats**: Pretty-print, TSV, and JSON

Quick Example
-------------

.. code-block:: bash

   # Standard MK test
   mkado test alignment.fa -i "dmel" -o "dsim"

   # Asymptotic MK test
   mkado test alignment.fa -i "dmel" -o "dsim" -a

   # Batch process a directory
   mkado batch alignments/ -i "dmel" -o "dsim"

   # Directly from VCF files
   mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   batch-workflow
   vcf-input
   asymptotic
   imputed
   alpha-tg
   dos
   file-formats
   api

Citation
--------

If you use MKado in your research, please cite:

   Rivera-Colón, A. G., Rehmann, C. T., & Kern, A. D. (2026). MKado: a toolkit for McDonald-Kreitman tests of natural selection. *bioRxiv*. https://doi.org/10.64898/2026.03.02.709122

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Murga-Moreno et al. 2022: https://doi.org/10.1093/g3journal/jkac206
.. _Stoletzki & Eyre-Walker 2011: https://doi.org/10.1093/molbev/msq249
