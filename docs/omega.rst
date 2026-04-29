Omega Decomposition (ω, ω_a, ω_na)
===================================

MKado reports the rate of nonsynonymous substitution per site (ω), and where
appropriate, its decomposition into adaptive (ω_a) and non-adaptive (ω_na)
components following `Gossmann, Keightley & Eyre-Walker (2012)`_.

Background
----------

The classical McDonald-Kreitman summary statistics (α, NI, DoS) are
*proportions* — they describe what fraction of substitutions were adaptive,
but not the *rate* at which adaptation occurs. As `Gossmann et al. (2012)`_
argue, comparing α between species or genomic regions is misleading when the
underlying neutral substitution rate varies. The authors recommend reporting
ω_a — the rate of adaptive substitution per site, scaled by the synonymous
substitution rate — as a more appropriate cross-comparison statistic.

The decomposition
-----------------

.. math::

   \omega    &= \frac{D_n}{D_s} \cdot \frac{L_s}{L_n} \\
   \omega_a  &= \alpha \cdot \omega \\
   \omega_{na} &= (1 - \alpha) \cdot \omega = \omega - \omega_a

Where:

- :math:`D_n, D_s` — nonsynonymous and synonymous fixed differences
- :math:`L_n, L_s` — nonsynonymous and synonymous site totals (Nei-Gojobori
  averaged across analyzed codons in both groups)
- :math:`\alpha` — proportion of adaptive substitutions

ω = (Dn × Ls) / (Ds × Ln) is the standard dN/dS ratio with Nei-Gojobori
site weighting; ω_a is the adaptive-substitution rate per synonymous site;
ω_na is its complement.

Site counting
-------------

MKado computes :math:`L_n` and :math:`L_s` with the **Nei-Gojobori (1986)**
method: for each clean codon, count the fraction of single-nucleotide
substitutions at each site that would be synonymous, sum across the three
codon positions, and average across all clean codons in both ingroup and
outgroup at each analyzed position.

.. note::

   `Gossmann et al. (2012)`_ used the F3×4 codon-frequency model implemented in
   PAML, not Nei-Gojobori. The two methods give numerically different site
   counts, so MKado's ω, ω_a, and ω_na are **not strictly comparable** to
   numbers reported under the F3×4 convention. The decomposition itself
   (α × ω) is identical; the difference is in how :math:`L_n` and :math:`L_s`
   are estimated.

When ω_a / ω_na are reported (and when they are not)
-----------------------------------------------------

ω is well-defined and well-behaved per gene: it is the standard dN/dS ratio,
routinely reported per-gene by tools like PAML's codeml. MKado therefore
reports ω on every result type.

ω_a and ω_na, by contrast, **inherit the variance of α**. Per-gene Smith &
Eyre-Walker α is known to be noisy on the small Dn/Ds/Pn/Ps counts typical
of a single gene; multiplying by ω propagates that noise into ω_a. Following
the convention of `Gossmann et al. (2012)`_, `Galtier (2016)`_, and
`Coronado-Zamora et al. (2019)`_, MKado therefore restricts ω_a / ω_na to
result types whose α estimator is appropriate for the decomposition:

.. list-table::
   :header-rows: 1
   :widths: 35 15 15 35

   * - Result type
     - ω
     - ω_a / ω_na
     - α estimator feeding ω_a
   * - ``MKResult`` (per-gene MK)
     - ✓
     - ✗
     - per-gene Smith-Eyre-Walker α
   * - ``PolarizedMKResult`` (per-gene)
     - ✓
     - ✗
     - per-gene polarized α
   * - ``AsymptoticMKResult`` (per-gene)
     - ✓
     - ✓
     - asymptotic α (`Messer & Petrov 2013`_)
   * - ``AsymptoticMKResult`` (aggregated)
     - ✓
     - ✓
     - asymptotic α — **canonical Gossmann-style ω_a**
   * - ``ImputedMKResult``
     - ✓
     - ✓
     - imputed α (`Murga-Moreno et al. 2022`_)
   * - ``AlphaTGResult``
     - ✓
     - ✓
     - Tarone-Greenland weighted α

Rationale: per-gene Smith-Eyre-Walker α is noisy because Dn, Ds, Pn, Ps
counts on a single gene are typically small. Multiplying that noisy α by ω
produces a ω_a whose variance is dominated by α-variance, not by any
biological signal. The asymptotic, Tarone-Greenland, and imputed estimators
were designed for use with sparse counts and are explicitly intended for
this kind of decomposition.

If you need a per-gene ω_a estimate with proper uncertainty quantification,
consider the SnIPRE Bayesian hierarchical estimator (`Eilertson et al. 2012`_)
— not currently implemented in MKado.

Edge cases
----------

ω, ω_a, and ω_na return ``None`` (rendered as ``NA`` in TSV, ``null`` in JSON)
when:

- :math:`D_s = 0` — denominator undefined
- :math:`L_s = 0` or :math:`L_n = 0` — site counts undefined
- :math:`L_n` / :math:`L_s` were not computed (e.g. legacy
  ``mk_test_from_counts()`` calls without the ``ln`` / ``ls`` arguments)

Additionally, ω_a and ω_na are ``None`` when α is ``None`` — typically when
:math:`D_n = 0` or :math:`P_s = 0`.

References
----------

.. _Gossmann, Keightley & Eyre-Walker (2012): https://doi.org/10.1093/gbe/evs027
.. _Gossmann et al. (2012): https://doi.org/10.1093/gbe/evs027
.. _Galtier (2016): https://doi.org/10.1371/journal.pgen.1005774
.. _Coronado-Zamora et al. (2019): https://doi.org/10.1093/gbe/evz046
.. _Messer & Petrov 2013: https://doi.org/10.1073/pnas.1220835110
.. _Murga-Moreno et al. 2022: https://doi.org/10.1093/g3journal/jkac206
.. _Eilertson et al. 2012: https://doi.org/10.1371/journal.pgen.1002806

Gossmann TI, Keightley PD, Eyre-Walker A (2012). The effect of variation in
the effective population size on the rate of adaptive molecular evolution
in eukaryotes. *Genome Biology and Evolution* 4(5):658-667.
https://doi.org/10.1093/gbe/evs027

Coronado-Zamora M, Salvador-Martínez I, Castellano D, Barbadilla A, Salazar-Ciudad I (2019).
Adaptation and conservation throughout the *Drosophila melanogaster* life-cycle.
*Genome Biology and Evolution* 11(5):1463-1482.
https://doi.org/10.1093/gbe/evz046

Galtier N (2016). Adaptive protein evolution in animals and the effective
population size hypothesis. *PLoS Genetics* 12(1):e1005774.
https://doi.org/10.1371/journal.pgen.1005774

Nei M, Gojobori T (1986). Simple methods for estimating the numbers of
synonymous and nonsynonymous nucleotide substitutions. *Molecular Biology
and Evolution* 3(5):418-426.
