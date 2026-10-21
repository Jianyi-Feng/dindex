# Minimum reporting checklist for a D-index

Use this list when publishing a D-index for any taxon pair, so that values can be compared and reproduced across studies.

## Definition

- [ ] **Groups.** P and M named exactly, with the taxonomy (database, release) and rank used to assign features.
- [ ] **Features.** What was summed: marker gene (e.g. PolB), ASVs, OTUs, MAGs, genes or reads.
- [ ] **Abundance units and normalization:** counts, coverage, TPM or relative abundance. D is invariant to normalization
  applied equally to P and M, but not to gene-length or copy-number corrections that differ between the groups.
- [ ] **Pair selection.** How the pair was chosen, on which data, and whether it was fixed before the tests reported.

## Quality control

- [ ] **Zero handling.** Number of samples with D = +1 or −1 (one group absent), and with P + M = 0 (undefined, excluded).
- [ ] **Share.** Distribution of (P + M)/total, and any `min_share` or `min_total` thresholds applied.
- [ ] **Protocol consistency.** Size fraction, extraction and sequencing protocol for every comparison. Do not pool across
  protocols without a design covariate.
- [ ] **Log-ratio sensitivity.** Any pseudocount used for log-ratio analyses, with results under alternative pseudocounts.

## Uncertainty

- [ ] **Per-sample intervals.** Method used: Beta interval for read counts, or descriptive feature-composition bootstrap for
  coverage/TPM. Do not present the latter as a sampling-error confidence interval.
- [ ] **Dependence structure.** Stations, repeated sampling or time; the robust errors or resampling units used.

## Calibration and interpretation

- [ ] **Sensitivity.** Slope of D on each environmental driver, with CI, covariates and error type.
- [ ] **Zero point.** Driver value at D = 0 with a bootstrap CI, reported separately for each system and fraction.
- [ ] **Baselines.** Period, harmonics and baseline R² for anomaly analyses; event inference at the year or event level.
- [ ] **Transfer.** Whether sensitivity and zero point were tested on independent data (other sites, basins or years).
- [ ] **Caveats.** D is relative, taxon-specific biases change it, and associations are not causal attribution.

## Reproducibility

- [ ] **Software version** of `dindex`, or the exact formula and code.
- [ ] **Sample-level output.** A table with `plus`, `minus`, `total`, `D`, `share`, `qc` and any intervals.
