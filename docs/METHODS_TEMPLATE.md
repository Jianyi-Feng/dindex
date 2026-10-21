# Adaptable Methods paragraph

Replace the bracketed items.

> **D-index.** For each sample we calculated D = (P − M)/(P + M). P and M are the summed [abundance measure, e.g. read counts
> or coverage] of [group P, e.g. *Prochlorococcus* ASVs] and [group M] assigned with [classifier/database and release] at the
> [rank] level. D ranges from −1 (only M detected) to +1 (only P detected). For positive P and M it equals tanh[½ ln(P/M)], a
> bounded transform of the two-part log-ratio, and is therefore invariant to sequencing depth and to normalization applied
> equally to both groups. Samples with P + M = 0 were excluded ([n] samples). Samples in which one group was absent took
> D = ±1 ([n] samples). Samples in which P and M together contributed less than [threshold] of [total community measure] were
> flagged. Per-sample 95% intervals were obtained from a Beta(P + ½, M + ½) distribution for p = P/(P + M) [or: by
> bootstrapping features within each group, 2,000 replicates, as a descriptive feature-composition sensitivity analysis]. The sensitivity of D to [driver] was estimated by linear
> regression [with covariates] and [cluster-robust/Newey–West] errors, and the value of [driver] at which D = 0 was estimated
> with a [station-level] bootstrap. Calculations used dindex v[version] (ref.).
