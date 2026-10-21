# inStrain pN/pS output for Supplementary Fig. 9

Tables and statistics for the pN/pS paragraph of Supplementary Note 6.2
("pN/pS analyses") — the SI caption needs *"the final pathway-level sample sizes,
filtering thresholds and BH-adjusted P values from the archived inStrain-derived
output"*. That is exactly what this package provides.

Built by `D_index_analysis/build_instrain_pnps.js`
(log: `D_index_analysis/out/instrain_pnps_log.txt`).

---

## 1. Inputs

| input | role |
|---|---|
| `SPOT_GV_instrain_gene_info/<sample>_gene_info.tsv` (167 files) | inStrain per-gene output, one file per SPOT sample; 786,133 gene records in total |
| `SPOT_GV_genes_pathway_for_Instrain.csv` (3,166 rows) | the curated gene list: `MAG, Order, Family, Protein_id, Gene, Pathway, Instrain_gene_id` |
| `D_index_analysis/out/fig1_tables/SPOT_sample_Dindex.csv` (167 rows) | sample date, size fraction, season, heatwave phase, campaign |

Joins and checks (all `PASS`, see `inStrain_pNpS_QC.csv`):

* every one of the 167 gene_info files was read; all have the standard 20-column header;
* `gene` in the gene_info file matches `Instrain_gene_id` in the curated list — 35,101
  gene × sample rows are carried over (35,078 distinct observations; two
  `Instrain_gene_id` values are each shared by two curated proteins, so those 23
  observations are reported once per protein/pathway);
* all 167 sample tokens resolve to a row of `SPOT_sample_Dindex.csv` (0 unmatched);
* **`pNpS_variants` is defined only when `SNV_S_count >= 1`.** Of the 6,459 mapped
  observations with `SNV_count >= 1`, exactly 5,551 carry a pN/pS and exactly 908 have
  `SNV_S_count = 0`. So 5,551 = 6,459 − 908 exactly. A gene with no *synonymous* variant
  has no pS denominator: those rows are **undefined, not zero**, and are excluded.
* `pNpS_variants = 0` is retained and is meaningful — it occurs in exactly the 2,005
  observations (154 in the conservative tier) where `SNV_N_count = 0`, i.e. synonymous
  variants but no nonsynonymous variant.

**One curated gene = one protein = one `Instrain_gene_id`.** A genome can carry several
copies of the same named gene (e.g. three HSP20 paralogues); those are separate proteins
in the curated list, so every count and every statistical unit is keyed on
(`Instrain_gene_id`, `Protein_id`) — never on (`genome`, `Gene`), which would silently
merge paralogues. `inStrain_pNpS_gene_summary.csv` records
`n_copies_of_this_gene_name_in_genome` so the collapse can be seen and reversed.

## 2. What `pNpS_variants` is, and what `dNdS_substitutions` is

inStrain's `gene_info.tsv` carries two ratio columns and four count columns:

| column | meaning |
|---|---|
| `SNV_count`, `SNV_S_count`, `SNV_N_count` | within-sample variant sites (all / synonymous / nonsynonymous) |
| `SNS_count`, `SNS_S_count`, `SNS_N_count` | the substitution-based counterpart counts used for dN/dS |
| `pNpS_variants` | **the pN/pS used here** |
| `dNdS_substitutions` | inStrain's substitution-based dN/dS |

The two ratios share one unreported, gene-specific pair of site counts
(*N*<sub>sites</sub> / *S*<sub>sites</sub>), which is why

```
pNpS_variants = (SNV_N / N_sites) / (SNV_S / S_sites)
dNdS_substitutions = (SNS_N / N_sites) / (SNS_S / S_sites)

pNpS_variants == dNdS_substitutions * (SNV_N * SNS_S) / (SNV_S * SNS_N)
```

This identity was verified to machine precision on **every** row in which both ratios
and all four counts are present (893 of 893 rows exact, worst relative error 0;
`D_index_analysis/tools/probe_pnps_formula2.js`, log
`out/instrain_pnps_formula2_log.txt`). It confirms that `pNpS_variants` is built from the
`SNV_*` counts and that `pN/P_S` therefore means *within-sample* polymorphism, not
between-sample divergence.

> The full definitions of the inStrain columns are not part of the supplied package; cite
> the inStrain paper (Olm et al. 2021) for them. In particular the exact semantics of
> `breadth` and `breadth_minCov` should be taken from inStrain's documentation — both are
> carried through verbatim in the master table and are not interpreted here.

## 3. Thresholds actually applied

Three separate layers — full machine-readable version in
`inStrain_pNpS_thresholds_and_provenance.csv`.

### Layer 1 — inStrain run settings (the thresholds stated in the manuscript)

> "Reads were mapped to dereplicated GVMAGs with BWA-MEM and processed with inStrain using
> a minimum site coverage of **10×** and **minor-allele-frequency threshold of 0.05**."
> — manuscript Methods / Supplementary Note 6.1

Both are baked into the supplied values, so they are not re-applied. The MAF ≥ 0.05 rule
means every counted `SNV_*` site has a within-sample allele frequency of at least 5%; no
per-variant allele frequencies are exported in `gene_info.tsv`, so this threshold cannot
be re-derived from the file and is carried over from the run settings.

### Layer 2 — pN/pS definedness

`SNV_S_count >= 1` (see §1). This is the only filter that changes the *set of genes* for
which pN/pS exists, and it is imposed by inStrain itself.

### Layer 3 — reporting filters (five tiers)

| tier | filter | observations | genes | genomes | samples | pathways with data | median pN/pS (gene level) |
|---|---|---|---|---|---|---|---|
| **A** | coverage ≥ 10 (conservative subset) | 1,087 | 497 | 130 | 64 | 13 | 0.1162 |
| B | coverage ≥ 5 | 1,640 | 741 | 166 | 91 | 14 | 0.1100 |
| C | coverage ≥ 1 | 3,564 | 1,131 | 205 | 152 | 14 | 0.0906 |
| **D** | **no coverage filter — full analysis set** | **5,551** | **1,219** | **213** | **166** | **15** | **0.0829** |
| E | coverage ≥ 10 and breadth ≥ 0.5 | 1,001 | 487 | 125 | 17 | 13 | 0.1175 |

**Which tier to report.** The stated 10× rule is a *per-site* setting inside inStrain, so
it is already respected by every value in these files. Tier D — every observation in which
inStrain reported a pN/pS — is therefore the analysis set implied by the stated methods,
and is the recommended one. Requiring the *gene mean* coverage to be ≥ 10× (tier A) is a
much stricter, additional filter that keeps only 64 of 167 samples; it is reported as the
conservative subset. **Tier A and tier D agree qualitatively on every contrast tested and
differ only in power** (§5). Exactly one sample (JUN-2013) contributes no defined pN/pS at
any coverage.

## 4. Statistics

* **Unit of analysis = gene.** The single best-sampled gene occurs in 97 samples, so gene ×
  sample observations are heavily pseudoreplicated; testing them directly would inflate
  *n* roughly a hundredfold. Each gene contributes one value: the median of its pN/pS
  across the samples in which it was observed.
* **Pathway vs background**: Mann–Whitney *U* (tie-corrected, continuity-corrected normal
  approximation), two-sided, with Cliff's delta, a 2,000-replicate bootstrap CI on delta,
  and the Hodges–Lehmann shift. Background = all other curated genes with a defined pN/pS
  in the same tier.
* **Order contrast**: Imitervirales vs Algavirales, at **two units** — gene (recommended)
  and genome (per-genome median, as an independence check).
* **Secondary**: Fisher exact on the fraction of genes with pN/pS > 1.
* **Multiplicity**: Benjamini–Hochberg, computed within tier across the pathways actually
  tested, and within tier across scopes for the order tests.
* An observation-level *P* is also reported in `inStrain_pNpS_pathway_tests.csv`
  (`P_observation_level_pseudoreplicated`) and is **explicitly not** the reported statistic.

## 5. Results to quote

Full number pack: `inStrain_pNpS_SI_caption_numbers.csv`.

**Tier D — full analysis set (recommended).** 5,551 gene × sample observations, 1,219
curated genes (proteins), 213 GVMAGs, 166 samples, 15 pathways. Median pN/pS = 0.083
(gene level), IQR 0.00–0.18; 12 of 1,219 genes (1.0%) have pN/pS > 1. 13 pathways have
enough genes to test; **2 survive BH correction**:

| pathway | median pathway vs background | *P* | BH-*P* | Cliff's δ [95% CI] |
|---|---|---|---|---|
| TCA cycle (49 vs 1,170 genes) | 0.025 vs 0.088 | 1.0 × 10⁻⁶ | **1.3 × 10⁻⁵** | −0.408 [−0.514, −0.292] |
| Ubiquitin system (411 vs 808 genes) | 0.115 vs 0.077 | 2.0 × 10⁻⁴ | **1.3 × 10⁻³** | +0.129 [+0.056, +0.199] |

**Tier C (coverage ≥ 1×)** gives the same two contrasts: TCA cycle
(0.043 vs 0.095; *P* = 6.1 × 10⁻⁴, BH-*P* = **7.3 × 10⁻³**, δ = −0.303) and ubiquitin
system (0.116 vs 0.084; *P* = 2.9 × 10⁻³, BH-*P* = **1.7 × 10⁻²**, δ = +0.107).

**Tier A (conservative, coverage ≥ 10×).** 1,087 observations, 497 genes, 130 GVMAGs, 64
samples, 13 pathways. Median pN/pS = 0.116 (gene level), IQR 0.058–0.20; 7 of 497 genes
(1.4%) exceed 1. 11 pathways are testable and **none survives BH correction** (smallest raw
*P*: TCA cycle 0.029, glycolysis 0.043; both BH-*P* = 0.23). This is a power difference,
not a disagreement: the direction of effect is the same as in tier D (TCA lower,
ubiquitin higher); `inStrain_pNpS_pathway_level.csv` shows all five tiers side by side.

**Order contrast.** pN/pS does not differ between Imitervirales and Algavirales at any
tier or unit (tier A, gene unit: median 0.119 vs 0.105, *P* = 0.68, δ = +0.024
[−0.090, +0.143]; tier A, genome unit: 0.110 vs 0.110, *P* = 0.88; tier D, gene unit:
0.083 vs 0.089, *P* = 0.50). Within pathways, the only nominal signal is heat-shock
proteins (tier A: Imitervirales 0.110 vs Algavirales 0.066, *P* = 0.015, BH-*P* = 0.11),
the same direction at every tier, but it never survives correction.

**Interpretation wording (per SI Note 6.2).** pN/pS > 1 cannot by itself distinguish
positive selection from relaxed purifying constraint. Use *"elevated pN/pS compatible with
relaxed purifying selection or positive selection"*, never treat pN/pS > 1 as proof of
adaptation. This package's numbers also do not support any claim of *higher* pN/pS in
Imitervirales overall.

### Suggested caption wording

> **pN/pS across selected pathways.** Synonymous and nonsynonymous within-sample variant
> counts were obtained with inStrain (minimum site coverage 10×; minor-allele frequency
> ≥ 0.05) for the curated GVMAG gene set; pN/pS is `pNpS_variants`. pN/pS is undefined for
> a gene with no synonymous variant, so the analysis set is the 5,551 gene × sample
> observations (1,219 genes from 213 GVMAGs in 166 samples) in which inStrain reported a
> value. To avoid pseudoreplication, each gene contributes a single value — the median of
> its pN/pS across samples — and pathways were compared with the remaining curated genes
> by two-sided Mann–Whitney *U* tests with Benjamini–Hochberg correction across the 13
> testable pathways. TCA-cycle genes had lower pN/pS than the genomic background
> (median 0.025 versus 0.088; Cliff's δ = −0.41, 95% CI −0.51 to −0.29; *P* = 1.0 × 10⁻⁶,
> BH-adjusted *P* = 1.3 × 10⁻⁵) and ubiquitin-system genes higher (0.115 versus 0.077;
> δ = +0.13, 95% CI +0.06 to +0.20; *P* = 2.0 × 10⁻⁴, BH-adjusted *P* = 1.3 × 10⁻³); no
> other pathway was significant. Overall pN/pS did not differ between Imitervirales and
> Algavirales (gene-level medians 0.083 versus 0.089; *P* = 0.50; genome-level
> 0.069 versus 0.079; *P* = 0.16). Values above 1 — 1.0% of genes — are compatible with
> relaxed purifying selection or positive selection rather than constituting evidence for
> adaptation. Restricting to genes with mean coverage ≥ 10× retains the same direction of
> effect for both pathways but leaves 11 testable pathways and no BH-significant contrast
> (Table X).

Replace the *n*, pathway names and *P* values if a different tier is chosen — every variant
is in `inStrain_pNpS_SI_caption_numbers.csv` and
`inStrain_pNpS_threshold_sensitivity.csv`.

## 6. Files

All in `tables/instrain/` (mirrored in `D_index_analysis/out/instrain/`).

| file | rows | contents |
|---|---|---|
| `inStrain_pNpS_gene_by_sample_master.csv` | 35,101 | **the mapping**: every curated-gene × sample observation, with sample date/size fraction/season/phase, MAG/order/phylum/family, pathway/gene/protein id, all raw inStrain columns (coverage, breadth, breadth_minCov, nucl_diversity, SNV/SNS counts, both ratios) and the five tier flags |
| `inStrain_pNpS_gene_by_sample_tierA_coverage10.csv` | 1,087 | the conservative subset, same 48 columns |
| `inStrain_pNpS_gene_summary.csv` | 2,903 | per curated gene (protein): the copy number of that gene name in its genome, tier A and all-tier observation counts, median/mean/spread of pN/pS, % observations > 1, median dN/dS and coverage |
| `inStrain_pNpS_genome_level.csv` | 256 | per GVMAG: number of curated genes with data, observations, samples, median pN/pS (tier A and all), % > 1, median π |
| `inStrain_pNpS_pathway_level.csv` | 75 | pathway × tier (15 × 5): gene/observation/sample/genome counts, median and IQR of pN/pS at gene and observation level, unadjusted SNV_N:SNV_S count ratio, % of genes > 1, order composition |
| `inStrain_pNpS_pathway_tests.csv` | 75 | pathway vs background per tier: *n*, medians, *U*, *z*, *P*, **BH-*P***, Cliff's δ + bootstrap CI, Hodges–Lehmann, Fisher *P* on pN/pS > 1, and the flagged pseudoreplicated observation-level *P* |
| `inStrain_pNpS_order_tests.csv` | 85 | Imitervirales vs Algavirales per tier, at gene and genome unit, overall and within each pathway, with BH-*P* |
| `inStrain_pNpS_threshold_sensitivity.csv` | 23 | the five tiers plus an 18-cell coverage × breadth ladder, each with counts and medians |
| `inStrain_pNpS_thresholds_and_provenance.csv` | 12 | every threshold, where it comes from, and the evidence that it was applied |
| `inStrain_pNpS_SI_caption_numbers.csv` | 57 | ready-to-quote numbers for the Supplementary Fig. 9 caption |
| `inStrain_pNpS_QC.csv` | 25 | the checks listed in §1, including the three exact identities |

## 7. Caveats a reviewer may raise

1. **Only the curated gene list is covered** — 3,166 proteins in 15 pathways (~50 gene
   families), of which 2,901 were observed in at least one sample. This is a curated
   subset, not whole-genome content, so "pathway pN/pS" means the pN/pS of the curated
   genes assigned to that pathway. Pathway *n* therefore ranges from 1 (glyoxylate shunt)
   to 411 (ubiquitin system) genes; pathways with fewer than 5 genes with data are listed
   as `not tested` rather than tested.
2. **The 5,551 observations are 1,219 genes observed across samples, not independent
   replicates.** All tests use the gene (or genome) as the unit; the observation-level test
   is included only to show the size of the pseudoreplication effect.
3. **Genes within one genome are not independent either.** The genome-unit order test is
   the independence check; it agrees with the gene-unit test in every tier.
4. **The background is "all other curated genes"**, so a pathway contrast partly reflects
   the genome composition of the pathway's genes (Imitervirales-dominated vs
   Algavirales-dominated pathways). This is why the order tests are reported at both the
   overall and the within-pathway level.
5. **Coverage is uneven across the curated genes** — median gene mean coverage among
   mapped observations is 0.17×, which is why the tier ladder is provided rather than a
   single threshold.

## 8. Reproducing

```powershell
# main build (writes tables/instrain/ and D_index_analysis/out/instrain/)
node .\D_index_analysis\build_instrain_pnps.js

# derivation checks for the two ratio columns
node .\D_index_analysis\tools\probe_pnps_formula2.js

# column / join profiling
node .\D_index_analysis\tools\profile_instrain.js
node .\D_index_analysis\tools\profile_instrain_tiers.js

# inspect the finished tables (caption | tiers | tests | order | pathways)
node .\D_index_analysis\tools\show_instrain_tables.js caption
```

`tools/docx2txt.js` was added to read the SI docx (Supplementary Note 6.2 supplies the
requirement quoted at the top of this file); the extracted text is
`D_index_analysis/out/SI.txt`.
