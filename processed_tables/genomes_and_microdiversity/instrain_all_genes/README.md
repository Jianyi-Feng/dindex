# inStrain pN/pS for **all genes** of the 266 GVMAGs

Companion to `D_index_analysis/INSTRAIN_PNPS_README.md` (which covers only the ~3,166
curated pathway proteins). Same pN/pS definition, same three threshold layers, same
statistics — but the analysis universe is now **every gene of the 266 MAGs** listed in
`SPOT_NCLDV_tax.csv`, annotated with the KEGG KO and Pfam tables, so the results no longer
depend on the hand-curated pathway list.

Built by `D_index_analysis/build_instrain_allgenes.js`
(log: `D_index_analysis/out/allgenes_build_log.txt`).
All tables in `tables/instrain_all_genes/` (mirrored in `D_index_analysis/out/instrain_all_genes/`).

---

## 1. Scope

| | curated build | **this build** |
|---|---|---|
| genes in scope | 3,166 curated proteins (~50 gene families, 15 pathways) | **84,846 genes — every gene of the 266 MAGs** |
| gene × sample observations with a defined pN/pS | 5,551 | **124,071** |
| genes with ≥1 defined pN/pS | 1,219 | **33,483** |
| genomes | 213 | **264 (of 266)** |
| annotation | 15 hand-curated pathways | **KEGG KO + BRITE L1/L2 + KEGG pathways + Pfam families** (431 pathways and 1,524 Pfam families exist; 392 pathways / 1,613 KOs / 956 Pfam families have data in tier D) |

## 2. Inputs and how genes were tied to taxonomy

| input | role |
|---|---|
| `SPOT_GV_instrain_gene_info/<sample>_gene_info.tsv` (167 files) | 786,133 gene records |
| `SPOT_NCLDV_tax.csv` | **266 MAGs**, `genome, Order, Family` |
| `SPOT_GV_nucl_diversity_metadata.csv` | second taxonomy source; also the π values |
| `SPOT_GV_bins_KO_annotation_final.csv` | KEGG KO + BRITE hierarchy, 80,146 rows / 30,056 genes |
| `SPOT_GV_bin_pfam_annotation_final.csv` | Pfam, 34,140 rows / 34,140 genes / 1,524 families |
| `SPOT_GV_genes_pathway_for_Instrain.csv` | flags which genes the curated build already covers |
| `D_index_analysis/out/fig1_tables/SPOT_sample_Dindex.csv` | sample date, size fraction, season, phase |

**Gene → MAG.** Every gene id is `<MAG>_contig_<n>_<k>` and `scaffold` is
`<MAG>_contig_<n>`, so the MAG is the scaffold with `_contig_<n>` stripped. Result:
**786,133 / 786,133 records resolve, 0 by the gene-id fallback, 0 outside the taxonomy.**

**Taxonomy naming (important).** `SPOT_NCLDV_tax.csv` has 266 rows, but the 25
Mirusviricota MAGs are written there in an intermediate form —
`Mirus_201408_bin_4`, `Mirus_2015_07VIROME_bin_205` — while the gene ids and both
annotation files use the canonical form `Mirus_AE_201408_bin_4`,
`Mirus_virome_2015_07_bin_205`. Canonicalising those 25 names makes **SPOT_NCLDV_tax.csv
cover all 266 MAGs**, and the two taxonomy sources then agree completely (0 MAGs in one
source only). Without that step, 2,904 genes (all Mirusviricota) would have dropped out.
`Order` is identical in both sources; the only disagreements are family labels of the
`incertae_sedis` families (`incertae_sedis` vs `AG_incertae_sedis` / `IM_incertae_sedis`),
and both spellings are carried in the genome table.

## 3. pN/pS and the thresholds applied

Identical to the curated build; see that README for the full derivation.

* **inStrain run settings** (manuscript Methods / SI Note 6.1): minimum site coverage
  **10×**, minor-allele frequency **≥ 0.05**. Both are baked into the supplied values.
* **Definedness**: `pNpS_variants` exists only when `SNV_S_count ≥ 1`. Verified exactly on
  the full gene set: **124,071 observations have a pN/pS and exactly the same 124,071 have
  `SNV_S_count ≥ 1`**. The remaining 662,062 records split into 333,891 with an explicit
  `SNV_S_count = 0` (variants called, none synonymous → no pS denominator) and 328,171 with
  a completely blank SNV/SNS block (no variant called at all).
* **Derivation re-verified on all genes**, not just the curated list: the identity
  `pNpS_variants == dNdS_substitutions × (SNV_N × SNS_S) / (SNV_S × SNS_N)` holds exactly on
  **21,475 / 21,475** testable observations (worst relative error 0).
* **Reporting tiers**

  | tier | filter | observations | genes | genomes | median pN/pS (gene level) |
  |---|---|---|---|---|---|
  | A | coverage ≥ 10× (conservative) | 52,543 | 18,589 | 232 | 0.1341 |
  | B | coverage ≥ 5× | 71,180 | 24,586 | 251 | 0.1262 |
  | C | coverage ≥ 1× | 105,410 | 31,958 | 263 | 0.1151 |
  | **D** | **no coverage filter (full analysis set)** | **124,071** | **33,483** | **264** | **0.1124** |
  | E | coverage ≥ 10× and breadth ≥ 0.5 | 50,102 | 18,371 | 216 | 0.1340 |

  Plus an 18-cell coverage × breadth ladder. All 167 samples contribute.

**One curated gene = one protein.** Every count and every statistical unit is keyed on
(`Instrain_gene_id`, `Protein_id`), never on (`genome`, `Gene`), so paralogous copies of the
same gene name are separate units.

## 4. Statistics

* **Unit of analysis = gene** (per-gene median pN/pS across the samples in which it was
  observed). Observations are heavily pseudoreplicated — the best-sampled single gene occurs
  in 142 of the 167 samples — so testing rows directly would inflate *n* by ~3.7×
  (124,071 observations / 33,483 genes).
* **Mann–Whitney *U***, two-sided, tie- and continuity-corrected, Cliff's δ with an analytic
  (DeLong) 95% CI, Hodges–Lehmann shift (computed for the BH-significant groups), Fisher
  exact on the fraction of genes with pN/pS > 1, **Benjamini–Hochberg** within tier and
  background.
* **Two backgrounds, deliberately.** Functional groups are tested twice:
  1. `all other genes with a defined pN/pS` — the headline comparison; and
  2. `all other <KEGG- or Pfam->annotated genes` — a comparison *within* annotation status.
  This matters because annotation status is itself strongly associated with pN/pS (§5.3), so
  the second background separates "this pathway behaves differently" from "this gene is
  annotated and annotated genes differ".
* The rank index for every test is built once per tier (and per background), so the ~5,600
  group tests are exact and fast rather than bootstrapped.

## 5. Results

### 5.1 Universe
84,846 genes observed across the 266 MAGs; 33,483 (39.5%) have at least one defined pN/pS;
124,071 gene × sample observations. Only 2 MAGs have no defined pN/pS anywhere. Among the
33,483 genes with data: 13,002 have a KO, 6,904 sit on ≥1 KEGG pathway, 15,417 have a Pfam
domain, and 16,531 have neither. 1,218 are in the curated pathway list.

### 5.2 Order contrasts
Imitervirales vs Algavirales, gene and genome unit, all five tiers
(*P* values below 10⁻³⁰⁰ are shown as < 1 × 10⁻³⁰⁰):

| tier | unit | Imitervirales | Algavirales | *P* | BH-*P* | Cliff's δ [95% CI] |
|---|---|---|---|---|---|---|
| A (cov ≥ 10×) | gene | 0.153 (n = 10,457) | 0.115 (n = 8,025) | < 1 × 10⁻³⁰⁰ | < 1 × 10⁻³⁰⁰ | **+0.149** [+0.130, +0.168] |
| A | genome | 0.168 (n = 128) | 0.119 (n = 81) | 3.2 × 10⁻⁷ | 8.8 × 10⁻⁷ | **+0.420** [+0.262, +0.559] |
| B (cov ≥ 5×) | gene | 0.135 | 0.110 | < 1 × 10⁻³⁰⁰ | < 1 × 10⁻³⁰⁰ | +0.106 |
| C (cov ≥ 1×) | gene | 0.119 (n = 22,204) | 0.102 (n = 9,419) | < 1 × 10⁻³⁰⁰ | < 1 × 10⁻³⁰⁰ | +0.059 |
| C | genome | 0.132 (n = 139) | 0.106 (n = 81) | 1.7 × 10⁻⁴ | 5.2 × 10⁻⁴ | +0.304 |
| D (no filter) | gene | 0.116 (n = 23,573) | 0.100 (n = 9,488) | 2.0 × 10⁻⁸ | 1.9 × 10⁻⁷ | **+0.039** [+0.025, +0.053] |
| D | genome | 0.120 (n = 139) | 0.104 (n = 81) | 0.020 | 0.040 | +0.188 [+0.054, +0.343] |

So **across the whole gene set Imitervirales genes have consistently higher pN/pS than
Algavirales genes** — significant at every tier and at both units — but the effect size
depends strongly on the coverage tier: δ falls from ≈ 0.15 (gene) / 0.42 (genome) at
≥10× coverage to 0.04 / 0.19 in the unfiltered set. The noisy low-coverage genes dilute it.
**This differs from the curated-gene result**, where no difference was detectable
(tier D *P* = 0.50): the curated list is dominated by a few large pathways (ubiquitin system
alone is 411 of 1,219 genes) and is not representative of the genome.

**Pandoravirales** genes are markedly higher than all other genes (tier A: 0.397 vs 0.134,
δ = +0.446; genome level 1.26 vs 0.140, δ = +0.540; the same direction in every tier), and
**MR_01** (Mirusviricota) lower (tier D: 0.053 vs 0.113, δ = −0.163). Nucleocytoviricota vs
Mirusviricota overall is not significant at any tier.

### 5.3 The annotation-status confound (read this before quoting §5.4)
Tier D, each group against all other genes:

| group | genes | median pN/pS | background | Cliff's δ [95% CI] |
|---|---|---|---|---|
| KO + Pfam | 11,467 | 0.0865 | 0.127 | **−0.161** [−0.170, −0.151] |
| Pfam only | 3,950 | 0.101 | 0.114 | −0.041 [−0.058, −0.024] |
| KO only | 1,535 | 0.137 | 0.111 | +0.113 [+0.085, +0.142] |
| unannotated (no KO, no Pfam) | 16,531 | 0.133 | 0.0935 | **+0.142** [+0.133, +0.151] |

Annotated genes carry a substantially lower pN/pS than unannotated ones, in every tier
(tier A: KO + Pfam δ = −0.136, unannotated δ = +0.127). Any functional comparison against
the whole gene set therefore partly measures *annotation*, not function — which is exactly
why every family is also tested against an annotation-matched background (§5.4).

### 5.4 Functional groups
Number of groups surviving BH < 0.05, **tier D**:

| background | groups tested | BH < 0.05 |
|---|---|---|
| all other genes with a defined pN/pS | 1,242 | 277 |
| other KEGG-annotated genes (pathway / KO / BRITE) | 810 | 186 |
| other Pfam-annotated genes | 428 | 73 |

The surviving signal against the **annotation-matched** background is coherent and matches
what selection theory predicts. Tier A, gene unit:

* **lower pN/pS** (stronger purifying selection) — informational and housekeeping functions:
  nucleotide excision repair δ = −0.236, pyrimidine metabolism −0.219, mismatch repair
  −0.207, glutathione metabolism −0.205, DNA replication −0.173, base excision repair
  −0.155, purine metabolism −0.141, RNA polymerase −0.131; Pfam PCNA_C −0.390, ubiquitin
  −0.349, Ribonuc_red_sm −0.312, HSP70 −0.294, DNA_topoisoIV −0.255.
* **higher pN/pS** (relaxed constraint / possible positive selection) — accessory and
  lineage-specific families: DUF285 δ = +0.750, Peptidase_S74 +0.543, Pkinase +0.344,
  HNH_3 +0.643, DUF5773 +0.366; tier D also Ank +0.784 (n = 16), F5_F8_type_C_2 +0.789
  (n = 13), DUF5855 +0.526, Collagen +0.318.
* Note that AAA is *not* significant at tier A (δ = −0.038, *P* = 0.53) although it is at
  tier D (δ = −0.172) — a reminder that δ, not just *P*, must be read.
* The complete table — every group, both backgrounds, three tiers, with *n*, medians, *P*,
  BH-*P*, δ and its DeLong CI — is `inStrain_allgenes_pNpS_functional_groups.csv`.

Because with tens of thousands of genes even δ ≈ 0.1 is highly significant, **quote δ and
its CI, not just *P***. A group with BH < 0.05 and |δ| < 0.1 is statistically detectable but
biologically small, and the difference between the two backgrounds shows how much of any
raw contrast is annotation status.

### 5.5 Heatwave / size-fraction summary
`inStrain_allgenes_pNpS_sample_summary.csv` (1,035 rows: sample × size fraction × order) and
`inStrain_allgenes_pNpS_phase_fraction_summary.csv` (42 rows: phase × fraction × order) give
per-sample medians so the time series can be inspected directly. **These are descriptive
only** — the unit is the sample, and genes recur across samples, so no *P* values are
attached. For orientation: within the viral (0.02–0.2 µm) fraction the per-sample median
pN/pS is similar in pre and heatwave samples for both Imitervirales (0.229 → 0.252) and
Algavirales (0.140 → 0.126), while the micro/nanoplankton fraction sits much higher in both
orders; a formal test needs a paired, gene-matched design and has not been run here.

## 6. Files

All in `tables/instrain_all_genes/` (mirrored in `D_index_analysis/out/instrain_all_genes/`).

| file | rows | contents |
|---|---|---|
| `inStrain_allgenes_pNpS_gene_level.csv` | 84,846 | **the main table.** Every gene: genome, phylum, order, family, whether it is in the curated list, KO, KO symbol, KO description, KEGG L1/L2, all KEGG pathways, number of KOs and pathways, Pfam name/accession/description, annotation status, observation and sample counts, pN/pS stats at tier A and for the full set (median, mean, IQR, min, max, % > 1), median dN/dS, median coverage, median π |
| `inStrain_allgenes_pNpS_gene_by_sample_defined.csv` | 124,071 | every gene × sample observation with a defined pN/pS, with sample date/size fraction/season/phase, MAG/order/family, coverage, breadth, π, all SNV/SNS counts, both ratios, and tier flags — the raw backing for re-analysis |
| `inStrain_allgenes_pNpS_genome_level.csv` | 266 | per GVMAG: gene/observation/sample counts, median pN/pS (tier A and full set), median of gene medians, % > 1, median coverage, median π, taxonomy |
| `inStrain_allgenes_MAG_taxonomy.csv` | 266 | the taxonomy used: order, family, phylum, the family/order name in the second source when it differs, the name as written in `SPOT_NCLDV_tax.csv`, taxonomy source, and per-MAG gene/observation counts |
| `inStrain_allgenes_pNpS_order_stats.csv` | 110 | every order vs the rest, plus Imitervirales vs Algavirales and Nucleocytoviricota vs Mirusviricota, at gene and genome unit, all five tiers, with *P*, BH-*P*, δ + CI |
| `inStrain_allgenes_pNpS_functional_groups.csv` | 6,574 | KEGG pathway / KO / BRITE L1 / BRITE L2 / Pfam family / annotation status, in three tiers, against two backgrounds, with *n*, medians, *P*, BH-*P*, δ + DeLong CI, HL shift, Fisher *P* |
| `inStrain_allgenes_pNpS_threshold_sensitivity.csv` | 18 | the five named tiers + 13 ladder cells: observations, genes, genomes, how many KEGG pathways / KOs / Pfam families have data, medians and % > 1 |
| `inStrain_allgenes_pNpS_sample_summary.csv` | 1,035 | per sample × size fraction × order: observation/gene counts, median and IQR of pN/pS, % > 1 |
| `inStrain_allgenes_pNpS_phase_fraction_summary.csv` | 42 | phase × size fraction × order: median of per-sample medians (descriptive) |
| `inStrain_allgenes_pNpS_QC.csv` | 25 | every join and identity check from §3 |

## 7. Caveats

1. **Only 39.5% of genes have a defined pN/pS** and only 15.8% of gene × sample records do.
   A gene needs a synonymous variant at ≥10× site coverage to have a pS denominator, so the
   analysis set is biased towards genes with enough coverage. The tier ladder quantifies this:
   tier A keeps 52,543 of the 124,071 observations but only 18,589 of the 33,483 genes.
2. **Tier A and tier D disagree in magnitude, not in direction, for the order contrast**
   (δ = +0.149 vs +0.039 at the gene unit). The high-coverage subset gives the cleaner
   estimate; the unfiltered set is diluted by genes whose pN/pS rests on very shallow data.
3. **The order contrast is not gene-set matched.** Imitervirales genomes are larger and carry
   a different mix of Pfam families, so part of the order difference could be composition.
   A within-Pfam-family matched test (compare the two orders only among genes of the same
   family, then combine) is the natural next step and has **not** been run.
4. **The functional-group table tests ~1,200 groups per tier × 2 backgrounds.** BH controls
   the false discovery rate within each family and background, not across the whole table.
5. **Annotation coverage is low** (28,032 of 84,846 genes have a KO, 31,996 a Pfam hit), so
   the "unannotated" class is large and heterogeneous; §5.3 shows it is not neutral with
   respect to pN/pS, which is why the within-annotation background is the one to prefer for
   functional claims.
6. `gene_length` is blank in a substantial share of the inStrain records; it is reported as
   given and was not used as a filter.
7. **Tier D's larger δ values against the all-genes background are partly the confound of
   §5.3** — the background there includes 16,531 unannotated genes with higher pN/pS. Always
   state which background a δ comes from; the `background_definition` column does that.
8. **Most per-gene values rest on very few observations.** Per-gene counts:
   tier A — median **1** observation (p75 = 3, max = 109), and **68.2% of the 18,589 tier-A
   genes have ≤ 2 observations**; full set — median 2 (p75 = 4, max = 142), 61.8% with ≤ 2.
   A gene's "median pN/pS" is therefore often a single number. `n_obs_tierA` and `n_obs_all`
   are in the table so a minimum observation count can be imposed; the group tests do **not**
   impose one beyond the tier filter and the ≥5-gene minimum, which makes individual tier-A
   effect sizes noisier than the aggregate statistics suggest.

## 8. Reproducing

```powershell
# main build (writes tables/instrain_all_genes/ and D_index_analysis/out/instrain_all_genes/)
node --max-old-space-size=8192 .\D_index_analysis\build_instrain_allgenes.js

# scope profiling (taxonomy resolution, gene/annotation coverage)
node .\D_index_analysis\tools\profile_allgenes.js

# inspect the finished tables (groups | order | tax | phase | tiers)
node .\D_index_analysis\tools\show_allgenes_tables.js groups
```
