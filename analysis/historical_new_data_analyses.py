"""Audit the 21 September data deposit and build source-auditable  outputs.

This script restores only analyses supported by the deposited row-level or finalized
machine-readable results. It deliberately keeps CCA/Mantel and rhodopsin phylogeny/
structure claims excluded because their required matrices, sequences, alignments or
model files are not present.
"""
from pathlib import Path
import math

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REV = HERE.parent
DATA = REV.parents[1] / "raw_tables"
SUPP = REV / "supplementary"
FIG = SUPP / "figures"
FIG.mkdir(exist_ok=True)

BLUE = "#32658f"
ORANGE = "#d98221"
TEAL = "#69b3b7"
INK = "#252525"
GREY = "#8a8a8a"
LIGHT = "#dedede"


def require(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def read_csv(relative):
    return pd.read_csv(DATA / relative)


# --------------------------------------------------------------------- load inputs
polb = read_csv("SPOT_PolB_order_relative_abundance_long.csv")
spot = read_csv("SPOT_sample_Dindex.csv")
nodes = read_csv("SPOT_18S_GV_network_nodes.csv")
edges = read_csv("SPOT_18S_GV_network_edges.csv")
tree_tax = read_csv("SPOT_ref_NCLDV_taxonomy_for_iqtree.csv")
rhod = pd.read_excel(DATA / "SPOT_GV_and_ref_rhodopsin_metadata.xlsx", sheet_name="Final")

ag_order = read_csv("instrain_all_genes/inStrain_allgenes_pNpS_order_stats.csv")
ag_fun = read_csv("instrain_all_genes/inStrain_allgenes_pNpS_functional_groups.csv")
ag_thr = read_csv("instrain_all_genes/inStrain_allgenes_pNpS_threshold_sensitivity.csv")
ag_qc = read_csv("instrain_all_genes/inStrain_allgenes_pNpS_QC.csv")
ag_gene = read_csv("instrain_all_genes/inStrain_allgenes_pNpS_gene_level.csv")

cp_path = read_csv("instrain_focus_on_special_pathway/inStrain_pNpS_pathway_tests.csv")
cp_order = read_csv("instrain_focus_on_special_pathway/inStrain_pNpS_order_tests.csv")
cp_thr = read_csv("instrain_focus_on_special_pathway/inStrain_pNpS_threshold_sensitivity.csv")
cp_qc = read_csv("instrain_focus_on_special_pathway/inStrain_pNpS_QC.csv")


# --------------------------------------------------------------------- data quality
checks = []


def check(name, passed, evidence, severity_if_failed="High"):
    require(passed, f"{name}: {evidence}")
    checks.append({"check": name, "status": "PASS", "evidence": evidence,
                   "severity_if_failed": severity_if_failed})


check("PolB required columns",
      {"sample_id", "order", "n_contigs_in_order", "polb_TPM", "total_PolB_TPM",
       "relative_abundance"}.issubset(polb.columns),
      f"{len(polb.columns)} columns present")
check("PolB sample-by-order grain", len(polb) == 167 * 9 and
      polb["sample_id"].nunique() == 167 and polb["order"].nunique() == 9 and
      not polb.duplicated(["sample_id", "order"]).any(),
      f"{len(polb):,} rows = 167 samples × 9 orders; no duplicate sample/order keys")
check("PolB completeness", int(polb.isna().sum().sum()) == 0,
      "no missing cells")
catalogue = (polb.groupby("order", as_index=False)["n_contigs_in_order"]
             .agg(lambda x: x.iloc[0] if x.nunique() == 1 else np.nan)
             .sort_values("n_contigs_in_order", ascending=False))
check("PolB catalogue counts constant", catalogue["n_contigs_in_order"].notna().all(),
      "one invariant catalogue count per order")
check("PolB catalogue total", int(catalogue["n_contigs_in_order"].sum()) == 1095,
      "order counts sum to 1,095; Imitervirales 721 and Algavirales 241")
sample_sums = polb.groupby("sample_id")["relative_abundance"].sum()
check("PolB compositional closure", np.allclose(sample_sums, 1.0, atol=2e-10),
      f"per-sample sum range {sample_sums.min():.12f}–{sample_sums.max():.12f}")

pivot = polb.pivot(index="sample_id", columns="order", values="relative_abundance")
calc_d = ((pivot["Imitervirales"] - pivot["Algavirales"]) /
          (pivot["Imitervirales"] + pivot["Algavirales"]))
spot_key = "sample_id" if "sample_id" in spot.columns else "sample"
spot_d = spot.set_index(spot_key)["D_index"]
common = calc_d.index.intersection(spot_d.index)
check("PolB-to-D join coverage", len(common) == 167,
      f"{len(common)}/167 sample identifiers match SPOT_sample_Dindex.csv")
max_d_diff = float((calc_d.loc[common] - spot_d.loc[common]).abs().max())
check("PolB D-index reproduction", max_d_diff < 1e-9,
      f"maximum absolute difference {max_d_diff:.3g}")

check("pN/pS QC status", set(ag_qc.loc[ag_qc["result"].astype(str).str.upper() == "PASS", "result"]) == {"PASS"}
      and set(cp_qc.loc[cp_qc["result"].astype(str).str.upper() == "PASS", "result"]) == {"PASS"},
      f"all-gene QC rows {len(ag_qc)}; curated-pathway QC rows {len(cp_qc)}")
check("All-gene pN/pS grain", len(ag_gene) == 84846 and
      not ag_gene.duplicated(["genome", "gene"]).any(),
      "84,846 unique genome/gene records")

target = ag_order[(ag_order["scope"] == "Imitervirales vs Algavirales") &
                  (ag_order["unit_of_analysis"].str.startswith("gene"))]
check("All-gene order sensitivity rows", set(target["tier"]) == set("ABCDE"),
      "gene-level Imitervirales–Algavirales result present for tiers A–E")
target_g = ag_order[(ag_order["scope"] == "Imitervirales vs Algavirales") &
                    (ag_order["unit_of_analysis"].str.startswith("genome"))]
check("All-genome order sensitivity rows", set(target_g["tier"]) == set("ABCDE"),
      "genome-level Imitervirales–Algavirales result present for tiers A–E")

check("Network identifiers", nodes["id"].is_unique and edges["id"].is_unique,
      "31 unique nodes and 53 unique edges")
check("Network endpoint integrity", set(edges["source"]).union(edges["target"]) <= set(nodes["id"]),
      "every source and target resolves to the node table")
edge_degree = pd.concat([edges["source"], edges["target"]]).value_counts()
edge_weight = pd.concat([
    edges[["source", "weight"]].rename(columns={"source": "id"}),
    edges[["target", "weight"]].rename(columns={"target": "id"})
]).groupby("id")["weight"].sum()
node_idx = nodes.set_index("id")
check("Network degree reproduction",
      (node_idx.loc[edge_degree.index, "degree"].astype(int) == edge_degree.astype(int)).all(),
      "node degree equals incident edge count")
check("Network weighted-degree reproduction",
      (node_idx.loc[edge_weight.index, "weighted_degree"].astype(float) == edge_weight.astype(float)).all(),
      "node weighted_degree equals sum of incident archived edge weights")
check("Network inference limitation", set(edges["interaction"]) == {"co_occurrence"},
      "aggregate table has no P value, signed FlashWeave weight or known-host flag",
      severity_if_failed="Medium")

check("Rhodopsin metadata completeness", len(rhod) == 183 and rhod.isna().sum().sum() == 0,
      "183 complete metadata rows")
check("Rhodopsin identifiers unique", rhod["Rhodopsin_name"].is_unique,
      "all 183 sequence labels are unique")
check("Rhodopsin topology limitation", set(rhod.columns) ==
      {"Rhodopsin_name", "organism", "domain", "Clade", "TM_count"},
      "metadata include labels/classes/TM counts but no sequences, alignment, tree or structure",
      severity_if_failed="Medium")
check("NCLDV tree metadata", len(tree_tax) == 1620 and tree_tax["MAG"].is_unique and
      tree_tax.isna().sum().sum() == 0,
      "1,620 unique, complete phylogeny-tip annotations")

pd.DataFrame(checks).to_csv(HERE / "new_data_quality_checks.csv", index=False)


# ---------------------------------------------------------------------  workbooks
readme_13 = pd.DataFrame({"note": [
    "Supplementary Table 14. pN/pS results from the deposited inStrain-derived outputs.",
    "pN/pS is defined only when at least one synonymous variant supplies a pS denominator; undefined values are excluded, not set to zero.",
    "All hypothesis tests use per-gene or per-genome summaries, not repeated gene-by-sample rows.",
    "Tier A adds mean gene coverage >=10x; Tier D uses every pN/pS value emitted by the inStrain run, whose site-level settings were minimum coverage 10x and MAF >=0.05.",
    "pN/pS above 1 is compatible with relaxed purifying selection or positive selection and is not proof of adaptation.",
]})
with pd.ExcelWriter(SUPP / "Supplementary_Table_14_pNpS.xlsx", engine="openpyxl") as w:
    readme_13.to_excel(w, sheet_name="README", index=False)
    ag_order.to_excel(w, sheet_name="A_allgene_order_tests", index=False)
    cp_path.to_excel(w, sheet_name="B_curated_pathway_tests", index=False)
    cp_order.to_excel(w, sheet_name="C_curated_order_tests", index=False)
    ag_thr.to_excel(w, sheet_name="D_allgene_thresholds", index=False)
    cp_thr.to_excel(w, sheet_name="E_curated_thresholds", index=False)
    ag_qc.to_excel(w, sheet_name="F_allgene_QC", index=False)
    cp_qc.to_excel(w, sheet_name="G_curated_QC", index=False)

readme_14 = pd.DataFrame({"note": [
    "Supplementary Table 15. Aggregated virus-order/eukaryote-class association network.",
    "Each edge is one order-class pair; weight is the archived number of underlying associations.",
    "These aggregate outputs reproduce the displayed network but do not contain the original phylotype-ASV scores, signs or P values, so the FlashWeave inference itself cannot be rerun from this table.",
    "Edges are ecological associations and are not host assignments.",
]})
with pd.ExcelWriter(SUPP / "Supplementary_Table_15_association_network.xlsx", engine="openpyxl") as w:
    readme_14.to_excel(w, sheet_name="README", index=False)
    nodes.to_excel(w, sheet_name="Nodes", index=False)
    edges.to_excel(w, sheet_name="Edges", index=False)

spot_mask = rhod["Rhodopsin_name"].str.contains("NCLDV|Mirus", case=False, regex=True)
rhod_summary = (rhod.assign(source=np.where(spot_mask, "SPOT giant-virus", "reference"))
                .groupby(["source", "Clade", "TM_count"], as_index=False).size())
polb_qc = pd.DataFrame(checks[:7])
readme_16 = pd.DataFrame({"note": [
    "Supplementary Table 16. SPOT PolB order composition: order-level abundance and catalogue counts.",
    "One row is one sample/order combination. Relative abundances sum to one within numerical tolerance in each of 167 samples.",
    "The invariant order catalogue counts sum to 1,095 PolB phylotypes: 721 Imitervirales, 241 Algavirales and 133 in seven other orders.",
]})
with pd.ExcelWriter(SUPP / "Supplementary_Table_16_SPOT_PolB_order_composition.xlsx", engine="openpyxl") as w:
    readme_16.to_excel(w, sheet_name="README", index=False)
    catalogue.to_excel(w, sheet_name="Catalogue_counts", index=False)
    polb.to_excel(w, sheet_name="Sample_order_abundance", index=False)
    polb_qc.to_excel(w, sheet_name="QC", index=False)

readme_18 = pd.DataFrame({"note": [
    "Supplementary Table 17. Phylogeny metadata: taxonomy of the 1,620 tips and Mirusviricota MCP clade assignment.",
    "The deposited table contains one complete, unique row per tip with order, family and plotting colour.",
    "Tree topology and branch support are supplied separately in the archived Newick tree.",
]})
with pd.ExcelWriter(SUPP / "Supplementary_Table_17_phylogeny_metadata.xlsx", engine="openpyxl") as w:
    readme_18.to_excel(w, sheet_name="README", index=False)
    tree_tax.to_excel(w, sheet_name="NCLDV_taxonomy", index=False)


# --------------------------------------------------------------------- pN/pS figure
# Chart contract: compare order-level pN/pS robustness across reporting tiers and
# show which curated pathway contrasts survive multiplicity correction. Static PNG
# for the manuscript supplement; blue/orange plus marker/line differences remain
# interpretable without color.
order_rows = ag_order[ag_order["scope"] == "Imitervirales vs Algavirales"].copy()
order_rows["unit"] = np.where(order_rows["unit_of_analysis"].str.startswith("gene"), "Gene", "Genome")
order_rows["tier"] = pd.Categorical(order_rows["tier"], list("ABCDE"), ordered=True)
order_rows = order_rows.sort_values(["unit", "tier"])

path_d = cp_path[(cp_path["tier"] == "D") & cp_path["P"].notna()].copy()
path_d = path_d.sort_values("cliffs_delta")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
                     "axes.labelsize": 9, "legend.fontsize": 8})
fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.05, 1.05, 1.65]})

ax = axes[0]
for unit, ls, marker in [("Gene", "-", "o"), ("Genome", "--", "s")]:
    d = order_rows[order_rows["unit"] == unit]
    x = np.arange(len(d))
    ax.plot(x, d["median_a"], color=BLUE, marker=marker, ls=ls, lw=1.5,
            label=f"Imitervirales, {unit.lower()}")
    ax.plot(x, d["median_b"], color=ORANGE, marker=marker, ls=ls, lw=1.5,
            label=f"Algavirales, {unit.lower()}")
ax.set_xticks(range(5), list("ABCDE"))
ax.set_xlabel("Reporting tier")
ax.set_ylabel("Median pN/pS")
ax.set_title("a  Order medians across filters", loc="left", weight="bold")
ax.grid(axis="y", color=LIGHT, lw=0.7)
ax.legend(frameon=False, loc="best")

ax = axes[1]
for unit, offset, marker, fill in [("Gene", -0.08, "o", BLUE), ("Genome", 0.08, "s", "white")]:
    d = order_rows[order_rows["unit"] == unit]
    y = np.arange(len(d)) + offset
    lo = d["cliffs_delta"].to_numpy() - d["cliffs_delta_CI_low"].to_numpy()
    hi = d["cliffs_delta_CI_high"].to_numpy() - d["cliffs_delta"].to_numpy()
    ax.errorbar(d["cliffs_delta"], y, xerr=np.vstack([lo, hi]), fmt=marker,
                mfc=fill, mec=BLUE, ecolor=BLUE, color=BLUE, capsize=2.5, lw=1.2,
                label=unit)
ax.axvline(0, color=INK, lw=0.9)
ax.set_yticks(range(5), list("ABCDE"))
ax.invert_yaxis()
ax.set_xlabel("Cliff's δ (Imitervirales − Algavirales)")
ax.set_title("b  Order effect and 95% CI", loc="left", weight="bold")
ax.grid(axis="x", color=LIGHT, lw=0.7)
ax.legend(frameon=False, loc="lower right")

ax = axes[2]
y = np.arange(len(path_d))
sig = path_d["BH_P"] < 0.05
colors = np.where(sig, ORANGE, GREY)
for yi, (_, row), color in zip(y, path_d.iterrows(), colors):
    ax.plot([row["cliffs_delta_CI_low"], row["cliffs_delta_CI_high"]], [yi, yi], color=color, lw=1.4)
    ax.plot(row["cliffs_delta"], yi, "o", color=color, ms=4.5)
ax.axvline(0, color=INK, lw=0.9)
ax.set_yticks(y, path_d["Pathway"])
ax.set_xlabel("Cliff's δ (pathway − other curated genes)")
ax.set_title("c  Curated pathways, full defined set", loc="left", weight="bold")
ax.grid(axis="x", color=LIGHT, lw=0.7)
ax.text(0.98, 0.02, "orange: BH-adjusted P < 0.05", transform=ax.transAxes,
        ha="right", va="bottom", color=ORANGE, fontsize=8)

fig.suptitle("pN/pS contrasts across orders, reporting tiers and curated pathways", x=0.07, ha="left",
             y=1.01, fontsize=12, weight="bold")
fig.text(0.07, -0.005,
         "Tier A adds mean gene coverage ≥10×; tier D includes every pN/pS value defined by the inStrain run. "
         "Intervals are analytic (order contrasts) or bootstrap (curated pathways).",
         ha="left", va="top", fontsize=8, color="#555555")
fig.tight_layout(rect=[0, 0.035, 1, 0.97])
fig.savefig(FIG / "SupplementaryFig4_pNpS.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)


# --------------------------------------------------------------------- aggregate network figure
# The deposited graph is already aggregated to order/class pairs; the plot therefore
# visualizes that output and does not imply that the original inference was rerun.
G = nx.Graph()
for row in nodes.itertuples(index=False):
    G.add_node(row.id, node_class=row.node_class, weighted_degree=row.weighted_degree)
for row in edges.itertuples(index=False):
    G.add_edge(row.source, row.target, weight=row.weight)

euk = nodes.loc[nodes["node_class"] == "Eukaryote_18S"].sort_values(["weighted_degree", "id"], ascending=[False, True])
vir = nodes.loc[nodes["node_class"] == "GV_order"].sort_values(["weighted_degree", "id"], ascending=[False, True])
pos = {}
for i, name in enumerate(euk["id"]):
    pos[name] = (0.08, 1 - (i + 0.5) / len(euk))
for i, name in enumerate(vir["id"]):
    pos[name] = (0.88, 1 - (i + 0.5) / len(vir))

fig, ax = plt.subplots(figsize=(10.5, 9.2))
for u, v, dat in G.edges(data=True):
    x1, y1 = pos[u]; x2, y2 = pos[v]
    ax.plot([x1, x2], [y1, y2], color="#a2a2a2", alpha=0.45,
            lw=0.45 + 0.55 * math.sqrt(dat["weight"]), zorder=1)
for frame, color, marker in [(euk, TEAL, "o"), (vir, ORANGE, "D")]:
    sizes = 45 + 18 * np.sqrt(frame["weighted_degree"].to_numpy())
    ax.scatter([pos[n][0] for n in frame["id"]], [pos[n][1] for n in frame["id"]],
               s=sizes, c=color, marker=marker, edgecolors="white", linewidths=0.7, zorder=3)
for name in euk["id"]:
    ax.text(pos[name][0] - 0.015, pos[name][1], name, ha="right", va="center", fontsize=7.2)
for name in vir["id"]:
    ax.text(pos[name][0] + 0.015, pos[name][1], name, ha="left", va="center", fontsize=8, weight="bold")
ax.text(0.08, 1.025, "Eukaryotic 18S classes", ha="center", va="bottom", weight="bold", color=TEAL)
ax.text(0.88, 1.025, "Giant-virus orders", ha="center", va="bottom", weight="bold", color=ORANGE)
ax.set_title("Aggregated SPOT virus–eukaryote association network", loc="left", fontsize=12, weight="bold", pad=18)
ax.text(0, -0.035, "31 nodes; 53 order–class edges. Edge width scales with the archived number of underlying associations (1–41).\n"
        "Aggregate outputs do not retain original scores or P values; associations are not host assignments.",
        transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#555555")
ax.set_xlim(-0.22, 1.16); ax.set_ylim(-0.02, 1.05); ax.axis("off")
fig.tight_layout()
fig.savefig(FIG / "SupplementaryFig5_association_network.png", dpi=300,
            bbox_inches="tight", facecolor="white")
plt.close(fig)


# --------------------------------------------------------------------- audit report
allgene_a = target[target["tier"] == "A"].iloc[0]
allgene_d = target[target["tier"] == "D"].iloc[0]
genome_a = target_g[target_g["tier"] == "A"].iloc[0]
genome_d = target_g[target_g["tier"] == "D"].iloc[0]
curated_d = cp_order[(cp_order["tier"] == "D") & (cp_order["scope"] == "all curated genes") &
                     (cp_order["unit_of_analysis"].str.startswith("gene"))].iloc[0]

report = f"""#  new-data audit and restored analyses

Data assessed: `raw_tables/` as available on 21 September 2026.

## Overall assessment

**Share with caveats.** The deposit fully supports the PolB catalogue totals and a reproducible pN/pS analysis. It supports a descriptive aggregate virus–eukaryote network and rhodopsin metadata summary, but not re-estimation of the original network inference or restoration of rhodopsin tree/structure claims.

## Restored analyses

- **SPOT PolB catalogue:** 1,095 phylotypes across nine orders, including 721 Imitervirales and 241 Algavirales. The table is complete at 167 samples × 9 orders, relative abundances close to 1 in every sample, and reproduces sample D with maximum absolute difference {max_d_diff:.3g}.
- **Genome-wide pN/pS:** at tier A (additional mean gene coverage ≥10×), Imitervirales versus Algavirales gene medians are {allgene_a.median_a:.3f} versus {allgene_a.median_b:.3f}, Cliff's δ = {allgene_a.cliffs_delta:.3f} (95% CI {allgene_a.cliffs_delta_CI_low:.3f} to {allgene_a.cliffs_delta_CI_high:.3f}); genome medians are {genome_a.median_a:.3f} versus {genome_a.median_b:.3f}, P = {genome_a.P:.2g}. Without that extra coverage filter (tier D), the effects attenuate to δ = {allgene_d.cliffs_delta:.3f} at gene level and {genome_d.cliffs_delta:.3f} at genome level.
- **Curated functional genes:** the overall order contrast is absent in tier D (medians {curated_d.median_a:.3f} versus {curated_d.median_b:.3f}, P = {curated_d.P:.2f}). TCA-cycle genes have lower pN/pS and ubiquitin-system genes higher than the remaining curated genes after BH correction; no tier-A pathway survives BH correction.
- **Aggregate association network:** 31 nodes and 53 order–class edges; endpoint, degree and weighted-degree checks pass. The graph is descriptive because original phylotype/ASV scores and P values are not present.
- **Rhodopsin metadata:** 183 complete, unique records, of which {int(spot_mask.sum())} carry SPOT NCLDV/Mirus labels. Metadata support clade and transmembrane-count summaries only.

## Not restored

- CCA/Mantel analyses: no phylotype community matrix or 16S/18S abundance matrices were added.
- Rhodopsin phylogeny and AlphaFold structure comparisons: no sequences, alignment, Newick tree, bootstrap output or structure files were added.
- Original FlashWeave inference: the deposited network is already aggregated and has no signed weights or P values.

## Reproducible outputs

- `analysis/new_data_quality_checks.csv`
- `supplementary/Supplementary_Table_14_pNpS.xlsx`
- `supplementary/Supplementary_Table_15_association_network.xlsx`
- `supplementary/Supplementary_Table_16_rhodopsin_metadata.xlsx`
- `supplementary/Supplementary_Table_16_SPOT_PolB_order_composition.xlsx`
- `supplementary/Supplementary_Table_17_phylogeny_metadata.xlsx`
- `supplementary/figures/SupplementaryFig4_pNpS.png`
- `supplementary/figures/SupplementaryFig5_association_network.png`

The pN/pS interpretation is deliberately non-adaptive: values above 1 are compatible with relaxed purifying selection or positive selection and do not establish adaptation.
"""
(REV / "NEW_DATA_AUDIT.md").write_text(report, encoding="utf-8")

print("wrote", HERE / "new_data_quality_checks.csv")
print("wrote  Supplementary Tables 14–18")
print("wrote  Supplementary Figs. 4–5")
print("wrote", REV / "NEW_DATA_AUDIT.md")
