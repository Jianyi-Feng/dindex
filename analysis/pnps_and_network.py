""" correction of the pN/pS order contrasts, and rebuilt Supplementary Figs 4-5.

Two defects in the / pN/pS analysis are corrected here.

1. Tier A was described as "mean gene coverage >= 10x". In the deposited data it is a filter on individual
   gene-by-sample observations (coverage >= 10x); `passes_tier_A` equals `coverage >= 10` in every row.
2. The genome-level order tests used the median of all pooled gene-by-sample observations in a genome
   (`median_pNpS_*` in the genome table). That weights each gene by how often it was sampled, which is the
   repeated-measures weighting the Methods say was avoided. Here a genome is summarised as the median of its
   per-gene medians, as the Methods describe. The pooled variant is kept as a labelled sensitivity row.

Gene-level Cliff's delta intervals use a genome-cluster bootstrap, because genes are nested in genomes;
genome-level intervals use an ordinary bootstrap over genomes.
Figures follow the house palette: Imitervirales #0072B2, Algavirales #D55E00.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tools"))
import paths as P  # noqa: E402  (repository layout is defined in tools/paths.py)

FIG = P.outdir(P.OUT / "figures")   # was ROOT/"supplementary"/"figures"
rng = np.random.default_rng(20260922)

I_C, A_C, GREY, LGREY, INK = "#0072B2", "#D55E00", "#7F7F7F", "#D9D9D9", "#222222"
TEAL = "#5AA6A6"
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "savefig.dpi": 600,
})
TIERS = {"A": ("passes_tier_A", "observation coverage >= 10x"),
         "B": ("passes_tier_B", "observation coverage >= 5x"),
         "C": ("passes_tier_C", "observation coverage >= 1x"),
         "D": (None, "no coverage filter (every defined pN/pS)"),
         "E": ("passes_tier_E", "observation coverage >= 10x and breadth >= 0.5")}


def cliff(a, b):
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2 * u / (len(a) * len(b)) - 1


# --------------------------------------------------------------------- data and definition check
cols = ["genome", "order", "gene", "coverage", "pNpS_variants", "passes_tier_A", "passes_tier_B",
        "passes_tier_C", "passes_tier_E"]
obs = pd.read_csv(P.resolve("data/pnps_network/instrain_all_genes/inStrain_allgenes_pNpS_gene_by_sample_defined.csv"), usecols=cols)
obs = obs[obs.order.isin(["Imitervirales", "Algavirales"])]
tier_a_is_observation_filter = bool((obs.passes_tier_A == (obs.coverage >= 10)).all())
assert tier_a_is_observation_filter, "tier A is not a per-observation coverage filter"

rows = []
for tier, (flag, label) in TIERS.items():
    sub = obs if flag is None else obs[obs[flag]]
    gene = sub.groupby(["order", "genome", "gene"]).pNpS_variants.median().reset_index()
    genome = gene.groupby(["order", "genome"]).pNpS_variants.median().reset_index()
    pooled = sub.groupby(["order", "genome"]).pNpS_variants.median().reset_index()

    # gene unit, with a genome-cluster bootstrap for the interval
    gI, gA = (gene[gene.order == o] for o in ("Imitervirales", "Algavirales"))
    by_genome = {o: {g: d.pNpS_variants.to_numpy() for g, d in gene[gene.order == o].groupby("genome")}
                 for o in ("Imitervirales", "Algavirales")}
    keys = {o: list(v) for o, v in by_genome.items()}
    boot = []
    for _ in range(1000):
        a = np.concatenate([by_genome["Imitervirales"][g] for g in rng.choice(keys["Imitervirales"], len(keys["Imitervirales"]))])
        b = np.concatenate([by_genome["Algavirales"][g] for g in rng.choice(keys["Algavirales"], len(keys["Algavirales"]))])
        boot.append(cliff(a, b))
    rows.append(dict(tier=tier, filter=label, unit="gene (per-gene median across samples)",
                     n_I=len(gI), n_A=len(gA), median_I=gI.pNpS_variants.median(), median_A=gA.pNpS_variants.median(),
                     P=stats.mannwhitneyu(gI.pNpS_variants, gA.pNpS_variants).pvalue,
                     cliffs_delta=cliff(gI.pNpS_variants, gA.pNpS_variants),
                     CI_low=np.percentile(boot, 2.5), CI_high=np.percentile(boot, 97.5),
                     CI_method="genome-cluster bootstrap (1,000)", role="secondary: genes are nested in genomes"))

    for agg, table, role in [("genome (median of per-gene medians)", genome, "primary"),
                             ("genome (median of pooled observations)", pooled,
                              "sensitivity: weights genes by sampling frequency; gene-weighted sensitivity")]:
        a = table[table.order == "Imitervirales"].pNpS_variants.to_numpy()
        b = table[table.order == "Algavirales"].pNpS_variants.to_numpy()
        boot = [cliff(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(4000)]
        rows.append(dict(tier=tier, filter=label, unit=agg, n_I=len(a), n_A=len(b),
                         median_I=np.median(a), median_A=np.median(b), P=stats.mannwhitneyu(a, b).pvalue,
                         cliffs_delta=cliff(a, b), CI_low=np.percentile(boot, 2.5), CI_high=np.percentile(boot, 97.5),
                         CI_method="genome bootstrap (4,000)", role=role))

res = pd.DataFrame(rows)
res.to_csv(P.outdir(P.OUT / "tables") / "pnps_order_tests.csv", index=False)  # was HERE/"pnps_order_tests.csv"

# --------------------------------------------------------------------- Supplementary Table 14 ()
SUPP_XLS = pd.ExcelFile(P.resolve("archive/supplementary/Supplementary_Table_14_pNpS.xlsx"))  # input sheet
readme = pd.DataFrame({"note": [
    "Supplementary Table 14. pN/pS results from the deposited inStrain-derived outputs ( corrected order tests).",
    "pN/pS is defined only when at least one synonymous variant supplies a pS denominator; undefined values are excluded, not set to zero.",
    "Reporting tiers filter individual gene-by-sample observations: A, coverage >= 10x; B, >= 5x; C, >= 1x; D, no filter; "
    "E, coverage >= 10x and breadth >= 0.5. The inStrain run itself used minimum site coverage 10x and MAF >= 0.05.",
    "Sheet A (): Imitervirales versus Algavirales. The primary genome unit is the median of per-gene medians; the pooled-"
    "observation genome median (gene-weighted sensitivity) is retained as a sensitivity row because it weights genes by sampling frequency.",
    "Gene-level intervals use a genome-cluster bootstrap because genes are nested in genomes; genome-level intervals use a genome bootstrap.",
    "Sheets B-G are unchanged from  (curated pathway tests, thresholds and QC); pathway tests are gene-level and descriptive.",
    "pN/pS above 1 is compatible with relaxed purifying selection or positive selection and is not proof of adaptation."]})
with pd.ExcelWriter(P.outdir(P.OUT / "supplementary_tables") / "Supplementary_Table_14_pNpS.xlsx", engine="openpyxl") as w:
    readme.to_excel(w, sheet_name="README", index=False)
    res.to_excel(w, sheet_name="A_order_tests", index=False)
    for s in SUPP_XLS.sheet_names[2:]:
        SUPP_XLS.parse(s).to_excel(w, sheet_name=s, index=False)
    if "A_allgene_order_tests" in SUPP_XLS.sheet_names:  # sheet of the predecessor workbook
        SUPP_XLS.parse("A_allgene_order_tests").to_excel(w, sheet_name="supplementary_allgene_order_tests", index=False)
    else:
        print("note: the deposited Supplementary_Table_14_pNpS.xlsx has no 'A_allgene_order_tests' sheet "
              "(the superseded predecessor workbook is not part of the deposit); "
              "sheet 'supplementary_allgene_order_tests' is not rewritten.")

# --------------------------------------------------------------------- Supplementary Fig. 4 ()
cp = SUPP_XLS.parse("B_curated_pathway_tests")
path_d = cp[(cp.tier == "D") & cp.P.notna()].sort_values("cliffs_delta")
fig, axes = plt.subplots(1, 3, figsize=(240 * MM, 88 * MM), gridspec_kw=dict(width_ratios=[1, 1.1, 1.45], wspace=0.55))
x = np.arange(5)
prim = res[res.unit.str.startswith("genome (median of per")].set_index("tier").loc[list("ABCDE")]
gene = res[res.unit.str.startswith("gene")].set_index("tier").loc[list("ABCDE")]

ax = axes[0]
for d, ls, mk, lab in [(prim, "-", "s", "genome"), (gene, "--", "o", "gene")]:
    ax.plot(x, d.median_I, color=I_C, ls=ls, marker=mk, ms=3.5, lw=1.1, label=f"Imitervirales, {lab}")
    ax.plot(x, d.median_A, color=A_C, ls=ls, marker=mk, ms=3.5, lw=1.1, label=f"Algavirales, {lab}")
ax.set_xticks(x, list("ABCDE")); ax.set_xlabel("Reporting tier"); ax.set_ylabel("Median pN/pS")
ax.set_title("Order medians across filters", loc="left", fontweight="bold")
ax.text(-0.22, 1.06, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax.legend(loc="upper right", fontsize=5.8, labelspacing=0.25, handlelength=2.2)

ax = axes[1]
for d, dy, mk, face, lab in [(prim, -0.12, "s", I_C, "Genome (primary)"), (gene, 0.12, "o", "white", "Gene (nested)")]:
    y = x + dy
    ax.errorbar(d.cliffs_delta, y, xerr=[d.cliffs_delta - d.CI_low, d.CI_high - d.cliffs_delta], fmt=mk, ms=3.8,
                mfc=face, mec=I_C, ecolor=I_C, elinewidth=0.9, capsize=1.8, label=lab)
ax.axvline(0, color=INK, lw=0.7)
ax.set_yticks(x, list("ABCDE")); ax.invert_yaxis(); ax.set_ylabel("Reporting tier")
ax.set_xlabel("Cliff's δ, Imitervirales − Algavirales (95% CI)")
ax.set_title("Order effect by analysis unit", loc="left", fontweight="bold")
ax.text(-0.22, 1.06, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=6)

ax = axes[2]
y = np.arange(len(path_d))
for yi, (_, r) in zip(y, path_d.iterrows()):
    sig = r.BH_P < 0.05
    ax.plot([r.cliffs_delta_CI_low, r.cliffs_delta_CI_high], [yi, yi], color=INK if sig else GREY, lw=0.9)
    ax.plot(r.cliffs_delta, yi, "o", ms=4, mfc=INK if sig else "white", mec=INK if sig else GREY, mew=0.8)
ax.axvline(0, color=INK, lw=0.7)
ax.set_yticks(y, path_d.Pathway)
ax.set_xlabel("Cliff's δ, pathway − other curated genes (95% CI)")
ax.set_title("Curated pathways, tier D (gene level)", loc="left", fontweight="bold")
ax.text(-0.5, 1.06, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")
ax.legend(handles=[Line2D([], [], marker="o", color=INK, mfc=INK, lw=0.9, ms=4, label="BH-adjusted P < 0.05"),
                   Line2D([], [], marker="o", color=GREY, mfc="white", lw=0.9, ms=4, label="not significant")],
          loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=6)
fig.savefig(FIG / "SupplementaryFig4_pNpS.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "SupplementaryFig4_pNpS.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)

# --------------------------------------------------------------------- Supplementary Fig. 5 ()
nodes = pd.read_csv(P.resolve("data/SPOT_18S_GV_network_nodes.csv"))
# the deposited CSV names the eukaryote side of an edge 'eukaryote_class'; the published
# association network (Supplementary Table 15, sheet Edges) names the same column 'source'.
# The rows and values are identical (verified); only the header differs.
edges = pd.read_csv(P.resolve("data/SPOT_18S_GV_network_edges.csv")).rename(columns={"eukaryote_class": "source"})
euk = nodes[nodes.node_class == "Eukaryote_18S"].sort_values("weighted_degree", ascending=False)
gv = nodes[nodes.node_class == "GV_order"].sort_values("weighted_degree", ascending=False)
pos = {n: (0.1, 1 - i / (len(euk) - 1)) for i, n in enumerate(euk.id)}
pos.update({n: (0.9, 0.92 - i * 0.84 / (len(gv) - 1)) for i, n in enumerate(gv.id)})
GV_COL = {"Imitervirales": I_C, "Algavirales": A_C}
wmax = edges.weight.max()
fig, ax = plt.subplots(figsize=(180 * MM, 175 * MM))
for r in edges.itertuples(index=False):
    (x1, y1), (x2, y2) = pos[r.source], pos[r.target]
    ax.plot([x1, x2], [y1, y2], color=GREY, alpha=0.35, lw=0.4 + 3.2 * r.weight / wmax, zorder=1)
size = lambda wd: 18 + 260 * wd / nodes.weighted_degree.max()
for r in euk.itertuples(index=False):
    ax.scatter(*pos[r.id], s=size(r.weighted_degree), color=TEAL, edgecolor="white", lw=0.6, zorder=3)
    ax.text(pos[r.id][0] - 0.018, pos[r.id][1], r.id, ha="right", va="center", fontsize=6)
for r in gv.itertuples(index=False):
    ax.scatter(*pos[r.id], s=size(r.weighted_degree), marker="D", color=GV_COL.get(r.id, GREY),
               edgecolor="white", lw=0.6, zorder=3)
    gap = 0.012 + 0.0022 * np.sqrt(size(r.weighted_degree))   # clear the diamond, which scales with degree
    ax.text(pos[r.id][0] + gap, pos[r.id][1], r.id, ha="left", va="center", fontsize=7, fontweight="bold")
ax.text(0.1, 1.045, "Eukaryotic 18S classes", ha="center", fontsize=7.5, fontweight="bold", color=INK)
ax.text(0.9, 1.045, "Giant-virus orders", ha="center", fontsize=7.5, fontweight="bold", color=INK)
leg = [Line2D([], [], color=GREY, alpha=0.6, lw=0.4 + 3.2 * w / wmax, label=f"{w} association" + ("s" if w != 1 else "")) for w in (1, 10, wmax)]
leg += [Line2D([], [], marker="o", lw=0, color=TEAL, ms=np.sqrt(size(v)) * 0.8,
               label=f"weighted degree {v}") for v in (10, 100)]
ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.62, -0.02), ncol=2, fontsize=6, labelspacing=0.9)
ax.set_xlim(-0.12, 1.12); ax.set_ylim(-0.1, 1.08); ax.axis("off")
fig.savefig(FIG / "SupplementaryFig5_association_network.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "SupplementaryFig5_association_network.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)

pd.set_option("display.width", 200)
print("tier A is a per-observation coverage filter:", tier_a_is_observation_filter)
print(res[["tier", "unit", "n_I", "n_A", "median_I", "median_A", "P", "cliffs_delta", "CI_low", "CI_high"]]
      .round(4).to_string(index=False))
print("weighted degree range:", nodes.weighted_degree.min(), "-", nodes.weighted_degree.max())
