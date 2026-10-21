"""Supplementary Tables 6 and 8 ( audit).
Table 6: every Tara model x validation design (station-blocked partitions, pooled LOBO, basin-specific LOBO).
Table 8: (A) five displayed traits, Fisher + Wilson, BH within the five; (B) all binary traits tested, BH over the full family;
         (C) marker-recovery-adjusted Mantel-Haenszel results (); (D) prevalence by viral family (partial phylogenetic check).
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUP = HERE.parent / "supplementary"
pkg = ROOT / "marker_recovery_FIXED_package"
rows = []

# ---------------- Table 6
cv = pd.read_csv(HERE / "tara_cross_validation.csv")
LABEL = {"temperature-only D": "Linear D ~ temperature + size fraction (main-text temperature model; Fig. 3d)",
         "temperature-only log-ratio": "Linear ln(I/A) ~ temperature + size fraction (sensitivity model)",
         "ridge log-ratio (7 predictors)": "Ridge ln(I/A), 7 environmental predictors + size fraction",
         "random forest log-ratio (7 predictors)": "Random forest ln(I/A), 7 environmental predictors + size fraction",
         "random forest D (7 predictors)": "Random forest D, 7 environmental predictors + size fraction"}
cv["model_description"] = cv.model.map(LABEL)
sb = cv[cv.scheme == "5-fold station-blocked"].groupby("model_description")[["R2", "RMSE", "MAE"]].agg(["mean", "min", "max"])
sb.columns = [f"{a}_{b}" for a, b in sb.columns]
lobo = cv[cv.scheme == "leave-one-basin-out"].set_index("model_description")[["R2", "RMSE", "MAE"]].add_prefix("LOBO_pooled_")
t6 = sb.join(lobo).reset_index()
basin = cv[cv.scheme.str.startswith("LOBO held-out basin")].assign(basin=lambda d: d.scheme.str.replace("LOBO held-out basin: ", ""))
basin = basin[["model_description", "basin", "n", "R2", "RMSE", "MAE", "bias"]]
ensv = pd.read_csv(HERE / "ensemble_validation.csv").iloc[0]
ens = pd.DataFrame([dict(model_description="Final 14-model ensemble, 7 predictors (I and A modelled separately), nested station-grouped CV",
                         R2_mean=ensv.R2, RMSE_mean=ensv.RMSE, MAE_mean=ensv.MAE)])
t6 = pd.concat([ens, t6], ignore_index=True)
with pd.ExcelWriter(SUP / "Supplementary_Table_5_Tara_model_validation.xlsx") as xw:
    pd.DataFrame({"note": ["Station-blocked: fivefold cross-validation over 20 random partitions of the 83 stations (mean, min, max across partitions).",
                           "LOBO pooled: R2 computed over the pooled held-out predictions from all seven withheld basins.",
                           "Basin-specific LOBO metrics are descriptive; R2 within a basin can be negative when within-basin variance is small.",
                           "The field labelled CO3 was removed after the source-registry audit showed that it is 19'-butanoyloxyfucoxanthin, not carbonate ion.",
                           "n = 196 metagenomes, 83 stations."]}).to_excel(xw, sheet_name="README", index=False)
    t6.to_excel(xw, sheet_name="A_summary", index=False)
    basin.to_excel(xw, sheet_name="B_basin_specific_LOBO", index=False)
print(t6.round(3).to_string())

# ---------------- Table 8
g = pd.read_csv(pkg / "Supplementary_Data_7_GVMAG_traits_marker_recovery.csv")
g = g[g.order.isin(["Imitervirales", "Algavirales"])].copy()
nI, nA = (g.order == "Imitervirales").sum(), (g.order == "Algavirales").sum()
skip = {"has_any_gene_annotation", "genome_id_in_gene_file"}
binary = [c for c in g.columns if not c.endswith("_copies") and c not in skip and set(pd.unique(g[c].dropna())) <= {0, 1, True, False}
          and c not in {"marker_recovery_score"}]
binary = [c for c in binary if not c.startswith(("PolB", "A32", "VLTF3", "SFII", "TFIIB", "TopoII", "MCP", "RNAP", "TFIIS"))]
res = []
for c in binary:
    kI, kA = int(g.loc[g.order == "Imitervirales", c].sum()), int(g.loc[g.order == "Algavirales", c].sum())
    orr, p = stats.fisher_exact([[kI, nI - kI], [kA, nA - kA]])
    liI, hiI = proportion_confint(kI, nI, method="wilson"); liA, hiA = proportion_confint(kA, nA, method="wilson")
    res.append(dict(trait=c, Imitervirales_pos=kI, Imitervirales_n=nI, prev_I=kI / nI, prev_I_lo=liI, prev_I_hi=hiI,
                    Algavirales_pos=kA, Algavirales_n=nA, prev_A=kA / nA, prev_A_lo=liA, prev_A_hi=hiA, fisher_OR=orr, P=p))
allt = pd.DataFrame(res)
allt["P_BH_all_traits"] = multipletests(allt.P, method="fdr_bh")[1]
DISPLAY = {"any_HSP": "Heat-shock proteins", "pw_DNA_repair": "DNA repair", "type_I_rhodopsin": "Type-I rhodopsin",
           "pw_glycolysis": "Glycolysis", "pw_TCA": "TCA cycle"}
five = allt[allt.trait.isin(DISPLAY)].copy()
five["display_label"] = five.trait.map(DISPLAY)
five["P_BH_five_displayed"] = multipletests(five.P, method="fdr_bh")[1]
mh = pd.read_csv(pkg / "marker_recovery_adjusted_genome_comparisons.csv")
fam_rows = []
for fam, gf in g[g.order == "Imitervirales"].groupby("family"):
    r = dict(order="Imitervirales", family=fam, n=len(gf))
    r.update({DISPLAY[t]: gf[t].mean() for t in DISPLAY}); fam_rows.append(r)
for fam, gf in g[g.order == "Algavirales"].groupby("family"):
    r = dict(order="Algavirales", family=fam, n=len(gf))
    r.update({DISPLAY[t]: gf[t].mean() for t in DISPLAY}); fam_rows.append(r)
fam = pd.DataFrame(fam_rows)
with pd.ExcelWriter(SUP / "Supplementary_Table_7_GVMAG_functional_contrasts.xlsx") as xw:
    pd.DataFrame({"note": [f"Focal genomes: {nI} Imitervirales and {nA} Algavirales GVMAGs (Supplementary Data 7).",
                           "A: the five traits displayed in Fig. 4a; two-sided Fisher exact tests, Wilson 95% intervals, BH across the five displayed traits.",
                           f"B: all {len(allt)} binary functional traits and pathway categories in Supplementary Data 7, BH across the full tested family.",
                           "C: Mantel-Haenszel common odds ratios stratified by six-marker recovery score (strata 4-6 shared by both orders).",
                           "D: prevalence of the displayed traits within each viral family (families with n >= 1); descriptive partial check for clade dependence."]}
                 ).to_excel(xw, sheet_name="README", index=False)
    five.to_excel(xw, sheet_name="A_displayed_traits", index=False)
    allt.sort_values("P").to_excel(xw, sheet_name="B_all_traits_BH", index=False)
    mh.to_excel(xw, sheet_name="C_marker_recovery_adjusted", index=False)
    fam.to_excel(xw, sheet_name="D_by_family", index=False)
print("binary traits tested:", len(allt)); print(five[["display_label", "Imitervirales_pos", "Algavirales_pos", "P", "P_BH_five_displayed", "P_BH_all_traits"]].to_string())
print("max BH (all traits) among displayed:", five.P_BH_all_traits.max())
big = fam[fam.n >= 5]; print(big.round(2).to_string())
pd.DataFrame([dict(statistic="binary traits tested", value=len(allt)),
              dict(statistic="max BH-adjusted P among five displayed traits, BH over all tested traits", value=five.P_BH_all_traits.max()),
              dict(statistic="Imitervirales families with n>=5", value=int((big.order == "Imitervirales").sum())),
              dict(statistic="Imitervirales families (n>=5) with HSP prevalence > max Algavirales family prevalence",
                   value=int((big[big.order == "Imitervirales"]["Heat-shock proteins"] > big[big.order == "Algavirales"]["Heat-shock proteins"].max()).sum()))]
             ).to_csv(HERE / "genome_table_statistics.csv", index=False)
