"""Would the Imitervirales-Algavirales pair have been chosen from pre-heatwave data alone?

Addresses the discovery/test dependence disclosed in . Uses only the 89 SPOT metagenomes collected before
2014 (47 viral, 42 picoplankton; the micro/nanoplankton series begins in 2014), and the per-sample order table
deposited on 21 September 2026. Selection rule, fixed before inspecting results: the pair of giant-virus orders
that jointly summarises the largest share of PolB coverage, with its balance D_ab = (a - b)/(a + b) regressed on SST
with size fraction as a covariate (HC3 errors). All 36 pairs of the nine orders are reported.
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "raw_tables"
p = pd.read_csv(DATA / "SPOT_PolB_order_relative_abundance_long.csv")
sst = pd.read_csv(DATA / "SPOT_sample_Dindex.csv")[["sample_id", "SST"]]
w = (p.pivot_table(index=["sample_id", "year", "size_fraction_ascii"], columns="order", values="relative_abundance")
     .reset_index().merge(sst, on="sample_id", how="left"))
pre = w[w.year < 2014].copy()
assert len(pre) == 89 and pre.SST.notna().all()
orders = [c for c in w.columns if c not in ("sample_id", "year", "size_fraction_ascii", "SST")]

rows = []
for a, b in itertools.combinations(sorted(orders), 2):
    # orient every pair so that a positive slope means the first-named order gains with warming
    d = pre[["size_fraction_ascii", "SST", a, b]].copy()
    d = d[(d[a] + d[b]) > 0]
    d["D"] = (d[a] - d[b]) / (d[a] + d[b])
    row = dict(pair=f"{a} vs {b}", n=len(d), median_combined_share=float((pre[a] + pre[b]).median()))
    if len(d) >= 30 and d.D.std() > 0:
        m = smf.ols("D ~ SST + C(size_fraction_ascii)", data=d).fit(cov_type="HC3")
        ci = m.conf_int().loc["SST"]
        row.update(slope_D_per_degC=m.params["SST"], CI_low=ci[0], CI_high=ci[1], P=m.pvalues["SST"],
                   abs_t=abs(m.tvalues["SST"]))
    rows.append(row)
res = pd.DataFrame(rows).sort_values("median_combined_share", ascending=False).reset_index(drop=True)
res["rank_by_share"] = np.arange(1, len(res) + 1)
res["rank_by_abs_t"] = res.abs_t.rank(ascending=False, method="min")
res.to_csv(HERE / "pre2014_pair_derivation.csv", index=False)

top = res.iloc[0]
ia = res[res.pair == "Algavirales vs Imitervirales"].iloc[0]
# report in the paper's orientation, D = (I - A)/(I + A)
summary = pd.DataFrame([
    dict(statistic="pre-2014 samples (viral + picoplankton)", value=len(pre)),
    dict(statistic="pair with the largest combined PolB share", value=top.pair),
    dict(statistic="Imitervirales + Algavirales median share, pre-2014", value=ia.median_combined_share),
    dict(statistic="next-largest combined share (any other pair)", value=res.iloc[1].median_combined_share),
    dict(statistic="D = (I - A)/(I + A) vs SST slope, pre-2014, fraction-adjusted (D per degC)", value=-ia.slope_D_per_degC),
    dict(statistic="slope 95% CI low", value=-ia.CI_high), dict(statistic="slope 95% CI high", value=-ia.CI_low),
    dict(statistic="slope P (HC3)", value=ia.P),
    dict(statistic="rank of the I-A pair by |t| among 36 pairs", value=int(ia.rank_by_abs_t)),
    dict(statistic="combined share of the pair ranked first by |t|", value=res.sort_values("abs_t", ascending=False).iloc[0].median_combined_share),
])
summary.to_csv(HERE / "pre2014_pair_derivation_summary.csv", index=False)
pd.set_option("display.width", 200, "display.max_colwidth", 90)
print(summary.to_string(index=False))
