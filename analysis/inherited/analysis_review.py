""" review additions and corrections for the D-index analysis.

This script does not replace the full / pipeline. It corrects labels/estimates that
were stale in  and adds reviewer-facing robustness summaries used in the  text.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "raw_tables"

# Start from the freshly rerun  headline table if present.
base_path = HERE / "headline_statistics.csv"
hs = pd.read_csv(base_path)

# Remove stale -managed rows if re-running.
remove_patterns = [
    "micro/nano season-adjusted D",
    "Tara log-ratio slope pseudocount",
    "SPOT vs Tara overlapping-range slope difference",
]
mask = np.ones(len(hs), dtype=bool)
for pat in remove_patterns:
    mask &= ~hs.statistic.str.contains(pat, na=False)
hs = hs[mask].copy()

new = []
def rec(section, statistic, value, n=None, ci_low=None, ci_high=None, P=None, note=""):
    new.append(dict(section=section, statistic=statistic, value=value, n=n,
                    ci_low=ci_low, ci_high=ci_high, P=P, note=note))

# 1) Correct the micro/nano follow-up contrast and its label.
s = pd.read_csv(DATA / "SPOT_sample_Dindex.csv", parse_dates=["date"])
mi = s[(s.size_fraction_ascii == "1.2-80um") & s.SST.notna()].copy()
mi["post"] = (mi.heatwave_phase == "post").astype(int)
mi["sin12"] = np.sin(2*np.pi*mi.month/12)
mi["cos12"] = np.cos(2*np.pi*mi.month/12)
fm = smf.ols("D_index ~ post + sin12 + cos12", data=mi).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
ci = fm.conf_int().loc["post"]
g0 = mi[mi.post == 0]; g1 = mi[mi.post == 1]
rec("G", "micro/nano season-adjusted D, post phase minus heatwave phase",
    fm.params["post"], n=len(mi), ci_low=ci[0], ci_high=ci[1], P=fm.pvalues["post"],
    note=(f"heatwave-phase matched-SST samples {g0.date.min().date()}..{g0.date.max().date()} (n={len(g0)}), "
          f"post-phase {g1.date.min().date()}..{g1.date.max().date()} (n={len(g1)}); "
          f"mean SST anomaly {g0.SST_anomaly.mean():.2f} -> {g1.SST_anomaly.mean():.2f} degC"))

# 2) Pseudocount sensitivity for Tara log-ratio analyses involving zero order abundances.
t = pd.read_csv(DATA / "Tara_observed_Dindex.csv")
t = t[t.D_index_qc == "ok"].copy()
t["femto"] = (t.size_fraction == "Femto-size").astype(int)
for eps in (1e-4, 1e-3, 1e-2):
    tt = t.copy()
    tt["L"] = np.log((tt.I_relative + eps)/(tt.A_relative + eps))
    m = smf.ols("L ~ temperature + femto", data=tt).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(tt.station_id)[0]})
    ci = m.conf_int().loc["temperature"]
    rec("H", f"Tara log-ratio slope pseudocount {eps:g}", m.params["temperature"], n=len(tt),
        ci_low=ci[0], ci_high=ci[1], P=m.pvalues["temperature"],
        note="ln[(I+eps)/(A+eps)] ~ temperature + size fraction; station-clustered errors")

# 3) Approximate independent-slope comparison for the overlapping 13-24 C range.
#    This is a sensitivity check only; the manuscript reports overlapping intervals rather than treating this as a formal equivalence test.
spot_row = hs[hs.statistic == "SPOT dD/dSST, fraction-adjusted (HC3)"].iloc[0]
tara_row = hs[hs.statistic == "Tara dD/dT within SPOT SST range 13-24 degC"].iloc[0]
se_spot = (spot_row.ci_high - spot_row.ci_low)/(2*1.96)
se_tara = (tara_row.ci_high - tara_row.ci_low)/(2*1.96)
z = (spot_row.value - tara_row.value)/np.sqrt(se_spot**2 + se_tara**2)
p = 2*stats.norm.sf(abs(z))
rec("H", "SPOT vs Tara overlapping-range slope difference (approx z)",
    spot_row.value - tara_row.value, n=int(spot_row.n + tara_row.n), P=p,
    note="independent-estimate sensitivity check from reported 95% intervals; not an equivalence test")

out = pd.concat([hs, pd.DataFrame(new)], ignore_index=True)
out.to_csv(HERE / "headline_statistics.csv", index=False)

# 4) Basin-specific leave-one-basin-out metrics, retained to prevent pooled R2 from being overinterpreted.
cv = pd.read_csv(HERE / "tara_cross_validation.csv") if (HERE / "tara_cross_validation.csv").exists() else pd.read_csv(HERE / "tara_cross_validation.csv")
lobo = cv[cv.scheme.str.startswith("LOBO held-out basin:")].copy()
lobo["basin"] = lobo.scheme.str.replace("LOBO held-out basin: ", "", regex=False)
lobo[["basin","model","n","R2","RMSE","MAE","bias"]].to_csv(HERE / "tara_lobo_basin_metrics.csv", index=False)

# 5) A concise machine-readable review table.
review = pd.DataFrame([
    dict(item="SPOT fraction-slope heterogeneity", status="supported", detail="joint SST x fraction interaction P=0.73; use 'similar', not 'identical' slopes"),
    dict(item="heatwave event inference", status="revise emphasis", detail="year-level exact permutation P=0.048 should be primary inferential check; sample-level P=0.0003 is complementary"),
    dict(item="micro/nano relaxation", status="corrected", detail=f"post phase is 2017-2018 matched-SST samples vs 2014-2016 heatwave phase; beta={fm.params['post']:.3f}, P={fm.pvalues['post']:.3f}"),
    dict(item="Tara leave-one-basin-out", status="revise wording", detail="pooled LOBO R2 summarizes all held-out predictions; basin-specific R2 varies and is often negative, so do not claim uniform basin-level skill"),
    dict(item="permutation importance", status="wording", detail="report descriptive permutation distribution across held-out folds/partitions; do not label pooled permutation quantiles as a confidence interval"),
    dict(item="CMIP6 front uncertainty", status="revise emphasis", detail="use four-ESM model range as primary projection spread; fixed-effect CI is not a full climate-projection uncertainty interval"),
])
review.to_csv(HERE / "analysis_review_items.csv", index=False)
print(out.tail(8).to_string(index=False))
print("\nWrote", HERE / "headline_statistics.csv")
print("Wrote", HERE / "tara_lobo_basin_metrics.csv")
print("Wrote", HERE / "analysis_review_items.csv")
