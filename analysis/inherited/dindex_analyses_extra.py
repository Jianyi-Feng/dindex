"""Follow-up checks appended to headline_statistics.csv (run after dindex_analyses.py)."""
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
HERE = Path(__file__).resolve().parent; DATA = HERE.parents[1] / "raw_tables"
hs = pd.read_csv(HERE / "headline_statistics.csv"); hs = hs[hs.section != "G"]
rows = []
def rec(name, value, n=None, lo=None, hi=None, p=None, note=""):
    rows.append(dict(section="G", statistic=name, value=value, n=n, ci_low=lo, ci_high=hi, P=p, note=note))
v = pd.read_csv(HERE / "spot_viral_heatwave_residuals.csv", parse_dates=["date"])
pre = v[v.date < "2014-01-01"]; hw = v[(v.date >= "2014-01-01") & (v.date <= "2015-12-31")]
# 1. absolute SST calibration learned before the event, applied to heatwave SST
cal = smf.ols("D_index ~ SST", data=pre).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
c = cal.conf_int().loc["SST"]
rec("viral pre-event dD/dSST", cal.params.SST, len(pre), c[0], c[1], cal.pvalues.SST)
pred = cal.get_prediction(hw).summary_frame(); pred.index = hw.index
exp_shift = pred["mean"].mean() - cal.fittedvalues.mean()
obs_shift = hw.D_index.mean() - pre.D_index.mean()
rec("heatwave D shift expected from pre-event SST calibration", exp_shift, len(hw))
rec("heatwave D shift observed (raw means)", obs_shift, len(hw))
exceed = hw.D_index - pred["mean"]
w = stats.wilcoxon(exceed)
rec("heatwave D above SST-calibrated expectation, median", exceed.median(), len(hw), p=w.pvalue, note=f"{(exceed>0).sum()}/{len(exceed)} positive; Wilcoxon")
rec("heatwave mean SST minus pre-event mean SST", hw.SST.mean() - pre.SST.mean(), len(hw))
# 2. Tara calibration with temperature + chlorophyll transferred to SPOT
EPS = 1e-3
t = pd.read_csv(DATA / "Tara_observed_Dindex.csv"); t = t[t.D_index_qc == "ok"].copy()
t["femto"] = (t.size_fraction == "Femto-size").astype(int); t["L"] = np.log((t.I_relative + EPS) / (t.A_relative + EPS)); t["lchl"] = np.log10(t.chlorophyll_a)
s = pd.read_csv(DATA / "SPOT_sample_Dindex.csv"); s = s[s.size_fraction_ascii.isin(["0.02-0.2um", "0.2-1.2um"])].dropna(subset=["SST", "Chla"]).copy()
s["femto"] = (s.size_fraction_ascii == "0.02-0.2um").astype(int); s["temperature"] = s.SST; s["lchl"] = np.log10(s.Chla)
for form in ["L ~ temperature + femto", "L ~ temperature + lchl + femto"]:
    m = smf.ols(form, data=t).fit(); d = np.tanh(m.predict(s) / 2)
    r, p = stats.pearsonr(d, s.D_index)
    rec(f"Tara->SPOT transfer [{form}] r", r, len(s), p=p, note=f"bias {np.mean(d - s.D_index):+.3f}; RMSE {np.sqrt(np.mean((d - s.D_index)**2)):.3f}")
rec("SPOT chlorophyll a median (mg m-3)", s.Chla.median(), len(s), note=f"Tara median {t.chlorophyll_a.median():.3f}")
# 3. within-Tara: do slopes differ between warm/cold halves? (linearity of the scale)
for lab, sub in [("<15 degC", t[t.temperature < 15]), (">=15 degC", t[t.temperature >= 15])]:
    m = smf.ols("L ~ temperature + femto", data=sub).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(sub.station_id)[0]})
    c = m.conf_int().loc["temperature"]
    rec(f"Tara dln(I/A)/dT {lab}", m.params.temperature, len(sub), c[0], c[1], m.pvalues.temperature)
# 4. micro/nano fraction: season-adjusted D after the event (2016-2019) vs 2014-2015
mi = pd.read_csv(DATA / "SPOT_sample_Dindex.csv", parse_dates=["date"])
mi = mi[mi.size_fraction_ascii == "1.2-80um"].dropna(subset=["SST"]).copy()
mi["post"] = (mi.heatwave_phase == "post").astype(int); mi["s"] = np.sin(2*np.pi*mi.month/12); mi["c"] = np.cos(2*np.pi*mi.month/12)
fm = smf.ols("D_index ~ post + s + c", data=mi).fit(cov_type="HAC", cov_kwds={"maxlags": 3}); c = fm.conf_int().loc["post"]
rec("micro/nano season-adjusted D, post-event (2016-2019) minus 2014-2015", fm.params.post, len(mi), c[0], c[1], fm.pvalues.post,
    note=f"mean SST anomaly {mi[mi.post==0].SST_anomaly.mean():.2f} -> {mi[mi.post==1].SST_anomaly.mean():.2f} degC")
out = pd.concat([hs, pd.DataFrame(rows)]); out.to_csv(HERE / "headline_statistics.csv", index=False)
print(pd.DataFrame(rows).to_string(float_format=lambda x: f"{x:.4g}"))
