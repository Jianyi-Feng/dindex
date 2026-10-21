"""Independent re-analysis of the D-index source tables.

Reads only the nine CSV tables in ../raw_tables and writes every statistic used in the
 manuscript to headline_statistics.csv (one row per number) plus supporting tables.

Analyses
  A  SPOT calibration: D vs SST, seasonal vs anomaly components, SST at D = 0
  B  Compositional identity: D = tanh(ln(I/A)/2)
  C  Heatwave: harmonic residuals, year-block inference, matched months,
     out-of-sample prediction of the heatwave shift from the pre-event anomaly slope
  D  Tara: station-level calibration, station-blocked and basin-blocked CV,
     out-of-fold permutation importance
  E  Cross-system transfer: Tara temperature calibration applied to SPOT
  F  CMIP6 front: per-ESM rates, signal-to-noise per decade, time of emergence
"""
from pathlib import Path
import itertools, json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "raw_tables"
OUT = HERE
RNG = np.random.default_rng(20260916)
KM_PER_DEG = 111.32
EPS = 1e-3  # pseudocount (relative-abundance units) for log-ratios with zeros

rows = []


def rec(section, name, value, n=None, ci_low=None, ci_high=None, p=None, note=""):
    rows.append(dict(section=section, statistic=name, value=value, n=n,
                     ci_low=ci_low, ci_high=ci_high, P=p, note=note))


spot = pd.read_csv(DATA / "SPOT_sample_Dindex.csv", parse_dates=["date"])
tara = pd.read_csv(DATA / "Tara_observed_Dindex.csv")
val = pd.read_csv(DATA / "Dindex_model_validation.csv")
latyr = pd.read_csv(DATA / "CMIP6_Dindex_latitude_year.csv")
front = pd.read_csv(DATA / "Dindex_front_trajectory.csv")
heat = pd.read_csv(DATA / "SPOT_heatwave_Dindex_stats.csv")

FRAC = {"0.02-0.2um": "viral", "0.2-1.2um": "pico", "1.2-80um": "micro"}
spot["frac"] = spot.size_fraction_ascii.map(FRAC)
spot["I"] = spot.Imitervirales_relative_abundance
spot["A"] = spot.Algavirales_relative_abundance
spot["L"] = np.log(spot.I / spot.A)
spot["SST_clim"] = spot.SST - spot.SST_anomaly

# ---------------------------------------------------------------- A. SPOT calibration
rec("A", "SPOT metagenomes", len(spot), note=json.dumps(spot.frac.value_counts().to_dict()))
rec("A", "SPOT I+A share of PolB coverage, median", spot.I_plus_A.median(), n=len(spot),
    ci_low=spot.I_plus_A.quantile(.25), ci_high=spot.I_plus_A.quantile(.75), note="IQR")
for f, g in spot.groupby("frac"):
    rec("A", f"I+A share median ({f})", g.I_plus_A.median(), n=len(g))

g = spot.dropna(subset=["SST"]).copy()
r, p = stats.pearsonr(g.SST, g.D_index)
rec("A", "SPOT D vs SST Pearson r", r, n=len(g), p=p)
rho, prho = stats.spearmanr(g.SST, g.D_index)
rec("A", "SPOT D vs SST Spearman rho", rho, n=len(g), p=prho)
m1 = smf.ols("D_index ~ SST + C(frac)", data=g).fit(cov_type="HC3")
ci = m1.conf_int().loc["SST"]
rec("A", "SPOT dD/dSST, fraction-adjusted (HC3)", m1.params.SST, n=len(g), ci_low=ci[0], ci_high=ci[1],
    p=m1.pvalues.SST, note="D units per degC")

# formal test for heterogeneity of the SST slope among size fractions
m_int = smf.ols("D_index ~ SST * C(frac)", data=g).fit(cov_type="HC3")
wald_int = m_int.f_test("SST:C(frac)[T.pico] = 0, SST:C(frac)[T.viral] = 0")
rec("A", "SPOT SST slope heterogeneity across fractions P", float(wald_int.fvalue), n=len(g),
    p=float(wald_int.pvalue), note="joint Wald/F test of SST x fraction interactions; statistic stored in value")

# seasonal (climatological) and anomaly components of SST
m2 = smf.ols("D_index ~ SST_clim + SST_anomaly + C(frac)", data=g).fit(cov_type="HC3")
for term in ["SST_clim", "SST_anomaly"]:
    c = m2.conf_int().loc[term]
    rec("A", f"SPOT dD/d{term}, fraction-adjusted (HC3)", m2.params[term], n=len(g),
        ci_low=c[0], ci_high=c[1], p=m2.pvalues[term])
wald = m2.t_test("SST_clim - SST_anomaly = 0")
rec("A", "SPOT seasonal minus anomaly slope", float(wald.effect[0]), n=len(g),
    ci_low=float(wald.conf_int()[0][0]), ci_high=float(wald.conf_int()[0][1]), p=float(wald.pvalue),
    note="equal slopes => one temperature sensitivity for season and anomaly")

# per fraction, time-ordered HAC errors
for f, gg in g.sort_values("date").groupby("frac"):
    mf = smf.ols("D_index ~ SST", data=gg).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
    c = mf.conf_int().loc["SST"]
    sst0 = -mf.params.Intercept / mf.params.SST
    boot = []
    for _ in range(2000):
        bb = gg.sample(len(gg), replace=True, random_state=RNG.integers(1e9))
        pb = np.polyfit(bb.SST, bb.D_index, 1)
        boot.append(-pb[1] / pb[0])
    rec("A", f"dD/dSST ({f}, HAC lag 3)", mf.params.SST, n=len(gg), ci_low=c[0], ci_high=c[1], p=mf.pvalues.SST)
    rec("A", f"SST at D=0 ({f})", sst0, n=len(gg), ci_low=np.percentile(boot, 2.5),
        ci_high=np.percentile(boot, 97.5), note="bootstrap 2000")

# ---------------------------------------------------------------- B. compositional identity
ident = np.max(np.abs(spot.D_index - np.tanh(spot.L / 2)))
rec("B", "max |D - tanh(ln(I/A)/2)| (SPOT)", ident, n=len(spot), note="identity for I,A>0")
rl, pl = stats.pearsonr(g.SST, g.L)
rec("B", "SPOT ln(I/A) vs SST Pearson r", rl, n=len(g), p=pl)
ml = smf.ols("L ~ SST + C(frac)", data=g).fit(cov_type="HC3")
c = ml.conf_int().loc["SST"]
rec("B", "SPOT dln(I/A)/dSST, fraction-adjusted", ml.params.SST, n=len(g), ci_low=c[0], ci_high=c[1], p=ml.pvalues.SST)

# ---------------------------------------------------------------- C. heatwave (viral fraction)
v = spot[spot.frac == "viral"].sort_values("date").copy()
v["s"] = np.sin(2 * np.pi * v.month / 12)
v["c"] = np.cos(2 * np.pi * v.month / 12)
pre = v[v.date < "2014-01-01"]
hw = v[(v.date >= "2014-01-01") & (v.date <= "2015-12-31")]
harm = smf.ols("D_index ~ s + c", data=pre).fit()
v["D_exp"] = harm.predict(v)
v["resid"] = v.D_index - v.D_exp
pre = v[v.date < "2014-01-01"]
hw = v[(v.date >= "2014-01-01") & (v.date <= "2015-12-31")]
rec("C", "viral pre-event samples", len(pre), note="2009-07..2013-12")
rec("C", "viral 2014-2015 samples", len(hw))
rec("C", "pre-event harmonic R2", harm.rsquared, n=len(pre), p=harm.f_pvalue)
d_med = hw.resid.median() - pre.resid.median()
d_mean = hw.resid.mean() - pre.resid.mean()
rec("C", "season-adjusted median shift (this analysis)", d_med, n=len(pre) + len(hw))
rec("C", "season-adjusted mean shift", d_mean, n=len(pre) + len(hw))
rec("C", "SST anomaly mean pre / 2014-15", pre.SST_anomaly.mean(), note=f"heatwave {hw.SST_anomaly.mean():.2f}")
rec("C", "SST anomaly mean 2014-15", hw.SST_anomaly.mean(), n=len(hw))

# year-block bootstrap (resample calendar years within period) for the mean shift
yrs_pre = pre.year.unique(); yrs_hw = hw.year.unique()
bs = []
for _ in range(5000):
    bp = pd.concat([pre[pre.year == y] for y in RNG.choice(yrs_pre, len(yrs_pre))])
    bh = pd.concat([hw[hw.year == y] for y in RNG.choice(yrs_hw, len(yrs_hw))])
    bs.append(bh.resid.mean() - bp.resid.mean())
rec("C", "season-adjusted mean shift, year-block bootstrap CI", d_mean, ci_low=np.percentile(bs, 2.5),
    ci_high=np.percentile(bs, 97.5), note="5000 resamples of calendar years")
# exact year-label permutation over the 7 years 2009-2015
allyrs = sorted(set(yrs_pre) | set(yrs_hw))
vv = pd.concat([pre, hw])
obs = d_mean; perm = []
for combo in itertools.combinations(allyrs, len(yrs_hw)):
    m = vv.year.isin(combo)
    perm.append(vv[m].resid.mean() - vv[~m].resid.mean())
perm = np.array(perm)
p_perm = np.mean(perm >= obs - 1e-12)
rec("C", "year-label permutation P (one-sided, exact)", p_perm, n=len(perm),
    note=f"observed rank {int((perm >= obs - 1e-12).sum())} of {len(perm)}; minimum attainable P = {1/len(perm):.3f}")

# matched months: each heatwave sample vs pre-event mean of the same calendar month
clim = pre.groupby("month").D_index.mean()
diffs = (hw.D_index - hw.month.map(clim)).dropna()
w = stats.wilcoxon(diffs)
rec("C", "matched-month D difference, median", diffs.median(), n=len(diffs), p=w.pvalue, note="Wilcoxon signed-rank")
rec("C", "matched-month differences > 0", int((diffs > 0).sum()), n=len(diffs))

# anomaly slope learned before the event, used to predict the heatwave shift out of sample
ma = smf.ols("resid ~ SST_anomaly", data=pre).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
c = ma.conf_int().loc["SST_anomaly"]
rec("C", "pre-event dResid/dSST_anomaly", ma.params.SST_anomaly, n=len(pre), ci_low=c[0], ci_high=c[1], p=ma.pvalues.SST_anomaly)
pred_shift = ma.params.SST_anomaly * (hw.SST_anomaly.mean() - pre.SST_anomaly.mean())
se = ma.bse.SST_anomaly * (hw.SST_anomaly.mean() - pre.SST_anomaly.mean())
rec("C", "heatwave shift predicted from pre-event anomaly slope", pred_shift, ci_low=pred_shift - 1.96 * se,
    ci_high=pred_shift + 1.96 * se, note="compare with observed mean shift")
mm = smf.ols("resid ~ SST_anomaly + MEI", data=v[v.date <= "2015-12-31"]).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
for term in ["SST_anomaly", "MEI"]:
    c = mm.conf_int().loc[term]
    rec("C", f"viral 2009-2015 resid ~ SST_anomaly + MEI: {term}", mm.params[term], n=int(mm.nobs),
        ci_low=c[0], ci_high=c[1], p=mm.pvalues[term])
v.to_csv(OUT / "spot_viral_heatwave_residuals.csv", index=False)

# ---------------------------------------------------------------- D. Tara
t = tara[tara.D_index_qc == "ok"].copy()
t["femto"] = (t.size_fraction == "Femto-size").astype(int)
t["L"] = np.log((t.I_relative + EPS) / (t.A_relative + EPS))
rec("D", "Tara samples with valid D and predictors", len(t), note=f"{t.station_id.nunique()} stations")
rec("D", "Tara samples with D = +1 or -1 (one order absent)", int((t.D_index.abs() >= 0.9999).sum()), n=len(t))
st = t.groupby("station_id").agg(D=("D_index", "mean"), T=("temperature", "mean"), chl=("chlorophyll_a", "mean"),
                                 lat=("latitude", "mean"), lon=("longitude", "mean"), basin=("ocean_basin", "first")).reset_index()
for col, lab in [("T", "temperature"), ("chl", "chlorophyll a")]:
    x = st[col] if col == "T" else np.log10(st[col])
    r, p = stats.pearsonr(x, st.D)
    rho, prho = stats.spearmanr(x, st.D)
    rec("D", f"Tara station-mean D vs {lab} Pearson r", r, n=len(st), p=p)
    rec("D", f"Tara station-mean D vs {lab} Spearman rho", rho, n=len(st), p=prho)
rec("D", "Tara station-mean |latitude| vs D Spearman rho", *(lambda s: (s[0], len(st), None, None, s[1]))(stats.spearmanr(st.lat.abs(), st.D)))

mt = smf.ols("L ~ temperature + femto", data=t).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(t.station_id)[0]})
c = mt.conf_int().loc["temperature"]
rec("D", "Tara dln(I/A)/dT (station-clustered)", mt.params.temperature, n=len(t), ci_low=c[0], ci_high=c[1], p=mt.pvalues.temperature)
mtd = smf.ols("D_index ~ temperature + femto", data=t).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(t.station_id)[0]})
c = mtd.conf_int().loc["temperature"]
rec("D", "Tara dD/dT (station-clustered)", mtd.params.temperature, n=len(t), ci_low=c[0], ci_high=c[1], p=mtd.pvalues.temperature)
mid = t[(t.temperature >= 13) & (t.temperature <= 24)]
mmid = smf.ols("D_index ~ temperature + femto", data=mid).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(mid.station_id)[0]})
c = mmid.conf_int().loc["temperature"]
rec("D", "Tara dD/dT within SPOT SST range 13-24 degC", mmid.params.temperature, n=len(mid), ci_low=c[0], ci_high=c[1], p=mmid.pvalues.temperature)

# temperature at D = 0 (pico and femto) from the log-ratio fit, station bootstrap
def t0(df):
    fit = smf.ols("L ~ temperature + femto", data=df).fit()
    b = fit.params
    return -b.Intercept / b.temperature, -(b.Intercept + b.femto) / b.temperature
t0p, t0f = t0(t)
stations = t.station_id.unique(); bp = []
for _ in range(2000):
    pick = RNG.choice(stations, len(stations))
    bp.append(t0(pd.concat([t[t.station_id == s] for s in pick])))
bp = np.array(bp)
rec("D", "Tara temperature at D=0, pico", t0p, n=len(t), ci_low=np.percentile(bp[:, 0], 2.5), ci_high=np.percentile(bp[:, 0], 97.5), note="station bootstrap")
rec("D", "Tara temperature at D=0, femto", t0f, n=len(t), ci_low=np.percentile(bp[:, 1], 2.5), ci_high=np.percentile(bp[:, 1], 97.5), note="station bootstrap")

# station-blocked and basin-blocked cross-validation
t["log_chl"] = np.log10(t.chlorophyll_a)
for col in ["CO3", "PO4", "NO3", "Si", "Fe"]:
    t["log_" + col] = np.log10(t[col].clip(lower=t[col][t[col] > 0].min() / 2))
FEAT = ["temperature", "salinity", "log_chl", "log_CO3", "log_PO4", "log_NO3", "log_Si", "log_Fe", "femto"]
X = t[FEAT].to_numpy(); yD = t.D_index.to_numpy(); yL = t.L.to_numpy()
groups = t.station_id.to_numpy(); basins = t.ocean_basin.to_numpy()


def models():
    return {
        "temperature-only D": ("D", ["temperature", "femto"], lambda: LinearRegression()),
        "temperature-only log-ratio": ("L", ["temperature", "femto"], lambda: LinearRegression()),
        "ridge log-ratio (8 env + fraction)": ("L", FEAT, lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25)))),
        "random forest log-ratio (8 env + fraction)": ("L", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
        "random forest D (8 env + fraction)": ("D", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
    }


def to_D(pred, target):
    return np.tanh(pred / 2) if target == "L" else np.clip(pred, -1, 1)


def cv_run(split_iter, tag, perm_imp=False):
    out = {}
    imp = []
    for name, (target, feats, mk) in models().items():
        pred = np.full(len(t), np.nan)
        for tr, te in split_iter():
            Xtr = t.iloc[tr][feats].to_numpy(); Xte = t.iloc[te][feats].to_numpy()
            y = (yL if target == "L" else yD)[tr]
            mdl = mk().fit(Xtr, y)
            pred[te] = to_D(mdl.predict(Xte), target)
            if perm_imp and name == "random forest log-ratio (8 env + fraction)":
                base = r2_score(yD[te], pred[te])
                for j, fname in enumerate(feats):
                    for _ in range(10):
                        Xp = Xte.copy(); Xp[:, j] = RNG.permutation(Xp[:, j])
                        imp.append(dict(feature=fname, dR2=base - r2_score(yD[te], to_D(mdl.predict(Xp), target))))
        out[name] = pred
    return out, imp


def station_folds(seed, k=5):
    def it():
        rng = np.random.default_rng(seed)
        us = rng.permutation(np.unique(groups))
        fold_of = {s: i % k for i, s in enumerate(us)}
        f = np.array([fold_of[s] for s in groups])
        for i in range(k):
            yield np.where(f != i)[0], np.where(f == i)[0]
    return it


cv_rows = []; all_imp = []; oof_keep = None
for rep in range(20):
    res, imp = cv_run(station_folds(rep), "station", perm_imp=True)
    all_imp += imp
    for name, pred in res.items():
        cv_rows.append(dict(scheme="5-fold station-blocked", rep=rep, model=name, R2=r2_score(yD, pred),
                            RMSE=mean_squared_error(yD, pred) ** 0.5, MAE=mean_absolute_error(yD, pred)))
    if rep == 0:
        oof_keep = pd.DataFrame({k: v for k, v in res.items()})


def lobo():
    for b in np.unique(basins):
        yield np.where(basins != b)[0], np.where(basins == b)[0]


res_b, _ = cv_run(lambda: lobo(), "basin")
for name, pred in res_b.items():
    cv_rows.append(dict(scheme="leave-one-basin-out", rep=0, model=name, R2=r2_score(yD, pred),
                        RMSE=mean_squared_error(yD, pred) ** 0.5, MAE=mean_absolute_error(yD, pred)))
    for b in np.unique(basins):
        m = basins == b
        cv_rows.append(dict(scheme=f"LOBO held-out basin: {b}", rep=0, model=name, n=int(m.sum()),
                            R2=r2_score(yD[m], pred[m]) if m.sum() > 2 else np.nan,
                            RMSE=mean_squared_error(yD[m], pred[m]) ** 0.5, MAE=mean_absolute_error(yD[m], pred[m]),
                            bias=float(np.mean(pred[m] - yD[m]))))
cv = pd.DataFrame(cv_rows)
cv.to_csv(OUT / "tara_cross_validation.csv", index=False)
summ = cv[cv.scheme == "5-fold station-blocked"].groupby("model")[["R2", "RMSE", "MAE"]].agg(["mean", "std"])
summ.to_csv(OUT / "tara_cv_summary_station_blocked.csv")
for name in models():
    s = cv[(cv.scheme == "5-fold station-blocked") & (cv.model == name)]
    rec("D", f"station-blocked CV R2: {name}", s.R2.mean(), n=len(t), ci_low=s.R2.min(), ci_high=s.R2.max(), note="20 random station partitions; range")
    rec("D", f"station-blocked CV RMSE: {name}", s.RMSE.mean(), n=len(t))
    b = cv[(cv.scheme == "leave-one-basin-out") & (cv.model == name)].iloc[0]
    rec("D", f"leave-one-basin-out R2: {name}", b.R2, n=len(t))
    rec("D", f"leave-one-basin-out RMSE: {name}", b.RMSE, n=len(t))
rec("D", "supplied pipeline OOF R2 (MD5 env-hash groups)", r2_score(val.observed_D, val.predicted_D), n=len(val))
rec("D", "supplied pipeline OOF RMSE", mean_squared_error(val.observed_D, val.predicted_D) ** 0.5, n=len(val))
# do the MD5 groups coincide with stations?
gmap = val.groupby("group").station_id.nunique().max(), val.groupby("station_id").group.nunique().max()
rec("D", "supplied CV: max stations per group / groups per station", float(gmap[0]), note=f"groups per station max {gmap[1]}")

impdf = pd.DataFrame(all_imp).groupby("feature").dR2.agg(["mean", "std", "count"]).sort_values("mean", ascending=False)
impdf["ci_low"] = pd.DataFrame(all_imp).groupby("feature").dR2.quantile(.025)
impdf["ci_high"] = pd.DataFrame(all_imp).groupby("feature").dR2.quantile(.975)
impdf.to_csv(OUT / "tara_oof_permutation_importance.csv")
for f_, r_ in impdf.iterrows():
    rec("D", f"OOF permutation dR2: {f_}", r_["mean"], n=int(r_["count"]), ci_low=r_.ci_low, ci_high=r_.ci_high,
        note="RF log-ratio; permuted within held-out station folds; 20 partitions x 10 permutations")
oof_keep.insert(0, "observed_D", yD); oof_keep.insert(0, "station_id", groups); oof_keep.insert(1, "basin", basins)
oof_keep["temperature"] = t.temperature.values; oof_keep["size_fraction"] = t.size_fraction.values
oof_keep.to_csv(OUT / "tara_oof_predictions_station_blocked_rep0.csv", index=False)
st.to_csv(OUT / "tara_station_means.csv", index=False)

# ---------------------------------------------------------------- E. cross-system transfer
cal = smf.ols("L ~ temperature + femto", data=t).fit()
sp = spot[spot.frac.isin(["viral", "pico"])].dropna(subset=["SST"]).copy()
sp["femto"] = (sp.frac == "viral").astype(int)
sp["D_tara"] = np.tanh(cal.predict(sp.rename(columns={"SST": "temperature"})) / 2)
r, p = stats.pearsonr(sp.D_tara, sp.D_index)
rec("E", "SPOT observed D vs D predicted from Tara temperature calibration, r", r, n=len(sp), p=p, note="viral+pico; Tara femto/pico calibration")
rec("E", "transfer bias (predicted - observed D), mean", float((sp.D_tara - sp.D_index).mean()), n=len(sp))
rec("E", "transfer RMSE (D units)", float(mean_squared_error(sp.D_index, sp.D_tara) ** 0.5), n=len(sp))
# inverse: temperature read from D with the Tara calibration
Lsp = np.log((sp.I + EPS) / (sp.A + EPS))
T_hat = (Lsp - cal.params.Intercept - cal.params.femto * sp.femto) / cal.params.temperature
r, p = stats.pearsonr(T_hat, sp.SST)
rec("E", "SST read from D (Tara calibration) vs observed SST, r", r, n=len(sp), p=p)
rec("E", "SST read from D, median absolute error degC", float(np.median(np.abs(T_hat - sp.SST))), n=len(sp))
for f, gg in sp.groupby("frac"):
    r, p = stats.pearsonr(gg.D_tara, gg.D_index)
    rec("E", f"transfer r ({f})", r, n=len(gg), p=p)
sp[["sample_id", "date", "frac", "SST", "D_index", "D_tara"]].assign(T_from_D=T_hat).to_csv(OUT / "spot_transfer_from_tara.csv", index=False)

# ---------------------------------------------------------------- F. CMIP6 front
fr = front[front.ESM != "ENSEMBLE-mean"].copy()
fr["abslat"] = fr.D0_front_latitude.abs()
toe_rows = []
for hem, gh in fr.groupby("hemisphere"):
    lr = stats.linregress(gh.year, gh.abslat)
    rec("F", f"pooled front rate {hem} km/decade", lr.slope * KM_PER_DEG * 10, n=len(gh),
        ci_low=(lr.slope - 2.024 * lr.stderr) * KM_PER_DEG * 10, ci_high=(lr.slope + 2.024 * lr.stderr) * KM_PER_DEG * 10,
        p=lr.pvalue, note=f"R2={lr.rvalue**2:.3f}; 4 ESM x 10 decades")
    # ESM fixed effects: removes between-model offset in front position
    fe = smf.ols("abslat ~ year + C(ESM)", data=gh).fit()
    c = fe.conf_int().loc["year"]
    rec("F", f"front rate with ESM fixed effects {hem} km/decade", fe.params.year * KM_PER_DEG * 10, n=len(gh),
        ci_low=c[0] * KM_PER_DEG * 10, ci_high=c[1] * KM_PER_DEG * 10, p=fe.pvalues.year, note=f"R2={fe.rsquared:.3f}")
    for esm, ge in gh.groupby("ESM"):
        l2 = stats.linregress(ge.year, ge.abslat)
        resid_sd = np.std(ge.abslat - (l2.intercept + l2.slope * ge.year), ddof=2)
        step = l2.slope * 10  # degrees per decade
        # decades until trend displacement exceeds 2 x decadal residual noise
        toe = 2 * resid_sd / step if step > 0 else np.inf
        toe_rows.append(dict(hemisphere=hem, ESM=esm, km_per_decade=l2.slope * KM_PER_DEG * 10, R2=l2.rvalue ** 2,
                             P=l2.pvalue, resid_sd_km=resid_sd * KM_PER_DEG, snr_per_decade=step / resid_sd,
                             decades_to_2sigma=toe, lat_2000s=ge.abslat.iloc[0], lat_2090s=ge.abslat.iloc[-1]))
toe = pd.DataFrame(toe_rows)
toe.to_csv(OUT / "cmip6_front_per_esm.csv", index=False)
for hem, gt in toe.groupby("hemisphere"):
    rec("F", f"four-model mean front rate {hem} km/decade", gt.km_per_decade.mean(), n=len(gt),
        ci_low=gt.km_per_decade.min(), ci_high=gt.km_per_decade.max(),
        note="mean across four ESM slopes; ci_low/ci_high store model range, not a confidence interval")
    rec("F", f"per-ESM rate range {hem} km/decade", gt.km_per_decade.median(), ci_low=gt.km_per_decade.min(), ci_high=gt.km_per_decade.max(), note="median, range")
    rec("F", f"per-ESM R2 range {hem}", gt.R2.median(), ci_low=gt.R2.min(), ci_high=gt.R2.max())
    rec("F", f"per-ESM decadal noise around trend {hem} km", gt.resid_sd_km.median(), ci_low=gt.resid_sd_km.min(), ci_high=gt.resid_sd_km.max())
    rec("F", f"per-ESM signal-to-noise per decade {hem}", gt.snr_per_decade.median(), ci_low=gt.snr_per_decade.min(), ci_high=gt.snr_per_decade.max())
    rec("F", f"per-ESM decades to 2-sigma displacement {hem}", gt.decades_to_2sigma.median(), ci_low=gt.decades_to_2sigma.min(), ci_high=gt.decades_to_2sigma.max())
    rec("F", f"ensemble front displacement 2000s->2090s {hem} km", (gt.lat_2090s - gt.lat_2000s).mean() * KM_PER_DEG,
        ci_low=(gt.lat_2090s - gt.lat_2000s).min() * KM_PER_DEG, ci_high=(gt.lat_2090s - gt.lat_2000s).max() * KM_PER_DEG)
    rec("F", f"all ESMs move poleward {hem}", int((gt.km_per_decade > 0).sum()), n=len(gt))

y0, y1 = latyr.year.min(), latyr.year.max()
b = latyr[latyr.year == y0].groupby(["ESM", "latitude_mid"]).mean_D.mean()
e = latyr[latyr.year == y1].groupby(["ESM", "latitude_mid"]).mean_D.mean()
dd = (e - b).reset_index().rename(columns={"mean_D": "dD"})
agree = dd.groupby("latitude_mid").dD.agg(lambda x: max((x > 0).sum(), (x < 0).sum()))
rec("F", "latitude bins where all 4 ESMs agree on sign of dD", float((agree == 4).mean()), n=len(agree), note="fraction of 1-degree bins")
band = dd.assign(band=pd.cut(dd.latitude_mid, [-80, -60, -30, 30, 60, 80])).groupby(["band", "ESM"], observed=True).dD.mean().groupby("band", observed=True).agg(["mean", "min", "max"])
band.to_csv(OUT / "cmip6_dD_by_band.csv")
for bnd, r_ in band.iterrows():
    rec("F", f"dD 2090s-2000s band {bnd}", r_["mean"], ci_low=r_["min"], ci_high=r_["max"], note="4-ESM mean, ESM range")

hs = pd.DataFrame(rows)
hs.to_csv(OUT / "headline_statistics.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 70, "display.max_rows", 400)
print(hs.drop(columns=["section"]).to_string(float_format=lambda x: f"{x:.4g}"))
