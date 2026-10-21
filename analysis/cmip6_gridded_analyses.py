"""Gridded CMIP6 analyses (run after cmip6_extract.py).

1  Reproduce the supplied latitude profiles and D = 0 fronts from the gridded predictions.
2  2-degree maps: four-model mean D (2000s, 2090s), change, sign agreement; agreement by latitude band.
3  Area of applicability (Meyer & Pebesma 2021): dissimilarity of each grid cell to the Tara calibration data.
4  Sensitivity projections refitted on Tara and applied to the gridded predictors:
     a  temperature-only log-ratio model   ln(I/A) ~ T + size class
     b  random forest, 8 predictors        (reimplementation of the supplied design)
     c  random forest, 7 predictors        (carbonate removed: CMIP6 and Tara CO3 differ in units)
   For each: D per size class, averaged over size classes, zonal means, fronts and poleward rates per ESM.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor

HERE = Path(__file__).resolve().parent
DATA_ROOT = HERE.parent  # repo root: raw_tables/ + processed_tables/
RNG = np.random.default_rng(29)
KM = 111.32
ESMS = ["CESM2-WACCM", "IPSL-CM6A-LR", "NorESM2-MM", "UKESM1-0-LL"]
DECS = [f"{y}-{y + 9}" for y in range(2000, 2100, 10)]
YEAR = {d: int(d[:4]) + 4.5 for d in DECS}
EPS = 1e-3
z = np.load(HERE / "cmip6_grid.npz")
C = {k: i for i, k in enumerate(z["columns"])}
rows = []


def rec(name, value, n=None, lo=None, hi=None, p=None, note=""):
    rows.append(dict(section="CMIP6-gridded", statistic=name, value=value, n=n, ci_low=lo, ci_high=hi, P=p, note=note))


# ------------------------------------------------------------------ helpers
def zonal(lat, D, width=1.0):
    b = np.floor(lat / width) * width + width / 2
    s = pd.DataFrame({"lat": b, "D": D}).groupby("lat").D.agg(["mean", "size"])
    return s


def fronts(profile, lo=15, hi=65):
    """First positive-to-negative crossing moving poleward from `lo` degrees, per hemisphere, linear interpolation."""
    out = {}
    for hem, sign in [("North", 1), ("South", -1)]:
        p = profile.copy()
        p = p[(sign * p.index >= lo) & (sign * p.index <= hi)].sort_index(ascending=(sign < 0) is False)
        la, d = p.index.to_numpy() * sign, p["mean"].to_numpy()
        order = np.argsort(la); la, d = la[order], d[order]
        f = np.nan
        for i in range(len(la) - 1):
            if d[i] > 0 >= d[i + 1]:
                f = la[i] + (la[i + 1] - la[i]) * d[i] / (d[i] - d[i + 1])
                break
        out[hem] = f
    return out


def rate_table(front_df, label):
    res = []
    for (esm, hem), g in front_df.groupby(["ESM", "hemisphere"]):
        g = g.dropna(subset=["front"])
        lr = stats.linregress(g.year, g.front)
        res.append(dict(projection=label, ESM=esm, hemisphere=hem, n=len(g), km_per_decade=lr.slope * KM * 10,
                        R2=lr.rvalue ** 2, P=lr.pvalue, front_2000s=g.front.iloc[0], front_2090s=g.front.iloc[-1]))
    return pd.DataFrame(res)


# ------------------------------------------------------------------ Tara calibration data (same transforms as DATA_ROOT)
tara = pd.read_csv(DATA_ROOT / "data/Tara_observed_Dindex.csv")
tara = tara[tara.D_index_qc == "ok"].copy()
tara["femto"] = (tara.size_fraction == "Femto-size").astype(int)
tara["L"] = np.log((tara.I_relative + EPS) / (tara.A_relative + EPS))
ENV_T = ["temperature", "salinity", "chlorophyll_a", "CO3", "PO4", "NO3", "Si", "Fe"]
ENV_M = ["Temperature", "Salinity", "ChlorophyllA", "CO3", "PO4", "NO3", "Si", "Iron"]
LOGV = {"chlorophyll_a", "CO3", "PO4", "NO3", "Si", "Fe"}
floors = {v: tara.loc[tara[v] > 0, v].min() / 2 for v in LOGV}


def transform(df, names_src):
    X = np.empty((len(df), 8), np.float32)
    for j, (vt, vs) in enumerate(zip(ENV_T, names_src)):
        x = np.asarray(df[vs], np.float64)
        X[:, j] = np.log10(np.clip(x, floors[vt], None)) if vt in LOGV else x
    return X


Xt = transform(tara, ENV_T)
femto = tara.femto.to_numpy()

# ------------------------------------------------------------------ 1 + 2: supplied predictions on grids
supplied = pd.read_csv(DATA_ROOT / "data/CMIP6_Dindex_latitude_year.csv")
sup_front = pd.read_csv(DATA_ROOT / "data/Dindex_front_trajectory.csv")
prof_rows, front_rows, grid = [], [], {}
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        lat, lon, D = a[:, C["latitude"]], a[:, C["longitude"]], a[:, C["D_index_pred"]]
        zp = zonal(lat, D)
        for la, r in zp.iterrows():
            prof_rows.append(dict(ESM=esm, decade=dec, year=YEAR[dec], lat=la, D=r["mean"], n=r["size"]))
        fr = fronts(zp)
        for hem, f in fr.items():
            front_rows.append(dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=hem, front=f))
        gi = ((lat + 90) // 2).astype(int).clip(0, 89); gj = ((lon + 180) // 2).astype(int).clip(0, 179)
        s = np.zeros((90, 180)); n = np.zeros((90, 180))
        np.add.at(s, (gi, gj), D); np.add.at(n, (gi, gj), 1)
        with np.errstate(invalid="ignore"):
            grid[(esm, dec)] = np.where(n > 0, s / n, np.nan)
prof = pd.DataFrame(prof_rows); fr_g = pd.DataFrame(front_rows)
prof.to_csv(HERE / "cmip6_gridded_zonal_profiles.csv", index=False)

# provenance: compare with supplied latitude table and front trajectories
m = prof.merge(supplied.rename(columns={"period": "decade", "latitude_mid": "lat"}), on=["ESM", "decade", "lat"], how="inner")
rec("gridded zonal-mean D vs supplied latitude table: max |difference|", float(np.abs(m.D - m.mean_D).max()), n=len(m))
cmp = fr_g.merge(sup_front.rename(columns={"D0_front_latitude": "sup"}), on=["ESM", "year", "hemisphere"])
cmp["sup_abs"] = cmp.sup.abs()
rec("gridded front vs supplied front: max |difference| (deg)", float(np.abs(cmp.front - cmp.sup_abs).max()), n=len(cmp))
rt_sup = rate_table(fr_g, "supplied ensemble (gridded)")

# maps and agreement
G = {dec: np.stack([grid[(e, dec)] for e in ESMS]) for dec in ["2000-2009", "2090-2099"]}
dG = G["2090-2099"] - G["2000-2009"]
common = np.all(np.isfinite(dG), axis=0)
agree = common & (np.all(dG > 0, axis=0) | np.all(dG < 0, axis=0))
latc = -89 + 2 * np.arange(90)
coslat = np.cos(np.radians(latc))[:, None] * np.ones((1, 180))
rec("2-degree ocean cells common to all four ESMs", int(common.sum()))
rec("area-weighted share of cells with four-model sign agreement on dD", float((agree * coslat).sum() / (common * coslat).sum()))
rec("area-weighted share of cells where four-model mean dD > 0", float(((np.nanmean(dG, 0) > 0) & common) .__mul__(coslat).sum() / (common * coslat).sum()))
band_rows = []
for lo, hi, lab in [(0, 30, "|lat| < 30"), (30, 45, "30-45"), (45, 60, "45-60"), (60, 90, "> 60")]:
    for hem, sgn in [("both", None), ("North", 1), ("South", -1)]:
        sel = (np.abs(latc) >= lo) & (np.abs(latc) < hi)
        if sgn is not None:
            sel &= (np.sign(latc) == sgn)
        mask = common & sel[:, None]
        if mask.sum() == 0:
            continue
        share = float((agree & mask).__mul__(coslat).sum() / (mask * coslat).sum())
        mean_dD = float(np.nansum(np.nanmean(dG, 0) * coslat * mask) / (mask * coslat).sum())
        band_rows.append(dict(band=lab, hemisphere=hem, n_cells=int(mask.sum()), share_sign_agreement=share, mean_dD=mean_dD))
        if hem == "both":
            rec(f"sign agreement, {lab}", share, n=int(mask.sum()), note=f"area-weighted; mean dD {mean_dD:.3f}")
pd.DataFrame(band_rows).to_csv(HERE / "cmip6_agreement_by_band.csv", index=False)
np.savez_compressed(HERE / "cmip6_maps_2deg.npz", D2000=G["2000-2009"], D2090=G["2090-2099"], esms=np.array(ESMS))

# ------------------------------------------------------------------ 3: area of applicability
imp = pd.read_csv(DATA_ROOT / "analysis/tara_oof_permutation_importance.csv").set_index("feature")["mean"].clip(lower=0)
w_map = {"temperature": "temperature", "salinity": "salinity", "log_chl": "chlorophyll_a", "log_CO3": "CO3",
         "log_PO4": "PO4", "log_NO3": "NO3", "log_Si": "Si", "log_Fe": "Fe"}
W = np.array([imp.get(k, 0) for k in w_map], float)  # order follows ENV_T
mu, sd = Xt.mean(0), Xt.std(0)


def di_setup(weights, use):
    Z = (Xt[:, use] - mu[use]) / sd[use] * weights[use]
    dmean = np.mean(stats.distributions.norm.pdf(0) * 0 + np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))[np.triu_indices(len(Z), 1)])
    # station-blocked cross-validated training DI (five folds of stations); the threshold is the median over 50 random
    # station partitions (seeded) because a single partition makes it unstable
    st = tara.station_id.to_numpy(); thrs = []
    for rep in range(50):
        us = np.random.default_rng(1000 + rep).permutation(np.unique(st)); fold = {s: i % 5 for i, s in enumerate(us)}
        f = np.array([fold[s] for s in st]); tdi = np.empty(len(Z))
        for k in range(5):
            te, tr = f == k, f != k
            tdi[te] = np.sqrt(((Z[te][:, None, :] - Z[tr][None, :, :]) ** 2).sum(-1)).min(1) / dmean
        q1, q3 = np.percentile(tdi, [25, 75])
        thrs.append(tdi[tdi <= q3 + 1.5 * (q3 - q1)].max())
    global AOA_THRESHOLD_RANGE
    AOA_THRESHOLD_RANGE = (float(np.min(thrs)), float(np.median(thrs)), float(np.max(thrs)))
    return Z, dmean, float(np.median(thrs))


def di(X, weights, use, Z, dmean):
    Q = (X[:, use] - mu[use]) / sd[use] * weights[use]
    out = np.empty(len(Q), np.float32)
    for s in range(0, len(Q), 4000):
        out[s:s + 4000] = np.sqrt(((Q[s:s + 4000, None, :] - Z[None, :, :]) ** 2).sum(-1)).min(1) / dmean
    return out


AOA_SPECS = {"importance-weighted, 8 predictors": (W, np.arange(8)),
             "unweighted, 8 predictors": (np.ones(8), np.arange(8)),
             "unweighted, 7 predictors (no CO3)": (np.ones(8), np.array([0, 1, 2, 4, 5, 6, 7]))}
aoa_rows = []
aoa_maps = {}
for spec, (wts, use) in AOA_SPECS.items():
    Z, dmean, thr = di_setup(wts, use)
    for esm in ESMS:
        for dec in ["2000-2009", "2090-2099"]:
            a = z[f"{esm}__{dec}"]
            Xm = transform(pd.DataFrame(a[:, [C[v] for v in ENV_M]], columns=ENV_M), ENV_M)
            d = di(Xm, wts, use, Z, dmean)
            inside = d <= thr
            D = a[:, C["D_index_pred"]]; lat = a[:, C["latitude"]]
            wlat = np.cos(np.radians(lat))
            near = np.abs(D) < 0.15
            polar = np.abs(lat) >= 60
            aoa_rows.append(dict(spec=spec, ESM=esm, decade=dec, threshold=thr,
                                 share_all=float((inside * wlat).sum() / wlat.sum()),
                                 share_near_front=float((inside * wlat)[near].sum() / wlat[near].sum()),
                                 share_polar=float((inside * wlat)[polar].sum() / wlat[polar].sum())))
            if spec == "unweighted, 7 predictors (no CO3)":
                gi = ((lat + 90) // 2).astype(int).clip(0, 89); gj = ((a[:, C["longitude"]] + 180) // 2).astype(int).clip(0, 179)
                s_ = np.zeros((90, 180)); n_ = np.zeros((90, 180)); np.add.at(s_, (gi, gj), inside); np.add.at(n_, (gi, gj), 1)
                with np.errstate(invalid="ignore"):
                    aoa_maps[(esm, dec)] = np.where(n_ > 0, s_ / n_, np.nan)
aoa = pd.DataFrame(aoa_rows)
aoa.to_csv(HERE / "cmip6_area_of_applicability.csv", index=False)
for spec in AOA_SPECS:
    for dec in ["2000-2009", "2090-2099"]:
        s = aoa[(aoa.spec == spec) & (aoa.decade == dec)]
        for col, lab in [("share_all", "all ocean"), ("share_near_front", "near fronts |D|<0.15"), ("share_polar", "|lat| >= 60")]:
            rec(f"AOA {spec}, {dec}, {lab}", s[col].mean(), n=4, lo=s[col].min(), hi=s[col].max(), note="mean and range over ESMs, area-weighted")
np.savez_compressed(HERE / "cmip6_aoa_2deg.npz",
                    aoa2000=np.nanmean(np.stack([aoa_maps[(e, "2000-2009")] for e in ESMS]), 0),
                    aoa2090=np.nanmean(np.stack([aoa_maps[(e, "2090-2099")] for e in ESMS]), 0))

# ------------------------------------------------------------------ 4: sensitivity projections
lin = np.linalg.lstsq(np.c_[np.ones(len(tara)), tara.temperature, femto], tara.L.to_numpy(), rcond=None)[0]
rec("temperature-only calibration: ln(I/A) = b0 + bT*T + bF*femto", float(lin[1]), note=f"b0={lin[0]:.3f}, bF={lin[2]:.3f}")
rec("temperature at D=0 (pico, femto) in temperature-only calibration", float(-lin[0] / lin[1]), note=f"femto {-(lin[0] + lin[2]) / lin[1]:.2f} degC")


def rf_fit(use):
    ms = {}
    for tgt in ["I_relative", "A_relative"]:
        Xf = np.c_[Xt[:, use], femto]
        ms[tgt] = RandomForestRegressor(300, min_samples_leaf=3, random_state=1, n_jobs=2).fit(Xf, tara[tgt].to_numpy())
    return ms


RF = {"random forest, 8 predictors": (rf_fit(np.arange(8)), np.arange(8)),
      "random forest, 7 predictors (no CO3)": (rf_fit(np.array([0, 1, 2, 4, 5, 6, 7])), np.array([0, 1, 2, 4, 5, 6, 7]))}
sens_fronts = {k: [] for k in ["supplied ensemble, 2-degree cells", "temperature-only", *RF]}
agg_rows = []
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        lat, lon = a[:, C["latitude"]], a[:, C["longitude"]]
        # aggregate predictors to 2 x 2 degree cells to keep random-forest prediction tractable
        key = pd.DataFrame({"la": np.floor(lat / 2) * 2 + 1, "lo": np.floor(lon / 2) * 2 + 1})
        env = pd.DataFrame(a[:, [C[v] for v in ENV_M]], columns=ENV_M).assign(Dsup=a[:, C["D_index_pred"]])
        cell = pd.concat([key, env], axis=1).groupby(["la", "lo"]).mean().reset_index()
        Xm = transform(cell, ENV_M)
        T = cell.Temperature.to_numpy()
        Dt = 0.5 * (np.tanh((lin[0] + lin[1] * T) / 2) + np.tanh((lin[0] + lin[1] * T + lin[2]) / 2))
        for hem, f in fronts(zonal(cell.la.to_numpy(), cell.Dsup.to_numpy(), 2.0)).items():
            sens_fronts["supplied ensemble, 2-degree cells"].append(dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=hem, front=f))
        for hem, f in fronts(zonal(cell.la.to_numpy(), Dt, 2.0)).items():
            sens_fronts["temperature-only"].append(dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=hem, front=f))
        for name, (ms, use) in RF.items():
            Ds = []
            for fr in (0, 1):
                Xf = np.c_[Xm[:, use], np.full(len(Xm), fr)]
                I_, A_ = ms["I_relative"].predict(Xf), ms["A_relative"].predict(Xf)
                Ds.append((I_ - A_) / np.clip(I_ + A_, 1e-6, None))
            Dr = 0.5 * (Ds[0] + Ds[1])
            for hem, f in fronts(zonal(cell.la.to_numpy(), Dr, 2.0)).items():
                sens_fronts[name].append(dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=hem, front=f))
        print(esm, dec, flush=True)

rates = [rt_sup] + [rate_table(pd.DataFrame(v), k) for k, v in sens_fronts.items()]
rates = pd.concat(rates)
rates.to_csv(HERE / "cmip6_front_rates_by_projection.csv", index=False)
pd.concat([pd.DataFrame(v).assign(projection=k) for k, v in sens_fronts.items()] + [fr_g.assign(projection="supplied ensemble (gridded)")]).to_csv(
    HERE / "cmip6_front_trajectories_by_projection.csv", index=False)
for (proj, hem), g in rates.groupby(["projection", "hemisphere"]):
    rec(f"front rate {hem}, {proj} (km/decade)", g.km_per_decade.mean(), n=len(g), lo=g.km_per_decade.min(), hi=g.km_per_decade.max(),
        note=f"mean and ESM range; ESMs poleward {int((g.km_per_decade > 0).sum())}/4; min R2 {g.R2.min():.2f}")

out = pd.DataFrame(rows)
out.to_csv(HERE / "cmip6_gridded_statistics.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 90, "display.max_rows", 300)
print(out.drop(columns=["section"]).to_string(float_format=lambda x: f"{x:.4g}"))
print(rates.round(3).to_string())
