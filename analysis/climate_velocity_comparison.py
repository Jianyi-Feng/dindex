"""Physical comparators for the D = 0 dominance-front rates, from the same CMIP6 fields.

Two quantities are computed from analysis/cmip6_grid.npz, so that the models, grid, decades,
scenario and front-tracking algorithm are identical to the ecological projection:

  1. Isotherm displacement: the zonal-mean SST isotherm that coincides with the D = 0 front in the
     2000s is tracked by decade with the same crossing algorithm used for the fronts (km per decade).
  2. Classic climate velocity (Burrows et al. 2011): VoCC = (dSST/dt) / |grad SST|, on a 2 deg grid,
     summarized globally and in the latitude bands that contain the fronts.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
KM = 111.32
ESMS = ["CESM2-WACCM", "IPSL-CM6A-LR", "NorESM2-MM", "UKESM1-0-LL"]
DECS = [f"{y}-{y + 9}" for y in range(2000, 2100, 10)]
YEARS = [int(d[:4]) + 4.5 for d in DECS]

z = np.load(HERE / "cmip6_grid.npz", allow_pickle=True)
cols = {c: i for i, c in enumerate(map(str, z["columns"]))}


def zonal(lat, v, width=1.0):
    b = np.floor(lat / width) * width + width / 2
    return pd.DataFrame({"lat": b, "v": v}).groupby("lat").v.mean()


def crossing(profile, lo=15, hi=65):
    """First positive-to-negative crossing moving poleward from `lo`, per hemisphere (as in fronts())."""
    out = {}
    for hem, sign in [("North", 1), ("South", -1)]:
        p = profile[(sign * profile.index >= lo) & (sign * profile.index <= hi)]
        la, d = p.index.to_numpy() * sign, p.to_numpy()
        order = np.argsort(la); la, d = la[order], d[order]
        f = np.nan
        for i in range(len(la) - 1):
            if d[i] > 0 >= d[i + 1]:
                f = la[i] + (la[i + 1] - la[i]) * d[i] / (d[i] - d[i + 1])
                break
        out[hem] = f
    return out


# ---------------------------------------------------------------- 1. isotherm tracking
dfront = pd.read_csv(HERE / "cmip6_front_trajectories_by_projection.csv")
dfront = dfront[dfront.projection == "final ensemble (7 predictors)"]

sst_prof, rows = {}, []
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        sst_prof[(esm, dec)] = zonal(a[:, cols["latitude"]], a[:, cols["Temperature"]])

for esm in ESMS:
    p0 = sst_prof[(esm, DECS[0])]
    for hem in ("North", "South"):
        f0 = dfront[(dfront.ESM == esm) & (dfront.hemisphere == hem) & (dfront.year == YEARS[0])].front
        if not len(f0) or not np.isfinite(f0.iloc[0]):
            continue
        lat0 = f0.iloc[0] * (1 if hem == "North" else -1)
        t_ref = float(np.interp(lat0, p0.index.to_numpy(), p0.to_numpy()))   # SST at the 2000s front
        lats = []
        for dec in DECS:
            lats.append(crossing(sst_prof[(esm, dec)] - t_ref)[hem])
        g = pd.DataFrame({"year": YEARS, "lat": lats}).dropna()
        lr = stats.linregress(g.year, g.lat)
        rows.append(dict(ESM=esm, hemisphere=hem, isotherm_degC=round(t_ref, 2),
                         km_per_decade=lr.slope * KM * 10, R2=lr.rvalue ** 2,
                         lat_2000s=g.lat.iloc[0], lat_2090s=g.lat.iloc[-1]))
iso = pd.DataFrame(rows)

# ---------------------------------------------------------------- 2. climate velocity on a 2 deg grid
lon2 = -179 + 2 * np.arange(180); lat2 = -89 + 2 * np.arange(90)


def grid2(lon, lat, val):
    xi = np.clip(((lon + 180) / 2).astype(int), 0, 179)
    yi = np.clip(((lat + 90) / 2).astype(int), 0, 89)
    s = np.zeros((90, 180)); n = np.zeros((90, 180))
    np.add.at(s, (yi, xi), val); np.add.at(n, (yi, xi), 1)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan)


vocc_rows = []
for esm in ESMS:
    stack = []
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        stack.append(grid2(a[:, cols["longitude"]], a[:, cols["latitude"]], a[:, cols["Temperature"]]))
    T = np.stack(stack)                                            # decade x lat x lon
    yrs = np.array(YEARS, float)
    Tm = T.mean(0)
    # temporal trend, degC per decade
    ok = np.isfinite(T).all(0)
    trend = np.full(Tm.shape, np.nan)
    y = yrs - yrs.mean()
    trend[ok] = ((T[:, ok] - T[:, ok].mean(0)) * y[:, None]).sum(0) / (y ** 2).sum() * 10
    # spatial gradient, degC per km (central differences, longitude wrapped)
    dy_km = 2 * KM
    dTdy = np.full(Tm.shape, np.nan)
    dTdy[1:-1] = (Tm[2:] - Tm[:-2]) / (2 * dy_km)
    dx_km = 2 * KM * np.cos(np.deg2rad(lat2))[:, None]
    dTdx = (np.roll(Tm, -1, axis=1) - np.roll(Tm, 1, axis=1)) / (2 * dx_km)
    grad = np.sqrt(dTdx ** 2 + dTdy ** 2)
    v = np.abs(trend) / grad                                        # km per decade
    w = np.cos(np.deg2rad(lat2))[:, None] * np.ones_like(Tm)
    finite = np.isfinite(v) & np.isfinite(w)
    for label, sel in [("global", np.ones_like(v, bool)),
                       ("30-45 N", (lat2[:, None] >= 30) & (lat2[:, None] < 45) & np.ones_like(v, bool)),
                       ("30-45 S", (lat2[:, None] <= -30) & (lat2[:, None] > -45) & np.ones_like(v, bool))]:
        m = finite & sel
        vals, wt = v[m], w[m]
        order = np.argsort(vals)
        cw = np.cumsum(wt[order]) / wt.sum()
        vocc_rows.append(dict(ESM=esm, domain=label, n_cells=int(m.sum()),
                              median_km_per_decade=float(vals[order][np.searchsorted(cw, 0.5)]),
                              area_weighted_mean=float((vals * wt).sum() / wt.sum()),
                              p25=float(vals[order][np.searchsorted(cw, 0.25)]),
                              p75=float(vals[order][np.searchsorted(cw, 0.75)])))
vocc = pd.DataFrame(vocc_rows)

iso.to_csv(HERE / "cmip6_isotherm_displacement.csv", index=False)
vocc.to_csv(HERE / "cmip6_climate_velocity.csv", index=False)

pd.set_option("display.width", 200)
print("ISOTHERM DISPLACEMENT (same models, grid, decades and front algorithm as the D = 0 fronts)")
print(iso.round(2).to_string(index=False))
print("\nmean by hemisphere:")
print(iso.groupby("hemisphere").km_per_decade.agg(["mean", "min", "max"]).round(1).to_string())
print("\nCLIMATE VELOCITY |dT/dt| / |grad T| (km per decade)")
print(vocc.round(1).to_string(index=False))
print("\nmean across models by domain:")
print(vocc.groupby("domain")[["median_km_per_decade", "area_weighted_mean"]].mean().round(1).to_string())


# ---------------------------------------------------------------- headline statistics for the manuscript
rat = pd.read_csv(HERE / "cmip6_front_rates_by_projection.csv")
d = (rat[rat.projection == "final ensemble (7 predictors)"][["ESM", "hemisphere", "km_per_decade"]]
     .rename(columns={"km_per_decade": "D_front"}))
pair = d.merge(iso[["ESM", "hemisphere", "km_per_decade", "isotherm_degC"]]
               .rename(columns={"km_per_decade": "isotherm"}), on=["ESM", "hemisphere"])
pair["ratio"] = pair.D_front / pair.isotherm
band = vocc[vocc.domain != "global"].assign(
    hemisphere=lambda x: np.where(x.domain.str.endswith("N"), "North", "South"))

rows = []
for hem in ("North", "South"):
    g = pair[pair.hemisphere == hem]
    t = stats.ttest_rel(g.D_front, g.isotherm)
    rows += [dict(statistic=f"D = 0 front rate, {hem}", value=g.D_front.mean(), n=len(g),
                  ci_low=g.D_front.min(), ci_high=g.D_front.max(), P=np.nan, note="four-model mean and range, km per decade"),
             dict(statistic=f"isotherm displacement rate, {hem}", value=g.isotherm.mean(), n=len(g),
                  ci_low=g.isotherm.min(), ci_high=g.isotherm.max(), P=np.nan,
                  note="SST isotherm coinciding with the 2000s front, same grid and algorithm, km per decade"),
             dict(statistic=f"front-to-isotherm rate ratio, {hem}", value=g.ratio.mean(), n=len(g),
                  ci_low=g.ratio.min(), ci_high=g.ratio.max(), P=t.pvalue,
                  note=f"paired t-test over four Earth-system models; ratio < 1 in {(g.ratio < 1).sum()} of {len(g)} models"),
             dict(statistic=f"climate velocity, 30-45 deg {hem[0]}", value=band[band.hemisphere == hem].median_km_per_decade.mean(),
                  n=4, ci_low=band[band.hemisphere == hem].median_km_per_decade.min(),
                  ci_high=band[band.hemisphere == hem].median_km_per_decade.max(), P=np.nan,
                  note="four-model mean of the area-weighted median |dT/dt|/|grad T|, km per decade")]
r = stats.pearsonr(pair.D_front, pair.isotherm)
rows.append(dict(statistic="correlation of front rate with isotherm rate", value=r[0], n=len(pair),
                 ci_low=np.nan, ci_high=np.nan, P=r[1], note="Pearson r across the eight model-hemisphere pairs"))
gv = vocc[vocc.domain == "global"]
rows.append(dict(statistic="climate velocity, global", value=gv.median_km_per_decade.mean(), n=4,
                 ci_low=gv.median_km_per_decade.min(), ci_high=gv.median_km_per_decade.max(), P=np.nan,
                 note="four-model mean of the area-weighted median, km per decade"))
head = pd.DataFrame(rows)
head.to_csv(HERE / "climate_velocity_headline_statistics.csv", index=False)
pair.to_csv(HERE / "cmip6_front_vs_isotherm_by_model.csv", index=False)
print("\nHEADLINE STATISTICS")
print(head.round(3).to_string(index=False))
