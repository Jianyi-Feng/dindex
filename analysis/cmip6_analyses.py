"""Gridded analyses of the  (no-'CO3') projection: fronts and rates per ESM, 2-degree maps and sign agreement, area of
applicability (7 predictors, unchanged) near the  fronts, and fronts from applicable cells only. Carbonate-free designs from
cmip6_gridded_analyses.py (temperature-only, random forest without CO3) and the archived  projection are merged for comparison."""
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
src = (HERE / "cmip6_gridded_analyses.py").read_text()
head = src.split("# ------------------------------------------------------------------ 1 + 2")[0]
aoa_defs = src.split("# ------------------------------------------------------------------ 3: area of applicability")[1].split("AOA_SPECS")[0]
ns = {"__file__": str(HERE / "cmip6_gridded_analyses.py")}; exec(head + aoa_defs, ns)
zonal, fronts, rate_table, transform, di_setup, di = (ns[k] for k in ["zonal", "fronts", "rate_table", "transform", "di_setup", "di"])
ESMS, DECS, YEAR, ENV_M = ns["ESMS"], ns["DECS"], ns["YEAR"], ns["ENV_M"]
z = np.load(HERE / "cmip6_grid.npz"); C = {k: i for i, k in enumerate(z["columns"])}
rows = []
rec = lambda name, value, n=None, lo=None, hi=None, note="": rows.append(dict(section="CMIP6-", statistic=name, value=value, n=n, ci_low=lo, ci_high=hi, note=note))
use = np.array([0, 1, 2, 4, 5, 6, 7]); wts = np.ones(8)
Z, dmean, thr = di_setup(wts, use)
thr_rng = ns["AOA_THRESHOLD_RANGE"]
print("AOA threshold min/median/max over 50 partitions:", thr_rng, flush=True)
prof_rows, fr_rows, fr_aoa, aoa_rows, grid, aoa_maps = [], [], [], [], {}, {}
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]; ok = np.isfinite(a[:, C["D_index_pred"]]); a = a[ok]
        lat, lon, D = a[:, C["latitude"]], a[:, C["longitude"]], a[:, C["D_index_pred"]]
        zp = zonal(lat, D)
        prof_rows += [dict(ESM=esm, decade=dec, year=YEAR[dec], lat=la, D=r["mean"], n=r["size"]) for la, r in zp.iterrows()]
        fr_rows += [dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=h, front=f) for h, f in fronts(zp).items()]
        dval = di(transform(pd.DataFrame(a[:, [C[v] for v in ENV_M]], columns=ENV_M), ENV_M), wts, use, Z, dmean)
        inside = dval <= thr; inside_cons = dval <= thr_rng[0]
        pz = zonal(lat[inside], D[inside]); pz = pz[pz["size"] >= 20]
        fr_aoa += [dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=h, front=f) for h, f in fronts(pz).items()]
        w = np.cos(np.radians(lat)); near = np.abs(D) < 0.15; polar = np.abs(lat) >= 60
        aoa_rows.append(dict(spec="unweighted, 7 predictors (no CO3)", ESM=esm, decade=dec, threshold=thr, share_all=float((inside * w).sum() / w.sum()),
                             share_near_front=float((inside * w)[near].sum() / w[near].sum()), share_polar=float((inside * w)[polar].sum() / w[polar].sum()),
                             conservative_threshold=thr_rng[0], share_all_conservative=float((inside_cons * w).sum() / w.sum()),
                             share_near_front_conservative=float((inside_cons * w)[near].sum() / w[near].sum()),
                             share_polar_conservative=float((inside_cons * w)[polar].sum() / w[polar].sum())))
        gi = ((lat + 90) // 2).astype(int).clip(0, 89); gj = ((lon + 180) // 2).astype(int).clip(0, 179)
        s = np.zeros((90, 180)); n = np.zeros((90, 180)); si = np.zeros((90, 180))
        np.add.at(s, (gi, gj), D); np.add.at(n, (gi, gj), 1); np.add.at(si, (gi, gj), inside)
        with np.errstate(invalid="ignore"):
            grid[(esm, dec)] = np.where(n > 0, s / n, np.nan); aoa_maps[(esm, dec)] = np.where(n > 0, si / n, np.nan)
        print(esm, dec, flush=True)
pd.DataFrame(prof_rows).to_csv(HERE / "cmip6_zonal_profiles.csv", index=False)
fr = pd.DataFrame(fr_rows); fa = pd.DataFrame(fr_aoa)
rt = rate_table(fr, "final ensemble (7 predictors)"); ra = rate_table(fa, "final ensemble, AOA cells only")
keep = {"supplied ensemble (gridded)": "archived final ensemble (8 predictors incl. mislabelled CO3)",
        "temperature-only": "temperature-only", "random forest, 7 predictors (no CO3)": "random forest, 7 predictors"}
# 早期存档结果只作为一致性核对，缺失时跳过
_old_rates = HERE / "archived__cmip6_front_rates_by_projection.csv"
_old_traj = HERE / "archived__cmip6_front_trajectories.csv"
if _old_rates.exists():
    _o = pd.read_csv(_old_rates)
    _o = _o[_o.projection.isin(keep)].assign(projection=lambda d: d.projection.map(keep))
    rates = pd.concat([rt, ra, _o])
else:
    rates = pd.concat([rt, ra])
rates.to_csv(HERE / "cmip6_front_rates_by_projection.csv", index=False)
_parts = [fr.assign(projection="final ensemble (7 predictors)"),
          fa.assign(projection="final ensemble, AOA cells only")]
if _old_traj.exists():
    _ot = pd.read_csv(_old_traj)
    _parts.append(_ot[_ot.projection.isin(keep)].assign(projection=lambda d: d.projection.map(keep)))
pd.concat(_parts).to_csv(HERE / "cmip6_front_trajectories_by_projection.csv", index=False)
for (p, h), g in rates.groupby(["projection", "hemisphere"]):
    rec(f"front rate {h}, {p} (km/decade)", g.km_per_decade.mean(), 4, g.km_per_decade.min(), g.km_per_decade.max(),
        f"ESMs poleward {int((g.km_per_decade > 0).sum())}/4; min R2 {g.R2.min():.2f}")
for h, g in rt.groupby("hemisphere"):
    d = (g.front_2090s - g.front_2000s) * 111.32
    rec(f" front displacement 2000s->2090s {h} (km)", d.mean(), 4, d.min(), d.max())
    rec(f" front latitude 2000s {h}", g.front_2000s.mean(), 4, g.front_2000s.min(), g.front_2000s.max())
G0 = np.stack([grid[(e, "2000-2009")] for e in ESMS]); G1 = np.stack([grid[(e, "2090-2099")] for e in ESMS]); dG = G1 - G0
common = np.all(np.isfinite(dG), 0); agree = common & (np.all(dG > 0, 0) | np.all(dG < 0, 0))
latc = -89 + 2 * np.arange(90); cos = np.cos(np.radians(latc))[:, None] * np.ones((1, 180))
rec("area-weighted share with four-model sign agreement on dD", float((agree * cos).sum() / (common * cos).sum()))
with np.errstate(invalid="ignore"):
    rec("area-weighted share where four-model mean dD > 0", float((((np.nanmean(dG, 0) > 0) & common) * cos).sum() / (common * cos).sum()))
band_rows = []
for lo, hi, lab in [(0, 30, "|lat| < 30"), (30, 45, "30-45"), (45, 60, "45-60"), (60, 90, "> 60")]:
    for hem, sg in [("both", None), ("North", 1), ("South", -1)]:
        sel = (np.abs(latc) >= lo) & (np.abs(latc) < hi)
        if sg is not None:
            sel &= np.sign(latc) == sg
        mk = common & sel[:, None]
        if mk.sum() == 0:
            continue
        share = float(((agree & mk) * cos).sum() / (mk * cos).sum())
        with np.errstate(invalid="ignore"):
            mdd = float(np.nansum(np.nanmean(dG, 0) * cos * mk) / (mk * cos).sum())
        band_rows.append(dict(band=lab, hemisphere=hem, n_cells=int(mk.sum()), share_sign_agreement=share, mean_dD=mdd))
        if hem == "both":
            rec(f"sign agreement, {lab}", share, int(mk.sum()), note=f"mean dD {mdd:.3f}")
pd.DataFrame(band_rows).to_csv(HERE / "cmip6_agreement_by_band.csv", index=False)
aoa = pd.DataFrame(aoa_rows); aoa.to_csv(HERE / "cmip6_area_of_applicability.csv", index=False)
for dec in ["2000-2009", "2090-2099"]:
    s = aoa[aoa.decade == dec]
    for col, lab in [("share_all", "all ocean"), ("share_near_front", "near fronts |D|<0.15"), ("share_polar", "|lat| >= 60"),
                     ("share_all_conservative", "all ocean, conservative threshold"), ("share_near_front_conservative", "near fronts, conservative threshold"),
                     ("share_polar_conservative", "|lat| >= 60, conservative threshold")]:
        rec(f"AOA share {dec}, {lab}", s[col].mean(), 4, s[col].min(), s[col].max())
with np.errstate(invalid="ignore"):
    np.savez_compressed(HERE / "cmip6_maps_2deg.npz", D2000=G0, D2090=G1, esms=np.array(ESMS))
    np.savez_compressed(HERE / "cmip6_aoa_2deg.npz", aoa2000=np.nanmean(np.stack([aoa_maps[(e, "2000-2009")] for e in ESMS]), 0),
                        aoa2090=np.nanmean(np.stack([aoa_maps[(e, "2090-2099")] for e in ESMS]), 0))
_old_grid_path = HERE / "archived__cmip6_grid.npz"
old_grid = np.load(_old_grid_path) if _old_grid_path.exists() else None
cmp = []
for esm in ESMS:
    for dec in ["2000-2009", "2090-2099"]:
        a_new = z[f"{esm}__{dec}"][:, C["D_index_pred"]]; a_old = old_grid[f"{esm}__{dec}"][:, C["D_index_pred"]]; ok = np.isfinite(a_new) & np.isfinite(a_old)
        cmp.append(dict(ESM=esm, decade=dec, r=float(np.corrcoef(a_new[ok], a_old[ok])[0, 1]), mean_abs_diff=float(np.abs(a_new[ok] - a_old[ok]).mean())))
cmp = pd.DataFrame(cmp); cmp.to_csv(HERE / "cmip6_vs_cellwise.csv", index=False)
rec("AOA dissimilarity threshold (median of 50 station partitions)", thr, 50, thr_rng[0], thr_rng[2])
rec(" vs archived  grid-cell D: mean r", cmp.r.mean(), len(cmp), cmp.r.min(), cmp.r.max(), f"mean |dD| {cmp.mean_abs_diff.mean():.3f}")
out = pd.DataFrame(rows); out.to_csv(HERE / "cmip6_statistics.csv", index=False)
pd.set_option("display.width", 220, "display.max_colwidth", 80)
print(out.drop(columns="section").round(3).to_string())
