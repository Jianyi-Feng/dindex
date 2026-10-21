"""Fronts recomputed from grid cells inside the area of applicability only (unweighted 7-predictor DI, CO3 excluded
because CMIP6 and Tara carbonate differ in units). Appends to cmip6_gridded_statistics.csv."""
import runpy, numpy as np, pandas as pd
from pathlib import Path
HERE = Path(__file__).resolve().parent
# reuse helpers without rerunning the whole analysis
src = (HERE / "cmip6_gridded_analyses.py").read_text()
head = src.split("# ------------------------------------------------------------------ 1 + 2")[0]
aoa_defs = src.split("# ------------------------------------------------------------------ 3: area of applicability")[1].split("AOA_SPECS")[0]
ns = {"__file__": str(HERE / "cmip6_gridded_analyses.py")}
exec(head + aoa_defs, ns)
use = np.array([0, 1, 2, 4, 5, 6, 7]); wts = np.ones(8)
Z, dmean, thr = ns["di_setup"](wts, use)
z, C, ESMS, DECS, YEAR = ns["z"], ns["C"], ns["ESMS"], ns["DECS"], ns["YEAR"]
rows = []
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        Xm = ns["transform"](pd.DataFrame(a[:, [C[v] for v in ns["ENV_M"]]], columns=ns["ENV_M"]), ns["ENV_M"])
        inside = ns["di"](Xm, wts, use, Z, dmean) <= thr
        prof = ns["zonal"](a[inside, C["latitude"]], a[inside, C["D_index_pred"]])
        prof = prof[prof["size"] >= 20]
        for hem, f in ns["fronts"](prof).items():
            rows.append(dict(ESM=esm, decade=dec, year=YEAR[dec], hemisphere=hem, front=f, share_inside=float(inside.mean())))
fr = pd.DataFrame(rows)
rt = ns["rate_table"](fr, "supplied ensemble, AOA cells only")
rt.to_csv(HERE / "cmip6_front_rates_aoa_only.csv", index=False)
fr.to_csv(HERE / "cmip6_front_trajectories_aoa_only.csv", index=False)
st = pd.read_csv(HERE / "cmip6_gridded_statistics.csv")
st = st[~st.statistic.str.contains("AOA cells only")]
new = [dict(section="CMIP6-gridded", statistic=f"front rate {hem}, supplied ensemble, AOA cells only (km/decade)", value=g.km_per_decade.mean(),
            n=len(g), ci_low=g.km_per_decade.min(), ci_high=g.km_per_decade.max(),
            note=f"ESMs poleward {int((g.km_per_decade > 0).sum())}/4; min R2 {g.R2.min():.2f}; DI threshold {thr:.3f}") for hem, g in rt.groupby("hemisphere")]
pd.concat([st, pd.DataFrame(new)]).to_csv(HERE / "cmip6_gridded_statistics.csv", index=False)
print(rt.round(3).to_string()); print(pd.DataFrame(new).to_string())
