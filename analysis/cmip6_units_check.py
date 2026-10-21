"""Compare Tara calibration predictors with co-located 2000-2009 CMIP6 fields (nearest ocean cell)."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
HERE = Path(__file__).resolve().parent
tara = pd.read_csv(HERE.parents[1] / "data/Tara_observed_Dindex.csv")
tara = tara[tara.D_index_qc == "ok"]
st = tara.groupby("station_id")[["latitude", "longitude", "temperature", "salinity", "chlorophyll_a", "CO3", "PO4", "NO3", "Si", "Fe"]].mean()
z = np.load(HERE / "archived__cmip6_grid.npz"); cols = list(z["columns"]); c = {k: i for i, k in enumerate(cols)}
PAIRS = [("temperature", "Temperature"), ("salinity", "Salinity"), ("chlorophyll_a", "ChlorophyllA"), ("CO3", "CO3"),
         ("PO4", "PO4"), ("NO3", "NO3"), ("Si", "Si"), ("Fe", "Iron")]
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
rows = []
for esm in ["CESM2-WACCM", "IPSL-CM6A-LR", "NorESM2-MM", "UKESM1-0-LL"]:
    a = z[f"{esm}__2000-2009"]
    tree = cKDTree(xyz(a[:, c["latitude"]], a[:, c["longitude"]]))
    dist, idx = tree.query(xyz(st.latitude.values, st.longitude.values))
    km = dist * 6371
    for tcol, mcol in PAIRS:
        tv, mv = st[tcol].values, a[idx, c[mcol]]
        ok = (km < 300) & np.isfinite(tv) & np.isfinite(mv)
        pos = ok & (tv > 0) & (mv > 0)
        rows.append(dict(ESM=esm, variable=tcol, n_stations=int(ok.sum()),
                         tara_min=np.nanmin(tv[ok]), tara_median=np.nanmedian(tv[ok]), tara_max=np.nanmax(tv[ok]),
                         cmip6_median_at_stations=np.nanmedian(mv[ok]), cmip6_global_p01=np.nanpercentile(a[:, c[mcol]], 1),
                         cmip6_global_p99=np.nanpercentile(a[:, c[mcol]], 99),
                         median_ratio_cmip6_over_tara=float(np.median(mv[pos] / tv[pos])) if pos.sum() else np.nan,
                         spearman_rho=stats.spearmanr(tv[ok], mv[ok])[0],
                         frac_global_cells_above_tara_max=float(np.mean(a[:, c[mcol]] > np.nanmax(tv[ok]))),
                         frac_global_cells_below_tara_min=float(np.mean(a[:, c[mcol]] < np.nanmin(tv[ok])))))
out = pd.DataFrame(rows)
out.to_csv(HERE / "cmip6_vs_tara_predictor_units.csv", index=False)
pd.set_option("display.width", 250)
print(out.groupby("variable")[["n_stations", "tara_min", "tara_median", "tara_max", "cmip6_median_at_stations", "cmip6_global_p99",
                               "median_ratio_cmip6_over_tara", "spearman_rho", "frac_global_cells_above_tara_max"]].median().round(4).to_string())
