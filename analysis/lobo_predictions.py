"""Observation-level leave-one-basin-out predictions for Supplementary Data 6.

Re-runs the leave-one-basin-out block of tara_benchmarks.py with identical data preparation, models and
seeds, and exports every held-out prediction. The run is accepted only if the pooled R2 of each model
reproduces the value stored in tara_cross_validation.csv (which feeds Fig. 3d and Supplementary Table 5).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
EPS = 1e-3
t = pd.read_csv(HERE.parents[1] / "data/Tara_observed_Dindex.csv")
t = t[t.D_index_qc == "ok"].copy().reset_index(drop=True)
t["femto"] = (t.size_fraction == "Femto-size").astype(int)
t["L"] = np.log((t.I_relative + EPS) / (t.A_relative + EPS))
t["log_chl"] = np.log10(t.chlorophyll_a)
for col in ["PO4", "NO3", "Si", "Fe"]:
    t["log_" + col] = np.log10(t[col].clip(lower=t[col][t[col] > 0].min() / 2))
FEAT = ["temperature", "salinity", "log_chl", "log_PO4", "log_NO3", "log_Si", "log_Fe", "femto"]
yD = t.D_index.to_numpy(); yL = t.L.to_numpy(); basins = t.ocean_basin.to_numpy()
MODELS = {
    "temperature-only D": ("D", ["temperature", "femto"], lambda: LinearRegression()),
    "temperature-only log-ratio": ("L", ["temperature", "femto"], lambda: LinearRegression()),
    "ridge log-ratio (7 predictors)": ("L", FEAT, lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25)))),
    "random forest log-ratio (7 predictors)": ("L", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
    "random forest D (7 predictors)": ("D", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
}
to_D = lambda p, tg: np.tanh(p / 2) if tg == "L" else np.clip(p, -1, 1)

stored = pd.read_csv(HERE / "tara_cross_validation.csv")
stored = stored[stored.scheme == "leave-one-basin-out"].set_index("model").R2
out = t[["sample_id", "station_id", "ocean_basin", "size_fraction", "D_index"]].rename(columns={"D_index": "observed_D"})
for name, (tg, feats, mk) in MODELS.items():
    pred = np.full(len(t), np.nan)
    for b in np.unique(basins):
        tr, te = np.where(basins != b)[0], np.where(basins == b)[0]
        mdl = mk().fit(t.iloc[tr][feats].to_numpy(), (yL if tg == "L" else yD)[tr])
        pred[te] = to_D(mdl.predict(t.iloc[te][feats].to_numpy()), tg)
    r2 = r2_score(yD, pred)
    assert abs(r2 - stored[name]) < 1e-9, f"{name}: pooled R2 {r2:.6f} does not reproduce stored {stored[name]:.6f}"
    out[f"LOBO_pred_D__{name}"] = pred
    print(f"{name:42s} pooled LOBO R2 {r2:.4f} (stored {stored[name]:.4f}) reproduced")
out.to_csv(HERE / "lobo_predictions.csv", index=False)
