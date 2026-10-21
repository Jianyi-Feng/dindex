""" Tara benchmark validation without the mislabelled 'CO3' predictor (7 environmental predictors + size fraction).
Same design as  dindex_analyses.py section D: 20 random 5-fold station partitions, leave-one-basin-out, and out-of-fold
permutation importance (random-forest log-ratio). Also summarizes the final ensemble nested CV (nested_cv_without_pigment.py)."""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260916); EPS = 1e-3
t = pd.read_csv(HERE.parents[1] / "data/Tara_observed_Dindex.csv")
t = t[t.D_index_qc == "ok"].copy().reset_index(drop=True)
t["femto"] = (t.size_fraction == "Femto-size").astype(int)
t["L"] = np.log((t.I_relative + EPS) / (t.A_relative + EPS))
t["log_chl"] = np.log10(t.chlorophyll_a)
for col in ["PO4", "NO3", "Si", "Fe"]:
    t["log_" + col] = np.log10(t[col].clip(lower=t[col][t[col] > 0].min() / 2))
FEAT = ["temperature", "salinity", "log_chl", "log_PO4", "log_NO3", "log_Si", "log_Fe", "femto"]
yD = t.D_index.to_numpy(); yL = t.L.to_numpy(); groups = t.station_id.to_numpy(); basins = t.ocean_basin.to_numpy()
MODELS = {
    "temperature-only D": ("D", ["temperature", "femto"], lambda: LinearRegression()),
    "temperature-only log-ratio": ("L", ["temperature", "femto"], lambda: LinearRegression()),
    "ridge log-ratio (7 predictors)": ("L", FEAT, lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25)))),
    "random forest log-ratio (7 predictors)": ("L", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
    "random forest D (7 predictors)": ("D", FEAT, lambda: RandomForestRegressor(500, min_samples_leaf=3, random_state=1, n_jobs=-1)),
}
to_D = lambda p, tg: np.tanh(p / 2) if tg == "L" else np.clip(p, -1, 1)


def cv_run(splits, perm_imp=False):
    out, imp = {}, []
    for name, (tg, feats, mk) in MODELS.items():
        pred = np.full(len(t), np.nan)
        for tr, te in splits():
            Xtr, Xte = t.iloc[tr][feats].to_numpy(), t.iloc[te][feats].to_numpy()
            mdl = mk().fit(Xtr, (yL if tg == "L" else yD)[tr]); pred[te] = to_D(mdl.predict(Xte), tg)
            if perm_imp and name == "random forest log-ratio (7 predictors)":
                base = r2_score(yD[te], pred[te])
                for j, f in enumerate(feats):
                    for _ in range(10):
                        Xp = Xte.copy(); Xp[:, j] = RNG.permutation(Xp[:, j])
                        imp.append(dict(feature=f, dR2=base - r2_score(yD[te], to_D(mdl.predict(Xp), tg))))
        out[name] = pred
    return out, imp


def station_folds(seed, k=5):
    def it():
        us = np.random.default_rng(seed).permutation(np.unique(groups)); fo = {s: i % k for i, s in enumerate(us)}
        f = np.array([fo[s] for s in groups])
        for i in range(k):
            yield np.where(f != i)[0], np.where(f == i)[0]
    return it


rows, imps = [], []
for rep in range(20):
    res, imp = cv_run(station_folds(rep), perm_imp=True); imps += imp
    for n, p in res.items():
        rows.append(dict(scheme="5-fold station-blocked", rep=rep, model=n, R2=r2_score(yD, p), RMSE=mean_squared_error(yD, p) ** 0.5, MAE=mean_absolute_error(yD, p)))
    print("partition", rep, flush=True)
res_b, _ = cv_run(lambda: ((np.where(basins != b)[0], np.where(basins == b)[0]) for b in np.unique(basins)))
for n, p in res_b.items():
    rows.append(dict(scheme="leave-one-basin-out", rep=0, model=n, R2=r2_score(yD, p), RMSE=mean_squared_error(yD, p) ** 0.5, MAE=mean_absolute_error(yD, p)))
    for b in np.unique(basins):
        m = basins == b
        rows.append(dict(scheme=f"LOBO held-out basin: {b}", rep=0, model=n, n=int(m.sum()), R2=r2_score(yD[m], p[m]) if m.sum() > 2 else np.nan,
                         RMSE=mean_squared_error(yD[m], p[m]) ** 0.5, MAE=mean_absolute_error(yD[m], p[m]), bias=float(np.mean(p[m] - yD[m]))))
cv = pd.DataFrame(rows); cv.to_csv(HERE / "tara_cross_validation.csv", index=False)
im = pd.DataFrame(imps).groupby("feature").dR2.agg(["mean", "std", "count"]).sort_values("mean", ascending=False)
im["ci_low"] = pd.DataFrame(imps).groupby("feature").dR2.quantile(.025); im["ci_high"] = pd.DataFrame(imps).groupby("feature").dR2.quantile(.975)
im.to_csv(HERE / "tara_oof_permutation_importance.csv")
oof = pd.read_csv(HERE / "nestedcv_oof_7pred.csv")
ens = dict(model="final ensemble (14 models, I and A separately; 7 predictors), nested station-grouped CV",
           R2=r2_score(oof.observed_D, oof.oof_D), RMSE=mean_squared_error(oof.observed_D, oof.oof_D) ** 0.5, MAE=mean_absolute_error(oof.observed_D, oof.oof_D), n=len(oof))
pd.DataFrame([ens]).to_csv(HERE / "ensemble_validation.csv", index=False)
print(cv[cv.scheme.isin(["5-fold station-blocked", "leave-one-basin-out"])].groupby(["scheme", "model"])[["R2", "RMSE", "MAE"]].mean().round(3))
print(im.round(3)); print(ens)
