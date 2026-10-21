"""Carbonate sensitivity of the supplied  projection ensemble.

Re-implements Predict_D-index.py exactly (same 14 candidate models, parameter spaces, log10(max(x,0)+1e-6)
transform, StandardScaler, Femto=1 dummy, RandomizedSearchCV n_iter=40 with KFold(3), clipping, per-model D, uniform
average over all models = the selected 'Uniform-All' ensemble, equal pico/femto averaging) and fits the final I and A
models on the 196 Tara metagenomes. Fronts are then computed on 2-degree CMIP6 cells for:
  A  reproduction      : CMIP6 predictors as supplied
  B  CO3 neutralised   : same fitted models, CMIP6 CO3 set to the Tara median in every cell
  C  CO3 removed       : ensemble refitted without CO3 (7 predictors)
The supplied ensemble on the same 2-degree cells is the reference.
"""
from pathlib import Path
import warnings, time
import numpy as np, pandas as pd
from scipy import stats
from sklearn.model_selection import KFold, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.base import clone
from sklearn.ensemble import ExtraTreesRegressor, AdaBoostRegressor, RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, Ridge, BayesianRidge, ElasticNet, HuberRegressor, Lasso, LassoLars
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
RS = 42; NJ = 2; KM = 111.32
ESMS = ["CESM2-WACCM", "IPSL-CM6A-LR", "NorESM2-MM", "UKESM1-0-LL"]
DECS = [f"{y}-{y + 9}" for y in range(2000, 2100, 10)]
LIN = ["Temperature", "Salinity"]; LOG = ["ChlorophyllA", "CO3", "PO4", "NO3", "Si", "Iron"]
log10_safe = lambda x: np.log10(np.maximum(x, 0) + 1e-6)

t = pd.read_csv(HERE.parents[1] / "data/Tara_observed_Dindex.csv")
t = t[(t.I_relative + t.A_relative) > 1e-12].reset_index(drop=True)
t = t.rename(columns={"temperature": "Temperature", "salinity": "Salinity", "chlorophyll_a": "ChlorophyllA", "Fe": "Iron"})
I_vals, A_vals = t.I_relative.to_numpy(float), t.A_relative.to_numpy(float)
femto_tr = (t.size_fraction == "Femto-size").astype(float).to_numpy()
co3_median = float(t.CO3.median())


def models():
    return [("LinearReg", LinearRegression(), None),
            ("Ridge", Ridge(random_state=RS), {"alpha": [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]}),
            ("Lasso", Lasso(max_iter=20000, random_state=RS), {"alpha": [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1.0]}),
            ("ElasticNet", ElasticNet(max_iter=20000, random_state=RS), {"alpha": [1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1], "l1_ratio": [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]}),
            ("BayesianRidge", BayesianRidge(), {"alpha_1": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3], "alpha_2": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3], "lambda_1": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3], "lambda_2": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]}),
            ("Huber", HuberRegressor(max_iter=5000, alpha=1e-3, epsilon=1.35), {"epsilon": [1.1, 1.35, 1.5, 2.0], "alpha": [1e-4, 1e-3, 1e-2, 0.1, 1.0]}),
            ("PLS", PLSRegression(scale=False), {"n_components": [1, 2, 3, 4, 5, 6, 8]}),
            ("LassoLars", LassoLars(max_iter=20000), {"alpha": [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1.0]}),
            ("AdaBoost", AdaBoostRegressor(random_state=RS), {"n_estimators": [50, 100, 200], "learning_rate": [0.1, 0.5, 1.0]}),
            ("ExtraTrees", ExtraTreesRegressor(random_state=RS, n_jobs=1), {"n_estimators": [200, 300, 400, 500], "max_depth": [3, 4, 5, 6, 8], "min_samples_split": [12, 15, 20, 30, 40], "min_samples_leaf": [8, 10, 12, 16, 20], "max_features": [0.3, 0.5, "sqrt", "log2"], "bootstrap": [True, False]}),
            ("RandomForest", RandomForestRegressor(random_state=RS, n_jobs=1), {"n_estimators": [200, 300, 500], "max_depth": [3, 4, 5, 6, 8], "min_samples_split": [12, 15, 20, 30], "min_samples_leaf": [8, 10, 12, 16], "max_features": [0.3, 0.5, "sqrt"]}),
            ("GradientBoosting", GradientBoostingRegressor(random_state=RS), {"n_estimators": [100, 200, 300], "learning_rate": [0.05, 0.1, 0.2], "max_depth": [2, 3, 4], "subsample": [0.8, 1.0]}),
            ("SVR", SVR(), {"C": [0.1, 1.0, 10.0], "epsilon": [0.05, 0.1, 0.2]}),
            ("KNN", KNeighborsRegressor(), {"n_neighbors": [3, 5, 8, 12], "weights": ["uniform", "distance"]})]


def fit_ensemble(log_feats):
    feats = LIN + [f + "_log" for f in log_feats]
    X_num = np.column_stack([t[f].to_numpy(float) for f in LIN] + [log10_safe(t[f].to_numpy(float)) for f in log_feats])
    X = np.hstack([X_num, femto_tr[:, None]]); n_num = X_num.shape[1]
    out = []
    for name, est, ps in models():
        fitted = []
        for y in (I_vals, A_vals):
            pipe = Pipeline([("prep", ColumnTransformer([("num", StandardScaler(), list(range(n_num))), ("cat", "passthrough", [n_num])])),
                             ("model", clone(est))])
            if ps:
                pipe = RandomizedSearchCV(pipe, {"model__" + k: v for k, v in ps.items()}, n_iter=40, cv=KFold(3, shuffle=True, random_state=RS),
                                          scoring="r2", n_jobs=NJ, random_state=RS, refit=True).fit(X, y).best_estimator_
            else:
                pipe.fit(X, y)
            fitted.append(pipe)
        out.append((name, *fitted))
        print("  fitted", name, flush=True)
    return out, n_num


def predict_D(ens, Xn):
    Ds = []
    for fr in (0.0, 1.0):
        X = np.hstack([Xn, np.full((len(Xn), 1), fr)])
        D = np.zeros(len(Xn))
        for _, mi, ma in ens:
            Ip = np.clip(np.ravel(mi.predict(X)), 0, 1); Ap = np.clip(np.ravel(ma.predict(X)), 0, 1)
            den = Ip + Ap
            D += np.divide(Ip - Ap, den, out=np.zeros_like(Ip), where=den > 1e-12) / len(ens)
        Ds.append(np.clip(D, -1, 1))
    return (Ds[0] + Ds[1]) / 2


def fronts(lat, D, lo=15, hi=65):
    prof = pd.DataFrame({"lat": np.floor(lat / 2) * 2 + 1, "D": D}).groupby("lat").D.mean()
    res = {}
    for hem, s in [("North", 1), ("South", -1)]:
        p = prof[(s * prof.index >= lo) & (s * prof.index <= hi)]
        la = p.index.to_numpy() * s; d = p.to_numpy(); o = np.argsort(la); la, d = la[o], d[o]; f = np.nan
        for i in range(len(la) - 1):
            if d[i] > 0 >= d[i + 1]:
                f = la[i] + (la[i + 1] - la[i]) * d[i] / (d[i] - d[i + 1]); break
        res[hem] = f
    return res


z = np.load(HERE / "archived__cmip6_grid.npz"); C = {k: i for i, k in enumerate(z["columns"])}
t0 = time.time()
print("fitting 8-predictor ensemble", flush=True)
ens8, _ = fit_ensemble(LOG)
print("fitting 7-predictor ensemble (no CO3)", flush=True)
ens7, _ = fit_ensemble([f for f in LOG if f != "CO3"])
print(f"fits done in {time.time() - t0:.0f}s", flush=True)

rows, agree_rows = [], []
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        cell = pd.DataFrame(a[:, [C[c] for c in ["latitude", "longitude", *LIN, *LOG, "D_index_pred"]]],
                            columns=["latitude", "longitude", *LIN, *LOG, "Dsup"])
        cell["la"] = np.floor(cell.latitude / 2) * 2 + 1; cell["lo"] = np.floor(cell.longitude / 2) * 2 + 1
        cell = cell.groupby(["la", "lo"]).mean().reset_index()
        Xn8 = np.column_stack([cell[f].to_numpy(float) for f in LIN] + [log10_safe(cell[f].to_numpy(float)) for f in LOG])
        Xn8_neut = Xn8.copy(); Xn8_neut[:, 2 + LOG.index("CO3")] = log10_safe(np.array([co3_median]))[0]
        Xn7 = np.delete(Xn8, 2 + LOG.index("CO3"), axis=1)
        variants = {"supplied ensemble, 2-degree cells": cell.Dsup.to_numpy(),
                    " refit, CO3 as supplied": predict_D(ens8, Xn8),
                    " refit, CO3 set to Tara median": predict_D(ens8, Xn8_neut),
                    " refit without CO3": predict_D(ens7, Xn7)}
        for v, D in variants.items():
            for hem, f in fronts(cell.la.to_numpy(), D).items():
                rows.append(dict(variant=v, ESM=esm, decade=dec, year=int(dec[:4]) + 4.5, hemisphere=hem, front=f))
        ref = variants["supplied ensemble, 2-degree cells"]
        for v, D in variants.items():
            agree_rows.append(dict(variant=v, ESM=esm, decade=dec, r_with_supplied=float(np.corrcoef(ref, D)[0, 1]),
                                   mean_abs_diff=float(np.abs(ref - D).mean())))
        print(esm, dec, f"{time.time() - t0:.0f}s", flush=True)

fr = pd.DataFrame(rows); fr.to_csv(HERE / "carbonate_front_trajectories.csv", index=False)
pd.DataFrame(agree_rows).to_csv(HERE / "carbonate_cellwise_agreement.csv", index=False)
rates = []
for (v, esm, hem), g in fr.dropna().groupby(["variant", "ESM", "hemisphere"]):
    lr = stats.linregress(g.year, g.front)
    rates.append(dict(variant=v, ESM=esm, hemisphere=hem, km_per_decade=lr.slope * KM * 10, R2=lr.rvalue ** 2, front_2000s=g.front.iloc[0]))
rates = pd.DataFrame(rates); rates.to_csv(HERE / "carbonate_front_rates.csv", index=False)
pd.set_option("display.width", 200)
print(rates.groupby(["variant", "hemisphere"]).km_per_decade.agg(["mean", "min", "max", lambda x: int((x > 0).sum())]).round(1))
print(pd.DataFrame(agree_rows).groupby("variant")[["r_with_supplied", "mean_abs_diff"]].mean().round(3))
