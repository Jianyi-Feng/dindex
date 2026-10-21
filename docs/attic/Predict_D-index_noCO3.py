#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D-index Prediction  -  without the Tara "CO3" predictor

 change: the Tara "CO3" column (TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB_20170515, "Carbonate ion") is identical to the HPLC
median of 19'-butanoyloxyfucoxanthin in the same PANGAEA release (mislabelled columns), so it is removed from training and
projection. Everything else is unchanged from .

( header) Framework C (predict I and A separately) + ENSEMBLE + I/A abundance

Framework C: predict Imitervirales (I) and Algavirales (A) separately, then combine D = (I-A)/(I+A).
Follows the earlier conclusion but drops the Temp x log10(Chl) interaction (main-effects model):
  - 8 main-effect features: Temperature and Salinity (linear), plus ChlorophyllA, CO3, PO4, NO3, Si and Iron (log10);
  - model pool: LinearReg, Ridge, Lasso, ElasticNet, BayesianRidge, Huber, PLS, LassoLars, AdaBoost, ExtraTrees (KernelRidge removed);
  - ensemble: average the per-model D predictions (uniform / R2-weighted / top-3) and choose the best by GroupKFold CV.
  - outputs the D-index plus the predicted Imitervirales/Algavirales relative abundances.

CMIP6 scanning covers all decades automatically.
"""
import os, sys, warnings, time, hashlib, argparse, glob, re
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import GroupKFold, KFold, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import (ExtraTreesRegressor, AdaBoostRegressor,
                              RandomForestRegressor, GradientBoostingRegressor)
from sklearn.linear_model import (LinearRegression, Ridge, BayesianRidge,
                                  ElasticNet, HuberRegressor, Lasso, LassoLars)
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.base import clone

parser = argparse.ArgumentParser()
parser.add_argument("--base-dir", type=str, default=None)
parser.add_argument("--cmip6-dir", type=str, default=None)
parser.add_argument("--n-jobs", type=int, default=None)
parser.add_argument("--n-iter", type=int, default=50)
args = parser.parse_args()

BASE_DIR = args.base_dir or os.path.dirname(os.path.abspath(__file__))
TRAIN_FILE = os.path.join(BASE_DIR,
    "Final_SPOT_TARA_ocean_Imitervirales_and_environments.csv")
if args.cmip6_dir:
    CMIP6_DIR = args.cmip6_dir
else:
    prefer = os.path.join(BASE_DIR, "CMIP6_100km_2000-2099_surface_ocean_environmental_factors")
    if os.path.isdir(prefer):
        CMIP6_DIR = prefer
    else:
        cands = [c for c in sorted(glob.glob(os.path.join(BASE_DIR, "*CMIP6*"))) if os.path.isdir(c)]
        CMIP6_DIR = cands[0] if cands else prefer
OUTPUT_DIR = os.path.join(BASE_DIR, "D_index_predictions_output_noCO3")
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

RANDOM_STATE = 42
N_JOBS = args.n_jobs if args.n_jobs else max(1, min(32, os.cpu_count() or 4))
N_JOBS_SEARCH = min(N_JOBS, 48)
os.environ["OMP_NUM_THREADS"] = str(min(N_JOBS, 8))
os.environ["OPENBLAS_NUM_THREADS"] = str(min(N_JOBS, 8))
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300,
                     "font.size": 10, "font.family": "sans-serif"})

LINEAR_FEATURES = ["Temperature", "Salinity"]
LOG_FEATURES = ["ChlorophyllA", "PO4", "NO3", "Si", "Iron"]  # : CO3 removed (mislabelled pigment)
DERIVED_FEATURES = []  # main-effects model: no interaction term (Temp_x_Chl_log removed)
CAT_FEATURE = "Size_fraction_group"


def log10_safe(x):
    return np.log10(np.maximum(x, 0) + 1e-6)


print("=" * 70)
print("D-index Prediction  - Framework C (main effects, no interaction, CO3 removed)")
print("=" * 70)
print(f"N_JOBS={N_JOBS}  n_iter={args.n_iter}")
print(f"CMIP6 directory: {CMIP6_DIR}")

# ── data ───────────────────────────────────────────────────────────
train_df = pd.read_csv(TRAIN_FILE, encoding="utf-8-sig")
train_df = train_df.loc[:, ~train_df.columns.str.startswith("Unnamed")].copy()
train_df.columns = [str(c).strip().lstrip("\ufeff") for c in train_df.columns]


def _find_col(df, needles):
    if isinstance(needles, str):
        needles = [needles]
    for n in needles:
        m = [c for c in df.columns if n.lower() in str(c).lower()]
        if m:
            return m[0]
    raise KeyError(f"{needles} not found; columns={list(df.columns)}")


ALG_COL = _find_col(train_df, ["algavirales"])
IMI_COL = _find_col(train_df, "imitervirales")

I_raw = train_df[IMI_COL].values.astype(np.float64)
A_raw = train_df[ALG_COL].values.astype(np.float64)
denom = I_raw + A_raw
valid = denom > 1e-12
train_df["D_index"] = np.divide(I_raw - A_raw, denom,
                                out=np.full_like(I_raw, np.nan), where=valid)
train_df = train_df[valid].reset_index(drop=True)

for feat in LOG_FEATURES:
    train_df[feat + "_log"] = log10_safe(train_df[feat].values)
model_features = LINEAR_FEATURES + [f + "_log" for f in LOG_FEATURES] + DERIVED_FEATURES

env_cols = ["Temperature", "ChlorophyllA", "Salinity",
            "CO3", "PO4", "NO3", "Si", "Iron"]  # grouping key unchanged from  (maps 1:1 to stations)
env_str = train_df[env_cols].astype(str).agg("_".join, axis=1)
train_df["station_id"] = env_str.apply(lambda x: hashlib.md5(x.encode()).hexdigest()[:12])

I_vals = train_df[IMI_COL].values.astype(np.float64)
A_vals = train_df[ALG_COL].values.astype(np.float64)
y_D = train_df["D_index"].values.astype(np.float64)
groups = train_df["station_id"].values
n_stations = train_df["station_id"].nunique()
print(f"Stations={n_stations}  Samples={len(train_df)}  D mean={y_D.mean():.4f} std={y_D.std():.4f}")

X_num = train_df[model_features].values
# Manual binarisation of Size_fraction: Pico=0, Femto=1 (equivalent to OneHotEncoder drop=first,
# and avoids the isnan error of the string OneHotEncoder in some sklearn/numpy versions)
X_cat = (train_df[CAT_FEATURE].values == "Femto-size").astype(np.float64).reshape(-1, 1)
X_full = np.hstack([X_num, X_cat])
n_num_cols = X_num.shape[1]


def make_prep():
    return ColumnTransformer([
        ("num", StandardScaler(), list(range(n_num_cols))),
        ("cat", "passthrough", [n_num_cols]),
    ])


outer_cv = GroupKFold(n_splits=5)


def get_all_models():
    m = []  # (name, estimator, param_space)
    m.append(("LinearReg", LinearRegression(), None))
    m.append(("Ridge", Ridge(random_state=RANDOM_STATE),
        {"alpha": [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]}))
    m.append(("Lasso", Lasso(max_iter=20000, random_state=RANDOM_STATE),
        {"alpha": [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1.0]}))
    m.append(("ElasticNet", ElasticNet(max_iter=20000, random_state=RANDOM_STATE),
        {"alpha": [1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1],
         "l1_ratio": [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]}))
    m.append(("BayesianRidge", BayesianRidge(),
        {"alpha_1": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
         "alpha_2": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
         "lambda_1": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
         "lambda_2": [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]}))
    m.append(("Huber", HuberRegressor(max_iter=5000, alpha=1e-3, epsilon=1.35),
        {"epsilon": [1.1, 1.35, 1.5, 2.0], "alpha": [1e-4, 1e-3, 1e-2, 0.1, 1.0]}))
    m.append(("PLS", PLSRegression(scale=False),
        {"n_components": [1, 2, 3, 4, 5, 6, 8]}))
    m.append(("LassoLars", LassoLars(max_iter=20000),
        {"alpha": [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1.0]}))
    m.append(("AdaBoost", AdaBoostRegressor(random_state=RANDOM_STATE),
        {"n_estimators": [50, 100, 200], "learning_rate": [0.1, 0.5, 1.0]}))
    m.append(("ExtraTrees", ExtraTreesRegressor(random_state=RANDOM_STATE, n_jobs=1),
        {"n_estimators": [200, 300, 400, 500],
         "max_depth": [3, 4, 5, 6, 8],
         "min_samples_split": [12, 15, 20, 30, 40],
         "min_samples_leaf": [8, 10, 12, 16, 20],
         "max_features": [0.3, 0.5, "sqrt", "log2"],
         "bootstrap": [True, False]}))
    m.append(("RandomForest", RandomForestRegressor(random_state=RANDOM_STATE, n_jobs=1),
        {"n_estimators": [200, 300, 500],
         "max_depth": [3, 4, 5, 6, 8],
         "min_samples_split": [12, 15, 20, 30],
         "min_samples_leaf": [8, 10, 12, 16],
         "max_features": [0.3, 0.5, "sqrt"]}))
    m.append(("GradientBoosting", GradientBoostingRegressor(random_state=RANDOM_STATE),
        {"n_estimators": [100, 200, 300],
         "learning_rate": [0.05, 0.1, 0.2],
         "max_depth": [2, 3, 4],
         "subsample": [0.8, 1.0]}))
    m.append(("SVR", SVR(),
        {"C": [0.1, 1.0, 10.0],
         "epsilon": [0.05, 0.1, 0.2]}))
    m.append(("KNN", KNeighborsRegressor(),
        {"n_neighbors": [3, 5, 8, 12],
         "weights": ["uniform", "distance"]}))
    return m


all_param_spaces = get_all_models()
print(f"Candidates: {len(all_param_spaces)}")


# ── nested CV for one target (returns OOF preds + per-fold R2/gap) ──
def nested_cv_oof(name, mt, ps, y_target):
    oof = np.full(len(y_target), np.nan)
    r2s, gaps = [], []
    for fold_i, (tr_i, te_i) in enumerate(outer_cv.split(X_full, y_target, groups)):
        X_tr, X_te = X_full[tr_i], X_full[te_i]
        y_tr, y_te = y_target[tr_i], y_target[te_i]
        pipe = Pipeline([("prep", make_prep()), ("model", clone(mt))])
        if ps:
            pg = {"model__" + k: v for k, v in ps.items()}
            pipe = RandomizedSearchCV(
                pipe, pg, n_iter=min(args.n_iter, 40),
                cv=KFold(3, shuffle=True, random_state=RANDOM_STATE + fold_i),
                scoring="r2", n_jobs=N_JOBS_SEARCH, random_state=RANDOM_STATE, refit=True
            ).fit(X_tr, y_tr).best_estimator_
        else:
            pipe.fit(X_tr, y_tr)
        pred_tr = pipe.predict(X_tr); pred_te = pipe.predict(X_te)
        oof[te_i] = pred_te
        r2s.append(r2_score(y_te, pred_te))
        gaps.append(r2_score(y_tr, pred_tr) - r2_score(y_te, pred_te))
    return oof, float(np.mean(r2s)), float(np.mean(gaps))


def combine_D(I_pred, A_pred):
    Ip = np.clip(I_pred, 0.0, 1.0)
    Ap = np.clip(A_pred, 0.0, 1.0)
    den = Ip + Ap
    return np.divide(Ip - Ap, den, out=np.zeros_like(Ip), where=den > 1e-12)


def dspace_r2(D_pred):
    r2s = [r2_score(y_D[te], D_pred[te]) for _, te in outer_cv.split(X_full, y_D, groups)]
    return float(np.mean(r2s))


# ── Framework C: per-model I & A -> D ──────────────────────────────
print("\n" + "=" * 70)
print("Framework C: nested CV on I and A separately, combine to D")
print("=" * 70)
results = []
for name, mt, ps in all_param_spaces:
    t0 = time.time()
    oof_I, r2_I, gap_I = nested_cv_oof(name, mt, ps, I_vals)
    oof_A, r2_A, gap_A = nested_cv_oof(name, mt, ps, A_vals)
    D_oof = combine_D(oof_I, oof_A)
    r2_D = dspace_r2(D_oof)
    results.append(dict(name=name, oof_I=oof_I, oof_A=oof_A, D_oof=D_oof,
                        r2_D=r2_D, r2_I=r2_I, r2_A=r2_A,
                        gap=float(np.mean([gap_I, gap_A]))))
    print(f"  [{name:14s}] R2(D)={r2_D:.4f}  (R2_I={r2_I:.3f}, R2_A={r2_A:.3f})  [{time.time()-t0:.0f}s]")

results.sort(key=lambda x: x["r2_D"], reverse=True)
print(f"\n{'Rank':<5s} {'Model':16s} {'R2(D)':>9s}")
print("-" * 34)
for i, r in enumerate(results):
    print(f"{i+1:<5d} {r['name']:16s} {r['r2_D']:9.4f}")

best = results[0]
print(f"\nBest single (Framework C): {best['name']} (R2(D)={best['r2_D']:.4f})")
pd.DataFrame([{"name": r["name"], "r2_D": r["r2_D"], "r2_I": r["r2_I"],
               "r2_A": r["r2_A"], "gap": r["gap"]} for r in results]
             ).to_csv(os.path.join(OUTPUT_DIR, "model_comparison.csv"), index=False)

# ── ensemble on D predictions (average top-N OOF D) ────────────────
TOP_N = len(results)
top = results[:TOP_N]

def cv_ens_D(names_idx, weights=None):
    # average the OOF D predictions of the selected models
    if weights is None:
        w = np.ones(len(names_idx)) / len(names_idx)
    else:
        w = np.array(weights) / np.sum(weights)
    D_ens = np.zeros_like(y_D)
    for j, idx in enumerate(names_idx):
        D_ens += w[j] * top[idx]["D_oof"]
    return dspace_r2(D_ens)


print("\n" + "=" * 70)
print("Ensemble strategies (average of D predictions, GroupKFold CV)")
print("=" * 70)
strategies = []

# uniform over all
r2 = cv_ens_D(list(range(TOP_N)))
strategies.append(("Uniform-All", r2))
print(f"  Uniform-All:       R2(D)={r2:.4f}")

# R2-weighted over all
wr = np.array([max(0.01, r["r2_D"]) for r in top])
r2 = cv_ens_D(list(range(TOP_N)), weights=wr)
strategies.append(("R2-weighted", r2))
print(f"  R2-weighted:       R2(D)={r2:.4f}")

# Top3 uniform
r2 = cv_ens_D([0, 1, 2])
strategies.append(("Top3-Uniform", r2))
print(f"  Top3-Uniform:      R2(D)={r2:.4f}")

# Top5 uniform
r2 = cv_ens_D(list(range(5)))
strategies.append(("Top5-Uniform", r2))
print(f"  Top5-Uniform:      R2(D)={r2:.4f}")

strategies.sort(key=lambda s: s[1], reverse=True)
best_ens_name, best_ens_r2 = strategies[0]
print(f"\nSelected ensemble: {best_ens_name} (R2(D)={best_ens_r2:.4f})")

# determine which models to use for the final ensemble
if best_ens_name == "Uniform-All":
    final_idx = list(range(TOP_N)); final_w = None
elif best_ens_name == "R2-weighted":
    final_idx = list(range(TOP_N)); final_w = wr / wr.sum()
elif best_ens_name == "Top3-Uniform":
    final_idx = [0, 1, 2]; final_w = None
elif best_ens_name == "Top5-Uniform":
    final_idx = list(range(5)); final_w = None
else:
    final_idx = [0]; final_w = None

# ── train the final I/A models on full data ────────────────────────
print("\nTraining final I & A models on full data ...")
final_I_models, final_A_models = {}, {}
for idx in final_idx:
    r = top[idx]
    for name, mt, ps in all_param_spaces:
        if name == r["name"]:
            def fit_full(y_target):
                pipe = Pipeline([("prep", make_prep()), ("model", clone(mt))])
                if ps:
                    pg = {"model__" + k: v for k, v in ps.items()}
                    pipe = RandomizedSearchCV(
                        pipe, pg, n_iter=min(args.n_iter, 40),
                        cv=KFold(3, shuffle=True, random_state=RANDOM_STATE),
                        scoring="r2", n_jobs=N_JOBS_SEARCH, random_state=RANDOM_STATE, refit=True
                    ).fit(X_full, y_target).best_estimator_
                else:
                    pipe.fit(X_full, y_target)
                return pipe
            final_I_models[name] = fit_full(I_vals)
            final_A_models[name] = fit_full(A_vals)
            break

prep = final_I_models[top[final_idx[0]]["name"]].named_steps["prep"]


def model_D_pred(nm, X):
    # NOTE: X is already the output of prep.transform, so only the model step is called here
    I_p = np.clip(final_I_models[nm].named_steps["model"].predict(X), 0, 1)
    A_p = np.clip(final_A_models[nm].named_steps["model"].predict(X), 0, 1)
    den = I_p + A_p
    return np.divide(I_p - A_p, den, out=np.zeros_like(I_p), where=den > 1e-12)


# full-train R2 (in-sample, for reference) -- computed on the standardised X
X_full_pp = prep.transform(X_full)
D_full = np.zeros_like(y_D)
for j, idx in enumerate(final_idx):
    w = (1.0 / len(final_idx)) if final_w is None else final_w[j]
    D_full += w * model_D_pred(top[idx]["name"], X_full_pp)
D_full = np.clip(D_full, -1, 1)
print(f"Full-training R2(D) = {r2_score(y_D, D_full):.4f}")

# ── CMIP6 ──────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("CMIP6 Global Prediction (Framework C: I & A separate)")
print("=" * 70)
CMIP6_FILES = sorted(glob.glob(os.path.join(CMIP6_DIR, "*_data_transformed.csv")))
print(f"Found {len(CMIP6_FILES)} CMIP6 file(s)")

CMIP6_LINEAR = ["Temperature", "Salinity"]
CMIP6_LOG = ["ChlorophyllA", "PO4", "NO3", "Si", "Iron"]  # : CO3 removed
OUTPUT_COLS = ["longitude", "latitude",
               "Temperature", "ChlorophyllA", "Salinity",
               "CO3", "PO4", "NO3", "Si", "Iron",
               "D_index_pred_Pico", "D_index_pred_Femto", "D_index_pred",
               "Imitervirales_pred_Pico", "Imitervirales_pred_Femto", "Imitervirales_pred",
               "Algavirales_pred_Pico", "Algavirales_pred_Femto", "Algavirales_pred",
               "Model", "Period"]


def predict_D(X_prep):
    D_ens = np.zeros(X_prep.shape[0])
    for j, idx in enumerate(final_idx):
        w = (1.0 / len(final_idx)) if final_w is None else final_w[j]
        D_ens += w * model_D_pred(top[idx]["name"], X_prep)
    return np.clip(D_ens, -1, 1)


def predict_I(X_prep):
    # Imitervirales relative abundance (0-1 fraction), weighted average of the ensemble I models
    I_ens = np.zeros(X_prep.shape[0])
    for j, idx in enumerate(final_idx):
        w = (1.0 / len(final_idx)) if final_w is None else final_w[j]
        I_ens += w * np.clip(final_I_models[top[idx]["name"]].named_steps["model"].predict(X_prep), 0, 1)
    return I_ens


def predict_A(X_prep):
    # Algavirales relative abundance (0-1 fraction), weighted average of the ensemble A models
    A_ens = np.zeros(X_prep.shape[0])
    for j, idx in enumerate(final_idx):
        w = (1.0 / len(final_idx)) if final_w is None else final_w[j]
        A_ens += w * np.clip(final_A_models[top[idx]["name"]].named_steps["model"].predict(X_prep), 0, 1)
    return A_ens


all_preds = []
for fpath in CMIP6_FILES:
    fname = os.path.basename(fpath)
    print(f"  {fname}...", end=" ", flush=True)
    t0 = time.time()
    cmip = pd.read_csv(fpath)
    for feat in CMIP6_LOG:
        if feat in cmip.columns:
            cmip[feat + "_log"] = log10_safe(cmip[feat].values)
    cmip_mf = CMIP6_LINEAR + [f + "_log" for f in CMIP6_LOG] + DERIVED_FEATURES
    cmip_num = cmip[cmip_mf].values
    ocean_mask = np.isfinite(cmip_num).all(axis=1)
    n_ocean = int(ocean_mask.sum())
    cmip["D_index_pred_Pico"] = np.nan
    cmip["D_index_pred_Femto"] = np.nan
    cmip["D_index_pred"] = np.nan
    cmip["Imitervirales_pred_Pico"] = np.nan
    cmip["Imitervirales_pred_Femto"] = np.nan
    cmip["Imitervirales_pred"] = np.nan
    cmip["Algavirales_pred_Pico"] = np.nan
    cmip["Algavirales_pred_Femto"] = np.nan
    cmip["Algavirales_pred"] = np.nan
    if n_ocean > 0:
        Xn = cmip_num[ocean_mask]
        Xp = prep.transform(np.hstack([Xn, np.zeros((n_ocean, 1))]))   # Pico = 0
        Xf = prep.transform(np.hstack([Xn, np.ones((n_ocean, 1))]))    # Femto = 1
        cmip.loc[ocean_mask, "D_index_pred_Pico"] = predict_D(Xp)
        cmip.loc[ocean_mask, "D_index_pred_Femto"] = predict_D(Xf)
        cmip.loc[ocean_mask, "D_index_pred"] = (cmip.loc[ocean_mask, "D_index_pred_Pico"] +
                                                cmip.loc[ocean_mask, "D_index_pred_Femto"]) / 2.0
        cmip.loc[ocean_mask, "Imitervirales_pred_Pico"] = predict_I(Xp)
        cmip.loc[ocean_mask, "Imitervirales_pred_Femto"] = predict_I(Xf)
        cmip.loc[ocean_mask, "Imitervirales_pred"] = (cmip.loc[ocean_mask, "Imitervirales_pred_Pico"] +
                                                       cmip.loc[ocean_mask, "Imitervirales_pred_Femto"]) / 2.0
        cmip.loc[ocean_mask, "Algavirales_pred_Pico"] = predict_A(Xp)
        cmip.loc[ocean_mask, "Algavirales_pred_Femto"] = predict_A(Xf)
        cmip.loc[ocean_mask, "Algavirales_pred"] = (cmip.loc[ocean_mask, "Algavirales_pred_Pico"] +
                                                      cmip.loc[ocean_mask, "Algavirales_pred_Femto"]) / 2.0
    cmip["Model"] = fname.split("_")[0]
    m_period = re.search(r"\d{4}-\d{4}", fname)
    cmip["Period"] = m_period.group(0) if m_period else "unknown"
    out_cols = [c for c in OUTPUT_COLS if c in cmip.columns]
    cmip_out = cmip.loc[ocean_mask, out_cols].copy()
    out_fname = fname.replace("_data_transformed.csv", "_D_index_predicted_noCO3.csv")
    cmip_out.to_csv(os.path.join(OUTPUT_DIR, out_fname), index=False)
    print(f"Ocean={n_ocean} ({time.time()-t0:.1f}s)")
    all_preds.append(cmip_out)

print("\nSummary")
print("=" * 70)
for df in all_preds:
    mask = df["D_index_pred"].notna()
    if mask.sum() == 0:
        continue
    p = df.loc[mask, "D_index_pred"]
    i = df.loc[mask, "Imitervirales_pred"]
    print(f"  {df['Model'].iloc[0]:20s} {df['Period'].iloc[0]:15s}  "
          f"D mean={p.mean():.4f}  std={p.std():.4f}  "
          f"I mean={i.mean():.4f}")

pd.DataFrame([
    {"Model": df["Model"].iloc[0], "Period": df["Period"].iloc[0],
     "Ocean_Cells": int(df["D_index_pred"].notna().sum()),
     "D_Mean": round(float(df["D_index_pred"].mean()), 4),
     "D_Std": round(float(df["D_index_pred"].std()), 4),
     "Imitervirales_Mean": round(float(df["Imitervirales_pred"].mean()), 4)}
    for df in all_preds
]).to_csv(os.path.join(OUTPUT_DIR, "prediction_summary.csv"), index=False)

print(f"\n{'='*70}\nFINAL REPORT\n{'='*70}")
print(f"  Best single:   {best['name']} (R2(D)={best['r2_D']:.4f})")
print(f"  Best ensemble: {best_ens_name} (R2(D)={best_ens_r2:.4f})")
print(f"  Output:        {OUTPUT_DIR}")
