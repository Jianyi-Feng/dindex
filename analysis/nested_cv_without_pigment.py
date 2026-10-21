"""Nested CV of the final ensemble with and without the Tara 'CO3' predictor, which the PANGAEA carbonate context table
(TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB_20170515) populates with the HPLC median of 19'-butanoyloxyfucoxanthin.
Outer: GroupKFold(5) by station (equivalent to the supplied MD5 env-hash groups); inner: as in Predict_D-index.py."""
import runpy, numpy as np, pandas as pd
from pathlib import Path
from sklearn.model_selection import GroupKFold, KFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.base import clone
from sklearn.metrics import r2_score, mean_squared_error
HERE = Path(__file__).resolve().parent
src = (HERE / "carbonate_sensitivity.py").read_text().split("def fit_ensemble")[0]
ns = {"__file__": str(HERE / "carbonate_sensitivity.py")}; exec(src, ns)
t, I_vals, A_vals, femto_tr, models, log10_safe, LIN, LOG, RS, NJ = (ns[k] for k in ["t", "I_vals", "A_vals", "femto_tr", "models", "log10_safe", "LIN", "LOG", "RS", "NJ"])
yD = (I_vals - A_vals) / (I_vals + A_vals); groups = t.station_id.to_numpy()
out = []
for label, logs in [("8 predictors (incl. mislabelled 'CO3' = 19'-butanoyloxyfucoxanthin)", LOG), ("7 predictors (without 'CO3')", [f for f in LOG if f != "CO3"])]:
    Xn = np.column_stack([t[f].to_numpy(float) for f in LIN] + [log10_safe(t[f].to_numpy(float)) for f in logs])
    X = np.hstack([Xn, femto_tr[:, None]]); nn = Xn.shape[1]
    D_models = []
    for name, est, ps in models():
        oofI = np.full(len(yD), np.nan); oofA = np.full(len(yD), np.nan)
        for fi, (tr, te) in enumerate(GroupKFold(5).split(X, yD, groups)):
            for y, oof in ((I_vals, oofI), (A_vals, oofA)):
                pipe = Pipeline([("prep", ColumnTransformer([("num", StandardScaler(), list(range(nn))), ("cat", "passthrough", [nn])])), ("model", clone(est))])
                if ps:
                    pipe = RandomizedSearchCV(pipe, {"model__" + k: v for k, v in ps.items()}, n_iter=40, cv=KFold(3, shuffle=True, random_state=RS + fi),
                                              scoring="r2", n_jobs=NJ, random_state=RS, refit=True).fit(X[tr], y[tr]).best_estimator_
                else:
                    pipe.fit(X[tr], y[tr])
                oof[te] = np.ravel(pipe.predict(X[te]))
        Ip, Ap = np.clip(oofI, 0, 1), np.clip(oofA, 0, 1); den = Ip + Ap
        D_models.append(np.divide(Ip - Ap, den, out=np.zeros_like(Ip), where=den > 1e-12))
        print(label[:12], name, flush=True)
    Dens = np.clip(np.mean(D_models, 0), -1, 1)
    fold_r2 = [r2_score(yD[te], Dens[te]) for _, te in GroupKFold(5).split(X, yD, groups)]
    out.append(dict(model=label, pooled_R2=r2_score(yD, Dens), mean_fold_R2=float(np.mean(fold_r2)), RMSE=mean_squared_error(yD, Dens) ** 0.5, n=len(yD)))
    pd.DataFrame({"sample_id": t.sample_id, "station_id": groups, "observed_D": yD, "oof_D": Dens}).to_csv(HERE / f"nestedcv_oof_{len(logs)+2}pred.csv", index=False)
res = pd.DataFrame(out); res.to_csv(HERE / "nestedcv_with_without_pigment.csv", index=False)
print(res.round(4).to_string())
