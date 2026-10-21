"""Full-resolution rerun of the projection ensemble without the mislabelled Tara 'CO3' predictor ( design).
Fits the 14-model final ensemble (I and A separately; identical transforms, search and clipping) on the 196 Tara metagenomes with
7 predictors, predicts every archived CMIP6 grid cell (4 ESMs x 10 decades), and writes cmip6_grid.npz in the same layout as
archived__cmip6_grid.npz (environmental columns copied; D/I/A columns replaced by  predictions)."""
import time, numpy as np, pandas as pd
from pathlib import Path
HERE = Path(__file__).resolve().parent
src = (HERE / "carbonate_sensitivity.py").read_text().split("z = np.load")[0]
ns = {"__file__": str(HERE / "carbonate_sensitivity.py")}; exec(src, ns)
LIN, LOG, log10_safe = ns["LIN"], [f for f in ns["LOG"] if f != "CO3"], ns["log10_safe"]
t0 = time.time()
ens, _ = ns["fit_ensemble"](LOG)
print(f"ensemble fitted {time.time() - t0:.0f}s", flush=True)
z = np.load(HERE / "archived__cmip6_grid.npz"); cols = list(z["columns"]); C = {k: i for i, k in enumerate(cols)}
out = {}


def predict_all(Xn):
    res = {}
    for fr_name, fr in (("Pico", 0.0), ("Femto", 1.0)):
        X = np.hstack([Xn, np.full((len(Xn), 1), fr)])
        D = np.zeros(len(Xn)); I = np.zeros(len(Xn)); A = np.zeros(len(Xn))
        for _, mi, ma in ens:
            Ip = np.clip(np.ravel(mi.predict(X)), 0, 1); Ap = np.clip(np.ravel(ma.predict(X)), 0, 1); den = Ip + Ap
            D += np.divide(Ip - Ap, den, out=np.zeros_like(Ip), where=den > 1e-12) / len(ens)
            I += Ip / len(ens); A += Ap / len(ens)
        res[fr_name] = (np.clip(D, -1, 1), I, A)
    return res


for key in [k for k in z.files if k != "columns"]:
    a = z[key].copy()
    Xn = np.column_stack([a[:, C[f]].astype(float) for f in LIN] + [log10_safe(a[:, C[f]].astype(float)) for f in LOG])
    ok = np.isfinite(Xn).all(1)
    r = predict_all(Xn[ok])
    for fr in ("Pico", "Femto"):
        a[ok, C[f"D_index_pred_{fr}"]], a[ok, C[f"Imitervirales_pred_{fr}"]], a[ok, C[f"Algavirales_pred_{fr}"]] = r[fr]
    a[ok, C["D_index_pred"]] = (r["Pico"][0] + r["Femto"][0]) / 2
    if "Imitervirales_pred" in C:
        a[ok, C["Imitervirales_pred"]] = (r["Pico"][1] + r["Femto"][1]) / 2
        a[ok, C["Algavirales_pred"]] = (r["Pico"][2] + r["Femto"][2]) / 2
    a[~ok, C["D_index_pred"]] = np.nan
    out[key] = a
    print(key, int(ok.sum()), f"{time.time() - t0:.0f}s", flush=True)
np.savez_compressed(HERE / "cmip6_grid.npz", **out, columns=np.array(cols))
print("DONE", flush=True)
