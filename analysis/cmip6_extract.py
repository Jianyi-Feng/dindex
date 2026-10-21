"""Stream the 40 gridded CMIP6 prediction files into one compact float32 store (archived__cmip6_grid.npz)
and record provenance checks on how D_index_pred was assembled."""
from pathlib import Path
import numpy as np, pandas as pd
SRC = Path(__file__).resolve().parents[2] / "external/marker_recovery_source/4_CMIP6_model_predictions_of_D-index_during_2000-2099"
COLS = ["longitude", "latitude", "Temperature", "ChlorophyllA", "Salinity", "CO3", "PO4", "NO3", "Si", "Iron",
        "D_index_pred_Pico", "D_index_pred_Femto", "D_index_pred", "Imitervirales_pred_Pico", "Imitervirales_pred_Femto",
        "Algavirales_pred_Pico", "Algavirales_pred_Femto", "Imitervirales_pred", "Algavirales_pred"]
store, checks = {}, []
for f in sorted(SRC.glob("*.csv")):
    esm, dec = f.name.split("_merged_")[0], f.name.split("_merged_")[1][:9]
    d = pd.read_csv(f, usecols=COLS, dtype={c: np.float32 for c in COLS})
    Dp, Df, D = d.D_index_pred_Pico, d.D_index_pred_Femto, d.D_index_pred
    Dfrac = lambda I, A: (I - A) / (I + A)
    checks.append(dict(ESM=esm, decade=dec, n_cells=len(d),
        max_abs_D_minus_mean_DpDf=float(np.abs(D - (Dp + Df) / 2).max()),
        max_abs_Dpico_minus_fromIA=float(np.abs(Dp - Dfrac(d.Imitervirales_pred_Pico, d.Algavirales_pred_Pico)).max()),
        max_abs_Dfemto_minus_fromIA=float(np.abs(Df - Dfrac(d.Imitervirales_pred_Femto, d.Algavirales_pred_Femto)).max()),
        mean_abs_D_minus_from_meanIA=float(np.abs(D - Dfrac(d.Imitervirales_pred, d.Algavirales_pred)).mean()),
        n_negative_chl=int((d.ChlorophyllA <= 0).sum()), n_negative_NO3=int((d.NO3 <= 0).sum())))
    store[f"{esm}|{dec}"] = d[COLS[:-2]].to_numpy(np.float32)
    print(esm, dec, len(d), flush=True)
np.savez_compressed(Path(__file__).parent / "archived__cmip6_grid.npz", **{k.replace("|", "__"): v for k, v in store.items()}, columns=np.array(COLS[:-2]))
pd.DataFrame(checks).to_csv(Path(__file__).parent / "cmip6_provenance_checks.csv", index=False)
print(pd.DataFrame(checks).describe().T.to_string())
