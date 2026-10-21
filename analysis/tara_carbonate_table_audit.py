"""Audit of the Tara Oceans carbonate-chemistry context table used for the 'CO3' predictor.
Downloads PANGAEA registry tables (CARB, HPLC, NUT; release 20170515) and tests, column by column, whether carbonate-system
columns reproduce values from other context tables, and whether training predictors equal registry medians."""
from pathlib import Path
import io, zipfile, urllib.request
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
CACHE = HERE / "tara_registry_cache"; CACHE.mkdir(exist_ok=True)
BASE = "https://store.pangaea.de/Projects/TARA-OCEANS/Samples_Registry/TARA_SAMPLES_CONTEXT_ENV-DEPTH-{}_20170515.zip"


def load(tag):
    f = CACHE / f"{tag}.xlsx"
    if not f.exists():
        z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(BASE.format(tag), timeout=300).read()))
        f.write_bytes(z.read(z.namelist()[0]))
    raw = pd.read_excel(f, header=None)
    pr = [i for i in range(40) if str(raw.iloc[i, 0]) == "PARAMETER"][0]
    st = [i for i in range(40) if str(raw.iloc[i, 0]).startswith("Row#")][0] + 1
    par, com = raw.iloc[pr].astype(str).tolist(), raw.iloc[st - 2].astype(str).tolist()
    data = raw.iloc[st:].reset_index(drop=True); data.columns = [f"{tag}{c}" for c in data.columns]
    return par, com, data


cpar, ccom, cd = load("CARB"); hpar, hcom, hd = load("HPLC"); npar, ncom, nd = load("NUT")
m = cd.merge(hd, left_on="CARB1", right_on="HPLC1").merge(nd, left_on="CARB1", right_on="NUT1")
rows = []
for c in range(19, len(cpar)):
    x = pd.to_numeric(m[f"CARB{c}"], errors="coerce")
    for tag, par, com in (("HPLC", hpar, hcom), ("NUT", npar, ncom)):
        for k in range(14, len(par)):
            y = pd.to_numeric(m[f"{tag}{k}"], errors="coerce"); ok = x.notna() & y.notna()
            if ok.sum() < 100:
                continue
            share = float(np.mean(np.isclose(x[ok], y[ok], atol=1e-9)))
            if share > 0.99:
                rows.append(dict(carb_column=c, carb_parameter=cpar[c], carb_statistic=ccom[c][:40], identical_to_table=tag,
                                 identical_to_parameter=par[k], identical_to_statistic=com[k][:40], n_rows=int(ok.sum()), share_identical=share))
aud = pd.DataFrame(rows); aud.to_csv(HERE / "tara_carbonate_table_audit.csv", index=False)
print(aud[["carb_parameter", "carb_statistic", "identical_to_parameter", "identical_to_statistic", "share_identical"]].to_string())

t = pd.read_csv(HERE.parents[1] / "data/Tara_observed_Dindex.csv")
tm = t.merge(m, left_on="sample_id", right_on="CARB1")
med = lambda par, tag, name: [k for k in range(14, len(par)) if par[k] == name and "median" in (ccom if tag == "CARB" else hcom if tag == "HPLC" else ncom)[k]][0]
checks = {"CO3 = CARB 'Carbonate ion' median": ("CO3", f"CARB{med(cpar, 'CARB', 'Carbonate ion')}"),
          "CO3 = HPLC 19'-butanoyloxyfucoxanthin median": ("CO3", f"HPLC{med(hpar, 'HPLC', '19-Butanoyloxyfucoxanthin')}"),
          "PO4 = NUT phosphate median": ("PO4", f"NUT{med(npar, 'NUT', 'Phosphate')}"),
          "Si = NUT silicate median": ("Si", f"NUT{med(npar, 'NUT', 'Silicate')}")}
for lab, (a, b) in checks.items():
    y = pd.to_numeric(tm[b], errors="coerce"); ok = y.notna()
    print(f"{lab}: {np.mean(np.isclose(tm.loc[ok, a], y[ok], atol=1e-9)):.3f} of {int(ok.sum())} matched samples identical")
