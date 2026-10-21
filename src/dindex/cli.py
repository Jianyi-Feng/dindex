"""Command-line interface.

  dindex compute   --abundance table.tsv --taxonomy taxonomy.tsv --rank order --plus Imitervirales --minus Algavirales -o d.tsv
  dindex calibrate --table d_with_env.tsv --driver temperature [--covariate size_fraction] [--cluster station]
  dindex baseline  --table d_with_dates.tsv --date date --baseline-end 2013-12-31 -o residuals.tsv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .core import d_ci_bootstrap, d_ci_counts, d_from_table


def _read(path: str, index_col=None) -> pd.DataFrame:
    sep = "," if Path(path).suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep, index_col=index_col)


def _write(df: pd.DataFrame, path: str | None, index=True):
    if path is None or path == "-":
        df.to_csv(sys.stdout, sep="\t", index=index)
    else:
        df.to_csv(path, sep="," if path.endswith(".csv") else "\t", index=index)


def _split(values):
    out = []
    for v in values:
        out += [x.strip() for x in v.split(",") if x.strip()]
    return out


def cmd_compute(a):
    ab = _read(a.abundance, index_col=0)
    tx = _read(a.taxonomy)
    feat = a.feature_col or tx.columns[0]
    if a.rank not in tx.columns:
        sys.exit(f"--rank '{a.rank}' not a column of the taxonomy table; columns: {list(tx.columns)}")
    taxonomy = tx.set_index(feat)[a.rank].dropna()
    plus, minus = _split(a.plus), _split(a.minus)
    res = d_from_table(ab, taxonomy, plus, minus, features_as=a.features_as, min_share=a.min_share, min_total=a.min_total,
                       pseudocount=a.pseudocount)
    if a.counts:
        _, lo, hi = d_ci_counts(res["plus"], res["minus"])
        res["D_lower"], res["D_upper"] = lo, hi
    elif a.bootstrap:
        ci = d_ci_bootstrap(ab, taxonomy, plus, minus, features_as=a.features_as, n_boot=a.bootstrap, seed=a.seed)
        res["D_lower"], res["D_upper"] = ci["D_lower"], ci["D_upper"]
    res.index.name = "sample"
    _write(res, a.output)
    unmapped = res.attrs.get("unmapped_features", [])
    print(f"[dindex] {len(res)} samples; plus={plus}; minus={minus}; {len(unmapped)} features absent from taxonomy were ignored; "
          f"qc counts: {res['qc'].value_counts().to_dict()}", file=sys.stderr)


def cmd_calibrate(a):
    from .calibrate import sensitivity, zero_point
    df = _read(a.table)
    out = {"sensitivity": sensitivity(df, a.driver, a.d_col, covariates=a.covariate or (), cluster=a.cluster, hac_lags=a.hac_lags)}
    if not a.covariate:
        out["zero_point"] = zero_point(df, a.driver, a.d_col, group=a.cluster, n_boot=a.n_boot)
    else:
        out["zero_point_by_covariate"] = {str(k): zero_point(g, a.driver, a.d_col, group=a.cluster, n_boot=a.n_boot)
                                          for k, g in df.groupby(a.covariate[0])}
    print(json.dumps(out, indent=2))


def cmd_baseline(a):
    from .calibrate import seasonal_baseline
    df = _read(a.table)
    res = seasonal_baseline(df, a.date, a.d_col, baseline_end=a.baseline_end, harmonics=a.harmonics)
    _write(res, a.output, index=False)
    print(f"[dindex] baseline R2={res.attrs['baseline_r2']:.3f} (n={res.attrs['baseline_n']})", file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(prog="dindex", description="Two-group log-ratio index D = (P - M)/(P + M) for microbial community tables.")
    p.add_argument("--version", action="version", version=f"dindex {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("compute", help="compute D per sample from an abundance table and a taxonomy table")
    c.add_argument("--abundance", required=True, help="TSV/CSV; first column = sample IDs, other columns = features (or use --features-as rows)")
    c.add_argument("--features-as", choices=["columns", "rows"], default="columns")
    c.add_argument("--taxonomy", required=True, help="TSV/CSV with a feature ID column and one column per rank")
    c.add_argument("--feature-col", default=None, help="feature ID column in the taxonomy table (default: first column)")
    c.add_argument("--rank", required=True, help="taxonomy column holding the group labels, e.g. order")
    c.add_argument("--plus", nargs="+", required=True, help="group label(s) counted as P (comma- or space-separated)")
    c.add_argument("--minus", nargs="+", required=True, help="group label(s) counted as M")
    c.add_argument("--min-share", type=float, default=None, help="flag samples where (P+M)/total < this value")
    c.add_argument("--min-total", type=float, default=None, help="flag samples where P+M < this value (e.g. reads)")
    c.add_argument("--pseudocount", type=float, default=None, help="pseudocount for the log_ratio column only")
    ci = c.add_mutually_exclusive_group()
    ci.add_argument("--counts", action="store_true", help="inputs are read counts: add Beta (Jeffreys) 95%% interval for D")
    ci.add_argument("--bootstrap", type=int, default=0, help="add descriptive feature-composition 95%% interval with this many replicates")
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("-o", "--output", default=None)
    c.set_defaults(func=cmd_compute)

    k = sub.add_parser("calibrate", help="sensitivity of D to a driver and the driver value at D = 0")
    k.add_argument("--table", required=True); k.add_argument("--driver", required=True); k.add_argument("--d-col", default="D")
    k.add_argument("--covariate", action="append"); k.add_argument("--cluster", default=None); k.add_argument("--hac-lags", type=int, default=None)
    k.add_argument("--n-boot", type=int, default=2000)
    k.set_defaults(func=cmd_calibrate)

    b = sub.add_parser("baseline", help="seasonal harmonic baseline and residuals (e.g. for heatwave analyses)")
    b.add_argument("--table", required=True); b.add_argument("--date", required=True); b.add_argument("--d-col", default="D")
    b.add_argument("--baseline-end", default=None); b.add_argument("--harmonics", type=int, default=1); b.add_argument("-o", "--output", default=None)
    b.set_defaults(func=cmd_baseline)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":  # pragma: no cover
    main()
