# Run and self-check (RUN)

This repository ships **all input tables** (`raw_tables/`) and the **analysis/figure output
tables** (`processed_tables/`), so a fresh clone needs no extra download before the integrity
self-check.

## One-command self-check

```bat
python tools\check_data_layout.py
```

It performs four steps:

1. checks (and if needed installs) the Python dependencies listed in `requirements.txt`;
2. **static check**: every input table referenced by the scripts in `analysis/` and `figures/`
   must exist in `raw_tables/` or `processed_tables/` (`tools/check_data_layout.py`);
3. runs the `dindex` unit tests (`tests/test_dindex.py`);
4. runs a worked example (`examples/spot_tara_example.py`).

Any failing step prints its reason, and a summary is given at the end.

## Input requirements per script

| Group | Rerunnable as shipped | Notes |
|---|---|---|
| `src/`, `R/`, `tests/`, `examples/` | yes | only uses observation tables in `raw_tables/` |
| `analysis/inherited/dindex_analyses*.py` | yes | inputs = `raw_tables/` |
| `analysis/tara_*.py`, `nested_cv_without_pigment.py`, `lobo_predictions.py` | yes | inputs = `raw_tables/Tara_observed_Dindex.csv` etc. |
| `analysis/pre2014_pair_derivation.py`, `pnps_and_network.py`, `assemble_supplementary_data.py` | yes | inputs = `raw_tables/` + `processed_tables/` |
| `analysis/cmip6_extract.py` | needs external input | raw CMIP6 model output (not distributed; see file header) |
| `analysis/cmip6_analyses.py`, `cmip6_gridded_analyses.py`, `cmip6_aoa_fronts.py`, `climate_velocity_comparison.py` | needs the `*.npz` grids | grids are too large and live in the data deposit (see README) |
| `figures/*.py` | partly needs the grids | figures that need an `.npz` skip that panel and say so in the log |
| `checks/*.py` | not standalone | pre-submission package checks against the full submission tree; kept for provenance |

## Path conventions

Every script locates itself with `Path(__file__).resolve()` and then resolves paths from the
repository root:

```
PKG        = <repository root>
raw_tables = PKG / "raw_tables"            # input tables
proc       = PKG / "processed_tables"      # analysis and figure output tables
outputs    = PKG / "outputs"               # figures and source data written by the scripts (git-ignored)
```

To keep data elsewhere, edit those constants at the top of the script; do not commit absolute paths.

## Environment

Python >= 3.9, plus R (only for `R/dindex.R` and `tests/test_dindex.R`). The map panels need
`cartopy`; the network figure needs `networkx`. Random seeds are fixed inside the scripts
(random forest `random_state=1`, inner random search `seed=42`, bootstrap `seed=0`).
