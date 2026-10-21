# 运行与自检（RUN）

本仓库自带**全部输入表格**（`raw_tables/`）与**分析/绘图产出的表格**（`processed_tables/`），
所以在任何一台机器上 clone 之后不需要额外下载数据就能做完整性自检。

## 一键自检

```bash
bash run_checks.sh
```

它依次做四件事：

1. 检查/安装 Python 依赖（`requirements.txt`）
2. **静态检查**：扫描 `analysis/`、`figures/` 里脚本引用的每个输入表，确认它们在
   `raw_tables/` 或 `processed_tables/` 里确实存在（`tools/check_data_layout.py`）
3. 跑 `dindex` 包的单元测试（`tests/test_dindex.py`）
4. 跑一个可运行的示例（`examples/spot_tara_example.py`）

每一步的失败原因都会打印出来，最后给出总结。

## 各脚本的输入要求

| 类别 | 能否直接重跑 | 说明 |
|---|---|---|
| `src/`、`R/`、`tests/`、`examples/` | ✅ 可以 | 只依赖 `raw_tables/` 里的观测表 |
| `analysis/inherited/dindex_analyses*.py` | ✅ 可以 | 输入 = `raw_tables/` |
| `analysis/tara_*.py`、`nested_cv_without_pigment.py`、`lobo_predictions.py` | ✅ 可以 | 输入 = `raw_tables/Tara_observed_Dindex.csv` 等 |
| `analysis/pre2014_pair_derivation.py`、`pnps_and_network.py`、`assemble_supplementary_data.py` | ✅ 可以 | 输入 = `raw_tables/` + `processed_tables/`；`assemble_supplementary_data.py` 另需 §1 的 `cmip6_maps_2deg.npz`（Data 9b） |
| `analysis/cmip6_extract.py` | ⚠️ 需要外部输入 | 需要 CMIP6 原始模型输出（未随包分发，见文件头注释） |
| `analysis/cmip6_analyses.py`、`cmip6_gridded_analyses.py`、`cmip6_aoa_fronts.py`、`climate_velocity_comparison.py` | ⚠️ 需要 `*.npz` 网格 | 网格文件体积过大，放在数据存档（见 README 的 Data 一节），下载后放进 `processed_tables/cmip6_gridded_fields/` |
| `figures/*.py` | ⚠️ 部分需要网格 | 依赖上述 npz 的图会跳过并在日志里说明 |
| `checks/*.py` | ❌ 不作为独立脚本 | 它们是投稿前的包完整性检查，针对当时的整棵投稿目录；保留供追溯 |

## 路径约定

每个脚本用 `Path(__file__).resolve()` 定位自己，再按仓库根解析：

```python
PKG        = <仓库根>
raw_tables = PKG / "raw_tables"            # 原始表格
proc       = PKG / "processed_tables"      # 分析与绘图产出的表格
outputs    = PKG / "outputs"               # 脚本写出的图与源数据（已被 .gitignore 忽略）
```

如需把数据放在别处，改脚本顶部的那几个路径常量即可；不要改成绝对路径再提交。

## 环境

Python >= 3.9，R（仅 `R/dindex.R` 与 `tests/test_dindex.R` 需要）。
地图面板需要 `cartopy`；网络图需要 `networkx`。随机种子已固定在各脚本内
（森林 `random_state=1`、内层随机搜索 `seed=42`、bootstrap `seed=0`）。

## Reproducing the pipeline (measured end to end on Windows, Python 3.11.4)

All numbers and figures below were produced by actually running the scripts in a clean
copy of this repository; logs were captured for every step.

### 0. Dependencies (required)

```
python -m pip install pandas openpyxl statsmodels scikit-learn cartopy shapely pyproj ^
      pyshp networkx pillow matplotlib pytest
```
Tested combination: pandas 2.2.3 / numpy 1.26.4 / scipy 1.14.0 / matplotlib 3.9.2 /
statsmodels 0.15.0 / scikit-learn 1.9.1 / cartopy 0.26.0 / openpyxl 3.1.5.
The map panels (Fig. 5, Extended Data Fig. 5) need cartopy Natural Earth shapefiles; they
are downloaded on first use, so an offline machine needs a pre-populated cartopy data dir.

### 1. Gridded CMIP6 arrays (from the figshare deposit)

Download from figshare (DOI 10.6084/m9.figshare.34419060, `07_CMIP6_gridded_fields/`) into
`processed_tables/cmip6_gridded_fields/`:

- required by the figure scripts: `cmip6_maps_2deg.npz`, `cmip6_aoa_2deg.npz`, `cmip6_grid.npz`
- only for the consistency check in `analysis/cmip6_analyses.py`: the three `archived__cmip6_*.npz`

`tools/paths.py` (`require_gridded()`) reports exactly which files are missing.

### 2. Path conventions

All input/output locations live in `tools/paths.py` (`RAW`, `SUPP_DATA`, `SUPP_TABLES`,
`FIG_SOURCE`, `MODEL_STATS`, `PHYLO`, `CMIP6`, `GENOMES`, `OUT`), together with `resolve()`
(maps the historical `data/…`, `analysis/…`, `archive/supplementary/…` names onto the
current layout) and `outdir()`. Scripts no longer hard-code directories. Outputs default
to `<repo>/outputs/` and the directory is in `.gitignore`; set `DINDEX_OUT` to redirect.
`raw_tables/` and `processed_tables/` are never written by the scripts.

### 3. Core chain (measured: 4.5 s / 56 s / 52 s / 15 s)

```
python analysis\assemble_supplementary_data.py            -> outputs/supplementary_data/*.csv
python analysis\pnps_and_network.py                       -> outputs/tables/pnps_order_tests.csv,
                                                              outputs/supplementary_tables/Supplementary_Table_14_pNpS.xlsx,
                                                              outputs/figures/SupplementaryFig4_pNpS.*, SupplementaryFig5_*
python figures\make_figures.py                            -> outputs/Fig1..Fig5 + ExtendedDataFig1..5 (.pdf/.png/.svg)
                                                              + outputs/source_data/*_source_data.xlsx
python figures\make_extended_data_fig4.py                 -> outputs/ExtendedDataFig5_projection_robustness.*
python figures\make_supplementary_figs4_5_pnps_network.py -> optional; same numbers, adds SVG
```

Two independent runs gave byte-identical PNGs (12/12); PDF/SVG differ only in the embedded
creation timestamp, so the statistics and plotting steps are deterministic.

### 4. Known content differences from the deposited tables (not runtime errors)

`analysis/assemble_supplementary_data.py` reproduces the deposited **Supplementary Data 2-10 and
9b cell for cell** when its output is diffed against the read-only deposit
(`figshare/dindex_data/01_Supplementary_Data/`): same rows, same columns, same values, same
wording. Nothing below changes a statistic, threshold, model or selection rule -- the differences
that used to remain were label spelling and number rendering, and both are reconciled in the
script's display layer.

| Table | Differing cells before | now | What is reconciled |
|---|---|---|---|
| Data 2 | 51 | 0 | short campaign labels, and `heatwave_phase` re-derived from the declared event window |
| Data 9 | 8 | 0 | duplicated word in the projection label collapsed in both the trajectory and the rate part |
| Data 10 | 579 | 0 | deposited display precision, deposited statistic wording, integer `n` |
| Data 8 | 3085 | 0 | copied back as text instead of through a float round-trip |
| Data 3-7, 9b | 0 | 0 | unchanged |

Why the inputs needed it (all four are input-table defects or properties of the deposited file
format, not runtime errors):

* **Data 2 -- `campaign` and `heatwave_phase`.** The deposited table drops the duplicated
  `extraction_method` column and uses the short campaign names, which the script maps. The phase
  is a genuine input defect: `raw_tables/SPOT_heatwave_Dindex_stats.csv` declares the event window
  as `2014-01 .. 2015-12`, and the deposited Data 2/3 label the samples inside it
  `early_heatwave` and every later sample `post`, whereas the per-sample column of
  `raw_tables/SPOT_sample_Dindex.csv` marks 2014-01 ... 2016-12 as `heatwave` (12 samples of 2016
  in the wrong phase, short label). The author's own full-precision copy of the same table
  carries exactly the deposited phases (pre 89 / early_heatwave 39 / post 39). The script reads
  the window from the shipped statistics table and re-labels those samples; `pre` is untouched.
* **Data 9 -- `projection`.** `archived ` was prefixed to a label that already began with
  `archived`, in both
  `processed_tables/cmip6_gridded_fields/cmip6_front_trajectories_by_projection.csv` and
  `processed_tables/model_and_statistics/cmip6_front_rates_by_projection.csv`.
* **Data 10 -- number rendering and wording.** The deposited file carries a UTF-8 BOM and CRLF
  line endings, and each of its numeric cells is rounded the way a spreadsheet renders "General"
  numbers (as many decimals as fit 11 characters, trailing zeros dropped, scientific with two
  decimals below 1e-4). The full-precision values computed here agree with the author's
  full-precision copy of the same table in every numeric cell except one Mann-Whitney P in its
  last floating-point digit, so rounding to the deposited precision is a display step only. The
  deposited table carries no internal version tags: its statistic labels are the computed labels,
  and the pN/pS sensitivity note is rendered from the `role` column of
  `processed_tables/model_and_statistics/pnps_order_tests.csv`. The script reproduces the
  deposited rendering and wording; it does not edit the source tables.
* **Data 8 -- encoding.** `Supplementary_Data_8_GVMAG_microdiversity.csv` in
  `processed_tables/supplementary_data/` is the deposited table itself (see `tools/paths.py`), so
  it is copied back as text: a float round-trip rewrote the two values the deposit prints as
  `8.09E-05` and `8.68E-05`, and its exact zeros as `0.0`.

Reproducing Data 9b needs the `cmip6_maps_2deg.npz` array from the deposit (section 1); it is not
shipped in the repository.

Two differences remain outside the cell values, and both are file format rather than content:

* `README_Supplementary_Data.csv` differs: the deposit lists Supplementary Data 1 first (1095
  rows, per-phylotype PolB catalogue) and words the Data 4 and Data 10 descriptions differently.
  Data 1 cannot be assembled from this repository -- no per-phylotype catalogue is shipped -- so
  this index keeps the "author action required" line at the end.
* The assembled files are written by pandas (UTF-8 without BOM, CRLF on Windows) while the
  deposited tables carry a UTF-8 BOM and, except Data 10, LF line endings; the 9b archive differs
  only in its gzip framing. `cmp` therefore still reports a difference for cells that are
  identical.

### 5. Scripts that cannot run from this repository alone

`analysis/Predict_D-index_noCO3.py` (needs three author-local input files),
`analysis/cmip6_extract.py` (needs the raw CMIP6 model output) and
`analysis/cmip6_analyses.py` (needs `archived__cmip6_front_trajectories.csv`).
`python tools/check_data_layout.py` lists these five repository-external inputs explicitly.

Three scripts are kept under `analysis/` but cannot run from the deposited data alone (`Predict_D-index_noCO3.py`, `cmip6_extract.py`, `cmip6_analyses.py`); they are listed by `python tools/check_data_layout.py`. `Predict_D-index_noCO3.py` has been moved to `docs/attic/` together with a note on the three author-local inputs it requires.
