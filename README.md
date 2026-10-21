# dindex — a two-group log-ratio index for microbial community tables

`dindex` computes the **D-index**, a single number that summarises the balance between two groups
of taxa in a sample:

$$D = \frac{I - A}{I + A} = \tanh\!\left(\tfrac{1}{2}\ln\frac{I}{A}\right)$$

- **I** is the summed abundance of *Imitervirales* (the "plus" group in the giant-virus study).
- **A** is the summed abundance of *Algavirales* (the "minus" group).

D runs from **−1** (only A detected) through **0** (equal abundance) to **+1** (only I detected).
It depends only on the ratio I/A, so it is unchanged by sequencing depth, by normalisation to
relative abundance, or by any factor applied equally to both groups. It is undefined when I + A = 0.

Nothing in the calculation is specific to viruses. It applies to any two groups whose balance you
expect to track an environmental gradient — ecotypes or clades, warm- vs cold-associated lineages,
copiotrophs vs oligotrophs, or two functional gene families. The input can be read counts, coverage,
TPM, relative abundance or cell counts.

This repository also contains the **analysis code for the accompanying manuscript**:

> Feng, J., Hou, S. *et al.* Recurrent warming-associated shifts in marine giant-virus dominance
> across scales. (in submission).

---

## Repository layout

| Path | Contents |
|---|---|
| `src/dindex/` | Python package: `core.py` (D-index), `calibrate.py` (sensitivity, zero point, seasonal baseline, pair screening), `cli.py` (command line) |
| `R/dindex.R` | Base-R implementation of the same functions |
| `tests/` | `pytest` suite (Python) and `test_dindex.R` (R) |
| `docs/METHODS_TEMPLATE.md` | Adaptable Methods paragraph for a new application |
| `docs/REPORTING_CHECKLIST.md` | Minimum reporting standard |
| `examples/` | Worked SPOT + Tara example and a small toy count table for the CLI |
| `raw_tables/` | **Raw input tables** used by the manuscript analyses (sample-level D values, Tara predictors, CMIP6 latitude–year field, network edge/node lists, rhodopsin metadata) |
| `analysis/` | Analysis scripts: SPOT time series, heatwave inference, Tara calibration and validation, CMIP6 projections and dominance fronts, genome tables, pN/pS and the association network |
| `processed_tables/` | **Analysis and figure output tables**: Supplementary Data, Supplementary Tables, figure source data, model and statistics tables, genome microdiversity, phylogenies, gridded CMIP6 fields |
| `figures/` | Figure-generation scripts (main Figs. 1–5, Extended Data Figs. 1–5) |
| `checks/` | Pre-submission consistency, figure and package checks (kept for provenance; they expect the full submission tree) |
| `outputs/` | Written by the analysis and figure scripts (ignored by git) |
| `RUN.md`, `run_checks.sh` | How to verify and rerun the package on any machine |

## Install

```bash
git clone https://github.com/Jianyi-Feng/dindex.git
cd dindex
pip install -e .                 # from this folder; core needs numpy, pandas, scipy
pip install -e ".[calibrate]"    # adds statsmodels for the calibration helpers
```

For R, run `source("R/dindex.R")`. It uses base R only, with an optional `d_from_phyloseq()` helper.

Requires Python >= 3.9.

## Quick start

### Command line

```bash
# abundance: samples x features (first column = sample IDs); taxonomy: feature ID + one column per rank
dindex compute --abundance counts.tsv --taxonomy taxonomy.tsv --rank genus \
               --plus Prochlorococcus --minus Synechococcus --counts --min-share 0.2 -o D.tsv

dindex calibrate --table D_env.tsv --driver temperature   # sensitivity + temperature at D = 0
dindex baseline  --table D_dates.tsv --date date --baseline-end 2013-12-31 -o residuals.tsv
```

Several labels can be pooled into one side, e.g. `--plus Imitervirales Pimascovirales`.
Use `--features-as rows` when features are rows. Instead of `--counts`, use `--bootstrap 2000`
for coverage or TPM tables.

### Python

```python
import pandas as pd, dindex
from dindex import calibrate

counts = pd.read_csv("counts.tsv", sep="\t", index_col=0)                    # samples x features
tax = pd.read_csv("taxonomy.tsv", sep="\t").set_index("feature")["order"]     # feature -> group
d = dindex.d_from_table(counts, tax, plus="Imitervirales", minus="Algavirales", min_share=0.5)
d[["D", "share", "qc"]]

dindex.d_index(p_abundance, m_abundance)          # vectorised core formula
dindex.d_ci_counts(n_plus_reads, n_minus_reads)   # Beta (Jeffreys) interval for read counts
dindex.d_ci_bootstrap(coverage, tax, "A", "B")    # feature-composition sensitivity interval

calibrate.sensitivity(df, "temperature", covariates=["size_fraction"], cluster="station")
calibrate.zero_point(df, "temperature", group="station")
calibrate.seasonal_baseline(df, "date", baseline_end="2013-12-31")
```

### R

```r
source("R/dindex.R")
d <- d_from_table(counts, taxonomy_vector, plus = "Imitervirales", minus = "Algavirales", min_share = 0.5)
d_ci_counts(d$plus, d$minus)
d_from_phyloseq(ps, rank = "Order", plus = "Imitervirales", minus = "Algavirales")
d_zero_point(df, "temperature", group = "station")
```

The generic command-line and API use the neutral names `plus` (= I) and `minus` (= A).

## What the output columns mean

| Column | Meaning |
|---|---|
| `plus`, `minus`, `total` | Summed abundance of P, of M, and of all features mapped in the taxonomy |
| `D` | The index; `NaN` when P + M = 0 |
| `log_ratio` | ln(P/M); ±inf for zeros unless `--pseudocount` is set. D itself never needs a pseudocount |
| `share` | (P + M)/total: how much of the table the index summarises |
| `qc` | `ok`, `undefined` (P + M = 0), `low_share` or `low_total` |
| `D_lower`, `D_upper` | 95% interval: Beta credible interval for read counts (`--counts`), descriptive feature-composition sensitivity interval otherwise (`--bootstrap`) |
| `n_plus_features`, `n_minus_features` | Number of non-zero features on each side |

## Using D well

1. **Choose the pair before testing it.** A good pair is abundant (high `share`) and responds in
   opposite directions to the gradient of interest. Derive it on one dataset and test it on
   independent data. `calibrate.screen_pairs()` ranks candidate pairs, but screening many pairs
   inflates false positives, so treat its output as hypotheses.
2. **Read D as change against a local baseline.** In the giant-virus study, D's sensitivity to
   temperature was similar across size fractions and between a coastal time series and the open
   ocean. Its zero point was not: D = 0 fell between 10 and 23 °C depending on fraction and system.
3. **Keep the compositional caveats.** D is a relative measure. It does not give absolute abundance,
   and extraction or mapping biases that differ between the two groups change it.
4. **Handle zeros explicitly.** Samples where one group is absent have D = ±1, which is a valid value.
   Samples where both are absent are undefined.
5. **Respect dependence.** Use cluster-robust errors for repeated samples from a station (`cluster=`),
   autocorrelation-robust errors for time series (`hac_lags=`), and year- or event-level inference for
   anomalies from `seasonal_baseline()`.

See `docs/REPORTING_CHECKLIST.md` for a minimum reporting standard and `docs/METHODS_TEMPLATE.md`
for an adaptable Methods paragraph.

## Data

Data availability statements, following the manuscript:

- **Newly generated SPOT metagenomic reads** — NCBI Sequence Read Archive, BioProject
  [`PRJNA1536449`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1536449)
- **Previously released viral metagenomes** — JGI Genome Portal, proposal ID `2799`
- **Previously released 0.2–1.2 µm size-fraction samples** — BioProject
  [`PRJNA814250`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA814250)
- **Nucleocytoviricota and Mirusviricota MAGs** (nucleotide and protein sequences) — figshare,
  [`https://doi.org/10.6084/m9.figshare.34419060`](https://doi.org/10.6084/m9.figshare.34419060)
- **Tara Oceans metagenomes** — ENA `PRJEB1787`, `PRJEB1788`, `PRJEB9740`
- **Companion 16S and 18S rRNA amplicon series** from the same programme are described in our
  previous studies
- **CMIP6 fields** — the four Earth-system models used (CESM2-WACCM, IPSL-CM6A-LR, NorESM2-MM,
  UKESM1-0-LL) are public through the [ESGF archive](https://esgf-node.llnl.gov/)

The observation tables needed to rerun the analyses ship with this repository (`raw_tables/` and
`processed_tables/`), so no extra download is required to get started. The large gridded CMIP6
arrays are archived with the figshare deposit (`07_CMIP6_gridded_fields/`); `RUN.md` lists which
scripts need them.

Run the packaging and runnability check with:

```bat
python tools\check_data_layout.py
```

## Tests

```bash
pip install -e ".[test]"
set PYTHONPATH=src && python -m pytest -q    # Python suite (Windows)
Rscript tests/test_dindex.R             # R checks (needs Rscript)
```

Both suites verify that the functions reproduce the published SPOT and Tara D values exactly.

## Citation

If you use the D-index, please cite the manuscript above and the archived release of this repository.
A `CITATION.cff` file is included; GitHub renders it as "Cite this repository".

## License

MIT — see `LICENSE`.
