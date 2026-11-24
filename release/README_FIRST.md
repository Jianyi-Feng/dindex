# D-index standalone reviewer release 0.1.2

Extract `dindex-0.1.2.tar.gz` for the complete source, base-R functions, examples,
tests, documentation and citation. Install the wheel for Python library/CLI use.
Inside the extracted source directory:

```bash
python -m pip install ".[calibrate,test]"
python -m pytest tests -q
Rscript tests/test_dindex.R
python examples/spot_tara_example.py
```

## What changed from 0.1.1

Documentation only. The Python and R code, the tests and the example data are
identical to 0.1.1; only the version number and three passages of text differ:

- the manuscript's projection model is now described as an ensemble of a ridge
  regression and a random forest (0.1.1 said "14-model ensemble", which the
  manuscript no longer uses);
- the PolB phylotype catalogues are now supplied with the manuscript
  (Supplementary Data 1 and 1b), so they are no longer listed as missing;
- this note no longer refers to earlier manuscript drafts.

## Scope

The README and docs/OWN_DATA_WORKFLOW.md cover calculation from your own abundance
and taxonomy tables, QC, uncertainty, sensitivity, and saved local bounded-D
prediction. Prediction uses a new single-driver calibration fitted to user data;
it is distinct from the manuscript's projection ensemble. It produces point
predictions and range flags; independent validation is required for scientific use.

## Validation records (validation/)

They reflect this release environment (see python_environment.txt), not the
software versions used historically for the manuscript.

- Source and wheel were built from the 0.1.2 source; the 14 Python tests pass when
  run from the extracted source against the installed wheel.
- The worked SPOT/Tara example reproduces the 0.1.1 output to 11 significant digits
  (D matches the published values to 1.3e-10; SPOT slope 0.112 D per degC, n = 157).
- The toy command-line outputs agree with those of 0.1.1 to floating-point rounding.
- The two bundled SPOT/Tara tables are byte-identical to the study's data tables,
  and D in the SPOT table equals Supplementary Data 2 of the manuscript.
- **R tests were not re-run**: R was not installed on the build machine. The R
  source and its test file are byte-identical to 0.1.1, where the tests passed.
  Run `Rscript tests/test_dindex.R` once before publishing.

## Author input before submission

- The manuscript's Code availability statement points to
  https://github.com/Jianyi-Feng/dindex. Make this exact release available there,
  preferably with a persistent archive and DOI. Nothing was uploaded by this local
  packaging, and public availability has not been verified.
- The same statement says the repository also holds the code for the heatwave and
  calibration analyses, cross-validation, CMIP6 projections, dominance fronts and
  figure generation. This toolkit alone does not satisfy that; the study scripts
  and their inputs must be added to the repository.
- Confirm software citation authorship (CITATION.cff lists J. Feng and S. Hou) and
  the terms for redistributing the SPOT/Tara example tables.
