#!/usr/bin/env bash
# Integrity and runnability self-check. Usage: bash run_checks.sh
set -uo pipefail
cd "$(dirname "$0")"
fail=0

echo "== [1/4] Python dependencies =="
python - <<'PY' || { echo "  missing dependencies, trying requirements.txt"; pip install -r requirements.txt || fail=1; }
import numpy, pandas, scipy, matplotlib, statsmodels, sklearn, openpyxl
PY

echo; echo "== [2/4] Static check of the data layout =="
python tools/check_data_layout.py || fail=1

echo; echo "== [3/4] dindex unit tests =="
PYTHONPATH=src python -m pytest -q tests/test_dindex.py || fail=1

echo; echo "== [4/4] Worked example =="
python examples/spot_tara_example.py || fail=1

echo
if [ "$fail" = "0" ]; then echo "All checks passed"; else echo "Some steps failed; see the output above"; fi
exit $fail
