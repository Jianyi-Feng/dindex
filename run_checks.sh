#!/usr/bin/env bash
# 完整性 + 可运行性自检。用法： bash run_checks.sh
set -uo pipefail
cd "$(dirname "$0")"
fail=0

echo "== [1/4] Python 依赖 =="
python - <<'PY' || { echo "  缺少依赖，尝试安装 requirements.txt"; pip install -r requirements.txt || fail=1; }
import numpy, pandas, scipy, matplotlib, statsmodels, sklearn, openpyxl
PY

echo; echo "== [2/4] 数据布局静态检查 =="
python tools/check_data_layout.py || fail=1

echo; echo "== [3/4] dindex 单元测试 =="
PYTHONPATH=src python -m pytest -q tests/test_dindex.py || fail=1

echo; echo "== [4/4] 示例脚本 =="
python examples/spot_tara_example.py || fail=1

echo
if [ "$fail" = "0" ]; then echo "全部通过 ✓"; else echo "有步骤失败，见上方输出"; fi
exit $fail
