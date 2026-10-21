# -*- coding: utf-8 -*-
"""Dependency and font gate for the figure scripts (same checks as the original render gate).

    python figures/check_render_env.py [--strict]

--strict makes a missing Arial face a hard failure (the published figures require Arial).
"""
import importlib.util, sys

NEED = ["numpy", "pandas", "scipy", "matplotlib", "cartopy", "statsmodels", "openpyxl", "PIL"]
missing = [m for m in NEED if importlib.util.find_spec(m) is None]
if missing:
    print("MISSING PACKAGES:", ", ".join(missing))
    print("install with:  pip install -r requirements.txt")
    sys.exit(1)

import matplotlib
import matplotlib.font_manager as fm

print("matplotlib:", matplotlib.__version__)
if not matplotlib.__version__.startswith("3.11"):
    print("note: the manuscript PDFs were rendered with matplotlib 3.11.2; "
          "other versions may shift tick padding and layout spacing")

faces = sorted({f.name for f in fm.fontManager.ttflist if f.name.lower().startswith("arial")})
resolved = fm.findfont(fm.FontProperties(family=["Arial", "sans-serif"]), fallback_to_default=True)
family = fm.FontProperties(fname=resolved).get_name()
print("Arial faces visible :", faces or "NONE")
print("resolved for 'Arial': %s  (family name: %s)" % (resolved, family))

ok = family.lower().startswith("arial") or "arial" in str(resolved).lower()
if ok:
    print("font check: PASS")
    sys.exit(0)

print("font check: FAIL -- matplotlib would substitute another family, so figures would not be Arial.")
print("  sudo apt-get install ttf-mscorefonts-installer")
print("  # or copy arial*.ttf into ~/.fonts && fc-cache -f -v ~/.fonts")
print("  rm -rf \"$(python -c 'import matplotlib;print(matplotlib.get_cachedir())')\"")
print("See figures/README.md for details.")
sys.exit(1 if "--strict" in sys.argv else 0)
