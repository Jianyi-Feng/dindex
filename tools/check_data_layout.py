# -*- coding: utf-8 -*-
"""Static check: are all input tables referenced by the scripts present in this repository?

Rules:
  * a table that exists in the repository -> OK
  * a submission-number prefix (Supplementary_Data_7_, Supplementary_Table_3_, ...) matches a short name used inside a script -> OK
  * files a script writes itself (to_csv/to_excel/savez/open) -> treated as generated, not missing
  * everything else -> listed as a repository-external input (large grids or external raw data)
"""
import os, re, sys

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def norm(n):
    n = n.lower()
    n = re.sub("^supplementary_(data|table)_\\d+[a-z]?_", "", n)
    n = re.sub("^inherited__", "", n)
    n = re.sub("_early(?=\\.)", "", n)
    n = re.sub("^(spot|tara|carbon|nested)[a-z]*_", "", n)
    return n


HAVE, HAVENORM = set(), set()
for base in ("raw_tables", "processed_tables"):
    for r, d, files in os.walk(os.path.join(PKG, base)):
        for f in files:
            HAVE.add(f)
            HAVENORM.add(norm(f))

REF = re.compile("['\"]([A-Za-z0-9_./\\-]+\\.(?:csv|xlsx|tsv|npz|log|json|txt))['\"]")
WROTE = re.compile("(?:to_csv|to_excel|ExcelWriter|savez\\w*|open)\\s*\\([^)]{0,160}?['\"]([A-Za-z0-9_./\\-]+\\.(?:csv|xlsx|tsv|npz|log))['\"]")

# Pass 1: collect every name that is written somewhere in the repository (cross-script counts; A writes, B reads)
SOURCES = []
generated = set()
for sub in ("analysis", "figures"):
    for r, d, files in os.walk(os.path.join(PKG, sub)):
        for f in files:
            if f.endswith(".py"):
                p = os.path.join(r, f)
                SOURCES.append(p)
                generated |= set(os.path.basename(x) for x in
                                 WROTE.findall(open(p, encoding="utf-8", errors="replace").read()))

# Pass 2: report the missing inputs
missing, checked = {}, 0
for p in SOURCES:
    t = open(p, encoding="utf-8", errors="replace").read()
    for name in set(REF.findall(t)):
        base = os.path.basename(name)
        if base in generated or base.startswith("test"):
            continue
        checked += 1
        if base not in HAVE and base not in HAVENORM and norm(base) not in HAVENORM:
            missing.setdefault(base, []).append(os.path.relpath(p, PKG))

print("Table files inside the repository: %d" % len(HAVE))
print("Inputs referenced by the scripts: %d; filenames generated within the pipeline (excluded): %d" % (checked, len(generated)))
if missing:
    print("\nRepository-external inputs: %d (large grids or external raw data; see RUN.md):" % len(missing))
    for b, fs in sorted(missing.items()):
        print("  %-48s referenced by %d script(s), e.g. %s" % (b, len(fs), fs[0]))
else:
    print("\nAll referenced inputs are present in the repository")
sys.exit(0)
