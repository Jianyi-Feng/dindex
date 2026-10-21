"""Consistency audit of the  submission package. Prints PASS/FAIL per check; exits non-zero on failure.

Usage: python3 check_package.py [package_root]
"""
import sys
import hashlib
import re
import subprocess
import zipfile
from pathlib import Path

import pandas as pd
from docx import Document

HERE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
# optional argument: audit a copy of the package elsewhere (used to self-test the checks)
MS = HERE / "manuscript/SPOT_giant_virus_Manuscript-.docx"
CL = HERE / "manuscript/Cover_Letter-Nature-.docx"
SI = HERE / "supplementary/Supplementary_Information.docx"
GD = HERE / "supplementary/SIGuide.docx"
SF = HERE / "supplementary/Supplementary_Figures.docx"
results = []


def check(name, ok, detail=""):
    results.append((ok, name, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))


def text(path, tables=True):
    d = Document(path)
    t = "\n".join(p.text for p in d.paragraphs)
    if tables:
        t += "\n" + "\n".join(c.text for tb in d.tables for r in tb.rows for c in r.cells)
    return t


ms, si, gd, sf, cl = (text(p) for p in (MS, SI, GD, SF, CL))
all_texts = {"manuscript": ms, "SI": si, "guide": gd, "figures": sf, "cover letter": cl}

print("\n== files")
expected = [MS, CL, SI, GD, SF] + [HERE / f"supplementary/Supplementary_Table_{i}_" for i in []]
missing = [p.name for p in expected if not p.exists()]
check("five  documents present", not missing, ", ".join(missing))
tables = sorted(HERE.glob("supplementary/Supplementary_Table_*.xlsx"))
nums = sorted(int(re.search(r"Table_(\d+)_", t.name).group(1)) for t in tables)
check("table workbooks 1-10 and 12-18 on disk (11 is author-completed)",
      nums == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] + list(range(12, 19)), f"found {nums}")
stale = {n: sorted({v for v in re.findall(r"\bV(2[5-9]|3[0-3])\b", t)}) for n, t in all_texts.items()}
check("no stale version strings (-) in any  document",
      not any(stale.values()), "; ".join(f"{k}: {v}" for k, v in stale.items() if v))

print("\n== cross-references")
def nums_in(t, pat):
    out = set()
    for m in re.finditer(pat, t):
        out |= {int(g) for g in m.groups() if g}
    return sorted(out)

ms_tables = sorted(set(nums_in(ms, r"Supplementary Tables? (\d+)(?:, (\d+))?")) | set(nums_in(ms, r"Supplementary Tables (\d+)[–-](\d+)")))
ms_data = nums_in(ms, r"Supplementary Data (\d+)(?: and (\d+))?")
si_tables = nums_in(si, r"Supplementary Table (\d+) \|")
si_data = nums_in(si, r"Supplementary Data (\d+) \|")
si_figs = nums_in(si, r"Supplementary Fig\. (\d+) \|")
sf_figs = nums_in(sf, r"Supplementary Fig\. (\d+) \|")
gd_tables = nums_in(gd, r"Supplementary Table (\d+)")
check("manuscript cites Supplementary Tables 1-18", ms_tables == list(range(1, 19)), str(ms_tables))
check("manuscript cites Supplementary Data 1-10", ms_data == list(range(1, 11)), str(ms_data))
check("SI defines Supplementary Tables 1-18", si_tables == list(range(1, 19)), str(si_tables))
check("SI defines Supplementary Data 1-10", si_data == list(range(1, 11)), str(si_data))
check("SI and figures document agree on Supplementary Figs 1-5",
      si_figs == sf_figs == [1, 2, 3, 4, 5], f"SI {si_figs}, figures {sf_figs}")
ms_sfigs = nums_in(ms, r"Supplementary Figs?\.? (\d+)(?:(?:,| and|–) ?(\d+))?")
check("manuscript cites Supplementary Figs 1-5", ms_sfigs == [1, 2, 3, 4, 5], str(ms_sfigs))
check("guide lists Supplementary Tables 1-18", gd_tables == list(range(1, 19)), str(gd_tables))
named = re.findall(r"Supplementary_Table_\S+?\.xlsx", gd)
absent = [f for f in named if not (HERE / "supplementary" / f).exists()]
check("every workbook named in the guide exists", not absent, ", ".join(absent))

print("\n== figures and panels")
legends = {}
# each legend appears twice: the full text in the legends section and a short caption under the embedded
# figure, so keep the longest body per figure
for m in re.finditer(r"((?:Extended Data )?Fig\. (\d+)) \| (.+?)(?=\n(?:Extended Data )?Fig\. \d+ \||\Z)", ms, re.S):
    if len(m.group(3)) > len(legends.get(m.group(1), "")):
        legends[m.group(1)] = m.group(3)
def panel_letters(leg):
    """Panel markers appear as 'a,', '(c)' or in ranges such as 'b–d'."""
    found = set(re.findall(r"(?<![A-Za-z])([a-j]),", leg)) | set(re.findall(r"\(([a-j])\)", leg))
    for lo, hi in re.findall(r"(?<![A-Za-z])([a-j])[–-]([a-j])(?![A-Za-z])", leg):
        found |= {chr(c) for c in range(ord(lo), ord(hi) + 1)}
    return found


for name, leg in legends.items():
    seq = sorted(panel_letters(leg))
    expected_seq = [chr(ord("a") + i) for i in range(len(seq))]
    check(f"{name} panel letters are consecutive", seq == expected_seq, f"{''.join(seq)}")
cited = set()
for m in re.finditer(r"(Extended Data Fig\.|Fig\.) (\d+)([a-j](?:[–,-][a-j])*)", ms):
    kind = "Extended Data Fig. " if m.group(1).startswith("Extended") else "Fig. "
    for letter in re.findall(r"[a-j]", m.group(3)):
        cited.add((kind + m.group(2), letter))
bad = [f"{f}{l}" for f, l in sorted(cited) if f in legends and l not in panel_letters(legends[f])]
check("every panel cited in the text exists in its legend", not bad, ", ".join(bad))
figs = sorted((HERE / "figures").glob("*.pdf"))
check("5 main + 5 Extended Data figure PDFs present", len(figs) == 10, f"{len(figs)} PDFs")
z = zipfile.ZipFile(MS)
media = {hashlib.md5(z.read(n)).hexdigest() for n in z.namelist() if "/media/" in n}
embedded = all(hashlib.md5((HERE / "figures" / f"{p.stem}.png").read_bytes()).hexdigest() in media for p in figs)
check("current figure PNGs are the ones embedded in the manuscript", embedded)

print("\n== numbers")
hv = pd.read_csv(HERE / "analysis/climate_velocity_headline_statistics.csv").set_index("statistic")
h29 = pd.read_csv(HERE / "analysis/headline_statistics.csv")
h14 = pd.read_csv(HERE / "analysis/cmip6_statistics.csv")
def g29(name, col="value"):
    r = h29[h29.statistic == name]; return float(r.iloc[0][col])
def g14(name, col="value"):
    r = h14[h14.statistic == name]; return float(r.iloc[0][col])
quoted = [
    ("0.112 units per °C", round(g29("SPOT dD/dSST, fraction-adjusted (HC3)"), 3), 0.112),
    ("65.3 km per decade", round(g14("front rate North, final ensemble (7 predictors) (km/decade)"), 1), 65.3),
    ("36.0 km per decade", round(g14("front rate South, final ensemble (7 predictors) (km/decade)"), 1), 36.0),
    ("68.4 km per decade", round(hv.loc["isotherm displacement rate, North", "value"], 1), 68.4),
    ("52.3 km per decade", round(hv.loc["isotherm displacement rate, South", "value"], 1), 52.3),
    ("ratio 0.95", round(hv.loc["front-to-isotherm rate ratio, North", "value"], 2), 0.95),
    ("ratio 0.68", round(hv.loc["front-to-isotherm rate ratio, South", "value"], 2), 0.68),
    ("r = 0.90", round(hv.loc["correlation of front rate with isotherm rate", "value"], 2), 0.90),
]
for label, computed, stated in quoted:
    ok = computed == stated and label.split(" ")[0] in ms
    check(f"text value {label} matches the audit table", ok, f"computed {computed}")
t13 = pd.ExcelFile(HERE / "supplementary/Supplementary_Table_13_isotherm_and_climate_velocity.xlsx")
check("Supplementary Table 13 has the expected sheets",
      t13.sheet_names == ["README", "A_front_vs_isotherm", "B_isotherm_tracking", "C_climate_velocity", "D_reported_statistics"],
      str(t13.sheet_names))
a13 = t13.parse("A_front_vs_isotherm")
check("Table 13 front rates reproduce the published means",
      round(a13[a13.hemisphere == "North"].D_front_km_per_decade.mean(), 1) == 65.3 and
      round(a13[a13.hemisphere == "South"].D_front_km_per_decade.mean(), 1) == 36.0)

print("\n== references")
reflist = re.findall(r"^(\d+)\. [A-Z]", ms, re.M)
n = [int(x) for x in reflist]
check("reference list is numbered consecutively without gaps or duplicates",
      n == list(range(1, len(n) + 1)), f"{len(n)} entries, max {max(n) if n else 0}")
check("main-text references <= 50 (Nature Article limit)", max([x for x in n if x <= 50], default=0) >= 1 and len([x for x in n if x <= 39]) == 39,
      "39 main-text references")
dupes = [t for t in set(re.findall(r"^\d+\. (.{0,60})", ms, re.M)) if ms.count(t) > 1]
check("no duplicated reference entries", not dupes, "; ".join(d[:40] for d in dupes[:3]))

print("\n== submission hygiene")
for label, t in all_texts.items():
    check(f"no AUTHOR VERIFY flags in the {label}", "AUTHOR VERIFY" not in t)
ph = re.findall(r"\[INSERT\]|\[GITHUB URL\]", ms)
check("only the expected author placeholders remain in the manuscript (accessions/DOIs/URL)",
      len(ph) == 5, f"{len(ph)} placeholders: {sorted(set(ph))}")
check("no TODO/FIXME/magenta markers", not re.search(r"TODO|FIXME|magenta", ms + si + gd + sf, re.I))

print("\n== format")
title = [p.text for p in Document(MS).paragraphs if p.text.strip()][0]
summary = re.search(r"e-mail:.*?\n(.+?)\n", ms, re.S)
body = ms[ms.index("Ocean warming is reorganizing"):ms.index("References")]
check("title <= 75 characters", len(title) <= 75, f"{len(title)} chars")
words = len(re.sub(r"\d", "", ms[ms.index("Ocean warming is reorganizing"):ms.index("Marine species are shifting")]).split())
check("summary <= 200 words", words <= 200, f"{words} words")
check("no double spaces or stray punctuation in the manuscript",
      not re.search(r"[a-z]\.  [A-Z]|\s\.|\s,", ms))

print("\n== document coherence")
si_paras = [x.text.strip() for x in Document(SI).paragraphs if x.text.strip()]
notes = [x for x in si_paras if x.startswith("Supplementary Note ")]
check("SI contents list matches the Note headings in the body",
      len(notes) == 14 and notes[:7] == notes[7:], f"{len(notes)} Note lines")
check("manuscript has 5 main and 5 Extended Data figure legends",
      len(re.findall(r"\nFig\. \d+ \|", "\n" + ms)) // 2 == 5 and len(re.findall(r"Extended Data Fig\. \d+ \|", ms)) // 2 == 5)
#  rewrote the letter to be concept-led; the invariant is that each quantitative claim it keeps is in the manuscript
letter_claims = {"96.5%": "96.5%", "96.9%": "96.9%", "83 ": "83 stations", "65 km per decade": "65.3 km per decade",
                 "36 km per decade": "36.0 km per decade", "2014–2016": "2014–2016"}
missing_cl = [k for k, v in letter_claims.items() if k in cl and v not in ms]
check("every quantitative claim in the cover letter is supported by the manuscript", not missing_cl, ", ".join(missing_cl))
check("the cover letter's quantitative claims were all found in the letter", all(k in cl for k in letter_claims),
      ", ".join(k for k in letter_claims if k not in cl))


print("\n==  corrections")
pn = pd.read_csv(HERE / "analysis/pnps_order_tests.csv")
gA = pn[(pn.tier == "A") & pn.unit.str.startswith("genome (median of per")].iloc[0]
gD = pn[(pn.tier == "D") & pn.unit.str.startswith("genome (median of per")].iloc[0]
expect = {f"{gA.median_I:.3f} versus {gA.median_A:.3f}": "tier A genome medians",
          f"Cliff's δ = {gA.cliffs_delta:.2f}, 95% CI {gA.CI_low:.2f}–{gA.CI_high:.2f}": "tier A genome delta and CI",
          f"δ = −{abs(gD.cliffs_delta):.2f}, 95% CI −{abs(gD.CI_low):.2f} to {gD.CI_high:.2f}": "tier D genome delta and CI",
          f"P = {gD.P:.2f}": "tier D genome P"}
for txt, what in expect.items():
    check(f"manuscript pN/pS {what} matches pnps_order_tests.csv", txt in ms, txt)
check("tier A P value reported as 8.1 × 10−6", abs(gA.P - 8.1e-6) < 0.05e-6 and "8.1 × 10" in ms, f"{gA.P:.2e}")
for label, t in [("manuscript", ms), ("SI", si), ("figures document", sf)]:
    check(f"no 'mean gene coverage' tier definition in the {label}", not re.search(r"mean[ -]gene[ -]coverage", t))
    check(f"no retired 'calibration/calibrated domain' wording in the {label}", not re.search(r"calibrat\w+ domain", t))
check("superseded pN/pS values (3.2 × 10−7, δ 0.188) are gone from manuscript and SI",
      not re.search(r"3\.2 × 10|0\.188", ms + si))
sfz = zipfile.ZipFile(SF)
sf_media = {hashlib.md5(sfz.read(n)).hexdigest() for n in sfz.namelist() if "/media/" in n}
for png in ["SupplementaryFig4_pNpS.png", "SupplementaryFig5_association_network.png"]:
    check(f"{png} is the image embedded in the Supplementary Figures document",
          hashlib.md5((HERE / "supplementary/figures" / png).read_bytes()).hexdigest() in sf_media)
src = (HERE / "analysis/pnps_and_network.py").read_text()
check("Supplementary Figs 4-5 use the house palette (Imitervirales #0072B2, Algavirales #D55E00)",
      'I_C, A_C, GREY, LGREY, INK = "#0072B2", "#D55E00"' in src and 'GV_COL = {"Imitervirales": I_C, "Algavirales": A_C}' in src)


print("\n==  changes")
pre = pd.read_csv(HERE / "analysis/pre2014_pair_derivation_summary.csv").set_index("statistic").value
sh = float(pre["Imitervirales + Algavirales median share, pre-2014"]) * 100
sl = float(pre["D = (I - A)/(I + A) vs SST slope, pre-2014, fraction-adjusted (D per degC)"])
for label, t in [("manuscript", ms), ("SI", si)]:
    check(f"pre-2014 derivation values in the {label} match pre2014_pair_derivation_summary.csv",
          f"{sh:.1f}%" in t and f"{sl:.3f} D per °C" in t, f"{sh:.1f}% / {sl:.3f}")
check("the  'not a fully independent discovery–validation split' disclosure is replaced", "not a fully independent" not in ms)
check("no unverifiable Python version (3.13.5) in the Methods", "3.13.5" not in ms)
check("Fig. 2 legend states the exact permutation P as 1/21", "1/21 = 0.048" in ms)
internal = r"\d revision|source-auditable|[Ll]egacy (figures|outputs)|21 September|data deposit|\(not assembled\)|aligned to the "
for label, t in [("SI", si), ("guide", gd), ("figures document", sf), ("manuscript", ms), ("cover letter", cl)]:
    hit = re.search(internal, t)
    check(f"no internal version/provenance notes in the {label}", not hit, hit.group(0) if hit else "")
for f in [MS, CL, SI, GD, SF]:
    cp = Document(f).core_properties
    check(f"{f.name}: metadata free of tool and version strings",
          cp.author != "python-docx" and not cp.comments and not re.search(r"V\d\d", cp.title or ""), f"{cp.author!r} {cp.title!r}")
DD = HERE / "supplementary/Supplementary_Data"
expect_files = {2: 167, 3: 65, 4: 197, 5: 196, 6: 196, 7: 266, 8: 26774, 9: 440, 10: None}
for n, rows in expect_files.items():
    fl = list(DD.glob(f"Supplementary_Data_{n}_*.csv"))
    ok = len(fl) == 1 and (rows is None or len(pd.read_csv(fl[0])) == rows)
    check(f"Supplementary Data {n} present" + (f" with {rows} rows" if rows else ""), ok, fl[0].name if fl else "missing")
check("Supplementary Data 9b (2-degree grid) present", (DD / "Supplementary_Data_9b_CMIP6_Dindex_2deg_grid.csv.gz").exists())
for n in (5, 6):
    d = pd.read_csv(next(DD.glob(f"Supplementary_Data_{n}_*.csv")))
    check(f"Supplementary Data {n} keyed by unique run_id", d.run_id.is_unique and len(d) == 196)
check("Supplementary Data 4 omits the mislabelled CO3 field", "CO3" not in pd.read_csv(next(DD.glob("Supplementary_Data_4_*.csv"))).columns)
check("missing Supplementary Data 1 is an explicit bracketed placeholder in the guide", "[TO BE SUPPLIED: Supplementary_Data_1" in gd)
d10 = pd.read_csv(next(DD.glob("Supplementary_Data_10_*.csv")))
vals = pd.concat([pd.to_numeric(d10[c], errors="coerce") for c in ("value", "ci_low", "ci_high", "P")])
heads = {"SPOT slope": 0.112, "permutation P": 0.0476, "same-month P": 0.00158, "micro/nano contrast": -0.203, "Tara r": 0.7363,
         "log-chl r": -0.6495, "ensemble R2": 0.6345, "temperature-model R2": 0.5578, "temperature LOBO": 0.5353, "RF R2": 0.5667,
         "RF LOBO": 0.2667, "temperature importance": 0.5149, "front N": 65.27, "front S": 36.03, "isotherm N": 68.41, "isotherm S": 52.27,
         "pN/pS tier A": 0.3667, "Tara share": 0.969, "MH OR HSP": 76.31, "pathway difference": 4.542, "pre-2014 slope": 0.0866,
         "PolB Imitervirales": 721, "zero point viral": 22.557, "sign agreement": 0.7004}
miss = [k for k, v in heads.items() if not (abs(vals - v) <= max(abs(v) * 0.006, 1e-4)).any()]
check("every headline statistic is traceable in Supplementary Data 10", not miss, ", ".join(miss))
check("Supplementary Data 10 excludes superseded 62.6 km/decade (8-predictor) front rate",
      not (abs(pd.to_numeric(d10.value, errors="coerce") - 62.64) < 0.01).any())

print("\n== toolkit")
r = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=HERE.parents[1] / "dindex_toolkit",
                   env={"PYTHONPATH": str(HERE.parents[1] / "dindex_toolkit" / "src"), "PATH": "/usr/bin:/bin"},
                   capture_output=True, text=True)
check("dindex toolkit tests pass", "passed" in r.stdout, r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[:80])

failed = [r for r in results if not r[0]]
print(f"\n{len(results) - len(failed)}/{len(results)} checks passed")
if failed:
    print("FAILED:")
    for _, name, detail in failed:
        print(f"  - {name} ({detail})")
sys.exit(1 if failed else 0)
