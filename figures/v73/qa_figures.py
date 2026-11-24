"""Objective checks for every figure PDF in a folder: width, height, text size range, panel-label size, fonts.
usage: python figures/v73/qa_figures.py [FOLDER]
  FOLDER defaults to <repo>/outputs/v73 (the output of build_main.py and build_ed.py).
  Nature production spec: 89 or 183 mm wide, <= 170 mm tall, text 5-7 pt, panel labels 8 pt bold, Arial."""
import sys, glob, os
from pathlib import Path
import pymupdf as fitz

FOLDER = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[2] / 'outputs' / 'v73')
files = sorted(glob.glob(os.path.join(FOLDER, '*.pdf')), key=lambda x: ('Extended' in x, x))
if not files:
    sys.exit(f'no figure PDFs in {FOLDER}; run build_main.py and build_ed.py first')
for f in files:
    p = fitz.open(f)[0]; w, h = p.rect.width / 72 * 25.4, p.rect.height / 72 * 25.4
    sp = [(s['size'], s['text'].strip(), s['font']) for b in p.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans'] if s['text'].strip()]
    sz = [a for a, _, _ in sp]; lab = sorted({round(a, 1) for a, t, fo in sp if len(t) == 1 and t in 'abcdefgh' and 'Bold' in fo})
    bad = []
    if abs(w - 183) > 1 and abs(w - 89) > 1: bad.append(f'width {w:.0f} mm')
    if h > 170: bad.append(f'height {h:.0f} mm > 170')
    if min(sz) < 4.95: bad.append(f'{sum(x < 4.95 for x in sz)} text runs < 5 pt (min {min(sz):.1f})')
    if lab and lab != [8.0]: bad.append(f'panel labels {lab} pt, not 8')
    if any('Arial' not in fo and 'Helvetica' not in fo for _, _, fo in sp): bad.append('non-Arial font')
    print(f"{os.path.basename(f)[:-4]:40s} {w:4.0f} x{h:4.0f} mm  text {min(sz):.1f}-{max(sz):.1f} pt  {'OK' if not bad else 'FIX: ' + '; '.join(bad)}")
