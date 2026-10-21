# -*- coding: utf-8 -*-
"""静态检查：脚本引用的输入表是否都在本仓库里。

判定规则：
  * 仓库里已有的表格 -> OK
  * 投稿编号前缀（Supplementary_Data_7_、Supplementary_Table_3_ 等）与脚本内部短名互相匹配 -> OK
  * 脚本自己写出来的文件（to_csv/to_excel/savez/open）-> 视为可生成，不算缺失
  * 其余 -> 列为"仓库外输入"（大体积网格或外部原始数据）
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

# 第 1 遍：收集"本仓库某处会写出这个名字"的集合（跨脚本也算，A 写 B 读很常见）
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

# 第 2 遍：报告缺失的输入
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

print("仓库内表格文件: %d 个" % len(HAVE))
print("脚本引用的输入: %d 处；脚本管线自行生成的文件名 %d 个（已排除）" % (checked, len(generated)))
if missing:
    print("\n仓库外输入 %d 个（体积过大的网格或外部原始数据，见 RUN.md）:" % len(missing))
    for b, fs in sorted(missing.items()):
        print("  %-48s %d 个脚本引用，例: %s" % (b, len(fs), fs[0]))
else:
    print("\n所有被引用的输入都在仓库内 ✓")
sys.exit(0)
