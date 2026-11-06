#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
matchbias — 靶点类目关键词匹配偏倚自检工具
==========================================

检测化合物-靶点注释中，无词边界（unanchored）关键词匹配是否把
包含目标子串的其他基因误归入某类目（如 VEGFR⊃EGFR、PERK⊃ERK）。

用法
----
    python matchbias.py --input compounds.csv --target-col Target --name-col Drug
    python matchbias.py --input compounds.csv --target-col Target --rules rules.tsv

输入
----
CSV/TSV，至少含两列：化合物名、靶点注释（可含多个靶点，分号分隔）。
若提供 `--pathway-col`，可一并匹配通路名。

输出
----
1. 每个类目的两种计数（无锚定 vs 词边界）
2. **受影响的类目清单**（哪些类目因匹配方式而改变规模）
3. **被误归的化合物清单**（无锚定进了类目、锚定后移出）
4. 退出码：0 = 无偏倚；1 = 检出偏倚

作者：郑家惟 / 连彬
许可：MIT
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------- 默认规则
# 每条：(类目名, [基因符号...], 特殊标记)
# 词边界模式：词首 (?<![A-Z0-9])，词尾允许一个数字或单字母后缀
DEFAULT_RULES: list[tuple[str, list[str], str | None]] = [
    ("AURK/mitotic", ["AURK", "AURORA"], "SL0101"),
    ("Tubulin", ["TUBULIN", "MICROTUBULE", "VINCA", "TAXANE"], None),
    ("BRAF", ["BRAF"], None),
    ("CDK", ["CDK"], None),
    ("MEK/MAPK", ["MEK", "MAPK", "ERK"], None),
    ("PI3K/AKT/mTOR", ["PI3K", "AKT", "MTOR"], None),
    ("EGFR/ERBB", ["EGFR", "ERBB", "HER2"], None),
    ("WNT", ["WNT", "PORCN", "TNKS"], None),
    ("Src family", ["SRC", "LYN", "FYN"], None),
    ("HDAC", ["HDAC"], None),
    ("IGF1R", ["IGF1R", "INSR"], None),
    ("Apoptosis", ["BCL", "BCL2", "APOPTOSIS"], None),
]


def anchored_pattern(gene: str) -> re.Pattern:
    """词边界模式：词首非字母数字 + 允许一个数字/单字母后缀"""
    return re.compile(
        r"(?<![A-Z0-9])" + gene + r"(?:[0-9]|[A-Z](?![A-Z0-9]))?"
    )


def unanchored_pattern(gene: str) -> re.Pattern:
    return re.compile(gene)


def classify(
    text: str,
    rules: list[tuple[str, list[str], str | None]],
    drug: str,
    anchor: bool,
) -> str:
    """按顺序匹配，首个命中的类目即结果"""
    for name, genes, special in rules:
        if special and special.upper() in drug.upper():
            return name
        for g in genes:
            pat = anchored_pattern(g) if anchor else unanchored_pattern(g)
            if pat.search(text):
                return name
    return "Other"


def load_rules(path: Path | None) -> list[tuple[str, list[str], str | None]]:
    if path is None:
        return DEFAULT_RULES
    rules = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                name = parts[0].strip()
                genes = [g.strip() for g in parts[1].split(";") if g.strip()]
                special = parts[2].strip() if len(parts) > 2 else None
                rules.append((name, genes, special))
    return rules


def read_table(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delim)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description="靶点类目关键词匹配偏倚自检",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--input", required=True, type=Path, help="化合物-靶点表（CSV/TSV）")
    ap.add_argument("--name-col", default="Drug", help="化合物名列名（默认 Drug）")
    ap.add_argument("--target-col", default="Target", help="靶点注释列名（默认 Target）")
    ap.add_argument("--pathway-col", default=None, help="通路列名（可选，一并参与匹配）")
    ap.add_argument("--rules", type=Path, default=None, help="自定义规则文件（TSV）")
    ap.add_argument("--quiet", action="store_true", help="只输出摘要")
    args = ap.parse_args()

    rules = load_rules(args.rules)
    header, rows = read_table(args.input)
    for col in (args.name_col, args.target_col):
        if col not in header:
            print(f"✗ 缺少列 {col}（实际列：{header}）", file=sys.stderr)
            return 2

    members: dict[str, dict[bool, set[str]]] = defaultdict(lambda: {False: set(), True: set()})
    for r in rows:
        drug = (r.get(args.name_col) or "").strip()
        text = " ".join(
            filter(None, [
                (r.get(args.target_col) or ""),
                (r.get(args.pathway_col) or "") if args.pathway_col else "",
            ])
        ).upper()
        if not drug or not text.strip():
            continue
        for anchor in (False, True):
            members[classify(text, rules, drug, anchor)][anchor].add(drug)

    # ---------- 报告 ----------
    classes = sorted(members)
    affected = [c for c in classes if members[c][False] != members[c][True]]
    moved = {
        c: sorted(members[c][False] - members[c][True]) for c in affected
    }

    print("=" * 60)
    print("matchbias — 靶点类目关键词匹配偏倚自检")
    print("=" * 60)
    print(f"输入文件      : {args.input}")
    print(f"化合物总数    : {len({d for c in classes for d in members[c][False] | members[c][True]})}")
    print(f"类目数        : {len(classes)}")
    print(f"受影响类目    : {len(affected)}")
    print()

    if not affected:
        print("✓ 未检出匹配偏倚：两种匹配方式归类结果一致。")
        return 0

    print(f"{'类目':<20}{'无锚定':>8}{'词边界':>8}{'差异':>8}{'被误归':>8}")
    print("-" * 60)
    for c in classes:
        a, b = len(members[c][False]), len(members[c][True])
        flag = "  ←" if c in affected else ""
        print(f"{c:<20}{a:>8}{b:>8}{a - b:>+8}{len(moved.get(c, [])):>8}{flag}")

    print()
    print("=" * 60)
    print("被误归化合物（无词边界匹配时进入类目、加词边界后移出）")
    print("=" * 60)
    for c in affected:
        if moved[c]:
            print(f"\n【{c}】共 {len(moved[c])} 种：")
            for d in moved[c]:
                print(f"    - {d}")

    print()
    print("=" * 60)
    print("建议：在方法部分报告匹配规则，并在补充材料公开类目成员。")
    print("=" * 60)
    return 1


if __name__ == "__main__":
    sys.exit(main())
