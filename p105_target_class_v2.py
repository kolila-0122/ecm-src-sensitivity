# -*- coding: utf-8 -*-
"""
p105_target_class_v2.py
锚定（词边界）靶点类目分类器
对应稿件 1.4 节，表 1、表 S13 的 anchored 列

规则:
  词首: (?<![A-Z0-9])  非字母数字断言
  词尾: (?:[0-9]|[A-Z](?![A-Z0-9]))?  允许一个数字或单个字母后缀
  例: EGFR 匹配 EGFR / EGFR2 / EGFRB，但不匹配 VEGFR

用法:
  python p105_target_class_v2.py --input TableS1_drug_sensitivity_full.csv --output TableS13_match_sensitivity.csv
"""
import argparse
import re
import pandas as pd
from scipy import stats


ECM23 = ['COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2', 'COL6A1', 'COL6A2',
         'COL8A1', 'COL11A1', 'LOX', 'LOXL1', 'LOXL2', 'FN1', 'POSTN', 'TGFB1',
         'DCN', 'LUM', 'ACTA2', 'TAGLN', 'TIMP1', 'SPARC', 'CCN2', 'FBN1']


def token(w, anchor=True):
    """构造匹配模式: anchor=True 时加词边界"""
    if anchor:
        return r"(?<![A-Z0-9])" + w + r"(?:[0-9]|[A-Z](?![A-Z0-9]))?"
    return w


RULES_ANCHORED = [
    ("AURK/mitotic", [token("AURK"), token("AURORA")], True),
    ("Tubulin", [token("TUBULIN"), token("MICROTUBULE"), token("VINCA"), token("TAXANE")], False),
    ("BRAF", [token("BRAF")], False),
    ("CDK", [token("CDK")], False),
    ("MEK/MAPK", [token("MEK"), token("MAPK"), token("ERK")], False),
    ("PI3K/AKT/mTOR", [token("PI3K"), token("AKT"), token("MTOR")], False),
    ("EGFR/ERBB", [token("EGFR"), token("ERBB"), token("HER2")], False),
    ("WNT", [token("WNT"), token("PORCN"), token("TNKS")], False),
    ("Src family", [token("SRC"), token("LYN"), token("FYN")], False),
    ("HDAC", [token("HDAC")], False),
    ("IGF1R", [token("IGF1R"), token("INSR")], False),
    ("Apoptosis", [token("BCL"), token("BCL2"), token("APOPTOSIS")], False),
]


def classify(row, anchor=True):
    """按 1.4 节规则归类，anchor 控制词边界"""
    t = (str(row.get("Target", "")) + " " + str(row.get("Pathway", ""))).upper()
    for name, pats, extra in RULES_ANCHORED:
        if extra and "SL0101" in str(row.get("Drug", "")).upper():
            return name
        if anchor:
            if any(re.search(p, t) for p in pats):
                return name
        else:
            # 无锚定版本（对照）
            plain = [re.sub(r"\(\?<!\[A-Z0-9\]\)", "", p).replace(r"(?:[0-9]|[A-Z](?![A-Z0-9]))?", "", p) for p in pats]
            if any(re.search(p, t) for p in plain):
                return name
    return "Other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.input, comment="#")
    df["cluster_anchored"] = df.apply(lambda r: classify(r, True), axis=1)
    df["cluster_unanchored"] = df.apply(lambda r: classify(r, False), axis=1)

    # 显著性
    if "FDR (sci)" in df.columns:
        df["sig"] = df["FDR (sci)"].apply(lambda x: float(str(x).replace("<", "")) < 0.05)
    df["resistant"] = df["sig"] & (df["rho"] > 0)

    rows = []
    for cls in sorted(set(df.cluster_anchored) | set(df.cluster_unanchored)):
        for mode, col in (("unanchored", "cluster_unanchored"), ("anchored", "cluster_anchored")):
            m = df[col] == cls
            a = int((m & df.resistant).sum())
            b = int((m & ~df.resistant).sum())
            c = int((~m & df.resistant).sum())
            d = int((~m & ~df.resistant).sum())
            if a + b > 0:
                orr, p = stats.fisher_exact([[a, b], [c, d]])
            else:
                orr, p = float("nan"), float("nan")
            rows.append({"class": cls, "mode": mode, "n_class": int(m.sum()),
                         "n_resistant": a, "OR": orr, "p": p})

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False, encoding="utf-8")
    print(f"已写入 {args.output}（{len(out)} 行）")


if __name__ == "__main__":
    main()
