"""ECM 模块与新辅助化疗 pCR"""
import os, sys, gzip, re
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

sys.stdout.reconfigure(line_buffering=True)
D = r"D:\stiffness_r3\data\GEO_neoadj"
OUT = r"D:\stiffness_r3\results"

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

MAP = pd.read_csv(os.path.join(D, "GPL96_map.csv"))
p2s = dict(zip(MAP.probe, MAP.symbol))
# 每个基因取方差最大的探针
best = {}
for g in ECM:
    ps = MAP.loc[MAP.symbol == g, "probe"].tolist()
    if ps:
        best[g] = ps[0]

print(f"ECM 基因可用 {len(best)}/{len(ECM)}\n")

# 各队列的 pCR 字段
PCR_FIELD = {
    "GSE20194": "pcr_vs_rd",
    "GSE20271": "pcr or rd",
    "GSE22093": "pcr",
    "GSE25055": "pcr",
}


def read_matrix(gse):
    """读 series matrix，返回 (expr DataFrame[probe x sample], meta dict)"""
    meta = {}
    lines = []
    with gzip.open(os.path.join(D, f"{gse}.gz"), "rt", encoding="utf-8", errors="ignore") as f:
        in_tab = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_tab = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_tab:
                if line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_geo_accession"):
                    parts = line.rstrip("\n").split("\t")[1:]
                    parts = [p.strip('"') for p in parts]
                    key = parts[0].split(":")[0].lower() if ":" in parts[0] else line.split("\t")[0]
                    meta.setdefault(key, []).append(parts)
                continue
            lines.append(line.rstrip("\n"))
    if not lines:
        return None, meta
    hdr = [c.strip('"') for c in lines[0].split("\t")]
    rows = [l.split("\t") for l in lines[1:] if l.strip()]
    X = pd.DataFrame([[float(v) if v not in ("", "null", "NA") else np.nan for v in r[1:]] for r in rows],
                     index=[r[0].strip('"') for r in rows], columns=hdr[1:])
    return X, meta


def get_pcr(meta, field):
    """从 characteristics 里提 pCR"""
    vals = None
    for key, blocks in meta.items():
        for b in blocks:
            if b and field.lower() in b[0].lower():
                vals = b
                break
        if vals is not None:
            break
    if vals is None:
        return None
    out = []
    for v in vals:
        s = v.split(":")[-1].strip().upper() if ":" in v else v.strip().upper()
        if "PCR" in s and "RD" not in s and "RESIDUAL" not in s:
            out.append(1)
        elif "RD" in s or "RESIDUAL" in s or "NOT" in s:
            out.append(0)
        else:
            out.append(np.nan)
    return pd.Series(out, index=[c.strip('"') for c in meta.get("!sample_geo_accession", [[""]])[0]] if "!sample_geo_accession" in meta else None)


rows = []
for gse, field in PCR_FIELD.items():
    fp = os.path.join(D, f"{gse}.gz")
    if not os.path.exists(fp):
        continue
    X, meta = read_matrix(gse)
    if X is None or X.empty:
        print(f"{gse}: 矩阵为空"); continue
    # 样本名
    samp = list(X.columns)
    # pCR
    pcr = None
    for key, blocks in meta.items():
        for b in blocks:
            if b and field.lower() in str(b[0]).lower():
                ser = []
                for v in b:
                    s = str(v).split(":")[-1].strip().upper()
                    if "PCR" in s and "RD" not in s and "RESIDUAL" not in s:
                        ser.append(1)
                    elif "RD" in s or "RESIDUAL" in s:
                        ser.append(0)
                    else:
                        ser.append(np.nan)
                pcr = pd.Series(ser)
                break
        if pcr is not None and len(pcr) == len(samp):
            break
    if pcr is None or len(pcr) != len(samp):
        # 尝试其他字段
        for key, blocks in meta.items():
            for b in blocks:
                if b and ("pcr" in str(b[0]).lower()):
                    ser = []
                    for v in b:
                        s = str(v).split(":")[-1].strip().upper()
                        if "PCR" in s and "RD" not in s:
                            ser.append(1)
                        elif "RD" in s or "RESIDUAL" in s:
                            ser.append(0)
                        else:
                            ser.append(np.nan)
                    if len(ser) == len(samp):
                        pcr = pd.Series(ser); break
            if pcr is not None and len(pcr) == len(samp):
                break
    if pcr is None or len(pcr) != len(samp):
        print(f"{gse}: 未找到 pCR (样本 {len(samp)}, meta keys {list(meta)[:4]})"); continue
    pcr.index = samp

    # 表达
    Xl = X.copy()
    if Xl.values.max() > 50:      # 线性 → log2
        Xl = np.log2(Xl.clip(lower=1))
    sub = Xl.loc[[p for p in best.values() if p in Xl.index]]
    if sub.shape[0] < 15:
        print(f"{gse}: 探针不足 {sub.shape[0]}"); continue
    sub.index = [MAP.set_index("probe").symbol.get(p, p) for p in sub.index]
    score = sub.mean(0)

    d = pd.DataFrame({"ecm": score, "pcr": pcr}).dropna()
    if len(d) < 40 or d.pcr.nunique() < 2:
        print(f"{gse}: 有效 {len(d)} 样本不足"); continue
    n1 = int((d.pcr == 1).sum()); n0 = int((d.pcr == 0).sum())
    auc = roc_auc_score(d.pcr, d.ecm)
    mw = stats.mannwhitneyu(d.ecm[d.pcr == 1], d.ecm[d.pcr == 0], alternative="two-sided")
    rows.append({"study": gse, "n": len(d), "pCR": n1, "RD": n0, "AUC": auc,
                 "med_pCR": float(d.ecm[d.pcr == 1].median()),
                 "med_RD": float(d.ecm[d.pcr == 0].median()), "p": float(mw.pvalue)})
    print(f"{gse:10s} n={len(d):3d} (pCR={n1}, RD={n0})  AUC={auc:.3f}  "
          f"ECM 中位 pCR={d.ecm[d.pcr==1].median():.2f} vs RD={d.ecm[d.pcr==0].median():.2f}  p={mw.pvalue:.4f}")

R = pd.DataFrame(rows)
if len(R):
    R.to_csv(f"{OUT}/ecm_pcr.csv", index=False, encoding="utf-8-sig")
    print("\n" + R.round(3).to_string(index=False))
    w = R.n.values
    z = np.arctanh(np.clip(R.AUC.values, .01, .99))
    zb = (w * z).sum() / w.sum()
    se = 1 / np.sqrt((w - 3).sum())
    print(f"\n加权合并 AUC = {np.tanh(zb):.3f} "
          f"({np.tanh(zb - 1.96*se):.3f}-{np.tanh(zb + 1.96*se):.3f})")
    print(f"AUC<0.5 的队列 {int((R.AUC < .5).sum())}/{len(R)}（ECM 高→pCR 率低）")
    print(f"总 {int(R.n.sum())} 例，pCR {int(R.pCR.sum())}")
