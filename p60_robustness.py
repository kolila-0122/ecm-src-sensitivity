"""短板1：ECM 评分的稳健性分析（不同构建方式 + 分位分组）"""
import os, sys, glob, json
import numpy as np
import pandas as pd
import requests
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
OUT = r"D:\stiffness_r3\results"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}

# 不同基因子集（预先定义，非数据驱动）
SUBSETS = {
    "All 23 genes": None,
    "Collagens only": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL8A1", "COL11A1"],
    "Cross-linkers only": ["LOX", "LOXL1", "LOXL2"],
    "No myofibroblast markers": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
                                  "COL8A1", "COL11A1", "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "TGFB1",
                                  "DCN", "LUM", "TIMP1", "SPARC", "CCN2", "FBN1"],
    "Core 12 (collagen+LOX+FN1+POSTN)": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1",
                                          "COL6A2", "LOX", "LOXL2", "FN1", "POSTN", "TGFB1"],
}

ALL_GENES = set()
for v in SUBSETS.values():
    if v:
        ALL_GENES.update(v)

fc = os.path.join(CACHE, "robust_expr.csv")
if os.path.exists(fc):
    X_all = pd.read_csv(fc, index_col=0, header=[0, 1])
    print(f"载入缓存 {X_all.shape}")
else:
    X_all = None

print("=== 拉取各队列基因表达 ===")
frames = {}
for study, label in STUDIES.items():
    try:
        profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
        pid = None
        for p in profs:
            if p.get("molecularAlterationType") == "MRNA_EXPRESSION" and "zscore" in p["molecularProfileId"].lower() \
               and "mirna" not in p["molecularProfileId"].lower():
                pid = p["molecularProfileId"]; break
        if not pid:
            print(f"  {label}: 无谱"); continue
        g = S.get(f"{B}/genes", timeout=120).json()
        gmap = {x["hugoGeneSymbol"].upper(): x["entrezGeneId"] for x in g if x.get("hugoGeneSymbol")}
        ids = [gmap[x] for x in ALL_GENES if x in gmap]
        inv = {v: k for k, v in gmap.items()}
        body = {"entrezGeneIds": ids, "sampleListId": f"{study}_all"}
        r = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch", data=json.dumps(body), timeout=180)
        df = pd.DataFrame(r.json())
        X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
        X.index = [inv.get(i, str(i)) for i in X.index]
        # 临床
        cc = S.get(f"{B}/studies/{study}/clinical-data",
                   params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=90)
        C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                                values="value", aggfunc="first")
        keep = [c for c in ["OS_MONTHS", "OS_STATUS", "AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE"] if c in C.columns]
        sm = S.get(f"{B}/studies/{study}/samples", timeout=40).json()
        m = {x["sampleId"]: x["patientId"] for x in sm if x.get("sampleId") in set(X.columns)}
        X = X.rename(columns=m).groupby(level=0, axis=1).mean()
        common = X.columns.intersection(C.index)
        X = X[common]; C = C.loc[common]
        frames[label] = {"X": X, "C": C[keep]}
        print(f"  {label}: {X.shape[0]} genes × {X.shape[1]} samples")
    except Exception as e:
        print(f"  {label}: 失败 {type(e).__name__}")

# ============ 稳健性分析 ============
rows = []
print("\n=== 不同构建方式的 pooled HR ===")
for sname, genes in SUBSETS.items():
    use = genes if genes else ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
                                "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
                                "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
    hrs, ses = [], []
    for label, d in frames.items():
        X, C = d["X"], d["C"]
        gs = [g for g in use if g in X.index]
        if len(gs) < 3:
            continue
        score = X.loc[gs].mean(0)
        T = pd.DataFrame({"score": score}).join(C, how="inner")
        T["OS_MONTHS"] = pd.to_numeric(T.OS_MONTHS, errors="coerce")
        T["AGE"] = pd.to_numeric(T.AGE, errors="coerce")
        T = T.dropna(subset=["OS_MONTHS", "OS_STATUS", "score"])
        T = T[T.OS_MONTHS > 0]
        if len(T) < 60:
            continue
        T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
        T["z"] = (T.score - T.score.mean()) / (T.score.std() + 1e-9)
        cols = ["z"]
        if T.AGE.notna().sum() > len(T) * .6:
            T["AGE_z"] = (T.AGE - T.AGE.mean()) / (T.AGE.std() + 1e-9); cols.append("AGE_z")
        try:
            c = CoxPHFitter()
            c.fit(T[cols + ["OS_MONTHS", "event"]], duration_col="OS_MONTHS", event_col="event")
            beta = c.params_["z"]; se = c.standard_errors_["z"]
            hrs.append(beta); ses.append(se)
        except Exception:
            continue
    if len(hrs) >= 3:
        w = np.array([1 / s ** 2 for s in ses])
        pb = (w * np.array(hrs)).sum() / w.sum()
        ps = np.sqrt(1 / w.sum())
        Q = float((w * (np.array(hrs) - pb) ** 2).sum())
        I2 = max(0, (Q - (len(hrs) - 1)) / Q) * 100 if Q > 0 else 0
        hr = float(np.exp(pb))
        lo, hi = float(np.exp(pb - 1.96 * ps)), float(np.exp(pb + 1.96 * ps))
        z = pb / ps
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        rows.append({"Subset": sname, "k cohorts": len(hrs), "HR": hr, "lo": lo, "hi": hi,
                     "p": p, "I2": I2})
        print(f"  {sname:34s} k={len(hrs)}  HR={hr:.3f} ({lo:.3f}-{hi:.3f})  p={p:.4f}  I²={I2:.0f}%")

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "robustness_ecm_subsets.csv"), index=False, encoding="utf-8-sig")
print("\n[saved robustness_ecm_subsets.csv]")

# ============ 分位分组（展示剂量效应） ============
print("\n=== 三分位分组的 log-rank（合并六队列） ===")
allT = []
for label, d in frames.items():
    X, C = d["X"], d["C"]
    use = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
           "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
           "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
    gs = [g for g in use if g in X.index]
    score = X.loc[gs].mean(0)
    T = pd.DataFrame({"score": score}).join(C, how="inner")
    T["OS_MONTHS"] = pd.to_numeric(T.OS_MONTHS, errors="coerce")
    T = T.dropna(subset=["OS_MONTHS", "OS_STATUS", "score"])
    T = T[T.OS_MONTHS > 0]
    if len(T) < 60: continue
    T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
    T["tertile"] = pd.qcut(T.score, 3, labels=["Low", "Mid", "High"])
    T["cohort"] = label
    allT.append(T[["cohort", "OS_MONTHS", "event", "tertile", "score"]])

if allT:
    A = pd.concat(allT)
    print(f"  总样本 {len(A)}, 事件 {int(A.event.sum())}")
    for coh, g in A.groupby("cohort"):
        lr = multivariate_logrank_test(g.OS_MONTHS, g.tertile, g.event)
        med = g.groupby("tertile", observed=True).OS_MONTHS.median()
        print(f"  {coh:12s} n={len(g):4d}  中位OS L/M/H = {med['Low']:.0f}/{med['Mid']:.0f}/{med['High']:.0f}  "
              f"log-rank p={lr.p_value:.3f}")
    A.to_csv(os.path.join(OUT, "tertile_survival.csv"), index=False, encoding="utf-8-sig")
