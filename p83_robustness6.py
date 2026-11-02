"""稳健性分析（修正版）：全部 6 队列 × 5 个基因子集"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
import requests
from scipy import stats
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")
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

SUBSETS = {
    "All 23 genes": ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","LOX","LOXL1","LOXL2",
                     "FN1","POSTN","CCN2","ACTA2","TAGLN","TIMP1","SPARC","DCN","LUM","TGFB1","COL11A1","COL8A1","FBN1"],
    "Collagens only (9)": ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL8A1","COL11A1"],
    "Cross-linkers only (3)": ["LOX","LOXL1","LOXL2"],
    "23 minus ACTA2/TAGLN (21)": ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","LOX","LOXL1","LOXL2",
                                   "FN1","POSTN","CCN2","TIMP1","SPARC","DCN","LUM","TGFB1","COL11A1","COL8A1","FBN1"],
    "Core 12": ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","LOX","LOXL2","FN1","POSTN","TGFB1"],
}

ALL = sorted({g for v in SUBSETS.values() for g in v})
g = S.get(f"{B}/genes", timeout=120).json()
gmap = {x["hugoGeneSymbol"].upper(): x["entrezGeneId"] for x in g if x.get("hugoGeneSymbol")}
inv = {v: k for k, v in gmap.items()}
ids = [gmap[x] for x in ALL if x in gmap]

# ---------- 取数 ----------
cache_f = os.path.join(CACHE, "robust6.pkl")
if os.path.exists(cache_f):
    frames = pd.read_pickle(cache_f)
    print(f"载入缓存 {len(frames)} 队列")
else:
    frames = {}
    for study, label in STUDIES.items():
        profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
        pid = next((p["molecularProfileId"] for p in profs
                    if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
                    and "zscore" in p["molecularProfileId"].lower()
                    and "mirna" not in p["molecularProfileId"].lower()), None)
        r = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                   data=json.dumps({"entrezGeneIds": ids, "sampleListId": f"{study}_all"}), timeout=180)
        X = pd.DataFrame(r.json()).pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
        X.index = [inv.get(i, str(i)) for i in X.index]
        cc = S.get(f"{B}/studies/{study}/clinical-data",
                   params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=90)
        C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
        sm = S.get(f"{B}/studies/{study}/samples", timeout=40).json()
        mp = {x["sampleId"]: x["patientId"] for x in sm if x.get("sampleId") in set(X.columns)}
        X = X.rename(columns=mp).T.groupby(level=0).mean().T
        common = X.columns.intersection(C.index)
        frames[label] = {"X": X[common], "C": C.loc[common, [c for c in ["OS_MONTHS","OS_STATUS","AGE"] if c in C.columns]]}
        print(f"  {label}: {len(common)} samples")
    pd.to_pickle(frames, cache_f)

# ---------- 稳健性 ----------
rows = []
detail = []
print("\n=== 六队列 × 五子集 ===")
for sname, genes in SUBSETS.items():
    hrs, ses, ks = [], [], []
    for label, d in frames.items():
        X, C = d["X"], d["C"]
        gs = [x for x in genes if x in X.index]
        if len(gs) < 3:
            continue
        T = pd.DataFrame({"score": X.loc[gs].mean(0)}).join(C, how="inner")
        T["OS_MONTHS"] = pd.to_numeric(T.OS_MONTHS, errors="coerce")
        T = T.dropna(subset=["OS_MONTHS", "OS_STATUS", "score"])
        T = T[T.OS_MONTHS > 0]
        T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
        T["z"] = (T.score - T.score.mean()) / (T.score.std() + 1e-9)
        if len(T) < 50 or T.event.sum() < 10:
            continue
        try:
            c = CoxPHFitter()
            c.fit(T[["z", "OS_MONTHS", "event"]], duration_col="OS_MONTHS", event_col="event")
            hrs.append(float(c.params_["z"]))
            ses.append(float(c.standard_errors_["z"]))
            ks.append(label)
            detail.append({"Subset": sname, "Cohort": label, "n": len(T),
                           "events": int(T.event.sum()),
                           "HR": np.exp(float(c.params_["z"]))})
        except Exception as e:
            print(f"  {sname} / {label}: 失败 {type(e).__name__}")
    if len(hrs) >= 2:
        w = 1 / np.array(ses) ** 2
        pb = (w * np.array(hrs)).sum() / w.sum()
        ps = np.sqrt(1 / w.sum())
        Q = float((w * (np.array(hrs) - pb) ** 2).sum())
        I2 = max(0, (Q - (len(hrs) - 1)) / Q) * 100 if Q > 0 else 0
        p = 2 * (1 - stats.norm.cdf(abs(pb / ps)))
        rows.append({"Gene subset": sname, "k cohorts": len(hrs),
                     "HR per SD": round(float(np.exp(pb)), 3),
                     "95% CI lower": round(float(np.exp(pb - 1.96 * ps)), 3),
                     "95% CI upper": round(float(np.exp(pb + 1.96 * ps)), 3),
                     "p": round(float(p), 4), "I2 (%)": round(I2, 0)})
        print(f"  {sname:28s} k={len(hrs)}  HR={np.exp(pb):.3f} "
              f"({np.exp(pb-1.96*ps):.3f}-{np.exp(pb+1.96*ps):.3f})  p={p:.4f}  I²={I2:.0f}%")

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "robustness_6cohorts.csv"), index=False)
D = pd.DataFrame(detail)
D.to_csv(os.path.join(OUT, "robustness_per_cohort.csv"), index=False)
print(f"\n[saved] k 值: {list(R['k cohorts'])}")
print("\n=== 各队列单变量 HR ===")
for coh in D[D.Subset == "All 23 genes"].itertuples():
    print(f"  {coh.Cohort:12s} n={coh.n:4d} ev={coh.events:3d} HR={coh.HR:.3f}")
