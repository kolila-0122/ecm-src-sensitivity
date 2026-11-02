"""ECM 模块与生存：TCGA 多癌种 Cox 回归"""
import os, sys, json, time
import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
OUT = r"D:\stiffness_r3\results"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}


def get_genes(symbols):
    out = {}
    for s in symbols:
        try:
            r = S.get(f"{B}/genes/{s}", timeout=25)
            if r.status_code == 200:
                out[s] = r.json()["entrezGeneId"]
        except Exception:
            pass
        time.sleep(0.1)
    return out


def profiles(study):
    r = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40)
    return r.json() if r.status_code == 200 else []


def pick_mrna(profs):
    cands = [p["molecularProfileId"] for p in profs
             if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
             and "zscore" in p["molecularProfileId"].lower()
             and "mirna" not in p["molecularProfileId"].lower()]
    if cands:
        return cands[0], True
    cands = [p["molecularProfileId"] for p in profs
             if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
             and "mirna" not in p["molecularProfileId"].lower()]
    return (cands[0], False) if cands else (None, False)


def fetch_expr(pid, entrez, study):
    body = {"entrezGeneIds": list(entrez), "sampleListId": f"{study}_all"}
    r = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
               data=json.dumps(body), timeout=120)
    if r.status_code != 200:
        return None
    df = pd.DataFrame(r.json())
    if df.empty:
        return None
    return df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")


def fetch_clin(study):
    r = S.get(f"{B}/studies/{study}/clinical-data",
              params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=60)
    if r.status_code != 200:
        return None
    df = pd.DataFrame(r.json())
    if df.empty:
        return None
    return df.pivot_table(index="patientId", columns="clinicalAttributeId",
                          values="value", aggfunc="first")


e2s = get_genes(ECM)
sym = {v: k for k, v in e2s.items()}
print(f"ECM 基因映射 {len(e2s)}/{len(ECM)}\n")

rows = []
for study, label in STUDIES.items():
    profs = profiles(study)
    if not profs:
        print(f"{label}: 无法访问"); continue
    pid, is_z = pick_mrna(profs)
    if not pid:
        print(f"{label}: 无 mRNA 谱"); continue
    X = fetch_expr(pid, e2s.values(), study)
    C = fetch_clin(study)
    if X is None or C is None:
        print(f"{label}: 取数失败"); continue
    X.index = [sym.get(i, str(i)) for i in X.index]
    gs = [g for g in ECM if g in X.index]
    X = X.loc[gs]
    if len(gs) < 15:
        print(f"{label}: 基因不足 {len(gs)}"); continue

    score = X.mean(0)   # 若已是 z-score，均值即模块分
    # 样本→患者
    smap = pd.DataFrame({"sampleId": score.index})
    smap["patientId"] = smap["sampleId"].str[:12]
    sc = score.groupby(smap.set_index("sampleId")["patientId"]).mean()

    need = ["OS_MONTHS", "OS_STATUS"]
    if not all(c in C.columns for c in need):
        print(f"{label}: 缺生存字段"); continue
    d = pd.DataFrame({"score": sc}).join(C[need], how="inner")
    d["OS_MONTHS"] = pd.to_numeric(d["OS_MONTHS"], errors="coerce")
    d = d.dropna(subset=["OS_MONTHS", "OS_STATUS"])
    d = d[d.OS_MONTHS > 0]
    if len(d) < 60:
        print(f"{label}: 样本不足 {len(d)}"); continue
    # 事件：1=死亡
    ev = d.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)

    # 单因素 Cox（分数连续）
    try:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter()
        tmp = pd.DataFrame({"T": d.OS_MONTHS.values, "E": ev.values,
                            "score": d.score.values})
        cph.fit(tmp, duration_col="T", event_col="E")
        hr = float(np.exp(cph.params_["score"]))
        ci = cph.confidence_intervals_.loc["score"].values
        p = float(cph.summary.loc["score", "p"])
    except Exception as e:
        print(f"{label}: Cox 失败 {type(e).__name__} {e}"); continue

    # 中位分组
    hi = d.score > d.score.median()
    from lifelines.statistics import logrank_test
    lr = logrank_test(d.OS_MONTHS[hi], d.OS_MONTHS[~hi],
                      ev[hi], ev[~hi])

    rows.append({"study": label, "n": len(d), "events": int(ev.sum()),
                 "HR_per_SD": hr, "CI_low": float(np.exp(ci[0])),
                 "CI_high": float(np.exp(ci[1])), "p": p,
                 "logrank_p": float(lr.p_value),
                 "median_OS_high": float(d.OS_MONTHS[hi].median()),
                 "median_OS_low": float(d.OS_MONTHS[~hi].median())})
    print(f"{label:12s} n={len(d):4d} 事件={int(ev.sum()):3d}  "
          f"HR/SD={hr:.2f} ({np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f})  p={p:.1e}  "
          f"logrank p={lr.p_value:.1e}")
    time.sleep(0.4)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/ecm_survival_tcga.csv", index=False, encoding="utf-8-sig")
print(f"\n完成 {len(df)} 个队列")
if len(df):
    print(df.round(3).to_string(index=False))
