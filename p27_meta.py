"""ECM 与生存：Meta 分析 + 多因素 Cox + 三分位"""
import os, sys, json, time
import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
OUT = r"D:\stiffness_r3\results"

d = pd.read_csv(f"{OUT}/ecm_survival_tcga.csv")
print("=== 各队列 ===")
print(d.round(3).to_string(index=False))

# ---- 固定效应 Meta（逆方差法）----
print("\n=== Meta 分析（固定效应，逆方差）===")
loghr = np.log(d.HR_per_SD.values)
se = (np.log(d.CI_high.values) - np.log(d.CI_low.values)) / (2 * 1.96)
w = 1 / se ** 2
pooled = (w * loghr).sum() / w.sum()
pooled_se = np.sqrt(1 / w.sum())
z = pooled / pooled_se
p_pool = 2 * (1 - stats.norm.cdf(abs(z)))
print(f"合并 HR/SD = {np.exp(pooled):.3f}  (95% CI {np.exp(pooled - 1.96 * pooled_se):.3f}-{np.exp(pooled + 1.96 * pooled_se):.3f})")
print(f"z = {z:.2f}, p = {p_pool:.4f}")
print(f"异质性 Q = {float((w * (loghr - pooled) ** 2).sum()):.2f}, "
      f"df = {len(d) - 1}, "
      f"I² = {max(0, (1 - (len(d) - 1) / float((w * (loghr - pooled) ** 2).sum()))) * 100:.1f}%")

# 符号检验
k = int((d.HR_per_SD > 1).sum())
print(f"\n方向一致性: {k}/{len(d)} 个队列 HR>1，二项检验 p = {stats.binomtest(k, len(d), 0.5, alternative='greater').pvalue:.4f}")

# ---- 多因素 Cox ----
print("\n=== 多因素 Cox（控制年龄、性别、分期）===")
B = "https://www.cbioportal.org/api"
S = requests.Session(); S.headers.update({"accept": "application/json"})
from lifelines import CoxPHFitter

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}
ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

e2s = {}
for g in ECM:
    r = S.get(f"{B}/genes/{g}", timeout=20)
    if r.status_code == 200:
        e2s[g] = r.json()["entrezGeneId"]
sym = {v: k for k, v in e2s.items()}

res = []
for study, label in STUDIES.items():
    profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
    pid = None
    for p in profs:
        if p.get("molecularAlterationType") == "MRNA_EXPRESSION" and "zscore" in p["molecularProfileId"].lower() \
           and "mirna" not in p["molecularProfileId"].lower():
            pid = p["molecularProfileId"]; break
    if not pid:
        continue
    body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{study}_all"}
    rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                data=json.dumps(body), timeout=120)
    if rr.status_code != 200:
        continue
    df = pd.DataFrame(rr.json())
    X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
    X.index = [sym.get(i, str(i)) for i in X.index]
    score = X.loc[[g for g in ECM if g in X.index]].mean(0)

    C = S.get(f"{B}/studies/{study}/clinical-data",
              params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=60)
    C = pd.DataFrame(C.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                           values="value", aggfunc="first")
    smap = pd.Series(score.index.str[:12], index=score.index)
    sc = score.groupby(smap).mean()
    cols = [c for c in ["OS_MONTHS", "OS_STATUS", "AGE", "SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE"] if c in C.columns]
    T = pd.DataFrame({"score": sc}).join(C[cols], how="inner")
    T["OS_MONTHS"] = pd.to_numeric(T["OS_MONTHS"], errors="coerce")
    T["AGE"] = pd.to_numeric(T["AGE"], errors="coerce")
    T = T.dropna(subset=["OS_MONTHS", "OS_STATUS"])
    T = T[T.OS_MONTHS > 0]
    T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
    T["sex_m"] = (T.SEX.astype(str).str.upper().str[0] == "M").astype(int)
    if "AJCC_PATHOLOGIC_TUMOR_STAGE" in T.columns:
        st = T.AJCC_PATHOLOGIC_TUMOR_STAGE.astype(str).str.extract(r"(\d)")[0]
        T["stage"] = pd.to_numeric(st, errors="coerce")
    else:
        T["stage"] = np.nan
    T = T.dropna(subset=["AGE"])
    cov = ["score", "AGE", "sex_m"] + (["stage"] if T.stage.notna().sum() > 0.7 * len(T) else [])
    T2 = T[cov + ["OS_MONTHS", "event"]].dropna()
    if len(T2) < 80 or T2.event.sum() < 15:
        print(f"{label}: 样本不足 n={len(T2)}"); continue
    try:
        cph = CoxPHFitter()
        cph.fit(T2, duration_col="OS_MONTHS", event_col="event")
        hr = float(np.exp(cph.params_["score"]))
        ci = cph.confidence_intervals_.loc["score"].values
        p = float(cph.summary.loc["score", "p"])
        res.append({"study": label, "n": len(T2), "covariates": "+".join(cov[1:]),
                    "HR_adj": hr, "CI_low": float(np.exp(ci[0])), "CI_high": float(np.exp(ci[1])),
                    "p": p})
        print(f"{label:12s} n={len(T2):4d} 调整后 HR/SD={hr:.2f} "
              f"({np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f}) p={p:.3f}  [校正: {'+'.join(cov[1:])}]")
    except Exception as e:
        print(f"{label}: Cox 失败 {type(e).__name__}")

if res:
    R = pd.DataFrame(res)
    R.to_csv(f"{OUT}/ecm_survival_multivar.csv", index=False, encoding="utf-8-sig")
    print("\n多因素汇总:")
    print(R.round(3).to_string(index=False))
    lg = np.log(R.HR_adj.values)
    se2 = (np.log(R.CI_high.values) - np.log(R.CI_low.values)) / (2 * 1.96)
    w2 = 1 / se2 ** 2
    pl = (w2 * lg).sum() / w2.sum()
    ps = np.sqrt(1 / w2.sum())
    print(f"\n调整后合并 HR = {np.exp(pl):.3f} "
          f"({np.exp(pl - 1.96 * ps):.3f}-{np.exp(pl + 1.96 * ps):.3f}), "
          f"p = {2 * (1 - stats.norm.cdf(abs(pl / ps))):.4f}")
