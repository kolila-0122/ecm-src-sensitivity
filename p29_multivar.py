"""步骤 2：多因素 Cox（本地，读缓存）"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

sys.stdout.reconfigure(line_buffering=True)
CACHE = r"D:\stiffness_r3\cache"
OUT = r"D:\stiffness_r3\results"

LAB = {"coadread": "Colon", "brca": "Breast", "luad": "Lung ad.",
       "lusc": "Lung sq.", "paad": "Pancreas", "ucec": "Endometrium"}

rows = []
for f in sorted(glob.glob(os.path.join(CACHE, "*.csv"))):
    key = os.path.basename(f).split("_")[0]
    label = LAB.get(key, key)
    T = pd.read_csv(f)
    T["OS_MONTHS"] = pd.to_numeric(T["OS_MONTHS"], errors="coerce")
    T["AGE"] = pd.to_numeric(T["AGE"], errors="coerce")
    T = T.dropna(subset=["OS_MONTHS", "OS_STATUS", "ecm_score"])
    T = T[T.OS_MONTHS > 0]
    T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
    T["sex_m"] = (T.SEX.astype(str).str.upper().str[0] == "M").astype(int)
    T["ecm_score"] = (T["ecm_score"] - T["ecm_score"].mean()) / (T["ecm_score"].std() + 1e-9)
    T["AGE"] = (T["AGE"] - T["AGE"].mean()) / (T["AGE"].std() + 1e-9)
    # 分期
    stcol = None
    for c in ["AJCC_PATHOLOGIC_TUMOR_STAGE", "NEOPLASM_DISEASEPATHOLOGIC_STAGE",
              "PATHOLOGIC_STAGE_AJCC7_LOWERCASE"]:
        if c in T.columns and T[c].notna().sum() > 0.5 * len(T):
            stcol = c; break
    if stcol:
        st = T[stcol].astype(str).str.upper().str.extract(r"(IV|III|II|I)")[0]
        mp = {"I": 1, "II": 2, "III": 3, "IV": 4}
        T["stage"] = st.map(mp)
    else:
        T["stage"] = np.nan
    T = T.dropna(subset=["AGE"])
    if len(T) < 60:
        print(f"{label}: 样本不足 {len(T)}"); continue

    # 单因素
    try:
        c1 = CoxPHFitter()
        d1 = T[["OS_MONTHS", "event", "ecm_score"]].copy()
        c1.fit(d1, duration_col="OS_MONTHS", event_col="event")
    except Exception as e:
        print(f"{label}: 单因素失败 {type(e).__name__}"); continue
    hr1 = float(np.exp(c1.params_["ecm_score"]))
    ci1 = c1.confidence_intervals_.loc["ecm_score"].values
    p1 = float(c1.summary.loc["ecm_score", "p"])

    # 多因素
    cov = ["ecm_score", "AGE", "sex_m"]
    if T.stage.notna().sum() > 0.6 * len(T):
        cov.append("stage")
    d2 = T[cov + ["OS_MONTHS", "event"]].dropna().copy()
    if len(d2) < 60 or d2.event.sum() < 12:
        print(f"{label}: 多因素样本不足"); continue
    # 去掉常数协变量
    for c in list(cov):
        if c in d2.columns and d2[c].nunique() < 2:
            d2 = d2.drop(columns=[c]); cov = [x for x in cov if x != c]
    try:
        c2 = CoxPHFitter()
        c2.fit(d2, duration_col="OS_MONTHS", event_col="event")
    except Exception as e:
        print(f"{label}: 多因素失败 {type(e).__name__} n={len(d2)}"); continue
    hr2 = float(np.exp(c2.params_["ecm_score"]))
    ci2 = c2.confidence_intervals_.loc["ecm_score"].values
    p2 = float(c2.summary.loc["ecm_score", "p"])

    # 中位分组 logrank
    hi = T.ecm_score > T.ecm_score.median()
    lr = logrank_test(T.OS_MONTHS[hi], T.OS_MONTHS[~hi], T.event[hi], T.event[~hi])

    rows.append({"study": label, "n_uni": len(T), "ev_uni": int(T.event.sum()),
                 "HR_uni": hr1, "lo_uni": float(np.exp(ci1[0])), "hi_uni": float(np.exp(ci1[1])),
                 "p_uni": p1,
                 "n_mv": len(d2), "ev_mv": int(d2.event.sum()),
                 "HR_adj": hr2, "lo_adj": float(np.exp(ci2[0])), "hi_adj": float(np.exp(ci2[1])),
                 "p_adj": p2,
                 "covs": "+".join(cov[1:]), "logrank_p": float(lr.p_value)})
    print(f"{label:12s} n={len(T):4d}  单因素 HR={hr1:.2f} p={p1:.3f}  |  "
          f"校正 n={len(d2):4d} HR={hr2:.2f} ({np.exp(ci2[0]):.2f}-{np.exp(ci2[1]):.2f}) p={p2:.3f}  "
          f"[{'+'.join(cov[1:])}]")


def meta(df, hr, lo, hi, tag):
    lg = np.log(df[hr].values)
    se = (np.log(df[hi].values) - np.log(df[lo].values)) / (2 * 1.96)
    w = 1 / se ** 2
    p = (w * lg).sum() / w.sum()
    ps = np.sqrt(1 / w.sum())
    z = p / ps
    Q = float((w * (lg - p) ** 2).sum())
    I2 = max(0, (Q - (len(df) - 1)) / Q) * 100 if Q > 0 else 0
    print(f"{tag}: 合并 HR={np.exp(p):.3f} ({np.exp(p - 1.96 * ps):.3f}-{np.exp(p + 1.96 * ps):.3f})  "
          f"z={z:.2f}  p={2 * (1 - stats.norm.cdf(abs(z))):.4f}  I²={I2:.1f}%")
    return np.exp(p), p


R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/ecm_survival_multivar.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 70)
print(R.round(3).to_string(index=False))
print("\n=== Meta ===")
meta(R, "HR_uni", "lo_uni", "hi_uni", "单因素  ")
meta(R, "HR_adj", "lo_adj", "hi_adj", "校正后  ")
k = int((R.HR_adj > 1).sum())
print(f"校正后方向一致: {k}/{len(R)}，二项检验 p={stats.binomtest(k, len(R), .5, alternative='greater').pvalue:.4f}")
