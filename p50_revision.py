"""按审稿意见全面整改（v2.0）"""
import os, re
import numpy as np
import pandas as pd

SRC = r"D:\workbuddy内容\A1_力学代理"
RES = r"D:\stiffness_r3\results"

# ============ 1. 重建 Table S1（科学计数法） ============
T = pd.read_csv(os.path.join(SRC, "TableS1_drug_sensitivity_full.csv"))
T.columns = ["Drug", "Target", "Pathway", "n cell lines", "rho", "p", "fdr", "rho_auc", "p_auc"]

def sci(x):
    if pd.isna(x):
        return ""
    if x == 0:
        return "<1e-300"
    return f"{x:.2e}"

T["p (sci)"] = T.p.map(sci)
T["FDR (sci)"] = T.fdr.map(sci)
T["|rho| rank"] = T.rho.abs().rank(ascending=False).astype(int)
T = T.sort_values("p")
out1 = T[["Drug", "Target", "Pathway", "n cell lines", "rho", "p (sci)", "FDR (sci)", "|rho| rank"]].copy()
out1["rho"] = out1.rho.round(4)
out1.to_csv(os.path.join(SRC, "TableS1_drug_sensitivity_full.csv"), index=False, encoding="utf-8")
print(f"Table S1 rebuilt: {len(out1)} rows, sci notation applied")

# 显著子集
sig = out1[out1["FDR (sci)"].apply(lambda s: float(s) if s and not s.startswith("<") else 0) < 0.05]
sig.to_csv(os.path.join(SRC, "TableS1b_significant_87.csv"), index=False, encoding="utf-8")
print(f"  significant: {len(sig)} (rho>0: {(sig.rho>0).sum()}, rho<0: {(sig.rho<0).sum()})")

# ============ 2. 重建 Table S5（UCLA log2 + 终点标注） ============
A = pd.read_csv(os.path.join(RES, "ecm_ici_all.csv"))
A.loc[A.study.str.contains("UCLA"), "ecm"] = np.log2(A.loc[A.study.str.contains("UCLA"), "ecm"].clip(lower=1))

from sklearn.metrics import roc_auc_score
END = {"Riaz 2017 (nivo)": "RESPONDER (RECIST)", "Hugo 2016 (pembro)": "RESPONDER (RECIST)",
       "Liu 2019 (anti-PD1)": "RESPONDER (RECIST)", "UCLA 2016 (pembro)": "DURABLE_CLINICAL_BENEFIT"}
rows = []
for st, g in A.groupby("study"):
    if g.resp.nunique() < 2:
        continue
    rows.append({"Study": st, "Endpoint": END.get(st, ""), "n": len(g),
                 "Responders": int(g.resp.sum()),
                 "ECM median (responders)": round(float(g.ecm[g.resp == 1].median()), 2),
                 "ECM median (non-responders)": round(float(g.ecm[g.resp == 0].median()), 2),
                 "AUC": round(float(roc_auc_score(g.resp, g.ecm)), 3)})
T5 = pd.DataFrame(rows)
T5.to_csv(os.path.join(SRC, "TableS5_immunotherapy.csv"), index=False, encoding="utf-8")
print("\nTable S5 rebuilt (log2 scale, endpoint labelled):")
print(T5.to_string(index=False))

# 分终点汇总
from scipy.stats import norm
print("\n--- 分终点 pooled AUC ---")
for end in T5.Endpoint.unique():
    sub = T5[T5.Endpoint == end]
    w = sub.n.values
    z = np.arctanh(np.clip(sub.AUC.values, .01, .99))
    zb = (w * z).sum() / w.sum()
    se = 1 / np.sqrt((w - 3).sum())
    print(f"  {end}: n={int(sub.n.sum())}, pooled AUC={np.tanh(zb):.3f} "
          f"({np.tanh(zb-1.96*se):.3f}-{np.tanh(zb+1.96*se):.3f})")

# ============ 3. 重建 Table S6（补 p 值 + 显著性标记） ============
T6 = pd.read_csv(os.path.join(SRC, "TableS6_chemotherapy_pcr.csv"))
T6["significant (p<0.05)"] = np.where(T6.p < 0.05, "yes", "no")
T6["direction"] = np.where(T6.AUC < 0.5, "ECM high → lower pCR", "ECM high → higher pCR")
T6.to_csv(os.path.join(SRC, "TableS6_chemotherapy_pcr.csv"), index=False, encoding="utf-8")
print("\nTable S6 rebuilt with significance and direction flags:")
print(T6[["study", "n", "AUC", "p", "significant (p<0.05)", "direction"]].round(4).to_string(index=False))

# ============ 4. 重建 Table S4（补缺失率说明） ============
T4 = pd.read_csv(os.path.join(SRC, "TableS4_prognosis_cohorts.csv"))
T4["covariates_note"] = np.where(T4.covs == "AGE",
    "sex and stage missing >40%; adjusted for age only",
    "age, sex and AJCC stage")
T4.to_csv(os.path.join(SRC, "TableS4_prognosis_cohorts.csv"), index=False, encoding="utf-8")
print("\nTable S4 rebuilt with covariate note")
print(T4[["study", "n_uni", "ev_uni", "n_mv", "ev_mv", "covs", "covariates_note"]].to_string(index=False))
print(f"\nSUM  univariable: n={T4.n_uni.sum()}, events={T4.ev_uni.sum()}")
print(f"SUM  multivariable: n={T4.n_mv.sum()}, events={T4.ev_mv.sum()}")
