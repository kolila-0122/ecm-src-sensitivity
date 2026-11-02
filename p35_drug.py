"""ECM 评分与 GDSC 药物敏感性"""
import os, sys, json, time, re
import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
D = r"D:\stiffness_r3\data\GDSC"
OUT = r"D:\stiffness_r3\results"
CACHE = r"D:\stiffness_r3\cache"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

# ---- 1. CCLE ECM 评分 ----
fc = os.path.join(CACHE, "ccle_ecm.csv")
if os.path.exists(fc):
    ecm = pd.read_csv(fc, index_col=0)
    print(f"载入 CCLE ECM 评分 {ecm.shape}")
else:
    e2s = {}
    for g in ECM:
        r = S.get(f"{B}/genes/{g}", timeout=20)
        if r.status_code == 200:
            e2s[g] = r.json()["entrezGeneId"]
    sym = {v: k for k, v in e2s.items()}
    st = "ccle_broad_2019"
    pid = "ccle_broad_2019_rna_seq_mrna_median_all_sample_Zscores"
    body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{st}_all"}
    rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                data=json.dumps(body), timeout=180)
    if rr.status_code != 200:
        raise SystemExit(f"取数失败 {rr.status_code}")
    df = pd.DataFrame(rr.json())
    X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
    X.index = [sym.get(i, str(i)) for i in X.index]
    gs = [g for g in ECM if g in X.index]
    print(f"ECM 基因 {len(gs)}/{len(ECM)}, 细胞系 {X.shape[1]}")
    ecm = pd.DataFrame({"ecm": X.loc[gs].mean(0)})
    ecm.index.name = "cell"
    ecm.to_csv(fc)
    print(f"[saved] {ecm.shape}")

# ---- 2. GDSC 药敏 ----
dose = pd.read_csv(os.path.join(D, "GDSC2_dose.csv"), low_memory=False)
print(f"GDSC2 {dose.shape}，药物 {dose.DRUG_NAME.nunique()}，细胞系 {dose.CELL_LINE_NAME.nunique()}")

# ---- 3. 名称匹配 ----
def norm(s):
    s = str(s).upper()
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s

# CCLE 名形如 143B_BONE；去掉组织后缀
def strip_tissue(s):
    parts = str(s).split("_")
    return "_".join(parts[:-1]) if len(parts) > 1 else str(s)

ccle_idx = {}
for i in ecm.index:
    ccle_idx[norm(i)] = i
    base = strip_tissue(i)
    ccle_idx.setdefault(norm(base), i)
dose["key"] = dose.CELL_LINE_NAME.map(norm)
dose["ccle"] = dose.key.map(ccle_idx)
matched = dose.dropna(subset=["ccle"])
print(f"匹配到 {matched.ccle.nunique()} 个细胞系，{len(matched)} 条药敏记录")

if matched.ccle.nunique() < 50:
    print("匹配率过低，列出未匹配样例：")
    print(list(dose.CELL_LINE_NAME.unique()[:10]))
    print(list(ecm.index[:10]))
    raise SystemExit("匹配失败")

# ---- 4. 每个药物的相关 ----
matched = matched.merge(ecm.reset_index(), left_on="ccle", right_on="cell", how="inner")
rows = []
for drug, g in matched.groupby("DRUG_NAME"):
    # 每个细胞系取平均 lnIC50
    per_cell = g.groupby("cell")[["LN_IC50", "AUC", "ecm"]].mean().dropna()
    if len(per_cell) < 40:
        continue
    r_ic, p_ic = stats.spearmanr(per_cell.ecm, per_cell.LN_IC50)
    r_auc, p_auc = stats.spearmanr(per_cell.ecm, per_cell.AUC)
    rows.append({"drug": drug, "n": len(per_cell),
                 "rho_ic50": r_ic, "p_ic50": p_ic,
                 "rho_auc": r_auc, "p_auc": p_auc,
                 "target": g.PUTATIVE_TARGET.iloc[0],
                 "pathway": g.PATHWAY_NAME.iloc[0]})

R = pd.DataFrame(rows)
R["fdr_ic50"] = stats.false_discovery_control(R.p_ic50) if hasattr(stats, "false_discovery_control") \
    else R.p_ic50 * len(R) / R.p_ic50.rank()
R = R.sort_values("p_ic50")
R.to_csv(f"{OUT}/ecm_drug_sensitivity.csv", index=False, encoding="utf-8-sig")

print(f"\n可分析药物 {len(R)} 个")
print(f"IC50 正相关（ECM 高→耐药）: {int((R.rho_ic50 > 0).sum())}")
print(f"FDR<0.05 的显著药物: {int((R.fdr_ic50 < .05).sum())}")

print("\n=== ECM 高 → 更耐药（rho>0）TOP 15 ===")
top = R[R.rho_ic50 > 0].head(15)
print(top[["drug", "n", "rho_ic50", "p_ic50", "fdr_ic50", "pathway"]].round(4).to_string(index=False))

print("\n=== ECM 高 → 更敏感（rho<0）TOP 8 ===")
neg = R[R.rho_ic50 < 0].sort_values("rho_ic50").head(8)
print(neg[["drug", "n", "rho_ic50", "p_ic50", "pathway"]].round(4).to_string(index=False))
