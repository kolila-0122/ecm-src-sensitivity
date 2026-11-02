"""
检验 SL0101 的 AURKB/慢增殖解释是否成立
思路：如果"慢增殖 → 抗有丝分裂药耐药"成立，
      那么所有 AURK 抑制剂、微管抑制剂都应在 ECM 高时耐药。
      同时看 73 个耐药化合物的靶点是否系统性偏向某几类。
"""
import os, sys, re
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
SRC = r"D:\workbuddy内容\A1_力学代理"
OUT = r"D:\stiffness_r3\results"

S1 = pd.read_csv(os.path.join(SRC, "TableS1_drug_sensitivity_full.csv"), comment="#")
S1["is_sig"] = S1["FDR (sci)"].apply(
    lambda x: float(str(x).replace("<", "")) if pd.notna(x) else 9) < 0.05
S1["direction"] = np.where(S1.rho > 0, "resistant",
                    np.where(S1.rho < 0, "sensitive", "none"))
S1.loc[~S1.is_sig, "direction"] = "ns"

print(f"总化合物 {len(S1)}，显著 {S1.is_sig.sum()}")
print(S1.direction.value_counts().to_string())

# ---------- 按靶点归类 ----------
def classify(row):
    t = str(row.get("Target", "")) + " " + str(row.get("Pathway", ""))
    t = t.upper()
    if re.search(r"AURK|AURORA", t) or "SL0101" in str(row.Drug).upper():
        return "AURK/mitotic"
    if re.search(r"TUBULIN|MICROTUBULE|VINCA|TAXANE", t):
        return "Tubulin"
    if re.search(r"CDK", t):
        return "CDK"
    if re.search(r"MEK|MAPK|ERK", t):
        return "MEK/MAPK"
    if re.search(r"PI3K|AKT|MTOR", t):
        return "PI3K/AKT/mTOR"
    if re.search(r"EGFR|ERBB|HER2", t):
        return "EGFR/ERBB"
    if re.search(r"WNT|PORCN|TNKS", t):
        return "WNT"
    if re.search(r"SRC|LYN|FYN", t):
        return "Src family"
    if re.search(r"HDAC", t):
        return "HDAC"
    if re.search(r"BRAF", t):
        return "BRAF"
    if re.search(r"IGF1R|INSR", t):
        return "IGF1R"
    if re.search(r"BCL| apoptosis", t):
        return "Apoptosis"
    return "Other"

S1["class"] = S1.apply(classify, axis=1)

print("\n=== 各靶点类别的方向分布 ===")
tab = pd.crosstab(S1["class"], S1["direction"])
tab["n_total"] = tab.sum(1)
tab = tab.sort_values("resistant", ascending=False)
print(tab.to_string())

# ---------- 关键检验：AURK/有丝分裂类 ----------
print("\n" + "=" * 62)
print("检验 1：AURK / 有丝分裂类化合物")
print("=" * 62)
aurk = S1[S1["class"] == "AURK/mitotic"]
print(f"n = {len(aurk)}")
if len(aurk):
    print(aurk[["Drug", "Target", "n cell lines", "rho", "p (sci)",
                "FDR (sci)", "direction"]].to_string(index=False))

# ---------- 检查 1：所有 AURK 抑制剂是否都耐药 ----------
print("\n" + "=" * 62)
print("检验 2：微管类（Tubulin）化合物")
print("=" * 62)
tub = S1[S1["class"] == "Tubulin"]
print(f"n = {len(tub)}")
if len(tub):
    print(tub[["Drug", "Target", "n cell lines", "rho", "p (sci)",
               "FDR (sci)", "direction"]].sort_values("rho", ascending=False).to_string(index=False))

# ---------- 检验 3：靶点类别是否过富集 ----------
print("\n" + "=" * 62)
print("检验 3：耐药组 vs 非显著组 的靶点类别富集（Fisher）")
print("=" * 62)
res_tab = []
for cls in S1["class"].unique():
    a = ((S1["class"] == cls) & (S1.direction == "resistant")).sum()
    b = ((S1["class"] == cls) & (S1.direction != "resistant")).sum()
    c = ((S1["class"] != cls) & (S1.direction == "resistant")).sum()
    d = ((S1["class"] != cls) & (S1.direction != "resistant")).sum()
    if a + b < 5:
        continue
    odds, p = stats.fisher_exact([[a, b], [c, d]])
    res_tab.append({"Class": cls, "resistant": a, "other": b,
                    "OR": round(float(odds), 2), "p": f"{p:.2e}"})
R = pd.DataFrame(res_tab).sort_values("p")
print(R.to_string(index=False))

# ---------- 检验 4：慢增殖假说的可检验预测 ----------
print("\n" + "=" * 62)
print("检验 4：慢增殖假说要求「所有抗有丝分裂药都耐药」")
print("=" * 62)
mitotic = S1[S1["class"].isin(["AURK/mitotic", "Tubulin"])]
print(f"抗有丝分裂类共 {len(mitotic)} 个：")
for _, r in mitotic.sort_values("rho", ascending=False).iterrows():
    print(f"  {r.Drug:22s} rho={r.rho:+.3f}  p={str(r['p (sci)']):>10s}  {r.direction}")

n_res = (mitotic.direction == "resistant").sum()
n_total = len(mitotic)
print(f"\n耐药 {n_res}/{n_total}")
if n_total and n_res / n_total > 0.8:
    print("→ 绝大多数抗有丝分裂药耐药，与慢增殖假说一致")
else:
    print("→ 并非全部耐药，慢增殖假说不能单独解释")

# ---------- 保存 ----------
R.to_csv(os.path.join(SRC, "TableS11_target_class_enrichment.csv"),
         index=False, encoding="utf-8")
S1[["Drug", "Target", "class", "n cell lines", "rho", "p (sci)", "FDR (sci)", "direction"]].to_csv(
    os.path.join(OUT, "drug_classes.csv"), index=False)
print("\n[saved TableS11_target_class_enrichment.csv]")
