"""第三轮整改：修复 S2/S3 的 p=0.0、BOM、衰减列定义、GSEA 数量、稳健性表"""
import os, sys
import numpy as np
import pandas as pd

SRC = r"D:\workbuddy内容\A1_力学代理"
OUT = r"D:\stiffness_r3\results"

def sci(x):
    if pd.isna(x):
        return ""
    if x == 0:
        return "<1e-300"
    return f"{x:.2e}"

def strip_bom(path):
    """去掉 UTF-8 BOM"""
    with open(path, "rb") as f:
        raw = f.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        with open(path, "wb") as f:
            f.write(raw[3:])
        return True
    return False

# ============ M6 + M9 + M10: Table S2 ============
p2 = os.path.join(SRC, "TableS2_partial_correlation.csv")
T2 = pd.read_csv(p2)
print("S2 原列:", list(T2.columns))

# 重算精确的 attenuation 与科学计数法
if {"rho (raw)", "rho (partial, adj. proliferation)"}.issubset(T2.columns):
    raw = T2["rho (raw)"]
    par = T2["rho (partial, adj. proliferation)"]
    T2["Attenuation %"] = ((1 - par.abs() / raw.abs()) * 100).round(2)
    T2["p (raw, sci)"] = T2["p (raw)"].map(sci)
    T2["p (partial, sci)"] = T2["p (partial)"].map(sci)

cols = ["Drug", "n", "rho (raw)", "p (raw, sci)",
        "rho (partial, adj. proliferation)", "p (partial, sci)", "Attenuation %"]
cols = [c for c in cols if c in T2.columns]
T2 = T2[cols]
# 表头加公式说明
T2.to_csv(p2, index=False, encoding="utf-8")
# 追加说明行
with open(p2, "a", encoding="utf-8") as f:
    f.write("\n# Attenuation % = (1 - |rho_partial| / |rho_raw|) x 100; "
            "positive = association weakened after adjusting for proliferation, "
            "negative = slightly strengthened\n")
print(f"S2 rebuilt: {len(T2)} rows")
print(T2.round(3).to_string(index=False))

# ============ M6 + M9: Table S3 ============
p3 = os.path.join(SRC, "TableS3_by_cancer_type.csv")
T3 = pd.read_csv(p3)
print("\nS3 原列:", list(T3.columns))
if "p" in T3.columns:
    T3["p (sci)"] = T3["p"].map(sci)
    T3["significant (p<0.05)"] = np.where(T3["p"] < 0.05, "yes", "no")
    T3 = T3[["Drug", "TCGA type", "n", "Spearman rho", "p (sci)", "significant (p<0.05)"]]
T3.to_csv(p3, index=False, encoding="utf-8")
print(f"S3 rebuilt: {len(T3)} rows")

# ============ M9: 去 BOM（所有表） ============
print("\n=== BOM 检查 ===")
for f in sorted(os.listdir(SRC)):
    if f.startswith("TableS") and f.endswith(".csv"):
        path = os.path.join(SRC, f)
        had = strip_bom(path)
        print(f"  {f:42s} {'移除 BOM' if had else '无 BOM'}")

# ============ M8: Table S7 稳健性 ============
p7 = os.path.join(OUT, "robustness_ecm_subsets.csv")
if os.path.exists(p7):
    R = pd.read_csv(p7)
    R = R.round(4)
    R.columns = ["Gene subset", "k cohorts", "HR per SD", "95% CI lower", "95% CI upper", "p", "I2 (%)"]
    R.to_csv(os.path.join(SRC, "TableS7_robustness_gene_subsets.csv"),
             index=False, encoding="utf-8")
    print("\n=== Table S7 稳健性 ===")
    print(R.to_string(index=False))
else:
    print("\n[警告] 缺 robustness_ecm_subsets.csv")

# ============ M7: GSEA 数量核实 ============
g = os.path.join(OUT, "gsea_hallmark.csv")
G = pd.read_csv(g)
G["NES"] = pd.to_numeric(G.NES, errors="coerce")
fc = [c for c in G.columns if "FDR" in c][0]
G[fc] = pd.to_numeric(G[fc], errors="coerce")
n_sig = int((G[fc] < 0.05).sum())
print(f"\n=== GSEA 核实 ===")
print(f"  总基因集 {len(G)}")
print(f"  FDR < 0.05: {n_sig}")
print(f"  图中绘制: 20（取 |NES| 最大的 20 个）")
G_sorted = G.reindex(G.NES.abs().sort_values(ascending=False).index)
print(f"  |NES| 排序前 20 的 FDR 全 <0.05: {(G_sorted.head(20)[fc] < 0.05).all()}")
