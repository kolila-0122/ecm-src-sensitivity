"""M1 修复：用 t 近似重算全部 p/FDR，消除 <1e-300 硬编码"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
SRC = r"D:\workbuddy内容\A1_力学代理"


def spearman_p(rho, n):
    """标准 t 近似，不依赖 scipy.spearmanr（后者在大 n 时下溢）"""
    if not np.isfinite(rho) or n < 3:
        return np.nan
    if abs(rho) >= 1:
        return 0.0
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    return float(2 * stats.t.sf(abs(t), df=n - 2))


def fmt_p(p):
    """科学计数法；只有真正为 0 才写 <1e-300"""
    if p is None or not np.isfinite(p):
        return ""
    if p == 0:
        return "<1e-300"
    return f"{p:.2e}"


# ============ Table S1 ============
p1 = os.path.join(SRC, "TableS1_drug_sensitivity_full.csv")
T1 = pd.read_csv(p1)
print("S1 列:", list(T1.columns))

ncol = [c for c in T1.columns if "cell" in c.lower() or c == "n"][0]
rho_col = "rho"
T1["p (sci)"] = [fmt_p(spearman_p(r, n)) for r, n in zip(T1[rho_col], T1[ncol])]

# BH-FDR 重算（286 项）
pvals = np.array([spearman_p(r, n) for r, n in zip(T1[rho_col], T1[ncol])])
order = np.argsort(pvals)
m = len(pvals)
fdr = np.empty(m)
prev = 1.0
for rank in range(m - 1, -1, -1):
    i = order[rank]
    val = pvals[i] * m / (rank + 1)
    prev = min(prev, val)
    fdr[i] = min(prev, 1.0)
T1["FDR (sci)"] = [fmt_p(x) for x in fdr]

T1.to_csv(p1, index=False, encoding="utf-8")
print(f"\nS1 重算完成 {len(T1)} 行")
print(f"  仍为 <1e-300 的 p: {(T1['p (sci)']=='<1e-300').sum()} 行")
print(f"  仍为 <1e-300 的 FDR: {(T1['FDR (sci)']=='<1e-300').sum()} 行")

print("\n=== 抽样核对（审稿人点名 5 行）===")
for d in ["Dasatinib", "LGK974", "Oxaliplatin", "Dihydrorotenone", "Wnt-C59"]:
    r = T1[T1.Drug == d]
    if len(r):
        r = r.iloc[0]
        p_calc = spearman_p(float(r[rho_col]), int(r[ncol]))
        print(f"  {d:18s} rho={r[rho_col]:+.4f} n={int(r[ncol]):3d}  "
              f"表报={r['p (sci)']:10s}  复算={p_calc:.2e}  FDR={r['FDR (sci)']}")

# ============ Table S1b ============
sig = T1[T1["FDR (sci)"].apply(lambda s: float(s) if s and not s.startswith("<") else 0) < 0.05]
sig.to_csv(os.path.join(SRC, "TableS1b_significant_87.csv"), index=False, encoding="utf-8")
print(f"\nS1b 重建: {len(sig)} 行（rho>0: {(sig[rho_col]>0).sum()}, rho<0: {(sig[rho_col]<0).sum()}）")

# ============ Table S2 ============
p2 = os.path.join(SRC, "TableS2_partial_correlation.csv")
raw = open(p2, encoding="utf-8").read()
T2 = pd.read_csv(p2, comment="#")
print("\nS2 列:", list(T2.columns))

T2["p (raw, sci)"] = [fmt_p(spearman_p(r, n)) for r, n in zip(T2["rho (raw)"], T2["n"])]
T2["p (partial, sci)"] = [fmt_p(spearman_p(r, n)) for r, n in
                          zip(T2["rho (partial, adj. proliferation)"], T2["n"])]

cols = ["Drug", "n", "rho (raw)", "p (raw, sci)",
        "rho (partial, adj. proliferation)", "p (partial, sci)", "Attenuation %"]
cols = [c for c in cols if c in T2.columns]
T2 = T2[cols]
T2.to_csv(p2, index=False, encoding="utf-8")
with open(p2, "a", encoding="utf-8") as f:
    f.write("\n# p values recomputed with the standard t approximation "
            "t = rho*sqrt((n-2)/(1-rho^2)), two-sided, df = n-2. "
            "FDR by Benjamini-Hochberg across 286 tests in Table S1. "
            "Attenuation % = (1 - |rho_partial|/|rho_raw|) x 100.\n")
print(f"S2 重建: {len(T2)} 行")

print("\n=== S2 抽样 ===")
for d in ["Dasatinib", "Oxaliplatin"]:
    r = T2[T2.Drug == d]
    if len(r):
        r = r.iloc[0]
        print(f"  {d:14s} raw p={r['p (raw, sci)']:10s}  partial p={r['p (partial, sci)']:10s}  "
              f"attenuation={r['Attenuation %']:+.2f}%")

# ============ Table S3 ============
p3 = os.path.join(SRC, "TableS3_by_cancer_type.csv")
T3 = pd.read_csv(p3, comment="#")
print("\nS3 列:", list(T3.columns))
rc = "Spearman rho"
if "p" in T3.columns:
    T3["p (sci)"] = [fmt_p(spearman_p(r, n)) for r, n in zip(T3[rc], T3["n"])]
    T3["significant (p<0.05)"] = np.where(
        [spearman_p(r, n) < 0.05 for r, n in zip(T3[rc], T3["n"])], "yes", "no")
T3.to_csv(p3, index=False, encoding="utf-8")
print(f"S3 重建: {len(T3)} 行，仍为 <1e-300 的: {(T3['p (sci)']=='<1e-300').sum() if 'p (sci)' in T3 else 'n/a'}")
