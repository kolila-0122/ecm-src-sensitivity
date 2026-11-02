"""
M3：对全部 73 个耐药化合物做增殖校正
回答审稿意见："73 个耐药是否主要由慢增殖驱动？"

做法：对每个 FDR<0.05 且 ρ>0 的化合物，
      算 raw Spearman ρ(ECM, IC50) 与 partial ρ(ECM, IC50 | proliferation)，
      比较校正前后的衰减幅度。
"""
import os, sys, re, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(line_buffering=True)

CACHE = r"D:\stiffness_r3\cache"
DATA = r"D:\stiffness_r3\data"
OUT = r"D:\stiffness_r3\results"
SRC = r"D:\workbuddy内容\A1_力学代理"


def spearman_p(rho, n):
    if not np.isfinite(rho) or n < 3:
        return np.nan
    if abs(rho) >= 1:
        return 0.0
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    return float(2 * stats.t.sf(abs(t), df=n - 2))


def partial_spearman(x, y, z):
    """秩残差法偏相关"""
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    # 分别对 z 回归，取残差
    def resid(a):
        b = np.polyfit(rz, a, 1)
        return a - np.polyval(b, rz)
    ex, ey = resid(rx), resid(ry)
    return float(stats.pearsonr(ex, ey)[0])


# ---------- 数据 ----------
M = pd.read_csv(os.path.join(CACHE, "ccle_ecm_prolif.csv"), index_col=0)
print(f"CCLE: {M.shape}, 列 {list(M.columns)}")

dose = pd.read_csv(os.path.join(DATA, "GDSC", "GDSC2_dose.csv"), low_memory=False)
print(f"GDSC2 dose: {dose.shape}")


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def strip_tissue(s):
    p = str(s).split("_")
    return "_".join(p[:-1]) if len(p) > 1 else str(s)


idx = {}
for i in M.index:
    idx[norm(i)] = i
    idx.setdefault(norm(strip_tissue(i)), i)
dose["cell"] = dose.CELL_LINE_NAME.map(norm).map(idx)

# ---------- 显著化合物 ----------
S1 = pd.read_csv(os.path.join(SRC, "TableS1_drug_sensitivity_full.csv"), comment="#")
sig = S1[S1["FDR (sci)"].apply(lambda x: float(str(x).replace("<", "")) if pd.notna(x) else 9) < 0.05]
resistant = sig[sig.rho > 0]["Drug"].tolist()
sensitive = sig[sig.rho < 0]["Drug"].tolist()
print(f"\n显著 {len(sig)}：耐药 {len(resistant)}，敏感 {len(sensitive)}")

# ---------- 逐化合物算 ----------
rows = []
for drug in sig.Drug:
    g = dose[dose.DRUG_NAME == drug].dropna(subset=["cell"])
    if not len(g):
        continue
    per = g.groupby("cell")[["LN_IC50"]].mean().join(M[["ecm", "prolif"]], how="inner").dropna()
    if len(per) < 40:
        continue
    r_raw = stats.spearmanr(per.ecm, per.LN_IC50)[0]
    r_par = partial_spearman(per.ecm.values, per.LN_IC50.values, per.prolif.values)
    n = len(per)
    rows.append({
        "Drug": drug,
        "n": n,
        "direction": "resistant" if r_raw > 0 else "sensitive",
        "rho_raw": round(float(r_raw), 4),
        "rho_partial": round(float(r_par), 4),
        "attenuation_%": round(float((1 - abs(r_par) / abs(r_raw)) * 100), 2) if r_raw else np.nan,
        "p_raw": f"{spearman_p(r_raw, n):.2e}",
        "p_partial": f"{spearman_p(r_par, n):.2e}",
    })

P = pd.DataFrame(rows)
P.to_csv(os.path.join(OUT, "all_partial_corr.csv"), index=False)

print("\n=== 按方向汇总 ===")
for d in ["resistant", "sensitive"]:
    sub = P[P.direction == d]
    if not len(sub):
        continue
    print(f"\n{d} (n={len(sub)}):")
    print(f"  |rho_raw|     中位 {sub.rho_raw.abs().median():.3f}，范围 {sub.rho_raw.abs().min():.3f}–{sub.rho_raw.abs().max():.3f}")
    print(f"  |rho_partial| 中位 {sub.rho_partial.abs().median():.3f}，范围 {sub.rho_partial.abs().min():.3f}–{sub.rho_partial.abs().max():.3f}")
    print(f"  衰减 %        中位 {sub['attenuation_%'].median():.1f}，范围 {sub['attenuation_%'].min():.1f} 至 {sub['attenuation_%'].max():.1f}")
    # 校正后仍有 p<0.05 的
    still = sub[sub.p_partial.apply(lambda x: float(x) < 0.05)]
    print(f"  校正后仍 p<0.05: {len(still)}/{len(sub)}")
    # 符号翻转
    flip = sub[np.sign(sub.rho_raw) != np.sign(sub.rho_partial)]
    print(f"  校正后符号翻转: {len(flip)}" + (f" → {flip.Drug.tolist()[:5]}" if len(flip) else ""))

# ---------- 关键对照：衰减幅度是否与 |rho_raw| 相关 ----------
print("\n=== 衰减 vs 效应量 ===")
if len(P) > 10:
    r_att = stats.spearmanr(P.rho_raw.abs(), P["attenuation_%"].abs())
    print(f"  |rho_raw| 与 |attenuation| 的 Spearman: rho={r_att[0]:.3f}, p={r_att[1]:.2e}")
    print(f"  → {'衰减与初始效应量相关，校正可能是回归到均值' if r_att[1] < 0.05 else '无系统性关系'}")

# ---------- 保存给稿件 ----------
P_sorted = P.sort_values("attenuation_%", ascending=False)
P_sorted.to_csv(os.path.join(SRC, "TableS8_proliferation_adjustment_all.csv"),
                index=False, encoding="utf-8")
print(f"\n[saved] TableS8_proliferation_adjustment_all.csv ({len(P)} compounds)")
print("\n衰减最大的 8 个：")
print(P_sorted.head(8)[["Drug", "n", "direction", "rho_raw", "rho_partial",
                         "attenuation_%", "p_partial"]].to_string(index=False))
print("\n衰减最小（=最不受增殖影响）的 8 个：")
print(P_sorted.tail(8)[["Drug", "n", "direction", "rho_raw", "rho_partial",
                         "attenuation_%", "p_partial"]].to_string(index=False))
