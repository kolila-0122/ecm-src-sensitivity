"""补充检验：C-index 比较 + 交互项似然比检验"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

sys.stdout.reconfigure(line_buffering=True)
OUT = r"D:\stiffness_r3\results"
A = pd.read_csv(f"{OUT}/ecm_imm_all.csv")
print(f"合并数据 {A.shape}，{A.study.nunique()} 队列\n")

rows = []
for st in sorted(A.study.unique()):
    T = A[A.study == st].copy()
    d = T[["OS_MONTHS", "event", "ecm_z", "imm_z"]].dropna()
    if len(d) < 80 or d.event.sum() < 15:
        print(f"{st}: 样本不足 n={len(d)}"); continue
    d = d.copy()
    d["inter"] = d.ecm_z * d.imm_z
    d["joint"] = d.ecm_z - d.imm_z     # 联合评分：ECM 高、免疫低 → 分高

    def fit(cols, name):
        c = CoxPHFitter()
        c.fit(d[cols + ["OS_MONTHS", "event"]], duration_col="OS_MONTHS", event_col="event")
        lp = c.predict_log_partial_hazard(d[cols]).values
        ci = concordance_index(d.OS_MONTHS, -lp, d.event)
        return c, ci

    # 四个模型
    c_ecm, ci_ecm = fit(["ecm_z"], "ecm")
    c_imm, ci_imm = fit(["imm_z"], "imm")
    c_both, ci_both = fit(["ecm_z", "imm_z"], "both")
    c_int, ci_int = fit(["ecm_z", "imm_z", "inter"], "inter")
    c_joint, ci_joint = fit(["joint"], "joint")

    # 似然比检验：both vs both+inter
    lrt = 2 * (c_int.log_likelihood_ - c_both.log_likelihood_)
    p_lrt = 1 - stats.chi2.cdf(max(lrt, 0), df=1)

    # 联合评分 vs 两项相加
    d_joint_hr = float(np.exp(c_joint.params_["joint"]))
    d_joint_p = float(c_joint.summary.loc["joint", "p"])
    ci_j = c_joint.confidence_intervals_.loc["joint"].values

    print(f"{st:12s} n={len(d):4d}  "
          f"C-index: ECM {ci_ecm:.3f} | 免疫 {ci_imm:.3f} | 两者 {ci_both:.3f} | "
          f"联合评分 {ci_joint:.3f} | +交互 {ci_int:.3f}")
    print(f"             LRT(交互) χ²={lrt:.2f}, df=1, p={p_lrt:.3f}   "
          f"联合评分 HR={d_joint_hr:.2f} ({np.exp(ci_j[0]):.2f}-{np.exp(ci_j[1]):.2f}) p={d_joint_p:.1e}")

    rows.append({"study": st, "n": len(d),
                 "c_ecm": ci_ecm, "c_imm": ci_imm, "c_both": ci_both,
                 "c_joint": ci_joint, "c_inter": ci_int,
                 "lrt_chi2": lrt, "lrt_p": p_lrt,
                 "hr_joint": d_joint_hr,
                 "lo_joint": float(np.exp(ci_j[0])), "hi_joint": float(np.exp(ci_j[1])),
                 "p_joint": d_joint_p})

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/ecm_imm_supptests.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 92)
print(R.round(3).to_string(index=False))

print("\n=== 汇总 ===")
print(f"C-index 中位：ECM {R.c_ecm.median():.3f} | 免疫 {R.c_imm.median():.3f} | "
      f"两者 {R.c_both.median():.3f} | 联合评分 {R.c_joint.median():.3f}")
print(f"C-index 提升（联合评分 − ECM）：中位 {R.c_joint.median() - R.c_ecm.median():+.3f}")
print(f"C-index 提升（联合评分 − 免疫）：中位 {R.c_joint.median() - R.c_imm.median():+.3f}")
print(f"加入交互项后 C-index 变化：中位 {(R.c_inter - R.c_both).median():+.3f}")
print(f"LRT 交互项显著（p<0.05）: {int((R.lrt_p < .05).sum())}/{len(R)}")

# 联合评分 meta
lg = np.log(R.hr_joint.values)
se = (np.log(R.hi_joint.values) - np.log(R.lo_joint.values)) / (2 * 1.96)
w = 1 / se ** 2
pl = (w * lg).sum() / w.sum(); ps = np.sqrt(1 / w.sum())
Q = float((w * (lg - pl) ** 2).sum())
I2 = max(0, (Q - (len(R) - 1)) / Q) * 100 if Q > 0 else 0
from scipy.stats import norm
print(f"\n联合评分（ECM − 免疫）Meta：HR = {np.exp(pl):.2f} "
      f"({np.exp(pl - 1.96 * ps):.2f}-{np.exp(pl + 1.96 * ps):.2f}), "
      f"p = {2 * (1 - norm.cdf(abs(pl / ps))):.1e}, I² = {I2:.1f}%")
