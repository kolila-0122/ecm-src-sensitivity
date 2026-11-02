import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4
FIG = r"D:\stiffness_r3\figs\panels"
OUT = r"D:\stiffness_r3\results"
C1 = "#B2182B"

R = pd.read_csv(f"{OUT}/ecm_survival_multivar.csv")
R = R.iloc[::-1].reset_index(drop=True)   # 倒序，使第一个在最上

n = len(R)
fig, ax = plt.subplots(figsize=(110 * MM, 75 * MM))

y = np.arange(n)
# 校正后 HR 为主
hr, lo, hi = R.HR_adj.values, R.lo_adj.values, R.hi_adj.values
ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="s", color=C1,
            ms=4, lw=.9, capsize=2, capthick=.7, label="Adjusted")

# 单因素（灰）
hr1, lo1, hi1 = R.HR_uni.values, R.lo_uni.values, R.hi_uni.values
ax.errorbar(hr1, y + .22, xerr=[hr1 - lo1, hi1 - hr1], fmt="o", color="#888",
            ms=2.6, lw=.7, capsize=1.6, capthick=.5, alpha=.75, label="Unadjusted")

# Meta 合并
lg = np.log(hr); se = (np.log(hi) - np.log(lo)) / (2 * 1.96)
w = 1 / se ** 2
p = (w * lg).sum() / w.sum(); ps = np.sqrt(1 / w.sum())
mhr, mlo, mhi = np.exp(p), np.exp(p - 1.96 * ps), np.exp(p + 1.96 * ps)
ym = -1.0
ax.errorbar([mhr], [ym], xerr=[[mhr - mlo], [mhi - mhr]], fmt="D", color="#1A1A1A",
            ms=5, lw=1.1, capsize=2.4, capthick=.8)

ax.axvline(1, color="k", lw=.6, ls="--")
ax.set_yticks(list(y) + [ym])
ax.set_yticklabels(list(R.study) + ["Pooled"], fontsize=6.4)
ax.set_xscale("log")
ax.set_xlim(.78, 1.55)
ax.set_xticks([.8, 0.9, 1.0, 1.1, 1.3, 1.5])
ax.set_xticklabels(["0.8", "0.9", "1.0", "1.1", "1.3", "1.5"], fontsize=6)
ax.set_xlabel("Hazard ratio per SD of ECM score (95% CI)")
ax.set_ylim(-2.2, n - .3)

# 右侧标注 n 与 p
for i, r in R.iterrows():
    ax.text(1.52, i, f"{int(r.n_mv)}", ha="right", va="center", fontsize=5.4, color="#555")
ax.text(1.52, n - .6, "n", ha="right", fontsize=5.6, fontweight="bold", color="#333")
ax.text(1.52, ym, f"{int(R.n_mv.sum())}", ha="right", va="center", fontsize=5.6,
        fontweight="bold", color="#333")

p_pool = 2 * (1 - stats.norm.cdf(abs(p / ps)))
ax.text(0.795, ym,
        f"{mhr:.2f} ({mlo:.2f}\u2013{mhi:.2f})",
        fontsize=6.0, color="#1A1A1A", ha="left", va="center", fontweight="bold")
ax.text(0.795, ym - .55, f"$p$ = {p_pool:.3f}, $I^2$ = 0%",
        fontsize=5.6, color="#555", ha="left", va="center")

ax.legend(fontsize=5.6, frameon=False, loc="upper left", bbox_to_anchor=(.01, .99))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PE1_forest.pdf", bbox_inches="tight")
plt.close(fig)
print("saved forest plot")
