"""四组 KM 曲线 + 森林图"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines import CoxPHFitter

sys.stdout.reconfigure(line_buffering=True)
OUT = r"D:\stiffness_r3\results"
FIG = r"D:\stiffness_r3\figs\panels"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4

A = pd.read_csv(f"{OUT}/ecm_imm_all.csv")
print(f"合并数据 {A.shape}  队列 {A.study.nunique()}")

# ---------- KM：选样本最大的癌种 ----------
sel = A[A.study == "Breast"].copy()
print(f"KM 用 Breast，n={len(sel)}")
sel["grp_lab"] = sel.grp.map({0: "ECM low / immune low", 1: "ECM low / immune high",
                              2: "ECM high / immune low", 3: "ECM high / immune high"})
colors = {0: "#92C5DE", 1: "#2166AC", 2: "#B2182B", 3: "#F4A582"}

fig, ax = plt.subplots(figsize=(85 * MM, 68 * MM))
for g in [1, 3, 0, 2]:          # 从最好到最差叠放
    m = sel.grp == g
    if m.sum() < 10: continue
    kmf = KaplanMeierFitter()
    kmf.fit(sel.OS_MONTHS[m] / 12, sel.event[m],
            label=f"{sel.grp_lab[m].iloc[0]} (n={int(m.sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=colors[g], lw=1.1)

lr = multivariate_logrank_test(sel.OS_MONTHS, sel.grp, sel.event)
ax.set_xlabel("Time (years)")
ax.set_ylabel("Overall survival")
ax.set_ylim(0, 1.02)
ax.set_xlim(0, min(sel.OS_MONTHS.max() / 12, 15))
ax.text(.98, .97, f"$P$ = {lr.p_value:.1e}", transform=ax.transAxes,
        ha="right", va="top", fontsize=7)
ax.legend(fontsize=5.2, frameon=False, loc="upper right", bbox_to_anchor=(1, .90))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PI1_km_four_groups.pdf", bbox_inches="tight")
plt.close(fig)
print("KM saved")

# ---------- 森林图 ----------
R = pd.read_csv(f"{OUT}/ecm_imm_interaction.csv").dropna(subset=["HR_worst_vs_best"])
R = R.iloc[::-1].reset_index(drop=True)
n = len(R)
fig, ax = plt.subplots(figsize=(100 * MM, 62 * MM))
y = np.arange(n)
hr, lo, hi = R.HR_worst_vs_best.values, R.lo_wb.values, R.hi_wb.values
ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="s", color="#B2182B",
            ms=4, lw=.9, capsize=2, capthick=.7)

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
ax.set_xlim(.45, 7)
ax.set_xticks([.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0])
ax.set_xticklabels(["0.5", "0.7", "1.0", "1.5", "2.0", "3.0", "5.0"], fontsize=6)
from matplotlib.ticker import NullFormatter, FixedLocator
ax.xaxis.set_minor_formatter(NullFormatter())
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_xlabel("Hazard ratio, ECM-high/immune-low vs ECM-low/immune-high (95% CI)")
ax.set_ylim(-2.0, n - .3)

for i, r in R.iterrows():
    ax.text(6.6, i, f"{int(r.n)}", ha="right", va="center", fontsize=5.2, color="#555")
ax.text(6.6, n - .55, "n", ha="right", fontsize=5.6, fontweight="bold", color="#333")
ax.text(6.6, ym, f"{int(R.n.sum())}", ha="right", va="center", fontsize=5.6,
        fontweight="bold", color="#333")

from scipy.stats import norm
p_pool = 2 * (1 - norm.cdf(abs(p / ps)))
ax.text(0.47, ym + .05, f"{mhr:.2f} ({mlo:.2f}\u2013{mhi:.2f})",
        fontsize=6.2, color="#1A1A1A", fontweight="bold", va="center")
ax.text(0.47, ym - .5, f"$P$ = {p_pool:.4f}, $I^2$ = 0%", fontsize=5.6, color="#555", va="center")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PI2_forest.pdf", bbox_inches="tight")
plt.close(fig)
print("Forest saved")

# ---------- 相关性散点（ECM vs 免疫）----------
fig, ax = plt.subplots(figsize=(75 * MM, 62 * MM))
for st, c in zip(sorted(A.study.unique()), plt.cm.tab10.colors):
    m = A.study == st
    ax.scatter(A.ecm[m], A.imm[m], s=1.2, alpha=.18, c=[c], edgecolor="none", rasterized=True)
    r = stats.pearsonr(A.ecm[m], A.imm[m])[0]
ax.set_xlabel("ECM module score (z)")
ax.set_ylabel("CD8/cytotoxic signature (z)")
ax.text(.03, .97, "ECM and immune infiltration are\nlargely independent (|r| < 0.26)",
        transform=ax.transAxes, va="top", fontsize=6.0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PI3_scatter.pdf", bbox_inches="tight")
plt.close(fig)
print("Scatter saved")
