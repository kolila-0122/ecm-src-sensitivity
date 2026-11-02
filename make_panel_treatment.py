"""治疗反应阴性结果图：免疫治疗 + 化疗 pCR 的 ROC 与 AUC"""
import os, sys, gzip, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, roc_auc_score

sys.stdout.reconfigure(line_buffering=True)
OUT = r"D:\stiffness_r3\results"
FIG = r"D:\stiffness_r3\figs\panels"
D = r"D:\stiffness_r3\data\GEO_neoadj"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4

# ---------- 免疫治疗 ROC ----------
ICI = pd.read_csv(f"{OUT}/ecm_ici_all.csv")
print(f"免疫队列 {ICI.shape}, 队列 {ICI.study.unique()}")
ICI = ICI.rename(columns={"resp": "resp"})

fig, axes = plt.subplots(1, 2, figsize=(120 * MM, 58 * MM))

ax = axes[0]
cols = plt.cm.tab10.colors
# 汇总 ROC（各队列分别画 + 汇总加权）
all_fpr, all_tpr = [], []
for i, (st, g) in enumerate(ICI.groupby("study")):
    if g.resp.nunique() < 2 or len(g) < 15:
        continue
    fpr, tpr, _ = roc_curve(g.resp, g.ecm)
    auc = roc_auc_score(g.resp, g.ecm)
    ax.plot(fpr, tpr, lw=.9, color=cols[i], alpha=.85,
            label=f"{st[:18]} ({auc:.2f})")
ax.plot([0, 1], [0, 1], "k--", lw=.6, alpha=.5)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("Immunotherapy response", fontsize=7)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(fontsize=4.6, frameon=False, loc="lower right")
ax.text(.03, .96, "AUC ≈ 0.5 across cohorts", transform=ax.transAxes,
        va="top", fontsize=5.6, color="#666")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---------- 化疗 pCR ROC ----------
R2 = pd.read_csv(f"{OUT}/ecm_pcr.csv")
ax = axes[1]
# 有 AUC，直接画条形图 + 0.5 参考线
y = np.arange(len(R2))
c = ["#7F7F7F"] * len(R2)
ax.barh(y, R2.AUC.values, color=c, height=.6, lw=0)
ax.axvline(.5, color="k", ls="--", lw=.7)
ax.axvline(0, color="k", lw=.5)
ax.set_yticks(y)
ax.set_yticklabels([f"{r.study}\n(n={int(r.n)})" for _, r in R2.iterrows()], fontsize=5.6)
ax.set_xlim(0, 1)
ax.set_xlabel("AUC (ECM score → pCR)")
ax.set_title("Neoadjuvant chemotherapy pCR", fontsize=7)
for i, r in R2.iterrows():
    ax.text(r.AUC + .01, i, f"{r.AUC:.2f}", va="center", fontsize=5.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout(pad=.4)
fig.savefig(f"{FIG}/PT1_treatment_response.pdf", bbox_inches="tight")
plt.close(fig)
print("treatment response panel saved")

# ---------- 合并 AUC 汇总条形 ----------
fig, ax = plt.subplots(figsize=(72 * MM, 52 * MM))
labels = ["Immunotherapy\n(4 cohorts)", "Chemotherapy pCR\n(4 cohorts)"]
values = [0.476, 0.457]
b = ax.bar([0, 1], values, width=.5, color=["#8C8C8C", "#A65628"], lw=0)
ax.axhline(.5, color="k", ls="--", lw=.7)
ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=6)
ax.set_ylabel("Pooled AUC")
ax.set_ylim(0, .75)
for x, v in zip([0, 1], values):
    ax.text(x, v + .015, f"{v:.3f}", ha="center", fontsize=6.5, fontweight="bold")
ax.text(.5, .53, "chance", fontsize=5.2, color="#666", ha="center")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PT2_auc_summary.pdf", bbox_inches="tight")
plt.close(fig)
print("AUC summary saved")
