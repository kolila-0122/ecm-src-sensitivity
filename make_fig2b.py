"""重画 Figure 2B：达沙替尼在 GDSC2 与 GDSC1 中的独立复现"""
import os, sys, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

sys.stdout.reconfigure(line_buffering=True)
D = r"D:\stiffness_r3\data\GDSC"
CACHE = r"D:\stiffness_r3\cache"
FIG = r"D:\stiffness_r3\figs\panels"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4

M = pd.read_csv(f"{CACHE}/ccle_ecm_prolif.csv", index_col=0)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def strip_t(s):
    p = str(s).split("_")
    return "_".join(p[:-1]) if len(p) > 1 else str(s)


idx = {}
for i in M.index:
    idx[norm(i)] = i
    idx.setdefault(norm(strip_t(i)), i)


def get(drug, fname):
    d = pd.read_csv(f"{D}/{fname}", low_memory=False)
    d["cell"] = d.CELL_LINE_NAME.map(norm).map(idx)
    g = d[(d.DRUG_NAME == drug)].dropna(subset=["cell"])
    return g.groupby("cell")[["LN_IC50"]].mean().join(M[["ecm"]], how="inner").dropna()


g2 = get("Dasatinib", "GDSC2_dose.csv")
g1 = get("Dasatinib", "GDSC1_dose.csv")
print(f"GDSC2 n={len(g2)}, GDSC1 n={len(g1)}")

fig, axes = plt.subplots(1, 2, figsize=(110 * MM, 55 * MM), sharey=False)
for ax, g, lab, color in [(axes[0], g2, "GDSC2", "#2166AC"),
                          (axes[1], g1, "GDSC1", "#B2182B")]:
    r, p = stats.spearmanr(g.ecm, g.LN_IC50)
    ax.scatter(g.ecm, g.LN_IC50, s=7, alpha=.35, c=color, edgecolor="none", rasterized=True)
    b = np.polyfit(g.ecm, g.LN_IC50, 1)
    xx = np.linspace(g.ecm.quantile(.02), g.ecm.quantile(.98), 50)
    ax.plot(xx, np.polyval(b, xx), color=color, lw=1.3)
    ax.set_xlabel("ECM module score (z)")
    ax.set_ylabel("ln IC$_{50}$, dasatinib")
    ax.set_title(f"{lab}  (n = {len(g)})", fontsize=7)
    ax.text(.04, .06, f"$\\rho$ = {r:.2f}\n$P$ = {p:.1e}", transform=ax.transAxes,
            fontsize=6.4, color=color, va="bottom")
    ax.text(.96, .94, "lower IC$_{50}$ = more sensitive", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.0, color="#666", style="italic")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

fig.tight_layout(pad=.5)
fig.savefig(f"{FIG}/PD2b_dasatinib_two_resources.pdf", bbox_inches="tight")
plt.close(fig)
print("saved PD2b")
