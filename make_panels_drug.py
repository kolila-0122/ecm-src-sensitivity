"""药敏结果图：火山图 + 达沙替尼散点 + 癌种分层"""
import os, sys, re, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

sys.stdout.reconfigure(line_buffering=True)
OUT = r"D:\stiffness_r3\results"
FIG = r"D:\stiffness_r3\figs\panels"
CACHE = r"D:\stiffness_r3\cache"
D = r"D:\stiffness_r3\data\GDSC"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4

R = pd.read_csv(f"{OUT}/ecm_drug_sensitivity.csv")
R["neglog"] = -np.log10(R.fdr_ic50.clip(lower=1e-16))
print(f"药物 {len(R)}，FDR<0.05 的 {int((R.fdr_ic50<.05).sum())}")

LAB = {"Dasatinib": "Dasatinib (Src)", "LGK974": "LGK974 (WNT)",
       "Oxaliplatin": "Oxaliplatin", "Trametinib": "Trametinib (MEK)",
       "Vorinostat": "Vorinostat (HDAC)", "Sorafenib": "Sorafenib",
       "IWP-2": "IWP-2 (WNT)", "Wnt-C59": "Wnt-C59 (WNT)",
       "Staurosporine": "Staurosporine", "Refametinib": "Refametinib (MEK)",
       "Fulvestrant": "Fulvestrant", "TAF1_5496": "TAF1i"}

# ---------- 图 1：火山图 ----------
fig, ax = plt.subplots(figsize=(85 * MM, 68 * MM))
sig = R.fdr_ic50 < .05
neg = sig & (R.rho_ic50 < 0)
pos = sig & (R.rho_ic50 > 0)
ns = ~sig
ax.scatter(R.rho_ic50[ns], R.neglog[ns], s=6, c="#CCCCCC", alpha=.55, edgecolor="none", rasterized=True)
ax.scatter(R.rho_ic50[pos], R.neglog[pos], s=12, c="#B2182B", alpha=.75, edgecolor="none", rasterized=True)
ax.scatter(R.rho_ic50[neg], R.neglog[neg], s=12, c="#2166AC", alpha=.75, edgecolor="none", rasterized=True)
ax.axhline(-np.log10(.05), color="k", ls="--", lw=.6)
ax.axvline(0, color="k", lw=.5)

# 标注关键药物
# 只标 5 个关键药物，手动错开
# 只标在两个资源中都复现的药物（GDSC2 + GDSC1）
PICK = {"Dasatinib": (-26, 16), "Trametinib": (-46, -10),
        "Refametinib": (-52, -4), "Vorinostat": (8, 6)}
for _, r in R[R.drug.isin(PICK)].iterrows():
    short = LAB[r.drug].split(" (")[0]
    dx, dy = PICK[r.drug]
    c = "#2166AC" if r.rho_ic50 < 0 else "#B2182B"
    ax.annotate(short, xy=(r.rho_ic50, r.neglog), xytext=(dx, dy),
                textcoords="offset points", fontsize=5.2, color=c,
                arrowprops=dict(arrowstyle="-", color=c, lw=.4))
ax.set_xlabel("Spearman $\\rho$ (ECM score vs ln IC$_{50}$)")
ax.set_ylabel("$-\\log_{10}$ FDR")
ax.text(.03, .96, f"{int(sig.sum())}/{len(R)} drugs at FDR < 0.05",
        transform=ax.transAxes, va="top", fontsize=6.2)
ax.text(.03, .88, "positive: ECM high → resistant", transform=ax.transAxes,
        va="top", fontsize=5.4, color="#B2182B")
ax.text(.03, .82, "negative: ECM high → sensitive", transform=ax.transAxes,
        va="top", fontsize=5.4, color="#2166AC")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PD1_volcano.pdf", bbox_inches="tight")
plt.close(fig)
print("volcano saved")

# ---------- 图 2：达沙替尼散点 ----------
M = pd.read_csv(f"{CACHE}/ccle_ecm_prolif.csv", index_col=0)
dose = pd.read_csv(f"{D}/GDSC2_dose.csv", low_memory=False)


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

fig, ax = plt.subplots(figsize=(72 * MM, 62 * MM))
for drug, c in [("Dasatinib", "#2166AC"), ("LGK974", "#B2182B")]:
    g = dose[dose.DRUG_NAME == drug].dropna(subset=["cell"])
    per = g.groupby("cell")[["LN_IC50"]].mean().join(M[["ecm"]], how="inner").dropna()
    r, p = stats.spearmanr(per.ecm, per.LN_IC50)
    ax.scatter(per.ecm, per.LN_IC50, s=5, alpha=.32, c=c, edgecolor="none",
               rasterized=True, label=f"{drug} ($\\rho$={r:+.2f})")
    # 拟合线
    b = np.polyfit(per.ecm, per.LN_IC50, 1)
    xx = np.linspace(per.ecm.quantile(.02), per.ecm.quantile(.98), 50)
    ax.plot(xx, np.polyval(b, xx), color=c, lw=1.2)

ax.set_xlabel("ECM module score (z)")
ax.set_ylabel("ln IC$_{50}$")
ax.legend(fontsize=5.4, frameon=False, loc="upper left")
ax.text(.97, .04, "low IC$_{50}$ = sensitive", transform=ax.transAxes,
        ha="right", fontsize=5.2, color="#666", style="italic")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PD2_dasatinib_scatter.pdf", bbox_inches="tight")
plt.close(fig)
print("scatter saved")

# ---------- 图 3：癌种分层 ----------
rows = []
for drug in ["Dasatinib", "Oxaliplatin", "LGK974"]:
    g = dose[dose.DRUG_NAME == drug].dropna(subset=["cell"])
    for tc in g.TCGA_DESC.dropna().unique():
        gg = g[g.TCGA_DESC == tc]
        per = gg.groupby("cell")[["LN_IC50"]].mean().join(M[["ecm"]], how="inner").dropna()
        if len(per) < 15:
            continue
        r, p = stats.spearmanr(per.ecm, per.LN_IC50)
        rows.append({"drug": drug, "tcga": tc, "n": len(per), "rho": r, "p": p})
S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/ecm_drug_bycancer.csv", index=False, encoding="utf-8-sig")

tcs = sorted(S.tcga.unique())
drugs = ["Dasatinib", "Oxaliplatin", "LGK974"]
pivot = S.pivot_table(index="tcga", columns="drug", values="rho")
pivot = pivot.reindex(tcs)

fig, ax = plt.subplots(figsize=(72 * MM, 78 * MM))
im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-.7, vmax=.7, aspect="auto")
ax.set_xticks(range(len(drugs)))
ax.set_xticklabels([d for d in drugs], fontsize=6, rotation=30, ha="right")
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=5.6)
for i in range(len(pivot.index)):
    for j in range(len(drugs)):
        v = pivot.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=4.8,
                    color="white" if abs(v) > .4 else "#222")
cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03)
cb.set_label("Spearman $\\rho$", fontsize=5.6)
cb.ax.tick_params(labelsize=5)
ax.set_title("ECM–drug correlation by cancer type", fontsize=6.6, pad=4)
fig.tight_layout(pad=.3)
fig.savefig(f"{FIG}/PD3_bycancer_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("heatmap saved")
