import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch

rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 7
rcParams["axes.linewidth"] = .6
rcParams["xtick.major.width"] = .6
rcParams["ytick.major.width"] = .6
rcParams["pdf.fonttype"] = 42
MM = 1 / 25.4

R = pd.read_csv(r"D:\stiffness_r3\results\gsea_hallmark.csv")
R["NES"] = pd.to_numeric(R.NES, errors="coerce")
fc = [c for c in R.columns if "FDR" in c][0]
R[fc] = pd.to_numeric(R[fc], errors="coerce")
R = R[R[fc] < 0.05].sort_values("NES")
top = pd.concat([R.head(8), R.tail(12)])

fig, ax = plt.subplots(figsize=(88 * MM, 78 * MM))
y = np.arange(len(top))
cols = ["#2166AC" if v < 0 else "#B2182B" for v in top.NES]
ax.barh(y, top.NES.values, color=cols, height=.7, lw=0)
ax.axvline(0, color="k", lw=.6)
ax.set_yticks(y)
ax.set_yticklabels([t[:44] for t in top.Term], fontsize=5.4)
ax.set_xlabel("Normalized enrichment score (NES)")
ax.set_xlim(-4.6, 4.8)
ax.legend(handles=[Patch(color="#B2182B", label="enriched in ECM-high"),
                   Patch(color="#2166AC", label="enriched in ECM-low")],
          fontsize=5.4, frameon=False, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=.3)
fig.savefig(r"D:\stiffness_r3\figs\panels\PG1_gsea.pdf", bbox_inches="tight")
print("ok")
