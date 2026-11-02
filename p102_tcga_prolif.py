"""
M4：在 TCGA 肿瘤组织中验证"ECM 高 = 增殖低"
用 MKI67 mRNA 与经典增殖基因作为指标（回应"无 Ki-67 染色"的批评）
对比 CCLE 细胞系中的同一关系。
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
import requests
from scipy import stats

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(line_buffering=True)

B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
SRC = r"D:\workbuddy内容\A1_力学代理"
OUT = r"D:\stiffness_r3\results"

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}

ECM = ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","LOX","LOXL1","LOXL2",
       "FN1","POSTN","CCN2","ACTA2","TAGLN","TIMP1","SPARC","DCN","LUM","TGFB1","COL11A1","COL8A1","FBN1"]

# 经典增殖标志（含 MKI67）
PROLIF = ["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CDK1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
          "BUB1","BUB1B","AURKA","AURKB","PLK1","CDC20","CENPF","TYMS","RRM2","UBE2C"]

S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})
g = S.get(f"{B}/genes", timeout=180).json()
gmap = {x["hugoGeneSymbol"].upper(): x["entrezGeneId"] for x in g if x.get("hugoGeneSymbol")}
inv = {v: k for k, v in gmap.items()}
need = ECM + PROLIF
ids = [gmap[x] for x in need if x in gmap]
print(f"需要 {len(need)} 个基因，找到 {len(ids)} 个\n")

rows = []
for study, label in STUDIES.items():
    profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
    pid = next(p["molecularProfileId"] for p in profs
               if "zscore" in p["molecularProfileId"].lower()
               and p.get("molecularAlterationType") == "MRNA_EXPRESSION"
               and "mirna" not in p["molecularProfileId"].lower())
    r = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
               data=json.dumps({"entrezGeneIds": ids, "sampleListId": f"{study}_all"}),
               timeout=200)
    X = pd.DataFrame(r.json()).pivot_table(index="entrezGeneId", columns="sampleId",
                                           values="value", aggfunc="mean")
    X.index = [inv.get(i, str(i)) for i in X.index]

    ecm_g = [x for x in ECM if x in X.index]
    pro_g = [x for x in PROLIF if x in X.index]
    ecm_score = X.loc[ecm_g].mean(0)
    prolif_score = X.loc[pro_g].mean(0)

    rho, p = stats.spearmanr(ecm_score, prolif_score)
    mk = X.loc["MKI67"] if "MKI67" in X.index else None
    if mk is not None:
        rho_mk, p_mk = stats.spearmanr(ecm_score, mk.loc[ecm_score.index])
    else:
        rho_mk, p_mk = np.nan, np.nan

    rows.append({
        "Cohort": label, "n": len(ecm_score),
        "ECM vs proliferation": round(float(rho), 3), "p": f"{p:.1e}",
        "ECM vs MKI67": round(float(rho_mk), 3) if np.isfinite(rho_mk) else np.nan,
        "p (MKI67)": f"{p_mk:.1e}" if np.isfinite(p_mk) else "",
    })
    print(f"{label:12s} n={len(ecm_score):4d}  ECM vs 增殖评分 ρ={rho:+.3f} (p={p:.1e})"
          + (f"   ECM vs MKI67 ρ={rho_mk:+.3f} (p={p_mk:.1e})" if np.isfinite(rho_mk) else ""))

T = pd.DataFrame(rows)
print("\n=== 汇总 ===")
print(f"六队列 ECM vs 增殖评分：中位 ρ = {T['ECM vs proliferation'].median():+.3f}，"
      f"范围 {T['ECM vs proliferation'].min():+.3f} 至 {T['ECM vs proliferation'].max():+.3f}")
print(f"六队列 ECM vs MKI67：  中位 ρ = {T['ECM vs MKI67'].median():+.3f}，"
      f"范围 {T['ECM vs MKI67'].min():+.3f} 至 {T['ECM vs MKI67'].max():+.3f}")
n_neg = (T["ECM vs MKI67"] < 0).sum()
print(f"MKI67 负相关的队列：{n_neg}/6")

# 与 CCLE 对比
M = pd.read_csv(os.path.join(CACHE, "ccle_ecm_prolif.csv"), index_col=0)
rho_ccle = stats.spearmanr(M.ecm, M.prolif)[0]
print(f"\nCCLE 细胞系 ECM vs 增殖：ρ = {rho_ccle:+.3f} (n={len(M)})")

T.to_csv(os.path.join(SRC, "TableS9_ecm_proliferation_tcga.csv"),
         index=False, encoding="utf-8")
print("\n[saved TableS9_ecm_proliferation_tcga.csv]")

# 结论句
print("\n" + "=" * 60)
if n_neg >= 5:
    print("结论：TCGA 肿瘤组织中 ECM 与增殖呈负相关 → 支持 'slow-cycling' 可在组织层面成立")
else:
    print("结论：TCGA 中 ECM 与增殖的相关方向不一致 → 不能把 slow-cycling 从细胞系迁移到肿瘤")
print("=" * 60)
