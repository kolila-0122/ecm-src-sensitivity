"""验证：分层癌种 + 控制增殖评分"""
import os, sys, json, re
import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
D = r"D:\stiffness_r3\data\GDSC"
OUT = r"D:\stiffness_r3\results"
CACHE = r"D:\stiffness_r3\cache"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
PROLIF = ["MKI67", "PCNA", "CCNA2", "CCNB1", "CDK1", "TOP2A", "BUB1", "AURKA",
          "PLK1", "CDC20", "RRM2", "TYMS"]
KEY_DRUGS = ["Dasatinib", "LGK974", "Oxaliplatin", "Trametinib", "Vorinostat",
             "Sorafenib", "Staurosporine", "IWP-2", "Wnt-C59", "Refametinib",
             "Fulvestrant", "Entospletinib"]

# ---- CCLE 取 ECM + 增殖 ----
fc = os.path.join(CACHE, "ccle_ecm_prolif.csv")
if os.path.exists(fc):
    M = pd.read_csv(fc, index_col=0)
    print(f"载入 {M.shape}")
else:
    genes = ECM + PROLIF
    e2s = {}
    for g in genes:
        r = S.get(f"{B}/genes/{g}", timeout=20)
        if r.status_code == 200:
            e2s[g] = r.json()["entrezGeneId"]
    sym = {v: k for k, v in e2s.items()}
    st = "ccle_broad_2019"
    pid = "ccle_broad_2019_rna_seq_mrna_median_all_sample_Zscores"
    body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{st}_all"}
    rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                data=json.dumps(body), timeout=180)
    df = pd.DataFrame(rr.json())
    X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
    X.index = [sym.get(i, str(i)) for i in X.index]
    eg = [g for g in ECM if g in X.index]
    pg = [g for g in PROLIF if g in X.index]
    print(f"ECM {len(eg)}/{len(ECM)}, 增殖 {len(pg)}/{len(PROLIF)}")
    M = pd.DataFrame({"ecm": X.loc[eg].mean(0), "prolif": X.loc[pg].mean(0)})
    M.index.name = "cell"
    M.to_csv(fc)
    print(f"[saved] {M.shape}")

r_ep = stats.spearmanr(M.ecm, M.prolif)
print(f"\nECM 与增殖评分相关: rho={r_ep[0]:+.3f} (p={r_ep[1]:.1e})")
print("→ 若强正相关，则药物敏感性差异可能由增殖速率驱动\n")

# ---- GDSC ----
dose = pd.read_csv(os.path.join(D, "GDSC2_dose.csv"), low_memory=False)


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def strip_tissue(s):
    parts = str(s).split("_")
    return "_".join(parts[:-1]) if len(parts) > 1 else str(s)


idx = {}
for i in M.index:
    idx[norm(i)] = i
    idx.setdefault(norm(strip_tissue(i)), i)

dose["key"] = dose.CELL_LINE_NAME.map(norm)
dose["cell"] = dose.key.map(idx)
m = dose.dropna(subset=["cell"]).merge(M.reset_index(), on="cell", how="inner")
print(f"匹配 {m.cell.nunique()} 细胞系, {len(m)} 记录")

# ---- 1. 偏相关（控制增殖）----
print("\n=== 控制增殖后的偏相关 ===")
rows = []
for drug in KEY_DRUGS:
    g = m[m.DRUG_NAME == drug]
    if g.empty:
        continue
    per = g.groupby("cell")[["LN_IC50", "ecm", "prolif"]].mean().dropna()
    if len(per) < 40:
        continue
    r_raw, p_raw = stats.spearmanr(per.ecm, per.LN_IC50)
    # 偏相关（秩残差法）
    rk_e = stats.rankdata(per.ecm)
    rk_p = stats.rankdata(per.prolif)
    rk_i = stats.rankdata(per.LN_IC50)
    # ecm 残差
    b1 = np.polyfit(rk_p, rk_e, 1); res_e = rk_e - np.polyval(b1, rk_p)
    b2 = np.polyfit(rk_p, rk_i, 1); res_i = rk_i - np.polyval(b2, rk_p)
    r_p, p_p = stats.pearsonr(res_e, res_i)
    rows.append({"drug": drug, "n": len(per), "rho_raw": r_raw, "p_raw": p_raw,
                 "rho_partial": r_p, "p_partial": p_p,
                 "attenuation": 1 - abs(r_p) / abs(r_raw) if abs(r_raw) > .01 else np.nan})
    print(f"  {drug:16s} n={len(per):3d}  raw rho={r_raw:+.3f} (p={p_raw:.1e})  "
          f"→ 偏相关={r_p:+.3f} (p={p_p:.1e})  衰减 {rows[-1]['attenuation']*100:.0f}%")

# ---- 2. 按癌种分层 ----
print("\n=== 按癌种分层（达沙替尼 & 奥沙利铂）===")
for drug in ["Dasatinib", "Oxaliplatin", "LGK974"]:
    g = m[m.DRUG_NAME == drug]
    print(f"\n{drug}:")
    sub = []
    for tc in g.TCGA_DESC.dropna().unique():
        gg = g[g.TCGA_DESC == tc]
        per = gg.groupby("cell")[["LN_IC50", "ecm"]].mean().dropna()
        if len(per) < 15:
            continue
        r, p = stats.spearmanr(per.ecm, per.LN_IC50)
        sub.append((tc, len(per), r, p))
    sub.sort(key=lambda x: x[2])
    for tc, n, r, p in sub:
        flag = "✓" if p < .05 else " "
        print(f"   {flag} {tc:12s} n={n:3d}  rho={r:+.3f}  p={p:.3f}")
    rs = [x[2] for x in sub]
    if len(rs) > 2:
        print(f"   方向一致: {sum(1 for r in rs if np.sign(r) == np.sign(np.median(rs)))}/{len(rs)}；"
              f"中位 rho={np.median(rs):+.3f}")

pd.DataFrame(rows).to_csv(f"{OUT}/ecm_drug_partial.csv", index=False, encoding="utf-8-sig")
print("\n[saved]")
