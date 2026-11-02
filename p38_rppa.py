"""ECM 高 vs Src 通路活性（CCLE RPPA）"""
import os, sys, json
import numpy as np
import pandas as pd
import requests
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
OUT = r"D:\stiffness_r3\results"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

# 关注的通路蛋白（RPPA 抗体名可能带磷酸化位点）
TARGETS = ["SRC", "SRC_pY416", "PTK2", "PTK2_pY397", "PXN", "PXN_pY118",
           "CRK", "CRKL", "PDPK1", "AKT1", "AKT1_pS473", "AKT1_pT308",
           "MAPK1", "MAPK1_pT202_Y204", "MAP2K1", "STAT3", "STAT3_pY705",
           "CTNNB1", "MYC", "CCND1", "VIM", "CDH1", "FN1", "COL1A1"]

ecm = pd.read_csv(f"{CACHE}/ccle_ecm_prolif.csv", index_col=0)
print(f"ECM 评分 {ecm.shape}")

# 取所有 RPPA 的基因列表
r = S.get(f"{B}/molecular-profiles/ccle_broad_2019_protein_quantification_zscores/molecular-data"
          f"?sampleListId=ccle_broad_2019_all&projection=SUMMARY", timeout=120)
print(f"RPPA fetch {r.status_code}")

# 更稳：用 genes 端点查存在的基因
found = {}
for g in TARGETS:
    try:
        rr = S.get(f"{B}/genes/{g}", timeout=15)
        if rr.status_code == 200:
            found[g] = rr.json()["entrezGeneId"]
    except Exception:
        pass
print(f"基因解析 {len(found)}/{len(TARGETS)}: {list(found)}")

if not found:
    raise SystemExit("无可用蛋白")

pid = "ccle_broad_2019_protein_quantification_zscores"
body = {"entrezGeneIds": list(found.values()), "sampleListId": "ccle_broad_2019_all"}
rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
            data=json.dumps(body), timeout=180)
if rr.status_code != 200:
    raise SystemExit(f"取数失败 {rr.status_code}")
df = pd.DataFrame(rr.json())
print(f"取到 {len(df)} 条记录")

X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
sym = {v: k for k, v in found.items()}
X.index = [sym.get(i, str(i)) for i in X.index]
print(f"蛋白矩阵 {X.shape}: {list(X.index)}")

# 匹配细胞系（CCLE 名形如 143B_BONE）
def base(s):
    p = str(s).split("_")
    return "_".join(p[:-1]) if len(p) > 1 else str(s)

map_ = {}
for i in X.columns:
    map_[i] = i
    map_.setdefault(base(i), i)

rows = []
for prot in X.index:
    vals = X.loc[prot]
    # 对齐 ECM
    pairs = [(map_[c], ecm.loc[c, "ecm"]) for c in X.columns
             if c in ecm.index and pd.notna(vals[c])]
    if len(pairs) < 100:
        continue
    a = np.array([p[1] for p in pairs])
    b = np.array([vals[p[0]] for p in pairs])
    r_s, p_s = stats.spearmanr(a, b)
    rows.append({"protein": prot, "n": len(pairs), "rho": r_s, "p": p_s})
    print(f"  {prot:16s} n={len(pairs):4d}  rho={r_s:+.3f}  p={p_s:.1e}")

R = pd.DataFrame(rows).sort_values("rho", ascending=False)
R.to_csv(f"{OUT}/ecm_rppa_correlation.csv", index=False, encoding="utf-8-sig")
print("\n" + R.round(3).to_string(index=False))
