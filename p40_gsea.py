"""任务2：ECM 高 vs 低细胞系的通路富集（GSEA）"""
import os, sys, json
import numpy as np
import pandas as pd
import requests, time
from scipy import stats

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
OUT = r"D:\stiffness_r3\results"
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

import gseapy as gp

fc = os.path.join(CACHE, "ccle_hallmark_expr.csv")
if os.path.exists(fc):
    X = pd.read_csv(fc, index_col=0)
    print(f"载入 Hallmark 表达 {X.shape}")
else:
    # 先拿 Hallmark 基因集里的全部基因
    sets = gp.get_library("MSigDB_Hallmark_2020")
    genes = sorted({g.upper() for v in sets.values() for g in v})
    print(f"Hallmark 基因 {len(genes)}")
    st = "ccle_broad_2019"
    pid = "ccle_broad_2019_rna_seq_mrna"
    # 一次性取全基因表（44896 条），本地映射
    gr = S.get(f"{B}/genes", timeout=180)
    gmap = {g["hugoGeneSymbol"].upper(): g["entrezGeneId"]
            for g in gr.json() if g.get("hugoGeneSymbol")}
    print(f"全基因表 {len(gmap)}")
    ent, sym = [], {}
    for g in genes:
        e = gmap.get(g)
        if e:
            ent.append(int(e)); sym[str(e)] = g
    print(f"映射到 entrez {len(ent)}")

    frames = []
    for i in range(0, len(ent), 600):
        body = {"entrezGeneIds": ent[i:i+600], "sampleListId": f"{st}_all"}
        try:
            rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                        data=json.dumps(body), timeout=180)
            if rr.status_code == 200 and rr.json():
                frames.append(pd.DataFrame(rr.json()))
        except Exception as e:
            print("batch fail", type(e).__name__)
        time.sleep(0.3)
    if not frames:
        raise SystemExit("全部批次失败")
    df = pd.concat(frames, ignore_index=True)
    X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
    X.index = [sym.get(str(i), str(i)) for i in X.index]
    X = X[~X.index.duplicated(keep="first")]
    X.to_csv(fc)
    print(f"[saved] {X.shape}")

# ECM 评分
ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
gs = [g for g in ECM if g in X.index]
print(f"ECM 基因 {len(gs)}/{len(ECM)}")
score = X.loc[gs].mean(0)
print(f"细胞系 {len(score)}")

# 每个基因与 ECM 的 Spearman 相关
print("计算基因级相关…")
rows = []
for g in X.index:
    v = X.loc[g]
    if v.std() < 1e-6:
        continue
    r, p = stats.spearmanr(score.values, v.values)
    rows.append((g, r, p))
R = pd.DataFrame(rows, columns=["gene", "rho", "p"]).dropna()
print(f"有效基因 {len(R)}")

# GSEA 需要排序的基因列表
R = R.sort_values("rho", ascending=False)
R.to_csv(f"{OUT}/ccle_ecm_gene_correlation.csv", index=False, encoding="utf-8-sig")

# 用 gseapy 做 prerank
import gseapy as gp
rnk = R.set_index("gene")["rho"]
print("\n=== GSEA (Hallmark) ===")
try:
    pre = gp.prerank(rnk=rnk, gene_sets="MSigDB_Hallmark_2020",
                     min_size=10, max_size=500, permutation_num=1000,
                     outdir=None, seed=1, threads=4)
    res = pre.res2d
    res = res.sort_values("NES", ascending=False)
    res.to_csv(f"{OUT}/gsea_hallmark.csv", index=False, encoding="utf-8-sig")
    print(f"  基因集 {len(res)}")
    print("\n  NES 最高（ECM 高时上调）TOP 12:")
    for _, r in res.head(12).iterrows():
        print(f"    {r.Term[:48]:50s} NES={float(r.NES):+.2f}  FDR={float(r.FDR_qval):.3f}")
    print("\n  NES 最低（ECM 高时下调）TOP 6:")
    for _, r in res.tail(6).iterrows():
        print(f"    {r.Term[:48]:50s} NES={float(r.NES):+.2f}  FDR={float(r.FDR_qval):.3f}")
except Exception as e:
    print("GSEA 失败:", type(e).__name__, e)
