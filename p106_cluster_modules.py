# -*- coding: utf-8 -*-
"""
p106_cluster_modules.py
药敏模块聚类 + 置换检验
对应稿件 2.7 节、图 8、图 9、表 S14

输入: GDSC2_lnIC50_matrix.csv（药物 × 细胞系 ln IC50 矩阵）
      TableS1_drug_sensitivity_full.csv（ECM 与药物 rho）
输出: TableS14_cluster_k10_full.csv

用法:
  python p106_cluster_modules.py --matrix GDSC2_lnIC50_matrix.csv --ecm TableS1_drug_sensitivity_full.csv --out TableS14_cluster_k10_full.csv
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="药物 x 细胞系 ln IC50 矩阵")
    ap.add_argument("--ecm", required=True, help="TableS1，含 Drug/rho 列")
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--nperm", type=int, default=2000)
    args = ap.parse_args()

    M = pd.read_csv(args.matrix, index_col=0)
    S1 = pd.read_csv(args.ecm, comment="#")
    C = M.T.corr(method="spearman")
    D = 1 - C.values
    np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), method="average")
    lab = fcluster(Z, args.k, criterion="maxclust")
    mod = pd.Series(lab, index=C.index, name="cluster")

    rho_all = S1.rho.values
    rows = []
    for k in range(1, args.k + 1):
        drugs = mod[mod == k].index.tolist()
        sub = S1[S1.Drug.isin(drugs)]
        if len(sub) == 0:
            continue
        non_member = S1[~S1.Drug.isin(drugs)]
        if len(sub) >= 2 and len(non_member) >= 2:
            u, p = stats.mannwhitneyu(sub.rho.values, non_member.rho.values,
                                      alternative="two-sided")
        else:
            p = float("nan")
        # 置换检验
        rng = np.random.default_rng(42)
        cnt = 0
        direction = "less" if sub.rho.mean() < 0 else "greater"
        for _ in range(args.nperm):
            pm = rho_all[rng.choice(len(rho_all), len(sub), replace=False)].mean()
            if (direction == "less" and pm <= sub.rho.mean()) or \
               (direction == "greater" and pm >= sub.rho.mean()):
                cnt += 1
        p_perm = (cnt + 1) / (args.nperm + 1)

        rows.append({
            "cluster": k, "n_members": len(drugs),
            "rho_mean": round(sub.rho.mean(), 4),
            "rho_min": round(sub.rho.min(), 4),
            "rho_max": round(sub.rho.max(), 4),
            "MW_P_vs_nonmembers": f"{p:.4e}" if not np.isnan(p) else "NA",
            "perm_P": f"{p_perm:.4f}",
            "significant": "yes" if (not np.isnan(p) and p < 0.05) else "no",
            "members": "; ".join(drugs),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False, encoding="utf-8")
    print(f"已写入 {args.out}（{len(out)} 个簇，{(out.significant=='yes').sum()} 个显著）")


if __name__ == "__main__":
    main()
