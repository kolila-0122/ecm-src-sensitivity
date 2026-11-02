"""ECM × CD8 免疫浸润 交互作用分析"""
import os, sys, json, time
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
OUT = r"D:\stiffness_r3\results"
S = __import__("requests").Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
# CD8 / 细胞毒 签名（Rooney 2015 类）
IMM = ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7",
       "CD3D", "CD3E", "CXCL9", "CXCL10"]

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}

GENES = ECM + IMM
e2s = {}
for g in GENES:
    try:
        r = S.get(f"{B}/genes/{g}", timeout=20)
        if r.status_code == 200:
            e2s[g] = r.json()["entrezGeneId"]
    except Exception:
        pass
sym = {v: k for k, v in e2s.items()}
print(f"基因映射 {len(e2s)}/{len(GENES)}")

rows = []
allsub = []
for study, label in STUDIES.items():
    fc = os.path.join(CACHE, f"ecmimm_{study}.csv")
    if os.path.exists(fc):
        T = pd.read_csv(fc, index_col=0)
    else:
        profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
        pid = None
        for p in profs:
            if p.get("molecularAlterationType") == "MRNA_EXPRESSION" \
               and "zscore" in p["molecularProfileId"].lower() \
               and "mirna" not in p["molecularProfileId"].lower():
                pid = p["molecularProfileId"]; break
        if not pid:
            print(f"{label}: 无 mRNA"); continue
        body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{study}_all"}
        rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                    data=json.dumps(body), timeout=120)
        if rr.status_code != 200:
            print(f"{label}: 取数失败"); continue
        df = pd.DataFrame(rr.json())
        X = df.pivot_table(index="entrezGeneId", columns="sampleId",
                           values="value", aggfunc="mean")
        X.index = [sym.get(i, str(i)) for i in X.index]

        ecm_g = [g for g in ECM if g in X.index]
        imm_g = [g for g in IMM if g in X.index]
        sc = pd.DataFrame({"ecm": X.loc[ecm_g].mean(0),
                           "imm": X.loc[imm_g].mean(0)})
        sm = S.get(f"{B}/studies/{study}/samples", timeout=40).json()
        m = {x["sampleId"]: x["patientId"] for x in sm if x.get("sampleId") in set(sc.index)}
        sc.index = [m.get(i, i) for i in sc.index]
        sc = sc.groupby(level=0).mean()

        cc = S.get(f"{B}/studies/{study}/clinical-data",
                   params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=90)
        C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                                values="value", aggfunc="first")
        keep = [c for c in ["OS_MONTHS", "OS_STATUS", "AGE", "SEX",
                            "AJCC_PATHOLOGIC_TUMOR_STAGE"] if c in C.columns]
        T = sc.join(C[keep], how="inner")
        T.index.name = "patient"
        T.to_csv(fc)
        print(f"{label}: 缓存 {T.shape}  ECM基因{len(ecm_g)} 免疫基因{len(imm_g)}")

    T["OS_MONTHS"] = pd.to_numeric(T["OS_MONTHS"], errors="coerce")
    T["AGE"] = pd.to_numeric(T["AGE"], errors="coerce")
    T = T.dropna(subset=["OS_MONTHS", "OS_STATUS", "ecm", "imm"])
    T = T[T.OS_MONTHS > 0]
    T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)
    if len(T) < 60:
        print(f"{label}: 样本不足"); continue

    # 标准化
    for c in ["ecm", "imm"]:
        T[c + "_z"] = (T[c] - T[c].mean()) / (T[c].std() + 1e-9)

    # 中位分组
    T["ecm_hi"] = (T.ecm_z > 0).astype(int)
    T["imm_hi"] = (T.imm_z > 0).astype(int)
    T["grp"] = T.ecm_hi * 2 + T.imm_hi   # 0=lo/lo, 1=lo/hi, 2=hi/lo, 3=hi/hi

    r_pear = stats.pearsonr(T.ecm, T.imm)
    print(f"\n{label:12s} n={len(T):4d}  ECM-免疫相关 r={r_pear[0]:+.2f} (p={r_pear[1]:.1e})")

    # 四组生存
    grp_stat = {}
    for gv, gname in [(0, "ECM低/免疫低"), (1, "ECM低/免疫高"), (2, "ECM高/免疫低"), (3, "ECM高/免疫高")]:
        m = T.grp == gv
        if m.sum() >= 8:
            grp_stat[gname] = (m.sum(), float(T.OS_MONTHS[m].median()), float(T.event[m].mean()))
            print(f"    {gname}: n={m.sum():3d} 中位OS={T.OS_MONTHS[m].median():6.1f} 事件率={T.event[m].mean():.2f}")

    # Cox: ECM + 免疫 + 交互
    d = T[["OS_MONTHS", "event", "ecm_z", "imm_z"]].copy()
    d["inter"] = d.ecm_z * d.imm_z
    try:
        c = CoxPHFitter()
        c.fit(d, duration_col="OS_MONTHS", event_col="event")
        hr_i = float(np.exp(c.params_["inter"]))
        p_i = float(c.summary.loc["inter", "p"])
        ci_i = c.confidence_intervals_.loc["inter"].values
        hr_e = float(np.exp(c.params_["ecm_z"]))
        p_e = float(c.summary.loc["ecm_z", "p"])
        hr_m = float(np.exp(c.params_["imm_z"]))
        p_m = float(c.summary.loc["imm_z", "p"])
        print(f"    Cox: ECM HR={hr_e:.2f} (p={p_e:.3f}) | 免疫 HR={hr_m:.2f} (p={p_m:.3f}) | "
              f"交互 HR={hr_i:.2f} (p={p_i:.3f})")
    except Exception as e:
        print(f"    Cox 失败 {type(e).__name__}")
        continue

    # 最差组 vs 最好组
    worst = T.grp == 2   # ECM高/免疫低
    best = T.grp == 1    # ECM低/免疫高
    hr_wb, p_wb, lo_wb, hi_wb = np.nan, np.nan, np.nan, np.nan
    if worst.sum() >= 8 and best.sum() >= 8:
        lr = logrank_test(T.OS_MONTHS[worst], T.OS_MONTHS[best],
                          T.event[worst], T.event[best])
        # 简单 HR 估计（两样本 Cox）
        dd = pd.DataFrame({"T": pd.concat([T.OS_MONTHS[worst], T.OS_MONTHS[best]]),
                           "E": pd.concat([T.event[worst], T.event[best]]),
                           "w": [1] * int(worst.sum()) + [0] * int(best.sum())})
        c2 = CoxPHFitter()
        c2.fit(dd, duration_col="T", event_col="E")
        hr_wb = float(np.exp(c2.params_["w"]))
        p_wb = float(c2.summary.loc["w", "p"])
        ci_wb = c2.confidence_intervals_.loc["w"].values
        lo_wb, hi_wb = float(np.exp(ci_wb[0])), float(np.exp(ci_wb[1]))
        print(f"    最差(ECM高/免疫低) vs 最好(ECM低/免疫高): HR={hr_wb:.2f}  logrank p={lr.p_value:.1e}")

    rows.append({"study": label, "n": len(T), "r_ecm_imm": float(r_pear[0]),
                 "HR_ecm": hr_e, "p_ecm": p_e, "HR_imm": hr_m, "p_imm": p_m,
                 "HR_inter": hr_i, "p_inter": p_i,
                 "HR_worst_vs_best": hr_wb, "lo_wb": lo_wb, "hi_wb": hi_wb,
                 "p_worst": p_wb})
    T["study"] = label
    allsub.append(T.reset_index())

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/ecm_imm_interaction.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 100)
print(R.round(3).to_string(index=False))

if len(R):
    print(f"\n交互项显著（p<0.05）: {int((R.p_inter < .05).sum())}/{len(R)}")
    print(f"交互 HR>1（ECM 放大免疫的保护作用）: {int((R.HR_inter > 1).sum())}/{len(R)}")
    print(f"最差 vs 最好 HR 中位: {R.HR_worst_vs_best.median():.2f}")

if len(R):
    sub = R.dropna(subset=["HR_worst_vs_best"])
    if len(sub) >= 3:
        lg = np.log(sub.HR_worst_vs_best.values)
        se = (np.log(sub.hi_wb.values) - np.log(sub.lo_wb.values)) / (2*1.96)
        w = 1/se**2
        pl = (w*lg).sum()/w.sum(); ps = np.sqrt(1/w.sum())
        Q = float((w*(lg-pl)**2).sum())
        I2 = max(0,(Q-(len(sub)-1))/Q)*100 if Q>0 else 0
        from scipy.stats import norm
        print(f"\n=== 最差 vs 最好 Meta ===")
        print(f"合并 HR = {np.exp(pl):.2f} ({np.exp(pl-1.96*ps):.2f}-{np.exp(pl+1.96*ps):.2f}), "
              f"z={pl/ps:.2f}, p={2*(1-norm.cdf(abs(pl/ps))):.4f}, I²={I2:.1f}%")

if allsub:
    pd.concat(allsub).to_csv(f"{OUT}/ecm_imm_all.csv", index=False)
    print("合并数据已存")
