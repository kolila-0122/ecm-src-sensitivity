"""ECM 模块与免疫治疗反应"""
import os, sys, json, time
import numpy as np
import pandas as pd
import requests
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
OUT = r"D:\stiffness_r3\results"
CACHE = r"D:\stiffness_r3\cache"
os.makedirs(CACHE, exist_ok=True)
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

STUDIES = {
    "mel_iatlas_riaz_nivolumab_2017": "Riaz 2017 (nivo)",
    "mel_iatlas_hugo_ucla_2016": "Hugo 2016 (pembro)",
    "mel_iatlas_liu_2019": "Liu 2019 (anti-PD1)",
    "nsclc_pd1_msk_2018": "NSCLC PD-1 (MSK)",
    "mel_ucla_2016": "UCLA 2016 (pembro)",
}
# 疗效字段优先级
RESP_KEYS = ["RESPONDER", "RESPONSE", "CLINICAL_BENEFIT", "DURABLE_CLINICAL_BENEFIT",
             "TREATMENT_RESPONSE"]

e2s = {}
for g in ECM:
    try:
        r = S.get(f"{B}/genes/{g}", timeout=20)
        if r.status_code == 200:
            e2s[g] = r.json()["entrezGeneId"]
    except Exception:
        pass
sym = {v: k for k, v in e2s.items()}
print(f"ECM 基因 {len(e2s)}\n")

rows = []
allscore = []
for study, label in STUDIES.items():
    fc = os.path.join(CACHE, f"ici2_{study}.csv")
    if os.path.exists(fc):
        T = pd.read_csv(fc, index_col=0)
    else:
        profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
        pid = None
        for p in profs:
            if p.get("molecularAlterationType") == "MRNA_EXPRESSION" \
               and "mirna" not in p["molecularProfileId"].lower():
                pid = p["molecularProfileId"]; break
        if not pid:
            print(f"{label}: 无 mRNA"); continue
        body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{study}_all"}
        rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                    data=json.dumps(body), timeout=120)
        if rr.status_code != 200:
            print(f"{label}: mRNA 失败"); continue
        df = pd.DataFrame(rr.json())
        X = df.pivot_table(index="entrezGeneId", columns="sampleId",
                           values="value", aggfunc="mean")
        X.index = [sym.get(i, str(i)) for i in X.index]
        gs = [g for g in ECM if g in X.index]
        score = X.loc[gs].mean(0)

        cc = S.get(f"{B}/studies/{study}/clinical-data",
                   params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=90)
        C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                                values="value", aggfunc="first")
        keep = [c for c in RESP_KEYS if c in C.columns]
        # 用官方 sample->patient 映射
        sm = S.get(f"{B}/studies/{study}/samples", timeout=40).json()
        m = {x["sampleId"]: x["patientId"] for x in sm if x.get("sampleId") in set(score.index)}
        sc = pd.Series({p: score[s].mean() for s, p in m.items()})
        T = pd.DataFrame({"ecm": sc}).join(C[keep], how="inner")
        T.index.name = "patient"
        T.to_csv(fc)
        print(f"{label}: 缓存 {T.shape} 疗效字段 {keep}")

    if T.empty:
        print(f"{label}: 空"); continue
    # 选响应字段
    rk = None
    for k in RESP_KEYS:
        if k in T.columns and T[k].notna().sum() >= 10:
            rk = k; break
    if not rk:
        print(f"{label}: 无可用疗效字段"); continue

    v = T[rk].astype(str).str.upper()
    # 解析响应：CR/PR/R = 1；PD/SD/NR = 0
    pos = v.str.contains("CR|PR|^R$|RESPOND|YES|TRUE|BENEFIT|COMPLETE|PARTIAL", regex=True)
    neg = v.str.contains("PD|SD|NR|NON|NO$|FALSE|PROGRESS", regex=True)
    resp = pd.Series(np.nan, index=T.index)
    resp[pos & ~neg] = 1
    resp[neg] = 0
    d = pd.DataFrame({"ecm": T.ecm, "resp": resp}).dropna()
    if len(d) < 20 or d.resp.nunique() < 2:
        print(f"{label}: 有效样本不足 (n={len(d)}, 字段={rk})"); continue

    r_idx = d.resp == 1
    mw = stats.mannwhitneyu(d.ecm[r_idx], d.ecm[~r_idx], alternative="two-sided")
    auc = roc_auc_score(d.resp, d.ecm)
    # 方向：若 AUC<0.5 则翻转，报告 AUC 与方向
    n1, n0 = int(r_idx.sum()), int((~r_idx).sum())
    rows.append({"study": label, "field": rk, "n": len(d), "resp": n1, "nonresp": n0,
                 "AUC": auc, "median_resp": float(d.ecm[r_idx].median()),
                 "median_nonresp": float(d.ecm[~r_idx].median()),
                 "p_mwu": float(mw.pvalue),
                 "diff": float(d.ecm[r_idx].median() - d.ecm[~r_idx].median())})
    print(f"{label:24s} n={len(d):3d} (R={n1}, NR={n0})  AUC={auc:.3f}  "
          f"ECM 中位 R={d.ecm[r_idx].median():+.2f} vs NR={d.ecm[~r_idx].median():+.2f}  "
          f"p={mw.pvalue:.3f}")
    d["study"] = label
    allscore.append(d)
    time.sleep(0.4)

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/ecm_ici_auc.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 78)
print(R.round(3).to_string(index=False))

if len(R):
    # 合并 AUC（加权平均，Fisher z）
    w = R.n.values
    z = np.arctanh(np.clip(R.AUC.values, .01, .99))
    zbar = (w * z).sum() / w.sum()
    se = 1 / np.sqrt((w - 3).sum())
    pooled = np.tanh(zbar)
    print(f"\n加权合并 AUC = {pooled:.3f} "
          f"(95% CI {np.tanh(zbar - 1.96 * se):.3f}-{np.tanh(zbar + 1.96 * se):.3f})")
    print(f"方向: {int((R.AUC > 0.5).sum())}/{len(R)} 个队列 AUC>0.5（ECM 高→反应差，应为 AUC<0.5）")
    print(f"总样本 {int(R.n.sum())}，响应 {int(R.resp.sum())}")

if allscore:
    A = pd.concat(allscore)
    A.to_csv(f"{OUT}/ecm_ici_all.csv")
    print(f"\n合并数据已存 ({len(A)} 例)")
