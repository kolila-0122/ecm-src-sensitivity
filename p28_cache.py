"""步骤 1：逐队列缓存 mRNA + 临床到本地 CSV"""
import os, sys, json, time
import numpy as np
import pandas as pd
import requests

sys.stdout.reconfigure(line_buffering=True)
B = "https://www.cbioportal.org/api"
CACHE = r"D:\stiffness_r3\cache"
os.makedirs(CACHE, exist_ok=True)
S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}

for study, label in STUDIES.items():
    f_out = os.path.join(CACHE, f"{study}.csv")
    if os.path.exists(f_out):
        print(f"{label}: 已缓存"); continue

    profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=40).json()
    pid = None
    for p in profs:
        if p.get("molecularAlterationType") == "MRNA_EXPRESSION" \
           and "zscore" in p["molecularProfileId"].lower() \
           and "mirna" not in p["molecularProfileId"].lower():
            pid = p["molecularProfileId"]; break
    if not pid:
        print(f"{label}: 无 mRNA 谱"); continue

    # 基因映射
    e2s = {}
    for g in ECM:
        try:
            r = S.get(f"{B}/genes/{g}", timeout=20)
            if r.status_code == 200:
                e2s[g] = r.json()["entrezGeneId"]
        except Exception:
            pass
    sym = {v: k for k, v in e2s.items()}

    # mRNA
    body = {"entrezGeneIds": list(e2s.values()), "sampleListId": f"{study}_all"}
    rr = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
                data=json.dumps(body), timeout=120)
    if rr.status_code != 200:
        print(f"{label}: mRNA 失败 {rr.status_code}"); continue
    df = pd.DataFrame(rr.json())
    X = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value", aggfunc="mean")
    X.index = [sym.get(i, str(i)) for i in X.index]
    gs = [g for g in ECM if g in X.index]
    score = X.loc[gs].mean(0)

    # 临床
    cc = S.get(f"{B}/studies/{study}/clinical-data",
               params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=90)
    if cc.status_code != 200:
        print(f"{label}: 临床失败 {cc.status_code}"); continue
    C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                            values="value", aggfunc="first")

    smap = pd.Series(score.index.str[:12], index=score.index)
    sc = score.groupby(smap).mean()
    keep = [c for c in ["OS_MONTHS", "OS_STATUS", "AGE", "SEX",
                        "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATHOLOGIC_STAGE_AJCC7_LOWERCASE",
                        "NEOPLASM_DISEASEPATHOLOGIC_STAGE"] if c in C.columns]
    T = pd.DataFrame({"ecm_score": sc}).join(C[keep], how="inner")
    T.index.name = "patient"
    T.to_csv(f_out)
    print(f"{label}: 缓存 {T.shape}  字段 {keep}  → {os.path.basename(f_out)}")
    time.sleep(0.5)

print("\n[缓存完成]")
