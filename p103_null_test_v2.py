"""
M2（重写）：阴性对照基因集分析
分队列缓存，可断点续跑。
"""
import os, sys, json, random, warnings, glob
import numpy as np
import pandas as pd
import requests
from scipy import stats
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(line_buffering=True)

B = "https://www.cbioportal.org/api"
CACHE = "cache"
POOL = os.path.join(CACHE, "randpool")
OUT = "results"
os.makedirs(POOL, exist_ok=True)

STUDIES = {
    "coadread_tcga_pan_can_atlas_2018": "Colon",
    "brca_tcga_pan_can_atlas_2018": "Breast",
    "luad_tcga_pan_can_atlas_2018": "Lung ad.",
    "lusc_tcga_pan_can_atlas_2018": "Lung sq.",
    "paad_tcga_pan_can_atlas_2018": "Pancreas",
    "ucec_tcga_pan_can_atlas_2018": "Endometrium",
}
N_POOL, N_SETS, SET_SIZE, SEED = 800, 500, 23, 42
ECM23 = ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","LOX","LOXL1","LOXL2",
         "FN1","POSTN","CCN2","ACTA2","TAGLN","TIMP1","SPARC","DCN","LUM","TGFB1","COL11A1","COL8A1","FBN1"]

S = requests.Session()
S.headers.update({"accept": "application/json", "Content-Type": "application/json"})


def load_or_fetch(study, label):
    f = os.path.join(POOL, f"{label.replace(' ','_').replace('.','')}.pkl")
    if os.path.exists(f):
        d = pd.read_pickle(f)
        print(f"  {label}: 缓存 {d['X'].shape[0]}g × {d['X'].shape[1]}s")
        return d

    random.seed(SEED)
    g = S.get(f"{B}/genes", timeout=180).json()
    allids = [x["entrezGeneId"] for x in g if x.get("entrezGeneId")]
    inv = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in g if x.get("entrezGeneId")}
    pool = random.sample(allids, N_POOL)
    # 确保 ECM23 在内
    for gene in ECM23:
        gid = next((x["entrezGeneId"] for x in g if x.get("hugoGeneSymbol", "").upper() == gene), None)
        if gid and gid not in pool:
            pool.append(gid)

    profs = S.get(f"{B}/studies/{study}/molecular-profiles", timeout=60).json()
    pid = next(p["molecularProfileId"] for p in profs
               if "zscore" in p["molecularProfileId"].lower()
               and p.get("molecularAlterationType") == "MRNA_EXPRESSION"
               and "mirna" not in p["molecularProfileId"].lower())
    r = S.post(f"{B}/molecular-profiles/{pid}/molecular-data/fetch",
               data=json.dumps({"entrezGeneIds": pool, "sampleListId": f"{study}_all"}),
               timeout=240)
    if r.status_code != 200:
        raise RuntimeError(f"{label} fetch {r.status_code}")
    X = pd.DataFrame(r.json()).pivot_table(index="entrezGeneId", columns="sampleId",
                                           values="value", aggfunc="mean")
    X.index = [inv.get(i, str(i)) for i in X.index]

    cc = S.get(f"{B}/studies/{study}/clinical-data",
               params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=120)
    C = pd.DataFrame(cc.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                            values="value", aggfunc="first")
    sm = S.get(f"{B}/studies/{study}/samples", timeout=60).json()
    mp = {x["sampleId"]: x["patientId"] for x in sm if x.get("sampleId") in set(X.columns)}
    X = X.rename(columns=mp).T.groupby(level=0).mean().T
    common = X.columns.intersection(C.index)

    T = C.loc[common].copy()
    T["OS_MONTHS"] = pd.to_numeric(T.OS_MONTHS, errors="coerce")
    T["AGE"] = pd.to_numeric(T.AGE, errors="coerce")
    T = T.dropna(subset=["OS_MONTHS", "OS_STATUS"])
    T = T[T.OS_MONTHS > 0]
    T["event"] = T.OS_STATUS.astype(str).str.upper().str.contains("1|DECEASED|DIED").astype(int)

    X = X[T.index]
    X = X[X.std(axis=1) > 0.05]
    d = {"X": X, "T": T}
    pd.to_pickle(d, f)
    print(f"  {label}: {X.shape[0]}g × {X.shape[1]}s, {int(T.event.sum())} events  [saved]")
    return d


print("=== 阶段 1：取数 ===")
data = {}
for study, label in STUDIES.items():
    try:
        data[label] = load_or_fetch(study, label)
    except Exception as e:
        print(f"  ✗ {label}: {type(e).__name__} {e}")

print(f"\n成功 {len(data)}/6 队列")
if len(data) < 4:
    sys.exit(1)


def run_meta(gene_sets, data):
    out = []
    for genes in gene_sets:
        hrs, ses = [], []
        for label, d in data.items():
            X, T = d["X"], d["T"]
            gs = [x for x in genes if x in X.index]
            if len(gs) < 5:
                continue
            score = X.loc[gs].mean(0)
            TT = T.copy()
            TT["z"] = (score - score.mean()) / (score.std() + 1e-9)
            TT = TT.dropna(subset=["z", "OS_MONTHS", "event"])
            if len(TT) < 50:
                continue
            try:
                c = CoxPHFitter()
                c.fit(TT[["z", "OS_MONTHS", "event"]], duration_col="OS_MONTHS", event_col="event")
                hrs.append(float(c.params_["z"]))
                ses.append(float(c.standard_errors_["z"]))
            except Exception:
                continue
        if len(hrs) < 4:
            out.append(np.nan); continue
        w = 1 / np.array(ses) ** 2
        out.append(float(np.exp((w * np.array(hrs)).sum() / w.sum())))
    return np.array(out)


print("\n=== 阶段 2：ECM 参照 ===")
ecm_in = [g for g in ECM23 if all(g in d["X"].index for d in data.values())]
print(f"ECM 基因可用 {len(ecm_in)}/{len(ECM23)}")
ecm_hr = run_meta([ecm_in], data)[0]
print(f"ECM pooled HR = {ecm_hr:.4f}")

print("\n=== 阶段 3：随机集 ===")
shared = sorted(set.intersection(*[set(d["X"].index) for d in data.values()]))
shared = [g for g in shared if g not in ECM23]
print(f"共有非 ECM 基因 {len(shared)} 个")

rng = np.random.default_rng(SEED)
sets = [list(rng.choice(shared, SET_SIZE, replace=False)) for _ in range(N_SETS)]
null = []
for i in range(0, len(sets), 25):
    null.extend(run_meta(sets[i:i + 25], data).tolist())
    if (i // 25) % 4 == 0:
        print(f"  {i + 25}/{N_SETS}，中位 HR {np.nanmedian(null):.4f}")
        sys.stdout.flush()

null = np.array([x for x in null if np.isfinite(x)])
p_emp = float((null >= ecm_hr).mean())

print("\n" + "=" * 58)
print(f"队列数             {len(data)}")
print(f"ECM23 pooled HR    {ecm_hr:.4f}")
print(f"随机集 n           {len(null)}")
print(f"随机 HR 中位数      {np.median(null):.4f}")
print(f"随机 HR 均值±SD     {null.mean():.4f} ± {null.std():.4f}")
print(f"随机 HR 5%/95%分位  {np.percentile(null,5):.4f} / {np.percentile(null,95):.4f}")
print(f"随机 HR 99%分位     {np.percentile(null,99):.4f}")
print(f"经验 p（随机≥ECM）  {p_emp:.4f}")
print(f"ECM 百分位          {(null < ecm_hr).mean()*100:.1f}%")
print("=" * 58)

pd.DataFrame({"null_HR": null}).to_csv(os.path.join(OUT, "null_distribution.csv"), index=False)
pd.DataFrame([{"ECM_HR": ecm_hr, "k_cohorts": len(data), "null_median": float(np.median(null)),
               "null_sd": float(null.std()), "null_p95": float(np.percentile(null, 95)),
               "null_p99": float(np.percentile(null, 99)), "empirical_p": p_emp,
               "ecm_percentile": float((null < ecm_hr).mean() * 100),
               "n_random": len(null)}]).to_csv(os.path.join(OUT, "null_test_summary.csv"), index=False)
print("\n[saved]")
