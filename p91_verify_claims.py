"""逐条验证 Response v5 中的每一项声明"""
import os, re, sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(line_buffering=True)
SRC = r"D:\workbuddy内容\A1_力学代理"
MS = os.path.join(SRC, "manuscript_ECM_v5.0.md")
s = open(MS, encoding="utf-8").read()
body = s[:s.index("## Figure legends")]
legends = s[s.index("## Figure legends"):s.index("## References")]

ok = fail = 0
def chk(label, cond, detail=""):
    global ok, fail
    if cond:
        ok += 1
        print(f"  ✓ {label}  {detail}")
    else:
        fail += 1
        print(f"  ✗ {label}  {detail}")

print("=" * 66)
print("Response v5 声明逐条验证")
print("=" * 66)

print("\n【M1】p 值重算")
for tbl, col in [("TableS1_drug_sensitivity_full.csv", "p (sci)"),
                 ("TableS1b_significant_87.csv", "p (sci)"),
                 ("TableS2_partial_correlation.csv", "p (raw, sci)"),
                 ("TableS3_by_cancer_type.csv", "p (sci)")]:
    T = pd.read_csv(os.path.join(SRC, tbl), comment="#")
    if col in T.columns:
        n_cens = int((T[col].astype(str) == "<1e-300").sum())
        chk(f"{tbl[:20]} 无 <1e-300", n_cens == 0, f"{n_cens} 行")
T1 = pd.read_csv(os.path.join(SRC, "TableS1_drug_sensitivity_full.csv"), comment="#")
d = T1[T1.Drug == "Dasatinib"].iloc[0]
chk("dasatinib p = 1.88e-15", str(d["p (sci)"]) == "1.88e-15", f"实际 {d['p (sci)']}")
sig = pd.read_csv(os.path.join(SRC, "TableS1b_significant_87.csv"), comment="#")
chk("87 个显著", len(sig) == 87, f"{len(sig)} 个")
chk("73 耐药 / 14 敏感", (sig.rho > 0).sum() == 73 and (sig.rho < 0).sum() == 14,
    f"{int((sig.rho>0).sum())} / {int((sig.rho<0).sum())}")

print("\n【M2】假不等式已移除")
chk("无 'p < 1×10⁻¹⁵'", "p < 1×10⁻¹⁵" not in s)
chk("无 'P = 1.9e-15'", "P = 1.9e-15" not in s)
chk("正文有 'p = 1.9×10⁻¹⁵'", "p = 1.9×10⁻¹⁵" in s)
chk("Methods 声明无删失", "none is censored" in s)

print("\n【M3】Figure S1 口径")
s1 = re.search(r"\*\*Figure S1\..*?(?=\n\n\*\*Figure S2)", legends, re.S)
s1t = s1.group(0) if s1 else ""
chk("图注给出 3,310", "3,310" in s1t)
chk("图注列出三个不一致队列", "497 vs 487" in s1t and "525 vs 523" in s1t)
chk("图注说明 age 要求", "requires age" in s1t or "additionally" in s1t)

print("\n【M4】衰减值统一")
att = re.findall(r"attenuation ([−-]?[0-9.]+)%", s)
chk("两处均为 −0.21%", len(att) >= 2 and all(a in ("−0.21", "-0.21") for a in att), str(att))

print("\n【m7】Fig1 图注双 n")
f1 = re.search(r"\*\*Figure 1\..*?(?=\n\n\*\*Figure 2)", legends, re.S)
f1t = f1.group(0) if f1 else ""
chk("图注含 3,291", "3,291" in f1t)
chk("图注含 3,255", "3,255" in f1t)

print("\n【m8】S2 注释在表头")
with open(os.path.join(SRC, "TableS2_partial_correlation.csv"), encoding="utf-8") as f:
    first = f.readline()
chk("首行为注释", first.startswith("#"))

print("\n【其他】")
chk("无 '828 deaths'", "828 deaths" not in s)
chk("无 'two small cohorts'", "two small cohorts" not in s)
chk("无 'five compounds'", "five compounds" not in s)
chk("无 'P = 3.2e-01'", "P = 3.2e-01" not in s)
chk("program 计数 31", len(re.findall(r"program", s)) == 31, str(len(re.findall(r"program", s))))
refs = re.findall(r"^(\d+)\.\s", s[s.index("## References"):], re.M)
cited = set()
for m in re.finditer(r"\[(\d+(?:,\s*\d+)*)\]", body):
    for x in m.group(1).split(","):
        cited.add(int(x))
chk("文献 14 条无孤儿", len(refs) == 14 and not (cited - {int(r) for r in refs}) and not ({int(r) for r in refs} - cited))

print("\n" + "=" * 66)
print(f"通过 {ok}，失败 {fail}")
print("=" * 66)
