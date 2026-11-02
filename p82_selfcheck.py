"""投稿前自检：逐项验证所有将在 Response 中声明的整改"""
import os, re, sys
import pandas as pd

sys.stdout.reconfigure(line_buffering=True)
SRC = r"D:\workbuddy内容\A1_力学代理"
MS = os.path.join(SRC, "manuscript_ECM_v3.0.md")
s = open(MS, encoding="utf-8").read()

body = s[:s.index("## Figure legends")]
legends = s[s.index("## Figure legends"):s.index("## References")]
refs_sec = s[s.index("## References"):]

print("=" * 60)
print("自检报告")
print("=" * 60)

# ---- 1. 图号引用计数 ----
print("\n[1] 图号引用计数（正文，不含图注）")
figs = ["Figure 1", "Figure 2", "Figure S1", "Figure S2", "Figure S3", "Figure S4", "Figure S5"]
counts = {}
for f in figs:
    # 匹配 Figure 1 / Figure S1，但要避免 Figure S1 被 Figure 1 计入
    pat = re.escape(f) + r"(?![0-9])"
    counts[f] = len(re.findall(pat, body))
for f in figs:
    flag = "✗ 零引用" if counts[f] == 0 else "✓"
    print(f"  {f:12s} {counts[f]:2d} 次  {flag}")

# ---- 2. 图注数量 ----
print("\n[2] 图注数量")
leg_figs = re.findall(r"\*\*Figure (S?\d+)\.", legends)
print(f"  图注: {leg_figs}")

# ---- 3. p 值写法一致性 ----
print("\n[3] p 值写法检查")
pats = {
    "科学计数法 1.9×10⁻¹⁵ 类": re.findall(r"1\.9×10[-−]15", s),
    "数学记法 < 1×10⁻¹⁵": re.findall(r"<\s*1×10[-−]15", s),
    "p = 0.0": re.findall(r"[pP]\s*=\s*0\.0(?![0-9])", s),
    "P = 1.9e-15 (计算机记法)": re.findall(r"[pP]\s*=\s*\d\.\de-\d+", s),
}
for k, v in pats.items():
    print(f"  {k:28s} {len(v)} 处  {'✓' if len(v)==0 else '需检查'}")

# ---- 4. 双句点 ----
print("\n[4] 双句点")
dbl = re.findall(r"\.\.(?!\.)", s)
print(f"  连续双句点: {len(dbl)} 处  {'✓' if len(dbl)==0 else '✗ '+str(dbl[:3])}")

# ---- 5. program / programme ----
print("\n[5] 用词统一")
print(f"  'programme' 出现: {s.count('programme')} 次")
print(f"  'program' 出现: {s.count('program')} 次")

# ---- 6. 关键数字一致性 ----
print("\n[6] 关键数字")
checks = {
    "809 deaths": s.count("809 deaths"),
    "828 deaths": s.count("828 deaths"),
    "3,255": s.count("3,255"),
    "3,310": s.count("3,310"),
    "AUC = 0.496": s.count("0.496"),
    "AUC = 0.476 (旧)": s.count("0.476"),
    "n = 205": s.count("n = 205"),
    "n = 232 (旧)": s.count("n = 232"),
}
for k, v in checks.items():
    print(f"  {k:20s} {v} 次")

# ---- 7. 参考文献 ----
print("\n[7] 参考文献")
listed = re.findall(r"^(\d+)\.\s", refs_sec, re.M)
cited = set()
for m in re.finditer(r"\[(\d+(?:,\s*\d+)*)\]", body):
    for x in m.group(1).split(","):
        cited.add(int(x))
L = set(int(x) for x in listed)
print(f"  文献条数: {len(L)}")
print(f"  正文引用: {sorted(cited)}")
print(f"  引用未列出: {sorted(cited - L)}  {'✓' if not (cited-L) else '✗'}")
print(f"  列出未引用: {sorted(L - cited)}  {'✓' if not (L-cited) else '✗'}")

# ---- 8. 补充表检查 ----
print("\n[8] 补充表")
for f in sorted(os.listdir(SRC)):
    if not (f.startswith("TableS") and f.endswith(".csv")):
        continue
    p = os.path.join(SRC, f)
    with open(p, "rb") as fh:
        head = fh.read(3)
    bom = head == b"\xef\xbb\xbf"
    T = pd.read_csv(p)
    zero_p = 0
    for c in T.columns:
        if "p" in str(c).lower() or "fdr" in str(c).lower():
            col = T[c].astype(str)
            zero_p += int((col == "0.0").sum())
    print(f"  {f:40s} BOM={'✗' if bom else '✓'}  0.0值={zero_p}  {'✓' if not bom and zero_p==0 else '✗'}")

# ---- 9. GSEA 数量 ----
print("\n[9] GSEA 数量口径")
if "25 at FDR < 0.05" in s and "20 with the largest" in s:
    print("  ✓ 正文与图注均说明 25 通过 / 20 展示")
else:
    print("  ✗ 未统一")

# ---- 10. 投稿包文件清单 ----
print("\n[10] 投稿包文件白名单")
pkg = os.path.join(SRC, "JENCI_v3")
allowed_ext = {".docx", ".md", ".tif", ".csv"}
bad = []
for root, dirs, files in os.walk(pkg):
    for f in files:
        ext = os.path.splitext(f)[1].lower()
        if ext not in allowed_ext:
            bad.append(f)
        if "TODO" in f.upper() or "ZENODO" in f.upper():
            bad.append(f)
print(f"  可疑文件: {bad if bad else '无  ✓'}")

print("\n" + "=" * 60)
