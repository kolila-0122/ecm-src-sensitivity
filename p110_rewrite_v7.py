"""v7：重新定位为药敏复现研究 —— 重写标题、摘要、Introduction 落点、Results 顺序、Discussion 重心"""

import re

SRC = r"D:\workbuddy内容\A1_力学代理\manuscript_ECM_v6.0.md"
OUT = r"D:\workbuddy内容\A1_力学代理\manuscript_ECM_v7.0.md"

s = open(SRC, encoding="utf-8").read()


def sec(name):
    k = f"## {name}"
    if k not in s:
        return ""
    i = s.index(k)
    rest = s[i + len(k):]
    m = re.search(r"^## ", rest, re.M)
    return rest[:m.start()] if m else rest


def replace_sec(name, new_body):
    global s
    k = f"## {name}"
    i = s.index(k)
    rest = s[i + len(k):]
    m = re.search(r"^## ", rest, re.M)
    end = i + len(k) + (m.start() if m else len(rest))
    s = s[:i] + k + "\n\n" + new_body.strip() + "\n\n" + s[end:]


# ============ 1. 标题 ============
s = re.sub(
    r"^# .*$",
    "# Cross-resource replication of extracellular-matrix–drug associations identifies "
    "EGFR/ERBB inhibitors as the dominant resistance class and Src inhibition as a "
    "reproducible sensitivity signal",
    s, count=1, flags=re.M)

# ============ 2. 摘要 ============
new_abs = """
**Background.** Drug-sensitivity screens such as GDSC2 generate hundreds of associations between a molecular feature and compound response. Most are reported once and never tested in an independent screen. We asked which extracellular matrix (ECM)-associated drug sensitivities survive cross-resource replication, and whether the surviving signals share a target class.

**Methods.** We defined a 23-gene ECM program score and correlated it with ln IC₅₀ for 286 compounds in GDSC2 (969 cell lines). Compounds reaching FDR < 0.05 were grouped by annotated target and tested for class enrichment. Eight representative compounds spanning both directions were re-analysed in GDSC1, an independent screen with a different library and protocol. As a secondary analysis, the same score was tested against overall survival in six TCGA PanCancer Atlas cohorts (n = 3,291) and against immunotherapy and chemotherapy response.

**Results.** Eighty-seven of 286 compounds correlated with the ECM score at FDR < 0.05, of which **73 indicated resistance and 14 sensitivity**. Target-class analysis showed that the only class significantly over-represented among the resistance associations was **EGFR/ERBB** (8 of 12 EGFR/ERBB compounds resistant versus 65 of 274 others; odds ratio 6.4, p = 0.0026), whereas mitotic inhibitors were not enriched (1 of 9 significant). Of the eight compounds tested in GDSC1, **four replicated** (dasatinib, trametinib, refametinib, vorinostat) and four did not. The Src inhibitor dasatinib was the only replicated sensitivity signal and the only one that strengthened on replication (GDSC2 ρ = −0.386; GDSC1 ρ = −0.586), consistent in eight of nine cancer types and unchanged after adjustment for proliferation. The ECM score did not predict immunotherapy (3 cohorts, n = 205; pooled AUC 0.496) or chemotherapy response (4 cohorts, n = 849; pooled AUC 0.457). In the secondary survival analysis the pooled hazard ratio was 1.085 per SD (95% CI 1.019–1.155), but a 500-set random-gene null placed this at only the 88th percentile (empirical p = 0.12), so the survival association cannot be separated from what any 23-gene average produces in these cohorts.

**Conclusions.** ECM-associated drug sensitivity is dominated by resistance in RTK signalling, and the reproducible sensitivity signal is confined to Src inhibition. Reporting a single-screen association without a class-level and cross-resource check risks presenting resource-specific noise as biology.
"""

replace_sec("Abstract", new_abs)

# ============ 3. Introduction 重写 ============
new_intro = """
Drug-sensitivity screens have made it routine to correlate a molecular feature with the response of hundreds of cell lines to hundreds of compounds. GDSC2 alone supports thousands of such correlations [5,11]. The usual workflow reports the strongest ones, discusses the biology they suggest, and moves on. What is rarely reported is how many of those associations survive when the same question is asked in a different screen.

This matters because the sources of resource-specific association are numerous and unremarkable: different compound libraries, different concentration ranges, different curve-fitting procedures, and different cell-line panels. An association that appears in one screen and not another may be biology, or may be the screen.

The extracellular matrix (ECM) is a useful test case. Matrix abundance is associated with poor outcome across solid tumours [1–4], and there is a plausible mechanistic story linking matrix density to drug penetration and to integrin-dependent survival signalling. A screen-wide correlation between an ECM expression score and drug response should therefore produce a coherent set of associations — and if it does not, that is informative about both the matrix biology and the screening literature.

We therefore took a deliberately narrow approach. We computed a 23-gene ECM program score, correlated it with all 286 compounds in GDSC2, grouped the significant associations by target class, and then re-tested representative compounds in GDSC1. We report three things: the direction asymmetry of the associations, the target class in which they concentrate, and which of them replicate. As a secondary analysis we asked whether the same score carries prognostic information, and report that result with the null comparison that determines how it should be read.
"""

replace_sec("Introduction", new_intro)

# ============ 4. Results 重排 ============
res_body = sec("Results")

# 提取各小节
def grab_sub(title):
    i = res_body.find(title)
    if i < 0:
        return ""
    rest = res_body[i:]
    m = re.search(r"\n### ", rest[1:])
    return rest[:m.start() + 1] if m else rest


sub_drug = grab_sub("### The ECM program systematically reshapes drug sensitivity")
sub_src = grab_sub("### Src inhibition is the best-supported signal, and it replicates")
sub_mech = grab_sub("### The ECM program is coupled to EMT, adhesion and a slow-cycling state")
sub_surv = grab_sub("### An ECM program shows a small survival association that does not exceed a random-gene-set null")
sub_immune = grab_sub("### The ECM program does not predict immunotherapy or chemotherapy response")
sub_inter = grab_sub("### A post-hoc subgroup analysis suggests synergy but does not survive formal testing")

# 新顺序：药敏 → Src → 机制 → 治疗反应 → 生存（降级） → 交互（合并到最后）
new_order = [
    sub_drug.replace("### The ECM program systematically reshapes drug sensitivity",
                     "### The ECM program is associated with drug sensitivity in 87 of 286 compounds, predominantly resistance"),
    sub_src.replace("### Src inhibition is the best-supported signal, and it replicates",
                    "### Four of eight tested associations replicate, and only Src inhibition strengthens"),
    sub_mech.replace("### The ECM program is coupled to EMT, adhesion and a slow-cycling state",
                     "### The programme is coupled to EMT, adhesion and a slow-cycling state"),
    sub_immune.replace("### The ECM program does not predict immunotherapy or chemotherapy response",
                       "### The programme does not predict immunotherapy or chemotherapy response"),
    sub_surv.replace("### An ECM program shows a small survival association that does not exceed a random-gene-set null",
                     "### Secondary analysis: the survival association does not exceed a random-gene-set null"),
    sub_inter,
]
replace_sec("Results", "\n\n".join(x.strip() for x in new_order if x.strip()))

# ============ 5. Discussion 重写开头 ============
disc = sec("Discussion")
m = re.search(r"Three findings emerge\.", disc)
if m:
    head = disc[:m.start()]
    tail = disc[m.start():]
    new_tail_head = ("Three findings emerge. First, ECM-associated drug sensitivity in GDSC2 is "
                     "overwhelmingly directional: 73 of 87 significant associations indicate resistance, "
                     "and the only target class over-represented among them is EGFR/ERBB (odds ratio 6.4). "
                     "Second, of the eight associations we could test in an independent screen, four "
                     "replicated — but only one, Src inhibition by dasatinib, strengthened on replication. "
                     "Third, the prognostic signal that the same score carries does not exceed a "
                     "random-gene-set null, and we therefore do not present it as a biomarker.")
    # 找到原来 "Three findings..." 那句的结束
    m2 = re.search(r"Three findings emerge\..*?(?=\n\n)", tail, re.S)
    if m2:
        tail = new_tail_head + tail[m2.end():]
    disc = head + tail
replace_sec("Discussion", disc)

open(OUT, "w", encoding="utf-8").write(s)
print(f"写入 {OUT}")

# 校验
v7 = open(OUT, encoding="utf-8").read()
print(f"\n标题: {v7.splitlines()[0][:110]}")
print(f"词数: {len(re.findall(r'[A-Za-z0-9][A-Za-z0-9-]*', v7))}")
print("\nResults 小节顺序:")
for m in re.finditer(r"^### (.+)$", v7[:v7.index("## Discussion")], re.M):
    print(f"  {m.group(1)}")
