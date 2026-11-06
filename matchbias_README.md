# matchbias

**靶点类目关键词匹配偏倚自检工具**

检测化合物-靶点注释中，无词边界的正则匹配是否把含目标子串的其他基因误归入某类目（如 `VEGFR` ⊃ `EGFR`、`PERK` ⊃ `ERK`），并列出被误归的化合物。

---

## 为什么需要它

药敏研究常用「按靶点归类 → 富集检验」的流程。归类多用正则匹配：

```python
re.search("EGFR|ERBB|HER2", target_text)      # 无词边界（默认写法）
```

这个写法会把 `VEGFR`、`VEGFR2` 一并命中——4 种 VEGFR 抑制剂被误归入 EGFR/ERBB 类。

**后果**：类目变大、富集变强。同一批数据，加一个 `\b` 就能让 P 值从 0.0026 变到 0.25。

本工具把这个问题变成可自动检测的检查项。

---

## 安装

单文件，无依赖（仅标准库）：

```bash
git clone https://github.com/kolila-0122/ecm-src-sensitivity.git
cd ecm-src-sensitivity/matchbias
python matchbias.py --help
```

---

## 用法

```bash
python matchbias.py --input compounds.csv --target-col Target --name-col Drug
```

带通路列：

```bash
python matchbias.py --input compounds.csv --target-col Target --pathway-col Pathway
```

自定义规则（TSV：`类目名 \t 基因1;基因2;...`）：

```bash
python matchbias.py --input compounds.csv --rules my_rules.tsv
```

---

## 输入格式

CSV 或 TSV，至少两列：

| Drug | Target | Pathway |
|---|---|---|
| Sorafenib | PDGFR, KIT, VEGFR, RAF | RTK signaling |
| Erlotinib | EGFR | EGFR signaling |

- `Drug`：化合物名
- `Target`：靶点注释（可含多个，逗号或分号分隔）
- `Pathway`：通路名（可选）

---

## 输出

```
============================================================
matchbias — 靶点类目关键词匹配偏倚自检
============================================================
化合物总数    : 378
类目数        : 13
受影响类目    : 3

类目                       无锚定     词边界      差异     被误归
------------------------------------------------------------
EGFR/ERBB                 24      13     +11      11  ←
Other                    219     229     -10       0  ←
Src family                 8       9      -1       0  ←

============================================================
被误归化合物
============================================================

【EGFR/ERBB】共 11 种：
    - Axitinib
    - Cabozantinib
    - Sorafenib
    - Sunitinib
    - Tivozanib
    ...
```

**退出码**：`0` = 无偏倚；`1` = 检出偏倚（可用于 CI 卡点）。

---

## 实测结果

| 数据集 | 化合物数 | 无锚定 EGFR/ERBB | 词边界 | 被误归 |
|---|---|---|---|---|
| **GDSC1** | 378 | 24 | 13 | **11** |
| **GDSC2** | 286 | 12 | 8 | **4** |

被误归的全部为 VEGFR 或多靶点酪氨酸激酶抑制剂。

---

## 默认类目规则

`AURK/mitotic`、`Tubulin`、`BRAF`、`CDK`、`MEK/MAPK`、`PI3K/AKT/mTOR`、`EGFR/ERBB`、`WNT`、`Src family`、`HDAC`、`IGF1R`、`Apoptosis`，其余归 `Other`。

规则按顺序匹配，首个命中即返回。可用 `--rules` 替换。

---

## 匹配模式

**无锚定**（问题写法）：

```python
re.compile("EGFR")
```

**词边界**（推荐）：

```python
re.compile(r"(?<![A-Z0-9])EGFR(?:[0-9]|[A-Z](?![A-Z0-9]))?")
```

- `(?<![A-Z0-9])`：词首不能是字母或数字 → 排除 `VEGFR`
- `(?:[0-9]|[A-Z](?![A-Z0-9]))?`：允许一个数字或单字母后缀 → 保留 `EGFR2`、`EGFRB`

---

## 建议

在方法部分报告匹配规则，并在补充材料公开类目成员。若本工具检出偏倚，做一次双口径敏感性分析。

---

## 引用

郑家惟, 连彬, 黄诗茵. 关键词匹配边界效应高估靶点类目富集[J]. 生物工程学报, 2026（投稿中）.

---

## 许可

MIT
