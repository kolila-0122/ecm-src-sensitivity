# 关键词匹配边界效应高估靶点类目富集

**Keyword-matching boundary effects inflate target-class enrichment in drug-sensitivity analysis**

《生物工程学报》投稿代码仓库 / Analysis code for the manuscript submitted to *Chinese Journal of Biotechnology*.

本仓库重现稿件正文、表 1–3、图 1–9 与补充表 S1–S16 的全部分析。所有数据从公开源下载，脚本不含中间数据文件。

---

## 研究做了什么

以 23 基因细胞外基质（ECM）评分为载体，检验药敏关联研究中**靶点类目关键词匹配的词边界敏感性**：

1. **双口径类目富集** — 13 个靶点类目，无词边界 vs 词边界两种归类方式（表 1、表 S13）。
2. **跨资源复测** — 8 种化合物在 GDSC1 中复测（表 3、图 3）。
3. **药敏模块聚类** — 286 种化合物按 IC₅₀ 相似性聚类，3 个模块与 ECM 显著相关（图 8、图 9、表 S14）。
4. **随机基因集零分布** — 500 组 23 基因集，ECM 评分位于第 88 百分位（表 S10、表 S15）。

**核心发现**：无词边界匹配 `EGFR|ERBB|HER2` 会把 4 种 VEGFR 抑制剂（索拉非尼、莫特塞尼、阿昔替尼、西地尼布）扫入 EGFR/ERBB 类，使一个经 Bonferroni 校正后不显著的富集（P = 0.25）变成校正后显著（P = 0.023）——而底层药理学数据完全未变。

---

## 脚本清单

| 脚本 | 作用 | 产出 |
|---|---|---|
| `p105_target_class_v2.py` | **锚定/无锚定分类器**（1.4 节） | 表 1、表 S13 |
| `p105_target_class.py` | 无锚定分类器（初版） | — |
| `p106_cluster_modules.py` | **层次聚类 + 置换检验**（2.7 节） | 表 S14、图 8、图 9 |
| `p103_null_test_v2.py` | **随机基因集零分布**（1.6 节） | 表 S10、表 S15 |
| `p101_prolif_all.py` | 增殖校正 | 表 S8 |
| `p102_tcga_prolif.py` | TCGA 增殖评分 | 表 S9 |
| `make_*.py` | 图件生成 | 图 1–9 |

---

## 运行环境

```
python >= 3.10
numpy, pandas, scipy, scikit-learn, lifelines, matplotlib
```

```bash
pip install -r requirements.txt
```

---

## 数据来源

- **GDSC2 / GDSC1**：Sanger FTP（https://www.sanger.ac.uk/）
- **CCLE**：cBioPortal `ccle_broad_2019`
- **TCGA PanCancer Atlas**：cBioPortal API

脚本首次运行会从上述源下载并缓存到 `cache/`。

---

## 引用

郑家惟, 连彬, 黄诗茵. 关键词匹配边界效应高估靶点类目富集[J]. 生物工程学报, 2026（投稿中）.

---

## 许可

MIT License，见 `LICENSE`。
