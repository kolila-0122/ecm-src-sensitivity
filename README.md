# Code for: An extracellular matrix gene program is independently prognostic and marks a Src-inhibitor-sensitive state across cancer types

Analysis code accompanying the manuscript submitted to the *Journal of the Egyptian National Cancer Institute*.

This repository reproduces every number, table and figure in the manuscript directly from public source datasets. No intermediate data files are included; all data are downloaded by the scripts.

---

## What the study does

A 23-gene extracellular matrix (ECM) transcriptional programme is evaluated in three settings:

1. **Prognosis** — overall survival across six TCGA PanCancer Atlas cohorts (n = 3,255), by Cox regression and inverse-variance meta-analysis.
2. **Drug sensitivity** — 286 compounds in GDSC2 (969 cell lines), with independent replication in GDSC1 (378 compounds), and mechanistic follow-up (proliferation adjustment, cancer-type stratification, RPPA, GSEA).
3. **Treatment response** — four immunotherapy cohorts (n = 232) and four neoadjuvant chemotherapy cohorts (n = 849).

**Principal findings.** The ECM programme is independently prognostic (pooled adjusted HR 1.085 per SD, p = 0.011, I² = 0%). Of 286 compounds, 87 reach FDR < 0.05, of which **73 indicate resistance and 14 sensitivity**. The best-supported sensitivity signal is dasatinib (ρ = −0.386, n = 395), which replicates and strengthens in GDSC1 (ρ = −0.586). The programme does not predict immunotherapy or chemotherapy response.

---

## Requirements

```
python >= 3.10
numpy, pandas, scipy, scikit-learn
lifelines
gseapy            # GSEA
matplotlib        # figures
pymupdf           # figure assembly and PDF → TIFF
Pillow            # CMYK conversion
requests
```

Install:

```bash
pip install numpy pandas scipy scikit-learn lifelines gseapy matplotlib pymupdf pillow requests
```

---

## Data sources

All data are public and are fetched by the scripts; nothing is bundled.

| Resource | Access | Used by |
|---|---|---|
| TCGA PanCancer Atlas (clinical + mRNA z-scores) | cBioPortal REST API, `https://www.cbioportal.org/api` | `p26`, `p27`, `p28`, `p29`, `p30`, `p31`, `p33` |
| CCLE transcriptomes (1,156 lines) | cBioPortal, study `ccle_broad_2019` | `p31`, `p35`, `p36`, `p38`, `p40` |
| CCLE RPPA (protein) | cBioPortal, profile `ccle_broad_2019_protein_quantification_zscores` | `p38` |
| GDSC1 / GDSC2 fitted dose-response (release 8.4) | `https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/` | `p35`, `p36` |
| Immunotherapy cohorts (Riaz, Hugo, Liu, UCLA) | cBioPortal, melanoma studies | `p30` |
| Neoadjuvant chemotherapy cohorts (GSE20194, GSE20271, GSE22093, GSE25055) | NCBI GEO FTP, GPL96 platform | `p31`, `p32` |

**Note on GEO access.** During development the NCBI eutils endpoint intermittently returned HTTP 500. All GEO series matrices were therefore retrieved directly from the FTP mirror:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{first2}nnn/GSE#####/matrix/GSE#####_series_matrix.txt.gz
```

Platform annotation for GPL96 was obtained from the SOFT file at

```
https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/soft/GPL96_family.soft.gz
```

(the `.annot.gz` file is no longer served).

---

## Run order

Scripts write intermediate files to `cache/` and results to `results/`. Run in this order:

```bash
# 1. Prognosis
python p26_survival.py        # per-cohort univariable Cox
python p29_multivar.py        # adjusted Cox + meta-analysis  → Table S4, Figure 1
python p27_meta.py            # pooled estimates
python make_forest.py         # Figure 1

# 2. Drug sensitivity
python p35_drug.py            # GDSC2, 286 compounds          → Table S1
python p36_validate.py        # proliferation adjustment + cancer-type stratification
                              #                                → Tables S2, S3
python p50_revision.py        # rebuild tables in scientific notation
python make_panels_drug.py    # Figure 2A
python make_fig2b.py          # Figure 2B (GDSC2 vs GDSC1)

# 3. Mechanism
python p38_rppa.py            # CCLE RPPA: paxillin, FAK, Src
python p40_gsea.py            # Hallmark GSEA                  → Figure S4
python make_gsea_panel.py

# 4. Treatment response
python p30_ici.py             # immunotherapy                  → Table S5, Figure S3A
python p32_pcr.py             # neoadjuvant chemotherapy       → Table S6, Figure S3B
python make_panel_treatment.py

# 5. Interaction analysis (reported as negative)
python p33_interaction.py     # ECM × immune quadrants         → Figure S1
python p34_supptests.py       # LRT + C-index: interaction does not hold

# 6. Robustness
python p60_robustness.py      # gene-subset sensitivity        → Supplementary

# 7. Figure assembly
python make_assembled_figs.py
```

---

## Key parameter choices

| Choice | Value | Rationale |
|---|---|---|
| ECM score | unweighted mean of 23 genes | Pre-specified; no training. Robustness to gene subset tested in `p60` (HR 1.13–1.15 across five subsets) |
| Proliferation score | unweighted mean of 12 genes, no overlap with ECM set | Used only as a covariate in partial correlation |
| Multiple testing | Benjamini–Hochberg FDR | Applied within the 286-compound panel |
| Partial correlation | rank-residual method | Both variables rank-transformed, regressed on ranked proliferation, residuals correlated |
| GSEA | prerank, Hallmark 2020, 1,000 permutations, **seed = 1** | Seed fixed for reproducibility |
| Meta-analysis | inverse-variance weighting of log HR | Heterogeneity by Cochran's Q and I² |
| Immunotherapy pooling | Fisher z of AUC weighted by n | Endpoint-stratified results also reported |
| FDR < 0.05 threshold | 87 of 286 compounds | Of these, 73 resistance, 14 sensitivity |

---

## Reproducibility notes

- The scripts download data at run time; total download is approximately 300 MB (GDSC1 is the largest at 52 MB compressed).
- Random seed for GSEA permutations is fixed at 1.
- All p and FDR values in the manuscript are reported in scientific notation. Values below floating-point resolution appear as `p < 1e-15` in the manuscript text and `<1e-300` in Table S1; they are never reported as zero.
- Figure colour: all figures are exported as 600 dpi CMYK TIFF for print production.

---

## Known limitations

Stated plainly, because they bound what these results can support:

1. **The prognostic effect is small.** Pooled HR 1.085 per SD. Tertile stratification is not significant in any of the six cohorts individually, so the score is not usable for patient stratification.
2. **The drug screen measures 2D monolayer sensitivity.** CCLE lines have no stroma; the "ECM score" there reflects what cancer cells themselves transcribe, not stromal content. The dasatinib finding is a cell-autonomous pharmacological association, not evidence that stroma causes resistance.
3. **No in vivo validation.** The dasatinib prediction remains a hypothesis generated from cell lines.
4. **Most GDSC2 drug associations do not replicate.** Of the associations significant in GDSC2, most fail in GDSC1 and are withdrawn in the manuscript. Only the Src and MEK classes survive.
5. **The interaction analysis is negative.** The ECM × immune quadrant comparison (HR 1.50) does not survive formal interaction testing (LRT non-significant in 6/6 cohorts; C-index change 0.004) and is reported as a supplementary observation only.

---

## Contact

Bin Lian (corresponding author) — drmilo@gkd.edu.cn
School of Health Sciences, Guangzhou Vocational University of Science and Technology
