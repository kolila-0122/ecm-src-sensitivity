"""用 PyMuPDF 把 panel 拼成期刊版式的 figure"""
import os, sys
import pymupdf

sys.stdout.reconfigure(line_buffering=True)
FIG = r"D:\stiffness_r3\figs\panels"
OUTDIR = r"D:\stiffness_r3\figs\assembled"
os.makedirs(OUTDIR, exist_ok=True)

MM = 72 / 25.4          # pt per mm


def place(page, src_path, rect):
    """把 src PDF 的第 1 页按比例缩放进 rect（保持宽高比，居中）"""
    doc = pymupdf.open(src_path)
    sr = doc[0].rect
    scale = min(rect.width / sr.width, rect.height / sr.height)
    w, h = sr.width * scale, sr.height * scale
    x0 = rect.x0 + (rect.width - w) / 2
    y0 = rect.y0 + (rect.height - h) / 2
    page.show_pdf_page(pymupdf.Rect(x0, y0, x0 + w, y0 + h), doc, 0)
    doc.close()


def label(page, xy, text, size=8):
    page.insert_text(xy, text, fontsize=size, fontname="hebo", color=(0, 0, 0))


# ============ Figure 1：预后森林图（单栏） ============
doc = pymupdf.open()
p = doc.new_page(width=110 * MM, height=80 * MM)
place(p, os.path.join(FIG, "PE1_forest.pdf"), pymupdf.Rect(0, 0, 110 * MM, 80 * MM))
label(p, (4 * MM, 7 * MM), "A", 10)
doc.save(os.path.join(OUTDIR, "Figure1_prognosis.pdf"))
doc.close()
print("Figure1 saved")

# ============ Figure 2：药敏（火山图 + 散点），双栏 ============
W, H = 180 * MM, 78 * MM
doc = pymupdf.open()
p = doc.new_page(width=W, height=H)
place(p, os.path.join(FIG, "PD1_volcano.pdf"), pymupdf.Rect(0, 0, 88 * MM, H))
place(p, os.path.join(FIG, "PD2b_dasatinib_two_resources.pdf"), pymupdf.Rect(90 * MM, 0, 178 * MM, H))
label(p, (3 * MM, 6 * MM), "A", 10)
label(p, (93 * MM, 6 * MM), "B", 10)
doc.save(os.path.join(OUTDIR, "Figure2_drug_sensitivity.pdf"))
doc.close()
print("Figure2 saved (A volcano + B two-resource replication)")

# ============ Figure S2：癌种分层热图（单栏） ============
doc = pymupdf.open()
p = doc.new_page(width=90 * MM, height=95 * MM)
place(p, os.path.join(FIG, "PD3_bycancer_heatmap.pdf"), pymupdf.Rect(0, 0, 90 * MM, 95 * MM))
label(p, (3 * MM, 6 * MM), "A", 10)
doc.save(os.path.join(OUTDIR, "FigureS2_bycancer_heatmap.pdf"))
doc.close()
print("FigureS2 saved")

# ============ Figure S3：治疗反应阴性（双栏） ============
doc = pymupdf.open()
p = doc.new_page(width=180 * MM, height=88 * MM)
place(p, os.path.join(FIG, "PT1_treatment_response.pdf"), pymupdf.Rect(0, 0, 180 * MM, 88 * MM))
label(p, (3 * MM, 6 * MM), "A", 10)
doc.save(os.path.join(OUTDIR, "FigureS3_treatment_response.pdf"))
doc.close()
print("FigureS3 saved")

# ============ Figure 3：泛化失败（methodology），双栏 ============
W, H = 180 * MM, 90 * MM
doc = pymupdf.open()
p = doc.new_page(width=W, height=H)
place(p, os.path.join(FIG, "PF1_scale_concept.pdf"), pymupdf.Rect(0, 0, 88 * MM, H))
place(p, os.path.join(FIG, "PA1_proxy_comparison.pdf"), pymupdf.Rect(90 * MM, 0, 178 * MM, 45 * MM))
place(p, os.path.join(FIG, "PA1_generalization.pdf"), pymupdf.Rect(90 * MM, 45 * MM, 178 * MM, H))
label(p, (3 * MM, 6 * MM), "A", 10)
label(p, (93 * MM, 6 * MM), "B", 10)
label(p, (93 * MM, 51 * MM), "C", 10)
doc.save(os.path.join(OUTDIR, "FigureS4_generalization.pdf"))
doc.close()
print("Figure3 saved")

# ============ Figure 4：ECM×免疫（含伪影说明） ============
W, H = 180 * MM, 92 * MM
doc = pymupdf.open()
p = doc.new_page(width=W, height=H)
place(p, os.path.join(FIG, "PI1_km_four_groups.pdf"), pymupdf.Rect(0, 0, 88 * MM, H))
place(p, os.path.join(FIG, "PI2_forest.pdf"), pymupdf.Rect(90 * MM, 0, 178 * MM, H))
label(p, (3 * MM, 6 * MM), "A", 10)
label(p, (93 * MM, 6 * MM), "B", 10)
doc.save(os.path.join(OUTDIR, "FigureS1_ecm_immune.pdf"))
doc.close()
print("Figure4 saved")

print("\n=== 输出 ===")
for f in sorted(os.listdir(OUTDIR)):
    sz = os.path.getsize(os.path.join(OUTDIR, f)) / 1024
    print(f"  {f}  {sz:.0f} KB")
