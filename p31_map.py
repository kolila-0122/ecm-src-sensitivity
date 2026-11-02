"""从 GPL96 soft 提取 probe -> gene symbol 映射"""
import gzip, sys
import pandas as pd

sys.stdout.reconfigure(line_buffering=True)
P = r"D:\stiffness_r3\data\GEO_neoadj\GPL96.soft.gz"
OUT = r"D:\stiffness_r3\data\GEO_neoadj\GPL96_map.csv"

rows = []
in_tab = False
hdr = None
with gzip.open(P, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            in_tab = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if not in_tab:
            continue
        parts = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = parts
            print(f"表头 {len(hdr)} 列: {hdr[:6]}")
            try:
                i_id = hdr.index("ID")
                i_sym = hdr.index("Gene Symbol")
            except ValueError:
                i_id, i_sym = 0, 10
            continue
        if len(parts) <= max(i_id, i_sym):
            continue
        gid, sym = parts[i_id], parts[i_sym].strip()
        if sym and sym != "---":
            # 多基因用 /// 分隔，取第一个
            s = sym.split("///")[0].strip()
            rows.append((gid, s))

df = pd.DataFrame(rows, columns=["probe", "symbol"])
df.to_csv(OUT, index=False)
print(f"提取 {len(df)} 条映射，{df.symbol.nunique()} 个基因")
print(df.head())

ECM = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
       "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CCN2", "ACTA2", "TAGLN",
       "TIMP1", "SPARC", "DCN", "LUM", "TGFB1", "COL11A1", "COL8A1", "FBN1"]
found = df[df.symbol.isin(ECM)]
print(f"\nECM 基因命中 {found.symbol.nunique()}/{len(ECM)}")
print(found.groupby("symbol").probe.apply(list).to_string())
