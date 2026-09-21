# -*- coding: utf-8 -*-
"""Builds data/tcga_brca_counts.npz from the downloaded GDC TSVs."""
import os, sys, time
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
SRC = os.path.join(DATA, "gdc_brca")

files = sorted(f for f in os.listdir(SRC) if f.endswith(".tsv"))
print(f"{len(files)} amostras", flush=True)
t0 = time.time()

cols, barcodes, types = [], [], []
genes_ref = None
for i, fn in enumerate(files):
    df = pd.read_csv(os.path.join(SRC, fn), sep="\t", comment="#",
                     usecols=["gene_id", "unstranded"],
                     dtype={"unstranded": np.int32})
    df = df[df["gene_id"].str.startswith("ENSG")].set_index("gene_id")
    if genes_ref is None:
        genes_ref = df.index
        col = df["unstranded"].to_numpy()
    else:
        col = df.reindex(genes_ref)["unstranded"].to_numpy()
    cols.append(col)
    b = fn[:-4]
    barcodes.append(b)
    types.append("tumor" if "-01" in b[12:16] else "normal")
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(files)} ({time.time()-t0:.0f}s)", flush=True)

X = np.ascontiguousarray(np.stack(cols, axis=0).T, dtype=np.float32)
np.savez_compressed(os.path.join(DATA, "tcga_brca_counts.npz"),
                    counts=X, genes=genes_ref.to_numpy())
pd.DataFrame({"barcode": barcodes, "type": types}).to_csv(
    os.path.join(DATA, "tcga_brca_samples.tsv"), sep="\t", index=False)
print(f"matriz: {X.shape[0]} genes x {X.shape[1]} amostras -> "
      f"tcga_brca_counts.npz")
