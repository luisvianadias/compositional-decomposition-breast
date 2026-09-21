# -*- coding: utf-8 -*-
"""Section 2.4 — TCGA-trained panel applied to independent cohorts:
GTEx v10 normal breast (99%+ called normal) and METABRIC tumors (92%+).
Output: results/07_transfer.json
"""
import os, sys, json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DATA, load_tcga, load_gct, fetch_metabric, pscore,
                    PANEL_C)

X, y, pos, _ = load_tcga()
syms = [g for g in PANEL_C if g in pos]
if set(syms) != set(PANEL_C):
    raise SystemExit(f"E15-GUARD: ausentes no TCGA: "
                     f"{sorted(set(PANEL_C) - set(syms))}")
colC = [pos[g] for g in syms]
# per-sample z (platform-invariant) on BOTH sides of the transfer:
# train and test features must live in the same space
Ztc = pscore(X[:, colC])
lr = LogisticRegression(max_iter=2000).fit(Ztc, y)

# GTEx normal breast
Lb, gpos_b, n_b = load_gct(os.path.join(DATA, "gtex_breast_v10_reads.gct.gz"))
syms_b = [g for g in syms if g in gpos_b]
if syms_b != syms:
    raise SystemExit(f"E15-GUARD: genes ausentes no GTEx: "
                     f"{sorted(set(syms) - set(syms_b))}")
Fte_b = pscore(Lb[[gpos_b[g] for g in syms_b], :].T.astype(np.float32))
p_gtex = lr.predict_proba(Fte_b)[:, 1]
gtex_normal = float((p_gtex < 0.5).mean())
print(f"GTEx normal breast: {gtex_normal:.1%} called normal (n={n_b})")

# METABRIC tumors (per-sample z, same gene order; loud guard)
mm = fetch_metabric(syms)
missing = [g for g in syms if g not in mm.columns]
if missing:
    raise SystemExit(f"E15-GUARD: genes ausentes no METABRIC: {missing}")
idx_m = mm[syms].dropna().index
Zm = pscore(mm.loc[idx_m, syms].to_numpy(dtype=np.float32))
p_met = lr.predict_proba(Zm)[:, 1]
met_tumor = float((p_met > 0.5).mean())
print(f"METABRIC tumors: {met_tumor:.1%} called tumor (n={len(p_met)})")

res = {"gtex_normal_frac": gtex_normal, "gtex_n": n_b,
       "metabric_frac_tumor": met_tumor,
       "metabric_n": int(len(p_met)),
       "panel": PANEL_C}
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results", "07_transfer.json"), "w"),
          indent=2)
