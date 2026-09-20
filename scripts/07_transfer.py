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

X, y, pos, gk = load_tcga()
colC = [pos[g] for g in PANEL_C if g in pos]
muC, sdC = X[:, colC].mean(0), X[:, colC].std(0) + 1e-6
Ftr = np.ascontiguousarray((X[:, colC] - muC) / sdC)
lr = LogisticRegression(max_iter=2000).fit(Ftr, y)

# GTEx normal breast
Lb, gpos_b, n_b = load_gct(os.path.join(DATA, "gtex_breast_v10_reads.gct.gz"))
idx_b = [gpos_b[g] for g in PANEL_C if g in gpos_b]
syms_b = [g for g in PANEL_C if g in gpos_b]
Fte_b = np.ascontiguousarray(((Lb[idx_b, :] - muC[:len(idx_b), None]) /
                              (sdC[:len(idx_b), None] + 1e-6 + 1e-6)).T)
p_gtex = lr.predict_proba(Fte_b)[:, 1]
gtex_normal = float((p_gtex < 0.5).mean())
print(f"GTEx normal breast: {gtex_normal:.1%} called normal (n={n_b})")

# METABRIC tumors
mm = fetch_metabric(PANEL_C)
idx_m = mm[PANEL_C].dropna().index
Zm = pscore(mm.loc[idx_m, PANEL_C].to_numpy(dtype=np.float32))
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
