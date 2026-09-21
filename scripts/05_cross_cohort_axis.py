# -*- coding: utf-8 -*-
"""Section 2.4 — adipose-axis score across populations (pooled z, pure
16-gene adipocyte panel; positive control: GTEx subcutaneous adipose).
Output: results/05_cross_cohort_axis.json
"""
import os, sys, json
import numpy as np
from common import DATA, load_tcga, load_gct, PURE16

X, y, pos, _ = load_tcga()
gtex_breast = os.path.join(DATA, "gtex_breast_v10_reads.gct.gz")
gtex_adip = os.path.join(DATA, "gtex_adipose_subcut_v10_reads.gct.gz")

pure = [g for g in PURE16 if g in pos]
rows_t = X[:, [pos[g] for g in pure]]           # samples x panel (log2CPM)

def panel_matrix(path):
    L, gpos, n = load_gct(path)
    idx = [gpos[g] for g in pure if g in gpos]
    return L[idx, :], n

rows_b, n_b = panel_matrix(gtex_breast)
rows_a, n_a = panel_matrix(gtex_adip)

pooled = np.concatenate([rows_t.T[:, y == 1], rows_t.T[:, y == 0],
                         rows_b, rows_a], axis=1)
mu = pooled.mean(axis=1, keepdims=True)
sd = pooled.std(axis=1, keepdims=True) + 1e-6
zp = (pooled - mu) / sd
n1 = int((y == 1).sum()); n2 = int((y == 0).sum())
n3 = rows_b.shape[1]; n4 = rows_a.shape[1]

def med(seg):
    return float(np.median(zp[:, seg].mean(axis=0)))

res = {
    "tcga_tumor": med(slice(0, n1)),
    "tcga_adjacent_normal": med(slice(n1, n1 + n2)),
    "gtex_breast": med(slice(n1 + n2, n1 + n2 + n3)),
    "gtex_adipose": med(slice(n1 + n2 + n3, None)),
    "n": {"tcga_tumor": n1, "tcga_normal": n2, "gtex_breast": n3,
          "gtex_adipose": n4},
}
print(json.dumps({k: round(v, 3) for k, v in res.items() if k != "n"},
                 indent=1))
print("n:", res["n"])
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results",
                                 "05_cross_cohort_axis.json"), "w"), indent=2)
