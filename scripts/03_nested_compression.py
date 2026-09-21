# -*- coding: utf-8 -*-
"""Section 2.1-2.2 — nested compressibility, pre-specified panels,
effect-size-matched (randDE) and uniformly random (null) controls.

All estimation inside training folds. 15 paired splits (5-fold x 3).
Logistic regression only. Output: results/03_nested_compression.json
"""
import os, sys, json, time
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_tcga, pscore, PANEL_B, PANEL_C
from pipeline_utils import pick_nearest_non_panel

t0 = time.time()
DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "..", "results")
os.makedirs(OUT, exist_ok=True)

X, y, pos, gk = load_tcga()
N_GENES = X.shape[1]
abs_lfc = np.abs(X[y == 1].mean(0) - X[y == 0].mean(0))
rank_global = np.argsort(-abs_lfc)

PANEL_A = PANEL_B[:12]
KS = [1, 2, 3, 5, 10, 20]
N_RAND_MATCHED = 5
N_NULL = 500
NULL_K = 10

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
splits = list(rskf.split(X, y))
n_tr, n_te = len(splits[0][0]), len(splits[0][1])

def lr_acc(Ftr, ytr, Fte, yte):
    m = LogisticRegression(max_iter=2000).fit(Ftr, ytr)
    return float(m.score(Fte, yte)), \
           float(roc_auc_score(yte, m.predict_proba(Fte)[:, 1]))

arms = {n: {"acc": np.zeros(len(splits)), "auc": np.zeros(len(splits))}
        for n in ["full", "panelA", "panelB", "panelC", "randDE"]
        + [f"k{k}" for k in KS] + ["null"]}
null_draws_acc = np.zeros(N_NULL)
aucs_null = np.zeros(N_NULL)

for i, (tr, te) in enumerate(splits):
    Xtr, Xte, ytr, yte = X[tr], X[te], y[tr], y[te]
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    lfc = np.abs(Ztr[ytr == 1].mean(0) - Ztr[ytr == 0].mean(0))
    rk = np.argsort(-lfc)

    for k in KS:
        a, u = lr_acc(Ztr[:, rk[:k]], ytr, Zte[:, rk[:k]], yte)
        arms[f"k{k}"]["acc"][i], arms[f"k{k}"]["auc"][i] = a, u
    a, u = lr_acc(Ztr, ytr, Zte, yte)
    arms["full"]["acc"][i], arms["full"]["auc"][i] = a, u

    # panels (pre-specified gene sets, evaluated on raw panel expression)
    colB = np.array([pos[g] for g in PANEL_B])
    colC = np.array([pos[g] for g in PANEL_C])
    colA = colB[:12]
    for nm, cc in (("panelA", colA), ("panelB", colB), ("panelC", colC)):
        a, u = lr_acc(Xtr[:, cc], ytr, Xte[:, cc], yte)
        arms[nm]["acc"][i], arms[nm]["auc"][i] = a, u

    # randDE: non-adipocyte gene with nearest training-fold |log2FC|
    # (E15-GUARD: matching por símbolo — gk é Ensembl; ver pipeline_utils)
    fat = set(PANEL_B)
    sym_of = {}
    for s, ci in pos.items():
        sym_of.setdefault(ci, s)
    rd = pick_nearest_non_panel(lfc, colB, fat, sym_of)
    a, u = lr_acc(Xtr[:, rd], ytr, Xte[:, rd], yte)
    arms["randDE"]["acc"][i] = a
    arms["randDE"]["auc"][i] = u

    # null: uniform random 10-gene panels (median over B draws per split)
    rng = np.random.default_rng(42)
    draws = np.zeros(N_NULL)
    for b in range(N_NULL):
        cols = rng.choice(N_GENES, size=NULL_K, replace=False)
        m = LogisticRegression(max_iter=2000).fit(Ztr[:, cols], ytr)
        draws[b] = m.score(Zte[:, cols], yte)
        aucs_null[b] = roc_auc_score(yte, m.predict_proba(Zte[:, cols])[:, 1])
    arms["null"]["acc"][i] = float(draws.mean())
    arms["null"]["auc"][i] = float(aucs_null.mean())
    null_draws_acc += draws / len(splits)

    if (i + 1) % 5 == 0:
        print(f"  split {i+1}/{len(splits)}", flush=True)

print(f"\n=== nested results (15 paired splits) ===")
for n, v in arms.items():
    print(f"  {n:<10} acc={v['acc'].mean():.1%}+-{v['acc'].std():.1%} "
          f"auc={v['auc'].mean():.3f}")

res = {n: {"acc": float(v["acc"].mean()), "acc_std": float(v["acc"].std()),
           "auc": float(v["auc"].mean()),
           "per_split_acc": [float(x) for x in v["acc"]]}
       for n, v in arms.items()}
res["null_per_draw_overall_mean"] = float(null_draws_acc.mean())
res["null_acc_per_draw"] = [float(x) for x in null_draws_acc]
res["null_auc_per_draw_mean_iqr"] = [float(np.quantile(aucs_null, q))
                                     for q in (0.25, 0.5, 0.75)]
json.dump(res, open(os.path.join(OUT, "nested_compression.json"), "w"),
          indent=2)
print(f"\nsalvo nested_compression.json ({time.time()-t0:.0f}s)")
