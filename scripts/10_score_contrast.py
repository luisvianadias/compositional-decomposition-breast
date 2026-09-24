# -*- coding: utf-8 -*-
"""Section 2.3 — paired-score specificity contrast for the global
ablation. The adipose score A is diffuse (Test 2b: removing the input
dimensions most correlated with A costs nothing), so the correct
specificity instrument at the score level: global ablation with A vs.
global ablation with B randomized control scores (15 non-adipose genes
|Δ|-matched to the panel genes per training fold, randomized among the
five nearest neighbours, without replacement).

D(A) >> D(ctl)  => adipose-specific drop
D(A) ~ D(ctl)   => generic route dependence (the drop does not certify
                   adipose identity; specificity rests on Sections 2.4-2.6)
Output: results/10_score_contrast.json
"""
import os, sys, json
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_tcga, PANEL_B

t0 = __import__("time").time()
DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "..", "results", "10_score_contrast.json")
B_RAND = 30
INPUT_TOP = 2000
NEIGHBORS = 5
FULL_ANCHOR = 0.99104   # 04_adipose_ablation.json::before/acc

X, y, pos, _ = load_tcga()
colB = np.array([pos[g] for g in PANEL_B if g in pos])
prev = float(np.mean(y))

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
splits = list(rskf.split(X, y))
n_te, n_tr = len(splits[0][1]), len(splits[0][0])

accA, accR = [], []
for i, (tr, te) in enumerate(splits):
    mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-6
    Ztr, Zte = (X[tr] - mu) / sd, (X[te] - mu) / sd
    muA, sdA = X[tr][:, colB].mean(0), X[tr][:, colB].std(0) + 1e-6
    A_tr = ((X[tr][:, colB] - muA) / sdA).mean(1)
    A_te = ((X[te][:, colB] - muA) / sdA).mean(1)
    dm = np.abs(Ztr[y[tr] == 1].mean(0) - Ztr[y[tr] == 0].mean(0))
    inp = np.argsort(-dm)[:INPUT_TOP]
    rng = np.random.default_rng(1000 + i)

    def abl(score_tr, score_te):
        Ac = score_tr - score_tr.mean()
        beta = (Ztr[:, inp].T @ Ac) / (Ac @ Ac)
        Zg_tr = Ztr[:, inp] - np.outer(Ac, beta)
        Zg_te = Zte[:, inp] - np.outer(score_te - score_tr.mean(), beta)
        m = LogisticRegression(max_iter=2000).fit(
            np.ascontiguousarray(Zg_tr), y[tr])
        return m.score(np.ascontiguousarray(Zg_te), y[te])

    accA.append(abl(A_tr, A_te))
    for b in range(B_RAND):
        alt = []
        pool = np.setdiff1d(np.arange(Ztr.shape[1]), colB)
        for g in colB:
            d = np.abs(dm[pool] - dm[g])
            k = int(rng.integers(0, min(NEIGHBORS, len(pool))))
            j = int(np.argsort(d, kind="stable")[k])
            alt.append(pool[j])
            pool = np.delete(pool, j)
        alt = np.array(alt)
        muR, sdR = X[tr][:, alt].mean(0), X[tr][:, alt].std(0) + 1e-6
        R_tr = ((X[tr][:, alt] - muR) / sdR).mean(1)
        R_te = ((X[te][:, alt] - muR) / sdR).mean(1)
        accR.append(abl(R_tr, R_te))
    print(f"  split {i+1:2d}/15: A={accA[-1]:.3f} "
          f"ctl_med={np.median(accR[-B_RAND:]):.3f}", flush=True)

accR_by_split = np.array(accR).reshape(len(splits), B_RAND)
dA = float(np.mean(accA))
dR_mean = float(np.mean(accR_by_split))
dR_dist = accR_by_split.mean(axis=0)
D_A = (FULL_ANCHOR - dA) / (FULL_ANCHOR - prev)
D_R = (FULL_ANCHOR - dR_mean) / (FULL_ANCHOR - prev)
d = accA - np.median(accR_by_split, axis=1)
t_nb = float(d.mean() / np.sqrt(d.var(ddof=1) * (1 / len(d) + n_te / n_tr)
                                + 1e-18))

res = {
    "acc_ablacao_A_adiposo": dA,
    "acc_ablacao_ctl_media": dR_mean,
    "acc_ablacao_ctl_por_draw_iqr": [float(np.percentile(dR_dist, 25)),
                                     float(np.percentile(dR_dist, 75))],
    "ctl_distintos": int(len(np.unique(np.round(dR_dist, 10)))),
    "prevalencia": prev,
    "D_A_adiposo": D_A, "D_ctl": D_R,
    "EXCESSO_especificidade": D_A - D_R,
    "pareado_NB_t_adiposo_vs_ctl": t_nb,
    "B_rand": B_RAND, "input_top": INPUT_TOP, "neighbors": NEIGHBORS,
    "nota": "Contraste no nivel do escore (regime difuso): ablaglobal "
            "com A adiposo vs. escores controle pareados por |Delta-media|, "
            "sorteados entre os 5 vizinhos mais proximos (sem reposicao).",
}
json.dump(res, open(OUT, "w", encoding="utf-8"), indent=2)
print(f"\nD_A={D_A:+.2f}  D_ctl={D_R:+.2f}  "
      f"EXCESSO={D_A - D_R:+.2f}  t={t_nb:+.2f}")
print(f"salvo 10_score_contrast.json ({__import__('time').time()-t0:.0f}s)")
