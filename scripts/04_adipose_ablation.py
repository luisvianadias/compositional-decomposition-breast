# -*- coding: utf-8 -*-
import os, sys, json
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_tcga, PANEL_B

X, y, pos, gk = load_tcga()
colB = np.array([pos[g] for g in PANEL_B if g in pos])

# exatamente como e13: per-gene z (train params), A do painel por-gene-z
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
splits = list(rskf.split(X, y))
acc_b, auc_b, acc_a, auc_a = [], [], [], []
for tr, te in splits:
    mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-6
    Ztr, Zte = (X[tr]-mu)/sd, (X[te]-mu)/sd
    muA, sdA = X[tr][:, colB].mean(0), X[tr][:, colB].std(0) + 1e-6
    A_tr = ((X[tr][:, colB]-muA)/sdA).mean(1)
    A_te = ((X[te][:, colB]-muA)/sdA).mean(1)
    Ac = A_tr - A_tr.mean()
    beta = (Ztr.T @ Ac) / (Ac @ Ac)
    Ztr_r = np.ascontiguousarray(Ztr - np.outer(Ac, beta))
    Zte_r = np.ascontiguousarray(Zte - np.outer(A_te - A_tr.mean(), beta))
    cc = np.argsort(np.abs(
        Ztr[y[tr] == 1].mean(0) - Ztr[y[tr] == 0].mean(0)))[::-1][:2000]
    m = LogisticRegression(max_iter=2000).fit(
        np.ascontiguousarray(Ztr_r[:, cc]), y[tr])
    acc_a.append(m.score(np.ascontiguousarray(Zte_r[:, cc]), y[te]))
    auc_a.append(roc_auc_score(y[te], m.predict_proba(
        np.ascontiguousarray(Zte_r[:, cc]))[:, 1]))
    m0 = LogisticRegression(max_iter=2000).fit(
        np.ascontiguousarray(Ztr[:, cc]), y[tr])
    acc_b.append(m0.score(np.ascontiguousarray(Zte[:, cc]), y[te]))
    auc_b.append(roc_auc_score(y[te], m0.predict_proba(
        np.ascontiguousarray(Zte[:, cc]))[:, 1]))

p1 = PCA(n_components=1, random_state=0).fit(Ztr)
pc1_full = p1.transform(Ztr)[:, 0]
# PC1 direction is sign-indeterminate (SVD): orient the AUC upward
auc_full = roc_auc_score(y[tr], pc1_full)
pc1_auc_full = max(auc_full, 1.0 - auc_full)
sep_full = abs(pc1_full[y[tr] == 1].mean() - pc1_full[y[tr] == 0].mean()) / \
           pc1_full.std()
p1r = PCA(n_components=1, random_state=0).fit(Ztr_r)
pc1_r = p1r.transform(Ztr_r)[:, 0]
auc_r = roc_auc_score(y[tr], pc1_r)
pc1_auc_r = max(auc_r, 1.0 - auc_r)
sep_r = abs(pc1_r[y[tr] == 1].mean() - pc1_r[y[tr] == 0].mean()) / pc1_r.std()

res = {
    "before": {"acc": float(np.mean(acc_b)), "auc": float(np.mean(auc_b))},
    "after": {"acc": float(np.mean(acc_a)), "auc": float(np.mean(auc_a))},
    "pc1_original": {"var": float(p1.explained_variance_ratio_[0]),
                     "auc": float(pc1_auc_full), "sep_sd": float(sep_full)},
    "pc1_residual": {"var": float(p1r.explained_variance_ratio_[0]),
                     "auc": float(pc1_auc_r), "sep_sd": float(sep_r)},
}
print(f"antes: {res['before']['acc']:.1%}/{res['before']['auc']:.3f} | "
      f"depois: {res['after']['acc']:.1%}/{res['after']['auc']:.3f}")
print(f"PC1 sep: {sep_full:.2f} -> {sep_r:.2f} sd")
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results",
                                 "04_adipose_ablation.json"), "w"), indent=2)
