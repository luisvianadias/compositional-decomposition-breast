# -*- coding: utf-8 -*-
"""Section 2.6 — ESTIMATE baseline (faithful port of the reference R
implementation; official gene sets) vs the adipose-associated panel.
Output: results/06_estimate_baseline.json
"""
import os, sys, json, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA, load_tcga, CORE10, pscore

t0 = time.time()
data_dir = DATA
pkg = os.path.join(data_dir, "estimate_pkg", "estimate", "inst", "extdata")
gmt = os.path.join(pkg, "SI_geneset.gmt")
signatures = {}
for line in open(gmt, encoding="utf-8", errors="replace"):
    p = line.rstrip("\n").split("\t")
    signatures[p[0]] = [g for g in p[2:] if g]
common = pd.read_csv(os.path.join(pkg, "common_genes.txt"), sep="\t")[
    "GeneSymbol"].dropna().tolist()

d = np.load(os.path.join(data_dir, "tcga_brca_counts.npz"), allow_pickle=True)
counts, genes = d["counts"], d["genes"].astype(str)
samples = pd.read_csv(os.path.join(data_dir, "tcga_brca_samples.tsv"),
                      sep="\t")
y = (samples["type"] == "tumor").to_numpy().astype(int)
lib = counts.sum(axis=0, keepdims=True)
lcpm = np.log2(counts / lib * 1e6 + 1.0)
gdc = os.path.join(data_dir, "gdc_brca", sorted(
    os.listdir(os.path.join(data_dir, "gdc_brca")))[0])
tsv = pd.read_csv(gdc, sep="\t", comment="#",
                  usecols=["gene_id", "gene_name"]).dropna()
tsv["gene_id"] = tsv["gene_id"].str.split(".").str[0]
smap = dict(tsv.drop_duplicates("gene_id").set_index("gene_id")["gene_name"])
sym_all = np.array([smap.get(g.split(".")[0], "") for g in genes])
sym_pos = {}
for i, s in enumerate(sym_all):
    if s and s not in sym_pos:
        sym_pos[s] = i

mask = np.isin(sym_all, common)
Xg = lcpm[mask]
gsyms = np.array(sorted(sym_all[mask]))
ng = Xg.shape[0]
ranks = np.argsort(np.argsort(Xg, axis=0, kind="stable"),
                   axis=0).astype(np.float64) + 1.0
m = 10000.0 * ranks / ng

scores = {}
for sname, gset in signatures.items():
    overlap = [g for g in gset if g in gsyms]
    sp = np.array([np.searchsorted(gsyms, g) for g in overlap])
    tag = np.zeros(ng, dtype=bool); tag[sp] = True
    Nh = len(sp); Nm = ng - Nh
    es = np.zeros(Xg.shape[1])
    order = np.argsort(-m, axis=0, kind="stable")
    for s in range(Xg.shape[1]):
        o = order[:, s]
        tg = tag[o]
        corr = np.abs(m[o, s]) ** 0.25
        Pn = np.where(tg, corr / corr[tg].sum(), 0.0)
        RES = np.cumsum(Pn) - np.cumsum((~tg) / Nm)
        es[s] = RES.sum()
    scores[sname] = es
    print(f"  {sname}: overlap {len(overlap)}/{len(gset)}")

stromal = scores["StromalSignature"]
estimate = stromal + scores["ImmuneSignature"]

# panel features (same pipeline as the nested analyses)
keep = counts.mean(axis=1) > 10
lk = lcpm[keep]
gk = np.array([g.split(".")[0] for g in genes[keep]])
sp = {}
for i, g in enumerate(gk):
    sp.setdefault(smap.get(g, ""), i)
Xp = lk[[sp[g] for g in CORE10], :].T
Zp = pscore(np.ascontiguousarray(Xp))

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
feats = {"stromal_only": stromal.reshape(-1, 1),
         "estimate_only": estimate.reshape(-1, 1),
         "panel": Zp,
         "stromal_plus_panel": np.column_stack([stromal, Zp])}
accs = {}
for name, F in feats.items():
    a = []
    for tr, te in skf.split(F, y):
        lr = LogisticRegression(max_iter=2000).fit(F[tr], y[tr])
        a.append(lr.score(F[te], y[te]))
    accs[name] = float(np.mean(a))
    print(f"  CV acc [{name}]: {np.mean(a):.1%}")

# additive test: panel vs stromal+panel, paired 5-fold
diffs = []
for tr, te in skf.split(Zp, y):
    lr1 = LogisticRegression(max_iter=2000).fit(Zp[tr], y[tr])
    F2 = np.column_stack([stromal, Zp])
    lr2 = LogisticRegression(max_iter=2000).fit(F2[tr], y[tr])
    diffs.append(lr1.score(Zp[te], y[te]) - lr2.score(F2[te], y[te]))

res = {"cv_acc": accs,
       "additive_delta_panel_minus_combo": float(np.mean(diffs)),
       "note": "ESTIMATE ported from the reference R implementation "
               "(stable-rank approximation)"}
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results",
                                 "06_estimate_baseline.json"), "w"),
          indent=2)
print(f"aditivo (painel - estromo+painel): {np.mean(diffs):+.1%}")
print(f"salvo 06_estimate_baseline.json ({time.time()-t0:.0f}s)")
