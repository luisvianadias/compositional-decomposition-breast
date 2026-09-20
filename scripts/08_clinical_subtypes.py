# -*- coding: utf-8 -*-
"""Section 2.5 — clinical stratifications (BH family) and PAM50 intrinsic
subtype capture. Logistic regression, CORE10 panel, per-sample z.
Output: results/08_clinical_subtypes.json
"""
import os, sys, json, math, time, requests
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, binomtest, norm
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DATA, load_tcga, pscore, fetch_metabric,
                    fetch_clinical, CORE10)

# TCGA training (CORE10, per-sample z)
X, y, pos, gk = load_tcga()
Ztc = pscore(X[:, [pos[g] for g in CORE10]])
lr = LogisticRegression(max_iter=2000).fit(Ztc, y)

# METABRIC capture
mm = fetch_metabric(CORE10)
idx_m = mm[CORE10].dropna().index
Zm = pscore(mm.loc[idx_m, CORE10].to_numpy(dtype=np.float32))
pr = lr.predict_proba(Zm)[:, 1]
calls = pd.DataFrame({"LR": pr > 0.5}, index=idx_m)

# clinical family (4 tests, BH)
clin = fetch_clinical(["ER_STATUS", "PR_STATUS", "GRADE",
                       "CANCER_TYPE_DETAILED"])
j = clin.join(calls, how="inner")

pvals = {}
g = j.dropna(subset=["GRADE"]).copy()
g["GRADE"] = g["GRADE"].astype(int)
grp = g.groupby("GRADE")["LR"].agg(["sum", "count"])
n_ = grp["count"].values.astype(float)
cpt = grp["sum"].values.astype(float)
sc = grp.index.values.astype(float)
N, R = n_.sum(), cpt.sum()
pbar = R / N
sbar = (n_ * sc).sum() / N
z_ca = ((sc * cpt).sum() - pbar * (n_ * sc).sum()) / \
       math.sqrt(pbar * (1 - pbar) * ((n_ * (sc - sbar) ** 2).sum()))
pvals["GRADE(Cochran-Armitage)"] = \
    2 * (1 - 0.5 * (1 + math.erf(abs(z_ca) / math.sqrt(2))))
for attr in ("ER_STATUS", "PR_STATUS"):
    pvals[attr] = chi2_contingency(
        pd.crosstab(j.dropna(subset=[attr])[attr],
                    j.dropna(subset=[attr])["LR"]))[1]
h = j["CANCER_TYPE_DETAILED"].map(
    lambda x: "mixed" if ("Ductal" in str(x) and "Lobular" in str(x))
    else "ductal" if "Ductal" in str(x)
    else "lobular" if "Lobular" in str(x) else "other")
j["hist"] = h
j2 = j[j["hist"] != "other"]
pvals["histology"] = chi2_contingency(
    pd.crosstab(j2["hist"], j2["LR"]))[1]

names = list(pvals)
pv = np.array([pvals[n] for n in names])
m_ = len(pv)
o = np.argsort(pv)
bh = np.empty(m_)
bh[o] = np.minimum.accumulate((pv[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
print("clinical family (BH):")
for n_, p0, b in zip(names, pv, bh):
    print(f"  {n_:<10} p={p0:.2e} BH={b:.2e}")

# PAM50 (genefu ssp2006 centroids; rdata parse)
import rdata
ssp_path = os.path.join(DATA, "ssp2006_genefu.rda")
ssp = rdata.read_rda(ssp_path)["ssp2006"]
cent = ssp["centroids"]
classes = [str(c) for c in cent.coords[cent.dims[1]].values]
cmap = ssp["centroids.map"]
entrez = []
for p_ in cmap["EntrezGene.ID"]:
    try:
        if p_ is None or (isinstance(p_, float) and np.isnan(p_)) or \
           str(type(p_)).find("NA") >= 0:
            entrez.append(None)
        else:
            entrez.append(int(p_))
    except (TypeError, ValueError):
        entrez.append(None)
CBIO = "https://www.cbioportal.org/api"
unique_ent = sorted({e for e in entrez if e is not None})
hugo_of = {}
for i in range(0, len(unique_ent), 100):
    chunk = unique_ent[i:i + 100]
    resp = requests.post(f"{CBIO}/genes/fetch", json=chunk,
                         params={"geneIdType": "ENTREZ_GENE_ID",
                                 "projection": "SUMMARY"}, timeout=120).json()
    for gg in resp:
        hugo_of[gg["entrezGeneId"]] = gg["hugoGeneSymbol"]
    time.sleep(0.2)
gene_cent = {}
for row, e in zip(range(cent.shape[0]), entrez):
    if e is None or e not in hugo_of:
        continue
    gene_cent.setdefault(hugo_of[e], []).append(
        np.asarray(cent.values[row, :], dtype=np.float64))
cent_genes = sorted(gene_cent)
CENT = np.stack([np.mean(gene_cent[g], axis=0) for g in cent_genes])

from scipy.stats import spearmanr
sym2ent = {g["hugoGeneSymbol"]: g["entrezGeneId"]
           for g in requests.post(f"{CBIO}/genes/fetch", json=cent_genes,
                                  params={"geneIdType": "HUGO_GENE_SYMBOL",
                                          "projection": "SUMMARY"},
                                  timeout=120).json()}
e2s = {v: k for k, v in sym2ent.items()}
md = requests.post(
    f"{CBIO}/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
    json={"entrezGeneIds": list(e2s.keys()),
          "sampleListId": "brca_metabric_all"},
    params={"projection": "SUMMARY"}, timeout=180).json()
rows = [{"s": d["sampleId"], "g": e2s[d["entrezGeneId"]], "v": d["value"]}
        for d in md
        if d.get("value") is not None and d.get("value") == d.get("value")]
mpam = pd.DataFrame(rows).pivot_table(index="s", columns="g", values="v")
gpos_p = {g: i for i, g in enumerate(mpam.columns)}
M = mpam.to_numpy(dtype=np.float64).T          # genes x samples
idx = [gpos_p[g] for g in cent_genes if g in gpos_p]
cgenes = [g for g in cent_genes if g in gpos_p]
sub = M[idx, :]
Cs = CENT[[cent_genes.index(g) for g in cgenes], :]
ranks = np.argsort(np.argsort(sub, axis=0), axis=0).astype(np.float64)
cranks = np.argsort(np.argsort(Cs, axis=0), axis=0).astype(np.float64)
subs = []
for s in range(sub.shape[1]):
    cors = [spearmanr(ranks[:, s], cranks[:, c]).statistic
            for c in range(Cs.shape[1])]
    subs.append(classes[int(np.argmax(cors))])
subs = np.array(subs)

sam = pd.Series(subs, index=mpam.index).reindex(idx_m)
Z10 = pscore(mm.loc[idx_m, CORE10].to_numpy(dtype=np.float32))
pr10 = lr.predict_proba(Z10)[:, 1]
df = pd.DataFrame({"PAM50": sam.values, "capt": pr10 > 0.5})
cap_sub = df.groupby("PAM50")["capt"].agg(["mean", "count"])
print("\nPAM50 capture (LR):")
print(cap_sub.round(3).to_string())

# PAM50 x capture chi-square, and without Normal-like
tab5 = pd.crosstab(df["PAM50"], df["capt"])
p5 = chi2_contingency(tab5)[1]
df4 = df[df["PAM50"] != "Normal"]
p4 = chi2_contingency(pd.crosstab(df4["PAM50"], df4["capt"]))[1]

# joint BH over the five-family tests + PAM50
import time as _t
all_p = dict(pvals)
all_p["PAM50(with Normal-like)"] = p5
nm_ = list(all_p)
pv = np.array([all_p[n] for n in nm_])
m_ = len(pv)
o = np.argsort(pv)
bhj = np.empty(m_)
bhj[o] = np.minimum.accumulate((pv[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
print("\njoint BH (5 clinical + PAM50):")
for n_, b in zip(nm_, bhj):
    print(f"  {n_:<12} BH={b:.2e}")

res = {
    "clinical_bh": {n_: float(b) for n_, b in zip(names, bh[:m_ - 1])},
    "pam50_chi2_with_normal": float(p5),
    "pam50_chi2_without_normal": float(p4),
    "joint_bh": {n_: float(b) for n_, b in zip(nm_, bhj)},
    "pam50_capture": {str(k): {"mean": float(v["mean"]), "n": int(v["count"])}
                      for k, v in cap_sub.iterrows()},
}
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results",
                                 "08_clinical_subtypes.json"), "w"),
          indent=2, default=float)
print("salvo 08_clinical_subtypes.json")
