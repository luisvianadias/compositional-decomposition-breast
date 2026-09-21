# -*- coding: utf-8 -*-
"""Section 2.4-2.5 — second-model analyses with a fully public model.

The second model is the adipose-axis score rule: A = mean over the 15-gene
adipose-associated panel of per-gene z-scores anchored on TCGA, thresholded
at tau* chosen on TCGA out-of-fold scores (balanced-accuracy maximizing).
No proprietary components. Outputs: results/09_clinical_models.json
"""
import os, sys, json, math, time, requests
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, binomtest, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DATA, load_tcga, pscore, fetch_metabric,
                    fetch_clinical, CORE10, PANEL_B)

t0 = time.time()
CBIO = "https://www.cbioportal.org/api"
X, y, pos, gk = load_tcga()

# ---- model 1: LR on CORE10 (per-sample z; identical to 08 protocol) ----
Ztc = pscore(X[:, [pos[g] for g in CORE10]])
lr = LogisticRegression(max_iter=2000).fit(Ztc, y)

# ---- model 2: adipose-axis score A (PANEL_B, TCGA-anchored z) ----
colB = [pos[g] for g in PANEL_B if g in pos]
XB = X[:, colB]
muB, sdB = XB.mean(0), XB.std(0) + 1e-6
A_full = ((XB - muB) / sdB).mean(1)

def tau_star_ba(scores, yy):
    taus = np.unique(np.quantile(scores, np.linspace(0.01, 0.99, 99)))
    best, btau = -1, 0.5
    for t in taus:
        p = (scores < t).astype(int)          # low A => tumor
        ba = (p[yy == 1].mean() + 1 - p[yy == 0].mean()) / 2
        if ba > best:
            best, btau = ba, t
    return float(btau)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
oof = np.zeros(len(y))
for tr, te in skf.split(XB, y):
    m_, s_ = XB[tr].mean(0), XB[tr].std(0) + 1e-6
    oof[te] = ((XB[te] - m_) / s_).mean(1)
t_star = tau_star_ba(oof, y)
print(f"tau* (TCGA OOF, low-A = tumor): {t_star:.3f}", flush=True)

# ---- METABRIC: panel for both models ----
mm = fetch_metabric(sorted(set(CORE10) | set(PANEL_B)))
idx_m = mm[sorted(set(CORE10) & set(mm.columns))].dropna().index
Zm = pscore(mm.loc[idx_m, CORE10].to_numpy(dtype=np.float32))
pr_lr = lr.predict_proba(Zm)[:, 1]
call_lr = pr_lr > 0.5
# A-score on METABRIC: per-gene z standardized within the platform
# (same construction as the TCGA OOF scores; platform offset removed)
genesB = [g for g in PANEL_B if g in mm.columns]
Bm = mm.loc[idx_m, genesB].to_numpy(dtype=np.float32)
A_met = ((Bm - Bm.mean(0)) / (Bm.std(0) + 1e-6)).mean(1)
call_ax = A_met < t_star
calls = pd.DataFrame({"LR": call_lr, "AXIS": call_ax}, index=idx_m)
print(f"METABRIC captures: LR {call_lr.mean():.1%} | "
      f"AXIS {call_ax.mean():.1%}", flush=True)

# ---- clinical stratifications: capture by attribute, BH per model ----
clin = fetch_clinical(["ER_STATUS", "PR_STATUS", "GRADE",
                       "CANCER_TYPE_DETAILED"])
j = clin.join(calls, how="inner")
j["hist"] = j["CANCER_TYPE_DETAILED"].map(
    lambda x: "mixed" if ("Ductal" in str(x) and "Lobular" in str(x))
    else "ductal" if "Ductal" in str(x)
    else "lobular" if "Lobular" in str(x) else "other")

def bh_family(pv):
    pv = np.asarray(pv, dtype=float)
    m_ = len(pv)
    o = np.argsort(pv)
    b = np.empty(m_)
    b[o] = np.minimum.accumulate(
        (pv[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
    return b

def model_pvalues(col):
    pv = {}
    g = j.dropna(subset=["GRADE"]).copy()
    g["GRADE"] = g["GRADE"].astype(int)
    grp = g.groupby("GRADE")[col].agg(["sum", "count"])
    n_ = grp["count"].values.astype(float)
    cpt = grp["sum"].values.astype(float)
    sc = grp.index.values.astype(float)
    N, R = n_.sum(), cpt.sum()
    pbar = R / N
    sbar = (n_ * sc).sum() / N
    z_ca = ((sc * cpt).sum() - pbar * (n_ * sc).sum()) / \
           math.sqrt(pbar * (1 - pbar) * ((n_ * (sc - sbar) ** 2).sum()))
    pv["GRADE"] = 2 * (1 - 0.5 * (1 + math.erf(abs(z_ca) / math.sqrt(2))))
    for attr in ("ER_STATUS", "PR_STATUS"):
        d = j.dropna(subset=[attr])
        pv[attr] = chi2_contingency(pd.crosstab(d[attr], d[col]))[1]
    d2 = j[j["hist"] != "other"]
    pv["histology"] = chi2_contingency(pd.crosstab(d2["hist"], d2[col]))[1]
    b = bh_family(list(pv.values()))
    return dict(zip(pv, b))

bh_models = {c: model_pvalues(c) for c in ("LR", "AXIS")}
print("histology capture by model:")
print(j.groupby("hist")[["LR", "AXIS"]].agg(["mean", "count"]).round(3))

# ---- PAM50 (genefu ssp2006 centroids, Spearman rule; 08 machinery) ----
import rdata
ssp = rdata.read_rda(os.path.join(DATA, "ssp2006_genefu.rda"))["ssp2006"]
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
unique_ent = sorted({e for e in entrez if e is not None})
hugo_of = {}
for i in range(0, len(unique_ent), 100):
    chunk = unique_ent[i:i + 100]
    resp = requests.post(f"{CBIO}/genes/fetch", json=chunk,
                         params={"geneIdType": "ENTREZ_GENE_ID",
                                 "projection": "SUMMARY"},
                         timeout=120).json()
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
cgenes = [g for g in cent_genes if g in gpos_p]
sub = mpam[cgenes].to_numpy(dtype=np.float64).T
Cs = CENT[[cent_genes.index(g) for g in cgenes], :]
ranks = np.argsort(np.argsort(sub, axis=0), axis=0).astype(np.float64)
cranks = np.argsort(np.argsort(Cs, axis=0), axis=0).astype(np.float64)
subs, margins = [], []
for s_ in range(sub.shape[1]):
    cors = np.array([spearmanr(ranks[:, s_], cranks[:, c]).statistic
                     for c in range(Cs.shape[1])])
    subs.append(classes[int(np.argmax(cors))])
    o = np.sort(cors)[::-1]
    margins.append(o[0] - o[1])
pam = pd.DataFrame({"PAM50": subs, "margin": margins},
                   index=mpam.index).reindex(idx_m)

# ---- confidence stratum (margin >= p75): captures, chi2, retention ----
q75 = pam["margin"].quantile(0.75)
strat = pam["margin"] >= q75
df = pd.DataFrame({"PAM50": pam["PAM50"], "LR": calls["LR"],
                   "AXIS": calls["AXIS"], "strat": strat})
cap_full = {m_: df.groupby("PAM50")[m_].agg(["mean", "count"])
            for m_ in ("LR", "AXIS")}
chi_full = {m_: chi2_contingency(pd.crosstab(df["PAM50"], df[m_]))[1]
            for m_ in ("LR", "AXIS")}
chi_full4 = {m_: chi2_contingency(
    pd.crosstab(df[df["PAM50"] != "Normal"]["PAM50"],
                df[df["PAM50"] != "Normal"][m_]))[1] for m_ in ("LR", "AXIS")}
st = df[df["strat"]]
cap_strat = {m_: st.groupby("PAM50")[m_].agg(["mean", "count"])
             for m_ in ("LR", "AXIS")}
p_strat = {m_: chi2_contingency(pd.crosstab(st["PAM50"], st[m_]))[1]
           for m_ in ("LR", "AXIS")}
retention = df.groupby("PAM50")["strat"].agg(["mean", "count"])

# ---- ER validation of PAM50 calls ----
er = fetch_clinical(["ER_STATUS"])
ve = er.join(pam["PAM50"], how="inner").dropna()
er_basal = (ve[ve["PAM50"] == "Basal"]["ER_STATUS"] == "Positive").mean()
er_luma = (ve[ve["PAM50"] == "LumA"]["ER_STATUS"] == "Positive").mean()

# ---- split-half: threshold stability + model concordance (held-out) ----
rng = np.random.default_rng(42)
perm = rng.permutation(len(idx_m))
A_, B_ = perm[:len(perm) // 2], perm[len(perm) // 2:]
res_split = {"LR@0.5": {"half_A": float(call_lr[A_].mean()),
                        "half_B": float(call_lr[B_].mean())},
             "AXIS@t*": {"half_A": float(call_ax[A_].mean()),
                         "half_B": float(call_ax[B_].mean())}}
# equalize the AXIS threshold on half A to the LR capture rate, then
# measure AXIS-vs-LR concordance on the untouched half B
rate_lr_A = float(call_lr[A_].mean())
t_eq = float(np.quantile(A_met[A_], rate_lr_A))
call_ax_eq = A_met < t_eq
agree_B = float((call_ax_eq[B_] == call_lr[B_]).mean())
print(f"split-half: LR halfB {res_split['LR@0.5']['half_B']:.1%} | "
      f"AXIS halfB {res_split['AXIS@t*']['half_B']:.1%} | "
      f"concordancia (held-out B, limiar igualado na A): {agree_B:.1%}")

out = {
    "tau_star": t_star,
    "overall_capture": {"LR": float(call_lr.mean()),
                        "AXIS": float(call_ax.mean())},
    "bh_by_model": bh_models,
    "hist_capture": {m_: j.groupby("hist")[m_].agg(["mean", "count"])
                     .apply(lambda r: {"mean": float(r["mean"]),
                                       "n": int(r["count"])}, axis=1)
                     .to_dict() for m_ in ("LR", "AXIS")},
    "stratum": {"threshold": float(q75), "n": int(strat.sum()),
                "chi2_lr": p_strat["LR"], "chi2_axis": p_strat["AXIS"],
                "capture": {m_: {str(k): {"mean": float(v["mean"]),
                                          "n": int(v["count"])}
                                 for k, v in cap_strat[m_].iterrows()}
                            for m_ in ("LR", "AXIS")}},
    "pam50_full_capture": {m_: {str(k): {"mean": float(v["mean"]),
                                         "n": int(v["count"])}
                                for k, v in cap_full[m_].iterrows()}
                           for m_ in ("LR", "AXIS")},
    "pam50_full_chi2": {"lr": chi_full["LR"], "axis": chi_full["AXIS"],
                        "lr_excl_normal": chi_full4["LR"],
                        "axis_excl_normal": chi_full4["AXIS"]},
    "retention": {str(k): {"mean": float(v["mean"]), "n": int(v["count"])}
                  for k, v in retention.iterrows()},
    "pam50_er_validation": {"er_positive_basal": float(er_basal),
                            "er_positive_lumA": float(er_luma)},
    "split_half": res_split,
    "concordance_heldout": agree_B,
    "axis_threshold_matched": {"half_A_lr_rate": rate_lr_A,
                               "t_eq": t_eq},
    "runtime_s": round(time.time() - t0, 1),
}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "results",
                                 "09_clinical_models.json"), "w"),
          indent=2, default=float)
print("salvo 09_clinical_models.json")
