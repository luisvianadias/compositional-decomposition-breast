# Scripts

Transcrição de 18 arquivos `.py` de `scripts/`.

---

## 01_download_data.py

```python
# -*- coding: utf-8 -*-
"""Downloads all public data (skips files already present).
- GDC API query -> manifest -> TCGA-BRCA STAR counts (open access),
  canonical cohort: all normals + seed-42 draw of 520 tumors (631 files)
- GTEx v10 gct.gz (breast, adipose subcutaneous) from public GCS bucket
- genefu ssp2006 centroids (GitHub, bhklab/genefu)
- ESTIMATE reference package (SourceForge)
"""
import os, sys, json, csv, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor

sys.stdout.reconfigure(encoding="utf-8")
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
os.makedirs(DATA, exist_ok=True)

def exists(p, min_size=1000):
    return os.path.exists(p) and os.path.getsize(p) > min_size

def fetch(url, dest):
    if exists(dest, 500):
        print(f"  [skip] {os.path.basename(dest)}")
        return
    print(f"  [get ] {os.path.basename(dest)}")
    urllib.request.urlretrieve(url, dest)

# ---- GTEx v10 ----
GCS = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/" \
      "counts-by-tissue/"
fetch(GCS + "gene_reads_v10_breast_mammary_tissue.gct.gz",
      os.path.join(DATA, "gtex_breast_v10_reads.gct.gz"))
fetch(GCS + "gene_reads_v10_adipose_subcutaneous.gct.gz",
      os.path.join(DATA, "gtex_adipose_subcut_v10_reads.gct.gz"))

# ---- genefu ssp2006 centroids ----
fetch("https://raw.githubusercontent.com/bhklab/genefu/master/data/"
      "ssp2006.rda", os.path.join(DATA, "ssp2006_genefu.rda"))

# ---- ESTIMATE reference package ----
est = os.path.join(DATA, "estimate_r.tar.gz")
if not exists(est, 1_000_000):
    print("  [get ] estimate_r.tar.gz (SourceForge)")
    req = urllib.request.Request(
        "https://sourceforge.net/projects/estimateproject/files/latest/"
        "download", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=300) as r, \
            open(est, "wb") as f:
        f.write(r.read())
import tarfile
pkg = os.path.join(DATA, "estimate_pkg", "estimate", "inst", "extdata")
if not os.path.exists(os.path.join(pkg, "SI_geneset.gmt")):
    print("  [extract] estimate_pkg")
    with tarfile.open(est) as t:
        t.extractall(os.path.join(DATA, "estimate_pkg"))

# ---- GDC: manifest query ----
man = os.path.join(DATA, "gdc_manifest.tsv")
if not exists(man, 100):
    print("  [get ] GDC manifest (API query)")
    filters = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": ["TCGA-BRCA"]}},
        {"op": "in", "content": {"field": "files.data_type",
                                 "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "files.analysis.workflow_type",
                                 "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "files.access",
                                 "value": ["open"]}}]}
    rows = []
    for offset in (0, 1000):
        params = urllib.parse.urlencode({
            "filters": json.dumps(filters),
            "fields": "file_id,cases.samples.submitter_id,"
                      "cases.samples.sample_type",
            "size": "1000", "from": str(offset)})
        req = urllib.request.Request(
            "https://api.gdc.cancer.gov/files?" + params)
        r = json.load(urllib.request.urlopen(req, timeout=120))
        for h in r["data"]["hits"]:
            if not h.get("cases"):
                continue
            s = h["cases"][0].get("samples", [{}])[0]
            rows.append((h["id"], s.get("submitter_id", ""),
                         s.get("sample_type", "")))
    with open(man, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["file_id", "sample_barcode", "sample_type"])
        w.writerows(rows)
    print(f"  manifesto: {len(rows)} arquivos")

# ---- GDC: download counts (cohort: all normals + seeded tumor draw) ----
# Canonical cohort (matches the manuscript: 518 tumors + 113 normals = 631):
# all "Solid Tissue Normal" files + a random 520-file draw of "Primary
# Tumor" (random.Random(42).shuffle), keeping the FIRST file per barcode.
# Two drawn barcodes have 2 GDC file entries (reprocessed workflow); the
# first in draw order is the canonical file (byte-verified, see
# data/PROVENANCE.md).
import random

gdc_dir = os.path.join(DATA, "gdc_brca")
os.makedirs(gdc_dir, exist_ok=True)
rows = list(csv.DictReader(open(man), delimiter="\t"))
normals = [r for r in rows if r["sample_type"] == "Solid Tissue Normal"]
tumors = [r for r in rows if r["sample_type"] == "Primary Tumor"]
random.Random(42).shuffle(tumors)
cohort, seen = [], set()
for r in normals + tumors[:520]:
    if r["sample_barcode"] not in seen:
        seen.add(r["sample_barcode"])
        cohort.append(r)
print(f"  coorte: {len(cohort)} amostras "
      f"({sum(1 for r in cohort if r['sample_type'] != 'Solid Tissue Normal')}"
      f" tumores + "
      f"{sum(1 for r in cohort if r['sample_type'] == 'Solid Tissue Normal')}"
      f" normais)")
todo = [r for r in cohort
        if not exists(os.path.join(gdc_dir, r["sample_barcode"] + ".tsv"),
                      1_000_000)]
print(f"  GDC counts: {len(cohort) - len(todo)} presentes, "
      f"{len(todo)} a baixar (~10 min na primeira vez)")
t0 = time.time()
done = 0
failed = []

def dl(row):
    dest = os.path.join(gdc_dir, row["sample_barcode"] + ".tsv")
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(
                f"https://api.gdc.cancer.gov/data/{row['file_id']}", dest)
            return True
        except Exception:
            time.sleep(2 * (attempt + 1))
    failed.append(row["sample_barcode"])
    return False

with ThreadPoolExecutor(max_workers=12) as ex:
    for ok in ex.map(dl, todo):
        done += 1
        if done % 50 == 0:
            print(f"    {done}/{len(todo)} ({time.time()-t0:.0f}s)",
                  flush=True)
if failed:
    raise SystemExit(f"FALHARAM {len(failed)} downloads "
                     f"(reexecute para retomar): {sorted(failed)[:10]} ...")
print("  GDC download concluido")
```

---

## 01b_download_estimate.py

```python
# -*- coding: utf-8 -*-
"""Downloads the ESTIMATE reference package (gene sets + common genes).
Needed by 06_estimate_baseline.py."""
import os, sys, tarfile, urllib.request

sys.stdout.reconfigure(encoding="utf-8")
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
DEST = os.path.join(DATA, "estimate_r.tar.gz")
PKG = os.path.join(DATA, "estimate_pkg")
os.makedirs(DATA, exist_ok=True)

ext_data = os.path.join(PKG, "estimate", "inst", "extdata")
if os.path.exists(os.path.join(ext_data, "SI_geneset.gmt")):
    print("[skip] estimate_pkg already extracted")
    sys.exit(0)

if not os.path.exists(DEST) or os.path.getsize(DEST) < 1_000_000:
    print("[get ] estimate_r.tar.gz (~3.7 MB, SourceForge)")
    req = urllib.request.Request(
        "https://sourceforge.net/projects/estimateproject/files/latest/"
        "download", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=300) as r, open(DEST, "wb") as f:
        f.write(r.read())

print("[extract] estimate_pkg...")
with tarfile.open(DEST) as t:
    t.extractall(os.path.join(DATA, "estimate_pkg"))

ext = os.path.join(PKG, "estimate", "inst", "extdata")
for fn in ("SI_geneset.gmt", "common_genes.txt"):
    p = os.path.join(ext, fn)
    assert os.path.exists(p), f"{fn} not found after extraction"
print(f"  SI_geneset.gmt: {os.path.getsize(os.path.join(ext, 'SI_geneset.gmt'))} bytes")
print(f"  common_genes.txt: {os.path.getsize(os.path.join(ext, 'common_genes.txt'))} bytes")
print("estimate_pkg ready")
```

---

## 02_build_matrices.py

```python
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
```

---

## 03_nested_compression.py

```python
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
from common import load_tcga, PANEL_B, PANEL_C
from pipeline_utils import pick_nearest_non_panel

t0 = time.time()
DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "..", "results")
os.makedirs(OUT, exist_ok=True)

X, y, pos, _ = load_tcga()
N_GENES = X.shape[1]

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
aucs_null_acc = np.zeros(N_NULL)

# null: 500 uniformly random 10-gene panels, drawn ONCE and paired across
# all splits (each panel evaluated on 15 test folds). The manuscript
# quotes the median of the pooled per-draw accuracies.
rng_null = np.random.default_rng(42)
null_panels = np.array([rng_null.choice(N_GENES, size=NULL_K, replace=False)
                        for _ in range(N_NULL)])

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

    # null: fixed panel set, per-split accuracies pooled across splits
    draws = np.zeros(N_NULL)
    aucs_null = np.zeros(N_NULL)
    for b in range(N_NULL):
        cols = null_panels[b]
        m = LogisticRegression(max_iter=2000).fit(Ztr[:, cols], ytr)
        draws[b] = m.score(Zte[:, cols], yte)
        aucs_null[b] = roc_auc_score(yte, m.predict_proba(Zte[:, cols])[:, 1])
    arms["null"]["acc"][i] = float(draws.mean())
    arms["null"]["auc"][i] = float(aucs_null.mean())
    null_draws_acc += draws / len(splits)
    aucs_null_acc += aucs_null / len(splits)

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
res["null_auc_per_draw_mean_iqr"] = [float(np.quantile(aucs_null_acc, q))
                                     for q in (0.25, 0.5, 0.75)]
json.dump(res, open(os.path.join(OUT, "nested_compression.json"), "w"),
          indent=2)
print(f"\nsalvo nested_compression.json ({time.time()-t0:.0f}s)")
```

---

## 04_adipose_ablation.py

```python
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
```

---

## 05_cross_cohort_axis.py

```python
# -*- coding: utf-8 -*-
"""Section 2.4 — adipose-axis score across populations (pooled z, pure
16-gene adipocyte panel; positive control: GTEx subcutaneous adipose).
Output: results/05_cross_cohort_axis.json
"""
import os, sys, json
import numpy as np
from common import DATA, load_tcga, load_gct, PURE16, strip_version

X, y, pos, gk = load_tcga()
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
```

---

## 06_estimate_baseline.py

```python
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
gdc_dir = os.path.join(data_dir, "gdc_brca")
gdc = os.path.join(gdc_dir, sorted(
    f for f in os.listdir(gdc_dir) if f.endswith(".tsv"))[0])
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
order_g = np.argsort(sym_all[mask], kind="stable")
Xg = lcpm[mask][order_g]
gsyms = np.array(sym_all[mask][order_g])   # sorted; rows aligned to gsyms
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
         "panel": Zp}
accs = {}
for name, F in feats.items():
    a = []
    for tr, te in skf.split(F, y):
        lr = LogisticRegression(max_iter=2000).fit(F[tr], y[tr])
        a.append(lr.score(F[te], y[te]))
    accs[name] = float(np.mean(a))
    print(f"  CV acc [{name}]: {np.mean(a):.1%}")

# combo: stromal z-scored with train-fold stats so the L2 penalty sees both
# features on comparable scales (panel is already per-gene z)
a = []
for tr, te in skf.split(Zp, y):
    mu_s, sd_s = stromal[tr].mean(), stromal[tr].std() + 1e-6
    F2 = np.column_stack([(stromal - mu_s) / sd_s, Zp])
    lr2 = LogisticRegression(max_iter=2000).fit(F2[tr], y[tr])
    a.append(lr2.score(F2[te], y[te]))
accs["stromal_plus_panel"] = float(np.mean(a))
print(f"  CV acc [stromal_plus_panel]: {np.mean(a):.1%}")

# additive test: panel vs stromal+panel, paired 5-fold (stromal z-scored
# with train-fold stats)
diffs = []
for tr, te in skf.split(Zp, y):
    lr1 = LogisticRegression(max_iter=2000).fit(Zp[tr], y[tr])
    mu_s, sd_s = stromal[tr].mean(), stromal[tr].std() + 1e-6
    F2 = np.column_stack([(stromal - mu_s) / sd_s, Zp])
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
```

---

## 07_transfer.py

```python
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
syms = [g for g in PANEL_C if g in pos]
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
```

---

## 08_clinical_subtypes.py

```python
# -*- coding: utf-8 -*-
"""Section 2.5 — clinical stratifications (BH family) and PAM50 intrinsic
subtype capture. Logistic regression, CORE10 panel, per-sample z.
Output: results/08_clinical_subtypes.json
"""
import os, sys, json, math, time, requests
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
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
    "clinical_bh": {n_: float(b) for n_, b in zip(names, bh)},
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
```

---

## 09_clinical_models.py

```python
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
```

---

## common.py

```python
# -*- coding: utf-8 -*-
"""Shared loaders and helpers — scikit-learn/SciPy/pandas only."""
import os, gzip, json
import numpy as np
import pandas as pd
import requests

DATA = os.environ.get("BCD_DATA_DIR", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "data"))

PANEL_A = ["ADIPOQ", "FABP4", "LEP", "PLIN1", "PPARG", "CD36", "LPL",
           "CIDEA", "LIPE", "DGAT2", "CFD", "ENPP2"]
PANEL_B = PANEL_A + ["ADH1B", "GPD1", "CIDEC"]
PANEL_C = PANEL_B + ["COL10A1", "MMP11"]
CORE10 = ["ADH1B", "FABP4", "LEP", "COL10A1", "PLIN1",
          "MMP11", "PLIN4", "ADIPOQ", "GPD1", "CIDEC"]
PURE16 = ["ADIPOQ", "FABP4", "LEP", "PLIN1", "PLIN4", "PPARG", "CD36",
          "LPL", "CIDEA", "CIDEC", "ADH1B", "GPD1", "DGAT2", "ENPP2",
          "CFD", "LIPE"]


def pscore(mat):
    """Per-sample z-score across the panel (row-wise)."""
    return (mat - mat.mean(axis=1, keepdims=True)) / \
           (mat.std(axis=1, keepdims=True) + 1e-6)


def strip_version(ensg):
    return ensg.split(".")[0]


def build_symbol_index(symbols_full, keep_mask):
    """symbol -> first position in the FILTERED axis."""
    kept = [strip_version(s) for s, k in zip(symbols_full, keep_mask) if k]
    idx = {}
    for i, s in enumerate(kept):
        idx.setdefault(s, i)
    return idx


def extract_gene_matrix(mat_full_filtered, symbol_pos, wanted):
    """mat (samples x filtered genes) -> samples x wanted. Raises on
    missing symbols (axis-mismatch guard, see README)."""
    missing = [g for g in wanted if g not in symbol_pos]
    if missing:
        raise KeyError(f"symbols absent from filtered axis: {missing}")
    return np.ascontiguousarray(mat_full_filtered[:, [symbol_pos[g]
                                                     for g in wanted]])


def load_tcga(data_dir=DATA):
    d = np.load(os.path.join(data_dir, "tcga_brca_counts.npz"),
                allow_pickle=True)
    counts, genes = d["counts"], d["genes"].astype(str)
    samples = pd.read_csv(os.path.join(data_dir, "tcga_brca_samples.tsv"),
                          sep="\t")
    y = (samples["type"] == "tumor").to_numpy().astype(int)
    lib = counts.sum(axis=0, keepdims=True)
    lcpm = np.log2(counts / lib * 1e6 + 1.0)
    keep = counts.mean(axis=1) > 10
    X = np.ascontiguousarray(lcpm[keep].T, dtype=np.float32)
    gk = np.array([g.split(".")[0] for g in genes[keep]])
    sym_file = _first_gdc_tsv(data_dir)
    tsv = pd.read_csv(sym_file, sep="\t", comment="#",
                      usecols=["gene_id", "gene_name"]).dropna()
    tsv["gene_id"] = tsv["gene_id"].str.split(".").str[0]
    smap = dict(tsv.drop_duplicates("gene_id").set_index("gene_id")["gene_name"])
    pos = {}
    for i, g in enumerate(gk):
        pos.setdefault(smap.get(g, ""), i)
    return X, y, pos, gk


def _first_gdc_tsv(data_dir):
    gdc = os.path.join(data_dir, "gdc_brca")
    return os.path.join(gdc, sorted(
        f for f in os.listdir(gdc) if f.endswith(".tsv"))[0])


def load_gct(path):
    with gzip.open(path, "rt") as f:
        f.readline(); f.readline()
        header = f.readline().rstrip("\n").split("\t")
        donors = header[2:]
        cols = {s: [] for s in donors}
        names = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            names.append(p[1])                    # Description = symbol
            for j, s in enumerate(donors):
                cols[s].append(p[2 + j])
    mat = np.array([cols[s] for s in donors], dtype=np.float32).T
    gpos = {g: i for i, g in enumerate(names)}
    lg = mat.sum(axis=0, keepdims=True)
    return np.log2(mat / lg * 1e6 + 1.0), gpos, len(donors)


def fetch_metabric(panel):
    CBIO = "https://www.cbioportal.org/api"
    sym2ent = {g["hugoGeneSymbol"]: g["entrezGeneId"]
               for g in requests.post(f"{CBIO}/genes/fetch", json=panel,
                                      params={"geneIdType":
                                              "HUGO_GENE_SYMBOL",
                                              "projection": "SUMMARY"},
                                      timeout=60).json()}
    e2s = {v: k for k, v in sym2ent.items()}
    md = requests.post(
        f"{CBIO}/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
        json={"entrezGeneIds": list(e2s.keys()),
              "sampleListId": "brca_metabric_all"},
        params={"projection": "SUMMARY"}, timeout=180).json()
    rows = [{"s": d["sampleId"], "g": e2s[d["entrezGeneId"]], "v": d["value"]}
            for d in md
            if d.get("value") is not None and d.get("value") == d.get("value")]
    return pd.DataFrame(rows).pivot_table(index="s", columns="g", values="v")


def fetch_clinical(attributes):
    cd = pd.DataFrame(requests.get(
        "https://www.cbioportal.org/api/studies/brca_metabric/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "DETAILED",
                "size": 100000, "direction": "ASC"}, timeout=120).json())
    return cd[cd["clinicalAttributeId"].isin(attributes)].pivot_table(
        index="sampleId", columns="clinicalAttributeId", values="value",
        aggfunc="first")
```

---

## figure2_pam50.py

```python
# -*- coding: utf-8 -*-
"""Figure 2 — Capture by intrinsic subtype (METABRIC, n=1,980).
Single classifier (logistic regression). Clean, no model comparison."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "results")
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "figures")
os.makedirs(FIG, exist_ok=True)

C_ORANGE = "#E69F00"

C_BAR = "#E69F00"
C_TEXT = "#B36B00"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

data_path = os.path.join(RES, "08_clinical_subtypes.json")
if os.path.exists(data_path):
    data = json.load(open(data_path))
    cap = data["pam50_capture"]
else:
    cap = {
        "Basal": {"mean": 0.90, "n": 220},
        "Her2": {"mean": 0.96, "n": 122},
        "LumB": {"mean": 1.00, "n": 52},
        "LumA": {"mean": 0.92, "n": 1577},
        "Normal": {"mean": 0.67, "n": 9},
    }

order = ["Basal", "Her2", "LumB", "LumA", "Normal"]
labels, n_str, vals = [], [], []

for st in order:
    v = cap.get(st, cap.get(st.replace("Basal", "Basal-like"), None))
    if v is None:
        continue
    labels.append(st)
    vals.append(v["mean"] * 100)
    n_str.append(f"n={int(v['n']):,}")

fig, ax = plt.subplots(figsize=(5.5, 3.4))
x = np.arange(len(vals))

ax.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)

bars = ax.bar(x, vals, 0.48, color=C_ORANGE, alpha=0.9,
              edgecolor="none", zorder=3)

for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 2.0, f"{v:.0f}%",
            ha="center", va="bottom", fontsize=7.5, fontweight="bold",
            color="#B36B00")

ax.axhline(50, color="#D95F02", ls="--", lw=0.9, alpha=0.6, zorder=2,
           label="50% (chance)")

ax.set_xticks(x)
ax.set_xticklabels(
    [f"{l}\n({n})" for l, n in zip(labels, n_str)],
    fontsize=7.8, color="#222222")
ax.set_ylabel("Capture (%)", fontsize=8.5, labelpad=6)
ax.set_ylim(0, 112)
ax.set_yticks([0, 20, 40, 60, 80, 100])

ax.legend(loc="upper right", fontsize=7.5, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG, "figure2_pam50.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "figure2_pam50.png"), bbox_inches="tight",
            dpi=300)
plt.close(fig)
print("Figure 2 salva — barras LogReg apenas, sem comparação de modelos")
```

---

## generate_paper_figures.py

```python
# -*- coding: utf-8 -*-
"""Gera Figura 1 e Figura S1 com proporções e contenção de texto rigorosas.
(versão do usuário, lendo results/ público — nested_compression.json +
05_cross_cohort_axis.json)"""

import json
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from matplotlib.lines import Line2D

matplotlib.use("Agg")
sys.stdout.reconfigure(encoding="utf-8")

DIR = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(DIR, "..", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 8,
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    }
)

C_NULL = "#999999"
C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_GREEN = "#009E73"
C_YELLOW = "#E69F00"
C_GRID = "#F0F0F0"

nested = json.load(open(os.path.join(DIR, "..", "results",
                                     "nested_compression.json")))
axis = json.load(open(os.path.join(DIR, "..", "results",
                                   "05_cross_cohort_axis.json")))
ablation = json.load(open(os.path.join(DIR, "..", "results",
                                       "04_adipose_ablation.json")))

K = [1, 2, 3, 5, 10, 20]
curve = [nested[f"k{k}"]["acc"] * 100 for k in K]
full_acc = nested["full"]["acc"] * 100
panel_accs = {n: nested[n]["acc"] * 100
              for n in ("panelA", "panelB", "panelC", "randDE")}
resid_acc = ablation["after"]["acc"] * 100
resid_auc = ablation["after"]["auc"]
null_draws = np.array(nested["null_acc_per_draw"]) * 100
null_med = float(np.median(null_draws))
null_q25, null_q75 = [float(q) for q in
                      np.quantile(null_draws, [0.25, 0.75])]

ladder_arms = [
    ("Null\n(B=500)", null_draws),
    ("Panel A\n(12 adip.)",
     np.array(nested["panelA"]["per_split_acc"]) * 100),
    ("Panel B\n(15 ad.-ass.)",
     np.array(nested["panelB"]["per_split_acc"]) * 100),
    ("randDE\n(matched)",
     np.array(nested["randDE"]["per_split_acc"]) * 100),
    ("Panel C\n(+stroma)",
     np.array(nested["panelC"]["per_split_acc"]) * 100),
    ("Full\n(23k genes)",
     np.array(nested["full"]["per_split_acc"]) * 100),
]


def panel_tag(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", va="top", ha="right")


fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(
    2, 3, height_ratios=[1.0, 1.05], width_ratios=[0.85, 0.95, 1.2],
    hspace=0.48, wspace=0.38, left=0.08, right=0.98, top=0.94, bottom=0.08,
)

axA = fig.add_subplot(gs[0, 0])
panel_tag(axA, "A")
axA.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
axA.axhspan(null_q25, null_q75, color=C_NULL, alpha=0.25, lw=0, zorder=1)
axA.axhline(null_med, color=C_NULL, lw=0.9, ls="--", zorder=2)
axA.axhline(full_acc, color="#333333", lw=0.8, ls=":", zorder=2)
axA.plot(K, curve, "o-", color=C_BLUE, lw=1.4, ms=4.2, mec="white",
         mew=0.5, zorder=3)
axA.set_xscale("log")
axA.set_xticks(K)
axA.set_xticklabels([str(k) for k in K])
axA.set_xlabel("Genes in nested panel ($k$)")
axA.set_ylabel("Accuracy (%)")
axA.set_ylim(80, 101)
axA.text(0.05, 0.08,
         f"Full = {full_acc:.1f}%\nNull median = {null_med:.1f}%",
         transform=axA.transAxes, fontsize=6.8, color="#444444")

axB = fig.add_subplot(gs[0, 1:])
panel_tag(axB, "B")
axB.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
colors = [C_NULL, C_BLUE, C_BLUE, C_YELLOW, C_GREEN, "#333333"]
pos = np.arange(len(ladder_arms))
data = [a * 100 if a.max() <= 1.01 else a for _, a in ladder_arms]
bp = axB.boxplot(data, positions=pos, widths=0.52, patch_artist=True,
                 showfliers=False,
                 medianprops=dict(color="black", lw=1.1),
                 whiskerprops=dict(lw=0.7, color="#444444"),
                 capprops=dict(lw=0.7, color="#444444"), zorder=2)
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.85 if c != C_NULL else 0.45)
    patch.set_edgecolor("#222222")
    patch.set_linewidth(0.7)
axB.set_xticks(pos)
axB.set_xticklabels([n for n, _ in ladder_arms], fontsize=6.8)
axB.set_ylabel("Accuracy (%)")
axB.set_ylim(78, 100.5)

axC = fig.add_subplot(gs[1, 0])
panel_tag(axC, "C")
axC.grid(axis="x", color=C_GRID, lw=0.6, zorder=0)
y_pos = [1, 0]
for y, val, col in zip(y_pos, [full_acc, resid_acc], [C_BLUE, C_YELLOW]):
    axC.hlines(y, xmin=82, xmax=val, color=col, lw=2.0, zorder=2)
    axC.plot(val, y, "o", color=col, ms=7, mec="black", mew=0.7, zorder=3)
    axC.text(val + 0.6, y, f"{val:.1f}%", va="center", ha="left",
             fontweight="bold", fontsize=7.2)
axC.set_yticks(y_pos)
axC.set_yticklabels(["Full model", "Residual model\n(non-adipose)"],
                    fontsize=7.2)
axC.set_xlim(82, 103)
axC.set_xlabel("Accuracy (%)")
axC.text(0.06, 0.14,
         f"AUC: {ablation['before']['auc']:.3f} -> {resid_auc:.3f}",
         transform=axC.transAxes, fontsize=6.8, color="#444444")

axD = fig.add_subplot(gs[1, 1])
panel_tag(axD, "D")
axD.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)
pops = ["TCGA\ntumor", "TCGA\nnormal", "GTEx\nbreast", "GTEx\nadipose"]
valsD = [axis["tcga_tumor"], axis["tcga_adjacent_normal"],
         axis["gtex_breast"], axis["gtex_adipose"]]
colsD = [C_ORANGE, C_BLUE, C_GREEN, C_YELLOW]
bars = axD.bar(pops, valsD, color=colsD, width=0.55, edgecolor="#222222",
               lw=0.6, zorder=2, alpha=0.9)
axD.axhline(0, color="#333333", lw=0.7, zorder=3)
axD.set_ylabel("Adipose score ($z$)")
axD.set_ylim(-1.6, 0.95)
for b, val in zip(bars, valsD):
    offset = 0.05 if val >= 0 else -0.12
    va = "bottom" if val >= 0 else "top"
    axD.text(b.get_x() + b.get_width() / 2, val + offset, f"{val:+.2f}",
             ha="center", va=va, fontsize=6.5, fontweight="bold")

axE = fig.add_subplot(gs[1, 2])
panel_tag(axE, "E")
axE.set_xlim(0, 10)
axE.set_ylim(0, 10)
axE.axis("off")


def node_box(x, y, w, h, text, fc):
    axE.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                                 fc=fc, ec="#333333", lw=0.7, zorder=2))
    axE.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=6.8, zorder=3, linespacing=1.2)


def conn_arrow(start, end):
    axE.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                  mutation_scale=8, lw=0.8,
                                  color="#333333", zorder=1))


node_box(0.8, 8.2, 8.4, 1.4,
         f"Full Transcriptome: {full_acc:.1f}% "
         f"(AUC {ablation['before']['auc']:.3f})", "#EAEAEA")
conn_arrow((3.0, 8.2), (2.5, 7.0))
conn_arrow((7.0, 8.2), (7.5, 7.0))
node_box(0.2, 4.8, 4.5, 2.2,
         f"Adipose-Stromal Axis\nPanel C: {panel_accs['panelC']:.1f}%\n"
         "- Adipocyte loss\n- Stroma exp.", "#D9EAF7")
node_box(5.3, 4.8, 4.5, 2.2,
         f"Residual Signal\n{resid_acc:.1f}% (AUC {resid_auc:.3f})\n"
         "- Tumor intrinsic\n- Non-adipose", "#FBE8D6")
conn_arrow((2.45, 4.8), (2.45, 3.2))
node_box(0.2, 1.6, 4.5, 1.6,
         "Tissue Architecture\nLobular < Ductal\n(p <= 1.7e-4)", "#EBF4FA")

fig.savefig(os.path.join(FIG, "figure1.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "figure1.png"), bbox_inches="tight", dpi=300)
plt.close(fig)

figS, axS = plt.subplots(figsize=(6.2, 3.5))
panel_tag(axS, "S1")
axS.grid(axis="y", color=C_GRID, lw=0.6, zorder=0)

n, bins, _ = axS.hist(
    null_draws,
    bins=30,
    color="#D8DCE3",
    edgecolor="#6B7280",
    lw=0.5,
    density=False,
    zorder=2,
)
max_h = max(n)
axS.set_ylim(0, max_h * 1.22)
axS.set_xlim(80, 101)

lines_data = [
    (null_med, C_NULL, "--", 1.1),
    (panel_accs["panelA"], C_BLUE, ":", 1.2),
    (panel_accs["panelB"], C_BLUE, "-.", 1.2),
    (panel_accs["randDE"], C_YELLOW, "--", 1.2),
    (panel_accs["panelC"], C_GREEN, "-", 1.4),
    (full_acc, "#111111", "-", 1.4),
]
for x, c, ls, lw in lines_data:
    axS.axvline(x, color=c, lw=lw, ls=ls, zorder=3)

legend_handles = [
    Line2D([0], [0], color=C_NULL, lw=1.2, ls="--",
           label=f"Null median ({null_med:.1f}%)"),
    Line2D([0], [0], color=C_BLUE, lw=1.2, ls=":",
           label=f"Panel A ({panel_accs['panelA']:.1f}%)"),
    Line2D([0], [0], color=C_BLUE, lw=1.2, ls="-.",
           label=f"Panel B ({panel_accs['panelB']:.1f}%)"),
    Line2D([0], [0], color=C_YELLOW, lw=1.2, ls="--",
           label=f"randDE ({panel_accs['randDE']:.1f}%)"),
    Line2D([0], [0], color=C_GREEN, lw=1.4, ls="-",
           label=f"Panel C ({panel_accs['panelC']:.1f}%)"),
    Line2D([0], [0], color="#111111", lw=1.4, ls="-",
           label=f"Full model ({full_acc:.1f}%)"),
]
axS.legend(
    handles=legend_handles, loc="upper left", frameon=True,
    facecolor="white", edgecolor="#D0D0D0", framealpha=0.95,
    fontsize=7.3, borderpad=0.6, labelspacing=0.45, handlelength=2.0,
)

null_max = float(null_draws.max())
n_ge = int((null_draws >= panel_accs["panelC"]).sum())
p_emp = (n_ge + 1) / (len(null_draws) + 1)
axS.annotate(
    f"Panel C & Full\n{n_ge}/{len(null_draws)} null >= "
    f"{panel_accs['panelC']:.1f}%\n"
    f"(empirical P = {p_emp:.3f}; null max = {null_max:.1f}%)",
    xy=(98.8, max_h * 0.75),
    xytext=(93.5, max_h * 0.95),
    fontsize=6.8, color="#222222", ha="center", va="center",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="#FAFAFA",
              edgecolor="#CCCCCC", lw=0.6),
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.15",
                    color="#444444", lw=0.8),
    zorder=4,
)

axS.set_xlabel("Accuracy (%) - 10-gene random panels (B = 500 draws)")
axS.set_ylabel("Count")
figS.tight_layout()
figS.savefig(os.path.join(FIG, "figureS1.pdf"), bbox_inches="tight")
figS.savefig(os.path.join(FIG, "figureS1.png"), bbox_inches="tight",
             dpi=300)
plt.close(figS)

print("Figuras 1 e S1 regeneradas com sucesso.")
```

---

## pipeline_utils.py

```python
# -*- coding: utf-8 -*-
"""E15-GUARD: utilitário para extração segura de painéis de genes.

Previne a classe de bug E15/E17: usar índices do eixo FILTRADO para
indexar a matriz FULL (ou vice-versa). Todos os scripts que extraem
subconjuntos de genes devem usar extract_gene_matrix() em vez de
indexação direta.
"""
import numpy as np


def strip_version(ensg):
    return ensg.split(".")[0]


def build_symbol_index(symbols_full, keep_mask):
    """symbol -> primeira posição no eixo FILTRADO."""
    kept = [strip_version(s) for s, k in zip(symbols_full, keep_mask) if k]
    idx = {}
    for i, s in enumerate(kept):
        idx.setdefault(s, i)
    return idx


def extract_gene_matrix(mat_full_filtered, symbol_pos, wanted):
    """Extrai amostras x genes(wanted) com guard de alinhamento.

    mat_full_filtered : matriz no eixo FILTRADO (amostras x genes)
    symbol_pos        : dict símbolo -> índice no eixo filtrado
    wanted            : lista de símbolos desejados

    Levanta KeyError se algum gene não existir no eixo filtrado.
    Levanta ValueError se houver descasamento de shape.
    """
    missing = [g for g in wanted if g not in symbol_pos]
    if missing:
        raise KeyError(
            f"E15-GUARD: símbolos ausentes no eixo filtrado: {missing}")
    cols = [symbol_pos[g] for g in wanted]
    if max(cols) >= mat_full_filtered.shape[1]:
        raise ValueError(
            f"E15-GUARD: índice {max(cols)} fora dos limites "
            f"(matriz tem {mat_full_filtered.shape[1]} colunas)")
    return np.ascontiguousarray(mat_full_filtered[:, cols])


def pick_nearest_non_panel(lfc, panel_cols, panel_symbols, sym_of):
    """Para cada coluna do painel, o gene mais próximo por |log2FC| que
    NÃO pertence ao painel.

    Exclusão por índice de coluna E por símbolo — guarda da classe de bug
    E15 (namespace errado): comparar IDs Ensembl contra um conjunto de
    símbolos é sempre True e faz o próprio gene do painel passar no
    filtro (o randDE degenerava no painel original — ver
    03_nested_compression.py antes do fix).

    sym_of : dict índice -> símbolo (ou array indexável por índice;
             entradas None/"" são tratadas como não-painel).
    """
    fat_cols = {int(c) for c in panel_cols}
    fat_syms = set(panel_symbols)
    rd = []
    for gp in panel_cols:
        for c in np.argsort(np.abs(lfc - lfc[gp])):
            c = int(c)
            if c in fat_cols:
                continue
            s = sym_of.get(c) if isinstance(sym_of, dict) else sym_of[c]
            if s is not None and s != "" and s in fat_syms:
                continue
            if c not in rd:
                rd.append(c)
                break
    return rd
```

---

## run_all.py

```python
# -*- coding: utf-8 -*-
"""Public reproduction orchestrator. Stages skip when outputs exist."""
import math
import os, sys, subprocess, argparse, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)

PIPELINE = ["01_download_data.py", "01b_download_estimate.py",
            "02_build_matrices.py",
            "03_nested_compression.py", "04_adipose_ablation.py",
            "05_cross_cohort_axis.py", "06_estimate_baseline.py",
            "07_transfer.py", "08_clinical_subtypes.py",
            "09_clinical_models.py"]

def get(d, path):
    cur = d
    for part in path.split("/"):
        cur = cur[part]
    return float(cur)

# verification manifest (path uses / to walk nested dicts)
EXPECTED = {
    "nested_compression.json": {
        "full/acc": 0.989, "panelA/acc": 0.947, "panelB/acc": 0.958,
        "panelC/acc": 0.987, "k1/acc": 0.970, "null/acc": 0.920,
        "randDE/acc": 0.981},
    "04_adipose_ablation.json": {
        "before/acc": 0.991, "after/acc": 0.909},
    "05_cross_cohort_axis.json": {
        "tcga_tumor": -1.374, "tcga_adjacent_normal": 0.400,
        "gtex_breast": 0.388, "gtex_adipose": 0.696},
    "06_estimate_baseline.json": {
        "cv_acc/stromal_only": 0.819, "cv_acc/panel": 0.984,
        "cv_acc/stromal_plus_panel": 0.981},
    "07_transfer.json": {
        "gtex_normal_frac": 0.996, "metabric_frac_tumor": 0.935},
    "08_clinical_subtypes.json": {
        "clinical_bh/histology": 1.34e-04, "clinical_bh/ER_STATUS": 1.61e-02,
        "clinical_bh/PR_STATUS": 3.97e-02,
        "pam50_chi2_with_normal": 2.25e-03,
        "pam50_chi2_without_normal": 4.04e-02},
    "09_clinical_models.json": {
        "overall_capture/LR": 0.924, "overall_capture/AXIS": 0.838,
        "hist_capture/LR/lobular/mean": 0.836,
        "hist_capture/AXIS/lobular/mean": 0.644,
        "pam50_full_chi2/lr": 2.25e-03,
        "pam50_er_validation/er_positive_basal": 0.027,
        "pam50_er_validation/er_positive_lumA": 0.907,
        "split_half/LR@0.5/half_B": 0.930,
        "split_half/AXIS@t*/half_B": 0.848,
        "concordance_heldout": 0.889},
}
TOL = 0.02

# script -> result filename (03 writes nested_compression.json, unprefixed)
OUT_NAMES = {"03_nested_compression.py": "nested_compression.json"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    for script in PIPELINE:
        out = os.path.join(RES, OUT_NAMES.get(script, script[:-3] + ".json"))
        if script.startswith(("01", "02")):
            print(f"[run ] {script} (idempotente: baixa apenas o que falta)",
                  flush=True)
        elif os.path.exists(out) and not args.force:
            print(f"[done] {script}")
            continue
        r = subprocess.run([sys.executable, os.path.join(HERE, script)],
                           cwd=HERE)
        if r.returncode != 0:
            raise SystemExit(f"FAILED: {script}")

    print("\n=== verification (headline numbers) ===")
    n_ok = n_bad = 0
    for jf, checks in EXPECTED.items():
        p = os.path.join(RES, jf)
        if not os.path.exists(p):
            print(f"  [skip] {jf}")
            continue
        d = json.load(open(p, encoding="utf-8"))
        for path, exp in checks.items():
            v = get(d, path)
            # accuracies: absolute tolerance; p-values/chi2/BH (<<1):
            # relative tolerance so the check stays meaningful at that scale
            if abs(exp) < 0.05:
                ok = math.isclose(v, exp, rel_tol=0.02)
            else:
                ok = abs(v - exp) <= TOL
            n_ok += ok
            n_bad += (not ok)
            print(f"  [{'PASS' if ok else 'FAIL'}] {jf}::{path} = "
                  f"{v:.6g} (esperado {exp})")
    print(f"\n{n_ok} PASS / {n_bad} FAIL "
          f"({time.time()-t0:.0f}s)")
    sys.exit(0 if n_bad == 0 else 1)

if __name__ == "__main__":
    main()
```

---

## test_guard.py

```python
# -*- coding: utf-8 -*-
"""E15-GUARD: extrai as colunas corretas da matriz."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_utils import extract_gene_matrix

# matriz: 4 amostras x 6 genes, gene j na coluna j tem valor j em todas as amostras
mat = np.tile(np.arange(6, dtype=np.float64), (4, 1))
symbols = [f"G{j}" for j in range(6)]
sym_pos = {s: i for i, s in enumerate(symbols)}

wanted = ["G1", "G3", "G5"]
result = extract_gene_matrix(mat, sym_pos, wanted)

print(f"result:\n{result}")
assert result.shape == (4, 3)
assert (result[:, 0] == 1.0).all(), "G1 incorreto"
assert (result[:, 1] == 3.0).all(), "G3 incorreto"
assert (result[:, 2] == 5.0).all(), "G5 incorreto"
print("E15-GUARD: PASS — colunas extraídas corretamente")
```

---

## test_namespace_guard.py

```python
# -*- coding: utf-8 -*-
"""E15-GUARD (namespace): regressão do bug do randDE.

O bug histórico: em 03_nested_compression.py, o filtro do randDE comparava
IDs Ensembl (gk) contra um conjunto de SÍMBOLOS (PANEL_B) — `gk[c] not in
fat` era sempre True, e o gene mais próximo por |log2FC| do próprio gene
do painel (distância 0) passava no filtro. Resultado: randDE ≡ panelB,
arrays per-split idênticos, controle de efeito de tamanho invalidado.

Este teste reproduz a armadilha sinteticamente: com os namespaces
deliberadamente trocados, pick_nearest_non_panel() deve EXCLUIR as colunas
do painel pela coluna (não pelo símbolo) e devolver os vizinhos corretos;
e com o símbolo correto, deve excluir também homólogos de outro índice.
"""
import numpy as np

from pipeline_utils import pick_nearest_non_panel

# 6 genes; colunas 0 e 1 são do painel ("P0", "P1").
# lfc[0] = 2.0 (o painel é o próprio alvo, distância 0),
# lfc[2] = 1.9 (vizinho não-painel mais próximo de 0),
# lfc[3] = 0.5, lfc[4] = 0.4, lfc[5] = 3.0.
lfc = np.array([2.0, 0.05, 1.9, 0.5, 0.4, 3.0])
panel_cols = np.array([0, 1])
panel_symbols = {"P0", "P1"}

# --- caso 1: namespaces trocados (o bug histórico) ---
# sym_of contém IDs Ensembl, que NUNCA casam com os símbolos do painel.
# A implementação com bug devolveria [0, 1] (o painel inteiro).
sym_of_broken = {i: f"ENSG{i:08d}" for i in range(6)}
rd = pick_nearest_non_panel(lfc, panel_cols, panel_symbols, sym_of_broken)
assert set(rd) == {2, 4}, (
    f"REGRESSÃO namespace: esperado {{2, 4}}, obtido {rd} "
    "(o painel vazou no filtro — classe E15)")
assert 0 not in rd and 1 not in rd, "gene do painel selecionado"

# --- caso 2: símbolo correto, homólogo em índice diferente ---
# aqui a coluna 3 (símbolo "P0", homólogo do painel em outro índice) é o
# vizinho MAIS PRÓXIMO de P1: só a exclusão por símbolo impede a escolha.
lfc_b = np.array([2.0, 0.05, 1.9, 0.045, 0.4, 3.0])
sym_of_dup = {0: "P0", 1: "P1", 2: "N2", 3: "P0", 4: "N4", 5: "N5"}
rd2 = pick_nearest_non_panel(lfc_b, panel_cols, panel_symbols, sym_of_dup)
assert 3 not in rd2, "homólogo por símbolo não excluído"
assert set(rd2) == {2, 4}, f"esperado {{2, 4}}, obtido {rd2}"

# --- caso 3: exclusão por símbolo NÃO deve excluir não-painel ---
sym_of_clean = {i: f"N{i}" for i in range(6)}
rd3 = pick_nearest_non_panel(lfc, panel_cols, panel_symbols, sym_of_clean)
assert set(rd3) == {2, 4}, f"esperado {{2, 4}}, obtido {rd3}"

# --- caso 4: sem símbolo (None), exclusão cai na coluna ---
sym_of_none = {i: None for i in range(6)}
rd4 = pick_nearest_non_panel(lfc, panel_cols, panel_symbols, sym_of_none)
assert set(rd4) == {2, 4}, f"esperado {{2, 4}}, obtido {rd4}"

print("E15-GUARD (namespace/randDE): PASS — painel nunca vaza no filtro, "
      "por coluna ou por símbolo")
```

---

## test_pipeline_utils.py

```python
# -*- coding: utf-8 -*-
"""E15-GUARD regression test: verifies that extract_gene_matrix raises
on missing symbols and returns correctly aligned columns."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_utils import extract_gene_matrix

# synthetic: 100 full genes, keep 50 (even indices), 20 samples
n_full, n_samples, n_keep = 100, 20, 50
symbols_full = np.array([f"GENE{i:04d}" for i in range(n_full)])
keep_mask = np.zeros(n_full, dtype=bool)
keep_mask[::2] = True  # keep even indices
mat = np.random.default_rng(0).random((n_samples, n_full))
sym_pos = {}
for i, s in enumerate(symbols_full[keep_mask]):
    sym_pos.setdefault(s, i)
mat_filtered = mat[:, keep_mask]

# test: extracts correct columns
wanted = ["GENE0010", "GENE0020", "GENE0030"]
cols = [sym_pos[g] for g in wanted]
expected = mat_filtered[:, cols]
result = extract_gene_matrix(mat_filtered, sym_pos, wanted)
assert result.shape == (n_samples, 3), f"shape {result.shape}"
assert (result == expected).all(), "values mismatch"

# test: raises KeyError on a symbol absent from the filtered axis
try:
    extract_gene_matrix(mat_filtered, sym_pos, ["GENE9999_MISSING"])
    raise AssertionError("expected KeyError for missing symbol")
except KeyError:
    pass
print(f"E15-GUARD regression test: PASS")
```

