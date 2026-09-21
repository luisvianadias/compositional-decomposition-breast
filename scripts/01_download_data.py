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
if not exists(est, 100000):
    print("  [get ] estimate_r.tar.gz (SourceForge)")
    urllib.request.urlretrieve(
        "https://sourceforge.net/projects/estimateproject/files/latest/"
        "download", est)
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
