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

if not os.path.exists(DEST):
    print("[get ] estimate_r.tar.gz (~3.7 MB, SourceForge)")
    urllib.request.urlretrieve(
        "https://sourceforge.net/projects/estimateproject/files/latest/"
        "download", DEST)

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
