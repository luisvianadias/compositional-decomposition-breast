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
