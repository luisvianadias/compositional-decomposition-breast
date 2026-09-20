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
