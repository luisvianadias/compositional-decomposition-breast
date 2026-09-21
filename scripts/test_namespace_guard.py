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
