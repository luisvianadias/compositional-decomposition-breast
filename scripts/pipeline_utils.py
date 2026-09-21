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
