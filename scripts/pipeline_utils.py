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
