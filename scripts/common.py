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
