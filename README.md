# Adipose-associated compositional axis in bulk breast tumor-normal separation

Reproduction code for the manuscript quantifying how much of the bulk
transcriptomic boundary between normal and tumor breast tissue is captured
by an adipose-associated transcriptional axis.

All analyses run on scikit-learn / SciPy / pandas only (no proprietary
components). Primary classifier: logistic regression, with Random Forest,
linear SVM and histogram gradient boosting as robustness checks.

## Data (all public)

| Cohort | Source | n |
|---|---|---|
| TCGA-BRCA STAR counts | GDC API (open access; GENCODE v36) | 631 (518 tumor / 113 adjacent normal) |
| GTEx v10 breast | GTEx Portal counts-by-tissue (public GCS bucket) | 514 donors |
| GTEx v10 adipose (subcutaneous) | same | positive compositional control |
| METABRIC | cBioPortal API, study `brca_metabric` (raw HT-12 v3 + clinical) | 1,980 tumors |

Access dates and versions: see `data/PROVENANCE.md`.

## Layout

```
scripts/
  01_download_data.py     # GDC manifest + counts + GTEx gct.gz (skip if present)
  02_build_matrices.py    # TCGA count matrix (npz) + symbol map
  03_nested_compression.py# Sec 2.1-2.2: compressibility curve, panels,
                          #            randDE control, null control (B=500)
  04_adipose_ablation.py  # Sec 2.3: residualization + PC1 (leakage-free)
  05_cross_cohort_axis.py # Sec 2.4: adipose-axis score, TCGA/GTEx/adipose
  06_estimate_baseline.py # Sec 2.6: ESTIMATE port + panel comparison
  07_transfer.py          # Sec 2.4: TCGA-trained panel applied to GTEx/METABRIC
  08_clinical_subtypes.py # Sec 2.5: clinical BH family + PAM50 subtypes
run_all.py
```

## Reproduction

```bash
python -m pip install -r requirements.txt
python scripts/run_all.py            # downloads data on first run (~10 min)
```

Every script writes JSON outputs to `results/` and uses fixed seeds
(42/123/777). The expected headline numbers (abstract) are embedded in
`run_all.py` as a verification manifest.

## Method summary

- Normalization: log2(CPM+1); genes with mean raw count > 10.
- Panels: 12 canonical adipocyte-associated genes (Panel A), +3
  (Panel B), +COL10A1/MMP11 (Panel C); randDE = non-adipocyte genes
  matched by training-fold |log2FC|; null = uniformly random 10-gene
  panels (B = 500).
- Nested evaluation: feature selection, z-scoring, residualization and
  randDE matching estimated exclusively within training folds
  (repeated stratified 5-fold, 5 x 3 = 15 paired splits).
- Residualization: per-gene linear coefficient against the adipose
  score, training-only.
- ESTIMATE: faithful port of the reference R implementation.
- Multiplicity: Benjamini-Hochberg over the clinical family; PAM50
  subtyping via genefu ssp2006 centroids (Spearman rule).

## License and provenance

Data: TCGA Research Network / GDC; GTEx Consortium v10; METABRIC via
cBioPortal (Curtis et al. 2012). Citation requirements in
`data/PROVENANCE.md`. Code: MIT.
