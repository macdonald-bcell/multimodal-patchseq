"""
01_h5ad_to_loom.py
==================
Convert a Seurat-exported AnnData (h5ad) object to loom format for use as
input to pySCENIC GRN inference.

Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
  "Gene regulatory network (GRN) inference was performed using pySCENIC on
   the quality-controlled scRNA-seq dataset."

Input:
  20260113_cells_plus_symbs.h5ad
    - AnnData object exported from R via sceasy::convertFormat()
    - Contains decontaminated raw counts (decontX-corrected)
    - Gene names are HGNC symbols (not Ensembl IDs) — required by pySCENIC
    - Cell metadata includes: CellID, Study, celltype, celltype_stage,
      Differentiation
    - Note: eGFP is retained as an explicit identity marker

Output:
  20260113_cells_plus_scenic.loom
    - Loom file formatted for pySCENIC input
    - Genes (rows) × Cells (columns), raw counts
    - Required by bash/01_pyscenic_grn_array.sh and bash/02_pyscenic_ctx_aucell.sh

Dependencies:
  numpy, scanpy, loompy

Environment:
  Conda environment: pyscenic2
  (pySCENIC requires its own environment due to dependency conflicts)

Usage:
  conda activate pyscenic2
  python 01_h5ad_to_loom.py

Notes:
  - pySCENIC requires raw (unnormalized) counts as input
  - Gene names must be HGNC symbols; Ensembl IDs are not recognized by the
    TF list or rankings databases
  - The h5ad was produced in R using sceasy::convertFormat() on the
    cells.plus Seurat object (see scripts/02_integrate_patchseq_datasets.R)
"""

import numpy as np
import scanpy as sc
import loompy as lp
import pandas as pd

# ── Parameters ────────────────────────────────────────────────────────────────

INPUT_H5AD  = "20260113_cells_plus_symbs.h5ad"
OUTPUT_LOOM = "20260113_cells_plus_scenic.loom"

# ── Load AnnData ──────────────────────────────────────────────────────────────

print(f"Loading {INPUT_H5AD}...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# ── QC: genes detected per cell ───────────────────────────────────────────────

# Print percentiles of gene detection as a sanity check before pySCENIC
# pySCENIC recommends cells with ≥200 detected genes; check distribution here
nGenesDetectedPerCell = np.array(np.sum(adata.X > 0, axis=1)).flatten()
percentile_values = np.percentile(nGenesDetectedPerCell, [1, 5, 10, 50, 100])
percentiles = pd.Series(percentile_values, index=["1%", "5%", "10%", "50%", "100%"])
print("\nGenes detected per cell (percentiles):")
print(percentiles)

# ── Build loom file ───────────────────────────────────────────────────────────

# pySCENIC expects:
#   - expression matrix: genes (rows) × cells (columns), dense
#   - row attributes: Gene (gene names)
#   - column attributes: CellID (cell barcodes), nGene, nUMI

# NOTE: loompy expects genes × cells, but AnnData stores cells × genes
# so we transpose here
expression_matrix_dense = adata.X.transpose().toarray()

row_attrs = {
    "Gene": np.array(adata.var.index),  # HGNC gene symbols
}

col_attrs = {
    "CellID": np.array(adata.obs.index),
    # nGene and nUMI are metadata pySCENIC may use for QC reporting
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI":  np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

print(f"\nWriting loom file: {OUTPUT_LOOM}")
lp.create(OUTPUT_LOOM, expression_matrix_dense, row_attrs, col_attrs)
print("Done.")
print(f"  Matrix shape: {expression_matrix_dense.shape[0]} genes × "
      f"{expression_matrix_dense.shape[1]} cells")
