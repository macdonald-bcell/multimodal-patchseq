"""
03_auc_extraction.py
====================
Extract regulon AUC scores from pySCENIC output loom file and save in
formats compatible with downstream R analysis (Seurat import via sceasy).

Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
  "...per-cell regulon activity scores (AUC values)."

Input:
  20260114_cells_plus_mi_pyscenic_output.loom
    - Output of bash/02_pyscenic_ctx_aucell.sh (pyscenic aucell step)
    - Contains regulon AUC scores per cell and regulon gene membership

Output:
  cells_plus_mi_auc_matrix.csv
    - Cells × regulons AUC matrix (primary output for Seurat import)
    - Loaded into Seurat as a separate "AUC" assay in R
  cells_plus_mi_regulon_thresholds.csv
    - Default AUC threshold per regulon from pySCENIC metadata
    - Used for binarizing AUC scores if needed
  20260115_adata_cells_plus_mi_SCENIC.h5ad
    - AnnData with AUC scores added to obsm['RegulonsAUC']
    - Useful for Python-side visualization

Dependencies:
  numpy, scanpy, loompy, pandas, json, zlib, base64, csv

Environment:
  Conda environment: pyscenic2

Usage:
  conda activate pyscenic2
  python 03_auc_extraction.py

Notes on regulon name formatting:
  - pySCENIC outputs regulon names as e.g. "SREBF1(+)"
  - These are reformatted to "SREBF1_(+)" by replacing "(" with "_("
  - This avoids parsing issues with parentheses in downstream R code
  - The reformatting is applied consistently to auc_mtx columns,
    the regulons structured array, and the regulonThresholds metadata

⚠ Bug fix applied:
  - Original script attempted `del adata.uns['RegulonSignatures']` at the
    end, but RegulonSignatures was never added to adata.uns (the line adding
    it was commented out). This caused a KeyError. The deletion has been
    removed.
"""

import numpy as np
import scanpy as sc
import loompy as lp
import pandas as pd
import json
import zlib
import base64
import csv

# ── Parameters ────────────────────────────────────────────────────────────────

INPUT_LOOM        = "20260114_cells_plus_mi_pyscenic_output.loom"
OUTPUT_AUC_CSV    = "cells_plus_mi_auc_matrix.csv"
OUTPUT_THRESH_CSV = "cells_plus_mi_regulon_thresholds.csv"
OUTPUT_H5AD       = "20260115_adata_cells_plus_mi_SCENIC.h5ad"

# ── Load pySCENIC output loom ─────────────────────────────────────────────────

print(f"Connecting to {INPUT_LOOM}...")
lf = lp.connect(INPUT_LOOM, mode='r+', validate=False)

# Extract metadata (contains regulon thresholds and embedding info)
meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))

# Extract AUC matrix: cells × regulons
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

# Extract regulon gene membership (structured array: genes × regulons)
regulons = lf.ra.Regulons

print(f"AUC matrix shape: {auc_mtx.shape} (cells × regulons)")

# ── Reformat regulon names ─────────────────────────────────────────────────────

# Replace "(" with "_(" in regulon names to avoid downstream parsing issues
# e.g. "SREBF1(+)" → "SREBF1_(+)"
# Applied to: AUC column names, regulons structured array, threshold metadata
auc_mtx.columns = auc_mtx.columns.str.replace('\\(', '_(', regex=True)
regulons.dtype.names = tuple([x.replace("(", "_(") for x in regulons.dtype.names])

# Update regulon names in threshold metadata to match
rt = meta['regulonThresholds']
for i, x in enumerate(rt):
    tmp = x.get('regulon').replace("(", "_(")
    x.update({'regulon': tmp})

# ── Save AUC matrix ───────────────────────────────────────────────────────────

# Primary output: cells × regulons AUC matrix
# Loaded into Seurat in R as:
#   auc_mtx <- read.csv("cells_plus_mi_auc_matrix.csv", row.names=1)
#   cells.plus[["AUC"]] <- CreateAssayObject(counts = t(auc_mtx))
print(f"Saving AUC matrix to {OUTPUT_AUC_CSV}...")
auc_mtx.to_csv(OUTPUT_AUC_CSV)

# ── Save regulon thresholds ───────────────────────────────────────────────────

# Default AUC thresholds from pySCENIC (for optional binarization)
print(f"Saving regulon thresholds to {OUTPUT_THRESH_CSV}...")
regulon_thresholds = meta.get('regulonThresholds', [])
with open(OUTPUT_THRESH_CSV, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["regulon", "defaultThresholdValue"])
    for entry in regulon_thresholds:
        writer.writerow([entry.get('regulon'), entry.get('defaultThresholdValue')])

# ── Add AUC scores to AnnData and save h5ad ───────────────────────────────────

# Reload the original h5ad (without pySCENIC normalization applied)
# and attach the AUC scores for Python-side visualization
print("Loading original h5ad to attach AUC scores...")
adata = sc.read_h5ad("20260113_cells_plus_symbs.h5ad")

# Align AUC matrix to adata cell order (some cells may differ due to filtering)
common_cells = adata.obs.index.intersection(auc_mtx.index)
n_dropped = len(adata.obs.index) - len(common_cells)
if n_dropped > 0:
    print(f"  ⚠ Warning: {n_dropped} cells in adata not found in AUC matrix — "
          f"these will have NaN AUC values")

adata = adata[common_cells].copy()
adata.obsm['RegulonsAUC'] = auc_mtx.loc[common_cells]

print(f"Saving {OUTPUT_H5AD}...")
adata.write_h5ad(OUTPUT_H5AD)

lf.close()
print("Done.")
print(f"  Regulons retained: {auc_mtx.shape[1]}")
print(f"  Cells retained: {len(common_cells)}")
