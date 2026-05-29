"""
02_aggregate_consensus_grn.py
=============================
Aggregate 100 independent pySCENIC GRN runs into a single consensus GRN.

Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
  "...the inference step (pyscenic grn) was run 100 times with different random
   seeds... transcription factor-target edges were retained only if present in
   90 or more of 100 runs. Edge importance scores were averaged across runs to
   yield a consensus GRN..."

Input:
  ./grn_results/grn_results_cells_plus_mi_run_{1..100}.csv
    - Output of bash/01_pyscenic_grn_array.sh
    - Each file contains columns: TF, target, importance

Output:
  cells_plus_mi_consensus_grn_results.csv
    - Consensus GRN: TF-target edges present in ≥90/100 runs
    - Columns: TF, target, importance (mean across runs in which edge appeared)

Dependencies:
  pandas, dask

Usage:
  python 02_aggregate_consensus_grn.py

Notes:
  - Uses dask for memory-efficient reading of 100 CSV files in parallel
  - The 90% threshold (≥90/100 runs) is the primary stability filter
  - Average importance is computed only over runs in which the edge appeared
    (i.e. edges with count < 90 are filtered before the importance average
    is exported, so avg_importance reflects the mean over present-run estimates)
  - Output feeds into bash/02_pyscenic_ctx_aucell.sh (ctx step)
"""

import pandas as pd
import dask.dataframe as dd

# ── Parameters ────────────────────────────────────────────────────────────────

# Directory containing the 100 GRN run output CSVs
# (produced by bash/01_pyscenic_grn_array.sh)
results_path = "./grn_results"

# Filename pattern matching the output of pyscenic_grn_array.sh
grn_files = [
    f"{results_path}/grn_results_cells_plus_mi_run_{i}.csv"
    for i in range(1, 101)
]

# Minimum number of runs an edge must appear in to be retained
# 90/100 = 90% consensus threshold (see Methods)
threshold_runs = int(0.9 * 100)  # = 90

# ── Load all 100 GRN run files ─────────────────────────────────────────────────

# Use dask for memory-efficient parallel loading
ddf = dd.read_csv(grn_files)

# ── Aggregate across runs ──────────────────────────────────────────────────────

# Count how many runs each TF-target edge appears in,
# and sum importance scores across runs
aggregated_ddf = ddf.groupby(["TF", "target"]).agg(
    count=("TF", "size"),             # number of runs edge appeared in
    importance_sum=("importance", "sum")
).reset_index()

# Compute (triggers dask execution)
aggregated = aggregated_ddf.compute()

# Compute mean importance score across runs in which the edge appeared
aggregated["avg_importance"] = aggregated["importance_sum"] / aggregated["count"]

# ── Apply consensus threshold ──────────────────────────────────────────────────

# Retain only edges present in ≥90% of runs
filtered = aggregated[aggregated["count"] >= threshold_runs]

# Save full filtered output including count (useful for QC)
filtered.to_csv("cells_plus_mi_filtered_grn_results.csv", index=False)

# ── Save consensus GRN ────────────────────────────────────────────────────────

# Produce final consensus GRN with standard column names (TF, target, importance)
# This is the format expected by pyscenic ctx (bash/02_pyscenic_ctx_aucell.sh)
# and by the R co-activity network analysis (scripts/06_SCENIC_networks.R)
df = filtered[["TF", "target", "avg_importance"]].copy()
df["importance"] = df["avg_importance"]
df = df.drop(columns=["avg_importance"])

df.to_csv("cells_plus_mi_consensus_grn_results.csv", index=False)

print(f"Consensus GRN: {len(df)} edges retained from {threshold_runs}/100 run threshold")
print(f"Unique TFs: {df['TF'].nunique()}")
print(f"Unique targets: {df['target'].nunique()}")
