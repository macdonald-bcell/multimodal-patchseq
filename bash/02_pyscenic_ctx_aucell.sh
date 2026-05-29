#!/bin/bash
# ==============================================================================
# 02_pyscenic_ctx_aucell.sh
# ==============================================================================
# Run pySCENIC cisTarget (ctx) and AUCell steps on the consensus GRN to
# generate motif-validated regulon definitions and per-cell AUC activity scores.
#
# Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
#   "...which was used as input for the cisTarget (ctx) and AUCell steps to
#    generate regulon definitions and per-cell regulon activity scores (AUC
#    values)."
#
# This script runs AFTER:
#   bash/01_pyscenic_grn_array.sh   → 100 GRN runs
#   python/02_aggregate_consensus_grn.py → consensus GRN
#
# Inputs:
#   cells_plus_mi_consensus_grn_results.csv
#     - Consensus GRN from python/02_aggregate_consensus_grn.py
#     - Edges present in ≥90/100 GRN runs, averaged importance scores
#   20260113_cells_plus_scenic.loom
#     - Expression matrix (from python/01_h5ad_to_loom.py)
#     - Required by ctx for masking dropouts and by aucell for AUC scoring
#   ./pyscenic_auxiliary/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#   ./pyscenic_auxiliary/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#     - Cisregulatory rankings databases (hg38, v10; from pySCENIC resources)
#   ./pyscenic_auxiliary/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
#     - Motif annotation table (from pySCENIC resources)
#
# Outputs:
#   cells_plus_mi_reg.csv
#     - Motif-validated regulon definitions (TF → target gene sets)
#     - Intermediate file; input to aucell step
#   20260114_cells_plus_mi_pyscenic_output.loom
#     - Final pySCENIC output loom containing per-cell AUC scores
#     - Used by python/03_auc_extraction.py to extract AUC matrix for R
#
# Resource requirements:
#   CPUs: 20 (ctx uses 15, aucell uses 20)
#   Memory: 160G
#   Time: 1 day
#
# Usage:
#   sbatch 02_pyscenic_ctx_aucell.sh
# ==============================================================================

#SBATCH --job-name=pyscenic_ctx_aucell
#SBATCH --output=pyscenic_%j.out
#SBATCH --error=pyscenic_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL

echo "Running on hostname $(hostname)"
echo "Starting run at: $(date)"

module load apptainer/1.2.4

# ── Input files ───────────────────────────────────────────────────────────────

INPUT_LOOM="20260113_cells_plus_scenic.loom"
OUTPUT_ADJ_AGG="cells_plus_mi_consensus_grn_results.csv"  # from python/02_aggregate_consensus_grn.py

RANKINGS="./pyscenic_auxiliary/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
          ./pyscenic_auxiliary/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOTATIONS="./pyscenic_auxiliary/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# ── Output files ──────────────────────────────────────────────────────────────

OUTPUT_REG="cells_plus_mi_reg.csv"                             # regulon definitions (ctx output)
OUTPUT_LOOM="20260114_cells_plus_mi_pyscenic_output.loom"      # AUC scores (aucell output)

# ── Step 1: cisTarget (ctx) ───────────────────────────────────────────────────
echo "Running pyscenic ctx..."

apptainer exec pyscenic.sif pyscenic ctx \
    $OUTPUT_ADJ_AGG \
    $RANKINGS \
    --annotations_fname $ANNOTATIONS \
    --expression_mtx_fname $INPUT_LOOM \
    --output $OUTPUT_REG \
    --mask_dropouts \
    --num_workers 15

echo "ctx step complete at: $(date)"

# ── Step 2: AUCell ────────────────────────────────────────────────────────────
echo "Running pyscenic aucell..."

apptainer exec pyscenic.sif pyscenic aucell \
    $INPUT_LOOM \
    $OUTPUT_REG \
    --output $OUTPUT_LOOM \
    --num_workers 20

echo "aucell step complete at: $(date)"
echo "Job finished with exit code $? at: $(date)"
