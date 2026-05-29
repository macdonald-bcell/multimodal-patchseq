#!/bin/bash
# ==============================================================================
# 01_pyscenic_grn_array.sh
# ==============================================================================
# Run pySCENIC GRN inference 100 times with different random seeds using a
# SLURM job array. Each of 10 array tasks runs 10 sequential GRN iterations,
# yielding 100 independent GRN estimates total.
#
# Corresponds to Methods: "SCENIC regulon inference and co-activity networks"
#   "...the inference step (pyscenic grn) was run 100 times with different
#    random seeds using SLURM job arrays..."
#
# Input:
#   20260113_cells_plus_scenic.loom
#     - Produced by python/01_h5ad_to_loom.py
#   ./pyscenic_auxiliary/allTFs_hg38.txt
#     - Human TF list from pySCENIC resources (https://scenic.aertslab.org/)
#
# Output:
#   ./grn_results/grn_results_cells_plus_mi_run_{1..100}.csv
#     - One GRN per run; columns: TF, target, importance
#     - Used as input to python/02_aggregate_consensus_grn.py
#
# Resource requirements (per array task):
#   CPUs: 20 (15 used by pySCENIC per iteration, leaving headroom)
#   Memory: 160G
#   Time: 2 days (10 sequential GRN runs per task)
#
# Usage:
#   sbatch 01_pyscenic_grn_array.sh
#
# Notes:
#   - pySCENIC GRN inference is stochastic; running 100 times and taking the
#     consensus (edges in ≥90/100 runs) reduces sensitivity to any single run
#   - --num_workers 15 is set below the --cpus-per-task 20 to avoid overloading
#     the node while leaving headroom for the SLURM overhead
#   - Array index 0-9: task 0 runs iterations 1-10, task 1 runs 11-20, etc.
#   - pySCENIC is run via Apptainer (Singularity) container to ensure
#     reproducible environment; pyscenic.sif must be in the working directory
# ==============================================================================

#SBATCH --job-name=pyscenic_grn_array
#SBATCH --output=pyscenic_%A_%a.out       # %A=job ID, %a=array task ID
#SBATCH --error=pyscenic_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-9                       # 10 array tasks × 10 iterations = 100 runs

echo "Running on hostname $(hostname)"
echo "Starting run at: $(date)"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

module load apptainer/1.2.4

# ── Input files ───────────────────────────────────────────────────────────────

INPUT_LOOM="20260113_cells_plus_scenic.loom"    # from python/01_h5ad_to_loom.py
ALL_TFS="./pyscenic_auxiliary/allTFs_hg38.txt"  # human TF list

OUTPUT_DIR="./grn_results"
mkdir -p $OUTPUT_DIR

# ── Calculate iteration range for this array task ─────────────────────────────

# Array task 0 → iterations 1-10
# Array task 1 → iterations 11-20
# ...
# Array task 9 → iterations 91-100
START=$((SLURM_ARRAY_TASK_ID * 10 + 1))
END=$((START + 9))

echo "Running iterations i = $START to $END"

# ── Run GRN inference ─────────────────────────────────────────────────────────

for (( i=$START; i<=$END; i++ ))
do
    OUTPUT_ADJ="$OUTPUT_DIR/grn_results_cells_plus_mi_run_$i.csv"
    echo "Running GRN inference for i=$i -> $OUTPUT_ADJ"

    # Output: CSV with columns TF, target, importance
    apptainer exec pyscenic.sif pyscenic grn \
        $INPUT_LOOM \
        $ALL_TFS \
        -o $OUTPUT_ADJ \
        --num_workers 15
done

echo "Job finished with exit code $? at: $(date)"
