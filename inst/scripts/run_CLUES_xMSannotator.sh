#!/usr/bin/env bash
#===============================================================================
# xMSannotator CLUES Multi-Database Annotation Pipeline
# Version: 260309
# Author: CLUES Lab, Emory University
#
# Description: SLURM batch script for xMSannotator annotation across multiple
#              compound databases in dual polarity mode (positive/negative)
#
# Usage:
#   sbatch run_CLUES_xMSannotator.sh
#
# Array Tasks:
#   Tasks 1..N     = Positive mode (one per database)
#   Tasks N+1..2N  = Negative mode (one per database)
#   N = number of databases in db_mapfile
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM Configuration - Update email and array range as needed
#-------------------------------------------------------------------------------
#SBATCH --job-name=CLUES_xMSannot
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email@email.com>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --array=1-18
#SBATCH --error=CLUES_xMSannotator_%A_%a.err
#SBATCH --output=CLUES_xMSannotator_%A_%a.out
#SBATCH --partition=clues

set -euo pipefail

#===============================================================================
# USER CONFIGURATION - Modify these settings for each study
#===============================================================================

# Study identification
study_id="CLU0046.1"
lc_mode="C18"                       # C18 or HILIC

# File paths
output_dir="/clues/STUDY/~Data_Extraction/"
mapfile_dir="/clues/STUDY/2a_LC_Mapfiles/"
db_mapfile="/clues/.../Compound_Database_Mapfile_v3.xlsx"

# Processing parameters
n_cores=10                          # Must match --cpus-per-task
fc_threshold=5
blank_names="Water,Blank"           # Comma-separated

# R script location. Do not change.
r_script="/clues/~Scripts/xMSannotator_CLUES_MultiDB_Runscript.R"

#===============================================================================
# SCRIPT EXECUTION - No modifications needed below
#===============================================================================

# Create output directory
mkdir -p "$output_dir"

# Setup logging
task_id="${SLURM_ARRAY_TASK_ID:-1}"
log_file="${output_dir}/xMSannotator_${SLURM_ARRAY_JOB_ID:-local}_task${task_id}.log"
exec > >(tee -a "$log_file") 2>&1

# Print header
echo "==============================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] xMSannotator Annotation - Task ${task_id}"
echo "==============================================================================="
echo ""
echo "Study ID:        ${study_id}"
echo "LC Mode:         ${lc_mode}"
echo "Cores:           ${n_cores}"
echo "FC Threshold:    ${fc_threshold}"
echo "Blank Names:     ${blank_names}"
echo "Output:          ${output_dir}"
echo "DB Mapfile:      ${db_mapfile}"
echo ""

# Load R module
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Loading R module..."
module load R/4.4.0

# Execute R script
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting xMSannotator annotation..."
echo "==============================================================================="

Rscript "$r_script" \
    "$task_id" \
    "$lc_mode" \
    "$output_dir" \
    "$study_id" \
    "$mapfile_dir" \
    "Unadjusted" \
    "$fc_threshold" \
    "$n_cores" \
    "$db_mapfile" \
    "$blank_names"

exit_code=$?

# Print summary
echo ""
echo "==============================================================================="
if [[ $exit_code -eq 0 ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: Task ${task_id} complete"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAILED: Task ${task_id} (exit code: ${exit_code})"
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Runtime: ${SECONDS}s ($(( SECONDS / 60 ))m)"
echo "==============================================================================="

exit $exit_code
