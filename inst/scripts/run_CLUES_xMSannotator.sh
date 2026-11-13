#!/usr/bin/env bash
#===============================================================================
# xMSannotator CLUES Multi-Database Annotation Pipeline
# Version: 1.2.2 (CLUES.xMSannotator)
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
#
# Paths: study_dir as in the Step 1-3 job scripts. Step 1 output is read from
# <study_dir>/2-Data_Extraction/<lc_mode>/ and the per-polarity sample mapfiles
# from <study_dir>/1-mzML/<lc_mode>/; output goes to
# <study_dir>/2-Data_Extraction/<lc_mode>/<mode>/7-<mode>_<study>_<DB>_CLUES_Annotation/.
#
# Needs: R/4.5.2 with CLUES.xMSannotator, readxl, writexl, dplyr, data.table.
#
# Feature grouping: output_dir is the Step 1 output_dir. When Step 1 v1.7.0+ has
# written {mode}/3-{mode}_{study_id}_Feature_Grouping/{mode}-{study_id}_Feature_Groups.txt,
# every task reads those groups (no WGCNA clustering; fc_threshold not applied).
# Without the file (older Step 1 output) each task filters by fc_threshold and
# runs WGCNA. The log line "Feature groups:" says which path was taken.
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM Configuration - Update email and array range as needed
#-------------------------------------------------------------------------------
#SBATCH --job-name=CLUES_xMSannot
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your-email@email.com>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
# Peak memory per task on CLU0127 (about 10,000 features x 350 injections, confirmed DB):
# about 4 GB with the Step 1 feature groups, about 10 GB with WGCNA.
# CLU0162 C18 (1,058 injections; 20,944 / 11,122 features, Step 1 groups, confirmed DB, 10 workers):
# 13.7 GB / 5.6 GB and under 1 min per task (Mac, run 88).
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

# Study folder, as in Step 1: Step 1 output in <study_dir>/2-Data_Extraction/<lc_mode>/, mapfiles in <study_dir>/1-mzML/<lc_mode>/
study_dir="/clues/CLU0046.1_STUDY/"
output_dir=""                       # Optional override ("" = derived): the Step 1 output_dir, if Step 1 used an override
mapfile_dir=""                      # Optional override ("" = derived)
db_mapfile="/clues/.../Compound_Database_Mapfile_v3.xlsx"

# Processing parameters
n_cores=10                          # Must match --cpus-per-task
fc_threshold=5                      # Water_fc >= threshold; used only without Step 1 feature groups
blank_names="Water,Blank"           # Comma-separated

# R script location. Do not change.
r_script="/clues/~Data_Extraction/xMSannotator/xMSannotator_CLUES_MultiDB_Runscript.R"

#===============================================================================
# SCRIPT EXECUTION - No modifications needed below
#===============================================================================

# Data paths: derived from study_dir unless set above; stop if Step 1 output or the mapfile folder is missing
task_id="${SLURM_ARRAY_TASK_ID:-1}"
output_dir="${output_dir:-${study_dir%/}/2-Data_Extraction/${lc_mode}/}"
mapfile_dir="${mapfile_dir:-${study_dir%/}/1-mzML/${lc_mode}/}"
for d in "$output_dir" "$mapfile_dir" "$db_mapfile"; do
    if [[ ! -e "$d" ]]; then
        echo "ERROR: not found: $d (has Step 1 run? check study_dir / lc_mode / db_mapfile, or set the override)" >&2
        exit 1
    fi
done

# Setup logging
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
echo "Mapfiles:        ${mapfile_dir}"
echo "DB Mapfile:      ${db_mapfile}"
echo ""

# Load R module
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Loading R module..."
module load R/4.5.2

# Execute R script
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting xMSannotator annotation..."
echo "==============================================================================="

# set -e would end the script at a failing Rscript before the summary below; capture the exit code instead
exit_code=0
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
    "$blank_names" || exit_code=$?

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
