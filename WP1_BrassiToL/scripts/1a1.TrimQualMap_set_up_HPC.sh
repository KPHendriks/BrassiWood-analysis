#!/usr/bin/env bash
#SBATCH -J mapRaw
#SBATCH -o mapRaw.o%j
#SBATCH -e mapRaw.e%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=31

set -euo pipefail
IFS=$'\n\t'

###############################################################################
# STEP 0 — project settings
###############################################################################

PROJECT_NAME="2026-06-16_BrassiWood_study_v1"
SAMPLES_FILE="WP1_BrassiToL/data/SampleListStart.txt"
SUBSET_NAME="Cardamine_baldensis_study"

###############################################################################
# STEP 0a — target settings
###############################################################################

USE_B764=yes
USE_A353=yes
USE_PLASTOME=no
USE_MITOGENOME=no
USE_FUNCTIONAL_GENES=no

###############################################################################
# STEP 0b — parallelism and CPU settings
###############################################################################

# Number of samples to process simultaneously with GNU parallel
PARALLEL_JOBS=20

# Number of CPUs each per-sample job may use internally
CPU_LIMIT_PER_SAMPLE=2

# Safety check: keep this <= allocated CPUs unless you intentionally oversubscribe
TOTAL_CPU_REQUEST=$((PARALLEL_JOBS * CPU_LIMIT_PER_SAMPLE))
SLURM_CPU_BUDGET="${SLURM_CPUS_PER_TASK:-20}"

if (( TOTAL_CPU_REQUEST > SLURM_CPU_BUDGET )); then
  echo "[WARNING] PARALLEL_JOBS * CPU_LIMIT_PER_SAMPLE = ${TOTAL_CPU_REQUEST}, which exceeds SLURM_CPUS_PER_TASK = ${SLURM_CPU_BUDGET}" >&2
  echo "[WARNING] This may oversubscribe CPUs. Consider lowering PARALLEL_JOBS or CPU_LIMIT_PER_SAMPLE." >&2
fi

###############################################################################
# STEP 0c — core paths
###############################################################################

REPO_ROOT="."
BIGDATA_ROOT="../BrassiWood_big_data"
RAW_SRA_DIR="${BIGDATA_ROOT}/raw_data_sra"

RESULTS_INT="${REPO_ROOT}/WP1_BrassiToL/results_intermediate"
HBP_SAMPLE_DIR="${RESULTS_INT}/results_hybpiper_bySample"
HBPH_SAMPLE_DIR="${RESULTS_INT}/results_hybphaser_bySample"
SNP_ASSESS_DIR="${RESULTS_INT}/results_hybphaser_SNP_assessment/${PROJECT_NAME}"

###############################################################################
# STEP 0d — reference target files
###############################################################################

TARGET_B764="${REPO_ROOT}/WP1_BrassiToL/data/ref-at_orf_minus_trailing_STOP_codons_exons_binned_to_by_genes.fasta"
TARGET_A353="${REPO_ROOT}/WP1_BrassiToL/data/A353_NewTargets_refgenome_minus_trailing_STOP_codons.fasta"
TARGET_PLASTOME="${REPO_ROOT}/WP1_BrassiToL/data/NikHay_chloro_bait_brassicaceae_new.fasta"
TARGET_MITOGENOME="${REPO_ROOT}/WP1_BrassiToL/data/mitogenome_reference.fasta"
TARGET_FUNCTIONAL_GENES="${REPO_ROOT}/WP1_BrassiToL/data/AHL15.fasta"

###############################################################################
# STEP 0e — called scripts
###############################################################################

SCRIPT_TRIM_MAP_BY_SAMPLE="${REPO_ROOT}/WP1_BrassiToL/scripts/1a2.TrimQualMap_by_sample.sh"

SCRIPT_COUNT_SNPS="${REPO_ROOT}/WP1_BrassiToL/scripts/hybphaser_scripts_updated_for_BrassiWood_project/1a_count_snps.py"
SCRIPT_ASSESS_DATASET="${REPO_ROOT}/WP1_BrassiToL/scripts/hybphaser_scripts_updated_for_BrassiWood_project/1b_assess_dataset.py"
SCRIPT_GENERATE_SEQUENCE_LISTS="${REPO_ROOT}/WP1_BrassiToL/scripts/hybphaser_scripts_updated_for_BrassiWood_project/1c_generate_sequence_lists.py"

SCRIPT_CLEAN_STATS="${REPO_ROOT}/WP1_BrassiToL/scripts/custom_scripts/clean_stats_tsv.py"
SCRIPT_DEDUP_A353_B764="${REPO_ROOT}/WP1_BrassiToL/scripts/custom_scripts/deduplicate_A353_B764_overlap.py"

###############################################################################
# STEP 0f — private raw-data access disabled in public version
###############################################################################

# Raw sequencing data are not fetched by the public version of this script.
# Place required FASTQ files manually in:
#   ../BrassiWood_big_data/raw_data_sra
# Or, collect raw sample data via NCBI SRA accession codes instead of sample numbers;
# this script will do so automatically when providing any of SRR*, ERR*, or DRR* numbers.
# See the publication's sample sheet for all relevant accession codes.
#
# The original private workflow could optionally retrieve missing reads from
# institutional/cloud storage, but all remote identifiers and credentials have
# been removed from the public repository.

export USE_GDRIVE_FALLBACK="no"
export RCLONE_BIN=""
export GDRIVE_REMOTE=""
export GDRIVE_FOLDER_ID=""

###############################################################################
# STEP 0g — helpers
###############################################################################

log() {
  echo "[$(date '+%F %T')] [MASTER] $*"
}

die() {
  echo "[ERROR] [MASTER] $*" >&2
  exit 1
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || die "Required file not found: $f"
}

require_dir() {
  local d="$1"
  [[ -d "$d" ]] || die "Required directory not found: $d"
}

require_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || die "Required command not found in PATH: $c"
}

safe_mkdir() {
  mkdir -p "$1"
}

###############################################################################
# STEP 0h — sanity checks
###############################################################################

require_dir "${REPO_ROOT}"
require_file "${SAMPLES_FILE}"

require_cmd parallel
require_cmd python3
require_cmd bash

require_file "${SCRIPT_TRIM_MAP_BY_SAMPLE}"
require_file "${SCRIPT_COUNT_SNPS}"
require_file "${SCRIPT_ASSESS_DATASET}"
require_file "${SCRIPT_GENERATE_SEQUENCE_LISTS}"
require_file "${SCRIPT_CLEAN_STATS}"
require_file "${SCRIPT_DEDUP_A353_B764}"

if [[ "${USE_B764}" == "yes" ]]; then
  require_file "${TARGET_B764}"
fi

if [[ "${USE_A353}" == "yes" ]]; then
  require_file "${TARGET_A353}"
fi

if [[ "${USE_PLASTOME}" == "yes" ]]; then
  require_file "${TARGET_PLASTOME}"
fi

if [[ "${USE_MITOGENOME}" == "yes" ]]; then
  require_file "${TARGET_MITOGENOME}"
fi

if [[ "${USE_FUNCTIONAL_GENES}" == "yes" ]]; then
  require_file "${TARGET_FUNCTIONAL_GENES}"
fi

if [[ "${USE_B764}" != "yes" && "${USE_A353}" != "yes" && "${USE_PLASTOME}" != "yes" && "${USE_MITOGENOME}" != "yes" && "${USE_FUNCTIONAL_GENES}" != "yes" ]]; then
  die "No targets enabled. Set at least one USE_* variable to yes."
fi

###############################################################################
# STEP 0i — create directories
###############################################################################

safe_mkdir "${BIGDATA_ROOT}"
safe_mkdir "${RAW_SRA_DIR}"

safe_mkdir "${BIGDATA_ROOT}/data_mapped_B764"
safe_mkdir "${BIGDATA_ROOT}/data_mapped_A353"
safe_mkdir "${BIGDATA_ROOT}/WP1_BrassiToL"
safe_mkdir "${BIGDATA_ROOT}/WP1_BrassiToL/data_temp"

safe_mkdir "${RESULTS_INT}"
safe_mkdir "${RESULTS_INT}/trimmomatic_logfiles"
safe_mkdir "${RESULTS_INT}/per_sample_read_counts"

safe_mkdir "${HBP_SAMPLE_DIR}"
safe_mkdir "${HBPH_SAMPLE_DIR}"
safe_mkdir "${HBP_SAMPLE_DIR}/per_sample_stats"
safe_mkdir "${HBP_SAMPLE_DIR}/per_sample_seq_lengths"

safe_mkdir "${RESULTS_INT}/results_hybpiper_byGene"
safe_mkdir "${RESULTS_INT}/results_hybphaser_SNP_assessment"
safe_mkdir "${SNP_ASSESS_DIR}"

###############################################################################
# STEP 0j — build target arguments for downstream Python scripts
###############################################################################

TARGET_ARGS=()

if [[ "${USE_B764}" == "yes" ]]; then
  TARGET_ARGS+=(--targets "B764=${TARGET_B764}")
fi

if [[ "${USE_A353}" == "yes" ]]; then
  TARGET_ARGS+=(--targets "A353=${TARGET_A353}")
fi

if [[ "${USE_PLASTOME}" == "yes" ]]; then
  TARGET_ARGS+=(--targets "plastome=${TARGET_PLASTOME}")
fi

if [[ "${USE_MITOGENOME}" == "yes" ]]; then
  TARGET_ARGS+=(--targets "mitogenome=${TARGET_MITOGENOME}")
fi

if [[ "${USE_FUNCTIONAL_GENES}" == "yes" ]]; then
  TARGET_ARGS+=(--targets "functional-genes=${TARGET_FUNCTIONAL_GENES}")
fi

###############################################################################
# STEP 0k — activate conda environment
###############################################################################

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate hybpiper

command -v hybpiper >/dev/null 2>&1 || {
  echo "[ERROR] hybpiper not found after activating conda environment" >&2
  exit 1
}

###############################################################################
# STEP 0l — export settings for per-sample jobs
###############################################################################

export PROJECT_NAME
export SUBSET_NAME

export REPO_ROOT
export BIGDATA_ROOT
export RAW_SRA_DIR
export RESULTS_INT
export HBP_SAMPLE_DIR
export HBPH_SAMPLE_DIR
export SNP_ASSESS_DIR

export USE_B764
export USE_A353
export USE_PLASTOME
export USE_MITOGENOME
export USE_FUNCTIONAL_GENES

export TARGET_B764
export TARGET_A353
export TARGET_PLASTOME
export TARGET_MITOGENOME
export TARGET_FUNCTIONAL_GENES

export CPU_LIMIT_PER_SAMPLE

###############################################################################
# STEP 0m — run summary
###############################################################################

log "Starting mapping pipeline"
log "Project name: ${PROJECT_NAME}"
log "Subset name: ${SUBSET_NAME}"
log "Sample list: ${SAMPLES_FILE}"
log "Repository root: ${REPO_ROOT}"
log "Big-data root: ${BIGDATA_ROOT}"
log "Results dir: ${RESULTS_INT}"
log "Parallel jobs: ${PARALLEL_JOBS}"
log "CPU per sample: ${CPU_LIMIT_PER_SAMPLE}"
log "Total requested by parallel layer: ${TOTAL_CPU_REQUEST}"
log "SLURM cpus-per-task: ${SLURM_CPU_BUDGET}"
log "Enabled targets: B764=${USE_B764}, A353=${USE_A353}, plastome=${USE_PLASTOME}, mitogenome=${USE_MITOGENOME}, functional-genes=${USE_FUNCTIONAL_GENES}"

###############################################################################
# STEP 1 — run trimming + mapping per sample
###############################################################################

log "STEP 1: starting per-sample trimming and mapping with GNU parallel"

parallel -j "${PARALLEL_JOBS}" \
  bash "${SCRIPT_TRIM_MAP_BY_SAMPLE}" {} \
  :::: "${SAMPLES_FILE}"

log "STEP 1: per-sample trimming and mapping finished"

###############################################################################
# STEP 2 — aggregate per-sample outputs and clean summary tables
###############################################################################

log "STEP 2: aggregating per-sample outputs"

# Read counts
if compgen -G "${RESULTS_INT}/per_sample_read_counts/*.tsv" > /dev/null; then
  cat "${RESULTS_INT}"/per_sample_read_counts/*.tsv > \
    "${RESULTS_INT}/Trimmomatic_read_counts.txt"
fi

aggregate_target_tables() {
  local target_name="$1"

  # hybpiper_stats
  if compgen -G "${HBP_SAMPLE_DIR}/per_sample_stats/${target_name}"__*.tsv > /dev/null; then
    cat "${HBP_SAMPLE_DIR}"/per_sample_stats/"${target_name}"__*.tsv > \
      "${HBP_SAMPLE_DIR}/hybpiper_stats_${target_name}.tsv"

    python3 "${SCRIPT_CLEAN_STATS}" \
      --stats_input_file "${HBP_SAMPLE_DIR}/hybpiper_stats_${target_name}.tsv" \
      --stats_output_file "${HBP_SAMPLE_DIR}/hybpiper_stats_${target_name}_cleaned.tsv"
  fi

  # seq_lengths
  if compgen -G "${HBP_SAMPLE_DIR}/per_sample_seq_lengths/${target_name}"__*.tsv > /dev/null; then
    cat "${HBP_SAMPLE_DIR}"/per_sample_seq_lengths/"${target_name}"__*.tsv > \
      "${HBP_SAMPLE_DIR}/hybpiper_seq_lengths_${target_name}.tsv"
  fi
}

[[ "${USE_B764}" == "yes" ]] && aggregate_target_tables "B764"
[[ "${USE_A353}" == "yes" ]] && aggregate_target_tables "A353"
[[ "${USE_PLASTOME}" == "yes" ]] && aggregate_target_tables "plastome"
[[ "${USE_MITOGENOME}" == "yes" ]] && aggregate_target_tables "mitogenome"
[[ "${USE_FUNCTIONAL_GENES}" == "yes" ]] && aggregate_target_tables "functional-genes"

log "STEP 2: aggregation finished"

###############################################################################
# STEP 3 — HybPhaser SNP assessment
###############################################################################

log "STEP 3: running HybPhaser SNP assessment"

####################################################################################################
# STEP 3a — Count SNP proportions per locus × sample
####################################################################################################

python3 "${SCRIPT_COUNT_SNPS}" \
  --consensus_dir "${HBPH_SAMPLE_DIR}" \
  --samples_file "${SAMPLES_FILE}" \
  "${TARGET_ARGS[@]}" \
  --targets_format DNA \
  --target_id_mode dash \
  --allow_extra_loci no \
  --out_root "${RESULTS_INT}/results_hybphaser_SNP_assessment/${PROJECT_NAME}" \
  --subset_name "${SUBSET_NAME}"

####################################################################################################
# STEP 3b — Assess dataset
####################################################################################################

python3 "${SCRIPT_ASSESS_DATASET}" \
  --out_root "${RESULTS_INT}/results_hybphaser_SNP_assessment/${PROJECT_NAME}" \
  --subset_name "${SUBSET_NAME}" \
  "${TARGET_ARGS[@]}" \
  --targets_format DNA \
  --target_id_mode dash \
  --thr_samples_prop_loci 0.01 \
  --thr_samples_prop_target 0.01 \
  --thr_loci_prop_samples 0.2 \
  --thr_loci_prop_target 0.2 \
  --min_len_frac_cell 0.25 \
  --paralogs_global none \
  --paralogs_each yes \
  --include_hybpiper yes \
  --hybpiper_dir "${HBP_SAMPLE_DIR}" \
  --plots_pdf \
  --save_args_txt

####################################################################################################
# STEP 3c — Generate sequence lists
####################################################################################################

python3 "${SCRIPT_GENERATE_SEQUENCE_LISTS}" \
  --out_root "${RESULTS_INT}/results_hybphaser_SNP_assessment/${PROJECT_NAME}" \
  --subset_name "${SUBSET_NAME}" \
  --samples_file "${SAMPLES_FILE}" \
  --layout flat \
  --consensus_dir "${HBPH_SAMPLE_DIR}" \
  --contig_dir "${HBP_SAMPLE_DIR}" \
  --intronerated no \
  --remove_hybpiper yes \
  "${TARGET_ARGS[@]}" \
  --targets_format DNA \
  --target_id_mode dash \
  --log_every 250

####################################################################################################
# STEP 3d — Deduplicate loci present in both A353 and B764 bait sets
####################################################################################################

if [[ "${USE_B764}" == "yes" && "${USE_A353}" == "yes" ]]; then
  log "Both B764 and A353 are enabled: deduplicating overlapping loci"

  python3 "${SCRIPT_DEDUP_A353_B764}" \
    --seq_root "${RESULTS_INT}/results_hybphaser_SNP_assessment/${PROJECT_NAME}/03_sequence_lists_${SUBSET_NAME}" \
    --prefer_on_tie B764 \
    --min_diff_n_consensus 50 \
    --min_diff_n_contig 50
else
  log "Skipping A353/B764 overlap deduplication (both bait sets not jointly enabled)"
fi

###############################################################################
# Final note
###############################################################################

log "Pipeline finished successfully"
