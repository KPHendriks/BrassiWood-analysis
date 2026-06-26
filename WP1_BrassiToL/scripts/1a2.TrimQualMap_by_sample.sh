#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

sample="${1:?Usage: 1a2.TrimQualMap_by_sample.sh <sample_id>}"

###############################################################################
# Configuration
###############################################################################

# Main repository root (small/medium persistent files, scripts, final outputs)
REPO_ROOT="${REPO_ROOT:-/data/KasperHendriks/BrassiWood}"

# Big-data root (temporary heavy files, mappings, downloads, etc.)
BIGDATA_ROOT="${BIGDATA_ROOT:-/data/KasperHendriks/BrassiWood_big_data}"

# Main intermediate results directory inside repo
RESULTS_INT="${RESULTS_INT:-${REPO_ROOT}/WP1_BrassiToL/results_intermediate}"

# Persistent per-sample result directories
HBP_SAMPLE_DIR="${HBP_SAMPLE_DIR:-${RESULTS_INT}/results_hybpiper_bySample}"
HBPH_SAMPLE_DIR="${HBPH_SAMPLE_DIR:-${RESULTS_INT}/results_hybphaser_bySample}"

# Additional persistent output directories
TRIM_LOG_DIR="${TRIM_LOG_DIR:-${RESULTS_INT}/trimmomatic_logfiles}"
READ_COUNT_DIR="${READ_COUNT_DIR:-${RESULTS_INT}/per_sample_read_counts}"
HBP_STATS_DIR="${HBP_STATS_DIR:-${HBP_SAMPLE_DIR}/per_sample_stats}"
HBP_SEQLEN_DIR="${HBP_SEQLEN_DIR:-${HBP_SAMPLE_DIR}/per_sample_seq_lengths}"

CPU_LIMIT_PER_SAMPLE="${CPU_LIMIT_PER_SAMPLE:-4}"

USE_B764="${USE_B764:-yes}"
USE_A353="${USE_A353:-yes}"
USE_PLASTOME="${USE_PLASTOME:-yes}"
USE_MITOGENOME="${USE_MITOGENOME:-yes}"
USE_FUNCTIONAL_GENES="${USE_FUNCTIONAL_GENES:-no}"

TARGET_B764="${TARGET_B764:-${REPO_ROOT}/WP1_BrassiToL/data/ref-at_orf_minus_trailing_STOP_codons_exons_binned_to_by_genes.fasta}"
TARGET_A353="${TARGET_A353:-${REPO_ROOT}/WP1_BrassiToL/data/A353_NewTargets_refgenome_minus_trailing_STOP_codons.fasta}"
TARGET_PLASTOME="${TARGET_PLASTOME:-${REPO_ROOT}/WP1_BrassiToL/data/NikHay_chloro_bait_brassicaceae_new.fasta}"
TARGET_MITOGENOME="${TARGET_MITOGENOME:-${REPO_ROOT}/WP1_BrassiToL/data/mitogenome_reference.fasta}"
TARGET_FUNCTIONAL_GENES="${TARGET_FUNCTIONAL_GENES:-${REPO_ROOT}/WP1_BrassiToL/data/AHL15.fasta}"

SAMPLES_TO_HYBSEQ_CSV="${SAMPLES_TO_HYBSEQ_CSV:-${REPO_ROOT}/WP1_BrassiToL/data/samples_to_hybseq.csv}"
MERGE_HELPER_PY="${MERGE_HELPER_PY:-${REPO_ROOT}/WP1_BrassiToL/scripts/custom_scripts/get_hybseq_ids_for_sample.py}"

CONSENSUS_SCRIPT="${CONSENSUS_SCRIPT:-${REPO_ROOT}/WP1_BrassiToL/scripts/hybphaser_scripts_updated_for_BrassiWood_project/1_generate_consensus_sequences.sh}"
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-/usr/share/java/trimmomatic-0.39.jar}"
ADAPTERS_FA="${ADAPTERS_FA:-${REPO_ROOT}/WP1_BrassiToL/data/illumina_adapters_for_trimmomatic_normal_and_palindrome_mode.fasta}"

# Temporary per-sample workspace in big-data area
TMP_SAMPLE_DIR="${TMP_SAMPLE_DIR:-${BIGDATA_ROOT}/WP1_BrassiToL/data_temp/${sample}}"

###############################################################################
# Optional Google Drive fallback
###############################################################################

# Public version: private cloud/raw-data retrieval is disabled.
USE_GDRIVE_FALLBACK="${USE_GDRIVE_FALLBACK:-no}"
RCLONE_BIN=""
GDRIVE_REMOTE=""
GDRIVE_FOLDER_ID=""

###############################################################################
# Helpers
###############################################################################

log() {
  echo "[$(date '+%F %T')] [$sample] $*"
}

die() {
  echo "[ERROR] [$sample] $*" >&2
  exit 1
}

require_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || die "Required command not found: $c"
}

safe_mkdir() {
  mkdir -p "$1"
}

cleanup_sample_tmp() {
  if [[ -d "${TMP_SAMPLE_DIR}" ]]; then
    rm -rf "${TMP_SAMPLE_DIR}"
  fi
}

on_exit() {
  cleanup_sample_tmp
}
trap on_exit EXIT

consensus_path_for_target() {
  local target_name="$1"
  echo "${HBPH_SAMPLE_DIR}/${sample}_${target_name}_consensus.fasta"
}

consensus_done_for_target() {
  local target_name="$1"
  [[ -f "$(consensus_path_for_target "$target_name")" ]]
}

target_enabled() {
  local flag="$1"
  [[ "$flag" == "yes" ]]
}

any_enabled_targets() {
  target_enabled "${USE_B764}" && return 0
  target_enabled "${USE_A353}" && return 0
  target_enabled "${USE_PLASTOME}" && return 0
  target_enabled "${USE_MITOGENOME}" && return 0
  target_enabled "${USE_FUNCTIONAL_GENES}" && return 0
  return 1
}

all_enabled_targets_done() {
  local needed=0
  local done=0

  if target_enabled "${USE_B764}"; then
    ((needed+=1))
    consensus_done_for_target "B764" && ((done+=1))
  fi
  if target_enabled "${USE_A353}"; then
    ((needed+=1))
    consensus_done_for_target "A353" && ((done+=1))
  fi
  if target_enabled "${USE_PLASTOME}"; then
    ((needed+=1))
    consensus_done_for_target "plastome" && ((done+=1))
  fi
  if target_enabled "${USE_MITOGENOME}"; then
    ((needed+=1))
    consensus_done_for_target "mitogenome" && ((done+=1))
  fi
  if target_enabled "${USE_FUNCTIONAL_GENES}"; then
    ((needed+=1))
    consensus_done_for_target "functional-genes" && ((done+=1))
  fi

  [[ "$needed" -gt 0 && "$needed" -eq "$done" ]]
}

write_sample_namelist() {
  echo "$sample" > namelist_sample.txt
}

count_gz_reads() {
  local fq="$1"
  echo $(( $(zcat "$fq" | wc -l) / 4 ))
}

count_reads() {
  local fq="$1"
  echo $(( $(wc -l < "$fq") / 4 ))
}

combine_consensus_fastas() {
  local input_dir="$1"
  local output_file="$2"

  : > "$output_file"

  shopt -s nullglob
  local found=0
  local file
  for file in "${input_dir}"/*.fasta; do
    found=1
    while IFS= read -r line; do
      if [[ "$line" == ">"* ]]; then
        local gene="${line#*-}"
        echo ">$gene" >> "$output_file"
      else
        echo "$line" >> "$output_file"
      fi
    done < "$file"
  done
  shopt -u nullglob

  [[ "$found" -eq 1 ]] || die "No consensus fasta files found in ${input_dir}"
}

combine_retrieved_fna() {
  local output_file="$1"

  : > "$output_file"

  shopt -s nullglob
  local found=0
  local file
  for file in *.FNA; do
    found=1
    local sample_name="${file%.FNA}"
    sed "1s/^>.*/>${sample_name}/" "$file" >> "$output_file"
  done
  shopt -u nullglob

  [[ "$found" -eq 1 ]] || die "No *.FNA files found after hybpiper retrieve_sequences"
}

copy_raw_data_for_local_sample() {
  [[ -f "${SAMPLES_TO_HYBSEQ_CSV}" ]] || die "samples_to_hybseq.csv not found: ${SAMPLES_TO_HYBSEQ_CSV}"

  python3 - <<'PY' "${SAMPLES_TO_HYBSEQ_CSV}" "${sample}"
import csv, sys
csv_file, sample = sys.argv[1], sys.argv[2]
hyb_ids = []
with open(csv_file, newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        if row["Sample_ID"] == sample:
            hyb_ids.append(row["HybSeq_ID"])
if not hyb_ids:
    raise SystemExit(f"[ERROR] No HybSeq_ID entries found for sample {sample}")
print("\n".join(hyb_ids))
PY
}

minio_raw_object_r1() {
  die "Private MinIO access is disabled in the public version."
}

minio_raw_object_r2() {
  die "Private MinIO access is disabled in the public version."
}

download_one_from_gdrive() {
  die "Private Google Drive fallback is disabled in the public version."
}

download_hybid_pair_from_minio() {
  die "Private MinIO access is disabled in the public version."
}

download_hybid_pair_from_gdrive() {
  die "Private Google Drive fallback is disabled in the public version."
}

fetch_raw_data() {
  log "Checking raw data locally"

  local local_r1="${RAW_SRA_DIR}/${sample}_R1.fastq.gz"
  local local_r2="${RAW_SRA_DIR}/${sample}_R2.fastq.gz"

  if [[ "$sample" == SRR* || "$sample" == ERR* ]]; then
    require_cmd gzip
    require_cmd rm

    log "Downloading public SRA/ENA reads for ${sample}"

    ~/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump --split-files -f -O ./ "$sample"
    gzip -c "${sample}_1.fastq" > "${sample}_R1.fastq.gz"
    gzip -c "${sample}_2.fastq" > "${sample}_R2.fastq.gz"
    rm -f "${sample}_1.fastq" "${sample}_2.fastq"
    rm -f "${BIGDATA_ROOT}/raw_data_sra/sra/${sample}"* || true

  elif [[ -f "$local_r1" && -f "$local_r2" ]]; then
    log "Using local FASTQ files from ${RAW_SRA_DIR}"
    cp "$local_r1" "./${sample}_R1.fastq.gz"
    cp "$local_r2" "./${sample}_R2.fastq.gz"

  else
    die "Raw FASTQ files not found for ${sample}. Public version does not fetch private data. Expected:
  ${local_r1}
  ${local_r2}"
  fi

  log "Finished preparing raw data"
}

run_trimmomatic() {
  log "Start trimming with Trimmomatic"

  require_cmd java
  require_cmd gunzip
  require_cmd zcat

  [[ -f "${TRIMMOMATIC_JAR}" ]] || die "Trimmomatic jar not found: ${TRIMMOMATIC_JAR}"
  [[ -f "${ADAPTERS_FA}" ]] || die "Adapter fasta not found: ${ADAPTERS_FA}"

  raw_reads_count="$(count_gz_reads "${sample}_R1.fastq.gz")"

  java -jar "${TRIMMOMATIC_JAR}" PE -phred33 \
    -threads "${CPU_LIMIT_PER_SAMPLE}" \
    "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" \
    "${sample}_R1_paired.fastq.gz" "${sample}_R1_unpaired.fastq.gz" \
    "${sample}_R2_paired.fastq.gz" "${sample}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS_FA}":2:30:10:2:true \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:40 \
    2> "${TRIM_LOG_DIR}/log_trim_${sample}.txt"

  cat "${sample}"_R*_unpaired.fastq.gz > "${sample}_unpaired.fastq.gz"

  gunzip -c "${sample}_R1_paired.fastq.gz" > "${sample}_R1_paired.fastq"
  gunzip -c "${sample}_R2_paired.fastq.gz" > "${sample}_R2_paired.fastq"
  gunzip -c "${sample}_unpaired.fastq.gz" > "${sample}_unpaired.fastq"

  paired_reads_count="$(count_reads "${sample}_R1_paired.fastq")"
  unpaired_reads_count="$(count_reads "${sample}_unpaired.fastq")"

  printf "%s\t%s\t%s\t%s\n" \
    "$sample" "$raw_reads_count" "$paired_reads_count" "$unpaired_reads_count" \
    > "${READ_COUNT_DIR}/${sample}.tsv"

  rm -f *.fastq.gz
  log "Finished trimming with Trimmomatic"
}

backup_if_possible() {
  # Public version: remote backups disabled.
  return 0
}

run_target_mapping() {
  local target_name="$1"
  local target_fasta="$2"

  if consensus_done_for_target "$target_name"; then
    log "Skipping ${target_name}; consensus fasta already exists at $(consensus_path_for_target "$target_name")"
    return 0
  fi

  [[ -f "${target_fasta}" ]] || die "Target fasta not found for ${target_name}: ${target_fasta}"

  log "Start mapping against ${target_name}"

  mkdir -p "${target_name}"
  cd "${target_name}"

  mv ../*.fastq .

  hybpiper assemble \
    --readfiles "${sample}"_R*_paired.fastq \
    --unpaired "${sample}_unpaired.fastq" \
    --targetfile_dna "${target_fasta}" \
    --prefix "$sample" \
    --bwa \
    --no_intronerate \
    --cpu "${CPU_LIMIT_PER_SAMPLE}"

  bash "${CONSENSUS_SCRIPT}" \
    -s "$sample" \
    -p . \
    -o "${sample}_hybphaser_output" \
    -t "${CPU_LIMIT_PER_SAMPLE}"

  hybpiper stats \
    --targetfile_dna "${target_fasta}" \
    gene \
    ../namelist_sample.txt

  cp hybpiper_stats.tsv \
    "${HBP_STATS_DIR}/${target_name}__${sample}.tsv"

  cp seq_lengths.tsv \
    "${HBP_SEQLEN_DIR}/${target_name}__${sample}.tsv"

  hybpiper retrieve_sequences \
    --targetfile_dna "${target_fasta}" \
    dna \
    --sample_names ../namelist_sample.txt

  combine_retrieved_fna "${sample}_${target_name}.fasta"

  if [[ -f "${sample}/${sample}_genes_with_long_paralog_warnings.txt" ]]; then
    backup_if_possible \
      "${sample}/${sample}_genes_with_long_paralog_warnings.txt" \
      "*/brassiwood/WP1_BrassiToL/results_intermediate/results_hybpiper_bySample/${sample}_${target_name}_genes_with_long_paralog_warnings.txt"

    mv "${sample}/${sample}_genes_with_long_paralog_warnings.txt" \
      "${HBP_SAMPLE_DIR}/${sample}_${target_name}_genes_with_long_paralog_warnings.txt"
  fi

  backup_if_possible \
    "${sample}_${target_name}.fasta" \
    "*/brassiwood/WP1_BrassiToL/results_intermediate/results_hybpiper_bySample/${sample}_${target_name}.fasta"

  mv "${sample}_${target_name}.fasta" \
    "${HBP_SAMPLE_DIR}/"

  combine_consensus_fastas \
    "${sample}_hybphaser_output/01_data/${sample}/consensus" \
    "${sample}_${target_name}_consensus.fasta"

  backup_if_possible \
    "${sample}_${target_name}_consensus.fasta" \
    "*/brassiwood/WP1_BrassiToL/results_intermediate/results_hybphaser_bySample/${sample}_${target_name}_consensus.fasta"

  mv "${sample}_${target_name}_consensus.fasta" \
    "${HBPH_SAMPLE_DIR}/"

  mv *.fastq ../
  cd ../
  rm -rf "${target_name}"

  log "Finished mapping against ${target_name}"
}

###############################################################################
# Main
###############################################################################

require_cmd hybpiper
require_cmd python3
require_cmd bash

any_enabled_targets || die "No targets enabled; set at least one USE_* flag to yes"

safe_mkdir "${TMP_SAMPLE_DIR}"
safe_mkdir "${HBP_SAMPLE_DIR}"
safe_mkdir "${HBPH_SAMPLE_DIR}"
safe_mkdir "${TRIM_LOG_DIR}"
safe_mkdir "${READ_COUNT_DIR}"
safe_mkdir "${HBP_STATS_DIR}"
safe_mkdir "${HBP_SEQLEN_DIR}"

cd "${TMP_SAMPLE_DIR}"
write_sample_namelist

log "Checking existing consensus files in ${HBPH_SAMPLE_DIR}"

if target_enabled "${USE_B764}"; then
  if consensus_done_for_target "B764"; then
    log "Existing found: $(consensus_path_for_target "B764")"
  else
    log "Missing: $(consensus_path_for_target "B764")"
  fi
fi

if target_enabled "${USE_A353}"; then
  if consensus_done_for_target "A353"; then
    log "Existing found: $(consensus_path_for_target "A353")"
  else
    log "Missing: $(consensus_path_for_target "A353")"
  fi
fi

if target_enabled "${USE_PLASTOME}"; then
  if consensus_done_for_target "plastome"; then
    log "Existing found: $(consensus_path_for_target "plastome")"
  else
    log "Missing: $(consensus_path_for_target "plastome")"
  fi
fi

if target_enabled "${USE_MITOGENOME}"; then
  if consensus_done_for_target "mitogenome"; then
    log "Existing found: $(consensus_path_for_target "mitogenome")"
  else
    log "Missing: $(consensus_path_for_target "mitogenome")"
  fi
fi

if target_enabled "${USE_FUNCTIONAL_GENES}"; then
  if consensus_done_for_target "functional-genes"; then
    log "Existing found: $(consensus_path_for_target "functional-genes")"
  else
    log "Missing: $(consensus_path_for_target "functional-genes")"
  fi
fi

if all_enabled_targets_done; then
  log "Skipping sample; all enabled targets already completed"
  exit 0
fi

fetch_raw_data
run_trimmomatic

if target_enabled "${USE_B764}"; then
  run_target_mapping "B764" "${TARGET_B764}"
fi

if target_enabled "${USE_A353}"; then
  run_target_mapping "A353" "${TARGET_A353}"
fi

if target_enabled "${USE_PLASTOME}"; then
  run_target_mapping "plastome" "${TARGET_PLASTOME}"
fi

if target_enabled "${USE_MITOGENOME}"; then
  run_target_mapping "mitogenome" "${TARGET_MITOGENOME}"
fi

if target_enabled "${USE_FUNCTIONAL_GENES}"; then
  run_target_mapping "functional-genes" "${TARGET_FUNCTIONAL_GENES}"
fi

log "All done"
