#!/bin/bash
# Create MSAs (per gene) + trimAl + IQ-TREE scaffolding on MaaS36
# Runs directly (no SLURM). Uses project-specific HybPhaser step-1c output from 1a1.

# Public version: this script does not download private raw data.
# It expects HybPhaser step-1c consensus FASTA files to already exist under:
#   WP1_BrassiToL/results_intermediate/results_hybphaser_SNP_assessment/<PROJECT_NAME>/03_sequence_lists_<SUBSET_NAME>/loci_consensus

set -euo pipefail
IFS=$'\n\t'
trap 'echo "[ERR] at line $LINENO"' ERR

# ----------- CONFIG / ROOTS -----------
ROOT="."
BIG="${ROOT%/*}/BrassiWood_big_data"
RINT="$ROOT/WP1_BrassiToL/results_intermediate"
RFIN="$ROOT/WP1_BrassiToL/results_final"
DATA="$ROOT/WP1_BrassiToL/data"
SCRIPTS="$ROOT/WP1_BrassiToL/scripts"

# ----------- PROJECT SETTINGS -----------
PROJECT_NAME="2026-06-16_BrassiWood_study_v1"
SAMPLES_FILE="WP1_BrassiToL/data/SampleListStart.txt"
SUBSET_NAME="Cardamine_baldensis_study"

USE_B764=yes
USE_A353=yes
USE_PLASTOME=no
USE_MITOGENOME=no

# Reference target files
TARGET_B764="$DATA/ref-at_orf_minus_trailing_STOP_codons_exons_binned_to_by_genes.fasta"
TARGET_A353="$DATA/A353_NewTargets_refgenome_minus_trailing_STOP_codons.fasta"


# Optional
SAMPLE_OUTGROUP="S2362"
INHERIT_PROJECT="none"   # or set to an older project_name to reuse finished outputs

# ----------- ALIGNMENT / TREE SETTINGS -----------
METHOD_ALIGNMENT="macse"   # macse only here
CPU_LIMIT=2
MAX_PARALLEL_GENES=20
MACSE_MEM_GB=15
IQTREE_MEM_GB=6
IQTREE_THREADS=4

# Executables
TRIMAL_BIN="${TRIMAL_BIN:-$HOME/trimal/source/trimal}"
IQTREE_BIN="${IQTREE_BIN:-$HOME/iqtree-3.0.1-Linux/bin/iqtree3}"

# ----------- SANITY CHECKS -----------
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1"; exit 127; }; }
need parallel
need singularity
need python3

[[ -f "$SAMPLES_FILE" ]] || { echo "[ERR] Samples file not found: $SAMPLES_FILE"; exit 1; }

[[ -f "$TARGET_B764"       || "$USE_B764"       != "yes" ]] || { echo "[ERR] Missing TARGET_B764: $TARGET_B764"; exit 1; }
[[ -f "$TARGET_A353"       || "$USE_A353"       != "yes" ]] || { echo "[ERR] Missing TARGET_A353: $TARGET_A353"; exit 1; }

# ----------- EXPORT SHARED VARS FOR 1b2 -----------
export project_name="$PROJECT_NAME"
export sample_list="$(basename "$SAMPLES_FILE")"
export sample_outgroup="$SAMPLE_OUTGROUP"
export inherit_project="$INHERIT_PROJECT"
export method_alignment="$METHOD_ALIGNMENT"
export cpu_limit="$CPU_LIMIT"

export MACSE_MEM_GB
export IQTREE_MEM_GB
export IQTREE_THREADS
export MAX_PARALLEL_GENES
export TRIMAL_BIN
export IQTREE_BIN

# ----------- WHICH GENOMES TO RUN -----------
genome_list=()
[[ "$USE_B764" == "yes" ]]      && genome_list+=(B764)
[[ "$USE_A353" == "yes" ]]      && genome_list+=(A353)

if (( ${#genome_list[@]} == 0 )); then
  echo "[ERR] No enabled genomes selected."
  exit 1
fi

# ----------- PROJECT-SPECIFIC INPUT LOCI DIR -----------
CONS_LOCI_DIR="$RINT/results_hybphaser_SNP_assessment/$PROJECT_NAME/03_sequence_lists_${SUBSET_NAME}/loci_consensus"

if [[ ! -d "$CONS_LOCI_DIR" ]]; then
  echo "[ERR] Project-specific loci directory not found:"
  echo "      $CONS_LOCI_DIR"
  exit 1
fi

# ----------- LOGGING -----------
RUNLOG="$RFIN/$PROJECT_NAME/0_stats/genes_1b2_runlog.csv"
export RUNLOG
mkdir -p "$(dirname "$RUNLOG")"
mkdir -p "$RFIN/$PROJECT_NAME/0_stats/gene_lists"
mkdir -p "$RFIN/$PROJECT_NAME/0_stats/gene_logs"

# ----------- DIRECTORIES -----------
mkdir -p "$BIG" "$BIG/WP1_BrassiToL/scripts/macse_pipeline"
mkdir -p "$RFIN/$PROJECT_NAME/0_stats"

for genome in "${genome_list[@]}"; do
  mkdir -p "$RFIN/$PROJECT_NAME/1_results_gene_alignments_${METHOD_ALIGNMENT}/$genome"
  mkdir -p "$RFIN/$PROJECT_NAME/2_results_gene_alignments_${METHOD_ALIGNMENT}_trimal/$genome"
  mkdir -p "$RFIN/$PROJECT_NAME/3_results_iqtree_gene_trees/$genome"
  mkdir -p "$RFIN/$PROJECT_NAME/0_stats/gene_logs/$genome"
done

# ----------- DEFAULT GENE LISTS (FROM TARGET FASTAS) -----------
mapfile -t gene_list_B764 < <(
  grep -a '^>' "$TARGET_B764" \
  | sed 's/^>//' \
  | cut -d'-' -f2- \
  | tr -d '\r' \
  | sort -u
)

mapfile -t gene_list_A353 < <(
  grep -a '^>' "$TARGET_A353" \
  | sed 's/^>//' \
  | cut -d'-' -f2- \
  | tr -d '\r' \
  | sort -u
)


# ----------- WRITE AND CHECK GENE LISTS BEFORE RUNNING ALIGNMENTS -----------

for genome in "${genome_list[@]}"; do
  declare -n arr="gene_list_${genome}"

  parsed_file="$RFIN/$PROJECT_NAME/0_stats/gene_lists/genes_parsed_from_target_${genome}.txt"
  torun_file="$RFIN/$PROJECT_NAME/0_stats/gene_lists/genes_to_run_${genome}.txt"

  printf '%s\n' "${arr[@]}" > "$parsed_file"

  mapfile -t filtered < <(
    for g in "${arr[@]}"; do
      if [[ -s "$CONS_LOCI_DIR/${g}_consensus.fasta" || -s "$CONS_LOCI_DIR/${g}_intronerated_consensus.fasta" ]]; then
        printf '%s\n' "$g"
      fi
    done
  )

  printf '%s\n' "${filtered[@]}" > "$torun_file"

  echo "[INFO] $genome: ${#arr[@]} genes parsed from target FASTA"
  echo "[INFO] $genome: ${#filtered[@]} genes found in loci_consensus"

  if (( ${#arr[@]} == 0 )); then
    echo "[WARN] $genome: no genes parsed from target FASTA"
  fi

  if (( ${#filtered[@]} == 0 )); then
    echo "[WARN] $genome: no matching consensus FASTAs found in:"
    echo "       $CONS_LOCI_DIR"
  fi
done

echo "[INFO] Gene-list files written:"
ls -lh "$RFIN/$PROJECT_NAME/0_stats/gene_lists"


# ----------- OMM-MACSE CONTAINER -----------
MACSE_DIR="$BIG/WP1_BrassiToL/scripts/macse_pipeline"
if [[ ! -s "$MACSE_DIR/omm_macse_v10.02.sif" ]]; then
  singularity pull --arch amd64 --dir "$MACSE_DIR" library://vranwez/default/omm_macse:v10.02
fi

# ----------- RUN 1b2 PER GENOME (PARALLEL PER GENE) -----------
export ROOT BIG RINT RFIN DATA SCRIPTS MACSE_DIR CONS_LOCI_DIR
export project_name sample_list sample_outgroup inherit_project method_alignment cpu_limit
export MACSE_MEM_GB IQTREE_MEM_GB IQTREE_THREADS MAX_PARALLEL_GENES

for genome in "${genome_list[@]}"; do
  echo "[INFO] Start alignment for genome: $genome"
  export genome

  torun_file="$RFIN/$PROJECT_NAME/0_stats/gene_lists/genes_to_run_${genome}.txt"

  if [[ ! -s "$torun_file" ]]; then
    echo "[WARN] No surviving loci for $genome — skipping."
    continue
  fi

  mapfile -t filtered < "$torun_file"

  echo "[INFO] $genome: ${#filtered[@]} genes queued for 1b2"

  parallel --lb \
           --joblog "$RFIN/$PROJECT_NAME/0_stats/parallel_${genome}.log" \
           -j "$MAX_PARALLEL_GENES" \
           --env project_name,genome,inherit_project,method_alignment,cpu_limit,ROOT,BIG,RINT,RFIN,DATA,SCRIPTS,MACSE_DIR,CONS_LOCI_DIR,MACSE_MEM_GB,IQTREE_MEM_GB,IQTREE_THREADS,TRIMAL_BIN,IQTREE_BIN,RUNLOG \
           bash "$SCRIPTS/1b2.Create_MSAs_by_gene.sh" ::: "${filtered[@]}" || true
done

echo "[OK] Finished alignment + IQ-TREE jobs."

# ---------- SUMMARIZE IQ-TREE SUPPORTS & BRANCH LENGTHS ----------
python3 "$SCRIPTS/custom_scripts/summarize_iqtree_supports.py" \
  --input-glob "$RFIN/$PROJECT_NAME/3_results_iqtree_gene_trees/*/*_iqtree.treefile" \
  --outdir "$RFIN/$PROJECT_NAME/0_stats" \
  --alrt-threshold 70 \
  --ufboot-threshold 80 \
  --bl-threshold 1e-5 \
  --log10-bl-hist

# ---------- SIMPLE FAILURE SUMMARY ----------
if [[ -s "$RUNLOG" ]]; then
  awk -F, 'NR==1 || $10 != 0' "$RUNLOG" > "$RFIN/$PROJECT_NAME/0_stats/genes_1b2_runlog_failures_only.csv" || true
  awk -F, 'NR>1 && $10 != 0 {print $3}' "$RUNLOG" | sort -u > "$RFIN/$PROJECT_NAME/0_stats/failed_genes_unique.txt" || true
fi

echo "[OK] Pipeline finished. See $RFIN/$PROJECT_NAME/0_stats/"
