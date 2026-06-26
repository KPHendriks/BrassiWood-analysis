#!/bin/bash
# 1b2.Create_MSAs_by_gene.sh
# Per-gene worker: OMM-MACSE -> trimAl -> IQ-TREE 3
# For MaaS36 direct execution via GNU parallel

set -euo pipefail
IFS=$'\n\t'

gene="${1:?Usage: 1b2.Create_MSAs_by_gene.sh <gene>}"

: "${project_name:?}"
: "${genome:?}"
: "${method_alignment:?}"
: "${ROOT:?}"
: "${RINT:?}"
: "${RFIN:?}"
: "${DATA:?}"
: "${SCRIPTS:?}"
: "${MACSE_DIR:?}"
: "${CONS_LOCI_DIR:?}"
: "${cpu_limit:?}"
: "${MACSE_MEM_GB:?}"
: "${IQTREE_MEM_GB:?}"
: "${IQTREE_THREADS:?}"
: "${RUNLOG:?}"

TRIMAL_BIN="${TRIMAL_BIN:-$HOME/trimal/source/trimal}"
IQTREE_BIN="${IQTREE_BIN:-$HOME/iqtree-3.0.1-Linux/bin/iqtree3}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1"; exit 127; }
}
need_exec() {
  [[ -x "$1" ]] || { echo "Missing executable: $1"; exit 127; }
}

need_exec "$TRIMAL_BIN"
need_exec "$IQTREE_BIN"
need_cmd /usr/bin/time

echo "Project: $project_name | Genome: $genome | Gene: $gene | Method: $method_alignment"

out_align_dir="$RFIN/$project_name/1_results_gene_alignments_${method_alignment}/$genome"
out_trimal_dir="$RFIN/$project_name/2_results_gene_alignments_${method_alignment}_trimal/$genome"
out_tree_dir="$RFIN/$project_name/3_results_iqtree_gene_trees/$genome"
log_dir="$RFIN/$project_name/0_stats/gene_logs/$genome/$gene"
fail_dir="$RFIN/$project_name/0_stats/gene_logs/$genome/${gene}_FAILED"

mkdir -p "$out_align_dir" "$out_trimal_dir" "$out_tree_dir" "$log_dir"

# ----------- timing / logging accumulators -----------
start_epoch=$(date +%s)
start_iso=$(date -Iseconds)
cpu_user_sum=0
cpu_sys_sum=0
max_rss=0
host=$(hostname -s || echo unknown)
pid=$$

FAIL_STAGE="ok"
FAIL_REASON="ok"

sanitize_csv_field() {
  echo "$1" | tr '\n' ' ' | tr ',' ';'
}

append_runlog() {
  local exit_code="$1"
  local end_epoch end_iso wall
  end_epoch=$(date +%s)
  end_iso=$(date -Iseconds)
  wall=$(( end_epoch - start_epoch ))

  local fail_stage_clean fail_reason_clean
  fail_stage_clean=$(sanitize_csv_field "$FAIL_STAGE")
  fail_reason_clean=$(sanitize_csv_field "$FAIL_REASON")

  {
    flock 200 2>/dev/null || true
    if [[ ! -s "$RUNLOG" ]]; then
      echo "project,genome,gene,start_iso,end_iso,wall_sec,cpu_user_sec,cpu_sys_sec,max_rss_kb,exit_code,host,pid,fail_stage,fail_reason" >> "$RUNLOG"
    fi
    echo "${project_name},${genome},${gene},${start_iso},${end_iso},${wall},${cpu_user_sum},${cpu_sys_sum},${max_rss},${exit_code},${host},${pid},${fail_stage_clean},${fail_reason_clean}" >> "$RUNLOG"
  } 200>>"$RUNLOG"
}

preserve_failure_artifacts() {
  mkdir -p "$fail_dir"
  cp -f "$log_dir"/* "$fail_dir/" 2>/dev/null || true
  cp -f "${gene}.fasta" "$fail_dir/" 2>/dev/null || true
  cp -rf omm_macse_out "$fail_dir/" 2>/dev/null || true
  cp -f "${gene}_macse_aligned.fasta" "$fail_dir/" 2>/dev/null || true
  cp -f "${gene}_macse_aligned_trimal.fasta" "$fail_dir/" 2>/dev/null || true
  cp -f "${gene}_iqtree.log" "$fail_dir/" 2>/dev/null || true
  cp -f "${gene}_iqtree.iqtree" "$fail_dir/" 2>/dev/null || true
}

log_and_exit() {
  local code="$1"
  echo "[ERR] ${gene} | stage=${FAIL_STAGE} | reason=${FAIL_REASON}"
  preserve_failure_artifacts
  append_runlog "$code"
  exit "$code"
}

trap 'FAIL_STAGE="${FAIL_STAGE:-unknown}"; FAIL_REASON="${FAIL_REASON:-untrapped_shell_error_at_line_$LINENO}"; log_and_exit $?' ERR

run_timed() {
  local label="$1"; shift
  local stdout_file="$log_dir/${label}.stdout.log"
  local stderr_file="$log_dir/${label}.stderr.log"
  local time_file="$log_dir/${label}.time.log"

  if ! /usr/bin/time -v -o "$time_file" "$@" >"$stdout_file" 2>"$stderr_file"; then
    FAIL_STAGE="$label"
    FAIL_REASON="${label}_command_failed"
    return 1
  fi

  local u s rss
  u=$(grep -E '^User time' "$time_file" | awk '{print $4}' || echo 0)
  s=$(grep -E '^System time' "$time_file" | awk '{print $4}' || echo 0)
  rss=$(grep -E '^Maximum resident set size' "$time_file" | awk '{print $6}' || echo 0)

  cpu_user_sum=$(awk -v a="$cpu_user_sum" -v b="$u" 'BEGIN{printf("%.3f", a+b)}')
  cpu_sys_sum=$(awk -v a="$cpu_sys_sum" -v b="$s" 'BEGIN{printf("%.3f", a+b)}')

  if [[ "$rss" =~ ^[0-9]+$ ]]; then
    if (( rss > max_rss )); then
      max_rss=$rss
    fi
  fi
}

# ----------- FAST SKIP IF CURRENT OUTPUT ALREADY EXISTS -----------
if [[ -s "$out_tree_dir/${gene}_iqtree.treefile" ]]; then
  echo "[SKIP] ${gene}: IQ-TREE treefile already present"
  FAIL_STAGE="skip_existing"
  FAIL_REASON="treefile_already_present"
  append_runlog 0
  exit 0
fi

# ----------- OPTIONAL INHERIT FROM AN EARLIER PROJECT -----------
if [[ "${inherit_project:-none}" != "none" ]]; then
  inh_align="$RFIN/$inherit_project/1_results_gene_alignments_${method_alignment}/$genome/${gene}_macse_aligned.fasta"
  inh_trimal="$RFIN/$inherit_project/2_results_gene_alignments_${method_alignment}_trimal/$genome/${gene}_macse_aligned_trimal.fasta"
  inh_tree="$RFIN/$inherit_project/3_results_iqtree_gene_trees/$genome/${gene}_iqtree.treefile"
  inh_log="$RFIN/$inherit_project/3_results_iqtree_gene_trees/$genome/${gene}_iqtree.log"

  if [[ -s "$inh_tree" ]]; then
    echo "[INHERIT] ${gene}: copying finished outputs from project ${inherit_project}"
    [[ -s "$inh_align"  ]] && cp -f "$inh_align" "$out_align_dir/"
    [[ -s "$inh_trimal" ]] && cp -f "$inh_trimal" "$out_trimal_dir/"
    [[ -s "$inh_tree"   ]] && cp -f "$inh_tree" "$out_tree_dir/"
    [[ -s "$inh_log"    ]] && cp -f "$inh_log" "$out_tree_dir/"
    FAIL_STAGE="inherit"
    FAIL_REASON="copied_from_previous_project"
    append_runlog 0
    exit 0
  fi
fi

# ----------- WORKDIR -----------
workdir="/tmp/$project_name/$genome/$gene"
rm -rf "$workdir"
mkdir -p "$workdir"
cd "$workdir"

# ----------- STEP 1: input FASTA -----------
src=""
if [[ -s "$CONS_LOCI_DIR/${gene}_consensus.fasta" ]]; then
  src="$CONS_LOCI_DIR/${gene}_consensus.fasta"
elif [[ -s "$CONS_LOCI_DIR/${gene}_intronerated_consensus.fasta" ]]; then
  src="$CONS_LOCI_DIR/${gene}_intronerated_consensus.fasta"
else
  FAIL_STAGE="input"
  FAIL_REASON="no_per_locus_fasta_found"
  log_and_exit 10
fi

cp "$src" "${gene}.fasta"

nseq=$(grep -c '^>' "${gene}.fasta" || true)
nbp=$(awk 'BEGIN{n=0} !/^>/{gsub(/[[:space:]]/,""); n+=length($0)} END{print n+0}' "${gene}.fasta")

echo "gene=${gene}"              >  "$log_dir/input_summary.txt"
echo "source_fasta=${src}"       >> "$log_dir/input_summary.txt"
echo "n_sequences=${nseq}"       >> "$log_dir/input_summary.txt"
echo "total_bp=${nbp}"           >> "$log_dir/input_summary.txt"

# ----------- STEP 2: OMM-MACSE -----------
if [[ "$method_alignment" != "macse" ]]; then
  FAIL_STAGE="config"
  FAIL_REASON="unsupported_method_alignment_${method_alignment}"
  log_and_exit 11
fi

cp "$MACSE_DIR/omm_macse_v10.02.sif" .
chmod +x omm_macse_v10.02.sif || true

echo "[MACSE] ${gene}.fasta -> omm_macse_out/${gene}_final_align_NT.aln (mem ${MACSE_MEM_GB}g)"
run_timed macse ./omm_macse_v10.02.sif \
  --in_seq_file "${gene}.fasta" \
  --out_dir "omm_macse_out" \
  --out_file_prefix "${gene}" \
  --java_mem "${MACSE_MEM_GB}g"

aligned="omm_macse_out/${gene}_final_align_NT.aln"
if [[ ! -s "$aligned" ]]; then
  FAIL_STAGE="macse"
  FAIL_REASON="macse_no_alignment_output"
  log_and_exit 3
fi

cp "$aligned" "${gene}_macse_aligned.fasta"

# ----------- STEP 3: trimAl -----------
run_timed trimal "$TRIMAL_BIN" \
  -in "${gene}_macse_aligned.fasta" \
  -out "${gene}_macse_aligned_trimal.fasta" \
  -htmlout "${gene}_macse_aligned_trimal.html" \
  -resoverlap 0.40 \
  -seqoverlap 40 \
  -gt 0.40

if [[ ! -s "${gene}_macse_aligned_trimal.fasta" ]]; then
  FAIL_STAGE="trimal"
  FAIL_REASON="trimal_no_output"
  log_and_exit 20
fi

# ----------- STEP 4: IQ-TREE 3 -----------
run_timed iqtree "$IQTREE_BIN" \
  -s "${gene}_macse_aligned_trimal.fasta" \
  -m GTR+F+R3 \
  -nt "${IQTREE_THREADS}" \
  --mem "${IQTREE_MEM_GB}G" \
  -keep-ident \
  --prefix "${gene}_iqtree" \
  --alrt 1000 \
  -B 1000

if [[ ! -s "${gene}_iqtree.treefile" ]]; then
  FAIL_STAGE="iqtree"
  FAIL_REASON="iqtree_no_treefile"
  log_and_exit 30
fi

# ----------- STEP 5: copy out -----------
cp -f "${gene}_macse_aligned.fasta"        "$out_align_dir/"
cp -f "${gene}_macse_aligned_trimal.fasta" "$out_trimal_dir/"
cp -f "${gene}_iqtree.treefile"            "$out_tree_dir/"
cp -f "${gene}_iqtree.log"                 "$out_tree_dir/"

[[ -s "${gene}_iqtree.iqtree"  ]] && cp -f "${gene}_iqtree.iqtree"  "$out_tree_dir/"
[[ -s "${gene}_iqtree.contree" ]] && cp -f "${gene}_iqtree.contree" "$out_tree_dir/"
[[ -s "${gene}_macse_aligned_trimal.html" ]] && cp -f "${gene}_macse_aligned_trimal.html" "$out_trimal_dir/"

# ----------- STEP 6: cleanup -----------
FAIL_STAGE="ok"
FAIL_REASON="completed"
append_runlog 0

cd /tmp
rm -rf "$workdir"

echo "[DONE] $gene"
