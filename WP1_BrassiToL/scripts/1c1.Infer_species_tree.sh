#!/bin/bash
# Species tree pipeline (1c):
#  - Step 2: summarize per-gene trees (moved from 1b1)
#  - Step 3: collapse weak branches in gene trees (aLRT/UFboot/BL thresholds)
#  - Step 4–5: TreeShrink (per-locus outlier removal)
#  - Step 6: ASTRAL-IV species tree (coalescent topology)
#  - Step 7: IQ-TREE 3 concordance factors (gCF, sCF) on ASTRAL topology
#  - Step 8: MCMCtree divergence time estimation (SortaDate + baseml + multi-profile MCMC)

set -euo pipefail
IFS=$'\n\t'
trap 'echo "[ERR] at line $LINENO"' ERR

# =====================================================
# ----------------------- STEP 1 ----------------------
# PREPARE: SET VARIABLES, CREATE DIRECTORIES, ETC.
# =====================================================

# ----------- CONFIG / ROOTS -----------
ROOT="."
BIG="${ROOT%/*}/BrassiWood_big_data"
RINT="$ROOT/WP1_BrassiToL/results_intermediate"
RFIN="$ROOT/WP1_BrassiToL/results_final"
DATA="$ROOT/WP1_BrassiToL/data"
SCRIPTS="$ROOT/WP1_BrassiToL/scripts"

# Project choices
export project_name="2026-06-16_BrassiWood_study_v1"
export sample_list="WP1_BrassiToL/data/SampleListStart.txt"
export sample_outgroup="PAFTOL_014331"

# Genomes to use (must match what 1b2 produced)
genome_list=(B764 A353)

# Alignment engine label (to find folders)
export method_alignment="macse"      # macse | mafft
export cpu_limit=10

# Binaries
IQTREE_BIN="${IQTREE_BIN:-$HOME/iqtree-3.0.1-Linux/bin/iqtree3}"
ASTRAL_BIN="${ASTRAL_BIN:-$HOME/ASTER-MacOS/bin/astral4}"
TREEPL_BIN="${TREEPL_BIN:-$HOME/.homebrew/bin/treePL}"

# PAML and MCMCtree locations
PAML_BIN_DIR="${PAML_BIN_DIR:-$HOME/paml/bin}"
export PATH="$PAML_BIN_DIR:$PATH"

MCMCTREE_BIN="${MCMCTREE_BIN:-$PAML_BIN_DIR/mcmctree}"
BASEML_BIN="${BASEML_BIN:-$PAML_BIN_DIR/baseml}"

if [[ ! -x "$MCMCTREE_BIN" ]]; then
  echo "[ERR] mcmctree not found or not executable: $MCMCTREE_BIN" | tee -a "$RUNLOG"
  exit 1
fi

if ! command -v baseml >/dev/null 2>&1; then
  echo "[ERR] baseml not found in PATH. Current PATH: $PATH" | tee -a "$RUNLOG"
  exit 1
fi

# Resources / knobs
export IQTREE_THREADS="AUTO"
export UFB_REPS=1000                    # 0 = skip constrained UFBoot; e.g. set 1000 to run
export SCF_QUARTETS=100              # sCF quartets per internal branch

# ----------- LOGGING -----------
RUNLOG="$RFIN/$project_name/0_stats/1c_species_tree_runlog.txt"
mkdir -p "$RFIN/$project_name/0_stats"
echo "[START] $(date)" | tee "$RUNLOG"

# ----------- DEFAULT GENE LISTS -----------
mapfile -t gene_list_B764 < <(
  grep -a "^>" "$DATA/ref-at_orf_minus_trailing_STOP_codons_exons_binned_to_by_genes.fasta" \
    | sed -E 's/^>//; s/[[:space:]].*$//; s/^.*-//' \
    | sort -u
)

mapfile -t gene_list_A353 < <(
  grep -a "^>" "$DATA/A353_NewTargets_refgenome_minus_trailing_STOP_codons.fasta" \
    | sed -E 's/^>//; s/[[:space:]].*$//; s/^.*-//' \
    | tr -d '\r' \
    | sort -u
)

# ----------- DIRECTORIES (with numbering) -----------
COLLAPSE_DIR="$RFIN/$project_name/4_gene_trees_collapsed"
TS_DIR="$RFIN/$project_name/5_results_treeshrink"
ASTRAL_DIR="$RFIN/$project_name/6_results_astral_species_tree"
IQSP_DIR="$RFIN/$project_name/7_results_iqtree_species_tree"
CALIB_DIR="$RFIN/$project_name/8_results_species_tree_calibrated"
PUBLISH_DIR="$RFIN/$project_name/9_results_species_tree_for_publication"
mkdir -p "$COLLAPSE_DIR" "$TS_DIR" "$ASTRAL_DIR" "$IQSP_DIR" "$CALIB_DIR" "$PUBLISH_DIR"


# =====================================================
# ----------------------- STEP 2 ----------------------
# SUMMARIZE IQ-TREE SUPPORTS & BRANCH LENGTHS
# =====================================================

# See script 1b1 and ideally use the same numbers to follow through on the stats output from that script.

python3 "$SCRIPTS/custom_scripts/summarize_iqtree_supports.py" \
  --input-glob "$RFIN/$project_name/3_results_iqtree_gene_trees/*/*_iqtree.treefile" \
  --outdir "$RFIN/$project_name/0_stats" \
  --alrt-threshold 70 \
  --ufboot-threshold 80 \
  --bl-threshold 1e-5 \
  --log10-bl-hist || true


# =====================================================
# ----------------------- STEP 3 ----------------------
# COLLAPSE GENE TREES BY THRESHOLDS
# =====================================================

python3 "$SCRIPTS/custom_scripts/collapse_gene_trees_by_support.py" \
  --input-glob "$RFIN/$project_name/3_results_iqtree_gene_trees/*/*_iqtree.treefile" \
  --outdir "$COLLAPSE_DIR" \
  --alrt-threshold 70 \
  --ufboot-threshold 80 \
  --bl-threshold 1e-5 \
  --logic any


# =====================================================
# ----------------------- STEP 4 ----------------------
# PREP FILES FOR TREESHRINK (use collapsed trees)
# =====================================================

echo "[INFO] Preparing TreeShrink inputs" | tee -a "$RUNLOG"
rm -rf "$TS_DIR"
mkdir -p "$TS_DIR/output_alignments" "$TS_DIR/output_trees" "$TS_DIR/output_logs"

n_ok=0
n_missing_aln=0
n_missing_tree=0

for genome in "${genome_list[@]}"; do
  echo "[INFO] Genome: $genome" | tee -a "$RUNLOG"
  declare -n glist="gene_list_${genome}"

  for gene in "${glist[@]}"; do
    aln_path="$RFIN/$project_name/2_results_gene_alignments_${method_alignment}_trimal/${genome}/${gene}_macse_aligned_trimal.fasta"
    tree_path="$COLLAPSE_DIR/${gene}_iqtree_collapsed.tree"

    if [[ ! -s "$aln_path" ]]; then
      ((n_missing_aln++))
      continue
    fi

    if [[ ! -s "$tree_path" ]]; then
      ((n_missing_tree++))
      continue
    fi

    gdir="$TS_DIR/$gene"
    mkdir -p "$gdir"
    cp "$aln_path" "$gdir/input.fasta"
    cp "$tree_path" "$gdir/input.tree"
    ((n_ok++))
  done
done

echo "[INFO] TreeShrink loci prepared: $n_ok" | tee -a "$RUNLOG"
echo "[INFO] Missing alignments: $n_missing_aln" | tee -a "$RUNLOG"
echo "[INFO] Missing trees: $n_missing_tree" | tee -a "$RUNLOG"
echo "[INFO] N loci for TreeShrink: $(find "$TS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name output_alignments ! -name output_trees ! -name output_logs | wc -l)" | tee -a "$RUNLOG"


# =====================================================
# ----------------------- STEP 5 ----------------------
# RUN TREESHRINK
# =====================================================

python3 ~/TreeShrink/run_treeshrink.py \
  -i "$TS_DIR" \
  -t input.tree \
  -a input.fasta \
  -f > "$TS_DIR/treeshrinklog.txt"

# Move outputs to flat folders and clean per-locus dirs
for d in "$TS_DIR"/*/; do
  gene="$(basename "$d")"
  [[ "$gene" == "output_alignments" || "$gene" == "output_trees" || "$gene" == "output_logs" ]] && continue
  [[ -f "$d/output.fasta" ]] && mv "$d/output.fasta" "$TS_DIR/output_alignments/${gene}_output_treeshrink.fasta"
  [[ -f "$d/output.tree"  ]] && mv "$d/output.tree"  "$TS_DIR/output_trees/${gene}_output_treeshrink.tree"
  [[ -f "$d/output.txt"   ]] && mv "$d/output.txt"   "$TS_DIR/output_logs/${gene}_output_treeshrink.txt"
  rm -rf "$d"
done


# =====================================================
# ----------------------- STEP 6 ----------------------
# ASTRAL-IV SPECIES TREE
# =====================================================

cat "$TS_DIR/output_trees/"*_output_treeshrink.tree > "$ASTRAL_DIR/gene_trees_combined.tree"
echo "[INFO] N gene trees for ASTRAL: $(wc -l < "$ASTRAL_DIR/gene_trees_combined.tree")" | tee -a "$RUNLOG"

# Increase stack size just in case (common for large trees in C++ binaries)
ulimit -s unlimited

# Run ASTRAL-IV
"$ASTRAL_BIN" \
  --input "$ASTRAL_DIR/gene_trees_combined.tree" \
  --output "$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.tree" \
  --thread 50 \
  --root "$sample_outgroup" \
  --support 2 \
  &> "$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.log"

#The above call is after increasing RAM, below is from before and can be deleted once it's confirmed the above works better  
#"$ASTRAL_BIN" \
#  --output "$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.tree" \
#  --support 2 \
#  --input "$ASTRAL_DIR/gene_trees_combined.tree" \
#  --root "$sample_outgroup" \
#  --thread 40 \
#  &> "$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.log"


# =====================================================
# ----------------------- STEP 7 ----------------------
# IQ-TREE CONCORDANCE FACTORS (gCF + sCF) ON ASTRAL TREE
# =====================================================

echo "[INFO] STEP 7: Calculating IQ-TREE concordance factors on ASTRAL species tree" | tee -a "$RUNLOG"

# ---- 7.0 Inputs / outputs ---------------------------------------------

CF_PREFIX="$IQSP_DIR/${project_name}_BrassiToL_concordance_factors"
CF_GENE_TREES="$IQSP_DIR/${project_name}_gene_trees_for_gcf.tree"
CF_REFERENCE_TREE="$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.tree"
CF_ALIGNMENT_DIR="$TS_DIR/output_alignments"

# ---- 7.1 Sanity checks ------------------------------------------------

if [[ ! -s "$CF_REFERENCE_TREE" ]]; then
  echo "[ERR] Reference ASTRAL tree not found or empty: $CF_REFERENCE_TREE" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -d "$CF_ALIGNMENT_DIR" ]]; then
  echo "[ERR] TreeShrink alignment directory not found: $CF_ALIGNMENT_DIR" | tee -a "$RUNLOG"
  exit 1
fi

N_ALN=$(find "$CF_ALIGNMENT_DIR" -maxdepth 1 -type f -name "*_output_treeshrink.fasta" | wc -l)
if (( N_ALN == 0 )); then
  echo "[ERR] No TreeShrink alignments found in: $CF_ALIGNMENT_DIR" | tee -a "$RUNLOG"
  exit 1
fi

N_COLLAPSED=$(find "$COLLAPSE_DIR" -maxdepth 1 -type f -name "*_iqtree_collapsed.tree" | wc -l)
if (( N_COLLAPSED == 0 )); then
  echo "[ERR] No collapsed gene trees found in: $COLLAPSE_DIR" | tee -a "$RUNLOG"
  exit 1
fi

echo "[INFO] ASTRAL reference tree: $CF_REFERENCE_TREE" | tee -a "$RUNLOG"
echo "[INFO] TreeShrink alignment directory: $CF_ALIGNMENT_DIR" | tee -a "$RUNLOG"
echo "[INFO] Collapsed gene trees found: $N_COLLAPSED" | tee -a "$RUNLOG"
echo "[INFO] TreeShrink alignments found: $N_ALN" | tee -a "$RUNLOG"
echo "[INFO] sCF quartets per branch: $SCF_QUARTETS" | tee -a "$RUNLOG"

# ---- 7.2 Combine collapsed gene trees for gCF -------------------------

find "$COLLAPSE_DIR" -maxdepth 1 -type f -name "*_iqtree_collapsed.tree" | sort > "$IQSP_DIR/collapsed_gene_tree_list.txt"

if [[ ! -s "$IQSP_DIR/collapsed_gene_tree_list.txt" ]]; then
  echo "[ERR] Collapsed gene tree list is empty." | tee -a "$RUNLOG"
  exit 1
fi

: > "$CF_GENE_TREES"
while read -r tf; do
  [[ -z "$tf" ]] && continue
  cat "$tf" >> "$CF_GENE_TREES"
  printf '\n' >> "$CF_GENE_TREES"
done < "$IQSP_DIR/collapsed_gene_tree_list.txt"

N_GCF_TREES=$(wc -l < "$CF_GENE_TREES")
echo "[INFO] Combined collapsed gene trees written to: $CF_GENE_TREES" | tee -a "$RUNLOG"
echo "[INFO] Number of gene trees in combined file: $N_GCF_TREES" | tee -a "$RUNLOG"

if (( N_GCF_TREES == 0 )); then
  echo "[ERR] Combined gene tree file is empty: $CF_GENE_TREES" | tee -a "$RUNLOG"
  exit 1
fi

# ---- 7.3 Run IQ-TREE concordance factor analysis ----------------------

echo "[INFO] Running IQ-TREE concordance factor analysis..." | tee -a "$RUNLOG"

"$IQTREE_BIN" \
  -s "$CF_ALIGNMENT_DIR/" \
  -t "$CF_REFERENCE_TREE" \
  --gcf "$CF_GENE_TREES" \
  --scf "$SCF_QUARTETS" \
  --prefix "$CF_PREFIX" \
  -n 0 \
  -redo \
  -nt 20 \
  &> "${CF_PREFIX}.log"

# ---- 7.4 Check expected outputs ---------------------------------------

if [[ -s "${CF_PREFIX}.cf.tree" ]]; then
  echo "[INFO] Concordance factor tree written to: ${CF_PREFIX}.cf.tree" | tee -a "$RUNLOG"
else
  echo "[ERR] Expected concordance factor tree not found: ${CF_PREFIX}.cf.tree" | tee -a "$RUNLOG"
  echo "[ERR] Check log: ${CF_PREFIX}.log" | tee -a "$RUNLOG"
  exit 1
fi

if [[ -s "${CF_PREFIX}.cf.stat" ]]; then
  echo "[INFO] Concordance factor statistics written to: ${CF_PREFIX}.cf.stat" | tee -a "$RUNLOG"
else
  echo "[WARN] Concordance factor statistics file not found: ${CF_PREFIX}.cf.stat" | tee -a "$RUNLOG"
fi

echo "[INFO] IQ-TREE concordance factor analysis completed." | tee -a "$RUNLOG"
echo "[INFO] Main CF outputs:" | tee -a "$RUNLOG"
echo "[INFO]   Tree: ${CF_PREFIX}.cf.tree" | tee -a "$RUNLOG"
echo "[INFO]   Stats: ${CF_PREFIX}.cf.stat" | tee -a "$RUNLOG"
echo "[INFO]   Log: ${CF_PREFIX}.log" | tee -a "$RUNLOG"




# =====================================================
# ----------------------- STEP 8 ----------------------
# SHARED LOCUS RANKING + SUBSET PREPARATION FOR DATING
# =====================================================

echo "[INFO] STEP 8: shared locus ranking and subset preparation for downstream dating" | tee -a "$RUNLOG"

# -------------------- USER SETTINGS -------------------

# Toggle downstream dating methods (used later in Steps 9–10)
RUN_MCMCTREE="${RUN_MCMCTREE:-yes}"
RUN_LSD2="${RUN_LSD2:-yes}"

# Reusable locus subset sizes to prepare once here
# Example: 20,50,100,200
DATE_SUBSET_SIZES_CSV="${DATE_SUBSET_SIZES_CSV:-20,50,100,200}"

# Method-specific defaults (used later; defined here so all dating settings live together)
MCT_N_LOCI="${MCT_N_LOCI:-20}"
LSD2_N_LOCI="${LSD2_N_LOCI:-200}"

# -------------------- PATHS ---------------------------

DATESEL_DIR="$CALIB_DIR/00_dating_locus_selection"
DATESEL_TREE_DIR="$DATESEL_DIR/01_treefiles"
DATESEL_CAND_ALIGN="$DATESEL_DIR/02_candidate_alignments"
DATESEL_CAND_TREES="$DATESEL_DIR/03_candidate_trees"
DATESEL_RESULTS="$DATESEL_DIR/04_results"
DATESEL_SUBSETS="$DATESEL_DIR/05_selected_subsets"

mkdir -p \
  "$DATESEL_TREE_DIR" \
  "$DATESEL_CAND_ALIGN" \
  "$DATESEL_CAND_TREES" \
  "$DATESEL_RESULTS" \
  "$DATESEL_SUBSETS"

# Main shared files
DATESEL_ASTRAL_TREE_RAW="$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.tree"
DATESEL_ASTRAL_TREE_CLEAN="$DATESEL_TREE_DIR/${project_name}_astral_clean_topology.nwk"
DATESEL_INGROUP_SPECIES_TREE="$DATESEL_RESULTS/${project_name}_astral_clean_topology_ingroup.nwk"

DATESEL_CANDIDATE_TSV="$DATESEL_RESULTS/candidate_loci.tsv"
DATESEL_FAIL_TSV="$DATESEL_RESULTS/candidate_loci_failures.tsv"
DATESEL_RANKED_FILE="$DATESEL_RESULTS/ranked_loci.txt"
DATESEL_MASTER_SELECTED_TSV="$DATESEL_RESULTS/selected_loci_master.tsv"
DATESEL_MASTER_SELECTED_LIST="$DATESEL_RESULTS/selected_loci_master_paths.txt"
DATESEL_COVERAGE_TSV="$DATESEL_RESULTS/selected_loci_coverage_progress.tsv"
DATESEL_UNCOVERED_TSV="$DATESEL_RESULTS/uncovered_taxa_after_selection.tsv"
DATESEL_SUBSET_SUMMARY_TSV="$DATESEL_RESULTS/subset_summary.tsv"

TS_ALN_DIR="$TS_DIR/output_alignments"
TS_TREE_DIR="$TS_DIR/output_trees"
OUTGROUP_TIP="${sample_outgroup}"

export DATESEL_ASTRAL_TREE_RAW
export DATESEL_ASTRAL_TREE_CLEAN
export DATESEL_INGROUP_SPECIES_TREE
export DATESEL_CAND_ALIGN
export DATESEL_CAND_TREES
export DATESEL_CANDIDATE_TSV
export DATESEL_FAIL_TSV
export DATESEL_RANKED_FILE
export DATESEL_MASTER_SELECTED_TSV
export DATESEL_MASTER_SELECTED_LIST
export DATESEL_COVERAGE_TSV
export DATESEL_UNCOVERED_TSV
export DATESEL_SUBSET_SUMMARY_TSV
export DATESEL_SUBSETS
export TS_ALN_DIR
export TS_TREE_DIR
export OUTGROUP_TIP
export DATE_SUBSET_SIZES_CSV
export MCT_N_LOCI
export LSD2_N_LOCI

# -------------------- SANITY CHECKS -------------------

if [[ ! -s "$DATESEL_ASTRAL_TREE_RAW" ]]; then
  echo "[ERR] ASTRAL species tree not found or empty: $DATESEL_ASTRAL_TREE_RAW" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -d "$TS_ALN_DIR" ]]; then
  echo "[ERR] TreeShrink alignment directory not found: $TS_ALN_DIR" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -d "$TS_TREE_DIR" ]]; then
  echo "[ERR] TreeShrink tree directory not found: $TS_TREE_DIR" | tee -a "$RUNLOG"
  exit 1
fi

N_TS_ALN=$(find "$TS_ALN_DIR" -maxdepth 1 -type f -name "*_output_treeshrink.fasta" | wc -l)
if (( N_TS_ALN == 0 )); then
  echo "[ERR] No TreeShrink alignments found in: $TS_ALN_DIR" | tee -a "$RUNLOG"
  exit 1
fi

if [[ -z "${OUTGROUP_TIP:-}" ]]; then
  echo "[WARN] sample_outgroup is empty; ingroup pruning will be skipped and full cleaned ASTRAL topology will be used." | tee -a "$RUNLOG"
fi

for cmd in pxrmt; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERR] Required command not found in PATH: $cmd" | tee -a "$RUNLOG"
    exit 1
  fi
done

echo "[INFO] Shared dating-selection folder: $DATESEL_DIR" | tee -a "$RUNLOG"
echo "[INFO] ASTRAL raw tree:                $DATESEL_ASTRAL_TREE_RAW" | tee -a "$RUNLOG"
echo "[INFO] TreeShrink alignments:          $TS_ALN_DIR" | tee -a "$RUNLOG"
echo "[INFO] TreeShrink trees:               $TS_TREE_DIR" | tee -a "$RUNLOG"
echo "[INFO] N TreeShrink loci:              $N_TS_ALN" | tee -a "$RUNLOG"
echo "[INFO] Date subset sizes requested:    $DATE_SUBSET_SIZES_CSV" | tee -a "$RUNLOG"
echo "[INFO] Default MCMCtree loci:          $MCT_N_LOCI" | tee -a "$RUNLOG"
echo "[INFO] Default LSD2 loci:              $LSD2_N_LOCI" | tee -a "$RUNLOG"

# -------------------- 8.1 Clean ASTRAL tree ----------------------------

echo "[INFO] STEP 8.1: cleaning ASTRAL topology" | tee -a "$RUNLOG"

python3 - <<'PY'
from pathlib import Path
import os
import re

inp = Path(os.environ["DATESEL_ASTRAL_TREE_RAW"])
outp = Path(os.environ["DATESEL_ASTRAL_TREE_CLEAN"])

if not inp.exists():
    raise SystemExit(f"[ERR] Input ASTRAL tree not found: {inp}")

txt = inp.read_text(errors="replace")

# Remove square-bracket annotations/comments
txt = re.sub(r"\[[^\]]*\]", "", txt)

# Remove single quotes
txt = txt.replace("'", "")

# Remove whitespace/newlines
txt = re.sub(r"\s+", "", txt)

if ";" not in txt:
    raise SystemExit("[ERR] No semicolon found in tree file.")

# Keep only the first tree if multiple are present
txt = txt.split(";", 1)[0] + ";"

# Remove internal node labels/support labels that appear immediately after ')'
txt = re.sub(r"\)([^:(),;\[\]]+)(?=[:),;])", ")", txt)

# Remove branch lengths
txt = re.sub(r":[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "", txt)

if not txt.endswith(";"):
    raise SystemExit("[ERR] Cleaned tree does not end with ';' as expected.")

outp.write_text(txt + "\n")

print(f"[INFO] Cleaned ASTRAL topology written to: {outp}")
PY

# -------------------- 8.2 Build ingroup species tree -------------------

echo "[INFO] STEP 8.2: building ingroup species tree for coverage-aware locus selection" | tee -a "$RUNLOG"

python3 - <<'PY'
from pathlib import Path
import os
import subprocess

inp = Path(os.environ["DATESEL_ASTRAL_TREE_CLEAN"])
outp = Path(os.environ["DATESEL_INGROUP_SPECIES_TREE"])
outgroup = os.environ.get("OUTGROUP_TIP", "").strip()

txt = inp.read_text().strip()
if not txt.endswith(";"):
    txt += ";"

if outgroup:
    res = subprocess.run(
        ["pxrmt", "-t", str(inp), "-n", outgroup],
        capture_output=True, text=True, check=False
    )
    if res.returncode == 0 and res.stdout.strip():
        txt2 = res.stdout.strip()
        if not txt2.endswith(";"):
            txt2 += ";"
        outp.write_text(txt2 + "\n")
        print(f"[INFO] Wrote ingroup-only species tree to: {outp}")
        raise SystemExit(0)

outp.write_text(txt + "\n")
print(f"[INFO] Could not prune outgroup cleanly (or no outgroup provided); using cleaned topology as-is: {outp}")
PY

# -------------------- 8.3 Collect candidate loci -----------------------

echo "[INFO] STEP 8.3: collecting candidate loci with relaxed validity checks" | tee -a "$RUNLOG"

rm -f \
  "$DATESEL_CANDIDATE_TSV" \
  "$DATESEL_FAIL_TSV" \
  "$DATESEL_RANKED_FILE" \
  "$DATESEL_MASTER_SELECTED_TSV" \
  "$DATESEL_MASTER_SELECTED_LIST" \
  "$DATESEL_COVERAGE_TSV" \
  "$DATESEL_UNCOVERED_TSV" \
  "$DATESEL_SUBSET_SUMMARY_TSV"

rm -rf "$DATESEL_CAND_ALIGN" "$DATESEL_CAND_TREES"
mkdir -p "$DATESEL_CAND_ALIGN" "$DATESEL_CAND_TREES"

printf "gene\talignment\ttree\tn_taxa_alignment\taln_len\tn_taxa_tree\tn_taxa_species_tree_overlap\tn_taxa_aln_tree_shared\tfraction_species_tree_covered\tfraction_alignment_in_tree\n" > "$DATESEL_CANDIDATE_TSV"
printf "gene\treason\tdetail\n" > "$DATESEL_FAIL_TSV"

python3 - <<'PY'
from pathlib import Path
import os
import re
import shutil

ts_aln_dir = Path(os.environ["TS_ALN_DIR"])
ts_tree_dir = Path(os.environ["TS_TREE_DIR"])
cand_align = Path(os.environ["DATESEL_CAND_ALIGN"])
cand_trees = Path(os.environ["DATESEL_CAND_TREES"])
cand_tsv = Path(os.environ["DATESEL_CANDIDATE_TSV"])
fail_tsv = Path(os.environ["DATESEL_FAIL_TSV"])
species_tree_fp = Path(os.environ["DATESEL_INGROUP_SPECIES_TREE"])

def read_fasta(fp):
    seqs = {}
    name = None
    chunks = []
    with open(fp) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.replace(" ", "").replace("\t", ""))
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    return seqs

def write_fasta(seqs, fp):
    with open(fp, "w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n{seq}\n")

def guess_tree(gene):
    candidates = [
        ts_tree_dir / f"{gene}_output_treeshrink.tree",
        ts_tree_dir / f"{gene}_output_treeshrink.tre",
        ts_tree_dir / f"{gene}_output_treeshrink.treefile",
        ts_tree_dir / f"{gene}.tree",
        ts_tree_dir / f"{gene}.tre",
        ts_tree_dir / f"{gene}.treefile",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None

def extract_leaf_names_from_newick(txt):
    txt = re.sub(r"\[[^\]]*\]", "", txt)
    txt = txt.strip()
    tokens = re.findall(r"(?<=[(,])\s*([^():;,]+)", txt)
    return {t.strip().strip("'") for t in tokens if t.strip()}

species_tree_taxa = extract_leaf_names_from_newick(species_tree_fp.read_text())
alns = sorted(ts_aln_dir.glob("*_output_treeshrink.fasta"))

n_total = len(alns)
n_with_tree = 0
n_rectangular = 0
n_min_taxa = 0
n_tree_taxa_ok = 0
n_species_overlap_ok = 0
kept = []
failed = []

for aln in alns:
    gene = aln.name.replace("_output_treeshrink.fasta", "")
    tree = guess_tree(gene)
    if tree is None:
        failed.append((gene, "missing_tree", ""))
        continue
    n_with_tree += 1

    seqs = read_fasta(aln)
    if not seqs:
        failed.append((gene, "empty_alignment", ""))
        continue

    lengths = {len(v) for v in seqs.values()}
    if len(lengths) != 1:
        failed.append((gene, "alignment_not_rectangular", ""))
        continue
    n_rectangular += 1

    aln_len = lengths.pop()
    aln_taxa_set = set(seqs.keys())

    if len(aln_taxa_set) < 3:
        failed.append((gene, "too_few_taxa_in_alignment", str(len(aln_taxa_set))))
        continue
    n_min_taxa += 1

    tree_txt = tree.read_text(errors="replace")
    tree_taxa = extract_leaf_names_from_newick(tree_txt)
    if len(tree_taxa) < 3:
        failed.append((gene, "too_few_taxa_in_tree", str(len(tree_taxa))))
        continue
    n_tree_taxa_ok += 1

    shared_aln_tree = aln_taxa_set & tree_taxa
    if len(shared_aln_tree) < 3:
        failed.append((gene, "too_few_shared_taxa_alignment_tree", str(len(shared_aln_tree))))
        continue

    species_overlap = aln_taxa_set & species_tree_taxa
    if len(species_overlap) < 3:
        failed.append((gene, "too_few_taxa_overlapping_species_tree", str(len(species_overlap))))
        continue
    n_species_overlap_ok += 1

    aln_dst = cand_align / f"{gene}.fasta"
    tree_dst = cand_trees / f"{gene}.tre"

    write_fasta(seqs, aln_dst)
    shutil.copyfile(tree, tree_dst)

    frac_species = len(species_overlap) / len(species_tree_taxa) if species_tree_taxa else 0.0
    frac_aln_tree = len(shared_aln_tree) / len(aln_taxa_set) if aln_taxa_set else 0.0

    kept.append({
        "gene": gene,
        "alignment": str(aln_dst),
        "tree": str(tree_dst),
        "n_taxa_alignment": len(aln_taxa_set),
        "aln_len": aln_len,
        "n_taxa_tree": len(tree_taxa),
        "n_taxa_species_tree_overlap": len(species_overlap),
        "n_taxa_aln_tree_shared": len(shared_aln_tree),
        "fraction_species_tree_covered": frac_species,
        "fraction_alignment_in_tree": frac_aln_tree,
    })

with open(cand_tsv, "a") as out:
    for row in kept:
        out.write(
            f"{row['gene']}\t{row['alignment']}\t{row['tree']}\t"
            f"{row['n_taxa_alignment']}\t{row['aln_len']}\t{row['n_taxa_tree']}\t"
            f"{row['n_taxa_species_tree_overlap']}\t{row['n_taxa_aln_tree_shared']}\t"
            f"{row['fraction_species_tree_covered']:.8f}\t{row['fraction_alignment_in_tree']:.8f}\n"
        )

with open(fail_tsv, "a") as out:
    for gene, reason, detail in failed:
        detail = str(detail).replace("\t", " ").replace("\n", " | ")
        out.write(f"{gene}\t{reason}\t{detail}\n")

print(f"[INFO] Total TreeShrink loci found:                  {n_total}")
print(f"[INFO] Loci with matching gene tree:                 {n_with_tree}")
print(f"[INFO] Loci with rectangular alignment:              {n_rectangular}")
print(f"[INFO] Loci with >=3 alignment taxa:                 {n_min_taxa}")
print(f"[INFO] Loci with >=3 taxa in tree:                   {n_tree_taxa_ok}")
print(f"[INFO] Loci with >=3 taxa overlapping species tree:  {n_species_overlap_ok}")
print(f"[INFO] Final relaxed candidate loci:                 {len(kept)}")
print(f"[INFO] Candidate failure table written to:           {fail_tsv}")
PY

N_CAND=$(( $(wc -l < "$DATESEL_CANDIDATE_TSV") - 1 ))
echo "[INFO] Candidate loci after relaxed filtering: $N_CAND" | tee -a "$RUNLOG"

if (( N_CAND == 0 )); then
  echo "[ERR] No candidate loci remain after relaxed filtering." | tee -a "$RUNLOG"
  exit 1
fi

# -------------------- 8.4 Rank candidate loci --------------------------

echo "[INFO] STEP 8.4: ranking candidate loci for coverage-aware downstream selection" | tee -a "$RUNLOG"

python3 - <<'PY'
from pathlib import Path
import os
import csv

cand_tsv = Path(os.environ["DATESEL_CANDIDATE_TSV"])
ranked_file = Path(os.environ["DATESEL_RANKED_FILE"])

rows = []
with open(cand_tsv, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        row["n_taxa_alignment"] = int(row["n_taxa_alignment"])
        row["aln_len"] = int(row["aln_len"])
        row["n_taxa_tree"] = int(row["n_taxa_tree"])
        row["n_taxa_species_tree_overlap"] = int(row["n_taxa_species_tree_overlap"])
        row["n_taxa_aln_tree_shared"] = int(row["n_taxa_aln_tree_shared"])
        row["fraction_species_tree_covered"] = float(row["fraction_species_tree_covered"])
        row["fraction_alignment_in_tree"] = float(row["fraction_alignment_in_tree"])
        rows.append(row)

rows.sort(
    key=lambda r: (
        -r["n_taxa_species_tree_overlap"],
        -r["n_taxa_aln_tree_shared"],
        -r["n_taxa_alignment"],
        -r["aln_len"],
        r["gene"],
    )
)

with open(ranked_file, "w") as out:
    for r in rows:
        out.write(
            f"{r['tree']}\t"
            f"species_overlap={r['n_taxa_species_tree_overlap']}\t"
            f"aln_tree_shared={r['n_taxa_aln_tree_shared']}\t"
            f"n_taxa={r['n_taxa_alignment']}\t"
            f"aln_len={r['aln_len']}\n"
        )

print(f"[INFO] Ranking file written to: {ranked_file}")
print(f"[INFO] Ranked candidate loci: {len(rows)}")
PY

# -------------------- 8.5 Build master coverage-aware ordered locus list

echo "[INFO] STEP 8.5: building master coverage-aware locus order" | tee -a "$RUNLOG"

python3 - <<'PY'
from pathlib import Path
import os
import re

ranked_file = Path(os.environ["DATESEL_RANKED_FILE"])
cand_tsv = Path(os.environ["DATESEL_CANDIDATE_TSV"])
selected_tsv = Path(os.environ["DATESEL_MASTER_SELECTED_TSV"])
selected_list = Path(os.environ["DATESEL_MASTER_SELECTED_LIST"])
species_tree = Path(os.environ["DATESEL_INGROUP_SPECIES_TREE"])
uncovered_tsv = Path(os.environ["DATESEL_UNCOVERED_TSV"])
coverage_tsv = Path(os.environ["DATESEL_COVERAGE_TSV"])

def read_fasta_names(fp):
    names = []
    with open(fp) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names

def extract_leaf_names_from_newick(txt):
    txt = re.sub(r"\[[^\]]*\]", "", txt)
    txt = txt.strip()
    tokens = re.findall(r"(?<=[(,])\s*([^():;,]+)", txt)
    return {t.strip().strip("'") for t in tokens if t.strip()}

full_taxa = extract_leaf_names_from_newick(species_tree.read_text())

gene_to_aln = {}
gene_to_taxa = {}
gene_to_stats = {}

with open(cand_tsv) as fh:
    next(fh)
    for line in fh:
        (
            gene, aln, tree, n_taxa_alignment, aln_len, n_taxa_tree,
            n_taxa_species_tree_overlap, n_taxa_aln_tree_shared,
            fraction_species_tree_covered, fraction_alignment_in_tree
        ) = line.rstrip("\n").split("\t")

        aln_path = Path(aln)
        taxa = set(read_fasta_names(aln_path))

        gene_to_aln[gene] = aln_path
        gene_to_taxa[gene] = taxa
        gene_to_stats[gene] = {
            "tree": tree,
            "n_taxa_alignment": int(n_taxa_alignment),
            "aln_len": int(aln_len),
            "n_taxa_tree": int(n_taxa_tree),
            "n_taxa_species_tree_overlap": int(n_taxa_species_tree_overlap),
            "n_taxa_aln_tree_shared": int(n_taxa_aln_tree_shared),
            "fraction_species_tree_covered": float(fraction_species_tree_covered),
            "fraction_alignment_in_tree": float(fraction_alignment_in_tree),
        }

ranked_genes = []
with open(ranked_file) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        tree_name = parts[0]
        gene = Path(tree_name).name.replace(".tre", "")
        if gene in gene_to_aln:
            ranked_genes.append((gene, line))

if not ranked_genes:
    raise SystemExit("[ERR] No ranked loci could be matched back to candidate alignments.")

all_ranked_taxa = set()
for gene, _ in ranked_genes:
    all_ranked_taxa |= gene_to_taxa[gene]

missing_even_with_all = sorted(full_taxa - all_ranked_taxa)
with open(uncovered_tsv, "w") as out:
    out.write("taxon\n")
    for t in missing_even_with_all:
        out.write(f"{t}\n")

if missing_even_with_all:
    raise SystemExit(
        f"[ERR] Even all ranked candidate loci together do not cover the full species tree. "
        f"Missing taxa: {len(missing_even_with_all)}. See {uncovered_tsv}"
    )

selected = []
selected_genes = set()
covered = set()

with open(coverage_tsv, "w") as out_cov:
    out_cov.write("selection_order\tgene\tnew_taxa_added\ttotal_taxa_covered\ttotal_taxa_full\n")

    remaining = ranked_genes[:]
    step = 0

    while remaining and covered != full_taxa:
        best_idx = None
        best_gain = -1
        best_tiebreak = None

        for i, (gene, raw_line) in enumerate(remaining):
            gain = len(gene_to_taxa[gene] - covered)
            tiebreak = (
                gene_to_stats[gene]["n_taxa_species_tree_overlap"],
                gene_to_stats[gene]["n_taxa_aln_tree_shared"],
                gene_to_stats[gene]["n_taxa_alignment"],
                gene_to_stats[gene]["aln_len"],
            )
            if (gain > best_gain) or (gain == best_gain and (best_tiebreak is None or tiebreak > best_tiebreak)):
                best_gain = gain
                best_idx = i
                best_tiebreak = tiebreak

        if best_idx is None or best_gain <= 0:
            break

        gene, raw_line = remaining.pop(best_idx)
        selected.append((gene, gene_to_aln[gene], raw_line))
        selected_genes.add(gene)

        new_taxa = gene_to_taxa[gene] - covered
        covered |= gene_to_taxa[gene]
        step += 1

        out_cov.write(f"{step}\t{gene}\t{len(new_taxa)}\t{len(covered)}\t{len(full_taxa)}\n")

    for gene, raw_line in ranked_genes:
        if gene in selected_genes:
            continue
        selected.append((gene, gene_to_aln[gene], raw_line))
        selected_genes.add(gene)
        step += 1
        out_cov.write(f"{step}\t{gene}\t0\t{len(covered)}\t{len(full_taxa)}\n")

with open(selected_tsv, "w") as out:
    out.write("selection_order\tgene\talignment\tranked_line\n")
    for i, (gene, aln, raw) in enumerate(selected, start=1):
        out.write(f"{i}\t{gene}\t{aln}\t{raw}\n")

with open(selected_list, "w") as out:
    for gene, aln, raw in selected:
        out.write(f"{aln}\n")

uncovered = sorted(full_taxa - covered)
with open(uncovered_tsv, "w") as out:
    out.write("taxon\n")
    for t in uncovered:
        out.write(f"{t}\n")

if uncovered:
    raise SystemExit(
        f"[ERR] Master selected locus order still does not cover the full species tree. "
        f"Uncovered taxa: {len(uncovered)}. See {uncovered_tsv}"
    )

print(f"[INFO] Full taxa in species tree:        {len(full_taxa)}")
print(f"[INFO] Taxa covered after master pass:  {len(covered)}")
print(f"[INFO] Total loci in master order:      {len(selected)}")
print(f"[INFO] Coverage progress table:         {coverage_tsv}")
print(f"[INFO] Master selected loci table:      {selected_tsv}")
print(f"[INFO] Master selected loci list:       {selected_list}")
PY

N_MASTER_SELECTED=$(( $(wc -l < "$DATESEL_MASTER_SELECTED_TSV") - 1 ))
echo "[INFO] Master ordered locus list size: $N_MASTER_SELECTED" | tee -a "$RUNLOG"

# -------------------- 8.6 Build reusable subset folders ----------------

echo "[INFO] STEP 8.6: building reusable top-N subset folders" | tee -a "$RUNLOG"

python3 - <<'PY'
from pathlib import Path
import os
import csv
import shutil

master_tsv = Path(os.environ["DATESEL_MASTER_SELECTED_TSV"])
subset_root = Path(os.environ["DATESEL_SUBSETS"])
subset_sizes_csv = os.environ["DATE_SUBSET_SIZES_CSV"]
summary_tsv = Path(os.environ["DATESEL_SUBSET_SUMMARY_TSV"])
mct_n = int(os.environ["MCT_N_LOCI"])
lsd2_n = int(os.environ["LSD2_N_LOCI"])

sizes = []
for x in subset_sizes_csv.split(","):
    x = x.strip()
    if not x:
        continue
    try:
        n = int(x)
    except Exception:
        raise SystemExit(f"[ERR] Invalid subset size in DATE_SUBSET_SIZES_CSV: {x}")
    if n <= 0:
        raise SystemExit(f"[ERR] Subset size must be > 0: {n}")
    sizes.append(n)

for extra in (mct_n, lsd2_n):
    if extra not in sizes:
        sizes.append(extra)

sizes = sorted(set(sizes))

rows = []
with open(master_tsv, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        row["selection_order"] = int(row["selection_order"])
        rows.append(row)

if not rows:
    raise SystemExit("[ERR] Master selected locus table is empty.")

n_available = len(rows)

with open(summary_tsv, "w") as out_sum:
    out_sum.write("subset_name\tn_requested\tn_written\tsubset_dir\talignments_dir\tloci_list\tselected_tsv\n")

    for n in sizes:
        n_use = min(n, n_available)
        subset_name = f"top_{n}"
        subset_dir = subset_root / subset_name
        align_dir = subset_dir / "alignments"
        loci_list = subset_dir / "loci_list.txt"
        subset_tsv = subset_dir / "selected_loci.tsv"

        if subset_dir.exists():
            shutil.rmtree(subset_dir)
        align_dir.mkdir(parents=True, exist_ok=True)

        chosen = rows[:n_use]

        with open(loci_list, "w") as out_list, open(subset_tsv, "w") as out_tsv:
            out_tsv.write("selection_order\tgene\talignment\tranked_line\n")
            for row in chosen:
                src = Path(row["alignment"])
                dst = align_dir / src.name
                shutil.copyfile(src, dst)
                out_list.write(f"{dst}\n")
                out_tsv.write(
                    f"{row['selection_order']}\t{row['gene']}\t{dst}\t{row['ranked_line']}\n"
                )

        out_sum.write(
            f"{subset_name}\t{n}\t{n_use}\t{subset_dir}\t{align_dir}\t{loci_list}\t{subset_tsv}\n"
        )

print(f"[INFO] Subset summary written to: {summary_tsv}")
print(f"[INFO] Available loci in master order: {n_available}")
print(f"[INFO] Requested subset sizes: {','.join(map(str, sizes))}")
PY

# -------------------- 8.7 Final reporting ------------------------------

echo "[INFO] STEP 8 completed." | tee -a "$RUNLOG"
echo "[INFO] Main Step 8 outputs:" | tee -a "$RUNLOG"
echo "[INFO]   Clean ASTRAL topology:         $DATESEL_ASTRAL_TREE_CLEAN" | tee -a "$RUNLOG"
echo "[INFO]   Ingroup species tree:          $DATESEL_INGROUP_SPECIES_TREE" | tee -a "$RUNLOG"
echo "[INFO]   Candidate loci table:          $DATESEL_CANDIDATE_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Candidate failure table:       $DATESEL_FAIL_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Ranked loci file:              $DATESEL_RANKED_FILE" | tee -a "$RUNLOG"
echo "[INFO]   Master selected loci table:    $DATESEL_MASTER_SELECTED_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Master selected loci list:     $DATESEL_MASTER_SELECTED_LIST" | tee -a "$RUNLOG"
echo "[INFO]   Coverage progress table:       $DATESEL_COVERAGE_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Uncovered taxa table:          $DATESEL_UNCOVERED_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Subset summary table:          $DATESEL_SUBSET_SUMMARY_TSV" | tee -a "$RUNLOG"
echo "[INFO]   Reusable subset root:          $DATESEL_SUBSETS" | tee -a "$RUNLOG"
echo "[INFO]   Default MCMCtree subset:       $DATESEL_SUBSETS/top_${MCT_N_LOCI}" | tee -a "$RUNLOG"
echo "[INFO]   Default LSD2 subset:           $DATESEL_SUBSETS/top_${LSD2_N_LOCI}" | tee -a "$RUNLOG"




# =====================================================
# ----------------------- STEP 10 ---------------------
# LSD2 DIVERGENCE TIME ESTIMATION ON FIXED ASTRAL TOPOLOGY
# USING A SHARED STEP-8 LOCUS SUBSET
# =====================================================

if [[ "${RUN_LSD2:-yes}" != "yes" ]]; then
  echo "[INFO] STEP 10 skipped because RUN_LSD2=${RUN_LSD2:-no}" | tee -a "$RUNLOG"
else

echo "[INFO] STEP 10: LSD2 dating on fixed ASTRAL topology using a selected Step-8 locus subset" | tee -a "$RUNLOG"

# -------------------- USER SETTINGS -------------------

# Number of loci to use from shared Step 8 subset preparation
LSD2_N_LOCI="${LSD2_N_LOCI:-200}"

# Folder tag for this specific LSD2 run
LSD2_RUN_TAG="lsd2_top${LSD2_N_LOCI}_fixedtopo"

# Fixed DNA model for partitioned branch-length estimation
LSD2_MODEL="${LSD2_MODEL:-GTR+F+R4}"

# Confidence interval replicates for LSD2
LSD2_DATE_CI="${LSD2_DATE_CI:-100}"

# Relaxed-clock SD for CI estimation in LSD2
LSD2_CLOCK_SD="${LSD2_CLOCK_SD:-0.2}"

# Optional outlier threshold for LSD2 (empty = off)
LSD2_DATE_OUTLIER="${LSD2_DATE_OUTLIER:-}"

# Minimum internal branch length in dated tree, in same units as dates.
# Here dates are in Ma, so 1e-3 = 0.001 Ma = 1000 years.
LSD2_MIN_INTERNAL_BLEN="${LSD2_MIN_INTERNAL_BLEN:-1e-2}"

# Rounding factor for minimum time-scaled branch length.
# With Ma units, -R 1000 gives resolution 1/1000 Ma = 0.001 Ma.
LSD2_ROUND_TIME="${LSD2_ROUND_TIME:-100}"

# Optional branch-collapse threshold for LSD2
# Examples:
#   empty = use LSD2 default
#   0     = collapse only null branches
#   -1    = do not collapse branches
LSD2_COLLAPSE_THRESHOLD="${LSD2_COLLAPSE_THRESHOLD:--1}"

# Threads for IQ-TREE branch-length estimation
LSD2_THREADS="${LSD2_THREADS:-20}"

# Fixed root age for the complete tree
# IMPORTANT: this must be NEGATIVE when tips are set to 0 with -z 0
LSD2_ROOT_AGE="${LSD2_ROOT_AGE:--90}"

# Standalone LSD2 binary
LSD2_BIN="${LSD2_BIN:-$HOME/lsd2/src/lsd2}"

echo "[INFO] LSD2 min internal branch (-u): $LSD2_MIN_INTERNAL_BLEN" | tee -a "$RUNLOG"
echo "[INFO] LSD2 round time (-R):           $LSD2_ROUND_TIME" | tee -a "$RUNLOG"

# -------------------- PATHS ---------------------------

LSD2_DIR="$CALIB_DIR/$LSD2_RUN_TAG"
LSD2_TREE_DIR="$LSD2_DIR/01_treefiles"
LSD2_ALIGN_SRC_DIR="$LSD2_DIR/02_selected_alignment_source"
LSD2_ALIGN_DIR="$LSD2_DIR/03_selected_alignments"
LSD2_BL_DIR="$LSD2_DIR/04_fixedtopo_branchlengths"
LSD2_DATE_DIR="$LSD2_DIR/05_lsd2_dating"

mkdir -p \
  "$LSD2_TREE_DIR" \
  "$LSD2_ALIGN_SRC_DIR" \
  "$LSD2_ALIGN_DIR" \
  "$LSD2_BL_DIR" \
  "$LSD2_DATE_DIR"

# Reuse Step-8 outputs
DATESEL_DIR="$CALIB_DIR/00_dating_locus_selection"
DATESEL_SUBSETS="$DATESEL_DIR/05_selected_subsets"

LSD2_SUBSET_DIR="$DATESEL_SUBSETS/top_${LSD2_N_LOCI}"
LSD2_SUBSET_ALIGN_DIR="$LSD2_SUBSET_DIR/alignments"
LSD2_SUBSET_LIST="$LSD2_SUBSET_DIR/loci_list.txt"
LSD2_SUBSET_TSV="$LSD2_SUBSET_DIR/selected_loci.tsv"

# Main tree inputs
LSD2_ASTRAL_TREE_RAW="$ASTRAL_DIR/${project_name}_BrassiToL_coalescent.tree"
LSD2_ASTRAL_TREE_STRIPPED="$LSD2_TREE_DIR/${project_name}_astral_stripped_topology.nwk"

# Outgroup file required by standalone LSD2
LSD2_OUTGROUP_FILE="$LSD2_TREE_DIR/lsd2_outgroups.txt"

# Main outputs
LSD2_BL_PREFIX="$LSD2_BL_DIR/${project_name}_lsd2_top${LSD2_N_LOCI}_fixedtopo_brlen"
LSD2_DATE_PREFIX="$LSD2_DATE_DIR/${project_name}_lsd2_top${LSD2_N_LOCI}_dated"
LSD2_ALIGNMENT_MANIFEST="$LSD2_ALIGN_SRC_DIR/copied_alignment_manifest.tsv"

# -------------------- SANITY CHECKS -------------------

if [[ ! -x "$LSD2_BIN" ]]; then
  echo "[ERR] LSD2 binary not found or not executable: $LSD2_BIN" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -s "$LSD2_ASTRAL_TREE_RAW" ]]; then
  echo "[ERR] ASTRAL species tree not found or empty: $LSD2_ASTRAL_TREE_RAW" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -d "$LSD2_SUBSET_ALIGN_DIR" ]]; then
  echo "[ERR] Requested Step-8 subset alignment directory not found: $LSD2_SUBSET_ALIGN_DIR" | tee -a "$RUNLOG"
  echo "[ERR] Check whether Step 8 created top_${LSD2_N_LOCI} successfully." | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -s "$LSD2_SUBSET_LIST" ]]; then
  echo "[ERR] Requested Step-8 subset loci list not found or empty: $LSD2_SUBSET_LIST" | tee -a "$RUNLOG"
  exit 1
fi

if [[ ! -s "$LSD2_SUBSET_TSV" ]]; then
  echo "[ERR] Requested Step-8 subset table not found or empty: $LSD2_SUBSET_TSV" | tee -a "$RUNLOG"
  exit 1
fi

N_LSD2_SUBSET=$(find "$LSD2_SUBSET_ALIGN_DIR" -maxdepth 1 -type f -name "*.fasta" | wc -l)
if (( N_LSD2_SUBSET == 0 )); then
  echo "[ERR] No FASTA alignments found in selected Step-8 subset: $LSD2_SUBSET_ALIGN_DIR" | tee -a "$RUNLOG"
  exit 1
fi

if [[ -z "${sample_outgroup:-}" ]]; then
  echo "[ERR] sample_outgroup is empty; LSD2 outgroup file cannot be created." | tee -a "$RUNLOG"
  exit 1
fi

if [[ -z "${LSD2_ROOT_AGE:-}" ]]; then
  echo "[ERR] LSD2_ROOT_AGE is empty; a fixed root age is required for this run." | tee -a "$RUNLOG"
  exit 1
fi

echo "[INFO] LSD2 run folder:              $LSD2_DIR" | tee -a "$RUNLOG"
echo "[INFO] Step-8 subset folder:         $LSD2_SUBSET_DIR" | tee -a "$RUNLOG"
echo "[INFO] Step-8 subset alignments:     $LSD2_SUBSET_ALIGN_DIR" | tee -a "$RUNLOG"
echo "[INFO] Requested top-N loci:         $LSD2_N_LOCI" | tee -a "$RUNLOG"
echo "[INFO] Loci available in subset:     $N_LSD2_SUBSET" | tee -a "$RUNLOG"
echo "[INFO] Outgroup tip:                 $sample_outgroup" | tee -a "$RUNLOG"
echo "[INFO] Fixed root age for LSD2:      $LSD2_ROOT_AGE" | tee -a "$RUNLOG"
echo "[INFO] IQ-TREE model for LSD2:       $LSD2_MODEL" | tee -a "$RUNLOG"
echo "[INFO] LSD2 CI replicates:           $LSD2_DATE_CI" | tee -a "$RUNLOG"
echo "[INFO] LSD2 clock SD:                $LSD2_CLOCK_SD" | tee -a "$RUNLOG"
echo "[INFO] LSD2 threads:                 $LSD2_THREADS" | tee -a "$RUNLOG"
echo "[INFO] LSD2 binary:                  $LSD2_BIN" | tee -a "$RUNLOG"

# -------------------- 10.1 Strip ASTRAL topology -----------------------

echo "[INFO] STEP 10.1: stripping branch lengths and node annotations from the ASTRAL tree, while keeping all taxa including the outgroup" | tee -a "$RUNLOG"

export LSD2_ASTRAL_TREE_RAW LSD2_ASTRAL_TREE_STRIPPED

python3 - <<'PY'
from pathlib import Path
import os
import re

inp = Path(os.environ["LSD2_ASTRAL_TREE_RAW"])
outp = Path(os.environ["LSD2_ASTRAL_TREE_STRIPPED"])

if not inp.exists():
    raise SystemExit(f"[ERR] Input ASTRAL tree not found: {inp}")

txt = inp.read_text(errors="replace")

# Remove square-bracket annotations/comments
txt = re.sub(r"\[[^\]]*\]", "", txt)

# Remove single quotes
txt = txt.replace("'", "")

# Remove whitespace/newlines
txt = re.sub(r"\s+", "", txt)

if ";" not in txt:
    raise SystemExit("[ERR] No semicolon found in tree file.")

# Keep only first tree if multiple are present
txt = txt.split(";", 1)[0] + ";"

# Remove internal node labels/support labels immediately after ')'
txt = re.sub(r"\)([^:(),;\[\]]+)(?=[:),;])", ")", txt)

# Remove all branch lengths
txt = re.sub(r":[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "", txt)

if not txt.endswith(";"):
    raise SystemExit("[ERR] Stripped topology does not end with ';'.")

outp.write_text(txt + "\n")
print(f"[INFO] Stripped ASTRAL topology written to: {outp}")
PY

if [[ ! -s "$LSD2_ASTRAL_TREE_STRIPPED" ]]; then
  echo "[ERR] Stripped ASTRAL topology was not created: $LSD2_ASTRAL_TREE_STRIPPED" | tee -a "$RUNLOG"
  exit 1
fi

# -------------------- 10.2 Write LSD2 outgroup file --------------------

echo "[INFO] STEP 10.2: writing LSD2 outgroup file" | tee -a "$RUNLOG"

{
  echo "1"
  echo "$sample_outgroup"
} > "$LSD2_OUTGROUP_FILE"

if [[ ! -s "$LSD2_OUTGROUP_FILE" ]]; then
  echo "[ERR] Failed to create LSD2 outgroup file: $LSD2_OUTGROUP_FILE" | tee -a "$RUNLOG"
  exit 1
fi

echo "[INFO] LSD2 outgroup file written to: $LSD2_OUTGROUP_FILE" | tee -a "$RUNLOG"

# -------------------- 10.3 Copy selected loci locally ------------------

echo "[INFO] STEP 10.3: copying selected Step-8 subset alignments into the LSD2 run folder" | tee -a "$RUNLOG"

rm -rf "$LSD2_ALIGN_DIR"
mkdir -p "$LSD2_ALIGN_DIR"

printf "source_alignment\tdestination_alignment\n" > "$LSD2_ALIGNMENT_MANIFEST"

find "$LSD2_SUBSET_ALIGN_DIR" -maxdepth 1 -type f -name "*.fasta" | sort | while read -r fp; do
  [[ -z "$fp" ]] && continue
  base="$(basename "$fp")"
  cp "$fp" "$LSD2_ALIGN_DIR/$base"
  printf "%s\t%s\n" "$fp" "$LSD2_ALIGN_DIR/$base" >> "$LSD2_ALIGNMENT_MANIFEST"
done

N_LSD2_ALN=$(find "$LSD2_ALIGN_DIR" -maxdepth 1 -type f -name "*.fasta" | wc -l)
echo "[INFO] Selected loci copied into LSD2 alignment folder: $N_LSD2_ALN" | tee -a "$RUNLOG"
echo "[INFO] Alignment manifest written to: $LSD2_ALIGNMENT_MANIFEST" | tee -a "$RUNLOG"

if (( N_LSD2_ALN == 0 )); then
  echo "[ERR] No alignments were copied into LSD2 alignment folder: $LSD2_ALIGN_DIR" | tee -a "$RUNLOG"
  exit 1
fi

# -------------------- 10.4 Estimate branch lengths ---------------------

echo "[INFO] STEP 10.4: estimating fixed-topology substitution branch lengths with IQ-TREE on the stripped ASTRAL topology" | tee -a "$RUNLOG"

"$IQTREE_BIN" \
  -p "$LSD2_ALIGN_DIR" \
  -te "$LSD2_ASTRAL_TREE_STRIPPED" \
  -m "$LSD2_MODEL" \
  -pre "$LSD2_BL_PREFIX" \
  -nt "$LSD2_THREADS" \
  &> "${LSD2_BL_PREFIX}.log"

if [[ ! -s "${LSD2_BL_PREFIX}.treefile" ]]; then
  echo "[ERR] Fixed-topology branch-length tree not found: ${LSD2_BL_PREFIX}.treefile" | tee -a "$RUNLOG"
  echo "[ERR] Check log: ${LSD2_BL_PREFIX}.log" | tee -a "$RUNLOG"
  exit 1
fi

echo "[INFO] Fixed-topology branch-length tree written to: ${LSD2_BL_PREFIX}.treefile" | tee -a "$RUNLOG"

# -------------------- 10.5 Calculate total sequence length -------------

echo "[INFO] STEP 10.5: calculating total aligned sequence length for standalone LSD2" | tee -a "$RUNLOG"

export LSD2_ALIGN_DIR
LSD2_SEQ_LEN=$(
python3 - <<'PY'
from pathlib import Path
import os

aln_dir = Path(os.environ["LSD2_ALIGN_DIR"])
total = 0

for fp in sorted(aln_dir.glob("*.fasta")):
    lengths = set()
    chunks = []

    with open(fp) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    lengths.add(len("".join(chunks)))
                chunks = []
            else:
                chunks.append(line.replace(" ", "").replace("\t", ""))
        if chunks:
            lengths.add(len("".join(chunks)))

    if not lengths:
        continue
    if len(lengths) != 1:
        raise SystemExit(f"[ERR] Alignment is not rectangular: {fp}")
    total += next(iter(lengths))

print(total)
PY
)

if [[ -z "$LSD2_SEQ_LEN" || "$LSD2_SEQ_LEN" -le 0 ]]; then
  echo "[ERR] Could not compute a valid total sequence length for LSD2." | tee -a "$RUNLOG"
  exit 1
fi

echo "[INFO] Total alignment length for -s: $LSD2_SEQ_LEN" | tee -a "$RUNLOG"

# -------------------- 10.6 Run standalone LSD2 dating ------------------

echo "[INFO] STEP 10.6: running standalone LSD2 on the fixed-topology branch-length tree" | tee -a "$RUNLOG"

echo "[INFO] LSD2 input tree:               ${LSD2_BL_PREFIX}.treefile" | tee -a "$RUNLOG"
echo "[INFO] LSD2 output prefix:            $LSD2_DATE_PREFIX" | tee -a "$RUNLOG"
echo "[INFO] LSD2 root age (-a):            $LSD2_ROOT_AGE" | tee -a "$RUNLOG"
echo "[INFO] LSD2 tips date (-z):           0" | tee -a "$RUNLOG"
echo "[INFO] LSD2 outgroup file (-g):       $LSD2_OUTGROUP_FILE" | tee -a "$RUNLOG"

LSD2_EXTRA_ARGS=()

if [[ -n "${LSD2_DATE_OUTLIER:-}" ]]; then
  LSD2_EXTRA_ARGS+=( -e "$LSD2_DATE_OUTLIER" )
fi

if [[ -n "${LSD2_COLLAPSE_THRESHOLD:-}" ]]; then
  LSD2_EXTRA_ARGS+=( -l "$LSD2_COLLAPSE_THRESHOLD" )
fi

if [[ -n "${LSD2_MIN_INTERNAL_BLEN:-}" ]]; then
  LSD2_EXTRA_ARGS+=( -u "$LSD2_MIN_INTERNAL_BLEN" )
fi

if [[ -n "${LSD2_ROUND_TIME:-}" ]]; then
  LSD2_EXTRA_ARGS+=( -R "$LSD2_ROUND_TIME" )
fi

if [[ "${LSD2_DATE_CI:-0}" =~ ^[0-9]+$ ]] && (( LSD2_DATE_CI > 0 )); then
  echo "[INFO] Running LSD2 with confidence intervals (n=${LSD2_DATE_CI}, q=${LSD2_CLOCK_SD})" | tee -a "$RUNLOG"

  "$LSD2_BIN" \
    -i "${LSD2_BL_PREFIX}.treefile" \
    -o "$LSD2_DATE_PREFIX" \
    -s "$LSD2_SEQ_LEN" \
    -a "$LSD2_ROOT_AGE" \
    -g "$LSD2_OUTGROUP_FILE" \
    -z 0 \
    -f "$LSD2_DATE_CI" \
    -q "$LSD2_CLOCK_SD" \
    "${LSD2_EXTRA_ARGS[@]}" \
    > "${LSD2_DATE_PREFIX}.stdout.log" \
    2> "${LSD2_DATE_PREFIX}.stderr.log"
else
  echo "[INFO] Running LSD2 without confidence intervals" | tee -a "$RUNLOG"

  "$LSD2_BIN" \
    -i "${LSD2_BL_PREFIX}.treefile" \
    -o "$LSD2_DATE_PREFIX" \
    -s "$LSD2_SEQ_LEN" \
    -a "$LSD2_ROOT_AGE" \
    -g "$LSD2_OUTGROUP_FILE" \
    -z 0 \
    "${LSD2_EXTRA_ARGS[@]}" \
    > "${LSD2_DATE_PREFIX}.stdout.log" \
    2> "${LSD2_DATE_PREFIX}.stderr.log"
fi

# -------------------- 10.7 LSD2 sensitivity analysis -------------------

echo "[INFO] STEP 10.7: LSD2 sensitivity analysis" | tee -a "$RUNLOG"

LSD2_DATE_DIR="$(dirname "$LSD2_DATE_PREFIX")"
LSD2_PARENT_DIR="$(dirname "$LSD2_DATE_DIR")"
LSD2_MAIN_LABEL="$(basename "$LSD2_DATE_PREFIX")"

LSD2_SENS_DIR="${LSD2_PARENT_DIR}/06_lsd2_sensitivity_analysis"
mkdir -p "$LSD2_SENS_DIR"

echo "[INFO] LSD2 sensitivity output dir: $LSD2_SENS_DIR" | tee -a "$RUNLOG"

score_lsd2_tree() {
  local label="$1"
  local collapse_threshold="$2"
  local min_internal_blen="$3"
  local round_time="$4"
  local date_ci="$5"
  local clock_sd="$6"
  local prefix="$7"
  local score_file="$8"
  
  if [[ ! -s "${prefix}.date.nexus" ]]; then
    echo "[WARN] Cannot score tree; missing: ${prefix}.date.nexus" | tee -a "$RUNLOG"
    return 0
  fi
  
  Rscript --vanilla - "$label" "$collapse_threshold" "$min_internal_blen" "$round_time" "$date_ci" "$clock_sd" "$prefix" "$score_file" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

label  <- args[1]
l_arg  <- args[2]
u_arg  <- args[3]
R_arg  <- args[4]
f_arg  <- args[5]
q_arg  <- args[6]
prefix <- args[7]
out    <- args[8]

suppressPackageStartupMessages(library(ape))

tree_file <- paste0(prefix, ".date.nexus")
tr <- read.nexus(tree_file)
if (inherits(tr, "multiPhylo")) tr <- tr[[1]]

el <- tr$edge.length

tol1 <- 1e-12
tol2 <- 1e-10
tol3 <- 1e-8
tol4 <- 1e-6

tip_depths <- node.depth.edgelength(tr)[seq_along(tr$tip.label)]
tree_height <- max(tip_depths, na.rm = TRUE)
tip_depth_span <- diff(range(tip_depths, na.rm = TRUE))

internal_edges <- tr$edge[, 2] > length(tr$tip.label)
terminal_edges <- tr$edge[, 2] <= length(tr$tip.label)

res <- data.frame(
  label = label,
  collapse_threshold = ifelse(nzchar(l_arg), l_arg, NA),
  min_internal_blen = ifelse(nzchar(u_arg), u_arg, NA),
  round_time = ifelse(nzchar(R_arg), R_arg, NA),
  date_ci = f_arg,
  clock_sd = ifelse(nzchar(q_arg), q_arg, NA),
  n_tips = length(tr$tip.label),
  n_edges = length(el),
  n_internal_edges = sum(internal_edges),
  n_terminal_edges = sum(terminal_edges),
  n_zero_exact = sum(el == 0, na.rm = TRUE),
  frac_zero_exact = mean(el == 0, na.rm = TRUE),
  n_near_zero_1e12 = sum(el <= tol1, na.rm = TRUE),
  frac_near_zero_1e12 = mean(el <= tol1, na.rm = TRUE),
  n_near_zero_1e10 = sum(el <= tol2, na.rm = TRUE),
  frac_near_zero_1e10 = mean(el <= tol2, na.rm = TRUE),
  n_near_zero_1e8 = sum(el <= tol3, na.rm = TRUE),
  frac_near_zero_1e8 = mean(el <= tol3, na.rm = TRUE),
  n_near_zero_1e6 = sum(el <= tol4, na.rm = TRUE),
  frac_near_zero_1e6 = mean(el <= tol4, na.rm = TRUE),
  n_zero_internal_exact = sum(el[internal_edges] == 0, na.rm = TRUE),
  n_zero_terminal_exact = sum(el[terminal_edges] == 0, na.rm = TRUE),
  min_edge = min(el, na.rm = TRUE),
  q001_edge = as.numeric(quantile(el, 0.001, na.rm = TRUE, names = FALSE)),
  q005_edge = as.numeric(quantile(el, 0.005, na.rm = TRUE, names = FALSE)),
  q01_edge = as.numeric(quantile(el, 0.01, na.rm = TRUE, names = FALSE)),
  median_edge = median(el, na.rm = TRUE),
  mean_edge = mean(el, na.rm = TRUE),
  max_edge = max(el, na.rm = TRUE),
  total_tree_length = sum(el, na.rm = TRUE),
  tree_height = tree_height,
  tip_depth_span = tip_depth_span,
  is_ultrametric_1e8 = is.ultrametric(tr, tol = 1e-8),
  is_ultrametric_1e6 = is.ultrametric(tr, tol = 1e-6)
)

write.table(
  res,
  file = out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
RSCRIPT
}

LSD2_SENS_CONFIGS=(
  "base_noCI|||0|"
  "u1e-6_R1000000_noCI||1e-6|1000000|0|"
  "u1e-4_R10000_noCI||1e-4|10000|0|"
  "u1e-3_R1000_noCI||1e-3|1000|0|"
  "l-1_u1e-3_R1000_noCI|-1|1e-3|1000|0|"
)

SENS_RUN_SUMMARY="${LSD2_SENS_DIR}/lsd2_sensitivity_runs.tsv"
SENS_SCORE_SUMMARY="${LSD2_SENS_DIR}/lsd2_sensitivity_tree_scores.tsv"

echo -e "label\tcollapse_threshold\tmin_internal_blen\tround_time\tdate_ci\tclock_sd\tstatus\tprefix\tresult_file\tnexus_file\tnwk_file\tstdout_log\tstderr_log" > "$SENS_RUN_SUMMARY"

for cfg in "${LSD2_SENS_CONFIGS[@]}"; do
  IFS='|' read -r sens_label sens_l sens_u sens_R sens_f sens_q <<< "$cfg"
  
  sens_prefix="${LSD2_SENS_DIR}/${LSD2_MAIN_LABEL}_${sens_label}"
  sens_args=()
  
if [[ -n "$sens_l" ]]; then
  sens_args+=( -l "$sens_l" )
fi

if [[ -n "$sens_u" ]]; then
  sens_args+=( -u "$sens_u" )
fi

if [[ -n "$sens_R" ]]; then
  sens_args+=( -R "$sens_R" )
fi
  
  echo "[INFO] LSD2 sensitivity run: $sens_label" | tee -a "$RUNLOG"
  echo "[INFO]   output prefix: ${sens_prefix}" | tee -a "$RUNLOG"
  echo "[INFO]   -l: ${sens_l:-<none>}" | tee -a "$RUNLOG"
  echo "[INFO]   -u: ${sens_u:-<none>}" | tee -a "$RUNLOG"
  echo "[INFO]   -R: ${sens_R:-<none>}" | tee -a "$RUNLOG"
  echo "[INFO]   -f: ${sens_f:-0}" | tee -a "$RUNLOG"
  echo "[INFO]   -q: ${sens_q:-<none>}" | tee -a "$RUNLOG"
  
  if [[ "$sens_f" =~ ^[0-9]+$ ]] && (( sens_f > 0 )); then
    q_use="${sens_q:-$LSD2_CLOCK_SD}"
    
    "$LSD2_BIN" \
      -i "${LSD2_BL_PREFIX}.treefile" \
      -o "$sens_prefix" \
      -s "$LSD2_SEQ_LEN" \
      -a "$LSD2_ROOT_AGE" \
      -g "$LSD2_OUTGROUP_FILE" \
      -z 0 \
      -f "$sens_f" \
      -q "$q_use" \
      "${sens_args[@]}" \
      > "${sens_prefix}.stdout.log" \
      2> "${sens_prefix}.stderr.log" || true
  else
    "$LSD2_BIN" \
      -i "${LSD2_BL_PREFIX}.treefile" \
      -o "$sens_prefix" \
      -s "$LSD2_SEQ_LEN" \
      -a "$LSD2_ROOT_AGE" \
      -g "$LSD2_OUTGROUP_FILE" \
      -z 0 \
      "${sens_args[@]}" \
      > "${sens_prefix}.stdout.log" \
      2> "${sens_prefix}.stderr.log" || true
  fi
  
  status="FAIL"
  if [[ -s "${sens_prefix}.result" ]]; then
    status="OK_RESULT"
  elif [[ -s "${sens_prefix}.date.nexus" || -s "${sens_prefix}.nwk" ]]; then
    status="OK_TREE"
  fi
  
  if [[ "$status" == "FAIL" ]]; then
    echo "[WARN] Sensitivity run produced no usable result/tree: $sens_label" | tee -a "$RUNLOG"
  else
    echo "[INFO] Sensitivity run status: $sens_label = $status" | tee -a "$RUNLOG"
  fi
  
  echo -e "${sens_label}\t${sens_l:-NA}\t${sens_u:-NA}\t${sens_R:-NA}\t${sens_f:-0}\t${sens_q:-NA}\t${status}\t${sens_prefix}\t${sens_prefix}.result\t${sens_prefix}.date.nexus\t${sens_prefix}.nwk\t${sens_prefix}.stdout.log\t${sens_prefix}.stderr.log" >> "$SENS_RUN_SUMMARY"
  
  score_file="${sens_prefix}.tree_score.tsv"
  score_lsd2_tree "$sens_label" "$sens_l" "$sens_u" "$sens_R" "$sens_f" "$sens_q" "$sens_prefix" "$score_file"
done

rm -f "$SENS_SCORE_SUMMARY"

first=1
for f in "${LSD2_SENS_DIR}"/*.tree_score.tsv; do
  [[ -s "$f" ]] || continue
  
  if [[ "$first" -eq 1 ]]; then
    cat "$f" > "$SENS_SCORE_SUMMARY"
    first=0
  else
    tail -n +2 "$f" >> "$SENS_SCORE_SUMMARY"
  fi
done

echo "[INFO] LSD2 sensitivity run table written to:   $SENS_RUN_SUMMARY" | tee -a "$RUNLOG"
echo "[INFO] LSD2 sensitivity score table written to: $SENS_SCORE_SUMMARY" | tee -a "$RUNLOG"


# -------------------- 10.8 Check and report LSD2 outputs ---------------------

echo "[INFO] STEP 10.8: checking and reporting LSD2 outputs" | tee -a "$RUNLOG"

if [[ -s "${LSD2_DATE_PREFIX}.result" ]]; then
  echo "[INFO] Main LSD2 result file written to: ${LSD2_DATE_PREFIX}.result" | tee -a "$RUNLOG"
else
  echo "[WARN] Main LSD2 .result file not found or empty: ${LSD2_DATE_PREFIX}.result" | tee -a "$RUNLOG"
fi

if [[ -s "${LSD2_DATE_PREFIX}.date.nexus" ]]; then
  echo "[INFO] Main LSD2 dated tree written to: ${LSD2_DATE_PREFIX}.date.nexus" | tee -a "$RUNLOG"
else
  echo "[ERR] Main LSD2 .date.nexus file not found: ${LSD2_DATE_PREFIX}.date.nexus" | tee -a "$RUNLOG"
  echo "[ERR] Check logs:" | tee -a "$RUNLOG"
  echo "[ERR]   ${LSD2_DATE_PREFIX}.stdout.log" | tee -a "$RUNLOG"
  echo "[ERR]   ${LSD2_DATE_PREFIX}.stderr.log" | tee -a "$RUNLOG"
  exit 1
fi

if [[ -s "${LSD2_DATE_PREFIX}.nwk" ]]; then
  echo "[INFO] Main LSD2 dated Newick written to: ${LSD2_DATE_PREFIX}.nwk" | tee -a "$RUNLOG"
else
  echo "[WARN] Main LSD2 .nwk file not found: ${LSD2_DATE_PREFIX}.nwk" | tee -a "$RUNLOG"
fi

MAIN_TREE_SCORE="${LSD2_SENS_DIR}/${LSD2_MAIN_LABEL}_main.tree_score.tsv"

score_lsd2_tree \
  "main" \
  "${LSD2_COLLAPSE_THRESHOLD:-}" \
  "${LSD2_MIN_INTERNAL_BLEN:-}" \
  "${LSD2_ROUND_TIME:-}" \
  "${LSD2_DATE_CI:-0}" \
  "${LSD2_CLOCK_SD:-}" \
  "$LSD2_DATE_PREFIX" \
  "$MAIN_TREE_SCORE"

echo "[INFO] Main LSD2 tree score written to: $MAIN_TREE_SCORE" | tee -a "$RUNLOG"

SENS_REPORT="${LSD2_SENS_DIR}/lsd2_sensitivity_report.tsv"

if [[ -s "$SENS_SCORE_SUMMARY" && -s "$MAIN_TREE_SCORE" ]]; then
  Rscript --vanilla - "$MAIN_TREE_SCORE" "$SENS_SCORE_SUMMARY" "$SENS_RUN_SUMMARY" "$SENS_REPORT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

main_file <- args[1]
sens_file <- args[2]
runs_file <- args[3]
out_file  <- args[4]

main <- read.delim(main_file, check.names = FALSE)
sens <- read.delim(sens_file, check.names = FALSE)
runs <- read.delim(runs_file, check.names = FALSE)

sens <- merge(
  sens,
  runs[, c("label", "status")],
  by = "label",
  all.x = TRUE,
  sort = FALSE
)

main$status <- "MAIN"
all <- rbind(main, sens)

main_height <- main$tree_height[1]
main_total  <- main$total_tree_length[1]
main_zero   <- main$n_zero_exact[1]
main_near10 <- main$n_near_zero_1e10[1]

report <- within(all, {
  delta_tree_height_vs_main <- tree_height - main_height
  pct_tree_height_vs_main <- 100 * (tree_height / main_height - 1)
  delta_total_tree_length_vs_main <- total_tree_length - main_total
  pct_total_tree_length_vs_main <- 100 * (total_tree_length / main_total - 1)
  delta_zero_exact_vs_main <- n_zero_exact - main_zero
  delta_near_zero_1e10_vs_main <- n_near_zero_1e10 - main_near10
})

report <- report[
  order(report$n_zero_exact, report$n_near_zero_1e10, report$tip_depth_span),
]

write.table(
  report,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[INFO] Sensitivity report, sorted by exact/near-zero branches:\n")
print(
  report[, c(
  "label",
  "status",
  "collapse_threshold",
  "min_internal_blen",
  "round_time",
  "n_zero_exact",
  "n_zero_internal_exact",
  "n_zero_terminal_exact",
  "n_near_zero_1e10",
  "frac_near_zero_1e10",
  "tree_height",
  "total_tree_length",
  "tip_depth_span",
  "delta_tree_height_vs_main",
  "delta_total_tree_length_vs_main"
)],
  row.names = FALSE
)
RSCRIPT
  
  echo "[INFO] LSD2 sensitivity report written to: $SENS_REPORT" | tee -a "$RUNLOG"
else
  echo "[WARN] Could not create sensitivity report; missing $SENS_SCORE_SUMMARY or $MAIN_TREE_SCORE" | tee -a "$RUNLOG"
fi

echo "[INFO] Standalone LSD2 dating, sensitivity analysis, and reporting completed." | tee -a "$RUNLOG"

