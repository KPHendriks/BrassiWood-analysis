#!/usr/bin/env Rscript

## BayesTraits: correlated evolution of woodiness (tip-level) and environment (drought + frost)
##
## Update (2026-02-23):
##   - Adds a mechanistic frost threshold for BIO6:
##       * thr_0C : frost = 1 if BIO6 <= 0
##   - Keeps percentile thresholds p40/p50/p60 (for both drought and frost)
##   - Plots can overlay p40/p50/p60 (and 0°C for frost) as before

## Outputs
##   WP2_BrassiNiche/results_final/5_bayestraits/
##     input/
##       brassicaceae_tree_bc.nex
##       bayestraits_wood_env__drought__thr_p40.txt  (and p50/p60)
##       bayestraits_wood_env__frost__thr_p40.txt    (and p50/p60 and thr_0C)
##       qc_tip_table_with_species_env_wood__*__*.csv
##       qc_bayestraits_thresholds__drought.csv
##       qc_bayestraits_thresholds__frost.csv
##       cmd_*_template.txt
##       bt_tasks.tsv
##       run_bt_parallel_all.sh
##     output/
##       runs_bt/drought/thr_p40/independent/chain_01/...
##       runs_bt/frost/thr_0C/dependent/chain_01/...
##     10_*.csv, 11_*.png, 12_*.csv, 13_*.csv, 13_*.txt, 13_*.png

## 0. Prepare ----

out_dir <- "WP2_BrassiNiche/results_final/5_bayestraits"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "input"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "output"), recursive = TRUE, showWarnings = FALSE)

message("[INFO] Output dir: ", normalizePath(out_dir, mustWork = FALSE))

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(coda)
  library(phytools)
})

stop_missing <- function(path) if (!file.exists(path)) stop("[ERROR] Missing: ", path)

std_species <- function(x) gsub(" ", "_", as.character(x))

bt_symbol <- function(x, missing_symbol = "-") {
  ifelse(is.na(x), missing_symbol, as.character(as.integer(x)))
}

fmt_int <- function(x) format(as.integer(round(x)), scientific = FALSE, trim = TRUE)

## 1. Inputs (paths) ----

path_2f <- "WP2_BrassiNiche/results_final/3_woodiness_evolution/8_objects.rds"

cli_niche <- list(
  ge1  = "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge1.csv",
  ge5  = "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge5.csv",
  ge10 = "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge10.csv"
)

stop_missing(path_2f)
invisible(lapply(cli_niche, stop_missing))

## Local publication input CSVs
cli <- list()
cli$species_csv  <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
cli$specimen_csv <- "WP2_BrassiNiche/data/brassiwood_specimens_publication.csv"

stop_missing(cli$species_csv)
stop_missing(cli$specimen_csv)


## 2. BayesTraits settings ----

tier_preference <- c("ge10", "ge5", "ge1")
mcwd_col        <- "MCWD_med"
bio6_col        <- "BIO6_med"
flip_sign_mcwd  <- FALSE
flip_sign_bio6  <- FALSE
missing_symbol  <- "-"

analyses <- c("drought", "frost")

## Percentile thresholds (both drought & frost)
#threshold_percentiles <- c(50)
threshold_percentiles <- c(40, 50, 60)

## Additional mechanistic frost threshold
use_frost_zeroC <- TRUE   # if TRUE: add thr_0C (BIO6 <= 0)

n_reps_each_model <- 20L #3L would be good for a test run
max_parallel_jobs <- 10L #50L whenever MaaS36 HPC is fully available for this analysis

threshold_percentiles <- c(40, 50, 60)
use_frost_zeroC <- TRUE

# # Use these numbers for a test run
# bt <- list(
#   burnin      = 200000L,
#   iterations  = 2000000L,
#   sample      = 2000L,
#   prior_mean  = 10,
#   stones      = 50L,
#   stone_iter  = 2000L
# )

# Use these numbers for a full run
bt <- list(
  burnin      = 1000000L,
  iterations  = 10000000L,
  sample      = 5000L,
  prior_mean  = 10,
  stones      = 100L,
  stone_iter  = 5000L
)

## 3. Load saved objects from 2f ----

obj2f <- readRDS(path_2f)

need_2f <- c("tree_bc", "woody_fac_bc")
miss_2f <- setdiff(need_2f, names(obj2f))

if (length(miss_2f) > 0) {
  stop(
    "[ERROR] 2f object lacks required object(s): ",
    paste(miss_2f, collapse = ", "),
    "\nExpected these in: ", path_2f,
    "\nRerun updated 2f.StudyWoodinessEvolution.R first."
  )
}

tree_bc <- obj2f$tree_bc
woody_fac_bc <- obj2f$woody_fac_bc

if (!inherits(tree_bc, "phylo")) {
  stop("[ERROR] obj2f$tree_bc is not a phylo object.")
}

if (is.null(names(woody_fac_bc))) {
  stop("[ERROR] obj2f$woody_fac_bc has no names.")
}

if (!all(tree_bc$tip.label %in% names(woody_fac_bc))) {
  stop("[ERROR] Not all tree_bc tips are present in woody_fac_bc.")
}

woody_fac_bc <- woody_fac_bc[tree_bc$tip.label]

if (any(is.na(woody_fac_bc))) {
  stop("[ERROR] woody_fac_bc contains NA after alignment to tree_bc.")
}

if (!all(as.character(woody_fac_bc) %in% c("0", "1"))) {
  stop("[ERROR] woody_fac_bc must contain only states '0' and '1'.")
}

if (any(is.na(tree_bc$edge.length))) {
  stop("[ERROR] tree_bc has NA edge lengths.")
}

if (any(tree_bc$edge.length < 0, na.rm = TRUE)) {
  stop("[ERROR] tree_bc has negative edge lengths.")
}

if (any(tree_bc$edge.length == 0, na.rm = TRUE)) {
  message("[INFO] Replacing zero-length edges in tree_bc with 1e-10.")
  tree_bc$edge.length[tree_bc$edge.length == 0] <- 1e-10
}

tip_depths <- ape::node.depth.edgelength(tree_bc)[seq_along(tree_bc$tip.label)]
tree_height <- max(tip_depths, na.rm = TRUE)
tip_depth_span <- diff(range(tip_depths, na.rm = TRUE))
ultra_tol_abs <- max(1e-6, tree_height * 1e-5)

message("[INFO] Loaded 2f BayesTraits tree object:")
message("[INFO]   source: ", normalizePath(path_2f))
message("[INFO]   tips: ", length(tree_bc$tip.label))
message("[INFO]   tree height: ", signif(tree_height, 6), " Ma")
message("[INFO]   tip-depth span: ", signif(tip_depth_span, 6), " Ma")
message("[INFO]   ultrametric tolerance: ", signif(ultra_tol_abs, 6), " Ma")

if (!ape::is.ultrametric(tree_bc, tol = ultra_tol_abs)) {
  stop(
    "[ERROR] tree_bc from 2f is not ultrametric within tolerance. ",
    "Rerun/check 2f before BayesTraits."
  )
}

message("[INFO] woody_fac_bc summary:")
print(table(woody_fac_bc, useNA = "ifany"))

## 4. Read local specimen_details + species_details CSVs ----

message("[INFO] Reading specimen metadata from local CSV: ", cli$specimen_csv)
specimen_details <- readr::read_csv(cli$specimen_csv, show_col_types = FALSE)

message("[INFO] Reading species metadata from local CSV: ", cli$species_csv)
species_details <- readr::read_csv(cli$species_csv, show_col_types = FALSE)

needed_specimen_cols <- c(
  "SAMPLE",
  "SPECIES_NAME_PRINT",
  "TRIBE",
  "SUPERTRIBE",
  "FAMILY",
  "GENUS",
  "loci_remaining"
)

needed_species_cols <- c(
  "SPECIES_NAME_PRINT",
  "GROWTH_FORM"
)

missing_specimen_cols <- setdiff(needed_specimen_cols, names(specimen_details))
missing_species_cols  <- setdiff(needed_species_cols, names(species_details))

if (length(missing_specimen_cols) > 0) {
  stop(
    "[ERROR] Missing columns in specimen CSV: ",
    paste(missing_specimen_cols, collapse = ", ")
  )
}

if (length(missing_species_cols) > 0) {
  stop(
    "[ERROR] Missing columns in species CSV: ",
    paste(missing_species_cols, collapse = ", ")
  )
}

specimen_details <- specimen_details %>%
  mutate(
    SAMPLE             = as.character(SAMPLE),
    SPECIES_NAME_PRINT = str_squish(as.character(SPECIES_NAME_PRINT)),
    TRIBE              = as.character(TRIBE),
    SUPERTRIBE         = as.character(SUPERTRIBE),
    FAMILY             = as.character(FAMILY),
    GENUS              = as.character(GENUS),
    loci_remaining     = suppressWarnings(as.numeric(loci_remaining)),
    tip_label          = paste0(TRIBE, "_", SPECIES_NAME_PRINT, " (", SAMPLE, ")")
  )
specimen_details$loci_remaining[is.na(specimen_details$loci_remaining)] <- 0

tip_info <- specimen_details %>%
  dplyr::select(SAMPLE, SPECIES_NAME_PRINT, FAMILY, TRIBE, SUPERTRIBE, GENUS) %>%
  dplyr::distinct() %>%
  dplyr::filter(SAMPLE %in% tree_bc$tip.label)

message("[INFO] tip_info rows (restricted to tree tips): ", nrow(tip_info))

## 5. Choose tier and extract species-level MCWD + BIO6 ----

avail_tiers <- names(cli_niche)

tier_use <- tier_preference[tier_preference %in% avail_tiers][1]

if (is.na(tier_use) || is.null(tier_use)) {
  stop("[ERROR] None of tier_preference found. Available: ", paste(avail_tiers, collapse = ", "))
}

message("[INFO] Using tier for climate: ", tier_use)

df_tier <- readr::read_csv(cli_niche[[tier_use]], show_col_types = FALSE)

if (is.null(df_tier) || nrow(df_tier) == 0) {
  stop("[ERROR] tier table is empty: ", cli_niche[[tier_use]])
}

if (!("species" %in% names(df_tier))) stop("[ERROR] tier table lacks column 'species'.")
if (!(mcwd_col %in% names(df_tier))) stop("[ERROR] tier table lacks column: ", mcwd_col)
if (!(bio6_col %in% names(df_tier))) stop("[ERROR] tier table lacks column: ", bio6_col)

if (is.null(df_tier) || nrow(df_tier) == 0) stop("[ERROR] tier model table is empty.")
if (!("species" %in% names(df_tier))) stop("[ERROR] tier model table lacks column 'species'.")
if (!(mcwd_col %in% names(df_tier))) stop("[ERROR] tier model table lacks column: ", mcwd_col)
if (!(bio6_col %in% names(df_tier))) stop("[ERROR] tier model table lacks column: ", bio6_col)

clim_species <- df_tier %>%
  dplyr::transmute(
    species = str_squish(as.character(species)),
    MCWD = suppressWarnings(as.numeric(.data[[mcwd_col]])),
    BIO6 = suppressWarnings(as.numeric(.data[[bio6_col]]))
  ) %>%
  dplyr::distinct(species, .keep_all = TRUE)

if (flip_sign_mcwd) clim_species$MCWD <- -clim_species$MCWD
if (flip_sign_bio6) clim_species$BIO6 <- -clim_species$BIO6

## 6. Build base tip table: SAMPLE -> species -> MCWD + BIO6 + wood ----

tip_tbl_base <- tibble(SAMPLE = tree_bc$tip.label) %>%
  dplyr::left_join(
    tip_info %>% transmute(
      SAMPLE = str_squish(as.character(SAMPLE)),
      SPECIES_NAME_PRINT = str_squish(as.character(SPECIES_NAME_PRINT)),
      species = std_species(SPECIES_NAME_PRINT)
    ),
    by = "SAMPLE"
  ) %>%
  dplyr::left_join(clim_species, by = "species")

wood_tip <- woody_fac_bc[tree_bc$tip.label]
wood_num <- suppressWarnings(as.numeric(as.character(wood_tip)))
wood_num[!wood_num %in% c(0, 1)] <- NA_real_

tip_tbl_base <- tip_tbl_base %>% mutate(wood = wood_num)

## 7. Build thresholds for drought + frost, write inputs ----

in_dir  <- file.path(out_dir, "input")
out_bt  <- file.path(out_dir, "output")
dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_bt, recursive = TRUE, showWarnings = FALSE)

tree_path <- file.path(in_dir, "brassicaceae_tree_bc.nex")
ape::write.nexus(tree_bc, file = tree_path)
message("[OK] Wrote tree: ", tree_path)

get_env_vec <- function(analysis, tip_tbl) {
  if (analysis == "drought") return(tip_tbl$MCWD)
  if (analysis == "frost")   return(tip_tbl$BIO6)
  stop("[ERROR] Unknown analysis: ", analysis)
}

compute_thr_tbl_percentiles <- function(analysis, tip_tbl, percentiles) {
  v <- get_env_vec(analysis, tip_tbl)
  v <- v[is.finite(v)]
  if (length(v) < 100) stop("[ERROR] Too few finite values for analysis '", analysis, "' to define percentile thresholds.")
  tibble(
    analysis = analysis,
    thr_type = "percentile",
    percentile = percentiles
  ) %>%
    mutate(
      thr_value = as.numeric(stats::quantile(v, probs = percentile / 100, names = FALSE, na.rm = TRUE)),
      thr_label = paste0("thr_p", percentile)
    )
}

compute_thr_tbl_fixed_frost0C <- function() {
  tibble(
    analysis = "frost",
    thr_type = "fixed",
    percentile = NA_real_,
    thr_value = 0,
    thr_label = "thr_0C"
  )
}

thr_tbl_all <- bind_rows(
  lapply(analyses, compute_thr_tbl_percentiles, tip_tbl = tip_tbl_base, percentiles = threshold_percentiles),
  if (isTRUE(use_frost_zeroC)) compute_thr_tbl_fixed_frost0C() else NULL
) %>%
  mutate(
    analysis = factor(analysis, levels = analyses)
  ) %>%
  arrange(analysis, dplyr::case_when(thr_label == "thr_0C" ~ 0L, TRUE ~ 1L), percentile)

write_csv(thr_tbl_all %>% filter(analysis == "drought"),
          file.path(in_dir, "qc_bayestraits_thresholds__drought.csv"))
write_csv(thr_tbl_all %>% filter(analysis == "frost"),
          file.path(in_dir, "qc_bayestraits_thresholds__frost.csv"))

message("[QC] Thresholds used:")
print(thr_tbl_all)

make_tip_tbl_for_threshold <- function(tip_tbl, analysis, thr_label, thr_value) {
  
  if (analysis == "drought") {
    tip_tbl %>%
      mutate(
        env_value = MCWD,
        env = case_when(
          is.na(MCWD) ~ NA_real_,
          MCWD <= thr_value ~ 1,
          TRUE ~ 0
        )
      )
    
  } else if (analysis == "frost") {
    
    ## NOTE: For frost, stress = frost. We keep convention:
    ##   env = 1 indicates "stress" (frost present/stronger)
    ## For thr_0C: frost if BIO6 <= 0
    tip_tbl %>%
      mutate(
        env_value = BIO6,
        env = case_when(
          is.na(BIO6) ~ NA_real_,
          BIO6 <= thr_value ~ 1,
          TRUE ~ 0
        )
      )
    
  } else {
    stop("[ERROR] Unknown analysis: ", analysis)
  }
}

write_bt_data_for_threshold <- function(analysis, thr_label, thr_value) {
  
  tip_tbl <- make_tip_tbl_for_threshold(tip_tbl_base, analysis, thr_label, thr_value)
  
  qc <- tibble(
    analysis = analysis,
    thr_label = thr_label,
    thr_percentile = suppressWarnings(as.numeric(sub("^thr_p", "", thr_label))),
    threshold_value = thr_value,
    n_tips_tree = length(tree_bc$tip.label),
    n_mapped_species = sum(!is.na(tip_tbl$species)),
    n_wood_present = sum(is.finite(tip_tbl$wood)),
    n_env_value_present = sum(is.finite(tip_tbl$env_value)),
    n_env_present = sum(is.finite(tip_tbl$env)),
    env_1_n = sum(tip_tbl$env == 1, na.rm = TRUE),
    env_0_n = sum(tip_tbl$env == 0, na.rm = TRUE)
  )
  
  write_csv(qc, file.path(in_dir, paste0("qc_input_coverage__", analysis, "__", thr_label, ".csv")))
  write_csv(tip_tbl, file.path(in_dir, paste0("qc_tip_table_with_species_env_wood__", analysis, "__", thr_label, ".csv")))
  
  bt_data <- tip_tbl %>%
    dplyr::transmute(
      Taxon = SAMPLE,
      wood = bt_symbol(wood, missing_symbol = missing_symbol),
      env  = bt_symbol(env,  missing_symbol = missing_symbol)
    )
  
  data_path <- file.path(in_dir, paste0("bayestraits_wood_env__", analysis, "__", thr_label, ".txt"))
  write_delim(bt_data, file = data_path, delim = "\t", col_names = TRUE)
  
  message("[OK] Wrote data: ", data_path,
          " | analysis=", analysis,
          " | cutoff=", signif(thr_value, 6),
          " | env1=", qc$env_1_n, " env0=", qc$env_0_n)
  
  list(qc = qc, data_path = data_path)
}

thr_outputs <- purrr::pmap(thr_tbl_all, \(analysis, thr_type, percentile, thr_value, thr_label) {
  write_bt_data_for_threshold(as.character(analysis), as.character(thr_label), thr_value)
})

qc_all <- dplyr::bind_rows(purrr::map(thr_outputs, "qc"))

write_csv(qc_all, file.path(out_dir, "07_qc_input_coverage_all_thresholds.csv"))

## Extra QC: compare frost thr_0C vs thr_p50 (coverage + balance)
qc_frost_compare <- qc_all %>%
  dplyr::filter(analysis == "frost", thr_label %in% c("thr_0C", "thr_p50")) %>%
  dplyr::mutate(
    env1_frac = if_else(n_env_present > 0, env_1_n / n_env_present, NA_real_),
    env0_frac = if_else(n_env_present > 0, env_0_n / n_env_present, NA_real_),
    env1_to_env0 = if_else(env_0_n > 0, env_1_n / env_0_n, NA_real_)
  ) %>%
  dplyr::select(
    analysis, thr_label, threshold_value,
    n_tips_tree, n_mapped_species, n_wood_present,
    n_env_value_present, n_env_present,
    env_1_n, env_0_n, env1_frac, env0_frac, env1_to_env0
  ) %>%
  dplyr::arrange(match(thr_label, c("thr_0C", "thr_p50")))

write_csv(qc_frost_compare, file.path(out_dir, "07_qc_frost_thr0C_vs_p50_coverage_balance.csv"))

message("[QC] Wrote frost threshold comparison: ",
        normalizePath(file.path(out_dir, "07_qc_frost_thr0C_vs_p50_coverage_balance.csv")))


## 8. Write BayesTraits command templates + task table + single runner ----

write_bt_cmd_template <- function(path,
                                  model = c("independent", "dependent"),
                                  use_rj = FALSE,
                                  logprefix_placeholder = "__LOGPREFIX__",
                                  seed_placeholder = "__SEED__") {
  model <- match.arg(model)
  lines <- c(
    if (model == "independent") "2" else "3",
    "2",
    paste("LogFile", logprefix_placeholder),
    paste("Seed", seed_placeholder),
    if (use_rj) paste("RevJump exp", bt$prior_mean) else paste("PriorAll exp", bt$prior_mean),
    paste("burnin", fmt_int(bt$burnin)),
    paste("iterations", fmt_int(bt$iterations)),
    paste("sample", fmt_int(bt$sample)),
    if (!is.null(bt$stones) && bt$stones > 0) paste("Stones", fmt_int(bt$stones), fmt_int(bt$stone_iter)),
    "Run",
    ""
  )
  writeLines(lines, con = path)
  invisible(path)
}

cmd_ind_tpl  <- file.path(in_dir, "cmd_ind_template.txt")
cmd_dep_tpl  <- file.path(in_dir, "cmd_dep_template.txt")
cmd_depR_tpl <- file.path(in_dir, "cmd_depR_template.txt")

write_bt_cmd_template(cmd_ind_tpl,  model = "independent", use_rj = FALSE)
write_bt_cmd_template(cmd_dep_tpl,  model = "dependent",   use_rj = FALSE)
write_bt_cmd_template(cmd_depR_tpl, model = "dependent",   use_rj = TRUE)

make_tasks <- function(thr_tbl_all,
                       n_rep = 20L,
                       models = c("independent", "dependent", "dependent_RJ")) {
  
  models <- match.arg(models, several.ok = TRUE)
  
  make_block <- function(analysis, thr_label, model, n, seed_base) {
    tibble(
      analysis = analysis,
      thr_label = thr_label,
      model = model,
      chain = sprintf("%02d", seq_len(n)),
      seed = as.integer(seed_base + seq_len(n))
    )
  }
  
  all_blocks <- list()
  for (i in seq_len(nrow(thr_tbl_all))) {
    a <- as.character(thr_tbl_all$analysis[i])
    t <- as.character(thr_tbl_all$thr_label[i])
    idx <- 0L
    for (m in models) {
      idx <- idx + 1L
      ## stable seed structure: depends on (analysis, thr_label, model, chain)
      seed_base <- 1000000L * i + 10000L * idx + 100L * 1L
      all_blocks[[length(all_blocks) + 1]] <- make_block(a, t, m, n_rep, seed_base)
    }
  }
  
  bind_rows(all_blocks)
}

tasks <- make_tasks(thr_tbl_all, n_rep = n_reps_each_model, models = c("independent", "dependent", "dependent_RJ"))

tasks_path <- file.path(in_dir, "bt_tasks.tsv")
write_tsv(tasks, file = tasks_path)

run_sh <- file.path(in_dir, "run_bt_parallel_all.sh")

run_lines <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# Run from the BayesTraits input directory on HPC:",
  "#   cd .../5_bayestraits/input",
  "#   module load parallel || true",
  "#   JOBS=50 ./run_bt_parallel_all.sh /path/to/BayesTraitsV5",
  "#",
  "# This runs ALL tasks in bt_tasks.tsv:",
  "#   analysis (drought,frost) x thresholds x models x chains",
  "# Outputs go to ../output/runs_bt/<analysis>/<thr_label>/<model>/chain_<chain>/",
  "",
  "BT_BIN=${1:-}",
  "if [[ -z \"$BT_BIN\" ]]; then",
  "  echo \"[ERROR] Provide BT binary path: ./run_bt_parallel_all.sh /path/to/BayesTraitsV5\" >&2",
  "  exit 1",
  "fi",
  "if [[ ! -x \"$BT_BIN\" ]]; then",
  "  echo \"[ERROR] BT binary not executable: $BT_BIN\" >&2",
  "  exit 1",
  "fi",
  "",
  "JOBS=${JOBS:-50}",
  "DRYRUN=${DRYRUN:-0}",
  "",
  "IN_DIR=$(pwd)",
  "OUT_DIR=$(cd .. && pwd)/output",
  "RUNS_DIR=\"$OUT_DIR/runs_bt\"",
  "mkdir -p \"$RUNS_DIR\"",
  "",
  "TREE=\"$IN_DIR/brassicaceae_tree_bc.nex\"",
  "TASKS=\"$IN_DIR/bt_tasks.tsv\"",
  "TPL_IND=\"$IN_DIR/cmd_ind_template.txt\"",
  "TPL_DEP=\"$IN_DIR/cmd_dep_template.txt\"",
  "TPL_DEPR=\"$IN_DIR/cmd_depR_template.txt\"",
  "",
  "for f in \"$TREE\" \"$TASKS\" \"$TPL_IND\" \"$TPL_DEP\" \"$TPL_DEPR\"; do",
  "  [[ -f \"$f\" ]] || { echo \"[ERROR] Missing required file: $f\" >&2; exit 1; }",
  "done",
  "",
  "if ! command -v parallel >/dev/null 2>&1; then",
  "  echo \"[ERROR] GNU parallel not found. Try: module load parallel\" >&2",
  "  exit 1",
  "fi",
  "",
  "data_for_task() {",
  "  local analysis=\"$1\"",
  "  local thr_label=\"$2\"",
  "  echo \"$IN_DIR/bayestraits_wood_env__${analysis}__${thr_label}.txt\"",
  "}",
  "",
  "run_one() {",
  "  local analysis=\"$1\"",
  "  local thr_label=\"$2\"",
  "  local model=\"$3\"",
  "  local chain=\"$4\"",
  "  local seed=\"$5\"",
  "",
  "  local data",
  "  data=$(data_for_task \"$analysis\" \"$thr_label\")",
  "  if [[ ! -f \"$data\" ]]; then",
  "    echo \"[ERROR] Missing data file: $data\" >&2",
  "    return 2",
  "  fi",
  "",
  "  local chain_dir=\"$RUNS_DIR/${analysis}/${thr_label}/${model}/chain_${chain}\"",
  "  mkdir -p \"$chain_dir\"",
  "",
  "  local tpl logprefix",
  "  if [[ \"$model\" == \"independent\" ]]; then",
  "    tpl=\"$TPL_IND\"",
  "    logprefix=\"Run_${analysis}_${thr_label}_independent_chain_${chain}\"",
  "  elif [[ \"$model\" == \"dependent\" ]]; then",
  "    tpl=\"$TPL_DEP\"",
  "    logprefix=\"Run_${analysis}_${thr_label}_dependent_chain_${chain}\"",
  "  elif [[ \"$model\" == \"dependent_RJ\" ]]; then",
  "    tpl=\"$TPL_DEPR\"",
  "    logprefix=\"Run_${analysis}_${thr_label}_dependent_RJ_chain_${chain}\"",
  "  else",
  "    echo \"[ERROR] Unknown model: $model\" >&2",
  "    return 2",
  "  fi",
  "",
  "  local cmdfile=\"$chain_dir/bayestraits_cmd.txt\"",
  "  sed -e \"s/__SEED__/${seed}/g\" -e \"s/__LOGPREFIX__/${logprefix}/g\" \"$tpl\" > \"$cmdfile\"",
  "",
  "  local stdout=\"$chain_dir/${logprefix}.stdout.txt\"",
  "  local stderr=\"$chain_dir/${logprefix}.stderr.txt\"",
  "",
  "  if [[ \"$DRYRUN\" == \"1\" ]]; then",
  "    echo \"[DRYRUN] (cd $chain_dir && $BT_BIN $TREE $data < $cmdfile > $stdout 2> $stderr)\"",
  "    return 0",
  "  fi",
  "",
  "  (",
  "    cd \"$chain_dir\"",
  "    echo \"[INFO] $(date) | start | analysis=$analysis thr=$thr_label model=$model chain=$chain seed=$seed\" >&2",
  "    \"$BT_BIN\" \"$TREE\" \"$data\" < \"$cmdfile\" > \"$stdout\" 2> \"$stderr\"",
  "    echo \"[INFO] $(date) | done  | analysis=$analysis thr=$thr_label model=$model chain=$chain\" >&2",
  "  )",
  "}",
  "",
  "export -f run_one data_for_task",
  "export BT_BIN TREE IN_DIR RUNS_DIR TPL_IND TPL_DEP TPL_DEPR DRYRUN",
  "",
  "echo \"[INFO] First tasks:\"",
  "awk -F'\\t' 'NR==1{next} {print \"  - analysis=\" $1 \" thr=\" $2 \" model=\" $3 \" chain=\" $4 \" seed=\" $5} NR==11{exit}' \"$TASKS\"",
  "n_tasks=$(awk 'NR>1{c++} END{print c+0}' \"$TASKS\")",
  "echo \"[INFO] Total tasks: $n_tasks | parallel jobs: $JOBS\"",
  "",
  "# analysis\\tthr_label\\tmodel\\tchain\\tseed",
  "tail -n +2 \"$TASKS\" | parallel --colsep '\\t' -j \"$JOBS\" run_one {1} {2} {3} {4} {5}",
  "",
  "echo \"[OK] All BayesTraits tasks finished. Outputs in: $RUNS_DIR\""
)

writeLines(run_lines, con = run_sh)
Sys.chmod(run_sh, mode = "0755")

message("[OK] Wrote templates + tasks + runner to input/:")
message("  - ", cmd_ind_tpl)
message("  - ", cmd_dep_tpl)
message("  - ", cmd_depR_tpl)
message("  - ", tasks_path)
message("  - ", run_sh)
message("[INFO] Total tasks: ", nrow(tasks), " = ", nrow(thr_tbl_all), " analysis-thresholds x 3 models x ", n_reps_each_model, " reps")
message("[INFO] Suggested HPC run: JOBS=", max_parallel_jobs, " ./run_bt_parallel_all.sh /path/to/BayesTraitsV5")



## 9 Now run BayesTraits on HPC! ----
message("[NEXT] Sync input/ to HPC, run the bash script, then sync output/ back.")

# Run on the PHC as follows; set parallel jobs to the number of CPUs available:
# cd /path/to/5_bayestraits/input
#JOBS=52 ./run_bt_parallel_all.sh /home/kasper.hendriks/BayesTraitsV5.0.3-Linux/BayesTraitsV5



## 10. Read BayesTraits outputs (all analyses/thresholds/models/chains) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(grid)   # unit()
})

runs_dir <- file.path(out_dir, "output", "runs_bt")
if (!dir.exists(runs_dir)) {
  stop("[ERROR] runs_dir does not exist: ", runs_dir, "\nHave you synced HPC results back?")
}

parse_path_meta <- function(path) {
  p <- gsub("\\\\", "/", path)
  m <- str_match(p, "runs_bt/([^/]+)/([^/]+)/([^/]+)/chain_([0-9]+)/")
  if (any(is.na(m))) {
    return(tibble(analysis = NA_character_, thr_label = NA_character_, model = NA_character_, chain = NA_character_))
  }
  tibble(analysis = m[, 2], thr_label = m[, 3], model = m[, 4], chain = m[, 5])
}

repair_names <- function(nms) {
  nms <- ifelse(is.na(nms) | nms == "", paste0("unnamed_", seq_along(nms)), nms)
  make.unique(nms, sep = "_")
}

q_params <- c("q12","q13","q21","q24","q31","q34","q42","q43")

thr_levels_all <- thr_tbl_all %>%
  arrange(analysis, dplyr::case_when(thr_label == "thr_0C" ~ 0L, TRUE ~ 1L), percentile) %>%
  pull(thr_label) %>%
  unique()

## 10.1 Stones
stones_files <- list.files(runs_dir, pattern = "\\.Stones\\.txt$", recursive = TRUE, full.names = TRUE)

parse_stones <- function(f) {
  x <- readLines(f, warn = FALSE)
  idx <- grep("^Log marginal likelihood:", x)
  logml <- NA_real_
  if (length(idx) > 0) {
    line <- x[idx[length(idx)]]
    logml <- suppressWarnings(as.numeric(str_trim(sub("^Log marginal likelihood:\\s*", "", line))))
  }
  bind_cols(tibble(file = f, log_marginal_likelihood = logml), parse_path_meta(f))
}

stones_tbl <- map_dfr(stones_files, parse_stones) %>%
  mutate(
    analysis = factor(analysis, levels = analyses),
    model = factor(model, levels = c("independent", "dependent", "dependent_RJ")),
    thr_label = factor(thr_label, levels = thr_levels_all),
    chain = suppressWarnings(as.integer(chain))
  )

write_csv(stones_tbl, file.path(out_dir, "10_stones_logml_by_chain.csv"))

## 10.2 Schedule
schedule_files <- list.files(runs_dir, pattern = "\\.Schedule\\.txt$", recursive = TRUE, full.names = TRUE)

parse_schedule <- function(f) {
  x <- readLines(f, warn = FALSE)
  hdr <- grep("^Rate Tried\\s+% Accepted", x)
  meta <- parse_path_meta(f)
  
  if (length(hdr) == 0) {
    return(bind_cols(tibble(file = f, sample_ave_accept = NA_real_, total_ave_accept = NA_real_), meta))
  }
  
  dat <- x[(hdr[length(hdr)] + 1):length(x)]
  dat <- dat[nzchar(trimws(dat))]
  dat <- dat[grepl("^[0-9]", dat)]
  if (length(dat) == 0) {
    return(bind_cols(tibble(file = f, sample_ave_accept = NA_real_, total_ave_accept = NA_real_), meta))
  }
  
  last <- dat[length(dat)]
  parts <- str_split(trimws(last), "\\s+")[[1]]
  sample_ave <- suppressWarnings(as.numeric(parts[length(parts) - 1]))
  total_ave  <- suppressWarnings(as.numeric(parts[length(parts)]))
  
  bind_cols(tibble(file = f, sample_ave_accept = sample_ave, total_ave_accept = total_ave), meta)
}

schedule_tbl <- map_dfr(schedule_files, parse_schedule) %>%
  mutate(
    analysis = factor(analysis, levels = analyses),
    model = factor(model, levels = c("independent", "dependent", "dependent_RJ")),
    thr_label = factor(thr_label, levels = thr_levels_all),
    chain = suppressWarnings(as.integer(chain))
  )

write_csv(schedule_tbl, file.path(out_dir, "10_schedule_accept_by_chain.csv"))

## 10.3 Log files (posterior draws)
log_files <- list.files(runs_dir, pattern = "\\.Log\\.txt$", recursive = TRUE, full.names = TRUE)

read_bayestraits_log_draws <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  i <- grep("^Iteration\\s", lines)
  if (length(i) == 0) {
    warning("[WARN] No 'Iteration' header found in: ", path)
    return(NULL)
  }
  
  tab_txt <- paste(lines[i[1]:length(lines)], collapse = "\n")
  
  df <- tryCatch(
    read.delim(
      text = tab_txt,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      fill = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      warning("[WARN] Failed to parse log: ", path, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  names(df) <- repair_names(names(df))
  
  meta <- parse_path_meta(path)
  meta$file <- path
  
  bind_cols(meta[rep(1, nrow(df)), , drop = FALSE], as_tibble(df))
}

log_draws_tbl <- map_dfr(log_files, read_bayestraits_log_draws) %>%
  dplyr::mutate(
    analysis = factor(analysis, levels = analyses),
    model = factor(model, levels = c("independent", "dependent", "dependent_RJ")),
    thr_label = factor(thr_label, levels = thr_levels_all),
    chain = suppressWarnings(as.integer(chain))
  )

if (nrow(log_draws_tbl) > 0 && all(c("Iteration", "Lh") %in% names(log_draws_tbl))) {
  write_csv(
    log_draws_tbl %>%
      dplyr::select(analysis, thr_label, model, chain, Iteration, Lh) %>%
      dplyr::arrange(analysis, thr_label, model, chain, Iteration),
    file.path(out_dir, "10_log_draws_index.csv")
  )
} else {
  warning("[WARN] log_draws_tbl is empty or missing Iteration/Lh; check that logs parsed correctly.")
  write_csv(
    tibble(analysis = character(), thr_label = character(), model = character(), chain = integer(), Iteration = numeric(), Lh = numeric()),
    file.path(out_dir, "10_log_draws_index.csv")
  )
}

## 10.4 RJ inclusion from Model string
rj_inclusion_from_model_string <- function(model_string, rate_names = q_params) {
  s <- gsub("'", "", model_string, fixed = TRUE)
  s <- trimws(s)
  nums <- suppressWarnings(as.integer(strsplit(s, "\\s+")[[1]]))
  
  out <- rep(NA_integer_, length(rate_names))
  names(out) <- rate_names
  if (length(nums) != length(rate_names)) return(out)
  
  out <- as.integer(nums > 0)
  names(out) <- rate_names
  out
}

if (!("Model string" %in% names(log_draws_tbl))) {
  warning("[WARN] No 'Model string' column in log_draws_tbl; RJ inclusion tables will be empty.")
  rj_inclusion_by_chain <- tibble()
  rj_inclusion_pooled   <- tibble()
} else {
  
  rj_long <- log_draws_tbl %>%
    dplyr::filter(model == "dependent_RJ") %>%
    dplyr::mutate(.model_string = as.character(.data[["Model string"]])) %>%
    dplyr::mutate(.incl = purrr::map(.model_string, rj_inclusion_from_model_string)) %>%
    tidyr::unnest_wider(.incl, names_sep = "_")
  
  ## By chain
  rj_inclusion_by_chain <- rj_long %>%
    dplyr::group_by(analysis, thr_label, chain) %>%
    dplyr::summarise(
      n_iter = dplyr::n(),
      dplyr::across(
        dplyr::all_of(paste0(".incl_", q_params)),
        ~ mean(.x == 1, na.rm = TRUE),
        .names = "incl{.col}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::rename_with(~ gsub("^incl\\.incl_", "incl_", .x), dplyr::starts_with("incl"))
  
  ## Pooled
  rj_inclusion_pooled <- rj_long %>%
    dplyr::group_by(analysis, thr_label) %>%
    dplyr::summarise(
      n_iter = dplyr::n(),
      dplyr::across(
        dplyr::all_of(paste0(".incl_", q_params)),
        ~ mean(.x == 1, na.rm = TRUE),
        .names = "incl{.col}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::rename_with(~ gsub("^incl\\.incl_", "incl_", .x), dplyr::starts_with("incl"))
}

readr::write_csv(rj_inclusion_by_chain, file.path(out_dir, "10_rj_inclusion_by_chain.csv"))
readr::write_csv(rj_inclusion_pooled,   file.path(out_dir, "10_rj_inclusion_pooled.csv"))

message("[OK] Parsed outputs:")
message("  stones rows:   ", nrow(stones_tbl))
message("  schedule rows: ", nrow(schedule_tbl))
message("  log draws rows:", nrow(log_draws_tbl))
message("  RJ incl rows (by chain): ", nrow(rj_inclusion_by_chain))
message("  RJ incl rows (pooled):   ", nrow(rj_inclusion_pooled))


## 11. Basic plots for SupMat-style QC ----

p_logml <- ggplot(stones_tbl, aes(x = model, y = log_marginal_likelihood)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 1) +
  facet_grid(analysis ~ thr_label, scales = "free_y") +
  labs(x = NULL, y = "Log marginal likelihood (stepping-stone)",
       title = "Stepping-stone log marginal likelihoods by analysis, model, and threshold")

ggsave(file.path(out_dir, "11A_logML_by_model_threshold_and_analysis.png"), p_logml, width = 13, height = 7.5, dpi = 220)

p_acc <- ggplot(schedule_tbl, aes(x = model, y = total_ave_accept)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(analysis ~ thr_label) +
  labs(x = NULL, y = "Total average acceptance rate",
       title = "Acceptance rates by analysis, model, and threshold")

ggsave(file.path(out_dir, "11B_acceptance_by_model_threshold_and_analysis.png"), p_acc, width = 13, height = 7.5, dpi = 220)
ggsave(file.path(out_dir, "11B_acceptance_by_model_threshold_and_analysis.pdf"), p_acc, width = 13, height = 7.5, dpi = 220)

if (requireNamespace("coda", quietly = TRUE)) {
  ess_params <- c("Lh", intersect(q_params, names(log_draws_tbl)))
  
  ess_tbl <- log_draws_tbl %>%
    dplyr::group_by(analysis, thr_label, model, chain) %>%
    group_modify(~{
      df <- .x %>% dplyr::select(any_of(ess_params))
      out <- lapply(names(df), function(p) {
        v <- df[[p]]
        v <- v[is.finite(v)]
        if (length(v) < 50) return(tibble(param = p, ess = NA_real_, n = length(v)))
        tibble(param = p, ess = as.numeric(coda::effectiveSize(v)), n = length(v))
      })
      dplyr::bind_rows(out)
    }) %>%
    ungroup()
  
  write_csv(ess_tbl, file.path(out_dir, "11D_ess_by_chain.csv"))
} else {
  warning("[WARN] Package 'coda' missing; skipping ESS table.")
  write_csv(tibble(), file.path(out_dir, "11D_ess_by_chain.csv"))
}

calc_logbf_pairs <- function(stones_tbl_one_group, modelA = "dependent", modelB = "independent") {
  if (!all(c(modelA, modelB) %in% unique(stones_tbl_one_group$model))) return(NULL)
  
  stones_tbl_one_group %>%
    filter(model %in% c(modelA, modelB)) %>%
    dplyr::select(model, chain, log_marginal_likelihood) %>%
    tidyr::pivot_wider(names_from = model, values_from = log_marginal_likelihood) %>%
    mutate(
      logBF = 2 * (.data[[modelA]] - .data[[modelB]]),
      modelA = modelA, modelB = modelB
    )
}

logbf_dep_ind <- stones_tbl %>%
  dplyr::group_by(analysis, thr_label) %>%
  dplyr::group_modify(~calc_logbf_pairs(.x, "dependent", "independent")) %>%
  dplyr::ungroup()

write_csv(logbf_dep_ind, file.path(out_dir, "11C_logBF_dep_vs_ind_by_threshold_and_analysis.csv"))

p_bf <- ggplot(logbf_dep_ind, aes(x = logBF)) +
  geom_histogram(bins = 25) +
  facet_grid(analysis ~ thr_label, scales = "free_y") +
  labs(x = "Log Bayes Factor: 2*(logML_dep - logML_ind)", y = "Number of paired chains",
       title = "Support for dependent vs independent across analyses and thresholds")

ggsave(file.path(out_dir, "11C_logBF_distribution_dep_vs_ind.png"), p_bf, width = 13, height = 7.5, dpi = 220)
ggsave(file.path(out_dir, "11C_logBF_distribution_dep_vs_ind.pdf"), p_bf, width = 13, height = 7.5, dpi = 220)

p_trace_lh <- ggplot(log_draws_tbl, aes(x = Iteration, y = Lh, group = chain)) +
  geom_line(alpha = 0.6) +
  facet_grid(model ~ analysis + thr_label, scales = "free_y") +
  labs(title = "Trace plots: log-likelihood", x = "Iteration", y = "Lh")

ggsave(file.path(out_dir, "11E_trace_Lh.png"), p_trace_lh, width = 13, height = 8.5, dpi = 220)

message("[OK] Step 11 wrote 11_*.png and 11C_*.csv to ", normalizePath(out_dir))



## 12. Strong sanity checks + absolute q summaries (directionality) ----
## Shared sign convention:
##   For woodiness transitions: x = ln(reduced-stress / stress)
##   x < 0  => transition faster under stronger stress
##   x > 0  => transition faster under reduced / absent stress
##
## Woodiness transitions:
##   Gain (H>W): ln(q13/q24)
##     q13 = gain under reduced stress
##     q24 = gain under stress
##
##   Loss (W>H): ln(q31/q42)
##     q31 = loss under reduced stress
##     q42 = loss under stress
##
## Environmental transitions:
##   Into stress: ln(q34/q12)
##     q12 = reduced stress -> stress in herbaceous lineages
##     q34 = reduced stress -> stress in woody lineages
##
##   Out of stress: ln(q43/q21)
##     q21 = stress -> reduced stress in herbaceous lineages
##     q43 = stress -> reduced stress in woody lineages

summarise_q_direction <- function(df) {
  need <- q_params
  present <- intersect(need, names(df))
  if (length(present) == 0) return(NULL)
  
  cover <- vapply(need, function(p) {
    if (!p %in% names(df)) return(0)
    mean(is.finite(df[[p]]), na.rm = TRUE)
  }, numeric(1))
  
  safe_prob <- function(a, b, op = `>`) {
    if (!(a %in% names(df) && b %in% names(df))) return(NA_real_)
    ok <- is.finite(df[[a]]) & is.finite(df[[b]])
    if (sum(ok) == 0) return(NA_real_)
    mean(op(df[[a]][ok], df[[b]][ok]))
  }
  
  safe_ratio <- function(a, b, probs = c(0.025, 0.975)) {
    if (!(a %in% names(df) && b %in% names(df))) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
    ok <- is.finite(df[[a]]) & is.finite(df[[b]]) & df[[b]] > 0
    if (sum(ok) == 0) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
    r <- df[[a]][ok] / df[[b]][ok]
    c(
      med = median(r),
      lo  = as.numeric(stats::quantile(r, probs[1], names = FALSE)),
      hi  = as.numeric(stats::quantile(r, probs[2], names = FALSE))
    )
  }
  
  safe_logratio <- function(a, b, probs = c(0.025, 0.975)) {
    if (!(a %in% names(df) && b %in% names(df))) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
    ok <- is.finite(df[[a]]) & is.finite(df[[b]]) & df[[a]] > 0 & df[[b]] > 0
    if (sum(ok) == 0) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
    lr <- log(df[[a]][ok] / df[[b]][ok])
    c(
      med = median(lr),
      lo  = as.numeric(stats::quantile(lr, probs[1], names = FALSE)),
      hi  = as.numeric(stats::quantile(lr, probs[2], names = FALSE))
    )
  }
  
  safe_abs <- function(x, probs = c(0.025, 0.975)) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
    c(
      med = median(x),
      lo  = as.numeric(stats::quantile(x, probs[1], names = FALSE)),
      hi  = as.numeric(stats::quantile(x, probs[2], names = FALSE))
    )
  }
  
  r13_24  <- safe_ratio("q13", "q24")
  lr13_24 <- safe_logratio("q13", "q24")
  
  r31_42  <- safe_ratio("q31", "q42")
  lr31_42 <- safe_logratio("q31", "q42")
  
  r34_12  <- safe_ratio("q34", "q12")
  lr34_12 <- safe_logratio("q34", "q12")
  
  r43_21  <- safe_ratio("q43", "q21")
  lr43_21 <- safe_logratio("q43", "q21")
  
  abs_q12 <- if ("q12" %in% names(df)) safe_abs(df$q12) else c(med = NA, lo = NA, hi = NA)
  abs_q13 <- if ("q13" %in% names(df)) safe_abs(df$q13) else c(med = NA, lo = NA, hi = NA)
  abs_q21 <- if ("q21" %in% names(df)) safe_abs(df$q21) else c(med = NA, lo = NA, hi = NA)
  abs_q24 <- if ("q24" %in% names(df)) safe_abs(df$q24) else c(med = NA, lo = NA, hi = NA)
  abs_q31 <- if ("q31" %in% names(df)) safe_abs(df$q31) else c(med = NA, lo = NA, hi = NA)
  abs_q34 <- if ("q34" %in% names(df)) safe_abs(df$q34) else c(med = NA, lo = NA, hi = NA)
  abs_q42 <- if ("q42" %in% names(df)) safe_abs(df$q42) else c(med = NA, lo = NA, hi = NA)
  abs_q43 <- if ("q43" %in% names(df)) safe_abs(df$q43) else c(med = NA, lo = NA, hi = NA)
  
  tibble(
    n_draws = nrow(df),
    
    cov_q12 = cover["q12"], cov_q13 = cover["q13"], cov_q21 = cover["q21"], cov_q24 = cover["q24"],
    cov_q31 = cover["q31"], cov_q34 = cover["q34"], cov_q42 = cover["q42"], cov_q43 = cover["q43"],
    
    P_q13_gt_q24 = safe_prob("q13", "q24"),
    P_q31_gt_q42 = safe_prob("q31", "q42"),
    P_q34_gt_q12 = safe_prob("q34", "q12"),
    P_q43_gt_q21 = safe_prob("q43", "q21"),
    
    median_ratio_q13_q24 = unname(r13_24["med"]),
    lo_ratio_q13_q24     = unname(r13_24["lo"]),
    hi_ratio_q13_q24     = unname(r13_24["hi"]),
    
    median_ratio_q31_q42 = unname(r31_42["med"]),
    lo_ratio_q31_q42     = unname(r31_42["lo"]),
    hi_ratio_q31_q42     = unname(r31_42["hi"]),
    
    median_ratio_q34_q12 = unname(r34_12["med"]),
    lo_ratio_q34_q12     = unname(r34_12["lo"]),
    hi_ratio_q34_q12     = unname(r34_12["hi"]),
    
    median_ratio_q43_q21 = unname(r43_21["med"]),
    lo_ratio_q43_q21     = unname(r43_21["lo"]),
    hi_ratio_q43_q21     = unname(r43_21["hi"]),
    
    median_logratio_q13_q24 = unname(lr13_24["med"]),
    lo_logratio_q13_q24     = unname(lr13_24["lo"]),
    hi_logratio_q13_q24     = unname(lr13_24["hi"]),
    
    median_logratio_q31_q42 = unname(lr31_42["med"]),
    lo_logratio_q31_q42     = unname(lr31_42["lo"]),
    hi_logratio_q31_q42     = unname(lr31_42["hi"]),
    
    median_logratio_q34_q12 = unname(lr34_12["med"]),
    lo_logratio_q34_q12     = unname(lr34_12["lo"]),
    hi_logratio_q34_q12     = unname(lr34_12["hi"]),
    
    median_logratio_q43_q21 = unname(lr43_21["med"]),
    lo_logratio_q43_q21     = unname(lr43_21["lo"]),
    hi_logratio_q43_q21     = unname(lr43_21["hi"]),
    
    median_q12 = unname(abs_q12["med"]), lo_q12 = unname(abs_q12["lo"]), hi_q12 = unname(abs_q12["hi"]),
    median_q13 = unname(abs_q13["med"]), lo_q13 = unname(abs_q13["lo"]), hi_q13 = unname(abs_q13["hi"]),
    median_q21 = unname(abs_q21["med"]), lo_q21 = unname(abs_q21["lo"]), hi_q21 = unname(abs_q21["hi"]),
    median_q24 = unname(abs_q24["med"]), lo_q24 = unname(abs_q24["lo"]), hi_q24 = unname(abs_q24["hi"]),
    median_q31 = unname(abs_q31["med"]), lo_q31 = unname(abs_q31["lo"]), hi_q31 = unname(abs_q31["hi"]),
    median_q34 = unname(abs_q34["med"]), lo_q34 = unname(abs_q34["lo"]), hi_q34 = unname(abs_q34["hi"]),
    median_q42 = unname(abs_q42["med"]), lo_q42 = unname(abs_q42["lo"]), hi_q42 = unname(abs_q42["hi"]),
    median_q43 = unname(abs_q43["med"]), lo_q43 = unname(abs_q43["lo"]), hi_q43 = unname(abs_q43["hi"])
  )
}

dir_tbl <- log_draws_tbl %>%
  group_by(analysis, thr_label, model) %>%
  group_modify(~{
    out <- summarise_q_direction(.x)
    if (is.null(out)) return(tibble())
    out
  }) %>%
  ungroup()

write_csv(dir_tbl, file.path(out_dir, "12_directionality_checks_by_analysis_threshold_and_model.csv"))

precedence_tbl <- dir_tbl %>%
  mutate(
    wood_gain_faster_under_stress =
      case_when(
        analysis %in% c("drought", "frost") ~ median_logratio_q13_q24 < 0,
        TRUE ~ NA
      ),
    wood_gain_faster_under_reduced_stress =
      case_when(
        analysis %in% c("drought", "frost") ~ median_logratio_q13_q24 > 0,
        TRUE ~ NA
      ),
    env_shift_more_herbaceous_than_woody =
      case_when(
        analysis == "drought" ~ median_logratio_q34_q12 < 0,
        analysis == "frost"   ~ median_logratio_q43_q21 < 0,
        TRUE ~ NA
      ),
    heuristic_supports_env_before_wood_gain =
      case_when(
        analysis == "drought" ~ (median_logratio_q13_q24 < 0) & (median_logratio_q34_q12 < 0),
        analysis == "frost"   ~ (median_logratio_q13_q24 > 0) & (median_logratio_q43_q21 < 0),
        TRUE ~ NA
      )
  )

write_csv(precedence_tbl, file.path(out_dir, "12_precedence_heuristic_summary.csv"))

thr_lookup <- thr_tbl_all %>%
  dplyr::select(analysis, thr_label, thr_type, percentile, thr_value) %>%
  mutate(thr_value = signif(thr_value, 6))

analysis_labels <- function(analysis) {
  if (analysis == "drought") {
    return(list(
      reduced_stress = "reduced drought stress",
      stress = "strong drought stress",
      var = "MCWD"
    ))
  }
  if (analysis == "frost") {
    return(list(
      reduced_stress = "frost-free conditions",
      stress = "frost",
      var = "BIO6"
    ))
  }
  list(reduced_stress = "reduced stress", stress = "stress", var = "env")
}

fmt_thr_label_pretty <- function(thr_label) {
  if (thr_label == "thr_0C") return("0°C")
  if (grepl("^thr_p", thr_label)) return(sub("^thr_p", "p", thr_label))
  thr_label
}


## 13. Summaries + main-text figure + SupMat histograms ----

## 13.1 Bayes factor summaries (Dep vs Ind; RJ vs Dep)
logbf_rj_dep <- stones_tbl %>%
  group_by(analysis, thr_label) %>%
  group_modify(~calc_logbf_pairs(.x, "dependent_RJ", "dependent")) %>%
  ungroup()

write_csv(logbf_rj_dep, file.path(out_dir, "13_logBF_depRJ_vs_dep_by_threshold_and_analysis.csv"))

bf_summary <- bind_rows(
  logbf_dep_ind %>% mutate(comparison = "dependent_vs_independent"),
  logbf_rj_dep  %>% mutate(comparison = "dependent_RJ_vs_dependent")
) %>%
  group_by(analysis, thr_label, comparison) %>%
  summarise(
    mean_logBF = mean(logBF, na.rm = TRUE),
    sd_logBF = sd(logBF, na.rm = TRUE),
    median_logBF = median(logBF, na.rm = TRUE),
    lo_logBF = quantile(logBF, 0.025, na.rm = TRUE),
    hi_logBF = quantile(logBF, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(bf_summary, file.path(out_dir, "13_bayesfactor_summary_by_threshold_and_analysis.csv"))

## 13.2 Full posterior summaries of all q’s (absolute med/CI) for SupMat tables
summarise_q_posteriors <- function(df) {
  need <- q_params
  present <- intersect(need, names(df))
  if (length(present) == 0) return(NULL)
  
  incl <- vapply(need, function(p) {
    if (!p %in% names(df)) return(0)
    mean(is.finite(df[[p]]), na.rm = TRUE)
  }, numeric(1))
  
  safe_med <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else median(x) }
  safe_q   <- function(x, pr) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else as.numeric(quantile(x, pr)) }
  
  tibble(
    param = need,
    inclusion = incl,
    median = vapply(need, \(p) if (p %in% names(df)) safe_med(df[[p]]) else NA_real_, numeric(1)),
    lo     = vapply(need, \(p) if (p %in% names(df)) safe_q(df[[p]], 0.025) else NA_real_, numeric(1)),
    hi     = vapply(need, \(p) if (p %in% names(df)) safe_q(df[[p]], 0.975) else NA_real_, numeric(1))
  )
}

q_post_tbl <- log_draws_tbl %>%
  filter(model %in% c("dependent", "dependent_RJ")) %>%
  group_by(analysis, thr_label, model) %>%
  group_modify(~{
    out <- summarise_q_posteriors(.x)
    if (is.null(out)) return(tibble())
    out
  }) %>%
  ungroup()

write_csv(q_post_tbl, file.path(out_dir, "13_q_posterior_summaries_all_q.csv"))

## 13.3 SupMat histogram grids for all q parameters (log10 rates for readability)
make_q_long <- function(df, models_keep = c("dependent", "dependent_RJ")) {
  df %>%
    filter(model %in% models_keep) %>%
    dplyr::select(analysis, thr_label, model, any_of(q_params)) %>%
    pivot_longer(cols = any_of(q_params), names_to = "param", values_to = "value") %>%
    filter(is.finite(value), value > 0) %>%
    mutate(value_log10 = log10(value))
}

q_long_all <- make_q_long(log_draws_tbl, models_keep = c("dependent", "dependent_RJ"))

p_q_hists_log <- ggplot(q_long_all, aes(x = value_log10)) +
  geom_histogram(bins = 45) +
  facet_grid(param ~ analysis + thr_label + model, scales = "free_x") +
  labs(
    x = "log10(posterior rate)",
    y = "Count",
    title = "Posterior distributions of all transition rates (q) across analyses, thresholds, and models",
    subtitle = "Shown on log10 scale for comparability; RJ kept as sensitivity analysis"
  )

ggsave(
  file.path(out_dir, "13C_supmat_all_q_log10_histograms.png"),
  p_q_hists_log,
  width = 15,
  height = 10.5,
  dpi = 220
)

ggsave(
  file.path(out_dir, "13C_supmat_all_q_log10_histograms.pdf"),
  p_q_hists_log,
  width = 15,
  height = 10.5,
  dpi = 220
)

## 13.4 Main-text figure: stacked drought + frost panels ----
## Shared sign convention:
##   x = ln(rate under reduced stress / rate under increased stress)
##   x < 0 => faster rates under increased stress
##   x > 0 => faster rates under reduced stress

if (!requireNamespace("scales", quietly = TRUE)) {
  stop("[ERROR] Package 'scales' is required.")
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("[ERROR] Package 'cowplot' is required.")
}

alpha_by_thr <- function(thr) {
  if (thr == "thr_p50") return(0.95)
  if (thr == "thr_0C")  return(0.85)
  if (thr == "thr_p60") return(0.35)
  return(0.18)
}

fmt_thr_label_pretty <- function(thr_label) {
  if (thr_label == "thr_0C") return("0°C")
  if (grepl("^thr_p", thr_label)) return(sub("^thr_p", "p", thr_label))
  thr_label
}

analysis_xlims <- function(analysis) {
  c(-3, 3.5)
}

make_main_hist_panel <- function(analysis_use) {
  
  xlim_use <- analysis_xlims(analysis_use)
  stress_word <- if (analysis_use == "drought") "drought" else "frost"
  
  thr_levels_this <- thr_tbl_all %>%
    dplyr::filter(analysis == analysis_use) %>%
    dplyr::arrange(
      dplyr::case_when(thr_label == "thr_0C" ~ 0L, TRUE ~ 1L),
      percentile
    ) %>%
    dplyr::pull(thr_label) %>%
    as.character()
  
  main_hr <- log_draws_tbl %>%
    dplyr::filter(model == "dependent", analysis == !!analysis_use) %>%
    dplyr::transmute(
      thr_label,
      `H→W shift rate`    = log(q13 / q24),
      `W→H reversal rate` = log(q31 / q42)
    ) %>%
    tidyr::pivot_longer(
      cols = c(`H→W shift rate`, `W→H reversal rate`),
      names_to = "direction",
      values_to = "logratio"
    ) %>%
    dplyr::filter(is.finite(logratio)) %>%
    dplyr::mutate(
      thr_label = factor(thr_label, levels = thr_levels_this),
      direction = factor(
        direction,
        levels = c("H→W shift rate", "W→H reversal rate")
      )
    )
  
  alpha_values <- setNames(
    vapply(thr_levels_this, alpha_by_thr, numeric(1)),
    thr_levels_this
  )
  
  alpha_labels <- setNames(
    vapply(thr_levels_this, fmt_thr_label_pretty, character(1)),
    thr_levels_this
  )
  
  arrow_df <- tibble::tibble(
    x0 = xlim_use[1] + 0.15,
    x1 = xlim_use[2] - 0.15,
    left_label  = paste0("faster rates under\nincreased ", stress_word, " stress"),
    right_label = paste0("faster rates under\nreduced ", stress_word, " stress")
  )
  
  ratio_box_df <- tibble::tibble(
    direction = c("H→W shift rate", "W→H reversal rate"),
    lab = c("q13/q24", "q31/q42"),
    x = xlim_use[2] - 0.15,
    y = c(Inf, Inf)
  )
  
  ggplot(
    main_hr,
    aes(x = logratio, fill = direction, alpha = thr_label)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_histogram(
      bins = 60,
      position = "identity",
      colour = NA
    ) +
    geom_segment(
      data = arrow_df,
      aes(x = x0, xend = x1, y = -Inf, yend = -Inf),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.10, "inches"), ends = "both"),
      linewidth = 0.35
    ) +
    geom_text(
      data = arrow_df,
      aes(x = x0 + 0.15, y = -Inf, label = left_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 2.0,
      size = 2.9
    ) +
    geom_text(
      data = arrow_df,
      aes(x = x1 - 0.15, y = -Inf, label = right_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 2.0,
      size = 2.9
    ) +
    geom_label(
      data = ratio_box_df,
      aes(x = x, y = y, label = lab, fill = direction),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = c(1.2, 2.7),
      size = 3.0,
      label.size = 0,
      alpha = 0.85
    ) +
    coord_cartesian(xlim = xlim_use, clip = "off") +
    scale_fill_manual(
      values = c(
        `H→W shift rate`    = "#892255",
        `W→H reversal rate` = "#44AA9A"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = seq(-3, 3, by = 1)
    ) +
    scale_alpha_manual(
      values = alpha_values,
      name = "Threshold",
      labels = alpha_labels
    ) +
    labs(
      title = stress_word,
      x = "ln(rate under reduced stress / rate under increased stress)",
      y = "Posterior draws"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      legend.background = element_rect(
        fill = scales::alpha("white", 0.75),
        colour = NA
      ),
      legend.key = element_rect(fill = NA, colour = NA),
      plot.margin = grid::unit(c(4, 8, 28, 4), "pt")
    )
}

p_drought <- make_main_hist_panel("drought")
p_frost   <- make_main_hist_panel("frost")

p_combined <- cowplot::plot_grid(
  p_drought,
  p_frost,
  ncol = 1,
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 16,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

ggsave(
  file.path(out_dir, "13A_maintext_two_panel_histograms_logratio_combined.png"),
  p_combined,
  width = 8.2,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(out_dir, "13A_maintext_two_panel_histograms_logratio_combined.pdf"),
  p_combined,
  width = 8.2,
  height = 7,
  dpi = 300
)

message("[OK] Step 13.4 wrote stacked Figure 6 to ", normalizePath(out_dir))



