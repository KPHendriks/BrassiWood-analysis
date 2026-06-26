  #!/usr/bin/env Rscript
  
  ## 2f.StudyWoodinessEvolution.R
  ## Phylogenetic woodiness evolution in Brassicaceae (+ Cleomaceae for context)
  ## using the LSD2-calibrated BrassiToL species tree.  ##
  ## Goals
  ##   Quantify and visualize growth-form shifts (H <-> W) on the dated tree:
  ##     2) Phylogenetic signal (Pagel’s lambda; Fritz & Purvis’ D)
  ##     3) ASR under Mk (ER vs ARD model test + posterior node probs)
  ##     4) SIMMAP posterior (event counts + timings)
  ##     5) Localization-aware synthesis + key figures (circular tree + timing plots)
  ##     6) Tempo + transition propensity through time (SIMMAP-derived)
  ##
  ## Notes on terminology
  ##   We use “shift” (not “switch”) for H<->W transitions.
  ##
  ## Outputs
  ##   WP2_BrassiNiche/results_final/3_woodiness_evolution/
  
  ## 0. Prepare -------------------------------------------------------------------

out_dir <- "WP2_BrassiNiche/results_final/3_woodiness_evolution"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(treeio)
  library(tidytree)
  library(phytools)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(readr)
  library(purrr)
  library(progressr)
})

message("[INFO] Output dir: ", normalizePath(out_dir, mustWork = FALSE))

## 0.1 Inputs: calibrated species tree -----------------------------------------

cli <- list()

cli$project <- "2026-03-09_BrassiWood_v2"
cli$calibration_source <- "lsd2"
cli$lsd2_run <- "lsd2_top200_fixedtopo"

cli$tree_lsd2_nexus <- file.path(
  "WP1_BrassiToL/results_final",
  cli$project,
  "8_results_species_tree_calibrated",
  cli$lsd2_run,
  "05_lsd2_dating",
  "2026-03-09_BrassiWood_v2_lsd2_top200_dated.date.nexus"
)

cli$species_csv <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
cli$specimen_csv <- "WP2_BrassiNiche/data/brassiwood_specimens_publication.csv"

stopifnot(file.exists(cli$species_csv))
stopifnot(file.exists(cli$specimen_csv))

if (!file.exists(cli$tree_lsd2_nexus)) {
  stop("[ERROR] LSD2 Nexus tree not found: ", cli$tree_lsd2_nexus)
}

## 0.2 Global colors -------------------------------------------------------------

clade_colors <- c(
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8",
  "Aethionemeae"      = "#ffed57",
  "Cleomaceae"        = "grey70"
)

col_hw <- "#892255"  # H>W
col_wh <- "#44AA9A"  # W>H

q025 <- function(x) {
  stats::quantile(x, 0.025, names = FALSE, na.rm = TRUE)
}
q975 <- function(x) {
  stats::quantile(x, 0.975, names = FALSE, na.rm = TRUE)
}

## 0.3 Helper functions -------------------------------------------------------------

safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) NA_real_ else x[1]
}

extract_lambda_estimate <- function(fit_lambda) {
  cand <- c(
    fit_lambda$opt$lambda,
    fit_lambda$opt$pars,
    fit_lambda$opt$par,
    fit_lambda$lambda
  )
  cand <- unlist(cand, recursive = TRUE, use.names = TRUE)
  cand_num <- suppressWarnings(as.numeric(cand))
  cand_num <- cand_num[is.finite(cand_num)]
  
  if (length(cand_num) == 0) return(NA_real_)
  
  ## Prefer values in the biologically sensible lambda range first
  in_range <- cand_num[cand_num >= 0 & cand_num <= 1.5]
  if (length(in_range) > 0) return(in_range[1])
  
  cand_num[1]
}

get_node_ages_bp <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)
  H <- max(depth)
  H - depth
}

interpret_lambda <- function(lambda_est, p_lrt) {
  if (!is.finite(lambda_est) || !is.finite(p_lrt)) {
    return("Could not determine support for phylogenetic signal.")
  }
  if (p_lrt >= 0.05) {
    return("Little support that the lambda-transformed model improves fit over the non-transformed model.")
  }
  if (lambda_est < 0.2) {
    return("Trait shows weak phylogenetic dependence.")
  }
  if (lambda_est < 0.8) {
    return("Trait shows moderate phylogenetic dependence.")
  }
  return("Trait shows strong phylogenetic dependence.")
}

interpret_D <- function(D_est, p_rand, p_bm) {
  if (!is.finite(D_est)) {
    return("D could not be estimated.")
  }
  
  pieces <- c()
  
  if (is.finite(p_rand)) {
    pieces <- c(
      pieces,
      if (p_rand < 0.05) {
        "distribution differs from random expectation (D = 1)"
      } else {
        "distribution is not distinguishable from random expectation (D = 1)"
      }
    )
  }
  
  if (is.finite(p_bm)) {
    pieces <- c(
      pieces,
      if (p_bm < 0.05) {
        "distribution differs from Brownian expectation (D = 0)"
      } else {
        "distribution is not distinguishable from Brownian expectation (D = 0)"
      }
    )
  }
  
  d_class <- dplyr::case_when(
    D_est < 0   ~ "more clumped than Brownian",
    D_est < 0.5 ~ "phylogenetically clumped",
    D_est < 1   ~ "intermediate",
    TRUE        ~ "overdispersed relative to Brownian"
  )
  
  paste0("D indicates a ", d_class, " distribution",
         if (length(pieces) > 0) paste0("; ", paste(pieces, collapse = "; "), ".") else ".")
}

make_ultrametric_by_extending_tips <- function(tree, rel_tol = 1e-5, min_tol = 1e-6) {
  tip_depths <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
  tree_height <- max(tip_depths, na.rm = TRUE)
  depth_span <- diff(range(tip_depths, na.rm = TRUE))
  tol_abs <- max(min_tol, tree_height * rel_tol)
  
  message("[INFO] Tree height: ", signif(tree_height, 6), " Ma")
  message("[INFO] Tip-depth span: ", signif(depth_span, 6), " Ma")
  message("[INFO] Ultrametric tolerance: ", signif(tol_abs, 6), " Ma")
  
  if (depth_span > tol_abs) {
    stop(
      "[ERROR] Tree is meaningfully non-ultrametric; tip-depth span = ",
      signif(depth_span, 6),
      " Ma, tolerance = ",
      signif(tol_abs, 6),
      " Ma."
    )
  }
  
  max_depth <- max(tip_depths, na.rm = TRUE)
  
  terminal_edges <- which(tree$edge[, 2] <= length(tree$tip.label))
  terminal_tips  <- tree$edge[terminal_edges, 2]
  
  add_len <- max_depth - tip_depths[terminal_tips]
  add_len[add_len < 0] <- 0
  
  if (sum(add_len > 0, na.rm = TRUE) > 0) {
    message("[INFO] Extending ",
            sum(add_len > 0, na.rm = TRUE),
            " terminal edges by <= ",
            signif(max(add_len, na.rm = TRUE), 6),
            " Ma.")
    
    tree$edge.length[terminal_edges] <-
      tree$edge.length[terminal_edges] + add_len
  } else {
    message("[INFO] No terminal-edge correction needed.")
  }
  
  tree
}


## 1.1 Read tree from LSD2 Nexus output -----------------------------------------
## Why: the plain .nwk may retain scaled branch lengths, whereas the Nexus tree
## contains the dated tree used for plotting.

if (!file.exists(cli$tree_lsd2_nexus)) {
  stop("[ERROR] LSD2 Nexus tree not found: ", cli$tree_lsd2_nexus)
}

tre <- ape::read.nexus(cli$tree_lsd2_nexus)

if (inherits(tre, "multiPhylo")) {
  message("[INFO] Nexus file contained multiple trees; using the first one.")
  tre <- tre[[1]]
}

message("[INFO] Read LSD2-calibrated tree from Nexus:")
message("[INFO]   tree: ", normalizePath(cli$tree_lsd2_nexus))
message("[INFO]   tips=", length(tre$tip.label), " ; Nnode=", tre$Nnode)

## Quick scale diagnostic
tip_depths_raw <- ape::node.depth.edgelength(tre)[seq_along(tre$tip.label)]
tree_height_raw <- max(tip_depths_raw, na.rm = TRUE)

message("[INFO] Raw tree height from Nexus: ",
        signif(tree_height_raw, 6), " Ma")

## 1.2 Read metadata from local publication CSVs
## Why: specimen table provides tip-level metadata; species table provides growth form.

message("[INFO] Reading specimen metadata from local CSV: ", cli$specimen_csv)
specimen_details <- readr::read_csv(cli$specimen_csv, show_col_types = FALSE)

message("[INFO] Reading species metadata from local CSV: ", cli$species_csv)
species_details <- readr::read_csv(cli$species_csv, show_col_types = FALSE)
         

## 1.4 Minimal cleaning (retain only what downstream steps need)
## Why: enforce types and create a consistent woody flag (H/W) for each sampled tip.
specimen_details <- specimen_details %>%
  dplyr::mutate(
    SAMPLE             = as.character(SAMPLE),
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    TRIBE              = as.character(TRIBE),
    SUPERTRIBE         = as.character(SUPERTRIBE),
    FAMILY             = as.character(FAMILY),
    GENUS              = as.character(GENUS),
    Type_status        = as.character(Type_status),
    loci_remaining     = suppressWarnings(as.numeric(loci_remaining)),
    tip_label          = paste0(TRIBE, "_", SPECIES_NAME_PRINT, " (", SAMPLE, ")")
  )
specimen_details$loci_remaining[is.na(specimen_details$loci_remaining)] <- 0

normalize_growth <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  x <- stringr::str_to_upper(x)
  dplyr::case_when(x %in% c("W", "W/L", "W/T") ~ "W",
                   x == "H"                   ~ "H",
                   TRUE                       ~ NA_character_)
}

## Lookup by print name
species_lu_print <- species_details %>%
  dplyr::transmute(
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    GROWTH_FORM        = normalize_growth(GROWTH_FORM)
  ) %>%
  dplyr::filter(!is.na(SPECIES_NAME_PRINT)) %>%
  dplyr::group_by(SPECIES_NAME_PRINT) %>%
  dplyr::summarise(
    has_W = any(GROWTH_FORM == "W", na.rm = TRUE),
    has_H = any(GROWTH_FORM == "H", na.rm = TRUE),
    WOODY_STATE_PRINT = dplyr::case_when(has_W & !has_H ~ "W", has_H &
                                           !has_W ~ "H", TRUE           ~ NA_character_),
    .groups = "drop"
  ) %>%
  dplyr::select(SPECIES_NAME_PRINT, WOODY_STATE_PRINT)

## Fallback lookup by binomial
species_lu_binom <- species_details %>%
  dplyr::transmute(
    BIONOMIAL   = stringr::str_squish(as.character(BIONOMIAL)),
    GROWTH_FORM = normalize_growth(GROWTH_FORM)
  ) %>%
  dplyr::filter(!is.na(BIONOMIAL)) %>%
  dplyr::group_by(BIONOMIAL) %>%
  dplyr::summarise(
    has_W = any(GROWTH_FORM == "W", na.rm = TRUE),
    has_H = any(GROWTH_FORM == "H", na.rm = TRUE),
    WOODY_STATE_BINOM = dplyr::case_when(has_W & !has_H ~ "W", has_H &
                                           !has_W ~ "H", TRUE           ~ NA_character_),
    .groups = "drop"
  ) %>%
  dplyr::select(BIONOMIAL, WOODY_STATE_BINOM)

specimen_details <- specimen_details %>%
  dplyr::mutate(species_binom = stringr::str_squish(
    sub("^([A-Za-z]+\\s+[A-Za-z-]+).*$", "\\1", SPECIES_NAME_PRINT)
  )) %>%
  dplyr::left_join(species_lu_print, by = "SPECIES_NAME_PRINT") %>%
  dplyr::left_join(species_lu_binom, by = c("species_binom" = "BIONOMIAL")) %>%
  dplyr::mutate(
    WOODY_STATE = dplyr::coalesce(WOODY_STATE_PRINT, WOODY_STATE_BINOM),
    woody_flag  = dplyr::case_when(
      WOODY_STATE == "W" ~ "W",
      WOODY_STATE == "H" ~ "H",
      TRUE               ~ NA_character_
    ),
    woody_source = dplyr::case_when(
      !is.na(WOODY_STATE_PRINT) ~ "print",
      is.na(WOODY_STATE_PRINT) &
        !is.na(WOODY_STATE_BINOM) ~ "binom",
      TRUE ~ "none"
    )
  )

message("[INFO] Woodiness source breakdown:")
print(table(specimen_details$woody_source, useNA = "ifany"))
message("[INFO] Still missing woody_flag: ", sum(is.na(specimen_details$woody_flag)))

## Tip info aligned to tree tips (used later)
tip_info <- specimen_details %>%
  dplyr::select(SAMPLE, SPECIES_NAME_PRINT, FAMILY, TRIBE, SUPERTRIBE, GENUS) %>%
  dplyr::distinct() %>%
  dplyr::filter(SAMPLE %in% tre$tip.label)


## 2. Phylogenetic signal -------------------------------------------------------

message("[INFO] Step 2: Phylogenetic signal & diagnostics (coded tree)")

## 2.1 Drop known errors/duplicates (but keep key woody island taxa for niche consistency)
## Why: ensure our trait inference is based on the intended final phylogenetic sampling.
errors <- c(
  )

## Keep overrides: ensure Lobularia canariensis and Brassica oleracea woody island taxa remain included.
## Why: the non-phylogenetic niche study includes these; we want consistency between pipelines.
keep_override_samples <- specimen_details %>%
  dplyr::filter(
    !is.na(SAMPLE),
    grepl("^Lobularia canariensis\\b", SPECIES_NAME_PRINT) |
      grepl("^Brassica oleracea\\b", SPECIES_NAME_PRINT)
  ) %>%
  dplyr::pull(SAMPLE) %>%
  unique()

errors_drop <- setdiff(errors, keep_override_samples)

if (length(intersect(errors, keep_override_samples)) > 0) {
  message(
    "[INFO] Override: retaining these samples even if listed as errors: ",
    paste(intersect(errors, keep_override_samples), collapse = ", ")
  )
}

tree0 <- ape::drop.tip(tre, intersect(tre$tip.label, errors_drop))
message("[INFO] After error-drop (with overrides): tips=",
        length(tree0$tip.label))

## 2.2 Check rooting, branch lengths and ultrametricity
## Why: LSD2 tree is already dated; edge lengths are time in Ma.
## We avoid force.ultrametric(), but allow tiny numerical deviations.

if (!ape::is.rooted(tree0)) {
  warning("[WARN] Tree is not rooted. Check LSD2 output/rooting before downstream analyses.")
} else {
  message("[INFO] Tree is rooted.")
}

if (any(is.na(tree0$edge.length))) {
  stop("[ERROR] Tree has NA edge lengths.")
}

if (any(tree0$edge.length < 0, na.rm = TRUE)) {
  neg_n <- sum(tree0$edge.length < 0, na.rm = TRUE)
  stop("[ERROR] Tree has ", neg_n, " negative edge lengths. Inspect LSD2 output.")
}

## Replace exact zero-length edges by an extremely small value.
## This avoids downstream problems while preserving dated branch lengths.
zero_edges <- sum(tree0$edge.length == 0, na.rm = TRUE)

if (zero_edges > 0) {
  message("[INFO] Replacing ", zero_edges, " zero-length edges with 1e-10 Ma.")
  tree0$edge.length[tree0$edge.length == 0] <- 1e-10
}

## Check how close the tree is to ultrametric
tip_depths <- ape::node.depth.edgelength(tree0)[seq_along(tree0$tip.label)]
depth_range <- range(tip_depths, na.rm = TRUE)
depth_span  <- diff(depth_range)
tree_height <- max(tip_depths, na.rm = TRUE)

message("[INFO] Tree height: ", signif(tree_height, 6), " Ma")
message("[INFO] Tip-depth range: ",
        signif(depth_range[1], 6), "–", signif(depth_range[2], 6), " Ma")
message("[INFO] Tip-depth span: ", signif(depth_span, 6), " Ma")

## Allow tiny absolute/relative deviations from ultrametricity.
## For a 90 Ma tree, 1e-5 of tree height = 0.0009 Ma.
ultra_tol_abs <- max(1e-6, tree_height * 1e-5)

message("[INFO] Ultrametric tolerance: ",
        signif(ultra_tol_abs, 6), " Ma")

if (!ape::is.ultrametric(tree0, tol = ultra_tol_abs)) {
  
  if (depth_span <= ultra_tol_abs) {
    message("[INFO] Tree is slightly non-ultrametric within tolerance; applying terminal-edge correction.")
    
    max_depth <- max(tip_depths, na.rm = TRUE)
    
    terminal_edges <- which(tree0$edge[, 2] <= length(tree0$tip.label))
    terminal_tips  <- tree0$edge[terminal_edges, 2]
    
    add_len <- max_depth - tip_depths[terminal_tips]
    add_len[add_len < 0] <- 0
    
    tree0$edge.length[terminal_edges] <-
      tree0$edge.length[terminal_edges] + add_len
    
  } else {
    stop(
      "[ERROR] Tree is meaningfully non-ultrametric; tip-depth span = ",
      signif(depth_span, 6),
      " Ma. Inspect the LSD2 output before correcting the tree."
    )
  }
  
} else {
  message("[INFO] Tree is ultrametric within tolerance.")
}

## Recheck after optional correction
tip_depths2 <- ape::node.depth.edgelength(tree0)[seq_along(tree0$tip.label)]
depth_range2 <- range(tip_depths2, na.rm = TRUE)
depth_span2  <- diff(depth_range2)

message("[INFO] Final tip-depth range: ",
        signif(depth_range2[1], 6), "–", signif(depth_range2[2], 6), " Ma")
message("[INFO] Final tip-depth span: ", signif(depth_span2, 6), " Ma")

if (depth_span2 > ultra_tol_abs) {
  stop("[ERROR] Tree still exceeds ultrametric tolerance after correction.")
}

message("[INFO] Tree passed ultrametric check.")

## 2.3 Build woody vector (0/1/NA) aligned to tips
## Why: consistent binary coding is required for fitDiscrete/fitMk/simmap.
woody_df <- specimen_details %>%
  dplyr::select(SAMPLE, woody_flag) %>%
  dplyr::filter(SAMPLE %in% tree0$tip.label)

woody_vec <- dplyr::case_when(woody_df$woody_flag == "W" ~ 1,
                              woody_df$woody_flag == "H" ~ 0,
                              TRUE                       ~ NA_real_)
names(woody_vec) <- woody_df$SAMPLE
woody_vec <- woody_vec[tree0$tip.label]

message("[INFO] woody_vec summary (0=H, 1=W):")
print(table(woody_vec, useNA = "ifany"))

missing_spp <- specimen_details %>%
  dplyr::filter(SAMPLE %in% names(woody_vec)[is.na(woody_vec)]) %>%
  dplyr::pull(SPECIES_NAME_PRINT) %>%
  unique()
if (length(missing_spp) > 0) {
  message("[INFO] Define growth form (W or H) for these species to increase coded sample size:")
  print(missing_spp)
}

## Keep coded tips for model-based steps
keep_tips_coded <- names(woody_vec)[!is.na(woody_vec)]
tree_coded <- ape::keep.tip(tree0, keep_tips_coded)

message("[INFO] Checking ultrametricity of coded tree after dropping uncoded tips.")
tree_coded <- make_ultrametric_by_extending_tips(tree_coded)

x_num <- woody_vec[tree_coded$tip.label]

x_fac <- factor(ifelse(x_num == 1, "1", "0"), levels = c("0", "1"))
names(x_fac) <- tree_coded$tip.label

## 2.4 Pagel’s lambda (discrete ER + lambda transform)
## Why: quantify phylogenetic dependence of growth form on the tree.
suppressPackageStartupMessages(library(geiger))

fit_er     <- geiger::fitDiscrete(tree_coded, x_fac, model = "ER", transform = "none")
fit_lambda <- geiger::fitDiscrete(tree_coded, x_fac, model = "ER", transform = "lambda")

LR <- 2 * (fit_lambda$opt$lnL - fit_er$opt$lnL)
p_lrt <- stats::pchisq(LR, df = 1, lower.tail = FALSE)

lambda_est <- extract_lambda_estimate(fit_lambda)

n_coded <- length(tree_coded$tip.label)

pagel_model_tbl <- tibble::tibble(
  analysis = "Pagel_lambda_binary_trait",
  model = c("ER_no_transform", "ER_lambda"),
  n_tips_coded = n_coded,
  logLik = c(
    safe_num(fit_er$opt$lnL),
    safe_num(fit_lambda$opt$lnL)
  ),
  k = c(1L, 2L),
  transform = c("none", "lambda"),
  lambda_estimate = c(NA_real_, lambda_est)
) %>%
  dplyr::mutate(
    AIC = 2 * k - 2 * logLik,
    deltaAIC = AIC - min(AIC, na.rm = TRUE),
    AIC_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
  )

pagel_test_tbl <- tibble::tibble(
  analysis = "Pagel_lambda_binary_trait",
  comparison = "ER_lambda_vs_ER_no_transform",
  n_tips_coded = n_coded,
  logLik_null = safe_num(fit_er$opt$lnL),
  logLik_alt = safe_num(fit_lambda$opt$lnL),
  lambda_estimate = lambda_est,
  LR = LR,
  df = 1L,
  p_value = p_lrt,
  interpretation = interpret_lambda(lambda_est, p_lrt)
)

write_csv(
  pagel_model_tbl,
  file.path(out_dir, "2.3_pagel_lambda_model_fits.csv")
)

write_csv(
  pagel_test_tbl,
  file.path(out_dir, "2.3_pagel_lambda_summary.csv")
)

message("[OK] Wrote Pagel lambda tables.")


## 2.5 Fritz & Purvis’ D (optional)
## Why: alternative phylogenetic signal summary for binary traits.
d_res <- NULL
if (requireNamespace("caper", quietly = TRUE)) {
  suppressPackageStartupMessages(library(caper))
  dat_fp <- data.frame(species = tree0$tip.label, woody   = woody_vec[tree0$tip.label])
  dat_fp2 <- dat_fp[!is.na(dat_fp$woody), , drop = FALSE]
  tree_d  <- ape::drop.tip(tree0, setdiff(tree0$tip.label, dat_fp2$species))
  
  comp_dat <- caper::comparative.data(
    phy       = tree_d,
    data      = dat_fp2,
    names.col = "species",
    vcv       = TRUE,
    na.omit   = TRUE
  )
  d_res <- caper::phylo.d(comp_dat, binvar = woody)
  
  ## Try to extract standard fields from phylo.d output
  D_est <- safe_num(d_res$DEstimate)
  
  ## In caper::phylo.d these are usually p-values against random and Brownian expectations.
  p_rand <- safe_num(d_res$Pval1)
  p_bm   <- safe_num(d_res$Pval0)
  
  d_summary_tbl <- tibble::tibble(
    analysis = "Fritz_Purvis_D_binary_trait",
    n_tips_coded = nrow(dat_fp2),
    D_estimate = D_est,
    p_value_vs_random_D1 = p_rand,
    p_value_vs_brownian_D0 = p_bm,
    interpretation = interpret_D(D_est, p_rand, p_bm)
  )
  
  write_csv(
    d_summary_tbl,
    file.path(out_dir, "2.4_fritz_purvis_D_summary.csv")
  )
  
  ## Optional: also save the raw printout for traceability
  write_lines(
    paste(capture.output(print(d_res)), collapse = "\n"),
    file.path(out_dir, "2.4_fritz_purvis_D_raw.txt")
  )
  
  message("[OK] Wrote Fritz & Purvis D tables.")
  message("[OK] Wrote Fritz & Purvis D results.")
} else {
  warning("[WARN] Package 'caper' not installed; skipping Fritz & Purvis D.")
}


# Write combined output table for publication
phylo_signal_summary_tbl <- tibble::tibble(
  metric = c("Pagel_lambda", "Fritz_Purvis_D"),
  n_tips_coded = c(
    n_coded,
    if (exists("dat_fp2")) nrow(dat_fp2) else NA_integer_
  ),
  estimate = c(
    lambda_est,
    if (exists("D_est")) D_est else NA_real_
  ),
  test_statistic = c(
    LR,
    NA_real_
  ),
  df = c(
    1L,
    NA_integer_
  ),
  p_value_primary = c(
    p_lrt,
    if (exists("p_rand")) p_rand else NA_real_
  ),
  p_value_secondary = c(
    NA_real_,
    if (exists("p_bm")) p_bm else NA_real_
  ),
  primary_label = c(
    "LRT_lambda_vs_no_transform",
    "P_vs_random_D1"
  ),
  secondary_label = c(
    NA_character_,
    "P_vs_brownian_D0"
  ),
  interpretation = c(
    interpret_lambda(lambda_est, p_lrt),
    if (exists("D_est")) interpret_D(D_est, p_rand, p_bm) else "D not calculated."
  )
)

write_csv(
  phylo_signal_summary_tbl,
  file.path(out_dir, "2.5_phylogenetic_signal_summary_for_supmat.csv")
)

message("[OK] Wrote combined phylogenetic signal summary table.")

## Save core objects (for reuse/debugging)
saveRDS(
  list(
    tree0 = tree0,
    tree_coded = tree_coded,
    woody_vec = woody_vec,
    x_fac = x_fac
  ),
  file.path(out_dir, "2_objects.rds")
)


## 3. Ancestral State Reconstruction --------------------------------------------

message("[INFO] Step 3: Mk ASR (ER vs ARD) + posterior nodes")

stopifnot(ape::is.ultrametric(tree_coded, tol = max(1e-6, max(ape::node.depth.edgelength(tree_coded)) * 1e-5)))

## 3.1 Fit Mk ER/ARD with a soft prior on a non-woody root
## Why: Mk inference can be sensitive to root state in highly asymmetric traits.
root_prior <- c("0" = 1 - 1e-8, "1" = 1e-8)

fit_mk_er  <- phytools::fitMk(
  tree = tree_coded,
  x = x_fac,
  model = "ER",
  pi = root_prior
)
fit_mk_ard <- phytools::fitMk(
  tree = tree_coded,
  x = x_fac,
  model = "ARD",
  pi = root_prior
)

assert_good_fit <- function(fit, label) {
  ll <- suppressWarnings(as.numeric(fit$logLik))
  if (!is.finite(ll) ||
      ll <= -1e40)
    stop("[ERROR] fitMk failed (", label, "): logLik=", ll)
  if (any(!is.finite(fit$rates)))
    stop("[ERROR] fitMk failed (", label, "): non-finite rates")
  invisible(TRUE)
}
assert_good_fit(fit_mk_er, "ER")
assert_good_fit(fit_mk_ard, "ARD")

## 3.2 Model comparison via IC (AIC/AICc)
get_ic <- function(fit, model, n_eff) {
  ll <- as.numeric(fit$logLik)
  k <- if (model == "ER")
    1L
  else
    2L
  AIC <- 2 * k - 2 * ll
  AICc <- if (n_eff <= (k + 1L))
    NA_real_
  else
    AIC + (2 * k * (k + 1)) / (n_eff - k - 1)
  tibble::tibble(
    logLik = ll,
    k = k,
    AIC = AIC,
    AICc = AICc,
    model = model
  )
}

mk_model_test <- dplyr::bind_rows(get_ic(fit_mk_er, "ER", n_coded),
                                  get_ic(fit_mk_ard, "ARD", n_coded)) %>%
  dplyr::mutate(
    IC      = dplyr::if_else(is.finite(AICc), AICc, AIC),
    deltaIC = IC - min(IC, na.rm = TRUE),
    weight  = exp(-0.5 * deltaIC) / sum(exp(-0.5 * deltaIC))
  ) %>%
  dplyr::arrange(IC) %>%
  dplyr::select(model, logLik, k, AIC, AICc, deltaIC, weight)

best_model <- mk_model_test$model[1]

fit_mk <- if (best_model == "ER") {
  fit_mk_er
} else {
  fit_mk_ard
}


write_csv(mk_model_test, file.path(out_dir, "3.0_mk_model_test.csv"))
message("[OK] Mk model test written; best_model=", best_model)

## 3.3 ASR posteriors at nodes (phytools::ancr)
## Why: node posteriors are used later for “localized” vs “diffuse” classification on the BC tree.
asr_res <- phytools::ancr(fit_mk, root.p = root_prior)
asr_mat <- asr_res$ace
if (!all(c("0", "1") %in% colnames(asr_mat)))
  colnames(asr_mat) <- c("0", "1")

asr_node_table <- tibble::tibble(
  node          = as.integer(rownames(asr_mat)),
  prob_nonwoody = asr_mat[, "0"],
  prob_woody    = asr_mat[, "1"],
  state_hat     = ifelse(asr_mat[, "1"] >= 0.5, 1, 0)
)

write_csv(asr_node_table, file.path(out_dir, "3.1_asr_node_table.csv"))
saveRDS(
  list(
    fit_mk_er = fit_mk_er,
    fit_mk_ard = fit_mk_ard,
    fit_mk = fit_mk,
    best_model = best_model,
    asr_node_table = asr_node_table
  ),
  file.path(out_dir, "3_objects.rds")
)

## 3.4 Extract transition rates from best model

extract_rates <- function(fit, model) {
  if (model == "ER") {
    q <- as.numeric(fit$rates)[1]
    return(
      tibble::tibble(
        model = "ER",
        q01 = q,
        q10 = q
      )
    )
  }
  
  ## ARD
  idx <- fit$index.matrix
  rates <- as.numeric(fit$rates)
  
  if (is.null(idx)) {
    stop("[ERROR] fit has no index.matrix; cannot extract ARD rates.")
  }
  
  ## Add default dimnames if missing
  if (is.null(rownames(idx))) rownames(idx) <- c("0", "1")
  if (is.null(colnames(idx))) colnames(idx) <- c("0", "1")
  
  q01_idx <- idx["0", "1"]
  q10_idx <- idx["1", "0"]
  
  q01 <- rates[q01_idx]
  q10 <- rates[q10_idx]
  
  tibble::tibble(
    model = "ARD",
    q01 = q01,
    q10 = q10
  )
}

rate_tbl <- extract_rates(fit_mk, best_model) %>%
  dplyr::mutate(
    rate_ratio_q01_q10 = q01 / q10,
    rate_ratio_q10_q01 = q10 / q01,
    interpretation = dplyr::case_when(
      abs(log(rate_ratio_q01_q10)) < 0.2 ~ "Rates approximately symmetric",
      rate_ratio_q01_q10 > 1 ~ "H→W transitions faster than W→H",
      rate_ratio_q01_q10 < 1 ~ "W→H transitions faster than H→W"
    )
  )

write_csv(
  rate_tbl,
  file.path(out_dir, "3.2_mk_rate_summary.csv")
)

message("[OK] Wrote Mk rate summary table.")


## 4. SIMMAP posterior -----------------------------------------------------------

message("[INFO] Step 4: SIMMAP posterior (events + timings)")

## 4.1 SIMMAP sampling parameters
## Why:
##   The SIMMAP posterior provides uncertainty-aware event counts and timing
##   summaries. We run in small batches for stability, then normalize the output
##   so downstream steps always see a clean list of simmap trees.

nsim <- 500
batch_size <- 10
set.seed(1)

get_Q <- function(fit) {
  ## Prefer an explicit Q matrix if present
  if (!is.null(fit$Q) && is.matrix(fit$Q)) {
    Q <- fit$Q
  } else if (!is.null(fit$index.matrix) && !is.null(fit$rates)) {
    idx <- fit$index.matrix
    rates <- as.numeric(fit$rates)
    
    if (is.null(rownames(idx))) rownames(idx) <- c("0", "1")
    if (is.null(colnames(idx))) colnames(idx) <- c("0", "1")
    
    Q <- matrix(
      0,
      nrow = nrow(idx),
      ncol = ncol(idx),
      dimnames = dimnames(idx)
    )
    
    for (i in seq_len(nrow(idx))) {
      for (j in seq_len(ncol(idx))) {
        kk <- idx[i, j]
        if (!is.na(kk) && kk > 0) {
          Q[i, j] <- rates[kk]
        }
      }
    }
    diag(Q) <- -rowSums(Q)
  } else {
    stop("[ERROR] Could not extract or reconstruct Q from fitMk.")
  }
  
  ## Final checks
  if (is.null(rownames(Q))) rownames(Q) <- c("0", "1")
  if (is.null(colnames(Q))) colnames(Q) <- c("0", "1")
  
  Q <- Q[c("0", "1"), c("0", "1"), drop = FALSE]
  
  if (any(!is.finite(Q))) {
    stop("[ERROR] Q contains non-finite values.")
  }
  if (any(Q[cbind(1:2, 1:2)] >= 0)) {
    stop("[ERROR] Q diagonal should be negative.")
  }
  if (any(Q[row(Q) != col(Q)] < 0)) {
    stop("[ERROR] Q off-diagonals should be non-negative.")
  }
  
  Q
}

Q <- get_Q(fit_mk)

pi_simmap <- "estimated"  # let simmap estimate pi

## Normalize make.simmap output so we ALWAYS get a plain list of simmap trees
as_simmap_list <- function(x) {
  if (inherits(x, "multiSimmap")) {
    x <- unclass(x)
  }
  
  ## If nsim = 1, some versions may return a single simmap object rather than a list
  if (inherits(x, "simmap")) {
    x <- list(x)
  }
  
  if (!is.list(x)) {
    stop("[ERROR] make.simmap did not return a list-like object.")
  }
  
  x
}

assert_simmap_tree <- function(tr, ref_tree = tree_coded) {
  if (!inherits(tr, "phylo")) {
    stop("[ERROR] A SIMMAP replicate is not a phylo object.")
  }
  if (is.null(tr$maps) || length(tr$maps) != nrow(tr$edge)) {
    stop("[ERROR] A SIMMAP replicate has missing or malformed $maps.")
  }
  if (!identical(sort(tr$tip.label), sort(ref_tree$tip.label))) {
    stop("[ERROR] A SIMMAP replicate has unexpected tip labels.")
  }
  TRUE
}

## 4.2 Run SIMMAP in batches with progress
handlers(global = TRUE)
handlers("txtprogressbar")

maps <- with_progress({
  p <- progressor(steps = nsim)
  out_maps <- vector("list", nsim)
  
  counter <- 1L
  while (counter <= nsim) {
    nsim_batch <- min(batch_size, nsim - counter + 1L)
    
    batch_maps <- phytools::make.simmap(
      tree  = tree_coded,
      x     = x_fac,
      model = best_model,
      Q     = Q,
      pi    = pi_simmap,
      nsim  = nsim_batch
    )
    
    batch_maps <- as_simmap_list(batch_maps)
    
    if (length(batch_maps) != nsim_batch) {
      stop(
        "[ERROR] make.simmap batch returned ",
        length(batch_maps),
        " trees, expected ",
        nsim_batch,
        "."
      )
    }
    
    ## Validate each replicate before storing
    for (i in seq_len(nsim_batch)) {
      assert_simmap_tree(batch_maps[[i]], ref_tree = tree_coded)
      out_maps[[counter + i - 1L]] <- batch_maps[[i]]
      p()
    }
    
    counter <- counter + nsim_batch
  }
  
  out_maps
})

if (length(maps) != nsim) {
  stop("[ERROR] SIMMAP run produced ", length(maps), " maps; expected ", nsim, ".")
}

message("[OK] SIMMAP finished: ", length(maps), " posterior maps.")

saveRDS(
  list(
    nsim = nsim,
    batch_size = batch_size,
    best_model = best_model,
    Q = Q,
    pi_simmap = pi_simmap,
    maps = maps
  ),
  file.path(out_dir, "4.1_simmap_maps.rds")
)

## 4.2b Extract SIMMAP events and per-edge posterior summaries
## Why:
##   Step 5 needs:
##     - change_tbl   : event-level transition table across all SIMMAP replicates
##     - edge_summary : per-edge posterior support for H→W / W→H / any shift
##
##   We compute ages as time before present on tree_coded.

message("[INFO] Step 4.2b: extracting SIMMAP event table and per-edge summaries")

norm_state <- function(s) {
  s <- as.character(s)
  dplyr::case_when(
    s %in% c("0", "H", "h", "herb", "herbaceous", "nonwoody") ~ "H",
    s %in% c("1", "W", "w", "woody")                          ~ "W",
    TRUE ~ NA_character_
  )
}

node_age_coded <- get_node_ages_bp(tree_coded)

empty_change_tbl <- tibble::tibble(
  sim        = integer(0),
  edge_id    = integer(0),
  parent     = integer(0),
  child      = integer(0),
  direction  = character(0),
  change_age = numeric(0)
)

extract_changes_one_sim <- function(sim_tree, sim_id, node_age_vec) {
  out <- vector("list", length(sim_tree$maps))
  k_out <- 0L
  
  for (e in seq_along(sim_tree$maps)) {
    segs <- sim_tree$maps[[e]]
    
    if (length(segs) <= 1) next
    
    parent <- sim_tree$edge[e, 1]
    child  <- sim_tree$edge[e, 2]
    
    parent_age <- node_age_vec[parent]
    if (!is.finite(parent_age)) next
    
    cur_age <- parent_age
    prev_state <- norm_state(names(segs)[1])
    
    if (is.na(prev_state)) next
    
    ## Walk from parent -> child along the mapped branch
    for (i in seq_along(segs)) {
      seg_len <- suppressWarnings(as.numeric(segs[i]))
      this_state <- norm_state(names(segs)[i])
      
      if (!is.finite(seg_len) || seg_len < 0 || is.na(this_state)) next
      
      if (i > 1) {
        ## Transition occurs at the start of this segment
        direction <- dplyr::case_when(
          prev_state == "H" & this_state == "W" ~ "H>W",
          prev_state == "W" & this_state == "H" ~ "W>H",
          TRUE ~ NA_character_
        )
        
        if (!is.na(direction)) {
          k_out <- k_out + 1L
          out[[k_out]] <- tibble::tibble(
            sim        = as.integer(sim_id),
            edge_id    = as.integer(e),
            parent     = as.integer(parent),
            child      = as.integer(child),
            direction  = direction,
            change_age = as.numeric(cur_age)
          )
        }
      }
      
      cur_age <- cur_age - seg_len
      prev_state <- this_state
    }
  }
  
  if (k_out == 0L) {
    return(empty_change_tbl)
  }
  
  dplyr::bind_rows(out[seq_len(k_out)])
}

change_tbl <- purrr::map2_dfr(
  maps,
  seq_along(maps),
  ~ extract_changes_one_sim(.x, .y, node_age_coded)
)

## Ensure stable column types even if empty
if (nrow(change_tbl) == 0) {
  change_tbl <- empty_change_tbl
} else {
  change_tbl <- change_tbl %>%
    dplyr::mutate(
      sim        = as.integer(sim),
      edge_id    = as.integer(edge_id),
      parent     = as.integer(parent),
      child      = as.integer(child),
      direction  = as.character(direction),
      change_age = as.numeric(change_age)
    )
}

write_csv(
  change_tbl,
  file.path(out_dir, "4.2_simmap_change_table_all_events.csv")
)

message("[INFO] Extracted ", nrow(change_tbl), " SIMMAP transition events.")

## Per-edge posterior summaries across replicates
all_edges_all_sims <- tidyr::expand_grid(
  sim = seq_along(maps),
  edge_id = seq_len(nrow(tree_coded$edge))
) %>%
  dplyr::mutate(
    parent = tree_coded$edge[edge_id, 1],
    child  = tree_coded$edge[edge_id, 2]
  )

if (nrow(change_tbl) == 0) {
  edge_event_counts <- tibble::tibble(
    sim = integer(0),
    edge_id = integer(0),
    `H>W` = integer(0),
    `W>H` = integer(0)
  )
} else {
  edge_event_counts <- change_tbl %>%
    dplyr::count(sim, edge_id, direction, name = "n_events") %>%
    tidyr::pivot_wider(
      names_from = direction,
      values_from = n_events,
      values_fill = 0
    )
  
  ## Guarantee both columns exist
  if (!"H>W" %in% names(edge_event_counts)) edge_event_counts$`H>W` <- 0L
  if (!"W>H" %in% names(edge_event_counts)) edge_event_counts$`W>H` <- 0L
  
  edge_event_counts <- edge_event_counts %>%
    dplyr::select(sim, edge_id, `H>W`, `W>H`)
}

edge_summary <- all_edges_all_sims %>%
  dplyr::left_join(edge_event_counts, by = c("sim", "edge_id")) %>%
  dplyr::mutate(
    `H>W` = dplyr::coalesce(`H>W`, 0L),
    `W>H` = dplyr::coalesce(`W>H`, 0L),
    any_shift = (`H>W` + `W>H`) > 0
  ) %>%
  dplyr::group_by(edge_id, parent, child) %>%
  dplyr::summarise(
    p_HW_any = mean(`H>W` > 0),
    p_WH_any = mean(`W>H` > 0),
    p_shift  = mean(any_shift),
    E_HW     = mean(`H>W`),
    E_WH     = mean(`W>H`),
    E_shift  = mean(`H>W` + `W>H`),
    .groups  = "drop"
  )

write_csv(
  edge_summary,
  file.path(out_dir, "4.2_simmap_edge_summary.csv")
)

message("[OK] Wrote SIMMAP event table and per-edge summary table.")

### 4.3 Combined SIMMAP summary (counts + timing)

message("[INFO] Step 4.3: combined SIMMAP summary")

## IMPORTANT:
##   We include ALL simulations here, including zero-event simulations.
##   Otherwise medians and proportions are biased upward and may break when a
##   direction is absent in some or all replicates.

sim_counts <- tidyr::expand_grid(sim = seq_along(maps)) %>%
  dplyr::left_join(
    change_tbl %>%
      dplyr::filter(direction %in% c("H>W", "W>H")) %>%
      dplyr::count(sim, direction, name = "n") %>%
      tidyr::pivot_wider(
        names_from = direction,
        values_from = n,
        values_fill = 0
      ),
    by = "sim"
  ) %>%
  dplyr::mutate(
    `H>W` = dplyr::coalesce(`H>W`, 0L),
    `W>H` = dplyr::coalesce(`W>H`, 0L),
    total = `H>W` + `W>H`
  )

count_summary <- sim_counts %>%
  dplyr::summarise(
    H_to_W_median = median(`H>W`, na.rm = TRUE),
    H_to_W_q025   = q025(`H>W`),
    H_to_W_q975   = q975(`H>W`),
    
    W_to_H_median = median(`W>H`, na.rm = TRUE),
    W_to_H_q025   = q025(`W>H`),
    W_to_H_q975   = q975(`W>H`),
    
    prop_H_to_W   = median(dplyr::if_else(total > 0, `H>W` / total, NA_real_), na.rm = TRUE),
    prop_W_to_H   = median(dplyr::if_else(total > 0, `W>H` / total, NA_real_), na.rm = TRUE)
  )

if (nrow(change_tbl) == 0) {
  timing_summary <- tibble::tibble(
    direction  = c("H>W", "W>H"),
    age_median = NA_real_,
    age_q025   = NA_real_,
    age_q975   = NA_real_
  )
} else {
  timing_summary <- change_tbl %>%
    dplyr::filter(direction %in% c("H>W", "W>H")) %>%
    dplyr::group_by(direction) %>%
    dplyr::summarise(
      age_median = median(change_age, na.rm = TRUE),
      age_q025   = q025(change_age),
      age_q975   = q975(change_age),
      .groups = "drop"
    )
  
  ## Guarantee both directions present
  timing_summary <- tibble::tibble(direction = c("H>W", "W>H")) %>%
    dplyr::left_join(timing_summary, by = "direction")
}

simmap_combined_tbl <- tibble::tibble(
  direction = c("H→W", "W→H"),
  n_events_median = c(count_summary$H_to_W_median, count_summary$W_to_H_median),
  n_events_q025   = c(count_summary$H_to_W_q025,   count_summary$W_to_H_q025),
  n_events_q975   = c(count_summary$H_to_W_q975,   count_summary$W_to_H_q975),
  proportion_median = c(count_summary$prop_H_to_W, count_summary$prop_W_to_H)
) %>%
  dplyr::left_join(
    timing_summary %>%
      dplyr::mutate(
        direction = dplyr::recode(direction, "H>W" = "H→W", "W>H" = "W→H")
      ),
    by = "direction"
  )

write_csv(
  simmap_combined_tbl,
  file.path(out_dir, "4.3_simmap_summary_counts_and_timing.csv")
)

saveRDS(
  list(
    Q = Q,
    maps = maps,
    change_tbl = change_tbl,
    edge_summary = edge_summary,
    sim_counts = sim_counts,
    simmap_combined_tbl = simmap_combined_tbl
  ),
  file.path(out_dir, "4_objects.rds")
)

message("[OK] Wrote combined SIMMAP summary table.")
message("[OK] Step 4 complete.")


## 5. Localization-aware synthesis + key figures ------------------------------
##
## Goal:
##   Summarize and visualize confident internal H↔W shifts on the
##   Brassicaceae + Cleomaceae (BC) tree.
##
## Definition of confident shift:
##   A BC-tree edge is considered a confident shift only if:
##     1) ASR implies a parent-child state flip on the BC tree, and
##     2) projected SIMMAP support is directionally concordant, and
##     3) direction-specific SIMMAP support >= p_event_min_confident.
##
## IMPORTANT:
##   Full-tree edge IDs and BC-tree edge IDs are NOT interchangeable after
##   pruning. Therefore, full-tree SIMMAP edge summaries and events are projected
##   onto the BC tree by descendant-tip mapping.

message("[INFO] Step 5: localization-aware synthesis + key figures")

if (!requireNamespace("phangorn", quietly = TRUE)) {
  stop("[ERROR] Package 'phangorn' is required for Step 5 descendant-tip mapping.")
}

### 5.0 Helper functions and BC-tree setup -------------------------------------

desc_tip_labels <- function(tree, node) {
  node <- as.integer(node)
  if (node <= length(tree$tip.label)) return(tree$tip.label[node])
  tips <- phangorn::Descendants(tree, node, type = "tips")[[1]]
  tree$tip.label[tips]
}

edge_index_tbl <- function(tree) {
  tibble::tibble(
    edge_id = seq_len(nrow(tree$edge)),
    parent  = tree$edge[, 1],
    child   = tree$edge[, 2]
  )
}

pick_majority_value <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

int_to_roman <- function(x) {
  x <- as.integer(x)
  if (length(x) == 0 || is.na(x) || x < 1) return(NA_character_)
  
  vals <- c(1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
  syms <- c("M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I")
  
  out <- character(1)
  n <- x
  
  for (i in seq_along(vals)) {
    while (n >= vals[i]) {
      out <- paste0(out, syms[i])
      n <- n - vals[i]
    }
  }
  
  out
}

norm_dir <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("H>W", "HW", "0>1", "H->W", "H→W") ~ "H>W",
    x %in% c("W>H", "WH", "1>0", "W->H", "W→H") ~ "W>H",
    TRUE ~ x
  )
}

## Keep Brassicaceae + Cleomaceae tips for context
keep_tbl <- specimen_details %>%
  dplyr::distinct(SAMPLE, SPECIES_NAME_PRINT, FAMILY, TRIBE, SUPERTRIBE) %>%
  dplyr::filter(FAMILY %in% c("Brassicaceae", "Cleomaceae"))

keep_tips <- intersect(tree_coded$tip.label, keep_tbl$SAMPLE)

if (length(keep_tips) < 10L) {
  stop("[ERROR] Too few tips after Brassicaceae + Cleomaceae filtering.")
}

tree_bc <- ape::keep.tip(tree_coded, keep_tips)
woody_fac_bc <- x_fac[tree_bc$tip.label]

keep_tbl_bc <- keep_tbl %>%
  dplyr::filter(SAMPLE %in% tree_bc$tip.label)

Ntip_bc  <- length(tree_bc$tip.label)
Nedge_bc <- nrow(tree_bc$edge)

message(
  "[INFO] BC tree tips=", Ntip_bc,
  " edges=", Nedge_bc,
  " internal=", tree_bc$Nnode
)

edge_index_full <- edge_index_tbl(tree_coded)
edge_index_bc   <- edge_index_tbl(tree_bc)

depth_bc <- ape::node.depth.edgelength(tree_bc)
tree_height_bc <- max(depth_bc)
node_age_bc <- tree_height_bc - depth_bc

edge_times_bc <- edge_index_bc %>%
  dplyr::mutate(
    age_start = node_age_bc[parent],
    age_end   = node_age_bc[child],
    age_mid   = (age_start + age_end) / 2
  )

## Brassicaceae stem and crown ages
brassi_tips <- keep_tbl_bc %>%
  dplyr::filter(FAMILY == "Brassicaceae") %>%
  dplyr::pull(SAMPLE) %>%
  intersect(tree_bc$tip.label)

if (length(brassi_tips) < 2) {
  stop("[ERROR] Too few Brassicaceae tips in tree_bc to compute crown/stem ages.")
}

brassi_crown_node <- ape::getMRCA(tree_bc, brassi_tips)
brassi_crown_age  <- as.numeric(node_age_bc[brassi_crown_node])

brassi_parent_edge <- which(tree_bc$edge[, 2] == brassi_crown_node)

brassi_stem_age <- if (length(brassi_parent_edge) == 1) {
  as.numeric(node_age_bc[tree_bc$edge[brassi_parent_edge, 1]])
} else {
  NA_real_
}

message("[INFO] Brassicaceae crown age: ", signif(brassi_crown_age, 6), " Ma")
message("[INFO] Brassicaceae stem age: ", signif(brassi_stem_age, 6), " Ma")

### 5.1 Project ASR node posteriors from full coded tree to BC nodes ------------

message("[INFO] Step 5.1: projecting ASR node posteriors to BC tree")

asr_full_map <- asr_node_table %>%
  dplyr::transmute(node_full = node, prob_woody_full = prob_woody) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    desc_bc = list(intersect(desc_tip_labels(tree_coded, node_full), tree_bc$tip.label))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(lengths(desc_bc) >= 2) %>%
  dplyr::mutate(
    node_bc = vapply(desc_bc, function(tips) {
      m <- ape::getMRCA(tree_bc, tips)
      if (is.null(m) || is.na(m)) NA_integer_ else as.integer(m)
    }, integer(1))
  ) %>%
  dplyr::filter(!is.na(node_bc)) %>%
  dplyr::group_by(node_bc) %>%
  dplyr::summarise(
    prob_woody = stats::median(prob_woody_full, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state_hat = dplyr::if_else(prob_woody >= 0.5, 1L, 0L)
  )

asr_node_table_bc <- asr_full_map %>%
  dplyr::transmute(node = node_bc, prob_woody, state_hat)

readr::write_csv(
  asr_node_table_bc,
  file.path(out_dir, "5.1_asr_node_table_bc_projected.csv")
)

### 5.2 Project SIMMAP per-edge posteriors from full tree to BC edges -----------

message("[INFO] Step 5.2: projecting SIMMAP edge summaries to BC tree")

edge_summary_full <- edge_summary %>%
  dplyr::select(
    edge_id,
    parent,
    child,
    p_HW_any,
    p_WH_any,
    p_shift,
    E_HW,
    E_WH,
    E_shift
  )

edge_summary_map <- edge_summary_full %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    desc_bc = list(intersect(desc_tip_labels(tree_coded, child), tree_bc$tip.label))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(lengths(desc_bc) >= 1) %>%
  dplyr::mutate(
    node_bc_child = vapply(desc_bc, function(tips) {
      if (length(tips) == 1) {
        return(match(tips, tree_bc$tip.label))
      }
      m <- ape::getMRCA(tree_bc, tips)
      if (is.null(m) || is.na(m)) NA_integer_ else as.integer(m)
    }, integer(1))
  ) %>%
  dplyr::filter(!is.na(node_bc_child)) %>%
  dplyr::left_join(
    edge_index_bc %>%
      dplyr::select(edge_id_bc = edge_id, child_bc = child),
    by = c("node_bc_child" = "child_bc")
  ) %>%
  dplyr::filter(!is.na(edge_id_bc)) %>%
  dplyr::group_by(edge_id_bc) %>%
  dplyr::summarise(
    p_HW_any = max(p_HW_any, na.rm = TRUE),
    p_WH_any = max(p_WH_any, na.rm = TRUE),
    p_shift  = max(p_shift,  na.rm = TRUE),
    E_HW     = sum(E_HW, na.rm = TRUE),
    E_WH     = sum(E_WH, na.rm = TRUE),
    E_shift  = sum(E_shift, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::rename(edge_id = edge_id_bc)

edge_summary_bc <- edge_index_bc %>%
  dplyr::select(edge_id) %>%
  dplyr::left_join(edge_summary_map, by = "edge_id") %>%
  dplyr::mutate(
    p_HW_any = dplyr::coalesce(p_HW_any, 0),
    p_WH_any = dplyr::coalesce(p_WH_any, 0),
    p_shift  = dplyr::coalesce(p_shift, 0),
    E_HW     = dplyr::coalesce(E_HW, 0),
    E_WH     = dplyr::coalesce(E_WH, 0),
    E_shift  = dplyr::coalesce(E_shift, 0)
  )

readr::write_csv(
  edge_summary_bc,
  file.path(out_dir, "5.2_edge_summary_bc_projected.csv")
)

### 5.3 ASR-derived edge flips on BC tree ---------------------------------------

message("[INFO] Step 5.3: identifying ASR-derived parent-child flips on BC tree")

pW_all_bc <- rep(NA_real_, Ntip_bc + tree_bc$Nnode)
pW_all_bc[seq_len(Ntip_bc)] <- as.numeric(woody_fac_bc == "1")
pW_all_bc[asr_node_table_bc$node] <- asr_node_table_bc$prob_woody

edge_node_asr_bc <- tibble::tibble(
  edge_id     = seq_len(Nedge_bc),
  parent      = tree_bc$edge[, 1],
  child       = tree_bc$edge[, 2],
  is_terminal = tree_bc$edge[, 2] <= Ntip_bc,
  pW_parent   = pW_all_bc[parent],
  pW_child    = pW_all_bc[child]
) %>%
  dplyr::mutate(
    indicate_ok  = !is.na(pW_parent) & !is.na(pW_child),
    state_parent = dplyr::if_else(pW_parent >= 0.5, "W", "H"),
    state_child  = dplyr::if_else(pW_child  >= 0.5, "W", "H"),
    node_flip    = indicate_ok & state_parent != state_child,
    node_dir     = dplyr::case_when(
      node_flip & state_parent == "H" & state_child == "W" ~ "H>W",
      node_flip & state_parent == "W" & state_child == "H" ~ "W>H",
      TRUE ~ NA_character_
    )
  )

readr::write_csv(
  edge_node_asr_bc,
  file.path(out_dir, "5.3_edge_node_asr_bc.csv")
)

### 5.4 Define confident shifts: ASR flip + SIMMAP support ----------------------

message("[INFO] Step 5.4: defining confident shifts")

edge_joined_bc <- edge_node_asr_bc %>%
  dplyr::left_join(edge_summary_bc, by = "edge_id") %>%
  dplyr::mutate(
    p_event_dir = dplyr::case_when(
      node_dir == "H>W" ~ p_HW_any,
      node_dir == "W>H" ~ p_WH_any,
      TRUE ~ NA_real_
    ),
    E_event_dir = dplyr::case_when(
      node_dir == "H>W" ~ E_HW,
      node_dir == "W>H" ~ E_WH,
      TRUE ~ NA_real_
    )
  )

## Restored original support threshold.
## A confident shift still requires an ASR flip; this threshold is only the
## additional directionally concordant SIMMAP support criterion.
p_event_min_confident <- 0.10

confident_shifts_all <- edge_joined_bc %>%
  dplyr::filter(
    node_flip,
    !is.na(node_dir),
    p_event_dir >= p_event_min_confident
  ) %>%
  dplyr::mutate(
    confidence_type = dplyr::if_else(
      is_terminal,
      "Confident terminal shift",
      "Confident internal shift"
    )
  )

diffuse_shift_signal_bc <- edge_joined_bc %>%
  dplyr::filter(
    !node_flip,
    p_shift >= p_event_min_confident
  )

readr::write_csv(
  confident_shifts_all,
  file.path(out_dir, "5.4_confident_shifts_all.csv")
)

confident_shift_summary <- confident_shifts_all %>%
  dplyr::mutate(
    subclass = dplyr::if_else(
      is_terminal,
      "Confident terminal shift",
      "Confident internal shift"
    )
  ) %>%
  dplyr::group_by(subclass, node_dir) %>%
  dplyr::summarise(
    n_edges = dplyr::n(),
    support_metric = "p_event_dir",
    support_median = stats::median(p_event_dir, na.rm = TRUE),
    support_q025   = q025(p_event_dir),
    support_q975   = q975(p_event_dir),
    .groups = "drop"
  )

readr::write_csv(
  confident_shift_summary,
  file.path(out_dir, "5.4_confident_shift_summary.csv")
)

message(
  "[OK] Confident shifts: ", nrow(confident_shifts_all),
  " (internal=", sum(!confident_shifts_all$is_terminal),
  "; terminal=", sum(confident_shifts_all$is_terminal), ")"
)

### 5.5 Circular BC tree with confident internal shifts -------------------------

message("[INFO] Step 5.5: circular BC tree with confident internal shifts")

tip_labels_bc <- tibble::tibble(label = tree_bc$tip.label) %>%
  dplyr::left_join(
    keep_tbl_bc %>%
      dplyr::select(SAMPLE, SPECIES_NAME_PRINT, SUPERTRIBE, TRIBE, FAMILY),
    by = c("label" = "SAMPLE")
  ) %>%
  dplyr::mutate(
    SPECIES_NAME_PRINT = dplyr::if_else(
      is.na(SPECIES_NAME_PRINT),
      label,
      SPECIES_NAME_PRINT
    )
  )

supertribe_cols2 <- clade_colors
supertribe_cols2["Aethionemeae"] <- "#ffed57"
supertribe_cols2["Cleomaceae"]   <- "grey70"

tip_annot_bc <- tibble::tibble(
  label = names(woody_fac_bc),
  tip_state = dplyr::if_else(as.character(woody_fac_bc) == "1", "W", "H")
)

node_annot_bc <- asr_node_table_bc %>%
  dplyr::mutate(
    node_state = dplyr::if_else(prob_woody >= 0.5, "W", "H")
  ) %>%
  dplyr::select(node, prob_woody, node_state)

edge_story_annot_bc <- tibble::tibble(
  edge_id = seq_len(Nedge_bc),
  dir_story = "none"
) %>%
  dplyr::left_join(
    confident_shifts_all %>%
      dplyr::transmute(edge_id, dir_story_conf = node_dir),
    by = "edge_id"
  ) %>%
  dplyr::mutate(
    dir_story = dplyr::coalesce(dir_story_conf, dir_story)
  ) %>%
  dplyr::select(edge_id, dir_story)

p_tree_bc <- ggtree::ggtree(
  tree_bc,
  layout = "circular",
  size   = 0.05,
  color  = "grey85"
) %<+% tip_labels_bc %<+% tip_annot_bc %<+% node_annot_bc

parent_pos_bc <- p_tree_bc$data %>%
  dplyr::select(parent_node = node, x_parent = x, y_parent = y)

p_tree_bc$data <- p_tree_bc$data %>%
  dplyr::left_join(
    edge_index_bc %>%
      dplyr::select(edge_id, child),
    by = c("node" = "child")
  ) %>%
  dplyr::left_join(edge_story_annot_bc, by = "edge_id") %>%
  dplyr::left_join(
    parent_pos_bc,
    by = c("parent" = "parent_node")
  ) %>%
  dplyr::mutate(
    x_mid_out = (x_parent + x) / 2,
    y_mid_out = y
  )

## Number confident internal shifts
confident_internal <- confident_shifts_all %>%
  dplyr::filter(!is_terminal)

switch_order_tbl <- edge_index_bc %>%
  dplyr::left_join(
    p_tree_bc$data %>%
      dplyr::select(node, y) %>%
      dplyr::distinct(),
    by = c("child" = "node")
  ) %>%
  dplyr::select(edge_id, y_child = y)

tip_to_species <- tip_labels_bc$SPECIES_NAME_PRINT
names(tip_to_species) <- tip_labels_bc$label

tip_to_family_vec <- keep_tbl_bc$FAMILY
names(tip_to_family_vec) <- keep_tbl_bc$SAMPLE

majority_family <- function(tips) {
  fam <- unname(tip_to_family_vec[tips])
  pick_majority_value(fam)
}

make_shift_name <- function(tips) {
  spp <- unname(tip_to_species[tips])
  spp <- spp[!is.na(spp)]
  if (length(spp) == 0) return(NA_character_)
  gen <- sub(" .*", "", spp)
  gen <- gen[!is.na(gen) & nzchar(gen)]
  if (length(gen) == 0) return(NA_character_)
  tab <- sort(table(gen), decreasing = TRUE)
  paste0(names(tab)[1], " clade")
}

shift_tbl_raw <- confident_internal %>%
  dplyr::left_join(switch_order_tbl, by = "edge_id") %>%
  dplyr::arrange(y_child) %>%
  dplyr::left_join(
    edge_index_bc %>%
      dplyr::select(edge_id, child_node = child),
    by = "edge_id"
  ) %>%
  dplyr::mutate(child_node = as.integer(child_node))

shift_tbl_raw$desc_tips <- lapply(
  shift_tbl_raw$child_node,
  function(n) desc_tip_labels(tree_bc, n)
)

shift_tbl_raw$major_family <- vapply(
  shift_tbl_raw$desc_tips,
  majority_family,
  character(1)
)

shift_tbl_raw$shift_name <- vapply(
  shift_tbl_raw$desc_tips,
  make_shift_name,
  character(1)
)

shift_tbl <- shift_tbl_raw %>%
  dplyr::arrange(y_child) %>%
  dplyr::mutate(
    brassicaceae_shift_id = dplyr::if_else(
      major_family == "Brassicaceae",
      cumsum(major_family == "Brassicaceae"),
      NA_integer_
    ),
    cleomaceae_shift_id = dplyr::if_else(
      major_family == "Cleomaceae",
      cumsum(major_family == "Cleomaceae"),
      NA_integer_
    ),
    shift_label = dplyr::case_when(
      major_family == "Brassicaceae" ~ as.character(brassicaceae_shift_id),
      major_family == "Cleomaceae"   ~ paste0("C", cleomaceae_shift_id),
      TRUE                           ~ NA_character_
    ),
    shift_id = dplyr::case_when(
      major_family == "Brassicaceae" ~ brassicaceae_shift_id,
      TRUE                           ~ NA_integer_
    )
  ) %>%
  dplyr::filter(!is.na(shift_label)) %>%
  dplyr::select(
    edge_id,
    shift_id,
    shift_label,
    shift_name,
    node_dir,
    child_node,
    major_family,
    p_event_dir,
    E_event_dir
  )

message(
  "[INFO] Labelled ",
  sum(shift_tbl$major_family == "Brassicaceae", na.rm = TRUE),
  " Brassicaceae internal confident shifts as 1, 2, ..."
)

message(
  "[INFO] Labelled ",
  sum(shift_tbl$major_family == "Cleomaceae", na.rm = TRUE),
  " Cleomaceae internal confident shifts as C1, C2, ..."
)

p_tree_bc$data <- p_tree_bc$data %>%
  dplyr::left_join(
    shift_tbl %>%
      dplyr::select(edge_id, shift_id, shift_label, shift_name, major_family),
    by = "edge_id"
  )

edge_pos_df <- p_tree_bc$data %>%
  dplyr::filter(!is.na(edge_id), !isTip) %>%
  dplyr::select(edge_id, x, y, x_parent, y_parent) %>%
  dplyr::distinct(edge_id, .keep_all = TRUE)

bubble_df <- shift_tbl %>%
  dplyr::left_join(edge_pos_df, by = "edge_id") %>%
  dplyr::mutate(
    x_mid_out = (x_parent + x) / 2,
    y_mid_out = y
  ) %>%
  dplyr::filter(
    is.finite(x_mid_out),
    is.finite(y_mid_out),
    !is.na(shift_label)
  )

message(
  "[CHECK 5.5] confident_internal=", nrow(confident_internal),
  "; shift_tbl=", nrow(shift_tbl),
  "; edge_pos_df=", nrow(edge_pos_df),
  "; bubble_df=", nrow(bubble_df)
)

fruit_df <- tip_labels_bc %>%
  dplyr::transmute(
    label,
    ring_rank = dplyr::case_when(
      FAMILY == "Cleomaceae"   ~ "Cleomaceae",
      TRIBE  == "Aethionemeae" ~ "Aethionemeae",
      TRUE                     ~ SUPERTRIBE
    )
  )

fruit_offset <- 1e-15
fruit_pwidth <- 0.03

p_tree_bc <- p_tree_bc +
  ggtree::geom_tree(size = 0.05, color = "grey85") +
  ggtree::geom_tree(
    aes(color = dplyr::na_if(dir_story, "none")),
    size = 0.15
  ) +
  ggplot2::scale_color_manual(
    name     = "Confident shifts",
    values   = c("H>W" = col_hw, "W>H" = col_wh),
    labels   = c("H>W" = "H→W", "W>H" = "W→H"),
    na.value = NA
  ) +
  ggplot2::geom_point(
    data = function(x) subset(x, !isTip),
    aes(x = x, y = y, fill = node_state),
    shape  = 21,
    size   = 0.3,
    colour = "transparent"
  ) +
  ggplot2::scale_fill_manual(
    name     = "Node state",
    values   = c("H" = col_wh, "W" = col_hw),
    na.value = "grey70"
  ) +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data    = fruit_df,
    geom    = geom_tile,
    mapping = aes(y = label, x = 1, fill = ring_rank),
    offset  = fruit_offset,
    pwidth  = fruit_pwidth
  ) +
  ggplot2::scale_fill_manual(
    name     = "Supertribe",
    values   = supertribe_cols2,
    na.value = "grey85"
  ) +
  ggnewscale::new_scale_color() +
  ggtree::geom_tippoint(aes(color = tip_state), size = 0.10) +
  ggplot2::scale_color_manual(
    name = "Tip state",
    values = c("H" = col_wh, "W" = col_hw)
  ) +
  ggtree::geom_tiplab(
    aes(label = SPECIES_NAME_PRINT),
    size   = 0.70,
    offset = 0.001
  ) +
  ggtree::theme_tree() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5, 120, 5, 5)
  ) +
  ggplot2::geom_point(
    data = bubble_df,
    aes(x = x_mid_out, y = y_mid_out),
    shape  = 21,
    size   = 1.6,
    stroke = 0.25,
    fill   = "white",
    colour = "black"
  ) +
  ggplot2::geom_text(
    data = bubble_df,
    aes(x = x_mid_out, y = y_mid_out, label = shift_label),
    size  = 0.8,
    vjust = 0.5
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "5.5_circular_tree_BrassiCleom_confident_internal.pdf"),
  plot     = p_tree_bc,
  width    = 30,
  height   = 30,
  dpi      = 300
)

ggplot2::ggsave(
  filename = file.path(out_dir, "5.5_circular_tree_BrassiCleom_confident_internal.png"),
  plot     = p_tree_bc,
  width    = 30,
  height   = 30,
  dpi      = 300
)

readr::write_csv(
  shift_tbl,
  file.path(out_dir, "5.5_confident_internal_shifts_lookup.csv")
)

message("[OK] Wrote Step 5.5 circular tree and lookup table.")

### 5.6 Project event table to BC edges and plot timing density -----------------

message("[INFO] Step 5.6: projecting SIMMAP event table to BC tree")

time_axis_lab <- "Time before present (Ma)"
shift_age_axis_lab <- "Estimated shift age—time before present (Ma)"

dir_col <- dplyr::case_when(
  "direction"  %in% names(change_tbl) ~ "direction",
  "node_dir"   %in% names(change_tbl) ~ "node_dir",
  "dir"        %in% names(change_tbl) ~ "dir",
  "transition" %in% names(change_tbl) ~ "transition",
  "change"     %in% names(change_tbl) ~ "change",
  TRUE ~ NA_character_
)

age_col <- dplyr::case_when(
  "change_age" %in% names(change_tbl) ~ "change_age",
  "age"        %in% names(change_tbl) ~ "age",
  TRUE ~ NA_character_
)

if (is.na(dir_col) || is.na(age_col)) {
  stop(
    "[ERROR] Could not find direction/age columns in change_tbl. Available columns: ",
    paste(names(change_tbl), collapse = ", ")
  )
}

change_tbl2 <- change_tbl %>%
  dplyr::transmute(
    sim          = as.integer(sim),
    edge_id_full = as.integer(edge_id),
    node_dir     = norm_dir(.data[[dir_col]]),
    age          = as.numeric(.data[[age_col]])
  ) %>%
  dplyr::filter(
    !is.na(sim),
    !is.na(edge_id_full),
    !is.na(node_dir),
    is.finite(age),
    node_dir %in% c("H>W", "W>H")
  )

edge_index_full2 <- tibble::tibble(
  edge_id_full = seq_len(nrow(tree_coded$edge)),
  child_full   = tree_coded$edge[, 2]
)

change_tbl2 <- change_tbl2 %>%
  dplyr::left_join(edge_index_full2, by = "edge_id_full")

unique_child_full <- sort(unique(change_tbl2$child_full))
child_map_tbl <- tibble::tibble(child_full = unique_child_full)

child_map_tbl$desc_bc <- lapply(child_map_tbl$child_full, function(child_node) {
  intersect(desc_tip_labels(tree_coded, child_node), tree_bc$tip.label)
})

child_map_tbl$child_bc <- vapply(child_map_tbl$desc_bc, function(tips) {
  if (length(tips) == 0) return(NA_integer_)
  if (length(tips) == 1) return(match(tips, tree_bc$tip.label))
  m <- ape::getMRCA(tree_bc, tips)
  if (is.null(m) || is.na(m)) NA_integer_ else as.integer(m)
}, integer(1))

edge_child_to_id_bc <- edge_index_bc %>%
  dplyr::transmute(edge_id_bc = edge_id, child_bc = child)

change_tbl_bc <- change_tbl2 %>%
  dplyr::left_join(
    child_map_tbl %>%
      dplyr::select(child_full, child_bc),
    by = "child_full"
  ) %>%
  dplyr::filter(!is.na(child_bc)) %>%
  dplyr::left_join(
    edge_child_to_id_bc,
    by = c("child_bc" = "child_bc")
  ) %>%
  dplyr::filter(!is.na(edge_id_bc)) %>%
  dplyr::transmute(
    sim,
    edge_id = as.integer(edge_id_bc),
    node_dir,
    age
  )

readr::write_csv(
  change_tbl_bc,
  file.path(out_dir, "5.6_change_tbl_projected_to_BC.csv")
)

confident_internal_edges <- confident_shifts_all %>%
  dplyr::filter(!is_terminal) %>%
  dplyr::transmute(edge_id, edge_dir = node_dir)

change_tbl_bc_confident <- change_tbl_bc %>%
  dplyr::inner_join(confident_internal_edges, by = "edge_id") %>%
  dplyr::filter(node_dir == edge_dir)

message(
  "[INFO] Events on BC confident internal edges, direction-consistent: ",
  nrow(change_tbl_bc_confident)
)

p_shift_density <- ggplot(
  change_tbl_bc_confident,
  aes(x = age, fill = node_dir, color = node_dir)
) +
  geom_density(alpha = 0.20, linewidth = 0.75) +
  scale_x_reverse() +
  scale_fill_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_color_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  labs(
    x = shift_age_axis_lab,
    y = "Density"
  ) +
  theme_classic()

ggplot2::ggsave(
  filename = file.path(out_dir, "5.6_confident_internal_shift_ages_density.png"),
  plot     = p_shift_density,
  width    = 6.5,
  height   = 4.0,
  dpi      = 350
)

ggplot2::ggsave(
  filename = file.path(out_dir, "5.6_confident_internal_shift_ages_density.pdf"),
  plot     = p_shift_density,
  width    = 6.5,
  height   = 4.0
)

message("[OK] Wrote Step 5.6 shift-age density plot.")

### 5.7 Per-shift timing plot with tribe-crown markers --------------------------

message("[INFO] Step 5.7: per-shift timing with tribe-crown markers")

if (!requireNamespace("scales", quietly = TRUE)) {
  stop("[ERROR] Package 'scales' is required for Step 5.7.")
}

tip_info_bc <- keep_tbl_bc %>%
  dplyr::select(SAMPLE, FAMILY, SUPERTRIBE, TRIBE) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    group = dplyr::case_when(
      FAMILY == "Cleomaceae"   ~ "Cleomaceae",
      TRIBE  == "Aethionemeae" ~ "Aethionemeae",
      TRUE                     ~ SUPERTRIBE
    )
  )

tip_lookup <- tip_info_bc %>%
  dplyr::select(SAMPLE, FAMILY, SUPERTRIBE, TRIBE, group)

shift_meta <- shift_tbl %>%
  dplyr::select(
    edge_id,
    shift_id,
    shift_label,
    shift_name,
    node_dir,
    child_node,
    major_family
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    child_node = as.integer(child_node),
    n_desc_tips = vapply(child_node, function(n) {
      length(desc_tip_labels(tree_bc, n))
    }, integer(1)),
    desc_tips = lapply(child_node, function(n) desc_tip_labels(tree_bc, n))
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    family = {
      df <- tip_lookup %>% dplyr::filter(SAMPLE %in% desc_tips)
      pick_majority_value(df$FAMILY)
    },
    tribe = {
      df <- tip_lookup %>% dplyr::filter(SAMPLE %in% desc_tips)
      pick_majority_value(df$TRIBE)
    },
    supertribe = {
      df <- tip_lookup %>% dplyr::filter(SAMPLE %in% desc_tips)
      pick_majority_value(df$SUPERTRIBE)
    },
    group = {
      df <- tip_lookup %>% dplyr::filter(SAMPLE %in% desc_tips)
      pick_majority_value(df$group)
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-desc_tips) %>%
  dplyr::mutate(
    tribe = dplyr::if_else(family == "Cleomaceae", NA_character_, tribe)
  )

tribe_tip_tbl <- keep_tbl_bc %>%
  dplyr::select(SAMPLE, FAMILY, TRIBE) %>%
  dplyr::distinct() %>%
  dplyr::filter(
    !is.na(TRIBE),
    TRIBE != "",
    SAMPLE %in% tree_bc$tip.label
  )

tribe_age_tbl <- tribe_tip_tbl %>%
  dplyr::group_by(TRIBE) %>%
  dplyr::summarise(
    tribe_family = dplyr::first(FAMILY),
    tribe_tips = list(SAMPLE),
    n_tribe_tips = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    tribe_mrca_node = vapply(tribe_tips, function(tips) {
      tips <- intersect(tips, tree_bc$tip.label)
      if (length(tips) == 0) return(NA_integer_)
      if (length(tips) == 1) return(match(tips, tree_bc$tip.label))
      m <- ape::getMRCA(tree_bc, tips)
      if (is.null(m) || is.na(m)) NA_integer_ else as.integer(m)
    }, integer(1)),
    tribe_crown_age = vapply(tribe_mrca_node, function(node) {
      if (is.na(node)) return(NA_real_)
      as.numeric(node_age_bc[node])
    }, numeric(1))
  ) %>%
  dplyr::transmute(
    tribe = TRIBE,
    n_tribe_tips,
    tribe_crown_age
  )

readr::write_csv(
  tribe_age_tbl,
  file.path(out_dir, "5.7_sampled_tribe_crown_ages_on_BC_tree.csv")
)

per_shift_age_raw <- change_tbl_bc_confident %>%
  dplyr::inner_join(
    shift_meta %>%
      dplyr::select(edge_id, shift_label, shift_name, node_dir),
    by = c("edge_id", "node_dir")
  ) %>%
  dplyr::group_by(edge_id, shift_label, shift_name, node_dir) %>%
  dplyr::summarise(
    n_events = dplyr::n(),
    age_med  = stats::median(age, na.rm = TRUE),
    age_lo   = q025(age),
    age_hi   = q975(age),
    .groups  = "drop"
  )

edge_age_tbl_bc <- edge_times_bc %>%
  dplyr::transmute(
    edge_id,
    edge_age_hi = age_start,
    edge_age_lo = age_end,
    edge_age_mid = age_mid
  )

per_shift_age <- shift_meta %>%
  dplyr::left_join(
    per_shift_age_raw,
    by = c("edge_id", "shift_label", "shift_name", "node_dir")
  ) %>%
  dplyr::left_join(tribe_age_tbl, by = "tribe") %>%
  dplyr::left_join(edge_age_tbl_bc, by = "edge_id") %>%
  dplyr::mutate(
    has_simmap_age = is.finite(age_med),
    tribe_crown_age = dplyr::if_else(
      node_dir == "H>W",
      tribe_crown_age,
      NA_real_
    ),
    lag_from_tribe_crown = dplyr::if_else(
      node_dir == "H>W" &
        is.finite(tribe_crown_age) &
        is.finite(age_med),
      tribe_crown_age - age_med,
      NA_real_
    ),
    tribe_shift_relationship = dplyr::case_when(
      node_dir != "H>W" ~ NA_character_,
      !is.finite(tribe_crown_age) ~ NA_character_,
      !is.finite(age_med) ~ NA_character_,
      tribe_crown_age >= age_lo & tribe_crown_age <= age_hi ~ "coincident with tribe crown",
      tribe_crown_age > age_hi ~ "subsequent to tribe crown",
      tribe_crown_age < age_lo ~ "check: shift appears older than sampled tribe crown",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(
    dplyr::case_when(
      major_family == "Cleomaceae" ~ 0L,
      major_family == "Brassicaceae" ~ 1L,
      TRUE ~ 2L
    ),
    suppressWarnings(as.numeric(shift_id)),
    shift_label
  ) %>%
  dplyr::group_by(shift_name) %>%
  dplyr::mutate(
    n_same_name = dplyr::n(),
    same_name_index = dplyr::row_number(),
    shift_name_display = dplyr::if_else(
      n_same_name == 1,
      shift_name,
      paste0(shift_name, " ", vapply(same_name_index, int_to_roman, character(1)))
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    shift_name_display_axis = stringr::str_remove_all(shift_name_display, "\\bclade\\s+|\\s+clade\\b"),
    shift_name_display_axis = stringr::str_squish(shift_name_display_axis),
    shift_label2 = paste0(shift_label, "  ", shift_name_display_axis),
    shift_label2 = factor(shift_label2, levels = rev(unique(shift_label2)))
  )

readr::write_csv(
  per_shift_age,
  file.path(out_dir, "5.7_confident_internal_shift_age_summary.csv")
)

## 5.7b Summary of lag times between tribe crown ages and H→W shifts

lag_tbl <- per_shift_age %>%
  dplyr::filter(
    node_dir == "H>W",
    is.finite(lag_from_tribe_crown)
  )

lag_summary_one <- function(df, subset_label) {
  
  zero_tol <- 0.25
  
  n_zero <- sum(abs(df$lag_from_tribe_crown) <= zero_tol, na.rm = TRUE)
  prop_zero <- n_zero / nrow(df)
  
  tibble::tibble(
    subset = subset_label,
    n_shifts = nrow(df),
    lag_min = min(df$lag_from_tribe_crown, na.rm = TRUE),
    lag_q025 = q025(df$lag_from_tribe_crown),
    lag_median = median(df$lag_from_tribe_crown, na.rm = TRUE),
    lag_mean = mean(df$lag_from_tribe_crown, na.rm = TRUE),
    lag_q975 = q975(df$lag_from_tribe_crown),
    lag_max = max(df$lag_from_tribe_crown, na.rm = TRUE),
    n_zero_lag = n_zero,
    prop_zero_lag = prop_zero,
    interpretation = paste0(
      "Median lag = ",
      sprintf("%.2f", median(df$lag_from_tribe_crown, na.rm = TRUE)),
      " Ma (range ",
      sprintf("%.2f", min(df$lag_from_tribe_crown, na.rm = TRUE)),
      "–",
      sprintf("%.2f", max(df$lag_from_tribe_crown, na.rm = TRUE)),
      " Ma); ",
      sprintf("%.1f", 100 * prop_zero),
      "% effectively zero lag (|lag| <= ",
      zero_tol,
      " Ma)."
    )
  )
}

lag_summary_tbl <- dplyr::bind_rows(
  
  lag_summary_one(
    lag_tbl,
    "All H→W shifts"
  ),
  
  lag_summary_one(
    lag_tbl %>%
      dplyr::filter(
        !tribe %in% c("Aethionemeae")
      ),
    "Excluding tribe Aethionemeae"
  )
)

readr::write_csv(
  lag_summary_tbl,
  file.path(out_dir, "5.7_lag_time_summary.csv")
)

message("[OK] Wrote lag-time summary table.")

message(
  "[CHECK] 5.7 plotted shifts total=", nrow(per_shift_age),
  "; with ages=", sum(is.finite(per_shift_age$age_med)),
  "; without ages=", sum(!is.finite(per_shift_age$age_med))
)

plot_df <- per_shift_age %>%
  dplyr::select(-dplyr::any_of(c("y_num", "y_lag", "y_num.x", "y_num.y"))) %>%
  dplyr::mutate(
    impact      = sqrt(pmax(n_desc_tips, 1)),
    impact_lw   = scales::rescale(impact, to = c(0.4, 2.0)),
    pt_size     = scales::rescale(impact, to = c(2.0, 4.0)),
    origin_size = scales::rescale(impact, to = c(2.0, 3.6))
  )

y_lu <- plot_df %>%
  dplyr::distinct(shift_label2) %>%
  dplyr::mutate(y_num = dplyr::row_number())

plot_df <- plot_df %>%
  dplyr::left_join(y_lu, by = "shift_label2") %>%
  dplyr::mutate(
    y_lag = y_num
  )

shade_df <- plot_df %>%
  dplyr::mutate(
    ymin = y_num - 0.5,
    ymax = y_num + 0.5,
    fill_group = dplyr::if_else(group == "Cleomaceae", NA_character_, group)
  ) %>%
  dplyr::select(ymin, ymax, fill_group) %>%
  dplyr::distinct()

missing_cols <- setdiff(unique(na.omit(shade_df$fill_group)), names(supertribe_cols2))

if (length(missing_cols) > 0) {
  warning(
    "[WARN] Some background groups have no color in supertribe_cols2; assigning grey85: ",
    paste(missing_cols, collapse = ", ")
  )
  for (nm in missing_cols) supertribe_cols2[nm] <- "grey85"
}

tribe_marker_df <- plot_df %>%
  dplyr::filter(
    node_dir == "H>W",
    is.finite(tribe_crown_age),
    is.finite(age_med)
  )

brassi_annot_y <- max(y_lu$y_num, na.rm = TRUE) + 0.35

p_per_shift <- ggplot() +
  geom_rect(
    data = shade_df,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax,
      fill = fill_group
    ),
    alpha = 0.10,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Supertribe",
    values = supertribe_cols2,
    na.value = NA,
    guide = guide_legend(override.aes = list(alpha = 0.3))
  ) +
  geom_vline(
    xintercept = brassi_crown_age,
    linewidth = 0.45,
    linetype = 3,
    color = "grey35"
  ) +
  annotate(
    "text",
    x = brassi_crown_age,
    y = brassi_annot_y,
    label = "Brassicaceae crown age",
    angle = 90,
    vjust = -0.35,
    hjust = 1,
    size = 2.8,
    color = "grey20"
  ) +
  geom_segment(
    data = tribe_marker_df,
    aes(
      x = tribe_crown_age,
      xend = age_med,
      y = y_lag,
      yend = y_lag,
      linetype = "Sampled tribe crown age + lag"
    ),
    linewidth = 0.3,
    color = "grey",
    alpha = 0.35,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  geom_errorbarh(
    data = plot_df,
    aes(
      y = y_num,
      xmin = age_lo,
      xmax = age_hi,
      color = node_dir,
      linewidth = impact_lw
    ),
    height = 0
  ) +
  geom_point(
    data = tribe_marker_df,
    aes(
      x = tribe_crown_age,
      y = y_num,
      shape = "Sampled tribe crown age + lag"
    ),
    size = 2.7,
    color = "black",
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  geom_point(
    data = plot_df %>% dplyr::filter(!has_simmap_age),
    aes(
      x = edge_age_mid,
      y = y_num
    ),
    shape = 4,
    size = 1.6,
    stroke = 0.35,
    color = "grey35",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = plot_df,
    aes(
      x = age_med,
      y = y_num,
      color = node_dir,
      size = pt_size
    ),
    stroke = 0.2
  ) +
  scale_color_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_linewidth_identity(guide = "none") +
  scale_size_identity(guide = "none") +
  scale_shape_manual(
    name = NULL,
    values = c("Sampled tribe crown age + lag" = 18)
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Sampled tribe crown age + lag" = 1)
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(linewidth = 1.0, shape = NA, size = NA, alpha = 1)
    ),
    fill = guide_legend(
      order = 2,
      override.aes = list(shape = 22, size = 5, colour = NA, alpha = 0.3)
    ),
    shape = guide_legend(
      order = 3,
      override.aes = list(color = "black", fill = "black", size = 3)
    ),
    linetype = guide_legend(
      order = 3,
      override.aes = list(color = "black", linewidth = 0.7, alpha = 0.35)
    )
  ) +
  scale_x_reverse() +
  scale_y_continuous(
    breaks = y_lu$y_num,
    labels = as.character(y_lu$shift_label2),
    expand = ggplot2::expansion(add = c(0.15, 0.55))
  ) +
  labs(
    x = "Estimated shift age—time before present (Ma)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = c(0.03, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "5.7_confident_internal_shift_ages_by_shift.png"),
  plot     = p_per_shift,
  width    = 10,
  height   = max(6, 0.18 * nrow(plot_df) + 2),
  dpi      = 300
)

ggplot2::ggsave(
  filename = file.path(out_dir, "5.7_confident_internal_shift_ages_by_shift.pdf"),
  plot     = p_per_shift,
  width    = 10,
  height   = max(6, 0.18 * nrow(plot_df) + 2)
)

message("[OK] Wrote Step 5.7 per-shift timing plot.")

message("[OK] Step 5 complete.")


message("[CHECK] Full-tree SIMMAP events: ", nrow(change_tbl))
message("[CHECK] BC-projected SIMMAP events: ", nrow(change_tbl_bc))

message("[CHECK] Full-tree event directions:")
print(table(change_tbl$direction, useNA = "ifany"))

message("[CHECK] BC-projected event directions:")
print(table(change_tbl_bc$node_dir, useNA = "ifany"))

message("[CHECK] BC tree total branch length: ", sum(tree_bc$edge.length, na.rm = TRUE))
message("[CHECK] tree_coded total branch length: ", sum(tree_coded$edge.length, na.rm = TRUE))
message("[CHECK] BC tree height: ", max(ape::node.depth.edgelength(tree_bc)))
message("[CHECK] coded tree height: ", max(ape::node.depth.edgelength(tree_coded)))



## 6. Process-level tempo, propensity, and persistence (SIMMAP + Mk) ----------
##
## Scope:
##   Step 6 is process-level. It uses all BC-projected SIMMAP events, not only
##   the confident internal shifts shown in Step 5.
##
## Outputs:
##   6.1_shift_rate_through_time_all_events.csv
##   6.2_shift_propensity_through_time_all_events.csv
##   6.3_spells_attained_with_censoring_all.rds
##   6.3_persistence_spells_with_censoring_summary.csv
##   6.4_mk_implied_sojourn_times.csv
##   6.1_6.2_tempo_propensity_2panel.(png/pdf)
##   6.4_mk_sojourn_times.(png/pdf)

message("[INFO] Step 6: process-level tempo + propensity + persistence (SIMMAP + Mk)")

required_change_cols <- c("sim", "edge_id", "node_dir", "age")

if (!exists("change_tbl_bc") || !all(required_change_cols %in% names(change_tbl_bc))) {
  stop(
    "[ERROR] Step 6 expects projected BC event table from corrected Step 5.6 with columns: ",
    paste(required_change_cols, collapse = ", ")
  )
}

RUN_EXPLORATORY_KM <- FALSE

bin_width <- 2
ymax_show <- 0.15

### 6.0 Helpers ----------------------------------------------------------------

norm_state <- function(s) {
  s <- as.character(s)
  if (s %in% c("0", "H", "h", "herb", "herbaceous", "nonwoody")) return("H")
  if (s %in% c("1", "W", "w", "woody"))                         return("W")
  NA_character_
}

norm_dir <- function(x) {
  dplyr::case_when(
    x %in% c("H>W", "HW", "0>1", "H->W", "H→W") ~ "H>W",
    x %in% c("W>H", "WH", "1>0", "W->H", "W→H") ~ "W>H",
    TRUE ~ as.character(x)
  )
}

add_brassi_origin_band <- function(p, stem_age, crown_age) {
  p +
    ggplot2::annotate(
      "rect",
      xmin = min(stem_age, crown_age),
      xmax = max(stem_age, crown_age),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.08,
      fill = "grey40"
    ) +
    ggplot2::geom_vline(
      xintercept = stem_age,
      linewidth = 0.4,
      linetype = 2,
      colour = "grey35"
    ) +
    ggplot2::geom_vline(
      xintercept = crown_age,
      linewidth = 0.4,
      linetype = 3,
      colour = "grey35"
    )
}

add_crown_label <- function(p, y = Inf) {
  p +
    ggplot2::annotate(
      "text",
      x = brassi_crown_age - 1.0,
      y = y,
      label = "Brassicaceae crown age",
      angle = 90,
      vjust = 0,
      hjust = 1,
      size = 2.6,
      color = "grey20"
    )
}

make_children_map <- function(tree) {
  split(tree$edge[, 2], tree$edge[, 1])
}

sim_ids <- seq_along(maps)

if (length(sim_ids) == 0) {
  stop("[ERROR] Step 6: no SIMMAP replicates found in 'maps'.")
}

depth_from_root_bc <- ape::node.depth.edgelength(tree_bc)
tree_height_bc     <- max(depth_from_root_bc)
node_age_bc2       <- tree_height_bc - depth_from_root_bc

edge_times_bc <- tibble::tibble(
  edge_id   = seq_len(nrow(tree_bc$edge)),
  parent    = tree_bc$edge[, 1],
  child     = tree_bc$edge[, 2],
  age_start = node_age_bc2[parent],
  age_end   = node_age_bc2[child]
)

stopifnot(all(edge_times_bc$age_start + 1e-10 >= edge_times_bc$age_end))

max_age <- max(edge_times_bc$age_start, na.rm = TRUE)
breaks  <- seq(0, ceiling(max_age / bin_width) * bin_width, by = bin_width)

bin_tbl <- tibble::tibble(
  bin    = seq_len(length(breaks) - 1),
  bin_lo = breaks[-length(breaks)],
  bin_hi = breaks[-1]
) %>%
  dplyr::mutate(
    bin_mid = (bin_lo + bin_hi) / 2
  )

events_binned <- function(events_df) {
  events_df %>%
    dplyr::mutate(
      node_dir = norm_dir(node_dir),
      bin = findInterval(age, vec = breaks, rightmost.closed = TRUE),
      bin = dplyr::if_else(
        bin < 1L | bin > nrow(bin_tbl),
        NA_integer_,
        as.integer(bin)
      )
    ) %>%
    dplyr::filter(
      !is.na(bin),
      node_dir %in% c("H>W", "W>H")
    )
}

exposure_tbl <- bin_tbl %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    exposure = sum(
      pmax(
        0,
        pmin(edge_times_bc$age_start, bin_hi) -
          pmax(edge_times_bc$age_end, bin_lo)
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(bin, bin_lo, bin_hi, bin_mid, exposure)

### 6.1 Tempo ------------------------------------------------------------------

message("[INFO] Step 6.1: tempo through time")

rate_summary <- {
  ev <- events_binned(change_tbl_bc)
  
  full_grid <- tidyr::expand_grid(
    sim = sim_ids,
    bin = bin_tbl$bin,
    node_dir = c("H>W", "W>H")
  )
  
  ev_counts <- ev %>%
    dplyr::count(sim, bin, node_dir, name = "n_events")
  
  per_sim <- full_grid %>%
    dplyr::left_join(ev_counts, by = c("sim", "bin", "node_dir")) %>%
    dplyr::mutate(
      n_events = dplyr::coalesce(n_events, 0L)
    ) %>%
    dplyr::left_join(exposure_tbl, by = "bin") %>%
    dplyr::mutate(
      rate = dplyr::if_else(exposure > 0, n_events / exposure, NA_real_)
    )
  
  per_sim %>%
    dplyr::group_by(bin, bin_mid, bin_lo, bin_hi, node_dir) %>%
    dplyr::summarise(
      rate_med = stats::median(rate, na.rm = TRUE),
      rate_lo  = q025(rate),
      rate_hi  = q975(rate),
      .groups = "drop"
    )
}

readr::write_csv(
  rate_summary,
  file.path(out_dir, "6.1_shift_rate_through_time_all_events.csv")
)

p_rate <- ggplot(
  rate_summary,
  aes(x = bin_mid, y = rate_med, color = node_dir, fill = node_dir)
) +
  geom_ribbon(
    aes(ymin = rate_lo, ymax = rate_hi),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_fill_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(0, ymax_show)) +
  labs(
    x = NULL,
    y = "Shift rate\n(events · Myr⁻¹ lineage time)",
    title = "A"
  ) +
  theme_classic()

p_rate <- add_brassi_origin_band(p_rate, brassi_stem_age, brassi_crown_age) +
  theme(
    legend.position = c(0.86, 0.83),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8)
  )

### 6.2 Propensity --------------------------------------------------------------

message("[INFO] Step 6.2: state-specific shift propensity through time")

exposure_from_sim_bc <- function(sim_bc) {
  node_age_sim <- get_node_ages_bp(sim_bc)
  
  out <- tidyr::expand_grid(
    bin = bin_tbl$bin,
    state = c("H", "W")
  ) %>%
    dplyr::mutate(exposure = 0)
  
  for (e in seq_len(nrow(sim_bc$edge))) {
    parent <- sim_bc$edge[e, 1]
    child  <- sim_bc$edge[e, 2]
    
    age_start <- node_age_sim[parent]
    age_end   <- node_age_sim[child]
    
    if (!is.finite(age_start) || !is.finite(age_end)) next
    
    segs <- sim_bc$maps[[e]]
    if (length(segs) == 0) next
    
    cur <- age_start
    
    for (k in seq_along(segs)) {
      st <- norm_state(names(segs)[k])
      L  <- as.numeric(segs[k])
      
      if (!st %in% c("H", "W") || !is.finite(L) || L <= 0) {
        cur <- cur - L
        next
      }
      
      seg_start <- cur
      seg_end   <- cur - L
      cur       <- seg_end
      
      overlaps <- pmax(
        0,
        pmin(seg_start, bin_tbl$bin_hi) -
          pmax(seg_end, bin_tbl$bin_lo)
      )
      
      if (any(overlaps > 0)) {
        idx <- which(overlaps > 0)
        out$exposure[out$state == st & out$bin %in% idx] <-
          out$exposure[out$state == st & out$bin %in% idx] + overlaps[idx]
      }
    }
  }
  
  out
}

ev_counts <- events_binned(change_tbl_bc) %>%
  dplyr::count(sim, bin, node_dir, name = "n_events")

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

expo_list <- progressr::with_progress({
  p <- progressr::progressor(steps = length(sim_ids))
  
  lapply(sim_ids, function(s) {
    p()
    
    sim_tree <- maps[[s]]
    drop <- setdiff(sim_tree$tip.label, tree_bc$tip.label)
    sim_bc <- if (length(drop) > 0) ape::drop.tip(sim_tree, drop) else sim_tree
    
    ex <- exposure_from_sim_bc(sim_bc)
    ex$sim <- s
    ex
  })
})

expo_tbl <- dplyr::bind_rows(expo_list)

prop_summary <- {
  full_grid <- tidyr::expand_grid(
    sim = sim_ids,
    bin = bin_tbl$bin,
    node_dir = c("H>W", "W>H")
  ) %>%
    dplyr::mutate(
      origin_state = dplyr::if_else(node_dir == "H>W", "H", "W")
    ) %>%
    dplyr::left_join(
      expo_tbl,
      by = c("sim", "bin", "origin_state" = "state")
    ) %>%
    dplyr::left_join(
      ev_counts,
      by = c("sim", "bin", "node_dir")
    ) %>%
    dplyr::mutate(
      n_events = dplyr::coalesce(n_events, 0L),
      prop = dplyr::if_else(exposure > 0, n_events / exposure, NA_real_)
    )
  
  full_grid %>%
    dplyr::left_join(bin_tbl, by = "bin") %>%
    dplyr::group_by(bin, bin_mid, bin_lo, bin_hi, node_dir) %>%
    dplyr::summarise(
      prop_med = stats::median(prop, na.rm = TRUE),
      prop_lo  = q025(prop),
      prop_hi  = q975(prop),
      .groups = "drop"
    )
}

readr::write_csv(
  prop_summary,
  file.path(out_dir, "6.2_shift_propensity_through_time_all_events.csv")
)

p_prop <- ggplot(
  prop_summary,
  aes(x = bin_mid, y = prop_med, color = node_dir, fill = node_dir)
) +
  geom_ribbon(
    aes(ymin = prop_lo, ymax = prop_hi),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_fill_manual(
    values = c("H>W" = col_hw, "W>H" = col_wh),
    name = "Direction",
    labels = c("H>W" = "H→W", "W>H" = "W→H")
  ) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(0, ymax_show)) +
  labs(
    x = time_axis_lab,
    y = "Shift propensity\n(events · Myr⁻¹ spent in origin state)",
    title = "B"
  ) +
  theme_classic() +
  theme(legend.position = "none")

p_prop <- add_brassi_origin_band(p_prop, brassi_stem_age, brassi_crown_age)

### 6.2b Save tempo + propensity 2-panel figure --------------------------------

p_rate2 <- p_rate + ggplot2::theme(plot.margin = ggplot2::margin(6, 18, 6, 6))
p_prop2 <- p_prop + ggplot2::theme(plot.margin = ggplot2::margin(6, 18, 6, 6))

out_png <- file.path(out_dir, "6.1_6.2_tempo_propensity_2panel.png")
out_pdf <- file.path(out_dir, "6.1_6.2_tempo_propensity_2panel.pdf")

combine_ok <- FALSE

if (requireNamespace("patchwork", quietly = TRUE)) {
  combine_ok <- TRUE
  p_combo_6 <- (p_rate2 / p_prop2) + patchwork::plot_layout(heights = c(1, 1))
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  combine_ok <- TRUE
  p_combo_6 <- cowplot::plot_grid(
    p_rate2,
    p_prop2,
    ncol = 1,
    rel_heights = c(1, 1),
    align = "v",
    axis = "lr"
  )
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  combine_ok <- TRUE
  p_combo_6 <- gridExtra::arrangeGrob(
    p_rate2,
    p_prop2,
    ncol = 1,
    heights = c(1, 1)
  )
}

if (!combine_ok) {
  warning("[WARN] Could not combine tempo/propensity panels; writing separate files.")
  
  ggsave(
    file.path(out_dir, "6.1_shift_rate_through_time_all_events.png"),
    plot = p_rate2,
    width = 4.2,
    height = 5.2,
    dpi = 400
  )
  
  ggsave(
    file.path(out_dir, "6.2_shift_propensity_through_time_all_events.png"),
    plot = p_prop2,
    width = 4.2,
    height = 5.2,
    dpi = 400
  )
  
} else {
  ggsave(out_png, plot = p_combo_6, width = 4.2, height = 7, dpi = 400)
  ggsave(out_pdf, plot = p_combo_6, width = 4.2, height = 7)
  message("[OK] Wrote tempo/propensity 2-panel figure.")
}

### 6.3 Descriptive SIMMAP spells with right-censoring --------------------------
# 
# message("[INFO] Step 6.3: descriptive SIMMAP spell summaries with right-censoring")
# 
# extract_spells_attained_with_censoring <- function(simmap_tree, tree_ref) {
#   Ntip <- length(tree_ref$tip.label)
#   root <- Ntip + 1L
#   node_age <- get_node_ages_bp(tree_ref)
#   children <- make_children_map(tree_ref)
#   
#   edge_child <- tree_ref$edge[, 2]
#   edge_index_by_child <- setNames(seq_len(nrow(tree_ref$edge)), edge_child)
#   
#   root_kids <- children[[as.character(root)]]
#   if (is.null(root_kids) || length(root_kids) == 0) {
#     stop("[ERROR] Could not find root children.")
#   }
#   
#   e0 <- edge_index_by_child[as.character(root_kids[1])]
#   seg0 <- simmap_tree$maps[[e0]]
#   root_state <- norm_state(names(seg0)[1])
#   
#   if (is.na(root_state)) {
#     stop("[ERROR] Could not infer root state from simmap.")
#   }
#   
#   stack <- list(list(
#     node = root,
#     state = root_state,
#     t = 0,
#     started = FALSE,
#     attain_age = NA_real_
#   ))
#   
#   out <- list()
#   k_out <- 0L
#   
#   while (length(stack) > 0) {
#     cur_item <- stack[[length(stack)]]
#     stack <- stack[-length(stack)]
#     
#     node <- cur_item$node
#     st   <- cur_item$state
#     tcur <- cur_item$t
#     started <- cur_item$started
#     attain_age <- cur_item$attain_age
#     
#     kids <- children[[as.character(node)]]
#     
#     if (is.null(kids) || length(kids) == 0) {
#       if (started && is.finite(tcur) && tcur > 0 && is.finite(attain_age)) {
#         k_out <- k_out + 1L
#         out[[k_out]] <- tibble::tibble(
#           state = st,
#           attain_age = attain_age,
#           duration = tcur,
#           event = 0L
#         )
#       }
#       next
#     }
#     
#     for (ch in kids) {
#       e_idx <- edge_index_by_child[as.character(ch)]
#       segs <- simmap_tree$maps[[e_idx]]
#       if (length(segs) == 0) next
#       
#       parent_age <- node_age[node]
#       if (!is.finite(parent_age)) next
#       
#       cur_age <- parent_age
#       st_edge <- norm_state(names(segs)[1])
#       if (is.na(st_edge)) next
#       
#       st_now <- st_edge
#       t_now  <- tcur
#       started_now <- started
#       attain_age_now <- attain_age
#       
#       if (!is.na(st) && !is.na(st_now) && st_now != st) {
#         if (started_now && is.finite(t_now) && t_now > 0 && is.finite(attain_age_now)) {
#           k_out <- k_out + 1L
#           out[[k_out]] <- tibble::tibble(
#             state = st,
#             attain_age = attain_age_now,
#             duration = t_now,
#             event = 1L
#           )
#         }
#         
#         started_now <- TRUE
#         t_now <- 0
#         attain_age_now <- parent_age
#       }
#       
#       for (i in seq_along(segs)) {
#         seg_state <- norm_state(names(segs)[i])
#         seg_len   <- as.numeric(segs[i])
#         
#         if (is.na(seg_state) || !is.finite(seg_len) || seg_len < 0) next
#         
#         if (i == 1) {
#           t_now <- t_now + seg_len
#           cur_age <- cur_age - seg_len
#         } else {
#           trans_age <- cur_age
#           
#           if (started_now && is.finite(t_now) && t_now > 0 && is.finite(attain_age_now)) {
#             k_out <- k_out + 1L
#             out[[k_out]] <- tibble::tibble(
#               state = st_now,
#               attain_age = attain_age_now,
#               duration = t_now,
#               event = 1L
#             )
#           }
#           
#           st_now <- seg_state
#           started_now <- TRUE
#           attain_age_now <- trans_age
#           t_now <- seg_len
#           cur_age <- cur_age - seg_len
#         }
#       }
#       
#       stack[[length(stack) + 1]] <- list(
#         node = ch,
#         state = st_now,
#         t = t_now,
#         started = started_now,
#         attain_age = attain_age_now
#       )
#     }
#   }
#   
#   if (k_out == 0L) {
#     return(tibble::tibble(
#       state = character(0),
#       attain_age = numeric(0),
#       duration = numeric(0),
#       event = integer(0)
#     ))
#   }
#   
#   dplyr::bind_rows(out)
# }
# 
# spell_tbl_all <- progressr::with_progress({
#   p <- progressr::progressor(steps = length(sim_ids))
#   
#   dplyr::bind_rows(lapply(sim_ids, function(s) {
#     p()
#     
#     sim_tree <- maps[[s]]
#     drop <- setdiff(sim_tree$tip.label, tree_bc$tip.label)
#     sim_bc <- if (length(drop) > 0) ape::drop.tip(sim_tree, drop) else sim_tree
#     
#     sp <- extract_spells_attained_with_censoring(sim_bc, tree_bc)
#     
#     if (nrow(sp) == 0) return(NULL)
#     
#     sp %>%
#       dplyr::mutate(sim = s) %>%
#       dplyr::filter(
#         is.finite(duration),
#         duration > 0,
#         is.finite(attain_age),
#         state %in% c("H", "W")
#       ) %>%
#       dplyr::mutate(
#         state = factor(state, levels = c("H", "W"))
#       )
#   }))
# })
# 
# if (is.null(spell_tbl_all) || nrow(spell_tbl_all) == 0) {
#   warning("[WARN] Step 6.3: no spells extracted; spell summaries not written.")
# } else {
#   saveRDS(
#     spell_tbl_all,
#     file.path(out_dir, "6.3_spells_attained_with_censoring_all.rds")
#   )
#   
#   persistence_summary <- spell_tbl_all %>%
#     dplyr::group_by(state) %>%
#     dplyr::summarise(
#       n_spells = dplyr::n(),
#       n_events = sum(event == 1L),
#       frac_censored = mean(event == 0L),
#       duration_median = stats::median(duration),
#       duration_q025 = q025(duration),
#       duration_q975 = q975(duration),
#       attain_age_median = stats::median(attain_age),
#       attain_age_q025 = q025(attain_age),
#       attain_age_q975 = q975(attain_age),
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       interpretation = dplyr::case_when(
#         state == "H" ~ "Descriptive SIMMAP spell summary for herbaceous-state attainments.",
#         state == "W" ~ "Descriptive SIMMAP spell summary for woody-state attainments."
#       )
#     )
#   
#   readr::write_csv(
#     persistence_summary,
#     file.path(out_dir, "6.3_persistence_spells_with_censoring_summary.csv")
#   )
#   
#   message("[OK] Wrote descriptive SIMMAP spell summary table.")
# }

### 6.4 Mk-implied sojourn times -----------------------------------------------

message("[INFO] Step 6.4: Mk-implied sojourn times")

if (!exists("Q")) {
  stop("[ERROR] Q not found; expected Q from earlier steps.")
}

if (is.null(rownames(Q)) || is.null(colnames(Q))) {
  rownames(Q) <- c("0", "1")
  colnames(Q) <- c("0", "1")
}

rate_out_H <- as.numeric(Q["0", "1"])
rate_out_W <- as.numeric(Q["1", "0"])

mk_sojourn_tbl <- tibble::tibble(
  state = c("H", "W"),
  rate_out = c(rate_out_H, rate_out_W),
  sojourn_mean   = 1 / c(rate_out_H, rate_out_W),
  sojourn_median = log(2) / c(rate_out_H, rate_out_W),
  sojourn_q025   = stats::qexp(0.025, rate = c(rate_out_H, rate_out_W)),
  sojourn_q975   = stats::qexp(0.975, rate = c(rate_out_H, rate_out_W))
) %>%
  dplyr::mutate(
    interpretation = dplyr::case_when(
      state == "H" ~ "Mk-implied waiting time until an herbaceous lineage shifts to woodiness.",
      state == "W" ~ "Mk-implied waiting time until a woody lineage shifts to herbaceousness."
    )
  )

readr::write_csv(
  mk_sojourn_tbl,
  file.path(out_dir, "6.4_mk_implied_sojourn_times.csv")
)

p_mk_sojourn <- ggplot(
  mk_sojourn_tbl,
  aes(x = sojourn_median, y = state, color = state)
) +
  geom_errorbarh(
    aes(xmin = sojourn_q025, xmax = sojourn_q975),
    height = 0,
    linewidth = 0.7
  ) +
  geom_point(size = 2.6) +
  scale_x_log10() +
  scale_color_manual(
    values = c("H" = col_wh, "W" = col_hw),
    guide = "none"
  ) +
  labs(
    x = "Mk-implied sojourn time (Myr)",
    y = NULL,
    title = paste0("Mk-implied persistence (best model = ", best_model, ")")
  ) +
  theme_classic()

ggsave(
  filename = file.path(out_dir, "6.4_mk_sojourn_times.png"),
  plot     = p_mk_sojourn,
  width    = 6.5,
  height   = 3.0,
  dpi      = 350
)

ggsave(
  filename = file.path(out_dir, "6.4_mk_sojourn_times.pdf"),
  plot     = p_mk_sojourn,
  width    = 6.5,
  height   = 3.0
)

message("[OK] Wrote Mk sojourn summary table and figure.")

### 6.5 Optional exploratory Kaplan–Meier diagnostics ---------------------------

if (isTRUE(RUN_EXPLORATORY_KM)) {
  message("[INFO] Step 6.5: optional exploratory Kaplan–Meier diagnostics")
  
  if (is.null(spell_tbl_all) || nrow(spell_tbl_all) == 0) {
    warning("[WARN] Step 6.5: spell_tbl_all is empty; skipping exploratory KM diagnostics.")
  } else if (!requireNamespace("survival", quietly = TRUE)) {
    warning("[WARN] Step 6.5: package 'survival' not installed; skipping exploratory KM diagnostics.")
  } else {
    suppressPackageStartupMessages(library(survival))
    
    sf_pool <- survival::survfit(
      survival::Surv(time = duration, event = event) ~ state,
      data = spell_tbl_all
    )
    
    km_df <- tibble::tibble(
      time = sf_pool$time,
      surv = sf_pool$surv,
      strata = rep(names(sf_pool$strata), sf_pool$strata)
    ) %>%
      dplyr::mutate(
        state = gsub("^state=", "", strata),
        state = factor(state, levels = c("H", "W"))
      )
    
    p_km_pool <- ggplot(km_df, aes(x = time, y = surv, color = state)) +
      geom_step(linewidth = 0.9) +
      scale_x_continuous(trans = "log10") +
      scale_color_manual(
        values = c("H" = col_wh, "W" = col_hw),
        name = "State"
      ) +
      labs(
        x = "Time since attainment (Myr)",
        y = "Exploratory persistence",
        title = "Exploratory Kaplan–Meier summary of SIMMAP spell durations"
      ) +
      theme_classic()
    
    ggsave(
      filename = file.path(out_dir, "6.5_exploratory_kaplan_meier_persistence_pooled.png"),
      plot     = p_km_pool,
      width    = 7.6,
      height   = 4.6,
      dpi      = 350
    )
  }
} else {
  message("[INFO] Step 6.5 skipped: exploratory Kaplan–Meier diagnostics are OFF by default.")
}

message("[OK] Step 6 complete.")


## 7. Combined figure: shift timing + tempo + propensity + sojourn -------------

message("[INFO] Step 7: combined shift timing + tempo + propensity + sojourn figure")

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("[ERROR] Package 'patchwork' is required for Step 7.")
}

stopifnot(exists("brassi_crown_age"))
stopifnot(exists("p_per_shift") || exists("p_per_switch"))
stopifnot(exists("p_rate"))
stopifnot(exists("p_prop"))
stopifnot(exists("p_mk_sojourn"))

if (!exists("p_per_shift") && exists("p_per_switch")) {
  p_per_shift <- p_per_switch
}

pad <- max(0.02 * brassi_crown_age, 0.25)
xlim_rev <- c(brassi_crown_age + pad, 0)

apply_common_time_scale <- function(p, ylim = NULL, clip = "on") {
  p +
    ggplot2::scale_x_reverse(
      limits = xlim_rev,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      ylim = ylim,
      clip = clip
    )
}

pA <- add_crown_label(
  apply_common_time_scale(p_per_shift, clip = "off"),
  y = max(y_lu$y_num) + 0.05
) +
  ggplot2::labs(title = "A") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.margin = ggplot2::margin(6, 10, 6, 6),
    legend.position = c(0.10, 0.10),
    legend.justification = c(0, 0),
    legend.background = ggplot2::element_rect(fill = "white", colour = NA),
    legend.title = ggplot2::element_text(size = 8),
    legend.text  = ggplot2::element_text(size = 8)
  )

pB <- add_crown_label(
  apply_common_time_scale(p_rate, ylim = c(0, ymax_show)),
  y = ymax_show
) +
  ggplot2::labs(title = "B") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.margin = ggplot2::margin(6, 6, 6, 6),
    axis.title.x = ggplot2::element_blank()
  )

pC <- add_crown_label(
  apply_common_time_scale(p_prop, ylim = c(0, ymax_show)),
  y = ymax_show
) +
  ggplot2::labs(title = "C", x = time_axis_lab) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.margin = ggplot2::margin(6, 6, 6, 6)
  )

pD <- p_mk_sojourn +
  ggplot2::labs(title = "D") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.margin = ggplot2::margin(6, 6, 6, 6)
  )

right_stack <- (pB / pC / pD) +
  patchwork::plot_layout(heights = c(1, 1, 0.55))

p_combo_7 <- (pA | right_stack) +
  patchwork::plot_layout(widths = c(1.45, 1.0))

ggplot2::ggsave(
  filename = file.path(out_dir, "7_combined_shift_timing_tempo_propensity_sojourn.png"),
  plot     = p_combo_7,
  width    = 12.5,
  height   = 10.2,
  dpi      = 400
)

ggplot2::ggsave(
  filename = file.path(out_dir, "7_combined_shift_timing_tempo_propensity_sojourn.pdf"),
  plot     = p_combo_7,
  width    = 12.5,
  height   = 10.2
)

message("[OK] Wrote combined Step 7 figure.")


## 8. Save objects for niche-coupling scripts ----------------------------------

saveRDS(
  list(
    tree_bc = tree_bc,
    tree_coded = tree_coded,
    woody_fac_bc = woody_fac_bc,
    asr_node_table = asr_node_table,
    asr_node_table_bc = asr_node_table_bc,
    edge_summary = edge_summary,
    edge_summary_bc = edge_summary_bc,
    change_tbl = change_tbl,
    change_tbl_bc = change_tbl_bc,
    confident_shifts_all = confident_shifts_all,
    diffuse_shift_signal_bc = diffuse_shift_signal_bc,
    shift_tbl = shift_tbl,
    per_shift_age = if (exists("per_shift_age")) per_shift_age else NULL,
    rate_summary = rate_summary,
    prop_summary = prop_summary,
    mk_sojourn_tbl = mk_sojourn_tbl,
    brassi_stem_age = brassi_stem_age,
    brassi_crown_age = brassi_crown_age
  ),
  file.path(out_dir, "8_objects.rds")
)

message("[DONE] 2f.StudyWoodinessEvolution.R complete.")
message("[INFO] Key figures:")
message("  ", file.path(out_dir, "5.5_circular_tree_BrassiCleom_confident_internal.pdf"))
message("  ", file.path(out_dir, "5.6_confident_internal_shift_ages_density.png"))
message("  ", file.path(out_dir, "5.7_confident_internal_shift_ages_by_shift.png"))
message("  ", file.path(out_dir, "6.1_6.2_tempo_propensity_2panel.png"))
message("  ", file.path(out_dir, "6.4_mk_sojourn_times.png"))
message("  ", file.path(out_dir, "7_combined_shift_timing_tempo_propensity_sojourn.png"))

