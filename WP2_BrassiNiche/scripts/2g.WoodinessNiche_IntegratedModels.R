#!/usr/bin/env Rscript

## 2g.WoodinessNiche_IntegratedModels.R
## Integrated non-phylogenetic and phylogenetic models of Brassicaceae woodiness
## as a function of climatic niche.
##
## Main goals:
##   1. Fit the same core model family to:
##      - full niche datasets (glm_full)
##      - matched niche+tree subsets (glm_matched)
##      - matched niche+tree subsets with phylogeny (phyloglm_matched)
##   2. Compare:
##      - subset effects   : glm_matched - glm_full
##      - phylogeny effects: phyloglm_matched - glm_matched
##   3. Retain a separate non-phylogenetic supertribe heterogeneity analysis
##      (ge10 only), focused on MCWD, BIO6, and elevation.
##
## Core models:
##   - m_clim
##   - m_clim_breadth
##   - m_clim_island
##   - m_clim_mainland
##
## Notes:
##   - Tree tips are SAMPLE IDs; species-level niche data are joined via binomial key.
##   - The full GLM analyses are not reduced to the phylogenetic subset.
##   - GLM_matched and phyloGLM_matched are fitted on identical complete-case subsets
##     per tier/model, so the effect of phylogeny is directly interpretable.

## 0. Prepare ----

out_dir <- "WP2_BrassiNiche/results_final/4_woodiness_niche_integrated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("[INFO] Output dir: ", normalizePath(out_dir, mustWork = FALSE))

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(phytools)
  library(patchwork)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

need_phylolm <- function() {
  if (!has_pkg("phylolm")) {
    stop("[ERROR] Package 'phylolm' is required. Install: install.packages('phylolm')")
  }
  TRUE
}

need_phytools_if_needed <- function() {
  if (!has_pkg("phytools")) {
    stop("[ERROR] Tree is not ultrametric and 'phytools' not installed. Install phytools or provide an ultrametric tree.")
  }
  TRUE
}

save_plot <- function(p, fn, w = 12, h = 8, dpi = 300) {
  png_file <- file.path(out_dir, paste0("PLOT_", fn, ".png"))
  pdf_file <- file.path(out_dir, paste0("PLOT_", fn, ".pdf"))
  
  ggsave(filename = png_file, plot = p, width = w, height = h, dpi = dpi)
  ggsave(filename = pdf_file, plot = p, width = w, height = h)
  
  message("[save_plot] wrote: ", png_file)
  message("[save_plot] wrote: ", pdf_file)
}

save_table <- function(x, fn) {
  out_file <- file.path(out_dir, paste0("TABLE_", fn, ".csv"))
  readr::write_csv(x, out_file)
  message("[save_table] wrote: ", out_file)
}

save_text <- function(lines, fn) {
  out_file <- file.path(out_dir, paste0("TABLE_", fn, ".txt"))
  cat(lines, file = out_file)
  message("[save_text] wrote: ", out_file)
}

tier_levels <- c("ge1", "ge5", "ge10")
as_tier <- function(x) factor(as.character(x), levels = tier_levels)

col_framework <- c(
  glm_full = "#7f7f7f",
  glm_matched = "#4C78A8",
  phyloglm_matched = "#F58518"
)

col_delta <- c(
  subset_effect = "#4C78A8",
  phylogeny_effect = "#F58518"
)

st_cols <- c(
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8",
  "Aethionemeae"      = "#ffed57"
)

include_special_cases <- FALSE
include_step10_supertribe <- TRUE
step10_min_n_per_supertribe <- 15

## 0.1 Inputs ----

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

cli$niche_ge1  <- "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge1.csv"
cli$niche_ge5  <- "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge5.csv"
cli$niche_ge10 <- "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge10.csv"

if (!file.exists(cli$tree_lsd2_nexus)) {
  stop("[ERROR] LSD2 Nexus tree not found: ", cli$tree_lsd2_nexus)
}
stopifnot(file.exists(cli$niche_ge1), file.exists(cli$niche_ge5), file.exists(cli$niche_ge10))

## 0.2 Helpers ----

normalize_growth <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  x <- stringr::str_to_upper(x)
  dplyr::case_when(
    x %in% c("W", "W/L", "W/T") ~ "W",
    x == "H"                    ~ "H",
    TRUE                        ~ NA_character_
  )
}

as_binom <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  m <- stringr::str_match(x, "^([^\\s]+)\\s+([^\\s]+)")
  genus <- m[, 2]
  spp   <- m[, 3]
  ifelse(is.na(genus) | is.na(spp), NA_character_, paste(genus, spp))
}

key_binom <- function(x) stringr::str_to_lower(as_binom(x))

taxon_key <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  
  key0 <- key_binom(x)
  if (!isTRUE(include_special_cases)) return(key0)
  
  is_bo_kale <- !is.na(x) & x %in% c(
    "Brassica oleracea (Canary Island Kale)",
    "Brassica oleracea (Jersey Kale)"
  )
  key0[is_bo_kale] <- stringr::str_to_lower(x[is_bo_kale])
  
  is_lobo <- !is.na(x) & stringr::str_detect(x, "^Lobularia\\s+canariensis\\b")
  if (any(is_lobo)) {
    m <- stringr::str_match(x[is_lobo], "^(Lobularia\\s+canariensis)(?:\\s+subsp\\.?\\s+([^\\s]+))?")
    base <- m[, 2]
    sube <- m[, 3]
    out  <- ifelse(is.na(sube), base, paste(base, "subsp", sube))
    key0[is_lobo] <- stringr::str_to_lower(out)
  }
  
  key0
}

zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

safe_keep_tip <- function(tr, tips) {
  tips2 <- intersect(tr$tip.label, tips)
  ape::keep.tip(tr, tips2)
}

get_ll <- function(fit) as.numeric(stats::logLik(fit))

mcfadden_r2 <- function(ll_full, ll_null) {
  if (!is.finite(ll_full) || !is.finite(ll_null) || ll_null == 0) return(NA_real_)
  1 - (ll_full / ll_null)
}

add_wald_ci <- function(tbl) {
  tbl %>%
    mutate(
      conf.low = estimate - 1.96 * se,
      conf.high = estimate + 1.96 * se
    )
}

assert_cols <- function(df, needed, label) {
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) stop("[ERROR] Missing columns in ", label, ": ", paste(miss, collapse = ", "))
  TRUE
}

safe_glm <- function(form, dat) {
  tryCatch(
    stats::glm(form, data = dat, family = stats::binomial()),
    error = function(e) {
      message("[WARN] glm failed: ", conditionMessage(e))
      NULL
    }
  )
}

safe_phyloglm <- function(form, dat, phy) {
  tryCatch(
    phylolm::phyloglm(
      formula = form,
      data    = dat,
      phy     = phy,
      method  = "logistic_MPLE",
      boot    = 0
    ),
    error = function(e) {
      message("[WARN] phyloglm failed: ", conditionMessage(e))
      NULL
    }
  )
}

is_special_bo_kale <- function(x) {
  stringr::str_squish(as.character(x)) %in% c(
    "Brassica oleracea (Jersey Kale)",
    "Brassica oleracea (Canary Island Kale)"
  )
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(
      estimate = NA_real_,
      p = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = length(x)
    ))
  }
  
  ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  
  out <- tibble(
    estimate = unname(ct$estimate),
    p = ct$p.value,
    n = length(x)
  )
  
  if (!is.null(ct$conf.int) && length(ct$conf.int) == 2) {
    out <- out %>%
      mutate(
        conf.low = ct$conf.int[1],
        conf.high = ct$conf.int[2]
      )
  } else {
    out <- out %>%
      mutate(
        conf.low = NA_real_,
        conf.high = NA_real_
      )
  }
  
  out
}

calc_vif_tbl <- function(df, preds) {
  preds <- unique(preds)
  preds <- preds[preds %in% names(df)]
  
  if (length(preds) == 0) {
    return(tibble(
      term = character(),
      vif = numeric(),
      tolerance = numeric(),
      r2_against_others = numeric(),
      n = integer()
    ))
  }
  
  d <- df %>%
    dplyr::select(all_of(preds)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(stats::complete.cases(.))
  
  if (nrow(d) < 5) {
    return(tibble(
      term = preds,
      vif = NA_real_,
      tolerance = NA_real_,
      r2_against_others = NA_real_,
      n = nrow(d)
    ))
  }
  
  out <- lapply(preds, function(v) {
    others <- setdiff(preds, v)
    
    if (length(others) == 0) {
      tibble(
        term = v,
        vif = 1,
        tolerance = 1,
        r2_against_others = 0,
        n = nrow(d)
      )
    } else {
      form <- stats::as.formula(
        paste(v, "~", paste(others, collapse = " + "))
      )
      fit <- stats::lm(form, data = d)
      r2 <- summary(fit)$r.squared
      tol <- 1 - r2
      vif <- ifelse(isTRUE(all.equal(tol, 0)), Inf, 1 / tol)
      
      tibble(
        term = v,
        vif = vif,
        tolerance = tol,
        r2_against_others = r2,
        n = nrow(d)
      )
    }
  })
  
  bind_rows(out)
}


## 1. Read tree + local metadata CSVs ----

message("[INFO] Step 1: Reading tree + local metadata CSVs")

tre <- ape::read.nexus(cli$tree_lsd2_nexus)

if (inherits(tre, "multiPhylo")) {
  message("[INFO] Nexus file contained multiple trees; using the first one.")
  tre <- tre[[1]]
}

message("[INFO] Read LSD2-calibrated tree from Nexus:")
message("[INFO]   tree: ", normalizePath(cli$tree_lsd2_nexus))
message("[INFO]   tips=", length(tre$tip.label), " ; Nnode=", tre$Nnode)

if (!ape::is.rooted(tre)) {
  warning("[WARN] Tree is not rooted. Check LSD2 output/rooting before downstream analyses.")
}

if (any(is.na(tre$edge.length))) {
  stop("[ERROR] Tree has NA edge lengths.")
}

if (any(tre$edge.length < 0, na.rm = TRUE)) {
  stop("[ERROR] Tree has negative edge lengths. Inspect LSD2 output.")
}

zero_edges <- sum(tre$edge.length == 0, na.rm = TRUE)
if (zero_edges > 0) {
  message("[INFO] Replacing ", zero_edges, " zero-length edges with 1e-10 Ma.")
  tre$edge.length[tre$edge.length == 0] <- 1e-10
}

tip_depths_raw <- ape::node.depth.edgelength(tre)[seq_along(tre$tip.label)]
message("[INFO] Raw LSD2 tree height: ", signif(max(tip_depths_raw, na.rm = TRUE), 6), " Ma")

message("[INFO] Reading specimen metadata from local CSV: ", cli$specimen_csv)
specimen_details <- readr::read_csv(cli$specimen_csv, show_col_types = FALSE)

message("[INFO] Reading species metadata from local CSV: ", cli$species_csv)
species_details <- readr::read_csv(cli$species_csv, show_col_types = FALSE)

specimen_details <- specimen_details %>%
  mutate(
    SAMPLE             = as.character(SAMPLE),
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    FAMILY             = as.character(FAMILY),
    loci_remaining     = suppressWarnings(as.numeric(loci_remaining))
  )
specimen_details$loci_remaining[is.na(specimen_details$loci_remaining)] <- 0

## 2. Build lookup tables ----

message("[INFO] Step 2: Building lookup tables")

species_lu_print <- species_details %>%
  transmute(
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    GROWTH_FORM        = normalize_growth(GROWTH_FORM)
  ) %>%
  filter(!is.na(SPECIES_NAME_PRINT)) %>%
  group_by(SPECIES_NAME_PRINT) %>%
  summarise(
    has_W = any(GROWTH_FORM == "W", na.rm = TRUE),
    has_H = any(GROWTH_FORM == "H", na.rm = TRUE),
    WOODY_STATE_PRINT = case_when(
      has_W & !has_H ~ "W",
      has_H & !has_W ~ "H",
      TRUE           ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(SPECIES_NAME_PRINT, WOODY_STATE_PRINT)

species_lu_binom <- species_details %>%
  transmute(
    BIONOMIAL   = stringr::str_squish(as.character(BIONOMIAL)),
    GROWTH_FORM = normalize_growth(GROWTH_FORM)
  ) %>%
  filter(!is.na(BIONOMIAL)) %>%
  group_by(BIONOMIAL) %>%
  summarise(
    has_W = any(GROWTH_FORM == "W", na.rm = TRUE),
    has_H = any(GROWTH_FORM == "H", na.rm = TRUE),
    WOODY_STATE_BINOM = case_when(
      has_W & !has_H ~ "W",
      has_H & !has_W ~ "H",
      TRUE           ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(BIONOMIAL, WOODY_STATE_BINOM) %>%
  mutate(BIONOMIAL = stringr::str_squish(BIONOMIAL))

specimen_details <- specimen_details %>%
  mutate(species_binom = as_binom(SPECIES_NAME_PRINT)) %>%
  left_join(species_lu_print, by = "SPECIES_NAME_PRINT") %>%
  left_join(species_lu_binom, by = c("species_binom" = "BIONOMIAL")) %>%
  mutate(
    WOODY_STATE = coalesce(WOODY_STATE_PRINT, WOODY_STATE_BINOM),
    woody_flag  = case_when(
      WOODY_STATE == "W" ~ "W",
      WOODY_STATE == "H" ~ "H",
      TRUE               ~ NA_character_
    )
  )

tip_tbl <- specimen_details %>%
  dplyr::select(SAMPLE, SPECIES_NAME_PRINT, FAMILY, woody_flag, loci_remaining) %>%
  distinct() %>%
  filter(SAMPLE %in% tre$tip.label)

brassi_tips <- tip_tbl %>%
  filter(FAMILY == "Brassicaceae") %>%
  pull(SAMPLE) %>%
  unique()

tree_b <- ape::keep.tip(tre, intersect(tre$tip.label, brassi_tips))
message("[INFO] Brassicaceae-only tree: tips=", length(tree_b$tip.label))

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

tree_b$edge.length[tree_b$edge.length == 0] <- 1e-10
tree_b <- make_ultrametric_by_extending_tips(tree_b)

stopifnot(
  ape::is.ultrametric(
    tree_b,
    tol = max(1e-6, max(ape::node.depth.edgelength(tree_b)) * 1e-5)
  )
)

tip_key <- tip_tbl %>%
  filter(SAMPLE %in% tree_b$tip.label) %>%
  mutate(
    SPECIES_NAME_PRINT = stringr::str_replace_all(SPECIES_NAME_PRINT, "_", " "),
    key = taxon_key(SPECIES_NAME_PRINT),
    woody01 = case_when(
      woody_flag == "W" ~ 1L,
      woody_flag == "H" ~ 0L,
      TRUE              ~ NA_integer_
    )
  ) %>%
  dplyr::select(SAMPLE, SPECIES_NAME_PRINT, key, woody01, woody_flag, loci_remaining) %>%
  filter(!is.na(key)) %>%
  distinct()

message("[INFO] tip_key summary:")
message("  rows:        ", nrow(tip_key))
message("  unique keys: ", dplyr::n_distinct(tip_key$key))
message("  woody coded: ", sum(!is.na(tip_key$woody01)))
message("  missing:     ", sum(is.na(tip_key$woody01)))

## Species lookup for niche datasets, matching 2e logic
req_cols_species <- c(
  "SPECIES_NAME_PRINT", "GROWTH_FORM", "SUPERTRIBE", "TRIBE_FULL",
  "SUBSP_VAR", "ACCEPTED", "Island_endemic"
)
miss_cols_species <- setdiff(req_cols_species, names(species_details))
if (length(miss_cols_species) > 0) {
  stop("[ERROR] Missing required columns in species_details: ", paste(miss_cols_species, collapse = ", "))
}

is_empty_cell <- function(x) {
  x <- stringr::str_squish(as.character(x))
  is.na(x) | x == ""
}

EXC_BINOM <- "Lobularia canariensis"
EXC_SPECIES_PRINT <- c(
  "Brassica oleracea (Jersey Kale)",
  "Brassica oleracea (Canary Island Kale)"
)
EXC_SPECIES_PRINT_NORM <- stringr::str_squish(EXC_SPECIES_PRINT)

species_sheet_filt <- species_details %>%
  mutate(
    ACCEPTED   = stringr::str_squish(as.character(ACCEPTED)),
    SUBSP_VAR  = stringr::str_squish(as.character(SUBSP_VAR)),
    TRIBE_FULL = stringr::str_squish(as.character(TRIBE_FULL)),
    SUPERTRIBE = stringr::str_squish(as.character(SUPERTRIBE)),
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    binom = stringr::str_c(
      stringr::word(SPECIES_NAME_PRINT, 1),
      stringr::word(SPECIES_NAME_PRINT, 2),
      sep = " "
    ),
    is_exc_print = SPECIES_NAME_PRINT %in% EXC_SPECIES_PRINT_NORM,
    is_exc_binom = binom == EXC_BINOM
  ) %>%
  filter(
    is_exc_print |
      (ACCEPTED == "Y" & (is_empty_cell(SUBSP_VAR) | is_exc_binom))
  )

species_lu <- species_sheet_filt %>%
  transmute(
    species = gsub(" ", "_", as.character(SPECIES_NAME_PRINT)),
    key = taxon_key(SPECIES_NAME_PRINT),
    woody_raw = GROWTH_FORM,
    woody_state = normalize_growth(woody_raw),
    woody_bin = case_when(
      woody_state == "W" ~ 1L,
      woody_state == "H" ~ 0L,
      TRUE ~ NA_integer_
    ),
    tribe_full     = dplyr::na_if(stringr::str_squish(TRIBE_FULL), ""),
    supertribe_raw = dplyr::na_if(stringr::str_squish(SUPERTRIBE), ""),
    is_island = dplyr::if_else(
      stringr::str_to_upper(stringr::str_squish(as.character(Island_endemic))) == "Y",
      1L, 0L, missing = 0L
    )
  ) %>%
  mutate(
    supertribe = case_when(
      tribe_full == "Aethionemeae" ~ "Aethionemeae",
      !is.na(supertribe_raw)       ~ supertribe_raw,
      TRUE                         ~ NA_character_
    ),
    supertribe = factor(supertribe, levels = names(st_cols)),
    island_f = factor(is_island, levels = c(0, 1), labels = c("mainland", "island"))
  ) %>%
  dplyr::select(species, key, supertribe, woody_state, woody_bin, is_island, island_f) %>%
  distinct(species, .keep_all = TRUE)

if (n_distinct(species_lu$species) != nrow(species_lu)) {
  dup <- species_lu %>% count(species) %>% filter(n > 1)
  message("[ERROR] species_lu has duplicate species keys. First 20:")
  print(utils::head(dup, 20))
  stop("[ERROR] Fix duplicates in species_lu before continuing.")
}

## 3. Read niche tables ----

message("[INFO] Step 3: Reading niche tables")

niche_ge1  <- readr::read_csv(cli$niche_ge1,  show_col_types = FALSE)
niche_ge5  <- readr::read_csv(cli$niche_ge5,  show_col_types = FALSE)
niche_ge10 <- readr::read_csv(cli$niche_ge10, show_col_types = FALSE)

niche_species_col <- "species"
stopifnot(niche_species_col %in% names(niche_ge1),
          niche_species_col %in% names(niche_ge5),
          niche_species_col %in% names(niche_ge10))

make_niche_key <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()
}

niche_add_key <- function(df, tier_label) {
  df %>%
    mutate(
      tier = tier_label,
      key = make_niche_key(.data[[niche_species_col]]),
      species = as.character(.data[[niche_species_col]])
    )
}

niche_all <- bind_rows(
  niche_add_key(niche_ge1,  "ge1"),
  niche_add_key(niche_ge5,  "ge5"),
  niche_add_key(niche_ge10, "ge10")
) %>%
  mutate(tier = factor(tier, levels = tier_levels))

message("[INFO] Niche rows by tier:")
print(table(niche_all$tier, useNA = "ifany"))

## 4. Merge niche tables with species lookup (full-data base) ----

message("[INFO] Step 4: Building full-data tier tables")


# 4b. Island vs mainland summaries ----

message("[INFO] Step 4b: Island vs mainland summaries")

## Build full datasets per tier (same base as modelling)
tier_full_tbl <- niche_all %>%
  dplyr::left_join(
    species_lu %>% dplyr::select(-key),
    by = "species"
  ) %>%
  dplyr::filter(!is.na(woody_bin)) %>%
  dplyr::mutate(
    woody_state = case_when(
      woody_bin == 1L ~ "W",
      woody_bin == 0L ~ "H",
      TRUE ~ NA_character_
    ),
    island_f = factor(
      ifelse(is_island == 1L, "island_endemic", "mainland"),
      levels = c("mainland", "island_endemic")
    )
  )

## --- Table 1: Distribution of woodiness across island vs mainland

tbl_wood_island <- tier_full_tbl %>%
  group_by(tier, island_f, woody_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = woody_state,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total = H + W,
    p_woody = W / total
  ) %>%
  arrange(tier, island_f)

save_table(tbl_wood_island, "4b_01_distribution_woodiness_island_mainland")

## --- Table 2: Proportion of island endemics within woody vs herbaceous

tbl_island_within <- tier_full_tbl %>%
  group_by(tier, woody_state, island_f) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = island_f,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total = mainland + island_endemic,
    p_island = island_endemic / total
  ) %>%
  arrange(tier, woody_state)

save_table(tbl_island_within, "4b_02_proportion_island_within_woodiness")



## 5. Build matched tier datasets for tree-based analyses ----

message("[INFO] Step 5: Building matched tier datasets")

build_tier_matched <- function(tier) {
  df_n <- niche_all %>% filter(tier == !!tier)
  
  df0 <- df_n %>%
    dplyr::left_join(
      species_lu %>% dplyr::select(-key),
      by = "species"
    )
  
  if (!"key" %in% names(df0)) {
    stop("[ERROR] Column 'key' missing after join in build_tier_matched().")
  }
  
  if ("status" %in% names(df0)) {
    df0 <- df0 %>% filter(status == "ok")
  }
  
  df0 <- df0 %>%
    dplyr::filter(!is.na(woody_bin))
  
  df_joined <- tip_key %>%
    dplyr::inner_join(df0, by = "key") %>%
    filter(!is.na(woody01))
  
  df <- df_joined %>%
    group_by(key) %>%
    arrange(is_special_bo_kale(SPECIES_NAME_PRINT), desc(loci_remaining), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  matched_tips <- intersect(tree_b$tip.label, df$SAMPLE)
  tr <- safe_keep_tip(tree_b, matched_tips)
  
  message(
    "[INFO] ", tier,
    ": full_model_rows=", nrow(df0),
    " matched_species=", nrow(df),
    " pruned_tips=", length(tr$tip.label)
  )
  
  if (length(matched_tips) < 25) {
    stop("[ERROR] Too few matched tips for ", tier, ". Fix keying before modelling.")
  }
  
  list(
    full_model = df0,
    matched_model = df,
    tree = tr
  )
}

tier_objs <- list(
  ge1  = build_tier_matched("ge1"),
  ge5  = build_tier_matched("ge5"),
  ge10 = build_tier_matched("ge10")
)

message("[INFO] Checking species-vs-tip woodiness agreement in matched datasets")
for (nm in tier_levels) {
  dchk <- tier_objs[[nm]]$matched_model
  if (all(c("woody01", "woody_bin") %in% names(dchk))) {
    n_mismatch <- sum(!is.na(dchk$woody01) & !is.na(dchk$woody_bin) & dchk$woody01 != dchk$woody_bin)
    message("  ", nm, ": mismatches woody01 vs woody_bin = ", n_mismatch)
  }
}

## 6. Define predictors and model family ----

message("[INFO] Step 6: Defining predictors and model family")

cols <- list(
  y = "woody_bin",
  MCWD_med = "MCWD_med",
  AI_med   = "AI_med",
  elev_med = "elev_med",
  BIO10_med = "BIO10_med",
  BIO6_med = "BIO6_med",
  MCWD_mad = "MCWD_mad",
  BIO10_mad = "BIO10_mad",
  BIO6_mad = "BIO6_mad",
  elev_mad = "elev_mad",
  is_island = "is_island"
)

assert_cols(
  niche_all,
  c(
    "tier", "key",
    cols$MCWD_med, cols$AI_med, cols$elev_med, cols$BIO10_med, cols$BIO6_med,
    cols$MCWD_mad, cols$BIO10_mad, cols$BIO6_mad, cols$elev_mad
  ),
  "niche_all"
)

model_set <- list(
  list(
    name = "m_clim",
    terms = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med),
    mainland_only = FALSE
  ),
  list(
    name = "m_clim_breadth",
    terms = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med,
              cols$MCWD_mad, cols$BIO10_mad, cols$BIO6_mad, cols$elev_mad),
    mainland_only = FALSE
  ),
  list(
    name = "m_clim_island",
    terms = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med, cols$is_island),
    mainland_only = FALSE
  ),
  list(
    name = "m_clim_mainland",
    terms = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med),
    mainland_only = TRUE
  )
)

model_levels <- vapply(model_set, `[[`, character(1), "name")
framework_levels <- c("glm_full", "glm_matched", "phyloglm_matched")

message("[INFO] Models: ", paste(model_levels, collapse = ", "))

## 6b. Predictor dependence / collinearity checks ----

message("[INFO] Step 6b: Predictor dependence / collinearity checks")

predictor_sets <- list(
  m_clim = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med),
  m_clim_breadth = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med,
                     cols$MCWD_mad, cols$BIO10_mad, cols$BIO6_mad, cols$elev_mad),
  m_clim_island = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med, cols$is_island),
  m_clim_mainland = c(cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med)
)

## 6b.1 Pairwise predictor correlations across tier tables
pairwise_cor_tbl <- bind_rows(lapply(tier_levels, function(tt) {
  
  df_tt <- niche_all %>%
    dplyr::filter(tier == tt) %>%
    dplyr::left_join(species_lu %>% dplyr::select(-key), by = "species")
  
  vars_here <- unique(c(
    cols$MCWD_med, cols$BIO10_med, cols$BIO6_med, cols$elev_med,
    cols$MCWD_mad, cols$BIO10_mad, cols$BIO6_mad, cols$elev_mad,
    cols$is_island
  ))
  
  vars_here <- vars_here[vars_here %in% names(df_tt)]
  
  combs <- utils::combn(vars_here, 2, simplify = FALSE)
  
  bind_rows(lapply(combs, function(vv) {
    x <- suppressWarnings(as.numeric(df_tt[[vv[1]]]))
    y <- suppressWarnings(as.numeric(df_tt[[vv[2]]]))
    
    sp <- safe_cor_test(x, y, method = "spearman")
    pe <- safe_cor_test(x, y, method = "pearson")
    
    tibble(
      tier = tt,
      var1 = vv[1],
      var2 = vv[2],
      spearman_rho = sp$estimate,
      spearman_p = sp$p,
      spearman_n = sp$n,
      pearson_r = pe$estimate,
      pearson_p = pe$p,
      pearson_n = pe$n
    )
  }))
}))

save_table(pairwise_cor_tbl, "6b_01_predictor_pairwise_correlations")

## 6b.2 Focused MCWD vs BIO6 tests for full and matched datasets
focus_mcwd_bio6_tbl <- dplyr::bind_rows(lapply(tier_levels, function(tt) {
  
  obj_t <- tier_objs[[tt]]
  
  dat_list <- list(
    glm_full_input = obj_t$full_model,
    matched_input  = obj_t$matched_model
  )
  
  bind_rows(lapply(names(dat_list), function(ff) {
    dd <- dat_list[[ff]]
    
    if (!all(c(cols$MCWD_med, cols$BIO6_med) %in% names(dd))) {
      return(NULL)
    }
    
    x <- suppressWarnings(as.numeric(dd[[cols$MCWD_med]]))
    y <- suppressWarnings(as.numeric(dd[[cols$BIO6_med]]))
    
    sp <- safe_cor_test(x, y, method = "spearman")
    pe <- safe_cor_test(x, y, method = "pearson")
    
    tibble(
      tier = tt,
      dataset = ff,
      var1 = cols$MCWD_med,
      var2 = cols$BIO6_med,
      spearman_rho = sp$estimate,
      spearman_p = sp$p,
      spearman_n = sp$n,
      pearson_r = pe$estimate,
      pearson_p = pe$p,
      pearson_n = pe$n
    )
  }))
}))

save_table(focus_mcwd_bio6_tbl, "6b_02_focus_MCWD_BIO6_correlations")

## 6b.3 VIF for the exact predictor sets used in each model
vif_tbl <- dplyr::bind_rows(lapply(tier_levels, function(tt) {
  
  obj_t <- tier_objs[[tt]]
  
  dataset_list <- list(
    glm_full_input = obj_t$full_model,
    matched_input  = obj_t$matched_model
  )
  
  dplyr::bind_rows(lapply(names(predictor_sets), function(mm) {
    preds <- predictor_sets[[mm]]
    
    bind_rows(lapply(names(dataset_list), function(dd_name) {
      dd <- dataset_list[[dd_name]]
      
      if (mm == "m_clim_mainland" && "is_island" %in% names(dd)) {
        dd <- dd %>% dplyr::filter(is_island == 0L)
      }
      
      vtbl <- calc_vif_tbl(dd, preds)
      
      if (nrow(vtbl) == 0) return(NULL)
      
      vtbl %>%
        mutate(
          tier = tt,
          dataset = dd_name,
          model = mm,
          vif_flag = case_when(
            !is.finite(vif) ~ "non-finite",
            vif >= 10 ~ "high",
            vif >= 5  ~ "moderate",
            TRUE ~ "low"
          )
        ) %>%
        dplyr::select(tier, dataset, model, term, n, r2_against_others, tolerance, vif, vif_flag)
    }))
  }))
}))

save_table(vif_tbl, "6b_03_predictor_vif_by_tier_dataset_model")

## 6b.4 Compact summary text for quick interpretation
vif_summary_txt <- c(
  "Predictor dependence summary",
  "============================",
  "",
  paste0("Rows in pairwise correlation table: ", nrow(pairwise_cor_tbl)),
  paste0("Rows in focused MCWD-BIO6 table: ", nrow(focus_mcwd_bio6_tbl)),
  paste0("Rows in VIF table: ", nrow(vif_tbl)),
  "",
  "VIF interpretation used here:",
  "  < 5   : low collinearity",
  "  5-10  : moderate collinearity",
  "  >= 10 : high collinearity",
  "",
  paste0("Max VIF observed: ", round(max(vif_tbl$vif, na.rm = TRUE), 3)),
  paste0("Number of VIFs >= 5: ", sum(vif_tbl$vif >= 5, na.rm = TRUE)),
  paste0("Number of VIFs >= 10: ", sum(vif_tbl$vif >= 10, na.rm = TRUE))
)

save_text(vif_summary_txt, "6b_04_predictor_dependence_summary")

## 7. Fit integrated model family ----

message("[INFO] Step 7: Fitting integrated model family")
need_phylolm()

prep_predictors <- function(dsub, preds) {
  dsub <- dsub %>%
    mutate(across(all_of(preds), ~ suppressWarnings(as.numeric(.x))))
  
  for (cc in preds) {
    if (cc == "is_island") next
    dsub[[cc]] <- zscore(dsub[[cc]])
  }
  
  dsub
}

make_formula <- function(terms, response = "woody_bin") {
  stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

tidy_glm <- function(fit) {
  s <- summary(fit)
  coefs <- as.data.frame(s$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  
  coefs %>%
    rename(
      estimate = Estimate,
      se = `Std. Error`,
      z = `z value`,
      p = `Pr(>|z|)`
    ) %>%
    add_wald_ci() %>%
    mutate(
      OR = exp(estimate),
      OR_lo = exp(conf.low),
      OR_hi = exp(conf.high),
      model_type = "glm",
      logLik = as.numeric(stats::logLik(fit)),
      AIC = stats::AIC(fit)
    )
}

tidy_phy <- function(fit) {
  s <- summary(fit)
  coefs <- as.data.frame(s$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  
  names(coefs) <- sub("^StdErr$", "se", names(coefs))
  names(coefs) <- sub("^z\\.value$", "z", names(coefs))
  names(coefs) <- sub("^p\\.value$", "p", names(coefs))
  
  ll_obj <- tryCatch(stats::logLik(fit), error = function(e) NA)
  ll <- as.numeric(ll_obj)
  
  k <- suppressWarnings(as.numeric(attr(ll_obj, "df")))
  if (!is.finite(k)) k <- nrow(coefs) + 1
  
  coefs %>%
    rename(estimate = Estimate) %>%
    add_wald_ci() %>%
    mutate(
      OR = exp(estimate),
      OR_lo = exp(conf.low),
      OR_hi = exp(conf.high),
      model_type = "phyloglm",
      logLik = ll,
      AIC = -2 * ll + 2 * k
    )
}

fit_glm_framework <- function(df, model_def, framework_name) {
  terms <- model_def$terms
  
  missing_preds <- setdiff(c("woody_bin", terms), names(df))
  if (length(missing_preds) > 0) {
    return(list(ok = FALSE, reason = paste0("Missing predictor columns: ", paste(missing_preds, collapse = ", "))))
  }
  
  dsub <- df
  
  if (isTRUE(model_def$mainland_only)) {
    dsub <- dsub %>% filter(is_island == 0L)
  }
  
  dsub <- dsub %>%
    dplyr::select(any_of(unique(c("species", "key", "woody_bin", "is_island", terms)))) %>%
    filter(stats::complete.cases(.))
  
  if (nrow(dsub) < 25) {
    return(list(ok = FALSE, reason = "Too few complete cases (<25)."))
  }
  
  dsub <- prep_predictors(dsub, terms)
  form <- make_formula(terms, response = "woody_bin")
  fit <- safe_glm(form, dsub)
  if (is.null(fit)) {
    return(list(ok = FALSE, reason = "glm failed"))
  }
  
  list(
    ok = TRUE,
    framework = framework_name,
    n = nrow(dsub),
    formula = form,
    data = dsub,
    fit = fit
  )
}

fit_phyloglm_framework <- function(df, tree, model_def) {
  terms <- model_def$terms
  
  missing_preds <- setdiff(c("SAMPLE", "woody_bin", terms), names(df))
  if (length(missing_preds) > 0) {
    return(list(ok = FALSE, reason = paste0("Missing predictor columns: ", paste(missing_preds, collapse = ", "))))
  }
  
  dsub <- df
  
  if (isTRUE(model_def$mainland_only)) {
    dsub <- dsub %>% filter(is_island == 0L)
  }
  
  dsub <- dsub %>%
    dplyr::select(any_of(unique(c("SAMPLE", "species", "key", "woody_bin", "is_island", terms)))) %>%
    filter(stats::complete.cases(.))
  
  if (nrow(dsub) < 25) {
    return(list(ok = FALSE, reason = "Too few complete cases (<25)."))
  }
  
  keep_tips <- intersect(tree$tip.label, dsub$SAMPLE)
  tr <- safe_keep_tip(tree, keep_tips)
  
  dsub <- dsub %>% filter(SAMPLE %in% tr$tip.label)
  if (nrow(dsub) < 25 || length(tr$tip.label) < 25) {
    return(list(ok = FALSE, reason = "Too few tips after pruning (<25)."))
  }
  
  dsub <- prep_predictors(dsub, terms)
  
  dsub <- as.data.frame(dsub)
  rownames(dsub) <- dsub$SAMPLE
  dsub <- dsub[tr$tip.label, , drop = FALSE]
  
  form <- make_formula(terms, response = "woody_bin")
  
  fit <- safe_phyloglm(form, dsub, tr)
  if (is.null(fit)) {
    return(list(ok = FALSE, reason = "phyloglm failed"))
  }
  
  list(
    ok = TRUE,
    framework = "phyloglm_matched",
    n = nrow(dsub),
    formula = form,
    data = dsub,
    tree = tr,
    fit = fit
  )
}

all_results <- list()
coef_tbl <- list()
meta_tbl <- list()
fit_index <- tibble()

for (tier in tier_levels) {
  obj_t <- tier_objs[[tier]]
  
  for (m in model_set) {
    model_name <- m$name
    
    ## glm_full
    res_full <- fit_glm_framework(obj_t$full_model, m, "glm_full")
    key_full <- paste(tier, model_name, "glm_full", sep = "__")
    all_results[[key_full]] <- res_full
    
    if (isTRUE(res_full$ok)) {
      coef_tbl[[paste0(key_full, "__coef")]] <- tidy_glm(res_full$fit) %>%
        mutate(tier = tier, model = model_name, framework = "glm_full")
      
      meta_tbl[[paste0(key_full, "__meta")]] <- tibble(
        tier = tier,
        model = model_name,
        framework = "glm_full",
        n = res_full$n,
        formula = paste(deparse(res_full$formula), collapse = ""),
        logLik = as.numeric(logLik(res_full$fit)),
        AIC = AIC(res_full$fit),
        alpha = NA_real_
      )
    }
    
    fit_index <- dplyr::bind_rows(
      fit_index,
      tibble(
        tier = tier,
        model = model_name,
        framework = "glm_full",
        ok = isTRUE(res_full$ok),
        n = if (isTRUE(res_full$ok)) res_full$n else NA_integer_,
        reason = if (isTRUE(res_full$ok)) NA_character_ else res_full$reason
      )
    )
    
    ## glm_matched
    res_match <- fit_glm_framework(obj_t$matched_model, m, "glm_matched")
    key_match <- paste(tier, model_name, "glm_matched", sep = "__")
    all_results[[key_match]] <- res_match
    
    if (isTRUE(res_match$ok)) {
      coef_tbl[[paste0(key_match, "__coef")]] <- tidy_glm(res_match$fit) %>%
        mutate(tier = tier, model = model_name, framework = "glm_matched")
      
      meta_tbl[[paste0(key_match, "__meta")]] <- tibble(
        tier = tier,
        model = model_name,
        framework = "glm_matched",
        n = res_match$n,
        formula = paste(deparse(res_match$formula), collapse = ""),
        logLik = as.numeric(logLik(res_match$fit)),
        AIC = AIC(res_match$fit),
        alpha = NA_real_
      )
    }
    
    fit_index <- dplyr::bind_rows(
      fit_index,
      tibble(
        tier = tier,
        model = model_name,
        framework = "glm_matched",
        ok = isTRUE(res_match$ok),
        n = if (isTRUE(res_match$ok)) res_match$n else NA_integer_,
        reason = if (isTRUE(res_match$ok)) NA_character_ else res_match$reason
      )
    )
    
    ## phyloglm_matched
    res_phy <- fit_phyloglm_framework(obj_t$matched_model, obj_t$tree, m)
    key_phy <- paste(tier, model_name, "phyloglm_matched", sep = "__")
    all_results[[key_phy]] <- res_phy
    
    if (isTRUE(res_phy$ok)) {
      coef_tbl[[paste0(key_phy, "__coef")]] <- tidy_phy(res_phy$fit) %>%
        mutate(tier = tier, model = model_name, framework = "phyloglm_matched")
      
      alpha_now <- NA_real_
      if (!is.null(res_phy$fit$alpha)) alpha_now <- as.numeric(res_phy$fit$alpha)
      
      meta_tbl[[paste0(key_phy, "__meta")]] <- tibble(
        tier = tier,
        model = model_name,
        framework = "phyloglm_matched",
        n = res_phy$n,
        formula = paste(deparse(res_phy$formula), collapse = ""),
        logLik = as.numeric(logLik(res_phy$fit)),
        AIC = unique(tidy_phy(res_phy$fit)$AIC),
        alpha = alpha_now
      )
    }
    
    fit_index <- dplyr::bind_rows(
      fit_index,
      tibble(
        tier = tier,
        model = model_name,
        framework = "phyloglm_matched",
        ok = isTRUE(res_phy$ok),
        n = if (isTRUE(res_phy$ok)) res_phy$n else NA_integer_,
        reason = if (isTRUE(res_phy$ok)) NA_character_ else res_phy$reason
      )
    )
    
    message("[INFO] Done: ", tier, " | ", model_name)
  }
}

coef_tbl2 <- if (length(coef_tbl) == 0) tibble() else bind_rows(coef_tbl)
meta_tbl2 <- if (length(meta_tbl) == 0) tibble() else bind_rows(meta_tbl)

fit_index <- fit_index %>%
  mutate(
    tier = factor(tier, levels = tier_levels),
    model = factor(model, levels = model_levels),
    framework = factor(framework, levels = framework_levels)
  ) %>%
  arrange(tier, model, framework)

coef_tbl2 <- coef_tbl2 %>%
  mutate(
    tier = factor(tier, levels = tier_levels),
    model = factor(model, levels = model_levels),
    framework = factor(framework, levels = framework_levels)
  ) %>%
  arrange(tier, model, framework, term)

meta_tbl2 <- meta_tbl2 %>%
  mutate(
    tier = factor(tier, levels = tier_levels),
    model = factor(model, levels = model_levels),
    framework = factor(framework, levels = framework_levels)
  ) %>%
  arrange(tier, model, framework)

save_table(fit_index, "7_01_model_index")
save_table(coef_tbl2, "7_02_model_coefficients_master")
save_table(meta_tbl2, "7_03_model_meta_master")

## 8. Comparison tables ----

message("[INFO] Step 8: Building comparison tables")

coef_compare <- coef_tbl2 %>%
  filter(term != "(Intercept)") %>%
  dplyr::select(tier, model, framework, term, estimate, OR, conf.low, conf.high, AIC, p, se) %>%
  pivot_wider(
    names_from = framework,
    values_from = c(estimate, OR, conf.low, conf.high, AIC, p, se),
    names_sep = "__"
  )

subset_effect_tbl <- coef_compare %>%
  mutate(
    delta_beta_subset = estimate__glm_matched - estimate__glm_full,
    delta_OR_subset = OR__glm_matched / OR__glm_full
  ) %>%
  dplyr::select(
    tier, model, term,
    estimate__glm_full, estimate__glm_matched,
    conf.low__glm_full, conf.high__glm_full,
    conf.low__glm_matched, conf.high__glm_matched,
    OR__glm_full, OR__glm_matched,
    p__glm_full, p__glm_matched,
    delta_beta_subset, delta_OR_subset
  ) %>%
  arrange(tier, model, term)

phylogeny_effect_tbl <- coef_compare %>%
  mutate(
    delta_beta_phylogeny = estimate__phyloglm_matched - estimate__glm_matched,
    delta_OR_phylogeny = OR__phyloglm_matched / OR__glm_matched
  ) %>%
  dplyr::select(
    tier, model, term,
    estimate__glm_matched, estimate__phyloglm_matched,
    conf.low__glm_matched, conf.high__glm_matched,
    conf.low__phyloglm_matched, conf.high__phyloglm_matched,
    OR__glm_matched, OR__phyloglm_matched,
    p__glm_matched, p__phyloglm_matched,
    delta_beta_phylogeny, delta_OR_phylogeny
  ) %>%
  arrange(tier, model, term)

save_table(subset_effect_tbl, "8_01_model_comparison_subset_effect")
save_table(phylogeny_effect_tbl, "8_02_model_comparison_phylogeny_effect")

## 9. Plots of side-by-side coefficients ----

message("[INFO] Step 9: Plotting side-by-side coefficient summaries")

term_levels <- c(
  "MCWD_med", "MCWD_mad",
  "BIO6_med", "BIO6_mad",
  "BIO10_med", "BIO10_mad",
  "elev_med", "elev_mad",
  "is_island"
)

term_labels <- c(
  MCWD_med  = "Drought (median)",
  MCWD_mad  = "Drought niche breadth",
  BIO6_med  = "Frost (median)",
  BIO6_mad  = "Frost niche breadth",
  BIO10_med = "Heat (median)",
  BIO10_mad = "Heat niche breadth",
  elev_med  = "Elevation (median)",
  elev_mad  = "Elevation niche breadth",
  is_island = "Island endemic"
)

model_labels <- c(
  m_clim         = "Climate medians",
  m_clim_breadth = "Climate medians + breadth",
  m_clim_island = "Climate medians + island",
  m_clim_mainland = "Climate medians, mainland only"
)

framework_levels <- c("glm_full", "glm_matched", "phyloglm_matched")

framework_labels <- c(
  glm_full = "Full GLM",
  glm_matched = "Matched GLM",
  phyloglm_matched = "Matched phyloGLM"
)

if (nrow(coef_tbl2) > 0) {
  
  plot_df <- coef_tbl2 %>%
    dplyr::filter(term %in% term_levels) %>%
    dplyr::mutate(
      term = factor(term, levels = term_levels),
      term_label = dplyr::recode(as.character(term), !!!term_labels),
      term_label = factor(term_label, levels = rev(unname(term_labels[term_levels]))),
      framework = factor(framework, levels = framework_levels),
      model = factor(model, levels = model_levels)
    )
  
  n_all_tbl <- meta_tbl2 %>%
    dplyr::mutate(
      framework = factor(framework, levels = framework_levels),
      tier = factor(tier, levels = tier_levels),
      model = factor(model, levels = model_levels)
    ) %>%
    dplyr::arrange(tier, model, framework) %>%
    dplyr::group_by(tier, model) %>%
    dplyr::summarise(
      label_full  = paste0("n full = ", n[framework == "glm_full"]),
      label_match = paste0("n matched = ", n[framework == "glm_matched"]),
      label_phylo = paste0("n phylo = ", n[framework == "phyloglm_matched"]),
      .groups = "drop"
    )
  
  p_eff <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = estimate, y = term_label, color = framework)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      height = 0.18,
      linewidth = 0.8,
      position = ggplot2::position_dodge(width = 0.65, reverse = TRUE)
    ) +
    ggplot2::geom_point(
      size = 2.4,
      position = ggplot2::position_dodge(width = 0.65, reverse = TRUE)
    ) +
    ggplot2::facet_grid(tier ~ model, scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(
      values = col_framework,
      drop = FALSE,
      labels = framework_labels
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 13),
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      panel.border = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.45),
      panel.spacing = grid::unit(1.1, "lines")
    ) +
    ggplot2::labs(
      x = "Standardized effect size (β)",
      y = NULL,
      color = NULL
    ) +
    ggplot2::geom_text(
      data = n_all_tbl,
      ggplot2::aes(x = Inf, y = Inf, label = label_full),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.2,
      color = col_framework["glm_full"],
      size = 3.4
    ) +
    ggplot2::geom_text(
      data = n_all_tbl,
      ggplot2::aes(x = Inf, y = Inf, label = label_match),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 2.5,
      color = col_framework["glm_matched"],
      size = 3.4
    ) +
    ggplot2::geom_text(
      data = n_all_tbl,
      ggplot2::aes(x = Inf, y = Inf, label = label_phylo),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 3.8,
      color = col_framework["phyloglm_matched"],
      size = 3.4
    )
  
  save_plot(
    p_eff,
    "9_01_effect_sizes_all_frameworks_by_tier_model",
    w = 20,
    h = 13
  )
  
  ## ge10-only MAIN TEXT FIGURE
  ## Built as separate panels, then combined with patchwork.
  ## This gives reliable outside-panel A-D labels.
  
  plot_df_ge10 <- plot_df %>%
    dplyr::filter(tier == "ge10")
  
  n_ge10_tbl <- meta_tbl2 %>%
    dplyr::filter(tier == "ge10") %>%
    dplyr::mutate(
      framework = factor(framework, levels = framework_levels),
      model = factor(model, levels = model_levels)
    ) %>%
    dplyr::arrange(model, framework) %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      label_full  = paste0("n full = ", n[framework == "glm_full"]),
      label_match = paste0("n matched = ", n[framework == "glm_matched"]),
      label_phylo = paste0("n phylo = ", n[framework == "phyloglm_matched"]),
      .groups = "drop"
    )
  
  make_ge10_panel <- function(mm) {
    
    dd <- plot_df_ge10 %>%
      dplyr::filter(model == mm)
    
    nn <- n_ge10_tbl %>%
      dplyr::filter(model == mm)
    
    ggplot2::ggplot(
      dd,
      ggplot2::aes(x = estimate, y = term_label, color = framework)
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = conf.low, xmax = conf.high),
        height = 0.18,
        linewidth = 0.8,
        position = ggplot2::position_dodge(width = 0.65, reverse = TRUE)
      ) +
      ggplot2::geom_point(
        size = 2.6,
        position = ggplot2::position_dodge(width = 0.65, reverse = TRUE)
      ) +
      ggplot2::geom_text(
        data = nn,
        ggplot2::aes(x = Inf, y = Inf, label = label_full),
        inherit.aes = FALSE,
        hjust = 1.05, vjust = 1.8,
        color = col_framework["glm_full"],
        size = 3.6
      ) +
      ggplot2::geom_text(
        data = nn,
        ggplot2::aes(x = Inf, y = Inf, label = label_match),
        inherit.aes = FALSE,
        hjust = 1.05, vjust = 3.8,
        color = col_framework["glm_matched"],
        size = 3.6
      ) +
      ggplot2::geom_text(
        data = nn,
        ggplot2::aes(x = Inf, y = Inf, label = label_phylo),
        inherit.aes = FALSE,
        hjust = 1.05, vjust = 5.8,
        color = col_framework["phyloglm_matched"],
        size = 3.6
      ) +
      ggplot2::scale_color_manual(
        values = col_framework,
        drop = FALSE,
        labels = framework_labels
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 11),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.title.x = ggplot2::element_text(size = 13),
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        panel.border = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.45),
        plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      ) +
      ggplot2::labs(
        title = dplyr::recode(as.character(mm), !!!model_labels),
        x = "Standardized effect size (β)",
        y = NULL,
        color = NULL
      )
  }
  
  ge10_panels <- lapply(model_levels, make_ge10_panel)
  
  p_eff_ge10 <- patchwork::wrap_plots(
    ge10_panels,
    ncol = 2,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 18),
        plot.tag.position = c(0, 1),
        legend.position = "bottom"
      )
    ) +
    patchwork::plot_layout(guides = "collect")
  
  save_plot(
    p_eff_ge10,
    "9_01b_effect_sizes_ge10_maintext",
    w = 13,
    h = 8.8
  )
  
  ## Subset effect
  
  subset_plot_df <- subset_effect_tbl %>%
    dplyr::filter(term %in% term_levels) %>%
    dplyr::mutate(
      term_label = dplyr::recode(term, !!!term_labels),
      term_label = factor(term_label, levels = rev(unname(term_labels[term_levels])))
    )
  
  p_subset <- ggplot2::ggplot(
    subset_plot_df,
    ggplot2::aes(x = delta_beta_subset, y = term_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = delta_beta_subset, y = term_label, yend = term_label),
      linewidth = 0.7,
      color = col_delta["subset_effect"]
    ) +
    ggplot2::geom_point(
      size = 3.0,
      color = col_delta["subset_effect"]
    ) +
    ggplot2::facet_grid(tier ~ model, scales = "free_y", space = "free_y") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      panel.border = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.45)
    ) +
    ggplot2::labs(
      x = expression(Delta*beta~"(matched GLM - full GLM)"),
      y = NULL
    )
  
  save_plot(
    p_subset,
    "9_02_delta_beta_subset_effect",
    w = 11,
    h = 11
  )
  
  ## Phylogeny effect
  
  phy_plot_df <- phylogeny_effect_tbl %>%
    dplyr::filter(term %in% term_levels) %>%
    dplyr::mutate(
      term_label = dplyr::recode(term, !!!term_labels),
      term_label = factor(term_label, levels = rev(unname(term_labels[term_levels])))
    )
  
  p_phy <- ggplot2::ggplot(
    phy_plot_df,
    ggplot2::aes(x = delta_beta_phylogeny, y = term_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = delta_beta_phylogeny, y = term_label, yend = term_label),
      linewidth = 0.7,
      color = col_delta["phylogeny_effect"]
    ) +
    ggplot2::geom_point(
      size = 3.0,
      color = col_delta["phylogeny_effect"]
    ) +
    ggplot2::facet_grid(tier ~ model, scales = "free_y", space = "free_y") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      panel.border = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.45)
    ) +
    ggplot2::labs(
      x = expression(Delta*beta~"(matched phyloGLM - matched GLM)"),
      y = NULL
    )
  
  save_plot(
    p_phy,
    "9_03_delta_beta_phylogeny_effect",
    w = 11,
    h = 11
  )
}

## 10. AIC comparison ----

message("[INFO] Step 10: AIC comparison")

if (nrow(meta_tbl2) > 0) {
  
  aic_tbl <- meta_tbl2 %>%
    dplyr::group_by(tier, framework) %>%
    dplyr::mutate(
      deltaAIC_within_framework = AIC - min(AIC, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(tier, framework, deltaAIC_within_framework)
  
  save_table(aic_tbl, "10_01_AIC_comparison_by_framework")
  
  p_aic <- ggplot2::ggplot(
    aic_tbl,
    ggplot2::aes(
      x = model,
      y = deltaAIC_within_framework,
      color = framework,
      group = framework
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 2,
      linetype = 2,
      linewidth = 0.45
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::facet_wrap(~ tier, ncol = 1) +
    ggplot2::scale_color_manual(
      values = col_framework,
      drop = FALSE,
      labels = c(
        glm_full = "Full GLM",
        glm_matched = "Matched GLM",
        phyloglm_matched = "Matched phyloGLM"
      )
    ) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "grey40",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(
        colour = "grey60",
        fill = NA,
        linewidth = 0.45
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "ΔAIC within tier and framework",
      color = NULL,
      title = "Model comparison across frameworks"
    )
  
  save_plot(
    p_aic,
    "10_02_AIC_comparison_by_framework",
    w = 9,
    h = 9
  )
}

## 11. Likelihood / pseudo-R2 summaries for matched GLM vs phyloGLM ----

message("[INFO] Step 11: Likelihood and alpha summaries for matched datasets")

partition_rows <- lapply(names(all_results), function(key) {
  res <- all_results[[key]]
  if (!isTRUE(res$ok)) return(NULL)
  if (!grepl("glm_matched|phyloglm_matched", key)) return(NULL)
  
  parts <- strsplit(key, "__")[[1]]
  tier <- parts[1]
  model <- parts[2]
  framework <- parts[3]
  
  if (framework == "glm_matched") {
    dsub <- res$data
    null_fit <- stats::glm(woody_bin ~ 1, data = dsub, family = stats::binomial())
    
    tibble(
      tier = tier,
      model = model,
      framework = framework,
      n = res$n,
      logLik_full = get_ll(res$fit),
      logLik_null = get_ll(null_fit),
      R2 = mcfadden_r2(get_ll(res$fit), get_ll(null_fit)),
      alpha = NA_real_
    )
    
  } else if (framework == "phyloglm_matched") {
    dsub <- res$data
    tr <- res$tree
    
    null_fit <- phylolm::phyloglm(
      woody_bin ~ 1,
      data = dsub,
      phy = tr,
      method = "logistic_MPLE",
      boot = 0
    )
    
    alpha_now <- NA_real_
    if (!is.null(res$fit$alpha)) alpha_now <- as.numeric(res$fit$alpha)
    
    tibble(
      tier = tier,
      model = model,
      framework = framework,
      n = res$n,
      logLik_full = get_ll(res$fit),
      logLik_null = get_ll(null_fit),
      R2 = mcfadden_r2(get_ll(res$fit), get_ll(null_fit)),
      alpha = alpha_now
    )
  } else {
    NULL
  }
})

partition_tbl <- bind_rows(partition_rows) %>%
  mutate(
    tier = factor(tier, levels = tier_levels),
    model = factor(model, levels = model_levels),
    framework = factor(framework, levels = c("glm_matched", "phyloglm_matched"))
  ) %>%
  arrange(tier, model, framework)

save_table(partition_tbl, "11_01_likelihood_pseudoR2_alpha_matched")

## 12. Supertribe heterogeneity (ge10 only; non-phylogenetic) ----

message("[INFO] Step 12: Supertribe heterogeneity (ge10 only; GLM only)")

if (!isTRUE(include_step10_supertribe)) {
  
  message("[INFO] Step 12 disabled by include_step10_supertribe=FALSE.")
  
} else if (!all(c("BIONOMIAL", "SUPERTRIBE", "TRIBE_FULL") %in% names(species_details))) {
  
  message("[WARN] species_details missing one of: BIONOMIAL, SUPERTRIBE, TRIBE_FULL. Skipping Step 12.")
  
} else {
  
  st_lu <- species_details %>%
    dplyr::transmute(
      key = taxon_key(BIONOMIAL),
      SUPERTRIBE = stringr::str_squish(as.character(SUPERTRIBE)),
      TRIBE_FULL = stringr::str_squish(as.character(TRIBE_FULL))
    ) %>%
    dplyr::mutate(
      SUPERTRIBE = dplyr::na_if(SUPERTRIBE, ""),
      SUPERTRIBE2 = dplyr::case_when(
        !is.na(SUPERTRIBE) ~ SUPERTRIBE,
        is.na(SUPERTRIBE) & !is.na(TRIBE_FULL) &
          stringr::str_detect(
            TRIBE_FULL,
            stringr::regex("^Aethionemeae\\b", ignore_case = TRUE)
          ) ~ "Aethionemeae",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(key, SUPERTRIBE = SUPERTRIBE2) %>%
    dplyr::filter(!is.na(key)) %>%
    dplyr::distinct()
  
  df10 <- tier_objs[["ge10"]]$full_model %>%
    dplyr::left_join(st_lu, by = "key") %>%
    dplyr::filter(!is.na(SUPERTRIBE))
  
  st_counts <- df10 %>%
    dplyr::count(SUPERTRIBE, name = "n_st") %>%
    dplyr::arrange(dplyr::desc(n_st))
  
  keep_st <- st_counts %>%
    dplyr::filter(n_st >= step10_min_n_per_supertribe) %>%
    dplyr::pull(SUPERTRIBE)
  
  df10 <- df10 %>%
    dplyr::filter(SUPERTRIBE %in% keep_st) %>%
    dplyr::mutate(SUPERTRIBE = factor(SUPERTRIBE))
  
  message("[INFO] ge10 supertribe heterogeneity dataset:")
  message("  n species: ", nrow(df10))
  message("  n supertribes kept: ", nlevels(df10$SUPERTRIBE))
  
  if (nrow(df10) < 300 || nlevels(df10$SUPERTRIBE) < 2) {
    
    message("[WARN] Too few rows or too few supertribes after filtering. Skipping Step 12.")
    
  } else {
    
    ## Include BIO10 for consistency with the integrated analyses.
    preds <- c(cols$MCWD_med, cols$BIO6_med, cols$BIO10_med, cols$elev_med)
    
    dsub <- df10 %>%
      dplyr::select(species, woody_bin, SUPERTRIBE, dplyr::all_of(preds)) %>%
      dplyr::mutate(
        woody_bin = as.integer(woody_bin),
        dplyr::across(
          dplyr::all_of(preds),
          ~ suppressWarnings(as.numeric(.x))
        )
      ) %>%
      dplyr::filter(stats::complete.cases(.)) %>%
      dplyr::mutate(SUPERTRIBE = factor(SUPERTRIBE))
    
    contrasts(dsub$SUPERTRIBE) <- stats::contr.sum(nlevels(dsub$SUPERTRIBE))
    
    for (cc in preds) {
      dsub[[cc]] <- zscore(dsub[[cc]])
    }
    
    if (nrow(dsub) < 300) {
      
      message("[WARN] Too few rows after complete-case filtering. Skipping Step 12.")
      
    } else {
      
      f_base <- stats::as.formula(paste0(
        "woody_bin ~ ",
        cols$MCWD_med, " + ",
        cols$BIO6_med, " + ",
        cols$BIO10_med, " + ",
        cols$elev_med
      ))
      
      f_int <- stats::as.formula(paste0(
        "woody_bin ~ ",
        cols$MCWD_med, " + ",
        cols$BIO6_med, " + ",
        cols$BIO10_med, " + ",
        cols$elev_med, " + ",
        "SUPERTRIBE + ",
        cols$MCWD_med, ":SUPERTRIBE + ",
        cols$BIO6_med, ":SUPERTRIBE + ",
        cols$BIO10_med, ":SUPERTRIBE + ",
        cols$elev_med, ":SUPERTRIBE"
      ))
      
      glm_base <- stats::glm(f_base, data = dsub, family = stats::binomial())
      glm_int  <- stats::glm(f_int,  data = dsub, family = stats::binomial())
      glm_null <- stats::glm(woody_bin ~ 1, data = dsub, family = stats::binomial())
      lrt_glm  <- stats::anova(glm_base, glm_int, test = "Chisq")
      
      R2_glm_base <- mcfadden_r2(
        as.numeric(stats::logLik(glm_base)),
        as.numeric(stats::logLik(glm_null))
      )
      
      R2_glm_int <- mcfadden_r2(
        as.numeric(stats::logLik(glm_int)),
        as.numeric(stats::logLik(glm_null))
      )
      
      comp_tbl <- tibble::tibble(
        tier = "ge10",
        n = nrow(dsub),
        n_supertribe = nlevels(dsub$SUPERTRIBE),
        min_n_supertribe = min(table(dsub$SUPERTRIBE)),
        glm_AIC_base = stats::AIC(glm_base),
        glm_AIC_int  = stats::AIC(glm_int),
        glm_deltaAIC = stats::AIC(glm_int) - stats::AIC(glm_base),
        glm_LRT_df = lrt_glm$Df[2],
        glm_LRT_LR = lrt_glm$Deviance[2],
        glm_LRT_p  = lrt_glm$`Pr(>Chi)`[2],
        glm_R2_base = R2_glm_base,
        glm_R2_int  = R2_glm_int,
        glm_deltaR2 = R2_glm_int - R2_glm_base
      )
      
      save_table(comp_tbl, "12_01_ge10_supertribe_interaction_model_comparison")
      
      coef_base_glm <- tidy_glm(glm_base) %>%
        dplyr::mutate(tier = "ge10", fit = "base")
      
      coef_int_glm <- tidy_glm(glm_int) %>%
        dplyr::mutate(tier = "ge10", fit = "interaction")
      
      coef_12 <- dplyr::bind_rows(coef_base_glm, coef_int_glm)
      
      save_table(coef_12, "12_02_ge10_supertribe_interaction_coefficients")
      
      fit_gain <- tibble::tibble(
        framework = "GLM",
        deltaAIC = comp_tbl$glm_deltaAIC,
        deltaR2  = comp_tbl$glm_deltaR2,
        LRT_p    = comp_tbl$glm_LRT_p
      ) %>%
        dplyr::mutate(
          neglog10_p = -log10(pmax(LRT_p, 1e-300))
        )
      
      p_gain_aic <- ggplot2::ggplot(
        fit_gain,
        ggplot2::aes(x = framework, y = deltaAIC)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45) +
        ggplot2::geom_col(width = 0.55) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::labs(
          x = NULL,
          y = expression(Delta*AIC~"(interaction - base)"),
          title = "Does allowing supertribe-specific climate slopes improve fit?"
        )
      
      save_plot(
        p_gain_aic,
        "12_03_ge10_supertribe_interaction_deltaAIC",
        w = 7,
        h = 5
      )
      
      p_gain_r2 <- ggplot2::ggplot(
        fit_gain,
        ggplot2::aes(x = framework, y = deltaR2)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45) +
        ggplot2::geom_col(width = 0.55) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::labs(
          x = NULL,
          y = expression(Delta*R[McFadden]^2~"(interaction - base)"),
          title = "Change in explanatory power from adding supertribe-specific climate slopes"
        )
      
      save_plot(
        p_gain_r2,
        "12_04_ge10_supertribe_interaction_deltaR2",
        w = 7,
        h = 5
      )
      
      p_gain_p <- ggplot2::ggplot(
        fit_gain,
        ggplot2::aes(x = framework, y = neglog10_p, group = 1)
      ) +
        ggplot2::geom_hline(
          yintercept = -log10(0.05),
          linetype = 2,
          linewidth = 0.45
        ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 3.0) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::labs(
          x = NULL,
          y = expression(-log[10](p)~"(LRT interaction vs base)"),
          title = "Evidence for supertribe-specific climate slopes"
        )
      
      save_plot(
        p_gain_p,
        "12_05_ge10_supertribe_interaction_LRT_neglog10p",
        w = 7,
        h = 5
      )
      
      ## 12.2b Reader-friendly supertribe deviation plot.
      ## Use finite-difference predictions, not coefficient-name parsing.
      
      slope_terms <- c("MCWD_med", "BIO6_med", "BIO10_med", "elev_med")
      
      term_label_map <- c(
        MCWD_med  = "Drought (MCWD)",
        BIO6_med  = "Frost (BIO6)",
        BIO10_med = "Heat (BIO10)",
        elev_med  = "Elevation"
      )
      
      st_levels <- levels(dsub$SUPERTRIBE)
      
      coef_base_terms <- stats::coef(summary(glm_base))
      coef_base_df <- as.data.frame(coef_base_terms)
      coef_base_df$term <- rownames(coef_base_df)
      rownames(coef_base_df) <- NULL
      
      coef_base_df <- coef_base_df %>%
        dplyr::rename(
          estimate = Estimate,
          se = `Std. Error`,
          z = `z value`,
          p = `Pr(>|z|)`
        ) %>%
        dplyr::mutate(
          conf.low = estimate - 1.96 * se,
          conf.high = estimate + 1.96 * se
        )
      
      overall_slopes <- coef_base_df %>%
        dplyr::filter(term %in% slope_terms) %>%
        dplyr::select(term, overall = estimate)
      
      st_dev_tbl <- dplyr::bind_rows(lapply(slope_terms, function(vv) {
        dplyr::bind_rows(lapply(st_levels, function(st) {
          
          nd0 <- tibble::tibble(
            SUPERTRIBE = factor(st, levels = st_levels),
            MCWD_med = 0,
            BIO6_med = 0,
            BIO10_med = 0,
            elev_med = 0
          )
          
          nd1 <- nd0
          nd1[[vv]] <- 1
          
          eta0 <- as.numeric(stats::predict(glm_int, newdata = nd0, type = "link"))
          eta1 <- as.numeric(stats::predict(glm_int, newdata = nd1, type = "link"))
          
          tibble::tibble(
            SUPERTRIBE = st,
            term = vv,
            slope_supertribe = eta1 - eta0
          )
        }))
      })) %>%
        dplyr::left_join(overall_slopes, by = "term") %>%
        dplyr::mutate(
          deviation = slope_supertribe - overall
        )
      
      st_levels_alpha <- sort(unique(as.character(st_dev_tbl$SUPERTRIBE)))
      
      st_dev_tbl <- st_dev_tbl %>%
        dplyr::mutate(
          SUPERTRIBE = factor(SUPERTRIBE, levels = rev(st_levels_alpha)),
          term_label = factor(
            dplyr::recode(term, !!!term_label_map),
            levels = c("Drought (MCWD)", "Frost (BIO6)", "Heat (BIO10)", "Elevation")
          )
        )
      
      st_n_tbl <- dsub %>%
        dplyr::count(SUPERTRIBE, name = "n_supertribe") %>%
        dplyr::mutate(SUPERTRIBE_chr = as.character(SUPERTRIBE))
      
      st_lab_map <- st_n_tbl %>%
        dplyr::mutate(
          SUPERTRIBE_label = paste0(SUPERTRIBE_chr, " (n = ", n_supertribe, ")")
        ) %>%
        dplyr::select(SUPERTRIBE_chr, SUPERTRIBE_label)
      
      st_dev_tbl <- st_dev_tbl %>%
        dplyr::mutate(SUPERTRIBE_chr = as.character(SUPERTRIBE)) %>%
        dplyr::left_join(st_lab_map, by = "SUPERTRIBE_chr") %>%
        dplyr::mutate(
          SUPERTRIBE_label = factor(
            SUPERTRIBE_label,
            levels = rev(st_lab_map$SUPERTRIBE_label[order(st_lab_map$SUPERTRIBE_chr)])
          )
        )
      
      save_table(st_dev_tbl, "12_02b_ge10_supertribe_derived_slope_deviations")
      
      p_st_dev <- ggplot2::ggplot(
        st_dev_tbl,
        ggplot2::aes(x = deviation, y = SUPERTRIBE_label, color = SUPERTRIBE)
      ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = 2,
          linewidth = 0.45
        ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0,
            xend = deviation,
            y = SUPERTRIBE_label,
            yend = SUPERTRIBE_label
          ),
          linewidth = 0.8
        ) +
        ggplot2::geom_point(size = 3.2) +
        ggplot2::scale_color_manual(values = st_cols, drop = FALSE) +
        ggplot2::facet_wrap(
          ~ term_label,
          ncol = 2,
          scales = "free_x"
        ) +
        ggplot2::theme_classic(base_size = 15) +
        ggplot2::theme(
          legend.position = "none",
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "grey40",
            linewidth = 0.5
          ),
          strip.text = ggplot2::element_text(face = "bold"),
          panel.border = ggplot2::element_rect(
            colour = "grey60",
            fill = NA,
            linewidth = 0.45
          ),
          panel.spacing = grid::unit(1.2, "lines")
        ) +
        ggplot2::labs(
          x = "Deviation from overall slope",
          y = NULL,
          title = "Supertribe deviations from the overall ge10 climate–woodiness relationship",
          subtitle = "Values right of zero indicate stronger-than-overall slopes; values left of zero indicate weaker-than-overall slopes"
        )
      
      save_plot(
        p_st_dev,
        "12_06_ge10_supertribe_derived_slope_deviations",
        w = 9,
        h = 9
      )
    }
  }
}

## 13. Save objects ----

message("[INFO] Step 13: Saving objects")

saveRDS(
  list(
    tree_b = tree_b,
    tip_key = tip_key,
    species_lu = species_lu,
    tier_objs = tier_objs,
    model_set = model_set,
    model_levels = model_levels,
    framework_levels = framework_levels,
    fit_index = fit_index,
    coef_tbl = coef_tbl2,
    meta_tbl = meta_tbl2,
    subset_effect_tbl = subset_effect_tbl,
    phylogeny_effect_tbl = phylogeny_effect_tbl,
    partition_tbl = partition_tbl,
    all_results = all_results
  ),
  file.path(out_dir, "OBJECT_13_01_integrated_model_objects.rds")
)

message("[DONE] 2g.WoodinessNiche_IntegratedModels.R finished successfully.")
