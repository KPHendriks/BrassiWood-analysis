#!/usr/bin/env Rscript

## 1c2.Draw_species_tree.R
## Extended version:
##   - original nuclear workflow retained
##
## Supported plot types:
##   Nuclear tree:
##      - calibrated MCMCtree representative tree
##      - ASTRAL tree
##      - CF tree with ASTRAL-style quartet support

##
## Output:
##   - nuclear annotated tree (optional)

suppressPackageStartupMessages({
  library(ape)
  library(treeio)
  library(tidytree)
  library(phytools)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(readr)
  library(purrr)
  library(tibble)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
left_join <- dplyr::left_join
arrange <- dplyr::arrange
summarise <- dplyr::summarise
transmute <- dplyr::transmute

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}


# ---------------- USER SETTINGS -----------------------

interactive_setup <- FALSE

## Project / root
ROOT <- "."
RFIN <- file.path(ROOT, "WP1_BrassiToL/results_final")

species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
specimen_file <- "WP2_BrassiNiche/data/brassiwood_specimens_publication.csv"

## ---------------- WHAT TO DRAW ----------------
draw_nuclear_tree  <- TRUE

## ---------------- NUCLEAR SETTINGS ----------------

project_name <- "2026-03-09_BrassiWood_v2"

## Calibration source for the nuclear tree
calibration_source <- "lsd2"   # "mcmctree" or "lsd2"

## MCMCtree settings
mct_run_tag   <- "mcmctree_5loci_sortadate"
mct_profile   <- "02_screen"
mct_replicate <- "rep_02"
calibrated_tree_manual <- NULL

## LSD2 settings
lsd2_n_loci <- 200
lsd2_run_tag <- paste0("lsd2_top", lsd2_n_loci, "_fixedtopo")
lsd2_tree_manual <- NULL

## LSD2 main output from STEP 10.6/10.7
lsd2_tree <- file.path(
  RFIN, project_name,
  "8_results_species_tree_calibrated",
  lsd2_run_tag,
  "05_lsd2_dating",
  paste0(project_name, "_lsd2_top", lsd2_n_loci, "_dated.date.nexus")
)

astral_tree <- file.path(
  RFIN, project_name, "6_results_astral_species_tree",
  paste0(project_name, "_BrassiToL_coalescent.tree")
)

cf_tree <- file.path(
  RFIN, project_name, "7_results_iqtree_species_tree",
  paste0(project_name, "_BrassiToL_concordance_factors.cf.tree")
)

primary_tree_type <- "calibrated"   # one of: "calibrated", "astral", "cf"

root_method <- "keep"               # "keep", "outgroup", "midpoint"
outgroup_tips <- c("PAFTOL_014331")

make_primary_ultrametric <- FALSE
ultrametric_method <- "extend"      # "chronos" or "extend"

## Select samples to represent each of the main clades when drawing a small backbone version
backbone_representatives <- c(
  "Cleomaceae outgroup" = "1768",
  "Aethionemeae"        = "S0735",
  "Arabodae (IV)"       = "S1120",
  "Brassicodae (II)"    = "S1376",
  "Camelinodae (I)"     = "S0554",
  "Heliophilodae (V)"   = "S0892",
  "Hesperodae (III)"    = "S0603"
)

## LSD2 was run with the fixed root age defined directly in Ma (-a -90).
## Use factor 1. Use 100 only for trees scaled as e.g. 0.90 = 90 Ma.
calibrated_time_scaling_factor <- 1

q1_breaks  <- c(0.50, 0.75)
scf_breaks <- c(33, 66)
lpp_high_threshold <- 0.95

nuclear_out_dir <- file.path(
  RFIN, project_name,
  "9_results_species_tree_for_publication"
)

nuclear_out_pdf <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_species_tree.pdf")
)

nuclear_out_png <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_species_tree.png")
)

draw_tribe_cladelabs <- TRUE
tribe_cladelab_offset <- 4.5
tribe_cladelab_barsize <- 0.55
tribe_cladelab_textsize <- 1.6
tribe_cladelab_min_tips <- 2

tribe_component_min_tips <- 2
tribe_component_min_fraction <- 0.75

tribe_crown_icon_size <- 2.2
tribe_crown_icon_alpha <- 0.25

supertribe_bar_offset <- 10
supertribe_bar_barsize <- 1.1
supertribe_bar_textsize <- 2.0


## ---------------- SUPPORT-OVER-TIME SETTINGS ----------------

draw_nuclear_support_time <- TRUE

support_time_out_pdf <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_support_over_time.pdf")
)

support_time_out_png <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_support_over_time.png")
)

support_time_bin_width <- 2     # Ma
support_time_min_age <- 0
support_time_max_age <- 95

support_time_highlight_tribes <- TRUE
support_time_highlight_supertribes <- TRUE
support_time_highlight_family <- TRUE

support_time_family_name <- "Brassicaceae"
support_time_aethionemeae_name <- "Aethionemeae"



## ---------------- SHARED DISPLAY SETTINGS ----------------

tip_label_size <- 0.5

node_icon_size_square <- 0.4
node_icon_size_circle <- 0.2
node_icon_stroke <- 0.1


gene_bar_offset <- -0.35
gene_bar_pwidth <- 0.16


geo_heat_offset <- 20
geo_col_spacing <- 0.35
geo_tile_height <- 0.8
geo_label_size <- 0.65

time_axis_limits <- c(-95, 0)

plot_width  <- 18
plot_height <- 49
plot_dpi    <- 400

## Optional taxonomy check settings
target_rank <- "family"
target_clade <- "Brassicaceae"
taxonomy_check_threshold <- 0.49

## Tip-level metadata column used for "genes retained"
genes_retained_col <- "loci_remaining"

st_cols <- c(
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8",
  "Aethionemeae"      = "#ffed57"
)

support_time_node_table_out <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_family_supertribe_tribe_node_ages_support.csv")
)

# ---------------- HELPER FUNCTIONS --------------------

ask_yes_no <- function(prompt, default = TRUE) {
  default_txt <- if (default) "Y/n" else "y/N"
  ans <- readline(paste0(prompt, " [", default_txt, "]: "))
  ans <- trimws(tolower(ans))
  if (ans == "") return(default)
  ans %in% c("y", "yes")
}

ask_choice <- function(prompt, choices, default = NULL) {
  show <- paste0(choices, collapse = "/")
  if (!is.null(default)) {
    ans <- readline(paste0(prompt, " [", show, "] (default=", default, "): "))
    ans <- trimws(ans)
    if (ans == "") ans <- default
  } else {
    ans <- readline(paste0(prompt, " [", show, "]: "))
    ans <- trimws(ans)
  }
  if (!ans %in% choices) {
    stop("[ERR] Invalid choice for ", prompt, ": ", ans)
  }
  ans
}

ask_character <- function(prompt, default = NULL) {
  if (is.null(default)) {
    ans <- readline(paste0(prompt, ": "))
  } else {
    ans <- readline(paste0(prompt, " (default=", default, "): "))
    if (trimws(ans) == "") ans <- default
  }
  trimws(ans)
}

resolve_calibrated_tree_path <- function(RFIN, project_name, mct_run_tag, mct_profile,
                                         mct_replicate, calibrated_tree_manual = NULL) {
  candidates <- c(
    file.path(
      RFIN, project_name, "8_results_species_tree_calibrated", mct_run_tag,
      "06_usedata2_profiles", mct_profile,
      paste0(project_name, "_", sub("^\\d+_", "", mct_profile), "_representative_figtree.tre")
    ),
    file.path(
      RFIN, project_name, "8_results_species_tree_calibrated", mct_run_tag,
      "06_usedata2_profiles", mct_profile, "replicates", mct_replicate, "FigTree.tre"
    ),
    calibrated_tree_manual
  )
  
  candidates <- candidates[!is.na(candidates) & !vapply(candidates, is.null, logical(1)) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) == 0) return(NULL)
  existing[1]
}

resolve_lsd2_tree_path <- function(RFIN, project_name, lsd2_run_tag,
                                   lsd2_n_loci,
                                   lsd2_tree_manual = NULL) {
  candidates <- c(
    lsd2_tree_manual,
    file.path(
      RFIN, project_name,
      "8_results_species_tree_calibrated",
      lsd2_run_tag,
      "05_lsd2_dating",
      paste0(project_name, "_lsd2_top", lsd2_n_loci, "_dated.date.nexus")
    ),
    file.path(
      RFIN, project_name,
      "8_results_species_tree_calibrated",
      lsd2_run_tag,
      "05_lsd2_dating",
      paste0(project_name, "_lsd2_top", lsd2_n_loci, "_dated.nwk")
    )
  )
  
  candidates <- candidates[
    !is.na(candidates) &
      !vapply(candidates, is.null, logical(1)) &
      nzchar(candidates)
  ]
  
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

read_any_tree <- function(path, type = c("newick", "astral", "iqtree", "mcmctree", "lsd2"),
                          force_ultrametric_mcmctree = FALSE) {
  type <- match.arg(type)
  
  if (is.null(path) || !file.exists(path)) {
    stop("[ERR] Tree file not found: ", path %||% "NULL")
  }
  
  td <- switch(
    type,
    newick = tryCatch(
      treeio::read.newick(path),
      error = function(e) treeio::as.treedata(ape::read.tree(path))
    ),
    astral = treeio::read.astral(path),
    iqtree = treeio::read.iqtree(path),
    mcmctree = treeio::read.mcmctree(path, force.ultrametric = force_ultrametric_mcmctree),
    lsd2 = tryCatch(
      treeio::read.beast(path),
      error = function(e1) {
        message("[WARN] treeio::read.beast() failed for LSD2 Nexus; trying ape::read.nexus().")
        tryCatch(
          treeio::as.treedata(ape::read.nexus(path)),
          error = function(e2) {
            message("[WARN] ape::read.nexus() failed; trying ape::read.tree().")
            treeio::as.treedata(ape::read.tree(path))
          }
        )
      }
    )
  )
  
  phy <- tryCatch(
    as.phylo(td),
    error = function(e) {
      if (inherits(td, "phylo")) {
        td
      } else if (!is.null(td@phylo)) {
        td@phylo
      } else {
        stop("[ERR] Could not convert tree object to phylo for file: ", path)
      }
    }
  )
  
  list(td = td, phy = phy, path = path, type = type)
}

root_tree_safely <- function(phy, root_method = c("keep", "outgroup", "midpoint"),
                             outgroup_tips = NULL) {
  root_method <- match.arg(root_method)
  
  if (root_method == "keep") return(phy)
  
  if (root_method == "midpoint") {
    return(phytools::midpoint.root(phy))
  }
  
  if (root_method == "outgroup") {
    if (is.null(outgroup_tips) || length(outgroup_tips) == 0) {
      stop("[ERR] root_method='outgroup' but no outgroup_tips supplied.")
    }
    outgroup_tips <- intersect(outgroup_tips, phy$tip.label)
    if (length(outgroup_tips) == 0) {
      stop("[ERR] None of the specified outgroup tips are present in the tree.")
    }
    return(ape::root(phy, outgroup = outgroup_tips, resolve.root = TRUE))
  }
  
  phy
}

make_ultrametric_if_requested <- function(phy,
                                          do_it = FALSE,
                                          method = c("chronos", "extend"),
                                          lambda = 1,
                                          model = "correlated",
                                          quiet = FALSE,
                                          control = NULL) {
  method <- match.arg(method)
  
  if (!do_it) return(phy)
  if (is.ultrametric(phy)) return(phy)
  
  if (!ape::is.rooted(phy)) {
    stop("[ERR] Cannot ultrametricize with chronos(): tree is not rooted.")
  }
  
  if (method == "chronos") {
    if (is.null(control)) {
      control <- ape::chronos.control(
        iter.max = 1e5,
        eval.max = 1e5,
        dual.iter.max = 40,
        epsilon = 1e-8
      )
    }
    
    message("[INFO] Making tree ultrametric using ape::chronos() with model='",
            model, "', lambda=", lambda)
    
    out <- ape::chronos(
      phy,
      lambda = lambda,
      model = model,
      quiet = quiet,
      control = control
    )
    
    ## chronos returns class c("chronos", "phylo");
    ## strip extra class so treeio::as.treedata() accepts it cleanly
    class(out) <- "phylo"
    
    return(out)
  }
  
  message("[INFO] Making tree ultrametric using phytools::force.ultrametric(method='extend')")
  phytools::force.ultrametric(phy, method = "extend")
}

rescale_phylo_to_target_depth <- function(phy, target_depth = 40) {
  if (is.null(phy$edge.length)) return(phy)
  
  current_depth <- max(node.depth.edgelength(phy))
  if (is.na(current_depth) || current_depth <= 0) return(phy)
  
  scale_factor <- target_depth / current_depth
  phy$edge.length <- phy$edge.length * scale_factor
  phy
}

ensure_label_column <- function(td_merged) {
  td_tbl <- tidytree::as_tibble(td_merged) %>% tibble::as_tibble()
  
  ## If label is already present and usable, keep it
  if ("label" %in% names(td_tbl)) {
    td_merged@data <- td_tbl
    return(td_merged)
  }
  
  ## Otherwise reconstruct from phylo
  phy <- td_merged@phylo
  n_tip <- ape::Ntip(phy)
  n_node <- phy$Nnode
  
  label_vec <- c(phy$tip.label, rep(NA_character_, n_node))
  
  if (!"node" %in% names(td_tbl)) {
    stop("[ERR] td_merged@data has no 'node' column, so labels cannot be reconstructed.")
  }
  
  td_tbl <- td_tbl %>%
    mutate(
      label = label_vec[match(node, seq_len(n_tip + n_node))]
    )
  
  td_merged@data <- td_tbl
  td_merged
}

plot_support_over_time_basic <- function(support_time_tbl) {
  
  plot_one <- function(df, med, lo, hi, title, ylab) {
    ggplot(df, aes(x = -age_mid)) +
      geom_ribbon(
        aes(
          ymin = .data[[lo]],
          ymax = .data[[hi]]
        ),
        alpha = 0.25
      ) +
      geom_line(aes(y = .data[[med]]), linewidth = 0.4) +
      geom_point(aes(y = .data[[med]]), size = 0.8) +
      scale_x_continuous(
        limits = c(-support_time_max_age, 0),
        breaks = seq(-90, 0, by = 10),
        labels = abs
      ) +
      labs(
        title = title,
        x = "Age from present (Ma)",
        y = ylab
      ) +
      theme_bw(base_size = 9)
  }
  
  p1 <- plot_one(support_time_tbl, "q1_med",  "q1_lo",  "q1_hi",  "ASTRAL Q1 over time", "Q1")
  p2 <- plot_one(support_time_tbl, "pp1_med", "pp1_lo", "pp1_hi", "ASTRAL LPP over time", "LPP")
  p3 <- plot_one(support_time_tbl, "sCF_med", "sCF_lo", "sCF_hi", "IQ-TREE sCF over time", "sCF")
  p4 <- plot_one(support_time_tbl, "gCF_med", "gCF_lo", "gCF_hi", "IQ-TREE gCF over time", "gCF")
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("[ERR] Package 'patchwork' is needed for support-over-time plots.")
  }
  
  p1 / p2 / p3 / p4
}

make_epoch_bar <- function(xlim = c(-95, 0)) {
  epochs <- tibble::tribble(
    ~name,          ~start, ~end,  ~fill,
    "Cretaceous",   -95,   -66,   "#9ACD32",
    "Paleocene",    -66,   -56,   "#FDB462",
    "Eocene",       -56,   -33.9, "#B3DE69",
    "Oligocene",    -33.9, -23.0, "#80B1D3",
    "Miocene",      -23.0, -5.33, "#FB8072",
    "Pliocene",     -5.33, -2.58, "#BEBADA",
    "Pleistocene",  -2.58, -0.0117, "#FFFFB3",
    "Holocene",     -0.0117, 0, "#D9D9D9"
  )
  
  ggplot(epochs) +
    geom_rect(
      aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = name),
      colour = "white",
      linewidth = 0.2
    ) +
    geom_text(
      aes(x = (start + end) / 2, y = 0.5, label = name),
      size = 2.5
    ) +
    scale_fill_manual(values = setNames(epochs$fill, epochs$name), guide = "none") +
    scale_x_continuous(limits = xlim, breaks = seq(-90, 0, 10), labels = abs) +
    theme_void() +
    theme(
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    )
}

make_rank_mrca_table <- function(td_merged,
                                 specimen_details,
                                 ranks = c("FAMILY", "SUPERTRIBE", "TRIBE"),
                                 family_filter = "Brassicaceae") {
  
  tre <- td_merged@phylo
  Ntip <- ape::Ntip(tre)
  
  node_depth <- ape::node.depth.edgelength(tre)
  tree_height <- max(node_depth, na.rm = TRUE)
  node_age_ma <- tree_height - node_depth
  
  node_dat <- tidytree::as_tibble(td_merged) %>%
    tibble::as_tibble()
  
  get_num_col <- function(tbl, col) {
    if (col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[col]])) else rep(NA_real_, nrow(tbl))
  }
  
  node_dat <- node_dat %>%
    mutate(
      age_Ma = node_age_ma[node],
      q1_num  = get_num_col(., "q1"),
      pp1_num = get_num_col(., "pp1"),
      gCF_num = get_num_col(., "gCF"),
      sCF_num = get_num_col(., "sCF")
    )
  
  tip_info <- specimen_details %>%
    mutate(
      SAMPLE = as.character(SAMPLE),
      FAMILY = as.character(FAMILY),
      SUPERTRIBE = as.character(SUPERTRIBE),
      TRIBE = as.character(TRIBE)
    ) %>%
    select(SAMPLE, FAMILY, SUPERTRIBE, TRIBE)
  
  tip_info <- tip_info[match(tre$tip.label, tip_info$SAMPLE), ]
  
  out <- list()
  
  for (rank_col in ranks) {
    
    if (rank_col == "FAMILY") {
      
      vals <- family_filter
      
    } else {
      
      vals <- tip_info %>%
        filter(FAMILY == family_filter) %>%
        pull(all_of(rank_col)) %>%
        unique() %>%
        na.omit() %>%
        sort()
    }
    
    for (val in vals) {
      if (rank_col != "FAMILY") {
        idx <- which(tip_info[[rank_col]] == val & tip_info$FAMILY == family_filter)
      } else {
        idx <- which(tip_info[[rank_col]] == val)
      }
      
      if (length(idx) == 0) next
      
      crown_node <- if (length(idx) == 1) {
        idx[1]
      } else {
        tryCatch(ape::getMRCA(tre, idx), error = function(e) NA_integer_)
      }
      
      if (is.na(crown_node)) next
      
      parent <- tre$edge[tre$edge[, 2] == crown_node, 1]
      stem_node <- if (length(parent) == 0) NA_integer_ else parent[1]
      
      crown_row <- node_dat %>% filter(node == crown_node)
      stem_row  <- node_dat %>% filter(node == stem_node)
      
      out[[paste(rank_col, val, sep = "__")]] <- tibble(
        rank = rank_col,
        clade = val,
        n_tips = length(idx),
        crown_node = crown_node,
        crown_age_Ma = crown_row$age_Ma %||% NA_real_,
        crown_q1 = crown_row$q1_num %||% NA_real_,
        crown_pp1 = crown_row$pp1_num %||% NA_real_,
        crown_gCF = crown_row$gCF_num %||% NA_real_,
        crown_sCF = crown_row$sCF_num %||% NA_real_,
        stem_node = stem_node,
        stem_age_Ma = stem_row$age_Ma %||% NA_real_,
        stem_q1 = stem_row$q1_num %||% NA_real_,
        stem_pp1 = stem_row$pp1_num %||% NA_real_,
        stem_gCF = stem_row$gCF_num %||% NA_real_,
        stem_sCF = stem_row$sCF_num %||% NA_real_
      )
    }
  }
  
  bind_rows(out)
}

plot_support_over_time_annotated <- function(support_time_tbl,
                                             highlight_tbl,
                                             xlim = c(-95, 0)) {
  
  highlight_tbl <- highlight_tbl %>%
    mutate(
      clade_clean = stringr::str_squish(as.character(clade)),
      is_unplaced = stringr::str_to_lower(clade_clean) == "unplaced",
      x = -crown_age_Ma,
      
      highlight_type = case_when(
        rank == "FAMILY" & clade_clean == "Brassicaceae" ~ "main",
        clade_clean %in% names(st_cols) ~ "main",
        rank == "SUPERTRIBE" & !is_unplaced ~ "main",
        rank == "TRIBE" & clade_clean == "Aethionemeae" ~ "main",
        rank == "TRIBE" & !is_unplaced ~ "tribe",
        TRUE ~ "skip"
      ),
      
      colour_group = case_when(
        rank == "FAMILY" & clade_clean == "Brassicaceae" ~ "Family crown",
        clade_clean %in% names(st_cols) ~ clade_clean,
        rank == "SUPERTRIBE" & !is_unplaced ~ clade_clean,
        rank == "TRIBE" & clade_clean == "Aethionemeae" ~ "Aethionemeae",
        TRUE ~ NA_character_
      ),
      
      point_size = case_when(
        rank == "FAMILY" ~ 3.2,
        rank == "SUPERTRIBE" ~ 2.5,
        clade_clean == "Aethionemeae" ~ 2.8,
        rank == "TRIBE" ~ 1.0,
        TRUE ~ 1.2
      )
    ) %>%
    filter(!is.na(x), highlight_type != "skip")
  
  main_tbl <- highlight_tbl %>%
    filter(highlight_type == "main", !is.na(colour_group))
  
  tribe_tbl <- highlight_tbl %>%
    filter(highlight_type == "tribe")
  
  legend_levels <- c(
    "Family crown",
    names(st_cols)
  )
  
  legend_cols <- c(
    "Family crown" = "black",
    st_cols
  )
  
  plot_one <- function(df, med, lo, hi, metric_col, title, ylab, show_legend = FALSE) {
    
    ggplot(df, aes(x = -age_mid)) +
      geom_ribbon(
        aes(ymin = .data[[lo]], ymax = .data[[hi]]),
        alpha = 0.20
      ) +
      geom_line(aes(y = .data[[med]]), linewidth = 0.4) +
      geom_point(aes(y = .data[[med]]), size = 0.7) +
      
      ## Tribe crown nodes: thin grey vertical lines + grey dots
      geom_vline(
        data = tribe_tbl,
        aes(xintercept = x),
        colour = "grey75",
        linetype = "dashed",
        linewidth = 0.18,
        alpha = 0.45,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = tribe_tbl,
        aes(x = x, y = .data[[metric_col]]),
        colour = "grey55",
        size = 0.9,
        alpha = 0.75,
        inherit.aes = FALSE
      ) +
      
      ## Main clades: family + supertribes + Aethionemeae
      geom_vline(
        data = main_tbl,
        aes(xintercept = x, colour = colour_group),
        linetype = "dashed",
        linewidth = 0.35,
        alpha = 0.75,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = main_tbl,
        aes(x = x, y = .data[[metric_col]], colour = colour_group, size = point_size),
        inherit.aes = FALSE,
        alpha = 0.98
      ) +
      scale_size_identity() +
      scale_colour_manual(
        name = "Main clade",
        values = legend_cols,
        breaks = legend_levels,
        limits = legend_levels,
        drop = FALSE,
        guide = if (show_legend) {
          guide_legend(override.aes = list(size = 3, linewidth = 0.6))
        } else {
          "none"
        }
      ) +
      scale_x_continuous(
        limits = xlim,
        breaks = seq(-90, 0, by = 10),
        labels = abs
      ) +
      labs(title = title, x = NULL, y = ylab) +
      theme_bw(base_size = 9) +
      theme(
        legend.position = if (show_legend) "right" else "none",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 10)
      )
  }
  
  p1 <- plot_one(
    support_time_tbl,
    "q1_med", "q1_lo", "q1_hi", "crown_q1",
    "ASTRAL quartet support over time", "Q1",
    show_legend = FALSE
  )
  
  p2 <- plot_one(
    support_time_tbl,
    "pp1_med", "pp1_lo", "pp1_hi", "crown_pp1",
    "ASTRAL local posterior probability", "LPP",
    show_legend = TRUE
  )
  
  p3 <- plot_one(
    support_time_tbl,
    "gCF_med", "gCF_lo", "gCF_hi", "crown_gCF",
    "Gene concordance factor", "gCF",
    show_legend = FALSE
  )
  
  p4 <- plot_one(
    support_time_tbl,
    "sCF_med", "sCF_lo", "sCF_hi", "crown_sCF",
    "Site concordance factor", "sCF",
    show_legend = FALSE
  )
  
  p_nodes <- ggplot(support_time_tbl, aes(x = -age_mid, y = n_nodes)) +
    geom_col(width = support_time_bin_width * 0.85) +
    scale_y_log10() +
    scale_x_continuous(limits = xlim, breaks = seq(-90, 0, 10), labels = abs) +
    labs(x = "Age from present (Ma)", y = "Nodes per bin (log10)", title = "Node density") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  p_epoch <- make_epoch_bar(xlim = xlim)
  
  p1 <- p1 + labs(tag = "A")
  p2 <- p2 + labs(tag = "B")
  p3 <- p3 + labs(tag = "C")
  p4 <- p4 + labs(tag = "D")
  p_nodes <- p_nodes + labs(tag = "E")
  
  panel_tag_theme <- theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.01, 0.98)
  )
  
  p_main <- patchwork::wrap_plots(
    p1 + panel_tag_theme,
    p2 + panel_tag_theme,
    p3 + panel_tag_theme,
    p4 + panel_tag_theme,
    p_nodes + panel_tag_theme,
    ncol = 1,
    heights = c(1, 1, 1, 1, 0.55)
  )
  
  patchwork::wrap_plots(
    p_main,
    p_epoch,
    ncol = 1,
    heights = c(5.55, 0.18)
  )
  
  patchwork::wrap_plots(
    p_main,
    p_epoch,
    ncol = 1,
    heights = c(5.55, 0.18)
  )
}

get_tribe_mrca <- function(tribe_name,
                           tree,
                           tip_info,
                           node_age,
                           desc_list,
                           manual_crown_vec = NULL,
                           max_other_in_raw = 10L,
                           min_core_fraction = 0.5) {
  
  Ntip <- ape::Ntip(tree)
  
  all_tribe_idx <- which(
    tip_info$TRIBE == tribe_name &
      tip_info$FAMILY == "Brassicaceae"
  )
  
  n_tribe_tips_total <- length(all_tribe_idx)
  
  empty_return <- function() {
    tibble(
      TRIBE = tribe_name,
      n_tribe_tips_total = n_tribe_tips_total,
      raw_crown_node = NA_integer_,
      raw_crown_age_Ma = NA_real_,
      raw_n_crown_tips_total = NA_integer_,
      raw_n_other_tips = NA_integer_,
      raw_status_phyly = NA_character_,
      main_node = NA_integer_,
      main_n_tribe_tips = NA_integer_,
      main_n_total_tips = NA_integer_,
      main_fraction = NA_real_
    )
  }
  
  if (n_tribe_tips_total == 0L) return(empty_return())
  
  raw_crown_node <- if (n_tribe_tips_total == 1L) {
    all_tribe_idx[1]
  } else {
    tryCatch(ape::getMRCA(tree, all_tribe_idx), error = function(e) NA_integer_)
  }
  
  if (is.na(raw_crown_node)) return(empty_return())
  
  raw_desc_all <- desc_list[[raw_crown_node]]
  raw_desc_tips <- raw_desc_all[raw_desc_all <= Ntip]
  
  raw_n_crown_tips_total <- length(raw_desc_tips)
  raw_n_other_tips <- sum(
    tip_info$TRIBE[raw_desc_tips] != tribe_name |
      tip_info$FAMILY[raw_desc_tips] != "Brassicaceae",
    na.rm = TRUE
  )
  
  tribe_tips_outside_raw <- setdiff(all_tribe_idx, raw_desc_tips)
  
  raw_status_phyly <- dplyr::case_when(
    length(tribe_tips_outside_raw) > 0L || raw_n_other_tips > max_other_in_raw ~ "polyphyletic",
    raw_n_other_tips > 0L ~ "paraphyletic",
    TRUE ~ "monophyletic"
  )
  
  if (raw_n_other_tips <= max_other_in_raw) {
    main_node <- raw_crown_node
  } else {
    internal_nodes <- (Ntip + 1L):(Ntip + tree$Nnode)
    
    best_node <- NA_integer_
    best_n_focal <- 0L
    best_fraction <- -Inf
    best_n_total <- NA_integer_
    
    for (nd in internal_nodes) {
      d <- desc_list[[nd]]
      tips_nd <- d[d <= Ntip]
      if (length(tips_nd) == 0L) next
      
      n_focal <- sum(
        tip_info$TRIBE[tips_nd] == tribe_name &
          tip_info$FAMILY[tips_nd] == "Brassicaceae",
        na.rm = TRUE
      )
      if (n_focal == 0L) next
      
      n_total <- length(tips_nd)
      n_other <- n_total - n_focal
      frac <- n_focal / n_total
      
      if (
        n_other <= max_other_in_raw &&
        frac >= min_core_fraction &&
        (n_focal > best_n_focal ||
         (n_focal == best_n_focal && frac > best_fraction))
      ) {
        best_node <- nd
        best_n_focal <- n_focal
        best_fraction <- frac
        best_n_total <- n_total
      }
    }
    
    main_node <- best_node
  }
  
  if (is.na(main_node)) return(empty_return())
  
  main_desc_all <- desc_list[[main_node]]
  main_tip_idx <- main_desc_all[main_desc_all <= Ntip]
  
  main_n_tribe_tips <- sum(
    tip_info$TRIBE[main_tip_idx] == tribe_name &
      tip_info$FAMILY[main_tip_idx] == "Brassicaceae",
    na.rm = TRUE
  )
  
  main_n_total_tips <- length(main_tip_idx)
  main_fraction <- main_n_tribe_tips / main_n_total_tips
  
  tibble(
    TRIBE = tribe_name,
    n_tribe_tips_total = n_tribe_tips_total,
    raw_crown_node = raw_crown_node,
    raw_crown_age_Ma = node_age[raw_crown_node],
    raw_n_crown_tips_total = raw_n_crown_tips_total,
    raw_n_other_tips = raw_n_other_tips,
    raw_status_phyly = raw_status_phyly,
    main_node = main_node,
    main_n_tribe_tips = main_n_tribe_tips,
    main_n_total_tips = main_n_total_tips,
    main_fraction = main_fraction
  )
}

add_tribe_cladelabs <- function(p,
                                td_merged,
                                specimen_details,
                                min_tips = 20,
                                offset = 2.5,
                                barsize = 0.35,
                                textsize = 1.6,
                                min_fraction = 0.75,
                                component_min_tips = 3,
                                component_min_fraction = 0.75) {
  
  tre <- td_merged@phylo
  Ntip <- ape::Ntip(tre)
  Nnode_total <- Ntip + tre$Nnode
  
  desc_list <- lapply(seq_len(Nnode_total), function(nd) {
    tryCatch(phytools::getDescendants(tre, nd), error = function(e) integer(0))
  })
  
  tip_info <- specimen_details %>%
    transmute(
      SAMPLE = as.character(SAMPLE),
      TRIBE = stringr::str_squish(as.character(TRIBE)),
      SUPERTRIBE = stringr::str_squish(as.character(SUPERTRIBE)),
      FAMILY = stringr::str_squish(as.character(FAMILY))
    )
  
  tip_info <- tip_info[match(tre$tip.label, tip_info$SAMPLE), ]
  
  internal_nodes <- (Ntip + 1L):(Ntip + tre$Nnode)
  
  find_components <- function(rank_col, clade_name, min_component_tips, min_component_fraction) {
    
    focal_idx <- which(
      tip_info[[rank_col]] == clade_name &
        tip_info$FAMILY == "Brassicaceae"
    )
    
    if (length(focal_idx) == 0) return(tibble())
    
    candidates <- lapply(internal_nodes, function(nd) {
      tips_nd <- desc_list[[nd]]
      tips_nd <- tips_nd[tips_nd <= Ntip]
      if (length(tips_nd) == 0) return(NULL)
      
      n_focal <- sum(tips_nd %in% focal_idx)
      if (n_focal < min_component_tips) return(NULL)
      
      frac <- n_focal / length(tips_nd)
      if (is.na(frac) || frac < min_component_fraction) return(NULL)
      
      tibble(
        rank = rank_col,
        clade = clade_name,
        node = nd,
        n_focal_tips = n_focal,
        n_total_tips = length(tips_nd),
        fraction = frac,
        focal_samples = paste(tip_info$SAMPLE[intersect(tips_nd, focal_idx)], collapse = "; ")
      )
    })
    
    cand <- bind_rows(candidates)
    if (nrow(cand) == 0) return(tibble())
    
    ## Greedy non-overlapping selection: largest focal clades first
    cand <- cand %>%
      arrange(desc(n_focal_tips), desc(fraction), n_total_tips)
    
    selected <- list()
    used_tips <- integer(0)
    
    for (i in seq_len(nrow(cand))) {
      nd <- cand$node[i]
      tips_nd <- desc_list[[nd]]
      tips_nd <- tips_nd[tips_nd <= Ntip]
      focal_here <- intersect(tips_nd, focal_idx)
      
      if (length(intersect(focal_here, used_tips)) == 0) {
        selected[[length(selected) + 1L]] <- cand[i, ]
        used_tips <- union(used_tips, focal_here)
      }
    }
    
    out <- bind_rows(selected) %>%
      arrange(desc(n_focal_tips), desc(fraction))
    
    if (nrow(out) > 1) {
      out <- out %>%
        mutate(
          component_id = dplyr::row_number(),
          plot_label = paste0(clade, "_", LETTERS[component_id])
        )
    } else {
      out <- out %>%
        mutate(
          component_id = 1L,
          plot_label = clade
        )
    }
    
    out
  }
  
  tribes <- tip_info %>%
    filter(FAMILY == "Brassicaceae") %>%
    pull(TRIBE) %>%
    unique() %>%
    na.omit() %>%
    sort()
  
  tribes <- tribes[!tolower(tribes) %in% "unplaced"]
  
  tribe_components <- bind_rows(lapply(
    tribes,
    find_components,
    rank_col = "TRIBE",
    min_component_tips = component_min_tips,
    min_component_fraction = component_min_fraction
  )) %>%
    mutate(
      plot_this = n_focal_tips >= min_tips | component_id > 1
    )
  
  tribe_plot <- tribe_components %>%
    filter(plot_this, !is.na(node))
  
  message("[INFO] geom_cladelab(): plotting ", nrow(tribe_plot), " tribe/components.")
  
  ## Large pale crown-node icons behind tribe/component crown nodes
  layout_df <- as.data.frame(p$data)
  
  tribe_icon_df <- layout_df %>%
    dplyr::filter(node %in% tribe_plot$node) %>%
    dplyr::select(node, x, y) %>%
    dplyr::left_join(tribe_plot, by = "node")
  
  if (nrow(tribe_icon_df) > 0) {
    p <- p +
      geom_point(
        data = tribe_icon_df,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        shape = 21,
        fill = "grey55",
        colour = NA,
        alpha = tribe_crown_icon_alpha,
        size = tribe_crown_icon_size
      )
  }
  
  ## Tribe/component bars
  for (i in seq_len(nrow(tribe_plot))) {
    nd <- tribe_plot$node[i]
    if (is.na(nd)) next
    
    p <- p +
      ggtree::geom_cladelab(
        node = nd,
        label = tribe_plot$plot_label[i],
        align = TRUE,
        offset = offset,
        offset.text = 0.4,
        barsize = barsize,
        fontsize = textsize,
        angle = 0
      )
  }
  
  ## Supertribe bars, wider and further right
  supertribes <- tip_info %>%
    filter(FAMILY == "Brassicaceae") %>%
    pull(SUPERTRIBE) %>%
    unique() %>%
    na.omit() %>%
    sort()
  
  supertribes <- supertribes[!tolower(supertribes) %in% "unplaced"]
  
  supertribe_components <- bind_rows(lapply(
    supertribes,
    find_components,
    rank_col = "SUPERTRIBE",
    min_component_tips = 3,
    min_component_fraction = 0.75
  ))
  
  supertribe_plot <- supertribe_components %>%
    filter(!is.na(node), clade %in% names(st_cols))
  
  layout_df <- as.data.frame(p$data)
  tree_x_max <- max(layout_df$x, na.rm = TRUE)
  
  for (i in seq_len(nrow(supertribe_plot))) {
    
    nd <- supertribe_plot$node[i]
    if (is.na(nd)) next
    
    clade_name <- supertribe_plot$clade[i]
    clade_col <- unname(st_cols[clade_name])
    if (is.na(clade_col) || length(clade_col) == 0) clade_col <- "grey40"
    
    node_rows <- layout_df %>%
      filter(node == nd)
    
    if (nrow(node_rows) == 0) next
    
    x_node <- node_rows$x[1]
    
    desc <- phytools::getDescendants(tre, nd)
    desc_tips <- desc[desc <= Ntip]
    
    yvals <- layout_df$y[match(desc_tips, layout_df$node)]
    y_min <- min(yvals, na.rm = TRUE)
    y_max <- max(yvals, na.rm = TRUE)
    
    x_bar <- tree_x_max + supertribe_bar_offset
    
    ## coloured vertical bar
    p <- p +
      geom_segment(
        data = tibble(x_bar = x_bar, y_min = y_min, y_max = y_max),
        aes(x = x_bar, xend = x_bar, y = y_min, yend = y_max),
        inherit.aes = FALSE,
        linewidth = supertribe_bar_barsize,
        colour = clade_col,
        lineend = "round"
      )
    
    ## text label
    p <- p +
      geom_text(
        data = tibble(
          x_lab = x_bar + 0.35,
          y_lab = mean(c(y_min, y_max)),
          lab = clade_name
        ),
        aes(x = x_lab, y = y_lab, label = lab),
        inherit.aes = FALSE,
        hjust = 0,
        size = supertribe_bar_textsize,
        colour = clade_col
      )
  }
  
  attr(p, "tribe_nodes_all") <- tribe_components
  attr(p, "tribe_nodes_plot") <- tribe_plot
  attr(p, "supertribe_nodes_plot") <- supertribe_plot
  
  p
}

draw_nuclear_backbone_tree <- TRUE

backbone_out_pdf <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_backbone_main_lineages.pdf")
)

backbone_out_png <- file.path(
  nuclear_out_dir,
  paste0(project_name, "_", calibration_source, "_backbone_main_lineages.png")
)

plot_backbone_main_lineages <- function(td_merged, specimen_details,
                                        out_pdf, out_png,
                                        width = 5.5, height = 3.5, dpi = 400) {
  
  tre <- td_merged@phylo
  Ntip <- ape::Ntip(tre)
  
  td_tbl <- tidytree::as_tibble(td_merged) %>%
    tibble::as_tibble()
  
  tip_info <- specimen_details %>%
    mutate(
      SAMPLE = as.character(SAMPLE),
      FAMILY = stringr::str_squish(as.character(FAMILY)),
      TRIBE = stringr::str_squish(as.character(TRIBE)),
      SUPERTRIBE = stringr::str_squish(as.character(SUPERTRIBE)),
      loci_remaining = suppressWarnings(as.numeric(loci_remaining))
    ) %>%
    filter(SAMPLE %in% tre$tip.label)
  
  ## Pick one good representative per lineage.
  ## Cleomaceae = best available outgroup representative.
  rep_tbl <- bind_rows(
    tip_info %>%
      filter(FAMILY == "Cleomaceae") %>%
      arrange(desc(loci_remaining)) %>%
      slice(1) %>%
      mutate(backbone_label = "Cleomaceae outgroup",
             backbone_group = "Cleomaceae outgroup"),
    
    tip_info %>%
      filter(FAMILY == "Brassicaceae", TRIBE == "Aethionemeae") %>%
      arrange(desc(loci_remaining)) %>%
      slice(1) %>%
      mutate(backbone_label = "Aethionemeae",
             backbone_group = "Aethionemeae"),
    
    tip_info %>%
      filter(FAMILY == "Brassicaceae", SUPERTRIBE %in% names(st_cols)) %>%
      group_by(SUPERTRIBE) %>%
      arrange(desc(loci_remaining), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(backbone_label = SUPERTRIBE,
             backbone_group = SUPERTRIBE)
  ) %>%
    filter(!is.na(SAMPLE), SAMPLE %in% tre$tip.label)
  
  if (nrow(rep_tbl) < 4) {
    stop("[ERR] Too few backbone representatives found. Check FAMILY/TRIBE/SUPERTRIBE labels.")
  }
  
  message("[INFO] Backbone representatives:")
  print(rep_tbl %>% select(backbone_label, SAMPLE, loci_remaining))
  
  ## Prune to representatives.
  keep_tips <- rep_tbl$SAMPLE
  bb_phy <- ape::keep.tip(tre, keep_tips)
  
  ## Remove branch-length meaning: every edge gets equal length.
  bb_phy$edge.length <- rep(1, nrow(bb_phy$edge))
  
  ## Replace sample IDs by lineage names.
  bb_phy$tip.label <- rep_tbl$backbone_label[match(bb_phy$tip.label, rep_tbl$SAMPLE)]
  
  ## Plot base tree.
  p <- ggtree::ggtree(bb_phy, branch.length = "none", linewidth = 0.6)
  
  p <- p %<+% rep_tbl
  
  p <- p +
    geom_tiplab(
      aes(label = label, colour = label),
      size = 3,
      align = FALSE,
      offset = 0.25
    ) +
    scale_colour_manual(
      values = c("Cleomaceae outgroup" = "black", st_cols),
      guide = "none"
    ) +
    theme_tree2() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  ## Add support values to internal nodes.
  layout_df <- as.data.frame(p$data)
  
  get_tip_labels_below <- function(phy, node) {
    desc <- phytools::getDescendants(phy, node)
    tip_idx <- desc[desc <= ape::Ntip(phy)]
    phy$tip.label[tip_idx]
  }
  
  ## Need mapping from lineage labels back to original sample labels.
  label_to_sample <- setNames(rep_tbl$SAMPLE, rep_tbl$backbone_label)
  
  node_support_df <- layout_df %>%
    filter(!isTip) %>%
    rowwise() %>%
    mutate(
      desc_lineages = list(get_tip_labels_below(bb_phy, node)),
      desc_samples = list(unname(label_to_sample[desc_lineages])),
      full_node = if (length(desc_samples) >= 2) {
        ape::getMRCA(tre, desc_samples)
      } else {
        NA_integer_
      }
    ) %>%
    ungroup() %>%
    left_join(
      td_tbl %>%
        transmute(
          full_node = node,
          q1 = suppressWarnings(as.numeric(q1)),
          pp1 = suppressWarnings(as.numeric(pp1)),
          gCF = suppressWarnings(as.numeric(gCF)),
          sCF = suppressWarnings(as.numeric(sCF))
        ),
      by = "full_node"
    ) %>%
    mutate(
      support_label = case_when(
        !is.na(pp1) & !is.na(q1) & !is.na(gCF) & !is.na(sCF) ~
          paste0(
            "LPP=", round(pp1, 2),
            "\nQ1=", round(q1, 2),
            "\ngCF=", round(gCF, 1),
            "\nsCF=", round(sCF, 1)
          ),
        
        !is.na(pp1) & !is.na(q1) ~
          paste0("LPP=", round(pp1, 2), "\nQ1=", round(q1, 2)),
        
        !is.na(gCF) & !is.na(sCF) ~
          paste0("gCF=", round(gCF, 1), "\nsCF=", round(sCF, 1)),
        
        TRUE ~ ""
      )
    )
  
  p <- p +
    geom_text(
      data = node_support_df %>% filter(support_label != ""),
      aes(x = x, y = y, label = support_label),
      inherit.aes = FALSE,
      hjust = 1.1,
      size = 2.4
    )
  
  ggsave(out_pdf, p, width = width, height = height, limitsize = FALSE, bg = "white")
  ggsave(out_png, p, width = width, height = height, dpi = dpi, limitsize = FALSE, bg = "white")
  
  message("[OK] Backbone mini-tree written: ", out_pdf)
  invisible(p)
}

# ---------------- ROBUST INTERVAL SCALING ----------------

scale_numeric_vector_if_possible <- function(x, factor) {
  xx <- suppressWarnings(as.numeric(x))
  if (length(xx) == length(x) && any(!is.na(xx))) {
    return(xx * factor)
  }
  x
}

scale_interval_entry <- function(x, factor) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) return(x)
  
  if (is.list(x) && length(x) == 1) {
    x <- x[[1]]
  }
  
  if (is.numeric(x)) {
    return(x * factor)
  }
  
  if (is.character(x) && length(x) == 1 && !is.na(x)) {
    nums <- stringr::str_extract_all(
      x,
      "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
    )[[1]]
    nums <- suppressWarnings(as.numeric(nums))
    
    if (length(nums) >= 2 && all(!is.na(nums[1:2]))) {
      nums[1:2] <- nums[1:2] * factor
      return(nums[1:2])
    }
  }
  
  x
}

scale_interval_column <- function(col, factor) {
  if (is.numeric(col)) {
    return(col * factor)
  }
  
  if (is.list(col)) {
    return(lapply(col, scale_interval_entry, factor = factor))
  }
  
  if (is.character(col)) {
    out <- lapply(col, scale_interval_entry, factor = factor)
    out_chr <- vapply(
      out,
      function(z) {
        if (is.numeric(z) && length(z) >= 2 && all(!is.na(z[1:2]))) {
          paste(z[1:2], collapse = " ")
        } else if (is.character(z) && length(z) == 1) {
          z
        } else {
          NA_character_
        }
      },
      character(1)
    )
    return(out_chr)
  }
  
  col
}

rescale_tree_obj <- function(tree_obj, factor = 100) {
  if (is.null(tree_obj)) return(tree_obj)
  
  if (!is.null(tree_obj$phy$edge.length)) {
    tree_obj$phy$edge.length <- tree_obj$phy$edge.length * factor
  }
  
  if (!is.null(tree_obj$td@phylo$edge.length)) {
    tree_obj$td@phylo$edge.length <- tree_obj$td@phylo$edge.length * factor
  }
  
  if (is.null(tree_obj$td@data) || nrow(tree_obj$td@data) == 0) {
    return(tree_obj)
  }
  
  td <- tree_obj$td@data
  td_cols <- names(td)
  
  direct_time_cols <- intersect(
    c(
      "length", "height", "height_95%_HPD", "length_95%_HPD",
      "height_median", "height_mean", "branch.length"
    ),
    td_cols
  )
  
  for (cc in direct_time_cols) {
    if (is.numeric(td[[cc]])) {
      td[[cc]] <- td[[cc]] * factor
    }
  }
  
  interval_cols <- grep(
    "HPD|hpd|interval|Interval|CI|credible|Credible",
    td_cols,
    value = TRUE
  )
  
  if (length(interval_cols) > 0) {
    for (cc in interval_cols) {
      td[[cc]] <- scale_interval_column(td[[cc]], factor)
    }
  }
  
  other_numeric_time_cols <- grep(
    "age|Age|time|Time|height|Height|length|Length",
    td_cols,
    value = TRUE
  )
  
  other_numeric_time_cols <- setdiff(other_numeric_time_cols, interval_cols)
  
  for (cc in other_numeric_time_cols) {
    if (is.numeric(td[[cc]])) {
      td[[cc]] <- td[[cc]] * factor
    }
  }
  
  tree_obj$td@data <- td
  tree_obj
}

# ---------------- NODE MAPPING / MERGING ----------------

map_nodes_by_splits <- function(src_phy, dst_phy, unroot_before_match = TRUE,
                                label = "node mapping") {
  if (unroot_before_match) {
    uu <- function(x) if (inherits(x, "phylo") && ape::is.rooted(x)) ape::unroot(x) else x
    src_phy <- uu(src_phy)
    dst_phy <- uu(dst_phy)
  }
  
  m <- phytools::matchNodes(dst_phy, src_phy)
  
  out <- if (is.matrix(m) && ncol(m) >= 2) {
    tibble(dst = as.integer(m[, 1]), src = as.integer(m[, 2]))
  } else if (is.data.frame(m) && ncol(m) >= 2) {
    tibble(dst = as.integer(m[[1]]), src = as.integer(m[[2]]))
  } else {
    stop("[ERR] Unexpected return from phytools::matchNodes().")
  }
  
  n_dst_internal <- dst_phy$Nnode
  n_src_internal <- src_phy$Nnode
  
  message(
    "[INFO] ", label, ": matched ",
    nrow(out), " nodes; dst internal=", n_dst_internal,
    ", src internal=", n_src_internal
  )
  
  if (nrow(out) < 0.95 * min(n_dst_internal, n_src_internal)) {
    warning(
      "[WARN] Low node-mapping success for ", label,
      ". This may indicate topology differences, collapsed branches, or rooting/tip-label issues."
    )
  }
  
  out
}

pick_cols <- function(tbl, patterns) {
  cols <- names(tbl)
  keep <- unique(unlist(lapply(patterns, function(p) grep(p, cols, ignore.case = TRUE, value = TRUE))))
  keep <- setdiff(keep, c("parent", "label", "isTip"))
  unique(c("node", keep))
}

merge_safe_by_node <- function(x, y) {
  if (is.null(y) || nrow(y) == 0L) return(x)
  
  common <- intersect(names(x), names(y))
  common <- setdiff(common, "node")
  
  y2 <- y
  if (length(common) > 0) {
    names(y2)[match(common, names(y2))] <- paste0(common, ".y")
  }
  
  dplyr::left_join(x, y2, by = "node")
}

# ---------------- ASTRAL/CF READER (NUCLEAR) ----------------

read_cf_astral_tree <- function(path) {
  if (!file.exists(path)) {
    stop("[ERR] Tree file not found: ", path)
  }
  
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  txt <- trimws(txt)
  
  if (!grepl(";", txt, fixed = TRUE)) {
    stop("[ERR] No semicolon found in tree file.")
  }
  
  base_pat <- "'\\[[^]]*\\]'"
  m <- gregexpr(base_pat, txt, perl = TRUE)[[1]]
  
  if (length(m) == 1 && m[1] == -1) {
    stop("[ERR] No quoted ASTRAL blocks found in tree.")
  }
  
  match_len <- attr(m, "match.length")
  starts <- as.integer(m)
  ends <- starts + match_len - 1L
  
  extend_one <- function(start, end, txt) {
    n <- nchar(txt)
    pos <- end + 1L
    while (pos <= n) {
      ch <- substr(txt, pos, pos)
      if (ch %in% c(":", ",", ")", ";")) break
      pos <- pos + 1L
    }
    c(start = start, end = pos - 1L)
  }
  
  spans <- Map(extend_one, starts, ends, MoreArgs = list(txt = txt))
  spans <- do.call(rbind, spans)
  
  raw_labels <- vapply(
    seq_len(nrow(spans)),
    function(i) substr(txt, spans[i, "start"], spans[i, "end"]),
    character(1)
  )
  
  placeholders <- sprintf("NODELAB_%06d", seq_along(raw_labels))
  
  txt2 <- txt
  ord <- order(spans[, "start"], decreasing = TRUE)
  
  for (k in seq_along(ord)) {
    i <- ord[k]
    s <- spans[i, "start"]
    e <- spans[i, "end"]
    ph <- placeholders[i]
    
    left  <- if (s > 1) substr(txt2, 1, s - 1) else ""
    right <- if (e < nchar(txt2)) substr(txt2, e + 1, nchar(txt2)) else ""
    txt2 <- paste0(left, ph, right)
  }
  
  phy <- ape::read.tree(text = txt2)
  
  if (is.null(phy$node.label)) {
    stop("[ERR] Tree read succeeded, but no internal node labels were found.")
  }
  
  idx <- match(phy$node.label, placeholders)
  restored_labels <- rep(NA_character_, length(phy$node.label))
  restored_labels[!is.na(idx)] <- raw_labels[idx[!is.na(idx)]]
  phy$node.label <- restored_labels
  
  parse_one_label <- function(x) {
    if (is.na(x) || x == "") {
      return(tibble::tibble(
        raw_label = x,
        pp1 = NA_real_, pp2 = NA_real_, pp3 = NA_real_,
        f1 = NA_real_, f2 = NA_real_, f3 = NA_real_,
        q1 = NA_real_, q2 = NA_real_, q3 = NA_real_,
        gCF = NA_real_, sCF = NA_real_
      ))
    }
    
    astral_block <- stringr::str_extract(x, "\\[[^]]*\\]")
    
    get_key <- function(key) {
      if (is.na(astral_block)) return(NA_real_)
      mm <- stringr::str_match(astral_block, paste0(key, "=([-0-9.eE]+)"))
      suppressWarnings(as.numeric(mm[, 2]))
    }
    
    tail_part <- sub("^'\\[[^]]*\\]'", "", x)
    tail_part <- sub("^/", "", tail_part)
    
    parts <- strsplit(tail_part, "/", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    
    gcf <- if (length(parts) >= 1) suppressWarnings(as.numeric(parts[1])) else NA_real_
    scf <- if (length(parts) >= 2) suppressWarnings(as.numeric(parts[2])) else NA_real_
    
    tibble::tibble(
      raw_label = x,
      pp1 = get_key("pp1"),
      pp2 = get_key("pp2"),
      pp3 = get_key("pp3"),
      f1  = get_key("f1"),
      f2  = get_key("f2"),
      f3  = get_key("f3"),
      q1  = get_key("q1"),
      q2  = get_key("q2"),
      q3  = get_key("q3"),
      gCF = gcf,
      sCF = scf
    )
  }
  
  parsed <- dplyr::bind_rows(lapply(phy$node.label, parse_one_label)) %>%
    dplyr::mutate(node = ape::Ntip(phy) + seq_len(phy$Nnode)) %>%
    dplyr::select(node, dplyr::everything()) %>%
    tibble::as_tibble()
  
  td <- treeio::as.treedata(phy)
  base_tbl <- tidytree::as_tibble(td) %>% tibble::as_tibble()
  merged_tbl <- base_tbl %>% dplyr::left_join(parsed, by = "node")
  td@data <- merged_tbl
  
  list(
    td = td,
    phy = phy,
    path = path,
    type = "cf_astral_combo"
  )
}

extract_astral_support <- function(td_astral, map_astral_to_base) {
  tbl_astral <- as_tibble(td_astral)
  
  astral_cols <- pick_cols(tbl_astral, c("^q1$", "^pp1$", "quartet", "lpp"))
  tbl_astral_keep <- tbl_astral[, intersect(astral_cols, names(tbl_astral)), drop = FALSE] %>%
    dplyr::rename(src = node)
  
  map_astral_to_base %>%
    dplyr::inner_join(tbl_astral_keep, by = "src") %>%
    dplyr::rename(node = dst) %>%
    dplyr::select(node, dplyr::everything(), -src)
}

extract_cf_support <- function(td_cf, map_cf_to_base) {
  tbl_cf <- as_tibble(td_cf)
  
  if (!all(c("node", "gCF", "sCF") %in% names(tbl_cf))) {
    warning("[WARN] gCF and/or sCF columns not found in CF tree.")
    return(tibble(node = integer(), gCF = numeric(), sCF = numeric()))
  }
  
  df_cf_src <- tbl_cf %>%
    dplyr::filter(!is.na(gCF) | !is.na(sCF)) %>%
    dplyr::transmute(
      src = node,
      gCF = suppressWarnings(as.numeric(gCF)),
      sCF = suppressWarnings(as.numeric(sCF))
    )
  
  map_cf_to_base %>%
    dplyr::inner_join(df_cf_src, by = "src") %>%
    dplyr::transmute(node = dst, gCF, sCF)
}


# ---------------- GENERAL SMALL HELPERS ----------------

extract_binomial <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  out <- ifelse(
    is.na(x) | x == "",
    NA_character_,
    sub("^([A-Za-z][A-Za-z-]*\\s+[A-Za-z][A-Za-z-]*).*$", "\\1", x)
  )
  out
}

normalize_growth <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  x <- stringr::str_to_upper(x)
  dplyr::case_when(
    x %in% c("W", "W/L", "W/T") ~ "W",
    x == "H" ~ "H",
    TRUE ~ NA_character_
  )
}

growth_cols <- c(
  W = "#892255",
  H = "#44AA9A"
)

first_nonempty <- function(...) {
  xs <- list(...)
  n <- max(vapply(xs, length, integer(1)))
  xs <- lapply(xs, function(z) {
    z <- as.character(z)
    length(z) <- n
    z
  })
  out <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    vals <- vapply(xs, function(z) z[i], character(1))
    vals <- vals[!is.na(vals) & str_squish(vals) != ""]
    out[i] <- if (length(vals) > 0) vals[1] else NA_character_
  }
  out
}

bin_three_levels <- function(x, breaks = c(0.5, 0.75), percent_scale = FALSE) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  
  out[!is.na(x) & x <  breaks[1]] <- "low"
  out[!is.na(x) & x >= breaks[1] & x < breaks[2]] <- "medium"
  out[!is.na(x) & x >= breaks[2]] <- "high"
  
  factor(out, levels = c("low", "medium", "high"))
}


# ---------------- SHARED METADATA PREP ----------------

prepare_shared_metadata <- function(td_merged,
                                    target_rank,
                                    target_clade,
                                    taxonomy_check_threshold,
                                    genes_retained_col = "loci_remaining") {
  
  td_merged <- ensure_label_column(td_merged)
  
  specimen_details <- readr::read_csv(specimen_file, show_col_types = FALSE)
  species_details  <- readr::read_csv(species_file, show_col_types = FALSE)
  
  specimen_details_local <- specimen_details %>%
    mutate(
      SAMPLE = as.character(SAMPLE),
      SPECIES_NAME_PRINT = as.character(SPECIES_NAME_PRINT),
      GENUS = as.character(GENUS),
      TRIBE = as.character(TRIBE),
      SUPERTRIBE = as.character(SUPERTRIBE),
      FAMILY = as.character(FAMILY),
      Type_status = as.character(Type_status),
      loci_remaining = suppressWarnings(as.numeric(loci_remaining))
    )
  
  if (!genes_retained_col %in% names(specimen_details_local)) {
    warning("[WARN] genes_retained_col='", genes_retained_col, "' not found; using loci_remaining instead.")
    genes_retained_col <- "loci_remaining"
  }
  
  specimen_details_local[[genes_retained_col]] <- suppressWarnings(as.numeric(specimen_details_local[[genes_retained_col]]))
  specimen_details_local[[genes_retained_col]][is.na(specimen_details_local[[genes_retained_col]])] <- 0
  
  species_details_local <- species_details %>%
    mutate(
      SPECIES_NAME_PRINT = as.character(SPECIES_NAME_PRINT),
      BIONOMIAL = as.character(BIONOMIAL),
      TRIBE_FULL = if ("TRIBE_FULL" %in% names(.)) as.character(TRIBE_FULL) else NA_character_,
      SUPERTRIBE = if ("SUPERTRIBE" %in% names(.)) as.character(SUPERTRIBE) else NA_character_,
      FAMILY = if ("FAMILY" %in% names(.)) as.character(FAMILY) else NA_character_,
      GROWTH_FORM = if ("GROWTH_FORM" %in% names(.)) as.character(GROWTH_FORM) else NA_character_,
      growth_norm_species = normalize_growth(GROWTH_FORM)
    )
  
  level_col <- grep(
    "WGSRPD.*LEVEL[_ ]?1.*native|WCVP.*LEVEL[_ ]?1.*native|LEVEL[_ ]?1.*native",
    names(species_details_local),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(level_col) == 0) {
    warning("[WARN] Could not find WGSRPD Level-1 native column; geography heatmap will be empty.")
    species_details_local$WCVP_L1_native <- NA_character_
  } else {
    message("[INFO] Using geography column: ", level_col[1])
    species_details_local$WCVP_L1_native <- as.character(species_details_local[[level_col[1]]])
  }
  
  species_lookup <- species_details_local %>%
    transmute(
      SPECIES_NAME_PRINT = as.character(SPECIES_NAME_PRINT),
      BIONOMIAL = as.character(BIONOMIAL),
      WCVP_L1_native = as.character(WCVP_L1_native),
      TRIBE_species = as.character(TRIBE_FULL),
      SUPERTRIBE_species = as.character(SUPERTRIBE),
      FAMILY_species = as.character(FAMILY),
      growth_norm_species = growth_norm_species
    ) %>%
    distinct()
  
  specimen_details_local <- specimen_details_local %>%
    mutate(
      SPECIES_NAME_PRINT = as.character(SPECIES_NAME_PRINT)
    ) %>%
    left_join(
      species_lookup,
      by = "SPECIES_NAME_PRINT"
    ) %>%
    mutate(
      species_old = SPECIES_NAME_PRINT,
      growth_norm = growth_norm_species,
      species_for_checks = species_old,
      species_binom_old = extract_binomial(species_old),
      species_binom_for_checks = species_binom_old,
      genus_for_checks = word(species_binom_for_checks, 1),
      tribe_for_checks = coalesce(TRIBE_species, TRIBE),
      supertribe_for_checks = coalesce(SUPERTRIBE_species, SUPERTRIBE),
      family_for_checks = coalesce(FAMILY_species, FAMILY),
      label_color = ifelse(!is.na(Type_status) & Type_status != "", "red", "black"),
      previous_sample = FALSE,

      tip_label = if_else(
        previous_sample,
        paste0("* ", species_old, " (", SAMPLE, ")"),
        paste0(species_old, " (", SAMPLE, ")")
      )
    )
  
  td_merged@data$tip_label <- specimen_details_local$tip_label[match(td_merged@data$label, specimen_details_local$SAMPLE)]
  td_merged@data$label_color <- specimen_details_local$label_color[match(td_merged@data$label, specimen_details_local$SAMPLE)]
  td_merged@data$growth_norm <- specimen_details_local$growth_norm[match(td_merged@data$label, specimen_details_local$SAMPLE)]
  
  tip_order <- td_merged@phylo$tip.label
  
  taxonomy_check_warnings <- data.frame(
    sample = td_merged@phylo$tip.label,
    species = specimen_details_local$species_for_checks[match(td_merged@phylo$tip.label, specimen_details_local$SAMPLE)],
    genus = specimen_details_local$genus_for_checks[match(td_merged@phylo$tip.label, specimen_details_local$SAMPLE)],
    tribe = specimen_details_local$tribe_for_checks[match(td_merged@phylo$tip.label, specimen_details_local$SAMPLE)],
    supertribe = specimen_details_local$supertribe_for_checks[match(td_merged@phylo$tip.label, specimen_details_local$SAMPLE)],
    family = specimen_details_local$family_for_checks[match(td_merged@phylo$tip.label, specimen_details_local$SAMPLE)],
    Warning_same_species_tip_elsewhere = NA_character_,
    Warning_sister_is_different_genus = NA_character_,
    Warning_sister_is_different_tribe = NA_character_,
    stringsAsFactors = FALSE
  )
  
  column_name <- switch(
    target_rank,
    family = "family",
    tribe = "tribe",
    genus = "genus",
    supertribe = "supertribe",
    NULL
  )
  
  tips_target_clade <- character(0)
  
  ## For very large trees, the sister-clade taxonomy check can be slow.
  ## Set target_rank <- NULL to skip it entirely.
  if (!is.null(target_rank) && !is.null(column_name)) {
    tips_target_clade <- taxonomy_check_warnings$sample[
      taxonomy_check_warnings[[column_name]] == target_clade &
        !is.na(taxonomy_check_warnings[[column_name]])
    ]
  }
  
  tips_target_clade_all_species <- taxonomy_check_warnings$species[
    match(tips_target_clade, taxonomy_check_warnings$sample)
  ]
  
  for (tip in tips_target_clade) {
    tip_species <- taxonomy_check_warnings$species[taxonomy_check_warnings$sample == tip]
    tip_genus   <- taxonomy_check_warnings$genus[taxonomy_check_warnings$sample == tip]
    tip_tribe   <- taxonomy_check_warnings$tribe[taxonomy_check_warnings$sample == tip]
    
    tip_number <- which(td_merged@phylo$tip.label == tip)
    parent_node <- td_merged@phylo$edge[td_merged@phylo$edge[, 2] == tip_number, 1]
    
    subtree <- extract.clade(td_merged@phylo, parent_node)
    
    sister_species <- taxonomy_check_warnings$species[match(subtree$tip.label, taxonomy_check_warnings$sample)]
    sister_genera  <- taxonomy_check_warnings$genus[match(subtree$tip.label, taxonomy_check_warnings$sample)]
    sister_tribes  <- taxonomy_check_warnings$tribe[match(subtree$tip.label, taxonomy_check_warnings$sample)]
    
    if (sum(tips_target_clade_all_species %in% tip_species) > 1) {
      if (!sum(sister_species %in% tip_species) == length(sister_species)) {
        taxonomy_check_warnings$Warning_same_species_tip_elsewhere[
          taxonomy_check_warnings$sample == tip
        ] <- "Y"
      }
    }
    
    if (sum(!sister_genera %in% tip_genus, na.rm = TRUE) > 0) {
      fraction <- sum(!sister_genera %in% tip_genus, na.rm = TRUE) / length(sister_genera)
      if (fraction > taxonomy_check_threshold) {
        taxonomy_check_warnings$Warning_sister_is_different_genus[
          taxonomy_check_warnings$sample == tip
        ] <- "Y"
      }
    }
    
    if (sum(!sister_tribes %in% tip_tribe, na.rm = TRUE) > 0) {
      fraction <- sum(!sister_tribes %in% tip_tribe, na.rm = TRUE) / length(sister_tribes)
      if (fraction > taxonomy_check_threshold) {
        taxonomy_check_warnings$Warning_sister_is_different_tribe[
          taxonomy_check_warnings$sample == tip
        ] <- "Y"
      }
    }
  }
  
  td_merged@data$warn_same_species <- taxonomy_check_warnings$Warning_same_species_tip_elsewhere[
    match(td_merged@data$label, taxonomy_check_warnings$sample)
  ]
  td_merged@data$warn_sister_genus <- taxonomy_check_warnings$Warning_sister_is_different_genus[
    match(td_merged@data$label, taxonomy_check_warnings$sample)
  ]
  td_merged@data$warn_sister_tribe <- taxonomy_check_warnings$Warning_sister_is_different_tribe[
    match(td_merged@data$label, taxonomy_check_warnings$sample)
  ]
  
  dat_genes <- tibble(
  SAMPLE = tip_order,
  genes_used = specimen_details_local[[genes_retained_col]][match(tip_order, specimen_details_local$SAMPLE)]
)
  
  data_dist <- tibble(
    SAMPLE = tip_order,
    species_old = specimen_details_local$species_old[match(tip_order, specimen_details_local$SAMPLE)],
    species_binom_for_checks = specimen_details_local$species_binom_for_checks[match(tip_order, specimen_details_local$SAMPLE)],
    WCVP_L1_native_direct = specimen_details_local$WCVP_L1_native[match(tip_order, specimen_details_local$SAMPLE)]
  ) %>%
    mutate(
      WCVP_L1_native_fallback = species_details_local$WCVP_L1_native[
        match(species_binom_for_checks, species_details_local$BIONOMIAL)
      ],
      L1_native = first_nonempty(WCVP_L1_native_direct, WCVP_L1_native_fallback)
    )
  
  data_dist_long <- data_dist %>%
    mutate(row = dplyr::row_number()) %>%
    filter(!is.na(L1_native) & str_squish(L1_native) != "") %>%
    separate_rows(L1_native, sep = ",\\s*") %>%
    mutate(value = "Y")
  
  data_dist_wide <- data_dist_long %>%
    pivot_wider(
      id_cols = row,
      names_from = L1_native,
      values_from = value,
      values_fill = NA
    ) %>%
    arrange(row)
  
  data_dist_final <- data_dist %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::left_join(data_dist_wide, by = "row") %>%
    dplyr::select(-row) %>%
    dplyr::mutate(dplyr::across(
      -c(
        species_old,
        species_binom_for_checks,
        WCVP_L1_native_direct,
        WCVP_L1_native_fallback,
        L1_native
      ),
      ~ tidyr::replace_na(.x, "N")
    )) %>%
    dplyr::select(
      -species_old,
      -species_binom_for_checks,
      -WCVP_L1_native_direct,
      -WCVP_L1_native_fallback,
      -L1_native
    ) %>%
    tibble::column_to_rownames("SAMPLE")
  
  int_colnames <- colnames(data_dist_final)[
    suppressWarnings(!is.na(as.integer(colnames(data_dist_final))))
  ]
  
  if (length(int_colnames) > 0) {
    data_dist_final <- data_dist_final[, int_colnames, drop = FALSE]
  }
  
  valid_names <- suppressWarnings(as.integer(colnames(data_dist_final)))
  valid_cols  <- which(!is.na(valid_names))
  int_names   <- valid_names[valid_cols]
  full_range  <- 1:9
  
  data_dist_numeric_full <- as.data.frame(
    matrix("N", nrow = nrow(data_dist_final), ncol = length(full_range))
  )
  colnames(data_dist_numeric_full) <- as.character(full_range)
  for (i in int_names) {
    cn <- as.character(i)
    data_dist_numeric_full[[cn]] <- data_dist_final[[cn]]
  }
  rownames(data_dist_numeric_full) <- tip_order
  
  for (col in colnames(data_dist_numeric_full)) {
    data_dist_numeric_full[[col]] <- ifelse(
      data_dist_numeric_full[[col]] == "Y", col, "N"
    )
  }
  
  list(
    td_merged = td_merged,
    specimen_details = specimen_details_local,
    species_details = species_details_local,
    taxonomy_check_warnings = taxonomy_check_warnings,
    dat_genes = dat_genes,
    data_dist_numeric_full = data_dist_numeric_full
  )
}

# ---------------- SHARED PLOT COMPONENTS ----------------

continent_colors <- c(
  "1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3", "4" = "#e7298a",
  "5" = "#66a61e", "6" = "orange", "7" = "#a6761d", "8" = "#666666", "9" = "#1f78b4"
)
continent_names <- c(
  "1" = "Europe", "2" = "Africa", "3" = "Asia-Temperate", "4" = "Asia-Tropical",
  "5" = "Australasia", "6" = "Pacific", "7" = "Northern America",
  "8" = "Southern America", "9" = "Antarctic"
)

add_shared_tip_side_panels <- function(p, layout_df, td_merged, dat_genes, data_dist_numeric_full,
                                       use_time_axis = TRUE) {
  layout_df$warn_same_species <- td_merged@data$warn_same_species[match(layout_df$node, td_merged@data$node)]
  layout_df$warn_sister_genus <- td_merged@data$warn_sister_genus[match(layout_df$node, td_merged@data$node)]
  layout_df$warn_sister_tribe <- td_merged@data$warn_sister_tribe[match(layout_df$node, td_merged@data$node)]
  
  # p <- p +
  #   geom_point(
  #     data = subset(layout_df, isTip & warn_same_species == "Y"),
  #     aes(x = x + 0.00, y = y),
  #     shape = 16, color = "blue", fill = "blue",
  #     size = 0.05, inherit.aes = FALSE
  #   ) +
  #   geom_point(
  #     data = subset(layout_df, isTip & warn_sister_genus == "Y"),
  #     aes(x = x + 0.06, y = y),
  #     shape = 16, color = "orange", fill = "orange",
  #     size = 0.05, inherit.aes = FALSE
  #   ) +
  #   geom_point(
  #     data = subset(layout_df, isTip & warn_sister_tribe == "Y"),
  #     aes(x = x + 0.12, y = y),
  #     shape = 16, color = "red", fill = "red",
  #     size = 0.05, inherit.aes = FALSE
  #   )
  
  p <- p + new_scale_fill()
  
  p <- p +
    geom_fruit(
      data = dat_genes,
      geom = geom_bar,
      mapping = aes(y = SAMPLE, x = -genes_used, fill = genes_used),
      orientation = "y",
      stat = "identity",
      offset = gene_bar_offset,
      pwidth = gene_bar_pwidth,
      grid.params = list(lwd = 0)
    ) +
    scale_fill_gradient(
      low = "red",
      high = "green",
      na.value = "white",
      name = "Genes retained",
      guide = guide_colorbar(order = 10)
    )

  
  tree_x_max <- max(layout_df$x, na.rm = TRUE)
  
  tip_xy <- layout_df %>%
    filter(isTip) %>%
    select(SAMPLE = label, y)
  
  geo_all <- data_dist_numeric_full %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SAMPLE") %>%
    tidyr::pivot_longer(
      cols = -SAMPLE,
      names_to = "L1",
      values_to = "value"
    ) %>%
    mutate(
      L1 = as.character(L1)
    ) %>%
    left_join(tip_xy, by = "SAMPLE") %>%
    mutate(
      L1 = factor(L1, levels = names(continent_colors)),
      col_index = as.integer(as.character(L1)),
      x = tree_x_max + geo_heat_offset + (col_index - 0.5) * geo_col_spacing
    )
  
  geo_Y <- geo_all %>%
    filter(!is.na(value), value != "N") %>%
    mutate(
      value = factor(as.character(value), levels = names(continent_colors))
    )
  
  p <- p +
    ggplot2::geom_tile(
      data = geo_all,
      aes(x = x, y = y),
      width = geo_col_spacing,
      height = geo_tile_height,
      fill = "white",
      colour = "grey90",
      linewidth = 0.1,
      inherit.aes = FALSE
    )
  
  p <- p + ggnewscale::new_scale_fill()
  p <- p +
    ggplot2::geom_tile(
      data = geo_Y,
      aes(x = x, y = y, fill = value),
      width = geo_col_spacing,
      height = geo_tile_height,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name = "Native distribution",
      values = continent_colors,
      breaks = names(continent_colors),
      labels = continent_names,
      na.value = "white",
      drop = FALSE,
      guide = guide_legend(order = 11, ncol = 1)
    )
  
  cont_lab_df <- geo_all %>%
    group_by(L1) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(L1))
  
  y_top <- max(tip_xy$y, na.rm = TRUE)
  label_y <- y_top + 2
  cont_lab_df$label <- continent_names[as.character(cont_lab_df$L1)]
  
  p <- p +
    ggplot2::geom_text(
      data = cont_lab_df,
      aes(x = x, y = label_y, label = label),
      angle = 45,
      vjust = 0,
      hjust = 0,
      size = geo_label_size,
      inherit.aes = FALSE
    )
  
  geo_x_max <- if (nrow(geo_all) > 0) max(geo_all$x, na.rm = TRUE) else max(layout_df$x, na.rm = TRUE)
  extra_right_space <- 8
  
  tree_x_min <- min(layout_df$x, na.rm = TRUE)
  tree_x_max <- max(layout_df$x, na.rm = TRUE)
  tree_x_span <- tree_x_max - tree_x_min
  left_pad <- max(0.03 * tree_x_span, 0.05)
  
  if (isTRUE(use_time_axis)) {
    final_x_limits <- c(min(time_axis_limits), geo_x_max + extra_right_space)
  } else {
    final_x_limits <- c(tree_x_min - left_pad, geo_x_max + extra_right_space)
  }
  
  p <- p +
    coord_cartesian(xlim = final_x_limits, clip = "off") +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = "grey70"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 8, r = 120, b = 8, l = 8)
    )
  
  p
}

standardize_lsd2_ci_columns <- function(td_merged) {
  td <- td_merged@data
  cols <- names(td)
  
  ## Prefer explicit LSD2/FigTree-style CI columns
  lower_candidates <- grep("(^|_)(lower|lo|min|ci_min|ci_lower)($|_)", cols, ignore.case = TRUE, value = TRUE)
  upper_candidates <- grep("(^|_)(upper|hi|max|ci_max|ci_upper)($|_)", cols, ignore.case = TRUE, value = TRUE)
  
  lower_col <- lower_candidates[1]
  upper_col <- upper_candidates[1]
  
  if (!is.na(lower_col) && !is.na(upper_col) &&
      lower_col %in% cols && upper_col %in% cols) {
    
    lo <- suppressWarnings(as.numeric(td[[lower_col]]))
    hi <- suppressWarnings(as.numeric(td[[upper_col]]))
    
    if (any(!is.na(lo)) && any(!is.na(hi))) {
      td$hpd_95 <- Map(c, lo, hi)
      message("[INFO] Created hpd_95 from columns: ",
              lower_col, " + ", upper_col)
    }
  }
  
  td_merged@data <- td
  td_merged
}

prepare_hpd_column_if_present <- function(td_merged) {
  td_cols <- names(td_merged@data)
  
  hpd_col <- c(
    "hpd_95",
    "CI_date",
    "CI_height",
    "height_0.95_HPD",
    "length_0.95_HPD",
    grep("95.*HPD|HPD.*95|HPD|hpd|CI", td_cols, value = TRUE)
  ) %>%
    unique() %>%
    .[. %in% td_cols] %>%
    .[1]
  
  if (length(hpd_col) == 0 || is.na(hpd_col)) {
    return(list(td_merged = td_merged, hpd_col = NA_character_))
  }
  
  if (hpd_col %in% names(td_merged@data)) {
    names(td_merged@data)[names(td_merged@data) == hpd_col] <- "hpd_95"
    hpd_col <- "hpd_95"
  }
  
  list(td_merged = td_merged, hpd_col = hpd_col)
}

# ---------------- NUCLEAR PLOT BUILDER ----------------
# Less monolithic version:
#   1. prepare_nuclear_plot_inputs()
#   2. build_nuclear_plot()
#   3. save_nuclear_plot()
#   4. plot_nuclear_tree() wrapper retained for backward compatibility

prepare_nuclear_plot_inputs <- function(td_merged,
                                        q1_breaks,
                                        scf_breaks,
                                        lpp_high_threshold) {
  
  message("[INFO] Preparing nuclear plot inputs...")
  
  td_merged <- standardize_lsd2_ci_columns(td_merged)
  hpd_info <- prepare_hpd_column_if_present(td_merged)
  td_merged <- hpd_info$td_merged
  hpd_col <- hpd_info$hpd_col
  
  td_tbl <- tidytree::as_tibble(td_merged)
  td_cols <- names(td_tbl)
  
  col_lpp <- grep("^pp1$|lpp", td_cols, ignore.case = TRUE, value = TRUE)[1]
  col_q1  <- grep("^q1$", td_cols, ignore.case = TRUE, value = TRUE)[1]
  col_scf <- if ("sCF" %in% td_cols) "sCF" else NA_character_
  
  message("[INFO] Support columns detected:")
  message("  q1 : ", col_q1  %||% "NOT FOUND")
  message("  lpp: ", col_lpp %||% "NOT FOUND")
  message("  sCF: ", col_scf %||% "NOT FOUND")
  message("  CI : ", hpd_col %||% "NOT FOUND")
  
  node_data <- td_tbl %>%
    mutate(
      q1_num  = if (!is.na(col_q1)  && col_q1  %in% names(.)) suppressWarnings(as.numeric(.data[[col_q1]])) else NA_real_,
      lpp_num = if (!is.na(col_lpp) && col_lpp %in% names(.)) suppressWarnings(as.numeric(.data[[col_lpp]])) else NA_real_,
      scf_num = if (!is.na(col_scf) && col_scf %in% names(.)) suppressWarnings(as.numeric(.data[[col_scf]])) else NA_real_,
      q1_bin  = bin_three_levels(q1_num, breaks = q1_breaks, percent_scale = FALSE),
      sCF_bin = bin_three_levels(scf_num, breaks = scf_breaks, percent_scale = TRUE),
      lpp_bin = case_when(
        is.na(lpp_num) ~ NA_character_,
        lpp_num >= lpp_high_threshold ~ "LPP ≥ 0.95",
        TRUE ~ "LPP < 0.95"
      ),
      lpp_bin = factor(lpp_bin, levels = c("LPP < 0.95", "LPP ≥ 0.95"))
    )
  
  node_data2 <- tibble::as_tibble(node_data)
  class(node_data2) <- c("tbl_df", "tbl", "data.frame")
  td_merged@data <- node_data2
  
  prep <- prepare_shared_metadata(
    td_merged = td_merged,
    target_rank = target_rank,
    target_clade = target_clade,
    taxonomy_check_threshold = taxonomy_check_threshold,
    genes_retained_col = genes_retained_col
  )
  
  message("[INFO] Nuclear plot inputs prepared.")
  
  list(
    td_merged = prep$td_merged,
    hpd_col = hpd_col,
    dat_genes = prep$dat_genes,
    data_dist_numeric_full = prep$data_dist_numeric_full,
    taxonomy_check_warnings = prep$taxonomy_check_warnings,
    specimen_details = prep$specimen_details,
    species_details = prep$species_details,
    detected_columns = list(
      q1 = col_q1,
      lpp = col_lpp,
      sCF = col_scf,
      hpd = hpd_col
    )
  )
}

build_nuclear_plot <- function(nuc_inputs) {
  
  message("[INFO] Building nuclear plot object...")
  
  td_merged <- nuc_inputs$td_merged
  hpd_col <- nuc_inputs$hpd_col
  dat_genes <- nuc_inputs$dat_genes
  data_dist_numeric_full <- nuc_inputs$data_dist_numeric_full
  
  p <- ggtree(td_merged, size = 0.03)
  
  if (!is.na(hpd_col) && hpd_col %in% names(td_merged@data)) {
    p <- p +
      geom_range(
        range = hpd_col,
        color = "red",
        alpha = 0.55,
        linewidth = 0.15
      )
    
    ## Move most recently added layer to the back, so CI bars sit behind branches
    p$layers <- c(p$layers[length(p$layers)], p$layers[-length(p$layers)])
  }
  
  p <- p +
    geom_tiplab(
      aes(label = tip_label, color = label_color),
      size = tip_label_size,
      offset = 0.02,
      align = FALSE
    ) +
    scale_color_identity() +
    theme_tree2()
  
  p <- revts(p)
  
  if (isTRUE(draw_tribe_cladelabs)) {
    p <- add_tribe_cladelabs(
      p = p,
      td_merged = td_merged,
      specimen_details = nuc_inputs$specimen_details,
      min_tips = tribe_cladelab_min_tips,
      offset = tribe_cladelab_offset,
      barsize = tribe_cladelab_barsize,
      textsize = tribe_cladelab_textsize
    )
    
    tribe_nodes_all <- attr(p, "tribe_nodes_all")
    
    if (!is.null(tribe_nodes_all)) {
      
      readr::write_csv(
        tribe_nodes_all,
        file.path(
          nuclear_out_dir,
          paste0(project_name, "_tribe_cladelab_status.csv")
        )
      )
    }
    
  }
  
  if ("label.x" %in% names(p$data) && "label.y" %in% names(p$data)) {
    if (identical(p$data$label.x, p$data$label.y)) {
      p$data$label <- p$data$label.x
    }
  }
  
  layout_df <- as.data.frame(p$data)
  layout_df$isTip <- as.logical(layout_df$isTip)
  
  growth_df <- layout_df %>%
    dplyr::filter(isTip) %>%
    dplyr::mutate(
      growth_norm = td_merged@data$growth_norm[
        match(node, td_merged@data$node)
      ]
    ) %>%
    dplyr::filter(!is.na(growth_norm))
  
  node_df <- layout_df %>%
    dplyr::filter(!isTip) %>%
    dplyr::select(node, x, y) %>%
    dplyr::left_join(
      td_merged@data %>%
        dplyr::select(node, q1_num, lpp_num, scf_num, q1_bin, lpp_bin, sCF_bin),
      by = "node"
    )
  
  p <- p + ggnewscale::new_scale_fill() +
    geom_point(
      data = node_df %>% filter(!is.na(sCF_bin)),
      aes(x = x, y = y, fill = sCF_bin),
      inherit.aes = FALSE,
      shape = 22,
      size = node_icon_size_square,
      stroke = 0,
      alpha = 0.8
    ) +
    scale_fill_manual(
      name = "Site concordance factor (sCF)",
      values = c(low = "red", medium = "yellow", high = "darkgreen"),
      breaks = c("high", "medium", "low"),
      labels = c(
        high = "high (≥ 66)",
        medium = "medium (33–66)",
        low = "low (< 33)"
      ),
      drop = FALSE,
      guide = guide_legend(order = 1, override.aes = list(shape = 22, size = 8))
    )
  
  p <- p + ggnewscale::new_scale_fill() + ggnewscale::new_scale_color() +
    geom_point(
      data = node_df %>% filter(!is.na(q1_bin) | !is.na(lpp_bin)),
      aes(x = x, y = y, fill = q1_bin, colour = lpp_bin),
      inherit.aes = FALSE,
      shape = 21,
      colour = "black",
      size = node_icon_size_circle,
      stroke = node_icon_stroke
    ) +
    scale_fill_manual(
      name = "Main quartet support (Q1)",
      values = c(low = "red", medium = "yellow", high = "darkgreen"),
      breaks = c("high", "medium", "low"),
      labels = c(
        high = "high (≥ 0.75)",
        medium = "medium (0.50–0.75)",
        low = "low (< 0.50)"
      ),
      drop = FALSE,
      guide = guide_legend(order = 2, override.aes = list(shape = 21, size = 7, colour = "black"))
    ) +
    scale_colour_manual(
      name = "Local posterior probability (LPP)",
      values = c("LPP < 0.95" = "transparent", "LPP ≥ 0.95" = "black"),
      breaks = c("LPP ≥ 0.95"),
      labels = c("LPP ≥ 0.95" = "high (≥ 0.95)"),
      drop = FALSE,
      guide = guide_legend(order = 3, override.aes = list(shape = 21, fill = "white", size = 7, stroke = 1.2))
    )
  
  p <- p +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = growth_df,
      aes(x = x + 0.35, y = y, fill = growth_norm),
      inherit.aes = FALSE,
      shape = 22,
      size = 0.8,
      colour = "black",
      stroke = 0.1
    ) +
    scale_fill_manual(
      name = "Growth form",
      values = growth_cols,
      breaks = c("W", "H"),
      labels = c(W = "Woody", H = "Herbaceous"),
      drop = FALSE
    )
  
  p <- add_shared_tip_side_panels(
    p = p,
    layout_df = layout_df,
    td_merged = td_merged,
    dat_genes = dat_genes,
    data_dist_numeric_full = data_dist_numeric_full
  )
  
  message("[INFO] Nuclear plot object built.")
  
  p
}

save_nuclear_plot <- function(p,
                              out_pdf,
                              out_png,
                              plot_width,
                              plot_height,
                              plot_dpi) {
  
  message("[INFO] Saving nuclear plot...")
  
  ggsave(
    filename = out_pdf,
    plot = p,
    width = plot_width,
    height = plot_height,
    limitsize = FALSE,
    bg = "white"
  )
  
  ggsave(
    filename = out_png,
    plot = p,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi,
    limitsize = FALSE,
    bg = "white"
  )
  
  message("[OK] Nuclear plot written: ", out_pdf)
  message("[OK] Nuclear plot written: ", out_png)
  
  invisible(TRUE)
}

plot_nuclear_tree <- function(td_merged,
                              out_pdf,
                              out_png,
                              plot_width,
                              plot_height,
                              plot_dpi,
                              q1_breaks,
                              scf_breaks,
                              lpp_high_threshold,
                              debug_out_dir = NULL,
                              save_debug_rds = TRUE) {
  
  nuc_inputs <- prepare_nuclear_plot_inputs(
    td_merged = td_merged,
    q1_breaks = q1_breaks,
    scf_breaks = scf_breaks,
    lpp_high_threshold = lpp_high_threshold
  )
  
  if (isTRUE(save_debug_rds) && !is.null(debug_out_dir)) {
    dir.create(debug_out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      nuc_inputs,
      file.path(debug_out_dir, "DEBUG_01_nuclear_plot_inputs.rds")
    )
    message("[INFO] Debug RDS written: ",
            file.path(debug_out_dir, "DEBUG_01_nuclear_plot_inputs.rds"))
  }
  
  p <- build_nuclear_plot(nuc_inputs)
  
  if (isTRUE(save_debug_rds) && !is.null(debug_out_dir)) {
    saveRDS(
      p,
      file.path(debug_out_dir, "DEBUG_02_nuclear_plot_object.rds")
    )
    message("[INFO] Debug RDS written: ",
            file.path(debug_out_dir, "DEBUG_02_nuclear_plot_object.rds"))
  }
  
  save_nuclear_plot(
    p = p,
    out_pdf = out_pdf,
    out_png = out_png,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_dpi = plot_dpi
  )
  
  invisible(p)
}



# ---------------- NUCLEAR SUPPORT-OVER-TIME HELPERS ----------------

get_node_age_table <- function(td_merged) {
  tre <- td_merged@phylo
  
  node_depth <- ape::node.depth.edgelength(tre)
  tree_height <- max(node_depth, na.rm = TRUE)
  node_age_ma <- tree_height - node_depth
  
  Ntip <- ape::Ntip(tre)
  
  tbl <- tidytree::as_tibble(td_merged) %>%
    tibble::as_tibble()
  
  tbl %>%
    mutate(
      isTip = node <= Ntip,
      age_Ma = node_age_ma[node],
      q1_num  = if ("q1"  %in% names(.)) suppressWarnings(as.numeric(q1))  else NA_real_,
      pp1_num = if ("pp1" %in% names(.)) suppressWarnings(as.numeric(pp1)) else NA_real_,
      gCF_num = if ("gCF" %in% names(.)) suppressWarnings(as.numeric(gCF)) else NA_real_,
      sCF_num = if ("sCF" %in% names(.)) suppressWarnings(as.numeric(sCF)) else NA_real_
    )
}

summarise_support_over_time <- function(node_tbl, bin_width = 2,
                                        min_age = 0,
                                        max_age = 95) {
  node_tbl %>%
    filter(!isTip) %>%
    filter(!is.na(age_Ma), age_Ma >= min_age, age_Ma <= max_age) %>%
    mutate(
      age_bin = floor(age_Ma / bin_width) * bin_width,
      age_mid = age_bin + bin_width / 2
    ) %>%
    group_by(age_bin, age_mid) %>%
    summarise(
      n_nodes = n(),
      q1_med  = median(q1_num,  na.rm = TRUE),
      q1_lo   = quantile(q1_num,  0.025, na.rm = TRUE),
      q1_hi   = quantile(q1_num,  0.975, na.rm = TRUE),
      pp1_med = median(pp1_num, na.rm = TRUE),
      pp1_lo  = quantile(pp1_num, 0.025, na.rm = TRUE),
      pp1_hi  = quantile(pp1_num, 0.975, na.rm = TRUE),
      gCF_med = median(gCF_num, na.rm = TRUE),
      gCF_lo  = quantile(gCF_num, 0.025, na.rm = TRUE),
      gCF_hi  = quantile(gCF_num, 0.975, na.rm = TRUE),
      sCF_med = median(sCF_num, na.rm = TRUE),
      sCF_lo  = quantile(sCF_num, 0.025, na.rm = TRUE),
      sCF_hi  = quantile(sCF_num, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}




# -------------------- NUCLEAR RUN ---------------------

if (isTRUE(draw_nuclear_tree)) {
  message("[INFO] ---------------- NUCLEAR TREE ----------------")
  dir.create(nuclear_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  calibrated_tree <- switch(
    calibration_source,
    mcmctree = resolve_calibrated_tree_path(
      RFIN = RFIN,
      project_name = project_name,
      mct_run_tag = mct_run_tag,
      mct_profile = mct_profile,
      mct_replicate = mct_replicate,
      calibrated_tree_manual = calibrated_tree_manual
    ),
    lsd2 = resolve_lsd2_tree_path(
      RFIN = RFIN,
      project_name = project_name,
      lsd2_run_tag = lsd2_run_tag,
      lsd2_n_loci = lsd2_n_loci,
      lsd2_tree_manual = lsd2_tree_manual
    ),
    stop("[ERR] Unknown calibration_source: ", calibration_source)
  )
  
  message("[INFO] Nuclear project: ", project_name)
  message("[INFO] Nuclear primary tree type: ", primary_tree_type)
  message("[INFO] Nuclear calibration source: ", calibration_source)
  message("[INFO] Tree paths:")
  message("  calibrated: ", calibrated_tree %||% "NOT FOUND / NOT SET")
  message("  astral    : ", astral_tree)
  message("  cf        : ", if (file.exists(cf_tree)) cf_tree else "NOT FOUND")
  
  astral_obj <- read_any_tree(astral_tree, type = "astral")
  
  calibrated_obj <- NULL
  if (!is.null(calibrated_tree) && file.exists(calibrated_tree)) {
    calibrated_obj <- read_any_tree(
      calibrated_tree,
      type = calibration_source,
      force_ultrametric_mcmctree = FALSE
    )
    
    if (calibration_source == "lsd2") {
      message("[INFO] LSD2 annotation columns detected:")
      message("  ", paste(names(calibrated_obj$td@data), collapse = ", "))
    }
    calibrated_obj <- rescale_tree_obj(
      calibrated_obj,
      factor = calibrated_time_scaling_factor
    )
  }
  
  cf_obj <- NULL
  if (!is.null(cf_tree) && file.exists(cf_tree)) {
    cf_obj <- read_cf_astral_tree(cf_tree)
  } else {
    message("[INFO] CF tree not found; continuing without concordance-factor annotations.")
  }
  
  primary_obj <- switch(
    primary_tree_type,
    calibrated = {
      if (is.null(calibrated_obj)) {
        stop("[ERR] primary_tree_type='calibrated' but calibrated tree could not be found.")
      }
      calibrated_obj
    },
    astral = astral_obj,
    cf = {
      if (is.null(cf_obj)) stop("[ERR] primary_tree_type='cf' but CF tree could not be found.")
      cf_obj
    }
  )
  
  primary_phy <- primary_obj$phy
  primary_phy <- root_tree_safely(
    primary_phy,
    root_method = root_method,
    outgroup_tips = outgroup_tips
  )
  
  if (primary_tree_type %in% c("astral", "cf")) {
    primary_phy <- make_ultrametric_if_requested(
      primary_phy,
      do_it = make_primary_ultrametric,
      method = ultrametric_method
    )
  }
  
  if (primary_tree_type == "calibrated" && root_method == "keep") {
    base_phy <- primary_phy
    td_base  <- primary_obj$td
  } else {
    base_phy <- primary_phy
    td_base  <- treeio::as.treedata(base_phy)
  }
  
  tip_base   <- sort(base_phy$tip.label)
  tip_astral <- sort(astral_obj$phy$tip.label)
  
  if (!identical(tip_base, tip_astral)) {
    stop("[ERR] Tip sets differ between primary nuclear tree and ASTRAL tree.")
  }
  
  if (!is.null(cf_obj)) {
    tip_cf <- sort(cf_obj$phy$tip.label)
    if (!identical(tip_base, tip_cf)) {
      warning("[WARN] Tip sets differ between primary nuclear tree and CF tree; CF mapping may fail or be partial.")
    }
  }
  
  map_astral_to_base <- map_nodes_by_splits(
    src_phy = astral_obj$phy,
    dst_phy = base_phy,
    unroot_before_match = TRUE,
    label = "ASTRAL-to-calibrated"
  )
  
  astral_on_base <- extract_astral_support(
    td_astral = astral_obj$td,
    map_astral_to_base = map_astral_to_base
  )
  
  cf_on_base <- NULL
  if (!is.null(cf_obj)) {
    cf_tip_overlap <- intersect(base_phy$tip.label, cf_obj$phy$tip.label)
    if (length(cf_tip_overlap) == length(base_phy$tip.label)) {
      map_cf_to_base <- map_nodes_by_splits(
        src_phy = cf_obj$phy,
        dst_phy = base_phy,
        unroot_before_match = TRUE,
        label = "CF-to-calibrated"
      )
      cf_on_base <- extract_cf_support(
        td_cf = cf_obj$td,
        map_cf_to_base = map_cf_to_base
      )
    } else {
      warning("[WARN] CF tree does not have the same full tip set as primary tree. Skipping CF mapping.")
    }
  }
  
  td_tab <- as_tibble(td_base)
  td_tab <- td_tab %>% merge_safe_by_node(astral_on_base)
  if (!is.null(cf_on_base)) {
    td_tab <- td_tab %>% merge_safe_by_node(cf_on_base)
  }
  td_merged_nuclear <- as.treedata(td_tab)
  
  plot_nuclear_tree(
    td_merged = td_merged_nuclear,
    out_pdf = nuclear_out_pdf,
    out_png = nuclear_out_png,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_dpi = plot_dpi,
    q1_breaks = q1_breaks,
    scf_breaks = scf_breaks,
    lpp_high_threshold = lpp_high_threshold,
    debug_out_dir = nuclear_out_dir,
    save_debug_rds = TRUE
  )
  
  
  if (isTRUE(draw_nuclear_support_time)) {
    message("[INFO] ---------------- NUCLEAR SUPPORT OVER TIME ----------------")
    
    node_tbl <- get_node_age_table(td_merged_nuclear)
    
    support_time_tbl <- summarise_support_over_time(
      node_tbl,
      bin_width = support_time_bin_width,
      min_age = support_time_min_age,
      max_age = support_time_max_age
    )
    
    rank_node_tbl <- make_rank_mrca_table(
      td_merged = td_merged_nuclear,
      specimen_details = specimen_details,
      ranks = c("FAMILY", "SUPERTRIBE", "TRIBE"),
      family_filter = "Brassicaceae"
    )
    
    readr::write_csv(
      node_tbl,
      file.path(nuclear_out_dir, paste0(project_name, "_", calibration_source, "_node_support_ages.csv"))
    )
    
    readr::write_csv(
      support_time_tbl,
      file.path(nuclear_out_dir, paste0(project_name, "_", calibration_source, "_support_over_time_summary.csv"))
    )
    
    readr::write_csv(
      rank_node_tbl,
      support_time_node_table_out
    )
    
    p_support_time <- plot_support_over_time_annotated(
      support_time_tbl = support_time_tbl,
      highlight_tbl = rank_node_tbl,
      xlim = c(-support_time_max_age, 0)
    )
    
    ggsave(
      support_time_out_pdf,
      p_support_time,
      width = 16,
      height = 14,
      limitsize = FALSE,
      bg = "white"
    )
    
    ggsave(
      support_time_out_png,
      p_support_time,
      width = 16,
      height = 14,
      dpi = plot_dpi,
      limitsize = FALSE,
      bg = "white"
    )
    
    message("[OK] Support-over-time plot written: ", support_time_out_pdf)
    message("[OK] Support/age table written: ", support_time_node_table_out)
  }
  
  if (isTRUE(draw_nuclear_backbone_tree)) {
    plot_backbone_main_lineages(
      td_merged = td_merged_nuclear,
      specimen_details = specimen_details,
      out_pdf = backbone_out_pdf,
      out_png = backbone_out_png,
      dpi = plot_dpi
    )
  }
  
  message("[INFO] Nuclear tree finished.")
}
