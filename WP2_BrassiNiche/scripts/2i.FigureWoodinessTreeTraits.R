#!/usr/bin/env Rscript

## 2i.FigureWoodinessTree_Main.R
## Main circular Brassicaceae/Cleomaceae tree:
## - LSD2-calibrated tree from 2f/8_objects.rds
## - branch colours = major clades
## - outer ring = growth form
## - confident internal shifts highlighted + labelled bubbles
## - radial Ma grid

out_dir <- "WP2_BrassiNiche/results_final/6_BrassiToL_traits_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Figure_circular_BC_tree_growthform_shifts.pdf")
out_png <- file.path(out_dir, "Figure_circular_BC_tree_growthform_shifts.png")

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(phangorn)
  library(readr)
})

## Settings ----

OPEN_ANGLE_DEG <- 9

col_hw <- "#892255"  # H>W
col_wh <- "#44AA9A"  # W>H

lighten_col <- function(col, amount = 0.45) {
  rgb <- grDevices::col2rgb(col) / 255
  rgb_new <- rgb + (1 - rgb) * amount
  grDevices::rgb(rgb_new[1], rgb_new[2], rgb_new[3])
}

col_hw_bubble <- lighten_col(col_hw, 0.45)
col_wh_bubble <- lighten_col(col_wh, 0.45)

clade_cols <- c(
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8",
  "Aethionemeae"      = "#ffed57",
  "Cleomaceae"        = "grey70",
  "Unplaced"          = "grey75"
)

clade_cols_light <- vapply(clade_cols, lighten_col, character(1), amount = 0.35)

growth_cols <- c(
  H = col_wh,
  W = col_hw
)

## Load current 2f objects ----

src_path <- "WP2_BrassiNiche/results_final/3_woodiness_evolution/8_objects.rds"
stopifnot(file.exists(src_path))
src <- readRDS(src_path)

need <- c("tree_bc", "woody_fac_bc", "confident_shifts_all")
miss <- setdiff(need, names(src))
if (length(miss) > 0) {
  stop("[ERROR] Missing objects in 2f RDS: ", paste(miss, collapse = ", "))
}

tree_bc <- src$tree_bc
woody_fac_bc <- src$woody_fac_bc[tree_bc$tip.label]
confident_shifts_all <- src$confident_shifts_all

if (any(tree_bc$edge.length == 0, na.rm = TRUE)) {
  tree_bc$edge.length[tree_bc$edge.length == 0] <- 1e-10
}

tree_height <- max(
  ape::node.depth.edgelength(tree_bc)[seq_along(tree_bc$tip.label)],
  na.rm = TRUE
)

message("[INFO] Loaded tree: tips=", length(tree_bc$tip.label),
        " height=", signif(tree_height, 6), " Ma")

## Shift-label lookup: prefer 8_objects.rds, fallback to Step 5.7 CSV ----

shift_lookup <- NULL

if ("per_shift_age" %in% names(src) && !is.null(src$per_shift_age)) {
  shift_lookup <- src$per_shift_age %>%
    dplyr::filter(!is.na(edge_id), !is.na(shift_label)) %>%
    dplyr::transmute(
      edge_id = as.integer(edge_id),
      shift_label = as.character(shift_label)
    ) %>%
    dplyr::distinct(edge_id, .keep_all = TRUE)
  
  message("[INFO] Loaded shift labels from src$per_shift_age: ", nrow(shift_lookup))
}

if (is.null(shift_lookup) || nrow(shift_lookup) == 0) {
  shift_csv <- "WP2_BrassiNiche/results_final/3_woodiness_evolution/5.7_confident_internal_shift_age_summary.csv"
  
  if (!file.exists(shift_csv)) {
    stop("[ERROR] No shift labels found in src$per_shift_age and missing fallback CSV: ", shift_csv)
  }
  
  shift_lookup <- readr::read_csv(shift_csv, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(edge_id), !is.na(shift_label)) %>%
    dplyr::transmute(
      edge_id = as.integer(edge_id),
      shift_label = as.character(shift_label)
    ) %>%
    dplyr::distinct(edge_id, .keep_all = TRUE)
  
  message("[INFO] Loaded shift labels from fallback CSV: ", nrow(shift_lookup))
}

## Read metadata from local publication CSV ----

specimen_csv <- "WP2_BrassiNiche/data/brassiwood_specimens_publication.csv"
stopifnot(file.exists(specimen_csv))

specimen_details <- readr::read_csv(specimen_csv, show_col_types = FALSE)

needed_specimen_cols <- c(
  "SAMPLE",
  "SPECIES_NAME_PRINT",
  "FAMILY",
  "TRIBE",
  "SUPERTRIBE"
)

missing_specimen_cols <- setdiff(needed_specimen_cols, names(specimen_details))
if (length(missing_specimen_cols) > 0) {
  stop(
    "[ERROR] Missing columns in specimen CSV: ",
    paste(missing_specimen_cols, collapse = ", ")
  )
}

specimen_details <- specimen_details %>%
  dplyr::mutate(
    SAMPLE = as.character(SAMPLE),
    SPECIES_NAME_PRINT = stringr::str_squish(as.character(SPECIES_NAME_PRINT)),
    FAMILY = as.character(FAMILY),
    TRIBE = as.character(TRIBE),
    SUPERTRIBE = as.character(SUPERTRIBE)
  )

tip_labels_bc <- specimen_details %>%
  dplyr::select(SAMPLE, SPECIES_NAME_PRINT, FAMILY, TRIBE, SUPERTRIBE) %>%
  dplyr::distinct(SAMPLE, .keep_all = TRUE) %>%
  dplyr::filter(SAMPLE %in% tree_bc$tip.label) %>%
  dplyr::transmute(
    label = SAMPLE,
    SPECIES_NAME_PRINT,
    FAMILY,
    TRIBE,
    SUPERTRIBE,
    clade_group = dplyr::case_when(
      FAMILY == "Cleomaceae"   ~ "Cleomaceae",
      TRIBE  == "Aethionemeae" ~ "Aethionemeae",
      !is.na(SUPERTRIBE) & SUPERTRIBE != "" ~ SUPERTRIBE,
      TRUE ~ "Unplaced"
    )
  )

tip_state_tbl <- tibble::tibble(
  label = tree_bc$tip.label,
  growth_form = dplyr::case_when(
    as.character(woody_fac_bc) == "1" ~ "W",
    as.character(woody_fac_bc) == "0" ~ "H",
    TRUE ~ NA_character_
  )
) %>%
  dplyr::mutate(
    growth_form = dplyr::case_when(
      growth_form %in% c("H", "Herbaceous", "0") ~ "H",
      growth_form %in% c("W", "Woody", "1")      ~ "W",
      TRUE ~ NA_character_
    )
  )

## Edge annotations ----

edge_index_bc <- tibble::tibble(
  edge_id = seq_len(nrow(tree_bc$edge)),
  parent  = tree_bc$edge[, 1],
  child   = tree_bc$edge[, 2]
)

desc_tip_labels <- function(tree, node) {
  node <- as.integer(node)
  if (node <= length(tree$tip.label)) return(tree$tip.label[node])
  tips <- phangorn::Descendants(tree, node, type = "tips")[[1]]
  tree$tip.label[tips]
}

pick_majority <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

tip_group_vec <- tip_labels_bc$clade_group
names(tip_group_vec) <- tip_labels_bc$label

edge_clades <- edge_index_bc %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    desc_tips = list(desc_tip_labels(tree_bc, child)),
    clade_group = pick_majority(unname(tip_group_vec[desc_tips]))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(edge_id, clade_group)

shift_edges <- confident_shifts_all %>%
  dplyr::transmute(
    edge_id = as.integer(edge_id),
    shift_dir = as.character(node_dir),
    is_terminal = as.logical(is_terminal)
  ) %>%
  dplyr::filter(
    !is_terminal,
    shift_dir %in% c("H>W", "W>H")
  ) %>%
  dplyr::left_join(shift_lookup, by = "edge_id") %>%
  dplyr::distinct(edge_id, .keep_all = TRUE)

message("[INFO] Confident internal shifts: ", nrow(shift_edges))
message("[INFO] Confident internal shifts with labels: ", sum(!is.na(shift_edges$shift_label)))

## Build fan tree ----

p <- suppressWarnings(
  ggtree::ggtree(
    tree_bc,
    layout = "fan",
    open.angle = OPEN_ANGLE_DEG,
    size = 0.15
  )
)

p$layers <- list()

p$data <- p$data %>%
  dplyr::mutate(
    parent = as.integer(parent),
    node   = as.integer(node)
  ) %>%
  dplyr::left_join(
    edge_index_bc %>%
      dplyr::transmute(
        edge_id = as.integer(edge_id),
        parent  = as.integer(parent),
        node    = as.integer(child)
      ),
    by = c("parent", "node")
  ) %>%
  dplyr::left_join(edge_clades, by = "edge_id") %>%
  dplyr::left_join(shift_edges, by = "edge_id")

tree_plot_data <- p$data
tip_radius <- max(tree_plot_data$x[tree_plot_data$isTip], na.rm = TRUE)

## Ma grid rings ----

time_ticks <- seq(0, floor(tree_height / 10) * 10, by = 10)
time_ticks <- time_ticks[time_ticks <= tree_height]

time_grid <- tibble::tibble(
  age_ma = time_ticks,
  x_pos  = tree_height - age_ma,
  label  = paste0(age_ma, " Ma")
)

## Precompute bubble positions before adding geom_fruit ----

node_coords <- tree_plot_data %>%
  dplyr::select(node, x, y) %>%
  dplyr::distinct()

bubble_df <- tree_plot_data %>%
  dplyr::filter(!is.na(shift_dir), !is.na(parent), !is.na(edge_id)) %>%
  dplyr::select(
    edge_id, parent, node, shift_dir, shift_label,
    x_child = x, y_child = y
  ) %>%
  dplyr::left_join(node_coords, by = c("parent" = "node")) %>%
  dplyr::rename(x_parent = x, y_parent = y) %>%
  dplyr::mutate(
    x_mid = (x_parent + x_child) / 2,
    y_mid = y_child
  ) %>%
  dplyr::filter(!is.na(shift_label)) %>%
  dplyr::arrange(y_mid)

message("[INFO] Shift bubbles: ", nrow(bubble_df))

## Grid rings first ----

p <- p +
  ggplot2::geom_vline(
    data = time_grid,
    ggplot2::aes(xintercept = x_pos),
    inherit.aes = FALSE,
    colour = "grey82",
    linewidth = 0.35
  )

## Coloured clade branches ----

p <- p +
  ggtree::geom_tree(
    data = tree_plot_data,
    ggplot2::aes(color = clade_group),
    size = 0.18,
    alpha = 1,
    lineend = "round"
  ) +
  ggplot2::scale_color_manual(
    name = "Clade",
    values = clade_cols_light,
    na.value = "grey85"
  )

## Overlay confident shifts ----

p <- p +
  ggnewscale::new_scale_color() +
  ggtree::geom_tree(
    data = tree_plot_data,
    ggplot2::aes(color = shift_dir),
    size = 0.55,
    alpha = 1,
    lineend = "round"
  ) +
  ggplot2::scale_color_manual(
    name = "Confident shifts",
    values = c("H>W" = col_hw, "W>H" = col_wh),
    breaks = c("H>W", "W>H"),
    labels = c("H>W" = "H→W", "W>H" = "W→H"),
    na.value = NA,
    na.translate = FALSE
  )

## Outer growth-form ring ----

ring_offset <- -0.12
ring_width  <- 0.18
tile_width  <- 4

p <- p +
  ggtree::xlim_tree(tip_radius + 0.20) +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = tip_state_tbl,
    geom = geom_tile,
    mapping = ggplot2::aes(
      y = label,
      x = 1,
      fill = growth_form
    ),
    offset = ring_offset,
    pwidth = ring_width,
    width = tile_width,
    axis.params = list(axis = "none")
  ) +
  ggplot2::scale_fill_manual(
    name = "Growth form",
    values = growth_cols,
    breaks = c("H", "W"),
    labels = c(H = "Herbaceous", W = "Woody"),
    na.value = "grey92"
  )

## Numbered shift bubbles ----

p <- p +
  ggnewscale::new_scale_color() +
  ggplot2::geom_point(
    data = bubble_df,
    ggplot2::aes(x = x_mid, y = y_mid, color = shift_dir),
    inherit.aes = FALSE,
    shape = 16,
    size = 5.6,
    alpha = 1
  ) +
  ggplot2::scale_color_manual(
    name = "Shift labels",
    values = c("H>W" = col_hw_bubble, "W>H" = col_wh_bubble),
    breaks = c("H>W", "W>H"),
    labels = c("H>W" = "H→W", "W>H" = "W→H"),
    na.translate = FALSE
  ) +
  ggplot2::geom_text(
    data = bubble_df,
    ggplot2::aes(x = x_mid, y = y_mid, label = shift_label),
    inherit.aes = FALSE,
    size = 3.0,
    fontface = "bold",
    colour = "black",
    vjust = 0.35
  )

## Time labels LAST so they cannot be hidden by later layers ----

p <- p +
  ggplot2::geom_text(
    data = time_grid,
    ggplot2::aes(x = x_pos, y = 0, label = label),
    inherit.aes = FALSE,
    colour = "grey35",
    size = 3.4,
    hjust = -0.10,
    vjust = 0.5
  )

## Final formatting ----

p <- p +
  ggtree::theme_tree() +
  ggplot2::theme(
    legend.position = "left",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 11),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )

ggplot2::ggsave(out_pdf, p, width = 11.7, height = 11.7)
ggplot2::ggsave(out_png, p, width = 11.7, height = 11.7, dpi = 400)

message("[DONE] Wrote: ", out_pdf)
message("[DONE] Wrote: ", out_png)
