#!/usr/bin/env Rscript

## 2e.EnvironmentalContextFigure.R
##
## Goal:
##   Create a three-panel environmental-context figure for Brassicaceae woodiness:
##     A) global drought map (MCWD)
##     B) global frost map (BIO6 = minimum temperature of coldest month)
##     C) species niche-space hex plot, filled by woody fraction
##
## Inputs:
##   - species_niche_summary_ge10.csv from 2d.RecordsToNiche.R
##   - brassiwood_species_publication.csv
##   - environmental rasters already used in 2d:
##       * WorldClim monthly precipitation
##       * WorldClim monthly mean temperature
##       * Zomer v3 monthly ET0
##       * WorldClim elevation (used only as raster template)
##   - WGSRPD L3 shapefile for outlines
##
## Output:
##   WP2_BrassiNiche/results_final/2_environmental_context/
##     - Figure_5_environmental_context.pdf
##     - Figure_5_environmental_context.png
##
## Notes:
##   - Panel C uses species-level medians from the ge10 tier.
##   - Woody fraction in panel C is capped at FRACTION_MAX for visual contrast.
##   - Panels A and B use coarsened rasters for plotting speed.
##   - This is a first integrated version for visual inspection and refinement.


## STEP 1. PREPARE PATHS, OUTPUTS, PACKAGES ----

out_dir <- "WP2_BrassiNiche/results_final/2_environmental_context"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Figure_5_environmental_context.pdf")
out_png <- file.path(out_dir, "Figure_5_environmental_context.png")

out_sup_pdf <- file.path(
  out_dir,
  "Figure_Sx_main_clade_niche_space.pdf"
)

out_sup_png <- file.path(
  out_dir,
  "Figure_Sx_main_clade_niche_space.png"
)

out_tribe_pdf <- file.path(
  out_dir,
  "Figure_Sx_woody_tribe_niche_space.pdf"
)

out_tribe_png <- file.path(
  out_dir,
  "Figure_Sx_woody_tribe_niche_space.png"
)

in_niche <- "WP2_BrassiNiche/results_final/1_niche_modelling_preparation/species_niche_summary_ge10.csv"

env_dir   <- "WP2_BrassiNiche/data/environmental_data"
cache_dir <- "WP2_BrassiNiche/results_intermediate/niche_cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

wgsrpd_root <- "WP2_BrassiNiche/data/wgsrpd-master"
shp_L3_path <- file.path(wgsrpd_root, "level3/level3.shp")
stopifnot(file.exists(shp_L3_path))

species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
stopifnot(file.exists(species_file))

req <- c(
  "dplyr", "readr", "stringr", "tidyr", "ggplot2",
  "terra", "sf", "patchwork", "scales",
  "hexbin", "grid", "geodata"
)
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop(
    "[ERROR] Missing packages: ", paste(miss, collapse = ", "),
    "\nInstall with install.packages(c(",
    paste(sprintf('\"%s\"', miss), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(terra)
  library(sf)
  library(patchwork)
  library(scales)
  library(grid)
  library(geodata)
})

message("[INFO] Output dir: ", normalizePath(out_dir))


## STEP 2. SETTINGS ----


## Panel C colours used throughout the project
COL_WOODY <- "#892255"
COL_HERB  <- "#44AA9A"

## Panel C fill scale
FRACTION_MAX    <- 0.20
FRACTION_SQUISH <- TRUE
HEX_BINS        <- 28

## Map raster coarsening for plotting speed
MAP_AGG_FACTOR <- 12

## Map colours
COL_Mcwd_LOW  <- "red"
COL_Mcwd_HIGH <- "#f7f7f7"

COL_Bio6_LOW  <- "#2c7fb8"  # cold
COL_Bio6_HIGH <- "#f7f7f7"  # warm / frost-free

## Optional plotting limits for climate maps
## Set to NULL to use raster ranges after aggregation.
MCWD_LIMITS <- c(-3000, 0)
BIO6_LIMITS <- c(-4, 2)

## Global map extent
MAP_XLIM <- c(-180, 180)
MAP_YLIM <- c(-90, 90)

## Theme sizes
BASE_SIZE <- 11


## STEP 3. HELPER FUNCTIONS ----


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

find_monthly_wc <- function(dir, var) {
  if (!dir.exists(dir)) return(character(0))
  pat <- paste0("^wc2\\.1_30s_", var, "_(0[1-9]|1[0-2])\\.tif$")
  list.files(dir, pattern = pat, full.names = TRUE)
}

find_single_wc <- function(dir, var) {
  if (!dir.exists(dir)) return(character(0))
  pat <- paste0("^wc2\\.1_30s_", var, "\\.tif$")
  list.files(dir, pattern = pat, full.names = TRUE)
}

save_plot <- function(p, pdf_file, png_file, w, h) {
  ggsave(pdf_file, plot = p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(png_file, plot = p, width = w, height = h, units = "in", dpi = 300)
}

raster_to_df <- function(r, value_name) {
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", value_name)
  df
}

cap_values <- function(x, limits = NULL) {
  if (is.null(limits)) return(x)
  pmin(pmax(x, limits[1]), limits[2])
}


## STEP 4. LOAD SPECIES METADATA + NICHE SUMMARIES FOR PANEL C ----


stopifnot(file.exists(in_niche))

message("[INFO] Reading species metadata from local publication CSV.")

species_details <- readr::read_csv(
  species_file,
  show_col_types = FALSE
)

niche_tbl <- readr::read_csv(in_niche, show_col_types = FALSE)

needed_sheet <- c(
  "ACCEPTED",
  "FAMILY",
  "SPECIES_NAME_PRINT",
  "SUBSP_VAR",
  "GROWTH_FORM",
  "SUPERTRIBE",
  "TRIBE_FULL"
)
miss_sheet <- setdiff(needed_sheet, names(species_details))
if (length(miss_sheet) > 0) {
  stop("[ERROR] Missing columns in species sheet: ", paste(miss_sheet, collapse = ", "))
}

needed_niche <- c("species", "MCWD_med", "BIO6_med", "n_cells_kept")
miss_niche <- setdiff(needed_niche, names(niche_tbl))
if (length(miss_niche) > 0) {
  stop("[ERROR] Missing columns in niche summary table: ", paste(miss_niche, collapse = ", "))
}

main_clade_cols <- c(
  "Aethionemeae"      = "#ffed57",
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8"
)

species_meta <- species_details %>%
  transmute(
    species = gsub(" ", "_", stringr::str_squish(as.character(SPECIES_NAME_PRINT))),
    ACCEPTED = as.character(ACCEPTED),
    FAMILY = stringr::str_squish(as.character(FAMILY)),
    SUPERTRIBE = stringr::str_squish(as.character(SUPERTRIBE)),
    TRIBE_FULL = stringr::str_squish(as.character(TRIBE_FULL)),
    SUBSP_VAR = as.character(SUBSP_VAR),
    GROWTH_FORM = normalize_growth(GROWTH_FORM)
  ) %>%
  filter(
    ACCEPTED == "Y",
    FAMILY == "Brassicaceae",
    is.na(stringr::str_squish(SUBSP_VAR)) |
      stringr::str_squish(SUBSP_VAR) == ""
  ) %>%
  mutate(
    MAIN_CLADE = case_when(
      stringr::str_detect(TRIBE_FULL, stringr::regex("^Aethionemeae\\b", ignore_case = TRUE)) ~ "Aethionemeae",
      stringr::str_detect(SUPERTRIBE, stringr::regex("Arabodae", ignore_case = TRUE)) ~ "Arabodae (IV)",
      stringr::str_detect(SUPERTRIBE, stringr::regex("Brassicodae", ignore_case = TRUE)) ~ "Brassicodae (II)",
      stringr::str_detect(SUPERTRIBE, stringr::regex("Camelinodae", ignore_case = TRUE)) ~ "Camelinodae (I)",
      stringr::str_detect(SUPERTRIBE, stringr::regex("Heliophilodae", ignore_case = TRUE)) ~ "Heliophilodae (V)",
      stringr::str_detect(SUPERTRIBE, stringr::regex("Hesperodae", ignore_case = TRUE)) ~ "Hesperodae (III)",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(species, .keep_all = TRUE)


plot_tbl_c <- niche_tbl %>%
  left_join(
    species_meta %>%
      dplyr::select(
        species,
        GROWTH_FORM,
        MAIN_CLADE
      ),
    by = "species"
  ) %>%
  filter(!is.na(GROWTH_FORM), GROWTH_FORM %in% c("H", "W")) %>%
  filter(is.finite(MCWD_med), is.finite(BIO6_med)) %>%
  mutate(
    woody01 = if_else(GROWTH_FORM == "W", 1, 0),
    growth_form_label = if_else(GROWTH_FORM == "W", "Woody", "Herbaceous")
  )

message("[DEBUG] MAIN_CLADE values:")
print(sort(unique(plot_tbl_c$MAIN_CLADE)))
print(setdiff(unique(na.omit(plot_tbl_c$MAIN_CLADE)), names(main_clade_cols)))

message("[INFO] Panel C species count: ", nrow(plot_tbl_c))
message("[INFO] Panel C woody species: ", sum(plot_tbl_c$woody01 == 1, na.rm = TRUE))
message("[INFO] Panel C herbaceous species: ", sum(plot_tbl_c$woody01 == 0, na.rm = TRUE))
message("[INFO] Main clade counts in Panel C:")
print(plot_tbl_c %>% count(MAIN_CLADE, sort = TRUE))


## STEP 5. LOAD ENVIRONMENTAL RASTERS FOR PANELS A AND B ----


## 5.1 Elevation as template
elev_path_local <- find_single_wc(env_dir, "elev")
elev_path_cache <- file.path(cache_dir, "wc2.1_30s_elev.tif")

if (length(elev_path_local) == 1) {
  elev <- terra::rast(elev_path_local)
  message("[ENV] Elevation: using local ", elev_path_local)
} else if (file.exists(elev_path_cache)) {
  elev <- terra::rast(elev_path_cache)
  message("[ENV] Elevation: using cache ", elev_path_cache)
} else {
  message("[ENV] Elevation: downloading WorldClim via geodata...")
  elev <- geodata::worldclim_global(var = "elev", res = 0.5, path = cache_dir)
  terra::writeRaster(elev, elev_path_cache, overwrite = TRUE)
  elev <- terra::rast(elev_path_cache)
}

## 5.2 Monthly precipitation and temperature
prec_files <- find_monthly_wc(env_dir, "prec")
tavg_files <- find_monthly_wc(env_dir, "tavg")

prec_cache <- file.path(cache_dir, "wc2.1_30s_prec_monthly.tif")
tavg_cache <- file.path(cache_dir, "wc2.1_30s_tavg_monthly.tif")
et0_cache  <- file.path(cache_dir, "et0_v3_monthly.tif")

if (length(prec_files) == 12) {
  prec <- terra::rast(sort(prec_files))
  message("[ENV] Precip: using local monthly tifs (n=12)")
} else if (file.exists(prec_cache)) {
  prec <- terra::rast(prec_cache)
  message("[ENV] Precip: using cache ", prec_cache)
} else {
  message("[ENV] Precip: downloading WorldClim via geodata...")
  prec <- geodata::worldclim_global(var = "prec", res = 0.5, path = cache_dir)
  terra::writeRaster(prec, prec_cache, overwrite = TRUE)
  prec <- terra::rast(prec_cache)
}

if (length(tavg_files) == 12) {
  tavg <- terra::rast(sort(tavg_files))
  message("[ENV] TAVG: using local monthly tifs (n=12)")
} else if (file.exists(tavg_cache)) {
  tavg <- terra::rast(tavg_cache)
  message("[ENV] TAVG: using cache ", tavg_cache)
} else {
  message("[ENV] TAVG: downloading WorldClim via geodata...")
  tavg <- geodata::worldclim_global(var = "tavg", res = 0.5, path = cache_dir)
  terra::writeRaster(tavg, tavg_cache, overwrite = TRUE)
  tavg <- terra::rast(tavg_cache)
}

stopifnot(terra::nlyr(prec) == 12, terra::nlyr(tavg) == 12)

## WorldClim monthly mean temperature: °C * 10
tavg_c <- tavg / 10

## 5.3 Monthly ET0
if (file.exists(et0_cache)) {
  et0 <- terra::rast(et0_cache)
  message("[ENV] ET0: using cached multi-layer tif ", et0_cache)
} else {
  et0_v3_dir <- file.path(env_dir, "Global-ET0_v3_monthly")
  et0_files_all <- list.files(
    et0_v3_dir,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  et0_files <- et0_files_all[
    grepl("^et0_v3_([0-9]{1,2})\\.tif$", basename(et0_files_all), ignore.case = TRUE)
  ]
  
  if (length(et0_files) == 0) {
    stop("[ENV] No ET0 v3 GeoTIFFs found in: ", et0_v3_dir)
  }
  
  get_month <- function(x) {
    as.integer(sub("^.*et0_v3_([0-9]{1,2})\\.tif$", "\\1", basename(x), ignore.case = TRUE))
  }
  
  months <- get_month(et0_files)
  ord <- order(months)
  et0_files <- et0_files[ord]
  months <- months[ord]
  
  if (!identical(months, 1:12)) {
    stop("[ENV] ET0 monthly files incomplete. Found months: ", paste(months, collapse = ", "))
  }
  
  et0 <- terra::rast(et0_files)
  names(et0) <- paste0("ET0_", sprintf("%02d", 1:12))
  terra::writeRaster(et0, et0_cache, overwrite = TRUE)
  message("[ENV] ET0: loaded from ", et0_v3_dir)
}

stopifnot(terra::nlyr(et0) == 12)


## STEP 6. DERIVE PLOTTING RASTERS FOR MCWD AND BIO6 ----


## Coarsen first for plotting speed
message("[ENV] Aggregating rasters for plotting (factor = ", MAP_AGG_FACTOR, ") ...")

prec_plot <- terra::aggregate(prec,   fact = MAP_AGG_FACTOR, fun = mean, na.rm = TRUE)
tavg_plot <- terra::aggregate(tavg_c, fact = MAP_AGG_FACTOR, fun = mean, na.rm = TRUE)
et0_plot  <- terra::aggregate(et0,    fact = MAP_AGG_FACTOR, fun = mean, na.rm = TRUE)

## A. MCWD = minimum cumulative deficit of pmin(P - ET0, 0)
message("[ENV] Deriving MCWD raster...")
def_plot <- terra::app(prec_plot - et0_plot, fun = function(x) pmin(x, 0))
mcwd_r <- terra::app(def_plot, fun = function(v) {
  if (all(!is.finite(v))) return(NA_real_)
  min(cumsum(v))
})

names(mcwd_r) <- "MCWD"

## B. BIO6 = minimum temperature of coldest month
message("[ENV] Deriving BIO6 raster...")
bio6_r <- terra::app(tavg_plot, fun = function(v) {
  if (all(!is.finite(v))) return(NA_real_)
  min(v)
})

names(bio6_r) <- "BIO6"

## Cap map values for cleaner visual contrast
mcwd_df <- raster_to_df(mcwd_r, "value") %>%
  mutate(value = cap_values(value, MCWD_LIMITS))

bio6_df <- raster_to_df(bio6_r, "value") %>%
  mutate(value = cap_values(value, BIO6_LIMITS))


## STEP 7. LOAD OUTLINE GEOMETRY FOR MAP PANELS ----


wgsrpd_L3 <- sf::st_read(shp_L3_path, quiet = TRUE)
wgsrpd_L3 <- suppressWarnings(sf::st_make_valid(wgsrpd_L3))


## STEP 8. BUILD PANEL A (GLOBAL DROUGHT MAP: MCWD) ----


p_A <- ggplot() +
  geom_raster(
    data = mcwd_df,
    aes(x = x, y = y, fill = value)
  ) +
  geom_sf(
    data = wgsrpd_L3,
    fill = NA,
    color = "grey20",
    linewidth = 0.15
  ) +
  scale_fill_gradient(
    name = "Drought (MCWD)",
    low = COL_Mcwd_LOW,
    high = COL_Mcwd_HIGH,
    na.value = "white",
    limits = MCWD_LIMITS,
    oob = scales::squish
  ) +
  coord_sf(
    crs = "+proj=natearth",
    xlim = MAP_XLIM,
    ylim = MAP_YLIM,
    expand = FALSE
  ) +
  labs(
    title = "B. Global drought",
    subtitle = "Maximum Cumulative Water Deficit (MCWD)"
  ) +
  theme_minimal(base_size = BASE_SIZE) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(4, 4, 4, 4)
  )


## STEP 9. BUILD PANEL B (GLOBAL FROST MAP: BIO6) ----


p_B <- ggplot() +
  geom_raster(
    data = bio6_df,
    aes(x = x, y = y, fill = value)
  ) +
  geom_sf(
    data = wgsrpd_L3,
    fill = NA,
    color = "grey20",
    linewidth = 0.15
  ) +
  scale_fill_gradient(
    name = "Frost (BIO6; °C)",
    low = COL_Bio6_LOW,
    high = COL_Bio6_HIGH,
    na.value = "white",
    limits = BIO6_LIMITS,
    oob = scales::squish
  ) +
  coord_sf(
    crs = "+proj=natearth",
    xlim = MAP_XLIM,
    ylim = MAP_YLIM,
    expand = FALSE
  ) +
  labs(
    title = "C. Global frost",
    subtitle = "Minimum temperature of coldest month (BIO6)"
  ) +
  theme_minimal(base_size = BASE_SIZE) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(4, 4, 4, 4)
  )


## STEP 10. BUILD PANEL C (HEX-BINNED NICHE SPACE) ----


brks <- seq(0, FRACTION_MAX, by = 0.05)
lbls <- scales::percent(brks, accuracy = 1)
if (length(lbls) > 0) {
  lbls[length(lbls)] <- paste0(scales::percent(FRACTION_MAX, accuracy = 1), " (or more)")
}

p_C <- ggplot(plot_tbl_c, aes(x = MCWD_med, y = BIO6_med)) +
  stat_summary_hex(
    aes(z = woody01, fill = after_stat(value)),
    fun = mean,
    bins = HEX_BINS,
    colour = "grey80",
    linewidth = 0.25
  ) +
  geom_point(
    data = tibble::tibble(
      MAIN_CLADE = names(main_clade_cols),
      MCWD_med = NA_real_,
      BIO6_med = NA_real_
    ),
    aes(
      x = MCWD_med,
      y = BIO6_med,
      colour = MAIN_CLADE
    ),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_fill_gradient(
    name = "Woody fraction",
    low = COL_HERB,
    high = COL_WOODY,
    limits = c(0, FRACTION_MAX),
    oob = if (FRACTION_SQUISH) scales::squish else scales::censor,
    breaks = brks,
    labels = lbls,
    na.value = "white"
  ) +
  labs(
    title = "A. Brassicaceae niche space",
    x = "Drought (MCWD; mm water deficit)",
    y = "Frost (BIO6; °C)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = c(0.08, 0.10),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.key.height = unit(22, "pt"),
    legend.key.width  = unit(10, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 6, 6, 6)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt")
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(linewidth = 1.4),
      order = 1
    ),
    fill = guide_colorbar(
      title.position = "top",
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt"),
      order = 2
    )
  )

## Optional: uncomment if you prefer drier conditions on the right
# p_C <- p_C + scale_x_reverse()


## STEP 10B. BUILD SUPPLEMENTARY MAIN-CLADE HEX PANELS ----

main_clade_levels <- names(main_clade_cols)

plot_tbl_sup <- plot_tbl_c %>%
  filter(MAIN_CLADE %in% main_clade_levels) %>%
  mutate(
    MAIN_CLADE = factor(
      MAIN_CLADE,
      levels = main_clade_levels
    )
  )

p_sup <- ggplot(
  plot_tbl_sup,
  aes(x = MCWD_med, y = BIO6_med)
) +
  stat_summary_hex(
    aes(
      z = woody01,
      fill = after_stat(value)
    ),
    fun = mean,
    bins = HEX_BINS,
    colour = "grey85",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ MAIN_CLADE,
    ncol = 2
  ) +
  scale_fill_gradient(
    name = "Woody fraction",
    low = COL_HERB,
    high = COL_WOODY,
    limits = c(0, FRACTION_MAX),
    oob = if (FRACTION_SQUISH) scales::squish else scales::censor,
    breaks = brks,
    labels = lbls,
    na.value = "white"
  ) +
  labs(
    x = "Drought (MCWD; mm water deficit)",
    y = "Frost (BIO6; °C)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey80"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "right",
    panel.spacing = unit(1.0, "lines"),
    aspect.ratio = 1
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top"
    )
  )

save_plot(
  p = p_sup,
  pdf_file = out_sup_pdf,
  png_file = out_sup_png,
  w = 10,
  h = 12
)

message("[DONE] Saved supplementary main-clade niche figure:")
message("  ", out_sup_pdf)
message("  ", out_sup_png)


## STEP 10C. BUILD SUPPLEMENTARY WOODY-TRIBE HEX PANELS ----

woody_tribes <- sort(c(
  "Alysseae",
  "Brassiceae",
  "Anastaticeae",
  "Lepidieae",
  "Heliophileae",
  "Thelypodieae",
  "Aethionemeae",
  "Arabideae",
  "Erysimeae",
  "Microlepidieae"
))


plot_tbl_tribes <- niche_tbl %>%
  left_join(
    species_meta %>%
      dplyr::select(
        species,
        GROWTH_FORM,
        TRIBE_FULL
      ),
    by = "species"
  ) %>%
  filter(
    TRIBE_FULL %in% woody_tribes,
    !is.na(GROWTH_FORM),
    GROWTH_FORM %in% c("H", "W"),
    is.finite(MCWD_med),
    is.finite(BIO6_med)
  ) %>%
  mutate(
    woody01 = if_else(GROWTH_FORM == "W", 1, 0),
    TRIBE_FULL = factor(
      TRIBE_FULL,
      levels = woody_tribes
    )
  )

p_tribes <- ggplot(
  plot_tbl_tribes,
  aes(x = MCWD_med, y = BIO6_med)
) +
  stat_summary_hex(
    aes(
      z = woody01,
      fill = after_stat(value)
    ),
    fun = mean,
    bins = HEX_BINS,
    colour = "grey85",
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ TRIBE_FULL,
    ncol = 3
  ) +
  scale_fill_gradient(
    name = "Woody fraction",
    low = COL_HERB,
    high = COL_WOODY,
    limits = c(0, FRACTION_MAX),
    oob = if (FRACTION_SQUISH) scales::squish else scales::censor,
    breaks = brks,
    labels = lbls,
    na.value = "white"
  ) +
  labs(
    x = "Drought (MCWD; mm water deficit)",
    y = "Frost (BIO6; °C)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey80"
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.spacing = unit(0.9, "lines"),
    aspect.ratio = 1
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top"
    )
  )

save_plot(
  p = p_tribes,
  pdf_file = out_tribe_pdf,
  png_file = out_tribe_png,
  w = 8,
  h = 10
)

message("[DONE] Saved woody-tribe supplementary figure:")
message("  ", out_tribe_pdf)
message("  ", out_tribe_png)


## STEP 11. COMBINE PANELS ----

## Rename panels conceptually (no need to change earlier code)
## p_C = new Panel A (main)
## p_A = new Panel B (drought)
## p_B = new Panel C (frost)

right_col <- p_A / p_B + patchwork::plot_layout(heights = c(1, 1))

Figure5 <- p_C | right_col

## STEP 12. SAVE FIGURE ----


save_plot(
  p = Figure5,
  pdf_file = out_pdf,
  png_file = out_png,
  w = 13.5,
  h = 8.5
)

message("[DONE] Saved:")
message("  ", out_pdf)
message("  ", out_png)

