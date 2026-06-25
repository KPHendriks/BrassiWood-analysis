# 2a.IntroductionOverview.R
# Global Brassicaceae distribution with woody emphasis + continent waffle composition
#
# Output:
#   WP2_BrassiNiche/results_final/0_introduction_and_overview/Figure_1_Global_distribution_woody_emphasis.{pdf,png}

## STEP 0: PREPARE -----

wd <- "~/Google Drive/My Drive/Publications/2026_Hendriks_et_al_BrassiWood/BrassiWood/"
setwd(wd)

out_dir <- "WP2_BrassiNiche/results_final/0_introduction_and_overview"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Figure_1_Global_distribution_woody_emphasis.pdf")
out_png <- file.path(out_dir, "Figure_1_Global_distribution_woody_emphasis.png")

species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
stopifnot(file.exists(species_file))

wgsrpd_root <- "WP2_BrassiNiche/data/wgsrpd-master"
shp_L3_path <- file.path(wgsrpd_root, "level3/level3.shp")
stopifnot(file.exists(shp_L3_path))

## Figure options
USE_MOLLWEIDE <- TRUE

SHADING_MODE     <- "woody_fraction"  # "woody_richness" | "woody_fraction"
FRACTION_MAX     <- 0.20
FRACTION_SQUISH  <- TRUE

## Remove problematic continents for geometry artefacts
DROP_L1_NAMES_FROM_MAP <- c("Antarctic")

## Remove geometry artefacts (“fat dots”)
DROP_NONPOLYGON_GEOMS <- TRUE
DROP_TINY_GEOMS       <- TRUE
TINY_GEOM_MIN_AREA_M2 <- 5e7  # 50 km^2 (tweak if you want)

## Globe outline
ADD_GLOBE_OUTLINE   <- TRUE
GLOBE_OUTLINE_LWD   <- 0.2
GLOBE_OUTLINE_COL   <- "black"
GLOBE_OUTLINE_ALPHA <- 0.9

## L1 continent outlines (dissolved)
ADD_L1_CONTOURS  <- TRUE
USE_AUTOMATIC_CONTINENT_BORDERS <- FALSE
L1_CONTOUR_LWD   <- 0.45 #was 0.2
L1_CONTOUR_COL   <- "black"
L1_CONTOUR_ALPHA <- 0.9
L1_CONTOUR_EXCLUDE_NAMES <- c("Pacific", "Antarctic")

## Waffles (single row below map)
ADD_L1_WAFFLES <- TRUE
WAFFLE_USE_FRACTIONAL_TILES <- FALSE
WAFFLE_ROUNDING_MODE <- "ceiling"  # "ceiling" | "round" | "floor"

WAFFLE_TILE_VALUE_SPP       <- 20
WAFFLE_NCOL                 <- 5
WAFFLE_EMPTY_FILL           <- "white"
WAFFLE_TILE_GAP_LINEWIDTH   <- 0.15

WAFFLE_WOODY_MARK <- "dot"  # "none" | "outline" | "dot" | "overlay"
WAFFLE_WOODY_DOT_SIZE <- 0.8
WAFFLE_WOODY_DOT_ALPHA <- 0.9
WAFFLE_WOODY_OUTLINE_COL     <- "black"
WAFFLE_WOODY_OUTLINE_LWD     <- 0.35
WAFFLE_WOODY_OUTLINE_INSET   <- 0.78
WAFFLE_PARTIAL_FINAL_TILE <- TRUE
WAFFLE_WOODY_DOT_MODE <- "presence"

WAFFLE_TITLE_SIZE <- 4.5  # halved-ish

WAFFLE_L1_ORDER_NAMES <- c(
  "Pacific",
  "Northern America",
  "Southern America",
  "Europe",
  "Africa",
  "Asia-Temperate",
  "Asia-Tropical",
  "Australasia"
)

## Supertribe colours (NO Cleomaceae, NO Unplaced)
supertribe_cols2 <- c(
  "Arabodae (IV)"     = "#806CBD",
  "Brassicodae (II)"  = "#37ab76",
  "Camelinodae (I)"   = "#c65ead",
  "Heliophilodae (V)" = "#d3571e",
  "Hesperodae (III)"  = "#61b5e8",
  "Aethionemeae"      = "#ffed57"
)
supertribe_levels <- names(supertribe_cols2)

## Packages
req <- c(
  "dplyr", "tidyr", "stringr", "sf", "ggplot2",
  "readr", "rWCVP",
  "patchwork", "cowplot", "scales", "tibble"
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
  library(tidyr)
  library(stringr)
  library(sf)
  library(ggplot2)
  library(readr)
  library(rWCVP)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(tibble)
})

message("[INFO] Output dir: ", normalizePath(out_dir))



# Continent recoding helper

recode_L1_broad <- function(x) {
  x0 <- stringr::str_squish(as.character(x))
  xu <- stringr::str_to_upper(x0)
  dplyr::case_when(
    xu %in% c("ASIA-TEMPERATE", "ASIA-TROPICAL") ~ "Asia",
    TRUE ~ x0
  )
}

## STEP 1: READ + CLEAN SPECIES DATA -----

species_details <- readr::read_csv(
  species_file,
  show_col_types = FALSE
)

needed_cols <- c(
  "ACCEPTED", "FAMILY", "SPECIES_NAME_PRINT", "SUBSP_VAR",
  "SUPERTRIBE", "TRIBE_FULL",
  "WCVP_WGSRPD_LEVEL_1_native",
  "WCVP_WGSRPD_LEVEL_2_native",
  "WCVP_WGSRPD_LEVEL_3_native",
  "GROWTH_FORM"
)

miss_cols <- setdiff(needed_cols, colnames(species_details))

if (length(miss_cols) > 0) {
  stop(
    "[ERROR] Missing columns in species input file: ",
    paste(miss_cols, collapse = ", ")
  )
}

normalize_growth <- function(x) {
  x <- str_squish(as.character(x))
  x <- dplyr::na_if(x, "")
  x <- str_to_upper(x)
  dplyr::case_when(
    x %in% c("W", "W/L", "W/T") ~ "W",
    x == "H"                    ~ "H",
    TRUE                        ~ NA_character_
  )
}

is_emptyish <- function(x) {
  x <- as.character(x)
  x <- dplyr::na_if(str_squish(x), "")
  is.na(x)
}

species_distribution <- species_details %>%
  transmute(
    ACCEPTED = as.character(ACCEPTED),
    FAMILY = str_squish(as.character(FAMILY)),
    SPECIES_NAME_PRINT = str_squish(as.character(SPECIES_NAME_PRINT)),
    SUBSP_VAR = as.character(SUBSP_VAR),
    SUPERTRIBE = as.character(SUPERTRIBE),
    TRIBE_FULL = as.character(TRIBE_FULL),
    l1_native = as.character(WCVP_WGSRPD_LEVEL_1_native),
    l2_native = as.character(WCVP_WGSRPD_LEVEL_2_native),
    l3_native = as.character(WCVP_WGSRPD_LEVEL_3_native),
    GROWTH_FORM = normalize_growth(GROWTH_FORM)
  ) %>%
  filter(ACCEPTED == "Y") %>%
  filter(FAMILY == "Brassicaceae") %>%
  filter(is_emptyish(SUBSP_VAR)) %>%
  filter(!(is.na(l1_native) & is.na(l2_native) & is.na(l3_native))) %>%
  filter(!(str_to_lower(str_squish(TRIBE_FULL)) %in% c("unplaced"))) %>%
  mutate(
    woody = (GROWTH_FORM == "W"),
    species = SPECIES_NAME_PRINT,
    supertribe6 = case_when(
      str_to_lower(str_squish(TRIBE_FULL)) == "aethionemeae" ~ "Aethionemeae",
      TRUE ~ str_squish(as.character(SUPERTRIBE))
    )
  ) %>%
  filter(supertribe6 %in% supertribe_levels)

message(
  "[INFO] species_distribution: n=", nrow(species_distribution),
  " | woody=", sum(species_distribution$woody, na.rm = TRUE)
)


## STEP 2: EXPAND NATIVE WGSRPD CODES (L3) -----

split_codes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  str_split(x, "\\s*,\\s*", simplify = FALSE)
}

species_distribution_l3 <- species_distribution %>%
  mutate(l3_list = split_codes(l3_native)) %>%
  tidyr::unnest(l3_list) %>%
  rename(l3_code = l3_list) %>%
  mutate(l3_code = str_squish(as.character(l3_code))) %>%
  filter(l3_code != "")

## STEP 3: SUMMARISE COUNTS (L3) -----

l3_counts <- species_distribution_l3 %>%
  group_by(l3_code) %>%
  summarise(
    n_total = n_distinct(species),
    n_woody = n_distinct(species[woody]),
    .groups = "drop"
  ) %>%
  mutate(
    woody_fraction = dplyr::case_when(
      is.na(n_total) | n_total == 0 ~ NA_real_,
      TRUE ~ n_woody / n_total
    )
  )

## STEP 4: READ SHAPEFILE + DATELINE FIX -----

wgsrpd_L3 <- sf::st_read(shp_L3_path, quiet = TRUE)

wrap_dateline_safe <- function(x) {
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326 else x <- sf::st_transform(x, 4326)
  x <- suppressWarnings(sf::st_make_valid(x))
  if ("st_wrap_dateline" %in% getNamespaceExports("sf")) {
    return(suppressWarnings(sf::st_wrap_dateline(
      x, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE
    )))
  }
  if (requireNamespace("lwgeom", quietly = TRUE)) {
    return(suppressWarnings(lwgeom::st_wrap_dateline(
      x, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
    )))
  }
  warning("[WARN] Could not wrap dateline (sf::st_wrap_dateline not available and lwgeom not installed).")
  x
}

wgsrpd_L3 <- wrap_dateline_safe(wgsrpd_L3)

map_L3 <- wgsrpd_L3 %>% left_join(l3_counts, by = c("LEVEL3_COD" = "l3_code"))

map_L3 <- map_L3 %>%
  mutate(
    L1_name = rWCVP::wgsrpd_mapping$LEVEL1_NAM[match(LEVEL1_COD, rWCVP::wgsrpd_mapping$LEVEL1_COD)]
  )

map_L3 <- map_L3 %>%
  mutate(L1_broad = recode_L1_broad(L1_name))

if (length(DROP_L1_NAMES_FROM_MAP) > 0) {
  map_L3 <- map_L3 %>% filter(!L1_name %in% DROP_L1_NAMES_FROM_MAP)
}

# Drop non-polygon geoms (prevents “fat dots” that are POINT/MULTIPOINT artefacts)
if (DROP_NONPOLYGON_GEOMS) {
  gt <- as.character(sf::st_geometry_type(map_L3))
  keep <- gt %in% c("POLYGON","MULTIPOLYGON")
  map_L3 <- map_L3[keep, ]
}

# L1 contours by dissolving L3
l1_contours <- NULL
if (ADD_L1_CONTOURS) {
  l1_contours <- map_L3 %>%
    filter(!L1_name %in% L1_CONTOUR_EXCLUDE_NAMES) %>%
    group_by(L1_name) %>%
    summarise(geometry = suppressWarnings(sf::st_union(geometry)), .groups = "drop") %>%
    suppressWarnings(sf::st_make_valid())
}

make_globe_outline <- function(n = 720) {
  lon <- seq(-180, 180, length.out = n)
  lat_top <- rep(90, n)
  lat_bot <- rep(-90, n)
  top <- cbind(lon, lat_top)
  right <- cbind(rep(180, n), seq(90, -90, length.out = n))
  bot <- cbind(rev(lon), lat_bot)
  left <- cbind(rep(-180, n), seq(-90, 90, length.out = n))
  ring <- rbind(top, right, bot, left, top[1, , drop = FALSE])
  poly <- sf::st_polygon(list(ring))
  sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
}
globe_outline <- if (ADD_GLOBE_OUTLINE) make_globe_outline() else NULL

## STEP 5: BUILD MAP PANEL -----

# Projection for plotting + area-based tiny-geometry filtering
MAP_CRS <- "+proj=natearth"

map_L3 <- map_L3 %>%
  mutate(
    n_total_plot = n_total,
    n_woody_plot = n_woody,
    
    woody_fraction_plot = dplyr::case_when(
      is.na(n_total_plot) | n_total_plot == 0 ~ NA_real_,
      TRUE ~ n_woody_plot / n_total_plot
    ),
    
    fill_value_raw = dplyr::case_when(
      SHADING_MODE == "woody_richness" ~ as.numeric(n_woody_plot),
      SHADING_MODE == "woody_fraction" ~ as.numeric(woody_fraction_plot),
      TRUE ~ as.numeric(n_woody_plot)
    ),
    
    fill_value = if (SHADING_MODE == "woody_fraction" && FRACTION_SQUISH) {
      pmin(fill_value_raw, FRACTION_MAX)
    } else {
      fill_value_raw
    }
  )

# Rebuild dissolved L1 contours after filtering, so they match the plotted map
USE_AUTOMATIC_CONTINENT_BORDERS <- FALSE

manual_borders <- sf::st_sfc(
  # Northern America / Southern America: Mexico–Central America
  sf::st_linestring(matrix(c(
    -92.3, 14.5,
    -91.2, 17.2,
    -89.1, 17.8,
    -88.2, 18.5
  ), ncol = 2, byrow = TRUE)),
  
  # Europe / Asia: Greece–Turkey / Aegean transition
  sf::st_linestring(matrix(c(
    26.0, 41.7,
    26.4, 40.8,
    26.8, 39.8,
    27.3, 38.8,
    27.8, 37.8,
    28.3, 36.8
  ), ncol = 2, byrow = TRUE)),
  
  # Europe / Africa: Strait of Gibraltar
  sf::st_linestring(matrix(c(
    -6.0, 36.2,
    -5.6, 35.9,
    -5.2, 35.7
  ), ncol = 2, byrow = TRUE)),
  
  # Asia / Australasia: Wallacea–New Guinea / Australia transition
  sf::st_linestring(matrix(c(
    118.0, -8.5,
    123.0, -10.0,
    128.0, -10.5,
    134.0, -9.5,
    141.0, -9.0
  ), ncol = 2, byrow = TRUE)),
  
  # Europe / Asia: Urals + Caucasus approximation
  sf::st_linestring(matrix(c(
    60.0, 67.0,
    59.0, 62.0,
    58.0, 57.0,
    56.0, 52.0,
    53.0, 48.0,
    50.0, 45.5,
    46.5, 42.2,
    44.0, 42.8,
    41.5, 43.4,
    39.5, 43.6
  ), ncol = 2, byrow = TRUE)),
  
  # Africa / Asia: Sinai–Levant / Red Sea transition
  sf::st_linestring(matrix(c(
    32.3, 31.2,
    34.2, 29.5,
    35.0, 28.0,
    36.0, 25.0,
    38.0, 21.0,
    41.0, 15.0,
    43.0, 12.5
  ), ncol = 2, byrow = TRUE)),
  
  crs = 4326
)

continent_internal_borders <- sf::st_sf(
  border = c(
    "Northern America / Southern America",
    "Europe / Asia",
    "Africa / Asia",
    "Europe / Asia",
    "Europe / Africa",
    "Asia / Australasia"
  ),
  geometry = manual_borders
)

p_map_base <- ggplot() +
  geom_sf(data = map_L3, aes(fill = fill_value), color = NA)

if (ADD_L1_CONTOURS && !is.null(continent_internal_borders)) {
  p_map_base <- p_map_base +
    geom_sf(
      data = continent_internal_borders,
      color = "black",
      linewidth = 0.55,
      alpha = 1
    )
}

if (ADD_GLOBE_OUTLINE && !is.null(globe_outline)) {
  p_map_base <- p_map_base +
    geom_sf(
      data = globe_outline,
      fill = NA,
      color = GLOBE_OUTLINE_COL,
      linewidth = GLOBE_OUTLINE_LWD,
      alpha = GLOBE_OUTLINE_ALPHA
    )
}

if (SHADING_MODE == "woody_fraction") {
  brks <- seq(0, FRACTION_MAX, by = 0.05)
  lbls <- scales::label_percent(accuracy = 1)(brks)
  if (length(lbls) > 0) {
    lbls[length(lbls)] <- paste0(
      scales::label_percent(accuracy = 1)(FRACTION_MAX),
      " (or more)"
    )
  }
  
  p_map_base <- p_map_base +
    scale_fill_gradient(
      name = "Woody fraction",
      low = "#eadfe6",
      high = "#8a2154",
      na.value = "grey97",
      limits = c(0, FRACTION_MAX),
      oob = scales::squish,
      breaks = brks,
      labels = lbls
    )
} else {
  p_map_base <- p_map_base +
    scale_fill_gradient(
      name = "Woody species\n(count)",
      low = "#eadfe6",
      high = "#8a2154",
      na.value = "grey97"
    )
}

# Natural Earth projection
p_map_base <- p_map_base + coord_sf(crs = MAP_CRS)

# Extract the map legend and place it INSIDE the map (so waffles can use full width)
map_legend_grob <- cowplot::get_legend(
  p_map_base +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
)

p_map <- cowplot::ggdraw() +
  cowplot::draw_plot(
    p_map_base +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        legend.position = "none",
        plot.margin = margin(6, 6, 6, 6)
      ) +
      labs(title = NULL, subtitle = NULL, caption = NULL),
    x = 0, y = 0, width = 1, height = 1
  ) +
  cowplot::draw_grob(map_legend_grob, x = 0.80, y = 0.32, width = 0.18, height = 0.38)

## STEP 6: WAFFLES (L1; single row below map) -----

# 6.0 Ensure palette + levels match what we want (no Cleomaceae, no Unplaced)
# (This assumes you already changed supertribe_cols2 earlier; this is just a safety check.)
supertribe_levels <- names(supertribe_cols2)

# 6.1 L1 numeric code -> name mapping (from your WGSRPD explanation)
l1_name_map <- c(
  "1"="Europe","2"="Africa","3"="Asia-Temperate","4"="Asia-Tropical",
  "5"="Australasia","6"="Pacific","7"="Northern America","8"="Southern America","9"="Antarctic"
)

# Expand L1 codes from the SHEET column (WCVP_WGSRPD_LEVEL_1_native)
species_distribution_l1 <- species_distribution %>%
  mutate(l1_list = split_codes(l1_native)) %>%
  tidyr::unnest(l1_list) %>%
  rename(l1_code = l1_list) %>%
  mutate(l1_code = stringr::str_squish(as.character(l1_code))) %>%
  filter(l1_code != "") %>%
  mutate(l1_name = unname(l1_name_map[as.character(l1_code)])) %>%
  filter(!is.na(l1_name)) %>%
  filter(l1_name != "Antarctic") # no Antarctic waffle

species_distribution_l1 <- species_distribution_l1 %>%
  mutate(l1_name = recode_L1_broad(l1_name))

message("[INFO] L1 names present (from sheet): ",
        paste(sort(unique(species_distribution_l1$l1_name)), collapse = ", "))

# Totals per L1 (for waffle sizing)
l1_totals <- species_distribution_l1 %>%
  group_by(l1_name) %>%
  summarise(
    n_total = n_distinct(species),
    n_woody = n_distinct(species[woody]),
    .groups = "drop"
  )

# Supertribe composition per L1 (all species)
l1_supertribe_counts <- species_distribution_l1 %>%
  group_by(l1_name, supertribe6) %>%
  summarise(n = n_distinct(species), .groups = "drop") %>%
  # Safety: drop anything outside our desired supertribes
  filter(supertribe6 %in% supertribe_levels)

# Supertribe composition per L1 (woody only)
l1_supertribe_counts_woody <- species_distribution_l1 %>%
  filter(woody) %>%
  group_by(l1_name, supertribe6) %>%
  summarise(n_woody = n_distinct(species), .groups = "drop") %>%
  filter(supertribe6 %in% supertribe_levels)

# Requested waffle order (L1 names)
WAFFLE_L1_ORDER_NAMES <- c(
  "Pacific",
  "Northern America",
  "Southern America",
  "Europe",
  "Africa",
  "Asia",
  "Australasia"
)

present_l1 <- sort(unique(l1_totals$l1_name))
l1_show <- WAFFLE_L1_ORDER_NAMES[WAFFLE_L1_ORDER_NAMES %in% present_l1]
if (length(l1_show) == 0) {
  warning("[WARN] None of requested L1 names found; using present L1 names instead.")
  l1_show <- present_l1
}


# Shared grid size across shown L1 names (reuse helper)

compute_shared_grid_tiles <- function(counts_df, l1_names, tile_value_spp = 20, ncol = 5) {
  dat <- counts_df %>%
    filter(l1_name %in% l1_names) %>%
    group_by(l1_name, supertribe6) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(n_tiles = ceiling(n / tile_value_spp)) %>%
    group_by(l1_name) %>%
    summarise(total_tiles = sum(n_tiles), .groups = "drop")
  
  if (nrow(dat) == 0) return(ncol)
  
  ceiling(max(dat$total_tiles, na.rm = TRUE) / ncol) * ncol
}

shared_grid_tiles <- compute_shared_grid_tiles(
  counts_df        = l1_supertribe_counts,
  l1_names         = l1_show,
  tile_value_spp   = WAFFLE_TILE_VALUE_SPP,
  ncol             = WAFFLE_NCOL
)

# Allocate tiles in FIXED order (supertribe_levels), returning full-length named integer vector
allocate_tiles_fixed_order <- function(counts_named, n_tiles_total, levels_fixed) {
  v <- setNames(rep(0, length(levels_fixed)), levels_fixed)
  if (!is.null(counts_named) && length(counts_named) > 0) {
    common <- intersect(names(counts_named), names(v))
    v[common] <- counts_named[common]
  }
  if (sum(v) <= 0 || n_tiles_total <= 0) return(as.integer(v))
  p <- v / sum(v)
  tiles_raw <- p * n_tiles_total
  tiles <- floor(tiles_raw)
  rem <- tiles_raw - tiles
  missing <- n_tiles_total - sum(tiles)
  if (missing > 0) {
    add_idx <- order(rem, decreasing = TRUE)
    add_idx <- add_idx[seq_len(min(missing, length(add_idx)))]
    tiles[add_idx] <- tiles[add_idx] + 1
  }
  tiles <- as.integer(tiles)
  names(tiles) <- levels_fixed
  tiles
}

# Build waffle plot for one L1 (NA-safe woody outlines)
make_waffle_L1 <- function(l1_name,
                           counts_df, woody_counts_df, totals_df,
                           grid_tiles, tile_value_spp = 20, ncol = 5,
                           supertribe_levels, supertribe_cols2) {
  
  tot <- totals_df %>% dplyr::filter(L1_name == !!l1_name)
  if (nrow(tot) == 0 || is.na(tot$n_total[1]) || tot$n_total[1] == 0) {
    return(ggplot() + theme_void())
  }
  
  dat_all <- counts_df %>%
    dplyr::filter(L1_name == !!l1_name) %>%
    tidyr::complete(supertribe6 = supertribe_levels, fill = list(n = 0)) %>%
    dplyr::mutate(supertribe6 = factor(supertribe6, levels = supertribe_levels)) %>%
    dplyr::arrange(supertribe6)
  
  dat_woody <- woody_counts_df %>%
    dplyr::filter(L1_name == !!l1_name) %>%
    tidyr::complete(supertribe6 = supertribe_levels, fill = list(n_woody = 0)) %>%
    dplyr::mutate(supertribe6 = factor(supertribe6, levels = supertribe_levels)) %>%
    dplyr::arrange(supertribe6)
  
  tile_list <- list()
  k <- 1L
  
  for (ii in seq_len(nrow(dat_all))) {
    nm <- as.character(dat_all$supertribe6[ii])
    n_sp <- dat_all$n[ii]
    if (is.na(n_sp) || n_sp <= 0) next
    
    if (isTRUE(WAFFLE_USE_FRACTIONAL_TILES)) {
      n_full <- floor(n_sp / tile_value_spp)
      rem <- n_sp %% tile_value_spp
      
      if (n_full > 0) {
        tile_list[[k]] <- tibble::tibble(
          supertribe6 = rep(nm, n_full),
          tile_fraction = rep(1, n_full)
        )
        k <- k + 1L
      }
      
      if (rem > 0) {
        tile_list[[k]] <- tibble::tibble(
          supertribe6 = nm,
          tile_fraction = rem / tile_value_spp
        )
        k <- k + 1L
      }
      
    } else {
      n_tiles <- switch(
        WAFFLE_ROUNDING_MODE,
        ceiling = ceiling(n_sp / tile_value_spp),
        round   = round(n_sp / tile_value_spp),
        floor   = floor(n_sp / tile_value_spp),
        ceiling(n_sp / tile_value_spp)
      )
      
      if (n_sp > 0) n_tiles <- max(1L, as.integer(n_tiles))
      
      tile_list[[k]] <- tibble::tibble(
        supertribe6 = rep(nm, n_tiles),
        tile_fraction = rep(1, n_tiles)
      )
      k <- k + 1L
    }
  }
  
  used_tiles <- dplyr::bind_rows(tile_list)
  n_used <- nrow(used_tiles)
  
  if (n_used > grid_tiles) {
    warning(
      "[WARN] Waffle for ", l1_name, " needs ", n_used,
      " tiles but grid has only ", grid_tiles,
      ". Increase shared_grid_tiles or WAFFLE_NCOL."
    )
    used_tiles <- used_tiles[seq_len(grid_tiles), , drop = FALSE]
    n_used <- grid_tiles
  }
  
  waffle_df <- tibble::tibble(
    tile = seq_len(grid_tiles),
    x = (tile - 1) %% ncol + 1,
    y = (tile - 1) %/% ncol + 1,
    supertribe6 = NA_character_,
    tile_fraction = 1
  )
  
  if (n_used > 0) {
    waffle_df$supertribe6[seq_len(n_used)] <- used_tiles$supertribe6
    waffle_df$tile_fraction[seq_len(n_used)] <- used_tiles$tile_fraction
  }
  
  waffle_df <- waffle_df %>%
    dplyr::mutate(
      supertribe6 = factor(supertribe6, levels = supertribe_levels),
      tile_width  = sqrt(tile_fraction),
      tile_height = sqrt(tile_fraction),
      x_plot = x - (1 - tile_width) / 2,
      y_plot = y - (1 - tile_height) / 2
    )
  
  woody_counts_aligned <- setNames(rep(0, length(supertribe_levels)), supertribe_levels)
  if (nrow(dat_woody) > 0) {
    woody_counts_aligned[as.character(dat_woody$supertribe6)] <- dat_woody$n_woody
  }
  
  woody_df <- NULL
  if (WAFFLE_WOODY_MARK %in% c("dot", "outline", "overlay")) {
    idx <- integer(0)
    
    for (nm in supertribe_levels) {
      if (woody_counts_aligned[[nm]] > 0) {
        candidate <- which(as.character(waffle_df$supertribe6) == nm)
        if (length(candidate) > 0) idx <- c(idx, candidate[1])
      }
    }
    
    if (length(idx) > 0) {
      woody_df <- waffle_df[idx, , drop = FALSE]
    }
  }
  
  woody_pct <- 100 * tot$n_woody[1] / tot$n_total[1]
  
  waffle_title <- sprintf(
    "%s\n(n = %s; %.1f%% woody)",
    l1_name,
    scales::comma(tot$n_total[1]),
    woody_pct
  )
  
  gg <- ggplot(waffle_df, aes(x, y)) +
    geom_tile(
      fill = WAFFLE_EMPTY_FILL,
      color = "white",
      linewidth = WAFFLE_TILE_GAP_LINEWIDTH
    ) +
    geom_tile(
      aes(
        x = x_plot,
        y = y_plot,
        fill = supertribe6,
        width = tile_width,
        height = tile_height
      ),
      color = "white",
      linewidth = WAFFLE_TILE_GAP_LINEWIDTH
    ) +
    coord_equal(expand = FALSE) +
    scale_y_reverse() +
    scale_fill_manual(
      values = supertribe_cols2,
      drop = FALSE,
      na.value = WAFFLE_EMPTY_FILL
    ) +
    labs(title = waffle_title) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(
        size = WAFFLE_TITLE_SIZE,
        face = "bold",
        margin = margin(b = 2)
      ),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
  
  if (!is.null(woody_df) && WAFFLE_WOODY_MARK == "dot") {
    gg <- gg +
      geom_point(
        data = woody_df,
        aes(x = x_plot, y = y_plot),
        inherit.aes = FALSE,
        shape = 16,
        size = WAFFLE_WOODY_DOT_SIZE,
        alpha = WAFFLE_WOODY_DOT_ALPHA,
        color = "black"
      )
  }
  
  if (!is.null(woody_df) && WAFFLE_WOODY_MARK == "outline") {
    gg <- gg +
      geom_tile(
        data = woody_df,
        aes(
          x = x_plot,
          y = y_plot,
          width = tile_width * WAFFLE_WOODY_OUTLINE_INSET,
          height = tile_height * WAFFLE_WOODY_OUTLINE_INSET
        ),
        inherit.aes = FALSE,
        fill = NA,
        color = WAFFLE_WOODY_OUTLINE_COL,
        linewidth = WAFFLE_WOODY_OUTLINE_LWD
      )
  }
  
  gg
}

# Build waffles
waffle_plots <- lapply(
  l1_show,
  make_waffle_L1,
  counts_df         = l1_supertribe_counts %>% rename(L1_name = l1_name),
  woody_counts_df   = l1_supertribe_counts_woody %>% rename(L1_name = l1_name),
  totals_df         = l1_totals %>% rename(L1_name = l1_name),
  grid_tiles        = shared_grid_tiles,
  tile_value_spp    = WAFFLE_TILE_VALUE_SPP,
  ncol              = WAFFLE_NCOL,
  supertribe_levels = supertribe_levels,
  supertribe_cols2  = supertribe_cols2
)

stopifnot(all(vapply(waffle_plots, inherits, logical(1), "ggplot")))

# Legend column for waffles (title required)
legend_dummy <- ggplot(
  tibble::tibble(
    x = seq_along(supertribe_levels),
    y = 1,
    supertribe6 = factor(supertribe_levels, levels = supertribe_levels)
  ),
  aes(x, y, fill = supertribe6)
) +
  geom_tile() +
  scale_fill_manual(values = supertribe_cols2, drop = FALSE) +
  guides(fill = guide_legend(title = "Main lineages", ncol = 1, byrow = TRUE)) +
  theme_void(base_size = 9) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

legend_grob <- cowplot::get_legend(legend_dummy)
p_legend <- cowplot::ggdraw(legend_grob)

# Combine waffles + legend as one row (full-width row because map legend is inside map panel)
plots_row <- c(waffle_plots, list(p_legend))
p_waffle_row <- patchwork::wrap_plots(plots_row, nrow = 1, ncol = length(plots_row))


## STEP 7: ASSEMBLE FINAL FIGURE -----

Figure1 <- p_map
if (ADD_L1_WAFFLES && exists("p_waffle_row") && !is.null(p_waffle_row)) {
  Figure1 <- p_map / p_waffle_row +
    patchwork::plot_layout(heights = c(3.35, 1.15))
}

## STEP 8: SAVE -----

ggsave(out_pdf, plot = Figure1, width = 12.5, height = 7.6, units = "in", device = cairo_pdf)
ggsave(out_png, plot = Figure1, width = 12.5, height = 7.6, units = "in", dpi = 300)

message("[DONE] Saved:\n  ", out_pdf, "\n  ", out_png)


## STEP 9: HOTSPOT ANALYSIS OF WOODINESS (SES + P VALUES) -----

# Goal:
# Identify WGSRPD Level-3 botanical countries with significantly
# more or fewer woody Brassicaceae species than expected given
# their total Brassicaceae richness.

# SETTINGS

HOTSPOT_N_REPLICATES <- 1000
HOTSPOT_SEED <- 123

HOTSPOT_MIN_SPECIES <- 10   # avoid unstable SES in tiny floras

HOTSPOT_OUT_PREFIX <- file.path(
  out_dir,
  "Figure_SX_woodiness_hotspots"
)

set.seed(HOTSPOT_SEED)


# PREPARE SPECIES × L3 MATRIX


message("[INFO] Preparing hotspot analysis input...")

hotspot_df <- species_distribution_l3 %>%
  distinct(species, woody, l3_code)

# observed counts per botanical country
obs_tbl <- hotspot_df %>%
  group_by(l3_code) %>%
  summarise(
    n_total = n_distinct(species),
    n_woody_obs = n_distinct(species[woody]),
    woody_fraction_obs = n_woody_obs / n_total,
    .groups = "drop"
  ) %>%
  filter(n_total >= HOTSPOT_MIN_SPECIES)


# NULL MODEL


message("[INFO] Running null model randomisations (n = ",
        HOTSPOT_N_REPLICATES, ")...")

species_pool <- hotspot_df %>%
  distinct(species, woody)

n_woody_global <- sum(species_pool$woody, na.rm = TRUE)

species_ids <- species_pool$species

null_mat <- matrix(
  NA_real_,
  nrow = nrow(obs_tbl),
  ncol = HOTSPOT_N_REPLICATES
)

rownames(null_mat) <- obs_tbl$l3_code

# incidence list speeds things up
l3_species_list <- hotspot_df %>%
  group_by(l3_code) %>%
  summarise(species_vec = list(unique(species)), .groups = "drop")

l3_species_lookup <- setNames(
  l3_species_list$species_vec,
  l3_species_list$l3_code
)

for (rr in seq_len(HOTSPOT_N_REPLICATES)) {
  
  if (rr %% 100 == 0) {
    message("[INFO] Randomisation ", rr, "/", HOTSPOT_N_REPLICATES)
  }
  
  woody_random_species <- sample(
    species_ids,
    size = n_woody_global,
    replace = FALSE
  )
  
  woody_random_species <- unique(woody_random_species)
  
  for (ii in seq_len(nrow(obs_tbl))) {
    
    spp_vec <- l3_species_lookup[[obs_tbl$l3_code[ii]]]
    
    null_mat[ii, rr] <- sum(spp_vec %in% woody_random_species)
  }
}


# SES + P VALUES


message("[INFO] Calculating SES statistics...")

null_mean <- rowMeans(null_mat, na.rm = TRUE)
null_sd <- apply(null_mat, 1, sd, na.rm = TRUE)

ses <- (obs_tbl$n_woody_obs - null_mean) / null_sd

p_upper <- vapply(
  seq_len(nrow(obs_tbl)),
  function(i) {
    mean(null_mat[i, ] >= obs_tbl$n_woody_obs[i], na.rm = TRUE)
  },
  numeric(1)
)

p_lower <- vapply(
  seq_len(nrow(obs_tbl)),
  function(i) {
    mean(null_mat[i, ] <= obs_tbl$n_woody_obs[i], na.rm = TRUE)
  },
  numeric(1)
)

hotspot_results <- obs_tbl %>%
  mutate(
    null_mean = null_mean,
    null_sd = null_sd,
    SES = ses,
    p_upper = p_upper,
    p_lower = p_lower,
    
    hotspot_class = case_when(
      p_upper <= 0.05 & SES > 0 ~ "More woody than expected",
      p_lower <= 0.05 & SES < 0 ~ "Less woody than expected",
      TRUE ~ "Not significant"
    )
  )


# JOIN TO MAP


message("[INFO] Joining hotspot statistics to map...")

map_hotspot <- map_L3 %>%
  left_join(
    hotspot_results,
    by = c("LEVEL3_COD" = "l3_code")
  )

sig_pts <- map_hotspot %>%
  filter(
    !is.na(SES),
    hotspot_class %in% c("More woody than expected", "Less woody than expected")
  ) %>%
  suppressWarnings(sf::st_point_on_surface()) %>%
  mutate(sig_label = "*")


# HOTSPOT MAP


p_hotspot <- ggplot() +
  
  geom_sf(
    data = map_hotspot,
    aes(fill = SES),
    color = NA
  ) +
  
  geom_sf_text(
    data = sig_pts,
    aes(label = sig_label),
    size = 2.6,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_sf(
    data = continent_internal_borders,
    color = "black",
    linewidth = 0.55
  ) +
  
  geom_sf(
    data = globe_outline,
    fill = NA,
    color = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_gradient2(
    name = "SES\n(woodiness)",
    low = "#2b6cb0",
    mid = "grey97",
    high = "#b2182b",
    midpoint = 0,
    na.value = "grey92",
    limits = c(-4, 4),
    oob = scales::squish
  ) +
  
  coord_sf(crs = MAP_CRS) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  
  labs(
    title = "Global hotspots of derived woodiness in Brassicaceae",
    subtitle =
      paste0(
        "SES from null-model randomisations (n = ",
        HOTSPOT_N_REPLICATES,
        "); positive values indicate more woody species than expected"
      )
  )


# SAVE


hotspot_pdf <- paste0(HOTSPOT_OUT_PREFIX, ".pdf")
hotspot_png <- paste0(HOTSPOT_OUT_PREFIX, ".png")
hotspot_tsv <- paste0(HOTSPOT_OUT_PREFIX, ".tsv")

ggsave(
  hotspot_pdf,
  p_hotspot,
  width = 12.5,
  height = 6.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  hotspot_png,
  p_hotspot,
  width = 12.5,
  height = 6.8,
  units = "in",
  dpi = 300
)

readr::write_tsv(hotspot_results, hotspot_tsv)

message("[DONE] Hotspot outputs saved:")
message("  ", hotspot_pdf)
message("  ", hotspot_png)
message("  ", hotspot_tsv)

