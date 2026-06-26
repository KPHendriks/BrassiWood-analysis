###############################################################################
# PUBLICATION NOTE
#
# This script is included as supplementary source code to document the
# methodology used to construct the Brassicaceae occurrence-record curation
# application used in this study.
#
# To prevent accidental deployment or modification of the live curation
# platform, all functionality requiring authentication, external databases,
# or online write access has been removed or disabled in this public version.
#
# Consequently, this script serves as a transparent description of the
# workflow and application architecture, but is not intended to function as
# a deployable Shiny application.
###############################################################################


## STEP 0: PREPARE ----

# Set working directory
setwd("~/Google Drive/My Drive/Publications/2026_Hendriks_et_al_BrassiWood/BrassiWood")

# Load necessary libraries
library(shiny)
library(leaflet)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(sf)
library(rsconnect)
library(geosphere)


# Create Shiny app directory if not yet present
if (!dir.exists("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app")){
  dir.create("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app", recursive = TRUE, showWarnings = FALSE)
}

# Create a directory inside the Shiny app directory to contain the occurrence tables for all species
if (!dir.exists("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data")){
  dir.create("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data", recursive = TRUE, showWarnings = FALSE)
}

# Create a directory inside the Shiny app directory to contain a favicon
if (!dir.exists("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/www")){
  dir.create("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/www", recursive = TRUE, showWarnings = FALSE)
}



## STEP 1: PREPARE TAXONOMIC DATA FOR APP ----

# Load species details table from frozen publication input
species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
stopifnot(file.exists(species_file))

species_details_cleaned <- readr::read_csv(
  species_file,
  show_col_types = FALSE
) %>%
  dplyr::filter(ACCEPTED == "Y") %>%
  dplyr::filter(is.na(SUBSP_VAR) | stringr::str_squish(as.character(SUBSP_VAR)) == "")

# Prepare and export taxonomy table
taxonomy_export <- data.frame(
  subfamily = species_details_cleaned$SUBFAMILY,
  supertribe = species_details_cleaned$SUPERTRIBE,
  tribe = species_details_cleaned$TRIBE_FULL,
  genus = species_details_cleaned$GENUS,
  species = gsub(" ", "_", species_details_cleaned$SPECIES_NAME_PRINT),
  species_print = species_details_cleaned$SPECIES_NAME_PRINT,
  native_range = species_details_cleaned$WCVP_WGSRPD_LEVEL_3_native
)

write.csv(taxonomy_export, "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/taxonomy_table.csv", row.names = FALSE)



## STEP 2: PREPARE SHAPEFILES FOR APP ----

# Copy all level3 shapefiles (shp, dbf, shx, prj, etc.)
if (!dir.exists("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/level3_maps")){
  dir.create("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/level3_maps", recursive = TRUE, showWarnings = FALSE)
}
file.copy(
  from = list.files("WP2_BrassiNiche/data/wgsrpd-master/level3", pattern = "^level3\\.", full.names = TRUE),
  to = "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/level3_maps/",
  overwrite = TRUE
)



## STEP 3: PREPARE SPECIES OCCURRENCE DATA FOR APP ----

# In script 2b we collected species occurrence data, but we only set a target minimum number of records.
# We will now use these data to create the curation app, but will limit the number of records to 100 per species (if more records are available).
# This way, we can keep the expert curation manageable.

# Below, we will apply the following strategy (species by species):
# 1. Subset the data to only those records having passed all tests;
# 2. Based on the location coordinates, assign geographic cluster to each record, which will help keeping records from as much of the species distribution range as possible;
# 3. From each cluster, select the best occurrence records based on the lowest coordinate uncertainty.


# Set source and destination directories
source_dir <- "WP2_BrassiNiche/results_intermediate/species_occurrence_records_tables"
dest_dir   <- "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data"

# Get full list of files in source
files <- list.files(source_dir, full.names = TRUE)

# constants
k_default           <- 10             # target clusters when n > 100
target_per_cluster  <- 10             # used implicitly by round-robin target of 100 when n>100
max_records         <- 10000
coord_decimals      <- 4              # 6 ≈ 0.11 m, 5 ≈ 1.1 m, 4 ≈ 11 m

for (src_file in files) {
  
  file_name <- basename(src_file)
  species   <- tools::file_path_sans_ext(file_name)
  
  occurrence_data_all <- read.csv(src_file, header = TRUE)
  
  # 1) keep only records that passed all tests
  occurrence_data_valid <- occurrence_data_all[occurrence_data_all$passed_all_tests %in% TRUE, ]
  
  # 2) cap computational size if huge
  if (nrow(occurrence_data_valid) > max_records) {
    set.seed(42)
    occurrence_data_valid <- occurrence_data_valid[sample.int(nrow(occurrence_data_valid), max_records), ]
  }
  
  # 3) if trivially small, just export what we have
  if (nrow(occurrence_data_valid) <= 10) {
    occurrence_data_for_export <- occurrence_data_valid
  } else {
    
    # 4) drop rows without coordinates and add row_id
    df <- occurrence_data_valid %>%
      filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
      mutate(row_id = dplyr::row_number())
    
    if (nrow(df) <= 10) {
      # (after removing NA coords it might drop below threshold)
      occurrence_data_for_export <- df
    } else {
      n <- nrow(df)
      
      # --- choose number of clusters per your rules ---
      # if 11–100: clusters = ceiling(n/10); if >100: clusters = 10; if <=10: 1 (already handled)
      if (n > 100) {
        local_k <- k_default
      } else {
        local_k <- min(ceiling(n / 10), n)  # safety: k cannot exceed n
      }
      
      # --- clustering (only if local_k > 1) ---
      if (local_k > 1) {
        coords <- df %>% select(decimalLongitude, decimalLatitude)
        D  <- geosphere::distm(as.matrix(coords), fun = geosphere::distHaversine)
        hc <- hclust(as.dist(D), method = "ward.D2")
        df$cluster <- factor(cutree(hc, k = local_k))
      } else {
        df$cluster <- factor(1L)
      }
      
      # --- dedup "same location" within cluster by rounded coords ---
      df <- df %>%
        mutate(
          lon_key   = sprintf(paste0("%.", coord_decimals, "f"), decimalLongitude),
          lat_key   = sprintf(paste0("%.", coord_decimals, "f"), decimalLatitude),
          coord_key = paste0(lon_key, "_", lat_key)
        )
      
      # --- deterministic ranking by priorities (best first within cluster) ---
      ranked <- df %>%
        mutate(
          source_rank = dplyr::case_when(
            coordinate_source == "GBIF"      ~ 0L,
            coordinate_source == "geocoding" ~ 1L,
            TRUE                             ~ 2L
          ),
          uncertainty_rank = tidyr::replace_na(as.numeric(coordinateUncertaintyInMeters), Inf)
        ) %>%
        arrange(cluster, source_rank, uncertainty_rank, row_id)
      
      # --- de-duplicate exact (rounded) coordinates within each cluster ---
      dedup <- ranked %>%
        group_by(cluster) %>%
        distinct(coord_key, .keep_all = TRUE) %>%
        mutate(rank_within_cluster = dplyr::row_number()) %>%
        ungroup()
      
      # --- decide total target per your rule set ---
      # if original n > 100: target = min(100, available uniques after dedup)
      # else (<=100): keep all available uniques (balanced by round-robin)
      if (n > 100) {
        total_target <- min(100L, nrow(dedup))
      } else {
        total_target <- nrow(dedup)
      }
      
      # --- round-robin fill across clusters by 'rank_within_cluster' to maintain coverage ---
      kept <- dedup[0, ]
      r <- 1L
      while (nrow(kept) < total_target) {
        batch <- dedup %>% filter(rank_within_cluster == r)
        if (nrow(batch) == 0L) break
        # keep global priority within the same rank
        batch <- batch %>% arrange(cluster, source_rank, uncertainty_rank, row_id)
        slots_left <- total_target - nrow(kept)
        if (nrow(batch) > slots_left) batch <- dplyr::slice_head(batch, n = slots_left)
        kept <- bind_rows(kept, batch)
        r <- r + 1L
      }
      
      occurrence_data_for_export <- kept
    }
  }
  
  # write output
  dest_file <- file.path(dest_dir, paste0(species, ".csv"))
  write.csv(occurrence_data_for_export, file = dest_file, row.names = FALSE)
  
  # clean local objects (optional)
  rm(list = intersect(ls(), c(
    "coords","D","df","hc","kept","occurrence_data_for_export",
    "occurrence_data_all","occurrence_data_valid","ranked","dedup",
    "dest_file","file_name","species","local_k","k", "n","r","batch","total_target"
  )))
}



## STEP 4: WRITE SHINY APP SCRIPT ----

# Set working directory
setwd("~/Google Drive/My Drive/Publications/2026_Hendriks_et_al_BrassiWood/BrassiWood/WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/")

# Write app.R into Brassicaceae_taxonomic_curation_app folder
# Make sure to copy the rest of the script from STEP 3 into the file "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/app.R"

# app.R (icons reinstated + species & project dashboard) -- FIXED onclick bug
# Note: after publication, this Google Sheets access was removed from public/reviewer version.
# Ratings are read from a local static CSV if present.
# Google authentication and online rating storage have been removed
# from the public version of this application.
suppressPackageStartupMessages({
  library(shiny)
  library(leaflet)
  library(dplyr)
  library(readr)
  library(stringr)
  library(sf)
  library(tidyr)
  library(glue)
  library(shinymanager)
  library(lubridate)
  library(jsonlite)
})

# =========================
# ====== CONFIG/SETUP =====
# =========================
# Public/reviewer version: no live Google Sheet connection.
# Optional example ratings used only for demonstration.
RATINGS_FILE <- "ratings_example.csv"

# Demo password only; change or disable for any deployment.
SHARED_PW <- Sys.getenv("BRASSIWOOD_APP_PASSWORD", unset = "demo")

TIME_TZ <- "Europe/Amsterdam"
CONTACT_EM <- "kasper.hendriks@naturalis.nl"

# 👉 Priority lists: each as plain text, one species per line (Genus epithet with a space)
BRASSITOL_SPECIES <- tryCatch({
  x <- readLines("brassitol_species.txt", warn = FALSE, encoding = "UTF-8")
  x <- trimws(x); x[nzchar(x)]
}, error = function(e) character(0))

DERIVED_WOODY_SPECIES <- tryCatch({
  x <- readLines("derived_woody_species.txt", warn = FALSE, encoding = "UTF-8")
  x <- trimws(x); x[nzchar(x)]
}, error = function(e) character(0))

# 👉 Custom icons (served from /www). Provide either SVG or PNG; we autodetect.
detect_icon <- function(base){
  files <- list.files("www", pattern = paste0("^", base, "\\.(svg|png)$"), ignore.case = TRUE)
  if (length(files)) return(files[1])
  NULL
}
PHYLO_ICON_URL <- detect_icon("phylo")
WOODY_ICON_URL <- detect_icon("woody")
# default URLs, served from /www (NO leading slash so it works locally & on shinyapps)
PHYLO_ICON_URL <- if (!is.null(PHYLO_ICON_URL) && nzchar(PHYLO_ICON_URL)) { paste0(PHYLO_ICON_URL, "?v=7") } else { "phylo.svg?v=7" }
WOODY_ICON_URL <- if (!is.null(WOODY_ICON_URL) && nzchar(WOODY_ICON_URL)) { paste0(WOODY_ICON_URL, "?v=7") } else { "woody.svg?v=7" }

# Authenticate Google Sheets (requires service-account.json in app dir)
gs4_auth(
  path = "service-account.json",
  scopes = c("https://www.googleapis.com/auth/spreadsheets",
             "https://www.googleapis.com/auth/drive")
)
options(googlesheets4_quiet = FALSE, gargle_verbosity = "debug")
ss <- NULL
SS_OBJ <- tryCatch(gs4_get(SHEET_ID), error = function(e) NULL)
try({
  message("[GS] Auth identity: ", as.character(gs4_user()))
  if (!is.null(SS_OBJ)) message("[GS] Sheet OK. Tabs: ", paste(sheet_names(SS_OBJ), collapse = ", "))
}, silent = TRUE)

# -----------------
# Helper functions
# -----------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

natural_cluster_levels <- function(x) {
  ux <- unique(as.character(x[is.finite(x) | !is.na(x)]))
  suppressWarnings(num <- as.integer(ux))
  if (all(!is.na(num))) as.character(sort(num)) else sort(ux, method = "radix")
}

latest_ratings_for <- function(gbif_ids, ratings_df) {
  if (is.null(ratings_df) || nrow(ratings_df) == 0) {
    return(tibble::tibble(
      gbifID = as.character(gbif_ids),
      latest = NA_character_,
      last3 = I(vector("list", length(gbif_ids)))
    ))
  }
  r <- ratings_df |>
    mutate(timestamp = lubridate::ymd_hms(timestamp_iso, quiet = TRUE)) |>
    arrange(.by_group = FALSE, gbifID, desc(timestamp))
  
  latest <- r |>
    group_by(gbifID) |>
    summarise(
      latest = dplyr::first(rating),
      last3 = list(head(paste0(timestamp, " — ", user, ": ", rating), 3)),
      .groups = "drop"
    )
  
  tibble::tibble(gbifID = as.character(gbif_ids)) |>
    left_join(latest, by = "gbifID")
}

linkify <- function(x) {
  if (is.null(x)) return(x)
  x <- as.character(x); x[is.na(x)] <- ""
  stringr::str_replace_all(
    x,
    stringr::regex("(https?://[^\\s<>()]+)", ignore_case = TRUE),
    function(m) sprintf("<a href='%s' target='_blank'>%s</a>", m, m)
  )
}

REQUIRED_COLS <- c(
  "decimalLatitude","decimalLongitude","coordinate_source",
  "gbifID","SPECIES_NAME_PRINT","TRIBE_FULL",
  "higherGeography","countryCode","country_full_name",
  "county","locality","locationRemarks","stateProvince","municipality",
  "coordinateUncertaintyInMeters","elevation","elevationAccuracy",
  "coordinate_cleaner_.summary","geocoding_errors",
  "basisOfRecord","eventDate","recordedBy","identifiedBy",
  "recordNumber","datasetName","occurrenceID",
  "cluster"
)
ensure_required_cols <- function(df) {
  missing <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing)) for (col in missing) df[[col]] <- NA
  df
}

AM_COLORS <- c("orange","green","darkgreen","blue","darkblue","purple","cadetblue","pink","gray","yellow")
normalize_species <- function(s) { tolower(gsub("\\s+", " ", gsub("_", " ", trimws(as.character(s))))) }
pretty_species <- function(s) gsub("_", " ", as.character(s))

# ---------------------------
# Data: taxonomy + WGSRPD L3
# ---------------------------
taxonomy <- read_csv("taxonomy_table.csv", show_col_types = FALSE) |>
  arrange(tribe, species)

wgsrpd_L3 <- tryCatch({
  st_read("level3_maps/level3.shp", quiet = TRUE)
}, error = function(e) {
  message("[WARN] Could not read WGSRPD Level 3 shapefile: ", conditionMessage(e))
  NULL
})

all_species <- sort(unique(taxonomy$species))
species_gbif_index <- local({
  message("[INIT] Building species GBIF index (including zero-record species)...")
  files <- setNames(file.path("occurrence_data", paste0(all_species, ".csv")), all_species)
  idx <- lapply(names(files), function(sp) {
    f <- files[[sp]]
    if (!file.exists(f)) return(list(species = sp, gbif_ids = character(0), total = 0L))
    df <- tryCatch(
      readr::read_csv(
        f, show_col_types = FALSE,
        col_types = readr::cols(.default = readr::col_skip(), gbifID = readr::col_character())
      ),
      error = function(e) NULL
    )
    if (is.null(df) || !nrow(df)) return(list(species = sp, gbif_ids = character(0), total = 0L))
    ids <- unique(df$gbifID[!is.na(df$gbifID) & nzchar(df$gbifID)])
    list(species = sp, gbif_ids = ids, total = length(ids))
  })
  tibble::tibble(
    species = vapply(idx, `[[`, character(1), "species"),
    total   = vapply(idx, `[[`, integer(1), "total"),
    gbif_ids = I(lapply(idx, `[[`, "gbif_ids"))
  )
})

# ===================
# ======== UI =======
# ===================

# Normalized sets for JS rendering (dropdown icons)
BTOL_SET <- unique(normalize_species(BRASSITOL_SPECIES))
WOODY_SET <- unique(normalize_species(DERIVED_WOODY_SPECIES))

base_ui <- fluidPage(
  tags$head(
    tags$link(rel="icon", type="image/png", href="favicon.png?v=3"),
    tags$link(rel="icon", sizes="16x16", href="favicon-16.png?v=3"),
    tags$link(rel="icon", sizes="32x32", href="favicon-32.png?v=3"),
    tags$link(rel="apple-touch-icon", sizes="180x180", href="apple-touch-icon.png?v=3"),
    # Pass sets + icon URLs to JavaScript
    tags$script(HTML(sprintf("
      window.BTOL_SET = %s;
      window.WOODY_SET = %s;
      window.PHYLO_ICON = %s;
      window.WOODY_ICON = %s;
    ",
                             toJSON(unname(BTOL_SET), auto_unbox = TRUE),
                             toJSON(unname(WOODY_SET), auto_unbox = TRUE),
                             toJSON(PHYLO_ICON_URL %||% "", auto_unbox = TRUE),
                             toJSON(WOODY_ICON_URL %||% "", auto_unbox = TRUE)
    ))),
    tags$style(HTML("
      html, body { height: 100%; margin: 0; padding: 0; }
      .container-fluid { height: 100%; }
      .leaflet-container { height: calc(100vh - 80px) !important; }
      .form-group { max-width: 200px; }
      .rate-wrap { margin: 6px 0 10px; }
      .rate-btn { display:inline-block; padding:4px 8px; margin:0 6px 6px 0; border-radius:4px; border:1px solid rgba(0,0,0,0.2); text-decoration:none; cursor:pointer; font-weight:600; }
      .rate-good { background:#4CAF50; color:#fff; }
      .rate-neutral { background:#ADD8E6; color:#000; }
      .rate-bad { background:#e53935; color:#fff; }
      .rate-active { box-shadow:0 0 0 2px rgba(0,0,0,0.25) inset; }
      .bulk-rate-table { width:100%; border-collapse:collapse; margin-top:8px; font-size:90%; }
      .bulk-rate-table th, .bulk-rate-table td { padding:4px 4px; border-bottom:1px solid #eee; vertical-align:middle; }
      .bulk-rate-btns .rate-btn { padding:2px 6px; font-size:12px; margin-right:4px; }
      .bulk-rate-row-label { display:flex; align-items:center; gap:6px; }
      .cluster-chip { width:12px; height:12px; border-radius:2px; border:1px solid rgba(0,0,0,0.2); }
      .instr-box { background:#f7f9fc; border:1px solid #e2e6ef; border-radius:6px; padding:8px 10px; margin-bottom:10px; font-size:90%; line-height:1.3; }
      .instr-box b { display:block; margin-bottom:4px; }
      .selectize-dropdown .optgroup-header { font-weight:700; color:#2b4b6f; background:#f1f6fb; border-bottom:1px solid #e4edf6; }
      .selectize-dropdown .no-data, .selectize-input .no-data { color:#666; font-style:italic; }
      .notice-box { background: rgba(255,255,255,0.95); padding: 10px 12px; border-radius: 6px; border: 1px solid #d9e3f0; line-height: 1.25; box-shadow: 0 2px 6px rgba(0,0,0,0.06); max-width: 320px; margin-bottom: 6px; }
      /* Dropdown icon styling */
      .sel-ico { height: 16px; width: 16px; margin-right: 6px; vertical-align: text-bottom; }
      .sel-ico-fallback { margin-right: 6px; }
      .opt-row { display: flex; align-items: center; gap: 2px; }
      /* KPI cards */
      .kpi { background:#ffffff; border:1px solid #e6ecf5; border-radius:8px; padding:8px 10px; margin-bottom:8px; box-shadow:0 1px 2px rgba(0,0,0,0.04); }
      .kpi b { font-size: 95%; }
      .bar { height: 7px; background:#f0f3f7; border-radius: 4px; overflow: hidden; margin-top: 4px; }
      .bar > div { height: 100%; }
    "))
  ),
  titlePanel("Brassicaceae Occurrence Records Expert Curation Tool"),
  sidebarLayout(
    sidebarPanel(
      actionButton("help_btn", "How this tool works", class = "btn-info"),
      br(),
      # ===== Project dashboard (global overview) =====
      uiOutput("project_dashboard"),
      br(),
      tags$div(id = "tribe_wrap",
               selectInput("tribe", "Tribe", choices = c("All"))
      ),
      tags$div(id = "species_wrap",
               selectizeInput("species", "Species", choices = NULL, selected = NULL,
                              options = list(
                                placeholder = "Type or select a species",
                                render = I("{
                                  option: function(item, escape) {
                                    var raw = item.label || '';
                                    var text = escape(raw);
                                    var val = item.value || '';
                                    var norm = val.replace(/_/g,' ').trim().toLowerCase();
                                    var icons = '';
                                    function iconHTML(url, emoji, alt){
                                      if (url && url.length) {
                                        return '<span class=ico-wrap>' +
                                               '<img src=\"' + url + '\" class=\"sel-ico\" alt=\"' + alt + '\" ' +
                                               'onerror=\"this.style.display=\\'none\\'; this.nextSibling.style.display=\\'inline-block\\';\">' +
                                               '<span class=sel-ico-fallback style=\"display:none\">' + emoji + '</span>' +
                                               '</span>';
                                      } else {
                                        return '<span class=sel-ico-fallback>' + emoji + '</span>';
                                      }
                                    }
                                    if (Array.isArray(window.BTOL_SET) && window.BTOL_SET.indexOf(norm) >= 0)
                                      icons += iconHTML(window.PHYLO_ICON, '🌳', 'phylo');
                                    if (Array.isArray(window.WOODY_SET) && window.WOODY_SET.indexOf(norm) >= 0)
                                      icons += iconHTML(window.WOODY_ICON, '🌲', 'woody');
                                    return '<div class=opt-row>' + icons + '<span>' + text + '</span></div>';
                                  },
                                  item: function(item, escape) {
                                    var raw = item.label || '';
                                    var text = escape(raw);
                                    var val = item.value || '';
                                    var norm = val.replace(/_/g,' ').trim().toLowerCase();
                                    var icons = '';
                                    function iconHTML(url, emoji, alt){
                                      if (url && url.length) {
                                        return '<span class=ico-wrap>' +
                                               '<img src=\"' + url + '\" class=\"sel-ico\" alt=\"' + alt + '\" ' +
                                               'onerror=\"this.style.display=\\'none\\'; this.nextSibling.style.display=\\'inline-block\\';\">' +
                                               '<span class=sel-ico-fallback style=\"display:none\">' + emoji + '</span>' +
                                               '</span>';
                                      } else {
                                        return '<span class=sel-ico-fallback>' + emoji + '</span>';
                                      }
                                    }
                                    if (Array.isArray(window.BTOL_SET) && window.BTOL_SET.indexOf(norm) >= 0)
                                      icons += iconHTML(window.PHYLO_ICON, '🌳', 'phylo');
                                    if (Array.isArray(window.WOODY_SET) && window.WOODY_SET.indexOf(norm) >= 0)
                                      icons += iconHTML(window.WOODY_ICON, '🌲', 'woody');
                                    return '<div class=opt-row>' + icons + '<span>' + text + '</span></div>';
                                  }
                                }")
                              ))
      ),
      tags$div(id = "show_map_wrap",
               actionButton("show_map", "SHOW MAP")
      ),
      br(),
      uiOutput("cluster_controls"),
      # ===== Species summary (now includes rating counts) =====
      htmlOutput("species_summary"),
      br(),
      uiOutput("whoami"),
      textInput("display_name",
                "Name to record with ratings (if other than signin name)",
                value = "", placeholder = "e.g., K. Hendriks"),
      helpText("Default is your login name if captured; you can override it here."),
      br(),
      uiOutput("bulk_rate_panel")
    ),
    mainPanel(
      leafletOutput("map", height = "100%")
    )
  )
)

ui <- secure_app(base_ui)

# ==================================================
# ===== Helpers that must exist at top-level =======
# ==================================================

# species choices grouped
build_species_choices <- function(taxonomy_df, tribe_sel, ratings_df, species_index, tol_set_norm, woody_set_norm) {
  sp_vec <- if (is.null(tribe_sel) || identical(tribe_sel, "All")) {
    sort(unique(taxonomy_df$species))
  } else {
    taxonomy_df |>
      filter(tribe == tribe_sel) |>
      arrange(species) |>
      pull(species) |>
      unique()
  }
  
  rated_ids <- if (!is.null(ratings_df) && nrow(ratings_df) > 0) unique(as.character(ratings_df$gbifID)) else character(0)
  
  stats <- species_index |>
    filter(species %in% sp_vec) |>
    mutate(
      rated = vapply(gbif_ids, function(v) sum(!is.na(v) & nzchar(v) & v %in% rated_ids), integer(1)),
      unrated = pmax(total - rated, 0L),
      disp = pretty_species(species),
      tol = normalize_species(disp) %in% tol_set_norm,
      woody = normalize_species(disp) %in% woody_set_norm
    ) |>
    arrange(species) |>
    mutate(
      label = dplyr::case_when(
        total == 0L ~ sprintf("%s (no records)", disp),
        unrated > 0 ~ sprintf("%s (%d total; %d unrated)", disp, total, unrated),
        TRUE ~ sprintf("%s (%d total; 0 unrated)", disp, total)
      ),
      value = species
    )
  
  needs <- stats |>
    filter(total > 0L, unrated > 0L)
  full <- stats |>
    filter(total > 0L, unrated == 0L)
  nodata <- stats |>
    filter(total == 0L)
  
  out <- list()
  if (nrow(needs)) out[["Needs attention"]] <- setNames(as.list(needs$value), needs$label)
  if (nrow(full))  out[["Fully rated"]]    <- setNames(as.list(full$value),  full$label)
  if (nrow(nodata)) out[["No records"]]    <- setNames(as.list(nodata$value), nodata$label)
  out
}

# tribe choices with counts
build_tribe_choices <- function(taxonomy_df, ratings_df, species_index) {
  rated_ids <- if (!is.null(ratings_df) && nrow(ratings_df) > 0) unique(as.character(ratings_df$gbifID)) else character(0)
  
  per_species <- species_index |>
    mutate(
      rated = vapply(gbif_ids, function(v) sum(!is.na(v) & nzchar(v) & v %in% rated_ids), integer(1)),
      unrated = pmax(total - rated, 0L),
      needs = total > 0L & unrated > 0L
    ) |>
    left_join(taxonomy_df |> select(species, tribe), by = "species")
  
  tribe_sum <- per_species |>
    group_by(tribe) |>
    summarise(
      n_species = n_distinct(species),
      n_need = sum(needs, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(tribe)
  
  tribe_labels <- setNames(
    tribe_sum$tribe,
    sprintf("%s (%d spp; %d need review)", tribe_sum$tribe, tribe_sum$n_species, tribe_sum$n_need)
  )
  
  all_n_species <- n_distinct(taxonomy_df$species)
  all_n_need <- sum(per_species$needs, na.rm = TRUE)
  all_label <- setNames("All", sprintf("All (%d spp; %d need review)", all_n_species, all_n_need))
  
  c(all_label, tribe_labels)
}

# ======================
# ======= SERVER =======
# ======================

server <- function(input, output, session) {
  
  # ---- Auth ----
  my_check <- function(user, password) {
    ok <- nzchar(user) && identical(password, SHARED_PW)
    list(result = ok, user = user, admin = FALSE, expire = NA, roles = NA, user_info = list(user = user))
  }
  res_auth <- secure_server(check_credentials = my_check)
  current_user <- reactive(res_auth$user)
  
  observe({
    u <- current_user()
    if (!is.null(u) && nzchar(u) && !nzchar(isolate(input$display_name))) {
      updateTextInput(session, "display_name", value = u)
    }
  })
  
  output$whoami <- renderUI({
    nm <- if (nzchar(input$display_name)) input$display_name else current_user()
    HTML(sprintf("<b>Signed in as:</b> %s", if (!is.null(nm) && nzchar(nm)) nm else "(unknown)"))
  })
  
  # ---- Persist map view ----
  map_view <- reactiveValues(lat = NULL, lng = NULL, zoom = NULL)
  observeEvent(input$map_center, {
    if (!is.null(input$map_center)) {
      map_view$lat <- input$map_center$lat; map_view$lng <- input$map_center$lng
    }}, ignoreInit = TRUE)
  observeEvent(input$map_zoom, { if (!is.null(input$map_zoom)) map_view$zoom <- input$map_zoom }, ignoreInit = TRUE)
  last_drawn_species <- reactiveVal(NULL)
  
  # ---- Ratings cache + initial load ----
ratings_df <- reactiveVal(tibble::tibble(
  gbifID = character(),
  user = character(),
  rating = character(),
  timestamp_iso = character()
))

observe({
  if (file.exists(RATINGS_FILE)) {
    df <- readr::read_csv(RATINGS_FILE, show_col_types = FALSE) |>
      dplyr::transmute(
        gbifID = as.character(gbifID),
        user = as.character(user),
        rating = tolower(trimws(as.character(rating))),
        timestamp_iso = as.character(timestamp_iso)
      )
    ratings_df(df)
  }
})
  ratings_for_dropdown <- reactive(ratings_df()) %>% debounce(600)
  
  # ========= PROJECT DASHBOARD STATS =========
  project_stats <- reactive({
    idx <- species_gbif_index
    total_records <- sum(idx$total, na.rm = TRUE)
    total_species <- length(all_species)
    species_with_records <- sum(idx$total > 0, na.rm = TRUE)
    
    r <- ratings_df()
    latest <- NULL
    if (!is.null(r) && nrow(r)) {
      latest <- r |>
        mutate(ts = ymd_hms(timestamp_iso, quiet = TRUE)) |>
        arrange(gbifID, desc(ts)) |>
        group_by(gbifID) |>
        summarise(latest = dplyr::first(rating), .groups = "drop")
    } else {
      latest <- tibble::tibble(gbifID = character(0), latest = character(0))
    }
    
    n_good   <- sum(latest$latest == "good", na.rm = TRUE)
    n_neutral<- sum(latest$latest == "neutral", na.rm = TRUE)
    n_bad    <- sum(latest$latest == "bad", na.rm = TRUE)
    n_rated  <- n_good + n_neutral + n_bad
    n_unrated<- max(total_records - n_rated, 0L)
    
    # species fully rated
    rated_id_set <- unique(latest$gbifID)
    per_species <- idx |>
      mutate(
        rated = vapply(gbif_ids, function(v) sum(!is.na(v) & nzchar(v) & v %in% rated_id_set), integer(1)),
        unrated = pmax(total - rated, 0L)
      )
    species_fully_rated <- sum(per_species$total > 0 & per_species$unrated == 0, na.rm = TRUE)
    species_needing     <- sum(per_species$total > 0 & per_species$unrated > 0, na.rm = TRUE)
    
    list(
      total_records = total_records,
      total_species = total_species,
      species_with_records = species_with_records,
      n_good = n_good, n_neutral = n_neutral, n_bad = n_bad, n_unrated = n_unrated,
      n_rated = n_rated,
      species_fully_rated = species_fully_rated,
      species_needing = species_needing
    )
  })
  
  pct <- function(a,b) ifelse(b > 0, round(100 * a / b, 1), 0)
  
  output$project_dashboard <- renderUI({
    ps <- project_stats()
    tot <- ps$total_records
    rated <- ps$n_rated
    unrated <- ps$n_unrated
    sp_tot <- ps$total_species
    sp_with <- ps$species_with_records
    sp_full <- ps$species_fully_rated
    sp_need <- ps$species_needing
    
    # bar widths
    w_good <- pct(ps$n_good, tot)
    w_neu  <- pct(ps$n_neutral, tot)
    w_bad  <- pct(ps$n_bad, tot)
    
    tagList(
      tags$div(class = "kpi",
               HTML(sprintf("<b>Project overview</b><br>
                 <span>Total records:</span> <b>%s</b> &nbsp;·&nbsp; <span>Rated:</span> <b>%s</b> (%s%%) &nbsp;·&nbsp; <span>Unrated:</span> <b>%s</b> (%s%%)",
                            format(tot, big.mark=","), format(rated, big.mark=","), pct(rated, tot),
                            format(unrated, big.mark=","), pct(unrated, tot))),
               tags$div(class="bar",
                        tags$div(style = sprintf("width:%s%%; background:#4CAF50; float:left;", w_good)),
                        tags$div(style = sprintf("width:%s%%; background:#ADD8E6; float:left;", w_neu)),
                        tags$div(style = sprintf("width:%s%%; background:#e53935; float:left;", w_bad))
               ),
               tags$small(HTML("&nbsp;Green = good, Blue = neutral, Red = bad"))
      ),
      tags$div(class = "kpi",
               HTML(sprintf("<span>Total species:</span> <b>%s</b> &nbsp;·&nbsp; <span>With records:</span> <b>%s</b><br>
                             <span>Fully rated species:</span> <b>%s</b> (%s%% of with-records) &nbsp;·&nbsp;
                             <span>Still to rate:</span> <b>%s</b>",
                            format(sp_tot, big.mark=","), format(sp_with, big.mark=","),
                            format(sp_full, big.mark=","), pct(sp_full, sp_with),
                            format(sp_need, big.mark=",")))
      )
    )
  })
  
  # ---- Tribe dropdown (preserve selection) ----
  observeEvent(ratings_for_dropdown(), {
    choices <- build_tribe_choices(taxonomy, ratings_df(), species_gbif_index)
    current <- isolate(input$tribe)
    valid_values <- unname(choices)
    if (is.null(current) || !nzchar(current) || !(current %in% valid_values)) current <- "All"
    updateSelectInput(session, "tribe", choices = choices, selected = current)
  }, ignoreInit = FALSE)
  
  # ---- Species dropdown (JS-rendered icons) ----
  tol_set_norm <- unique(normalize_species(BRASSITOL_SPECIES))
  woody_set_norm <- unique(normalize_species(DERIVED_WOODY_SPECIES))
  observeEvent(list(input$tribe, ratings_for_dropdown()), {
    choices_named <- build_species_choices(
      taxonomy, input$tribe %||% "All", ratings_df(), species_gbif_index, tol_set_norm, woody_set_norm
    )
    valid_values <- unname(unlist(choices_named, use.names = FALSE))
    current <- isolate(input$species)
    if (is.null(current) || !nzchar(current) || !(current %in% valid_values)) current <- NULL
    
    updateSelectizeInput(
      session, "species",
      choices = choices_named, selected = current, server = TRUE,
      options = list(
        placeholder = if (length(choices_named)) "Type or select a species" else "No species for this tribe",
        render = I("{
          option: function(item, escape) {
            var raw = item.label || '';
            var text = escape(raw);
            var val = item.value || '';
            var norm = val.replace(/_/g,' ').trim().toLowerCase();
            var icons = '';
            function iconHTML(url, emoji, alt){
              if (url && url.length) {
                return '<span class=ico-wrap>' +
                       '<img src=\"' + url + '\" class=\"sel-ico\" alt=\"' + alt + '\" ' +
                       'onerror=\"this.style.display=\\'none\\'; this.nextSibling.style.display=\\'inline-block\\';\">' +
                       '<span class=sel-ico-fallback style=\"display:none\">' + emoji + '</span>' +
                       '</span>';
              } else {
                return '<span class=sel-ico-fallback>' + emoji + '</span>';
              }
            }
            if (Array.isArray(window.BTOL_SET) && window.BTOL_SET.indexOf(norm) >= 0)
              icons += iconHTML(window.PHYLO_ICON, '🌳', 'phylo');
            if (Array.isArray(window.WOODY_SET) && window.WOODY_SET.indexOf(norm) >= 0)
              icons += iconHTML(window.WOODY_ICON, '🌲', 'woody');
            return '<div class=opt-row>' + icons + '<span>' + text + '</span></div>';
          },
          item: function(item, escape) {
            var raw = item.label || '';
            var text = escape(raw);
            var val = item.value || '';
            var norm = val.replace(/_/g,' ').trim().toLowerCase();
            var icons = '';
            function iconHTML(url, emoji, alt){
              if (url && url.length) {
                return '<span class=ico-wrap>' +
                       '<img src=\"' + url + '\" class=\"sel-ico\" alt=\"' + alt + '\" ' +
                       'onerror=\"this.style.display=\\'none\\'; this.nextSibling.style.display=\\'inline-block\\';\">' +
                       '<span class=sel-ico-fallback style=\"display:none\">' + emoji + '</span>' +
                       '</span>';
              } else {
                return '<span class=sel-ico-fallback>' + emoji + '</span>';
              }
            }
            if (Array.isArray(window.BTOL_SET) && window.BTOL_SET.indexOf(norm) >= 0)
              icons += iconHTML(window.PHYLO_ICON, '🌳', 'phylo');
            if (Array.isArray(window.WOODY_SET) && window.WOODY_SET.indexOf(norm) >= 0)
              icons += iconHTML(window.WOODY_ICON, '🌲', 'woody');
            return '<div class=opt-row>' + icons + '<span>' + text + '</span></div>';
          }
        }")
      )
    )
  }, ignoreInit = FALSE)
  
  # ---- Selected species (latched on SHOW MAP) ----
  selected_species <- eventReactive(input$show_map, { req(input$species); input$species })
  
  # ---- Load selected species data ----
  species_data <- eventReactive(input$show_map, {
    req(input$species)
    f <- file.path("occurrence_data", paste0(input$species, ".csv"))
    if (!file.exists(f)) {
      df0 <- tibble::tibble()[0, ]
      return(ensure_required_cols(df0))
    }
    df <- read_csv(f, show_col_types = FALSE,
                   col_types = cols(.default = col_guess(), gbifID = col_character()))
    df <- ensure_required_cols(df)
    df |>
      mutate(
        decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
        decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude)),
        coordinate_source = tolower(trimws(as.character(coordinate_source))),
        coordinateUncertaintyInMeters = suppressWarnings(as.numeric(coordinateUncertaintyInMeters)),
        elevation = suppressWarnings(as.numeric(elevation)),
        elevationAccuracy= suppressWarnings(as.numeric(elevationAccuracy))
      )
  })
  
  # ---- Cluster UI ----
  clusters_available <- reactive({
    d <- species_data()
    !is.null(d) && "cluster" %in% names(d) && any(!is.na(d$cluster))
  })
  output$cluster_controls <- renderUI({
    req(input$show_map)
    if (!clusters_available()) return(NULL)
    radioButtons(
      "color_by", "Colour points by",
      choices = c("Rating" = "rating", "Cluster" = "cluster"),
      selected = "rating"
    )
  })
  
  # ---- Bulk-rate panel ----
  output$bulk_rate_panel <- renderUI({
    req(input$show_map)
    d <- species_data(); req(d)
    
    have_clusters <- clusters_available()
    cluster_rows <- NULL
    
    # Counts for geocoded / gbif
    n_all <- nrow(d)
    n_geocoded <- sum(tolower(d$coordinate_source) != "gbif", na.rm = TRUE)
    n_gbif <- sum(tolower(d$coordinate_source) == "gbif", na.rm = TRUE)
    
    if (have_clusters) {
      d$cluster <- as.factor(as.character(d$cluster))
      lvls <- natural_cluster_levels(d$cluster)
      d$cluster <- factor(d$cluster, levels = lvls)
      lvl_to_am <- setNames(AM_COLORS[(seq_along(lvls)-1) %% length(AM_COLORS) + 1], lvls)
      
      cs <- d |>
        filter(!is.na(cluster)) |>
        count(cluster, name = "n") |>
        mutate(col = unname(lvl_to_am[as.character(cluster)]))
      
      cluster_rows <- lapply(seq_len(nrow(cs)), function(i) {
        cl <- as.character(cs$cluster[i]); col <- cs$col[i]; n <- cs$n[i]
        htmltools::tags$tr(
          htmltools::tags$td(
            class = "bulk-rate-row-label",
            htmltools::tags$span(class = "cluster-chip", style = paste0("background:", col, ";")),
            htmltools::tags$span(sprintf("Cluster %s", cl)),
            htmltools::tags$span(sprintf("(n=%d)", n), style = "color:#666; font-size:90%;")
          ),
          htmltools::tags$td(
            class = "bulk-rate-btns",
            htmltools::tags$a(class = "rate-btn rate-good",
                              "Good",
                              onclick = sprintf("Shiny.setInputValue('bulk_rate_click', {scope:'cluster', cluster:'%s', rating:'good', nonce:Math.random()}, {priority:'event'})", htmltools::htmlEscape(cl))),
            htmltools::tags$a(class = "rate-btn rate-neutral",
                              "Neutral",
                              onclick = sprintf("Shiny.setInputValue('bulk_rate_click', {scope:'cluster', cluster:'%s', rating:'neutral', nonce:Math.random()}, {priority:'event'})", htmltools::htmlEscape(cl))),
            htmltools::tags$a(class = "rate-btn rate-bad",
                              "Bad",
                              onclick = sprintf("Shiny.setInputValue('bulk_rate_click', {scope:'cluster', cluster:'%s', rating:'bad', nonce:Math.random()}, {priority:'event'})", htmltools::htmlEscape(cl)))
          )
        )
      })
    }
    
    # All points row
    all_row <- htmltools::tags$tr(
      htmltools::tags$td(
        class = "bulk-rate-row-label",
        htmltools::tags$span("All points"),
        htmltools::tags$span(sprintf("(n=%d)", n_all), style = "color:#666; font-size:90%;")
      ),
      htmltools::tags$td(
        class = "bulk-rate-btns",
        htmltools::tags$a(class = "rate-btn rate-good",
                          "Good",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'all', rating:'good', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-neutral",
                          "Neutral",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'all', rating:'neutral', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-bad",
                          "Bad",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'all', rating:'bad', nonce:Math.random()}, {priority:'event'})")
      )
    )
    
    # GBIF-only row
    gbif_row <- htmltools::tags$tr(
      htmltools::tags$td(
        class = "bulk-rate-row-label",
        htmltools::tags$span("GBIF coordinate points"),
        htmltools::tags$span(sprintf("(n=%d)", n_gbif), style = "color:#666; font-size:90%;")
      ),
      htmltools::tags$td(
        class = "bulk-rate-btns",
        htmltools::tags$a(class = "rate-btn rate-good",
                          "Good",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'gbif', rating:'good', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-neutral",
                          "Neutral",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'gbif', rating:'neutral', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-bad",
                          "Bad",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'gbif', rating:'bad', nonce:Math.random()}, {priority:'event'})")
      )
    )
    
    # Geocoded-only row
    geocoded_row <- htmltools::tags$tr(
      htmltools::tags$td(
        class = "bulk-rate-row-label",
        htmltools::tags$span("Geocoded points"),
        htmltools::tags$span(sprintf("(n=%d)", n_geocoded), style = "color:#666; font-size:90%;")
      ),
      htmltools::tags$td(
        class = "bulk-rate-btns",
        htmltools::tags$a(class = "rate-btn rate-good",
                          "Good",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'geocoded', rating:'good', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-neutral",
                          "Neutral",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'geocoded', rating:'neutral', nonce:Math.random()}, {priority:'event'})"),
        htmltools::tags$a(class = "rate-btn rate-bad",
                          "Bad",
                          onclick = "Shiny.setInputValue('bulk_rate_click', {scope:'geocoded', rating:'bad', nonce:Math.random()}, {priority:'event'})")
      )
    )
    
    tbody_children <- htmltools::tagList(all_row, gbif_row, geocoded_row, cluster_rows %||% list())
    htmltools::tags$div(
      htmltools::tags$h5(htmltools::tags$strong("Bulk rate option:")),
      htmltools::tags$p("Apply a rating to all points, only GBIF points, only geocoded points, or to an entire cluster:"),
      htmltools::tags$table(
        class = "bulk-rate-table",
        htmltools::tags$thead(
          htmltools::tags$tr(
            htmltools::tags$th("Scope"),
            htmltools::tags$th("Set rating")
          )
        ),
        htmltools::tags$tbody(tbody_children)
      )
    )
  })
  
  # ---- Map base ----
  output$map <- renderLeaflet({
    req(input$show_map)
    sp <- selected_species()
    
    native_range_string <- taxonomy$native_range[taxonomy$species == sp]
    native_polygons <- NULL; native_polygons_buffered <- NULL
    if (length(native_range_string) == 1 && !is.na(native_range_string)) {
      native_range_codes <- unlist(strsplit(native_range_string, ",\\s*"))
      native_polygons <- tryCatch(
        wgsrpd_L3[wgsrpd_L3$LEVEL3_COD %in% native_range_codes, ],
        error = function(e) NULL
      )
      if (!is.null(native_polygons) && nrow(native_polygons) > 0) {
        native_polygons_buffered <- st_transform(native_polygons, 3857) |>
          st_buffer(dist = 100000) |>
          st_transform(4326)
      }
    }
    
    m <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
      addProviderTiles("OpenStreetMap", group = "Street Map") %>%
      addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
      addProviderTiles("Esri.WorldTopoMap", group = "Topographic") %>%
      addMapPane("native", zIndex = 300) %>%
      addMapPane("native_buffer", zIndex = 310) %>%
      addMapPane("hulls", zIndex = 420) %>%
      addMapPane("points", zIndex = 430)
    
    if (!is.null(native_polygons_buffered) && nrow(native_polygons_buffered) > 0) {
      m <- m %>% addPolygons(
        data = native_polygons_buffered, options = pathOptions(pane = "native_buffer"),
        fillColor = "lightgreen", fillOpacity = 0.2, color = NA, group = "Buffered Native Range"
      )
    }
    if (!is.null(native_polygons) && nrow(native_polygons) > 0) {
      m <- m %>% addPolygons(
        data = native_polygons, options = pathOptions(pane = "native"),
        fillColor = "lightblue", fillOpacity = 0.3, color = "blue", weight = 1, group = "Native Range"
      )
    }
    
    m %>% addLayersControl(
      baseGroups = c("Street Map", "Satellite", "Topographic"),
      overlayGroups = c("Native Range", "Buffered Native Range", "Cluster hulls", "Cluster labels", "Points"),
      options = layersControlOptions(collapsed = FALSE)
    )
  })
  
  # ---- Redraw points / hulls / labels + notices ----
  observeEvent(list(input$show_map, input$color_by, ratings_df()), {
    d <- species_data(); req(d)
    sp <- selected_species()
    color_by <- input$color_by %||% "rating"
    lp <- leafletProxy("map")
    
    cur_lat <- isolate(map_view$lat); cur_lng <- isolate(map_view$lng); cur_zoom <- isolate(map_view$zoom)
    active_overlays <- isolate(input$map_groups %||% character(0))
    
    lp <- lp %>%
      clearMarkerClusters() %>%
      clearGroup("Points") %>%
      clearGroup("Cluster hulls") %>%
      clearGroup("Cluster labels") %>%
      clearControls()
    
    fit_to <- NULL
    native_range_string <- taxonomy$native_range[taxonomy$species == sp]
    if (length(native_range_string) == 1 && !is.na(native_range_string) && !is.null(wgsrpd_L3)) {
      codes <- unlist(strsplit(native_range_string, ",\\s*"))
      polys <- tryCatch(wgsrpd_L3[wgsrpd_L3$LEVEL3_COD %in% codes, ], error = function(e) NULL)
      if (!is.null(polys) && nrow(polys) > 0) {
        bb <- sf::st_bbox(sf::st_transform(polys, 4326))
        fit_to <- list(west = bb["xmin"], south = bb["ymin"], east = bb["xmax"], north = bb["ymax"])
      }
    }
    
    d_clean0 <- d %>% filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
    if (nrow(d_clean0) == 0) {
      lp <- lp %>% addControl(
        html = htmltools::HTML(sprintf(
          "<div class='notice-box'>
             <b>No occurrence records found for this species.</b><br/>
             If you have coordinates/records to share, please email:
             <a href='mailto:%1$s'>%1$s</a>.
           </div>", CONTACT_EM)),
        position = "bottomleft", layerId = "empty_msg"
      )
      lp <- lp %>% addLayersControl(
        baseGroups = c("Street Map", "Satellite", "Topographic"),
        overlayGroups = c("Native Range", "Buffered Native Range", "Cluster hulls", "Cluster labels", "Points"),
        options = layersControlOptions(collapsed = FALSE)
      )
      if (!is.null(fit_to)) {
        lp <- lp %>% fitBounds(fit_to$west, fit_to$south, fit_to$east, fit_to$north)
        map_view$lat <- map_view$lng <- map_view$zoom <- NULL
      } else if (!is.null(cur_lat) && !is.null(cur_lng) && !is.null(cur_zoom)) {
        lp <- lp %>% setView(lng = cur_lng, lat = cur_lat, zoom = cur_zoom)
      }
      last_drawn_species(sp)
      return(invisible())
    }
    
    # Few-records notice (<20)
    n_total <- nrow(d)
    if (n_total < 20) {
      lp <- lp %>% addControl(
        html = htmltools::HTML(sprintf(
          "<div class='notice-box'>
             <b>We have few records for this species.</b><br/>
             We aim for at least <b>10</b> records rated \"good\" per species.<br/>
             If you have additional records, please email:
             <a href='mailto:%1$s'>%1$s</a>.
           </div>", CONTACT_EM)),
        position = "bottomleft", layerId = "few_msg"
      )
    }
    
    RATING2AMCOL <- c(good = "green", neutral = "lightblue", bad = "red", unrated = "lightgray")
    
    have_clusters <- clusters_available()
    if (have_clusters) {
      d$cluster <- as.factor(as.character(d$cluster))
      lvls <- natural_cluster_levels(d$cluster)
      d$cluster <- factor(d$cluster, levels = lvls)
      lvl_to_am <- setNames(AM_COLORS[(seq_along(lvls)-1) %% length(AM_COLORS) + 1], lvls)
    }
    
    d <- d |>
      mutate(across(where(is.character), linkify)) |>
      mutate(src_key = ifelse(coordinate_source == "gbif", "gbif", "other"),
             src_lbl = ifelse(src_key == "gbif", "GBIF", "Geocoded"))
    
    lr <- latest_ratings_for(unique(as.character(d$gbifID)), ratings_df())
    
    d <- d |>
      mutate(gbifID = as.character(gbifID)) |>
      left_join(lr, by = "gbifID") |>
      mutate(
        rating_display = if_else(!is.na(latest) & latest %in% c("good","neutral","bad"), latest, "unrated"),
        recent_html = vapply(
          last3,
          function(x) if (length(x) == 0 || all(is.na(x))) "(no recent ratings)"
          else paste0(paste0("&bull; ", x), collapse = "<br>"),
          character(1)
        ),
        btn_good_class   = ifelse(rating_display == "good", " rate-active", ""),
        btn_neutral_class= ifelse(rating_display == "neutral", " rate-active", ""),
        btn_bad_class    = ifelse(rating_display == "bad", " rate-active", "")
      ) |>
      mutate(
        popup_html = glue::glue("
          <b>RECORD DETAILS:</b><br>
          <b>GBIF ID:</b> <a href='https://www.gbif.org/occurrence/{gbifID}' target='_blank'>{gbifID}</a><br>
          <b>Species:</b> {SPECIES_NAME_PRINT}<br>
          <b>Tribe:</b> {TRIBE_FULL}<br><br>
          <b>LOCATION DETAILS:</b><br>
          <b>Higher Geography:</b> {higherGeography}<br>
          <b>Country Code:</b> {countryCode}<br>
          <b>Country (Full):</b> {country_full_name}<br>
          <b>County:</b> {county}<br>
          <b>Locality:</b> {locality}<br>
          <b>Remarks:</b> {locationRemarks}<br>
          <b>State/Province:</b> {stateProvince}<br>
          <b>Municipality:</b> {municipality}<br>
          <b>Latitude:</b> {decimalLatitude}<br>
          <b>Longitude:</b> {decimalLongitude}<br>
          <b>Uncertainty (m):</b> {coordinateUncertaintyInMeters}<br>
          <b>Elevation (m):</b> {elevation}<br><br>
          <b>LOCATION CLEANING:</b><br>
          <b>Coordinate Source:</b> {src_lbl}<br>
          <b>Passed automated location curation tests:</b> {coordinate_cleaner_.summary}<br>
          <b>Geocoding Errors:</b> {geocoding_errors}<br><br>
          <b>DATA SOURCE:</b><br>
          <b>Basis of Record:</b> {basisOfRecord}<br>
          <b>Date:</b> {eventDate}<br>
          <b>Recorded By:</b> {recordedBy}<br>
          <b>Identified By:</b> {identifiedBy}<br>
          <b>Record Number:</b> {recordNumber}<br>
          <b>Dataset:</b> {datasetName}<br>
          <b>Occurrence ID:</b> {occurrenceID}<br><br>
          <b>RATINGS:</b><br>
          <b>Latest rating:</b> {ifelse(is.na(latest), '(none yet)', latest)}<br>
          <b>Latest three ratings:</b><br>{recent_html}<br>
          <b>Your rating please:</b><br>
          <div class='rate-wrap'>
            <a class='rate-btn rate-good{btn_good_class}' onclick=\"Shiny.setInputValue('rate_click', {{gbifID: '{gbifID}', rating: 'good', nonce: Math.random()}}, {{priority: 'event'}})\">Good</a>
            <a class='rate-btn rate-neutral{btn_neutral_class}' onclick=\"Shiny.setInputValue('rate_click', {{gbifID: '{gbifID}', rating: 'neutral', nonce: Math.random()}}, {{priority: 'event'}})\">Neutral</a>
            <a class='rate-btn rate-bad{btn_bad_class}' onclick=\"Shiny.setInputValue('rate_click', {{gbifID: '{gbifID}', rating: 'bad', nonce: Math.random()}}, {{priority: 'event'}})\">Bad</a>
          </div>
        ")
      )
    
    d <- d |>
      mutate(gbifID = as.character(gbifID)) |>
      filter(!is.na(gbifID) & nzchar(gbifID))
    
    d_clean <- d |> filter(is.finite(decimalLongitude), is.finite(decimalLatitude))
    if (!nrow(d_clean)) return(invisible())
    
    is_new_species <- !identical(last_drawn_species(), sp)
    if (is_new_species && is.null(fit_to) && nrow(d_clean) > 0) {
      lon <- range(d_clean$decimalLongitude, na.rm = TRUE)
      lat <- range(d_clean$decimalLatitude,  na.rm = TRUE)
      pad <- 0.05
      if (diff(lon) == 0) lon <- lon + c(-0.5, 0.5) else lon <- lon + c(-1, 1) * pad * diff(lon)
      if (diff(lat) == 0) lat <- lat + c(-0.5, 0.5) else lat <- lat + c(-1, 1) * pad * diff(lat)
      fit_to <- list(west = lon[1], south = lat[1], east = lon[2], north = lat[2])
    }
    
    hulls <- NULL; labels_sf <- NULL
    if (clusters_available() && nrow(d_clean) > 0) {
      pts <- st_as_sf(d_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
      ptsm <- st_transform(pts, 3857)
      hull_list <- lapply(split(ptsm, ptsm$cluster), function(g) {
        coords <- st_coordinates(g)
        if (nrow(unique(coords[, c("X","Y"), drop = FALSE])) < 3) return(NULL)
        ch <- tryCatch(st_convex_hull(st_union(g)), error = function(e) NULL)
        if (is.null(ch) || length(ch) == 0) return(NULL)
        st_geometry(ch)[[1]]
      })
      keep <- !vapply(hull_list, is.null, logical(1))
      if (any(keep)) {
        hulls <- st_sfc(hull_list[keep], crs = 3857) |>
          st_transform(4326) |>
          st_as_sf() |>
          mutate(cluster = names(hull_list)[keep], .before = 1)
      }
      
      labels_sf <- ptsm |>
        group_by(cluster) |>
        summarise(n = dplyr::n(), .groups = "drop") |>
        st_centroid() |>
        st_transform(4326) |>
        mutate(label_html = glue("Cluster {cluster}<br/>(n = {n})"))
    }
    
    if (!is.null(hulls) && nrow(hulls) > 0) {
      hulls$col <- unname(lvl_to_am[as.character(hulls$cluster)])
      lp <- lp %>% addPolygons(
        data = hulls, options = pathOptions(pane = "hulls"),
        fillColor = ~col, fillOpacity = 0.2, color = ~col, weight = 1,
        group = "Cluster hulls"
      )
    }
    
    d$src_key <- ifelse(d$coordinate_source == "gbif", "gbif", "other")
    d$icon_shape <- ifelse(d$src_key == "gbif", "circle", "play")
    make_icons <- function(marker_color_vec)
      awesomeIcons(icon = d$icon_shape, library = "fa", iconColor = "white", markerColor = marker_color_vec)
    
    if (clusters_available() && identical(color_by, "cluster")) {
      marker_cols <- unname(lvl_to_am[as.character(d$cluster)])
      icons <- make_icons(marker_cols)
      lp <- lp %>% addAwesomeMarkers(
        data = d,
        lng = ~decimalLongitude, lat = ~decimalLatitude,
        popup = ~lapply(popup_html, htmltools::HTML),
        icon = icons, options = markerOptions(pane = "points"),
        group = "Points", layerId = ~gbifID,
        clusterOptions = markerClusterOptions(
          showCoverageOnHover = FALSE, zoomToBoundsOnClick = TRUE,
          spiderfyOnMaxZoom = TRUE, disableClusteringAtZoom = 18, maxClusterRadius = 35,
          spiderLegPolylineOptions = list(weight = 1, opacity = 0.6), spiderfyDistanceMultiplier = 2
        )
      )
      lp <- lp %>% addLegend(
        position = "bottomright", colors = unname(lvl_to_am[lvls]), labels = lvls,
        title = "Cluster (fill)", opacity = 1
      )
    } else {
      d$rating_display <- ifelse(d$rating_display %in% names(RATING2AMCOL), d$rating_display, "unrated")
      marker_cols <- unname(RATING2AMCOL[d$rating_display])
      icons <- make_icons(marker_cols)
      lp <- lp %>% addAwesomeMarkers(
        data = d,
        lng = ~decimalLongitude, lat = ~decimalLatitude,
        popup = ~lapply(popup_html, htmltools::HTML),
        icon = icons, options = markerOptions(pane = "points"),
        group = "Points", layerId = ~gbifID,
        clusterOptions = markerClusterOptions(
          showCoverageOnHover = FALSE, zoomToBoundsOnClick = TRUE,
          spiderfyOnMaxZoom = TRUE, disableClusteringAtZoom = 18, maxClusterRadius = 35,
          spiderLegPolylineOptions = list(weight = 1, opacity = 0.6), spiderfyDistanceMultiplier = 2
        )
      )
      lp <- lp %>% addLegend(
        position = "bottomright",
        colors = c("green","lightblue","red","lightgray"),
        labels = c("good","neutral","bad","unrated"),
        title = "Rating (fill)", opacity = 1
      )
    }
    
    lp <- lp %>% addControl(
      html = HTML('
        <div class="notice-box" style="margin-bottom:6px;">
          <b>Source (shape)</b><br/>
          <i class="fa fa-circle" style="margin-right:6px;"></i> GBIF<br/>
          <i class="fa fa-play" style="margin-right:6px;"></i> Geocoded
        </div>
      '),
      position = "bottomleft", layerId = "shape-legend"
    )
    
    if (!is.null(labels_sf) && nrow(labels_sf) > 0) {
      lp <- lp %>% addLabelOnlyMarkers(
        data = labels_sf,
        lng = ~sf::st_coordinates(geometry)[,1],
        lat = ~sf::st_coordinates(geometry)[,2],
        label = ~lapply(label_html, htmltools::HTML),
        labelOptions = labelOptions(noHide = TRUE, textOnly = TRUE, direction = "center",
                                    textsize = "16px", style = list("font-weight" = "bold")),
        group = "Cluster labels"
      )
    }
    
    overlay_groups <- c("Native Range", "Buffered Native Range", "Points")
    if (clusters_available()) overlay_groups <- c("Native Range", "Buffered Native Range", "Cluster hulls", "Cluster labels", "Points")
    lp <- lp %>% addLayersControl(
      baseGroups = c("Street Map", "Satellite", "Topographic"),
      overlayGroups = overlay_groups,
      options = layersControlOptions(collapsed = FALSE)
    )
    
    if (length(active_overlays)) for (g in active_overlays) lp <- lp %>% showGroup(g)
    lp <- lp %>% showGroup("Points")
    
    if (!is.null(fit_to) && !identical(last_drawn_species(), sp)) {
      lp <- lp %>% fitBounds(fit_to$west, fit_to$south, fit_to$east, fit_to$north)
      map_view$lat <- map_view$lng <- map_view$zoom <- NULL
    } else if (!is.null(cur_lat) && !is.null(cur_lng) && !is.null(cur_zoom)) {
      lp <- lp %>% setView(lng = cur_lng, lat = cur_lat, zoom = cur_zoom)
    }
    last_drawn_species(sp)
  })
  
  # Map clicks helper
  selected_gbif <- reactiveVal(NULL)
  observeEvent(input$map_marker_click, {
    click <- input$map_marker_click
    if (!is.null(click$id) && nzchar(click$id)) selected_gbif(click$id)
    else showNotification("Click a single point (not a cluster). Zoom in if needed.", type = "warning")
  })
  
  # ---- Save ratings (bulk + single) ----
  observeEvent(input$bulk_rate_click, {

    showNotification(
  "This public/reviewer version is read-only; ratings are not saved.",
  type = "message"
    )
    return(invisible())
    
    info <- input$bulk_rate_click; req(info)
    rate <- tolower(info$rating %||% "")
    if (!nzchar(rate) || !rate %in% c("good","neutral","bad")) {
      showNotification("Invalid bulk rating.", type = "error"); return(invisible())
    }
    usr <- if (nzchar(input$display_name)) input$display_name else (current_user() %||% "")
    if (!nzchar(usr)) {
      showNotification("Please enter your name before saving.", type = "error"); return(invisible())
    }
    
    d <- species_data(); req(d)
    d <- d |>
      mutate(gbifID = as.character(gbifID)) |>
      filter(!is.na(gbifID) & nzchar(gbifID))
    
    scope <- info$scope %||% "all"
    gbif_ids <- if (identical(scope, "all")) {
      d$gbifID
    } else if (identical(scope, "cluster")) {
      cl <- info$cluster %||% NA
      if (is.na(cl)) { showNotification("No cluster specified.", type = "error"); return(invisible()) }
      d |>
        mutate(cluster = as.character(cluster)) |>
        filter(cluster == as.character(cl)) |>
        pull(gbifID)
    } else if (identical(scope, "geocoded")) {
      d |>
        filter(tolower(coordinate_source) != "gbif") |>
        pull(gbifID)
    } else if (identical(scope, "gbif")) {
      d |>
        filter(tolower(coordinate_source) == "gbif") |>
        pull(gbifID)
    } else {
      showNotification("Unknown bulk scope.", type = "error"); return(invisible())
    }
    
    gbif_ids <- unique(gbif_ids)
    if (!length(gbif_ids)) {
      showNotification("No records found for this selection.", type = "warning"); return(invisible())
    }
    
    spec_val <- input$species %||% ""
    tribe_val <- taxonomy$tribe[taxonomy$species == spec_val][1] %||% ""
    tribe_val <- ifelse(is.na(tribe_val), "", tribe_val)
    ts_iso <- sub('(\\+|\\-)(\\d{2})(\\d{2})$', '\\1\\2:\\3',
                  format(lubridate::with_tz(lubridate::now(), TIME_TZ), '%Y-%m-%dT%H:%M:%S%z'))
    entry <- tibble::tibble(
      gbifID = gbif_ids, user = usr, rating = rate,
      species = spec_val, tribe = tribe_val, timestamp_iso = ts_iso
    )
    
    if (is.null(ss)) ss <<- tryCatch(gs4_get(SHEET_ID), error = function(e) NULL)
    tr <- try({
      if (!is.null(ss)) sheet_append(ss, entry, sheet = SHEET_TAB)
      else sheet_append(SHEET_ID, entry, sheet = SHEET_TAB)
    }, silent = TRUE)
    
    if (inherits(tr, "try-error")) {
      msg <- conditionMessage(attr(tr, "condition"))
      message("[GS] bulk sheet_append failed: ", msg)
      showNotification(paste("Could not save bulk rating:", msg), type = "error", duration = 10)
      return(invisible())
    }
    
    ratings_df(bind_rows(ratings_df(), entry))
    
    scope_msg <- if (identical(scope, "all")) {
      sprintf("all %d points", length(gbif_ids))
    } else if (identical(scope, "geocoded")) {
      sprintf("all %d geocoded points", length(gbif_ids))
    } else if (identical(scope, "gbif")) {
      sprintf("all %d GBIF points", length(gbif_ids))
    } else {
      sprintf("cluster %s (%d points)", info$cluster, length(gbif_ids))
    }
    showNotification(sprintf("Saved '%s' for %s.", rate, scope_msg), type = "message")
  }, ignoreInit = TRUE)
  
  observeEvent(input$rate_click, {


    showNotification(
  "This public/reviewer version is read-only; ratings are not saved.",
  type = "message"
    )
    return(invisible())

    
    info <- input$rate_click
    gid <- as.character(info$gbifID %||% "")
    rate <- tolower(info$rating %||% "good")
    usr <- if (nzchar(input$display_name)) input$display_name else (current_user() %||% "")
    
    if (!nzchar(gid)) { showNotification("No GBIF ID captured.", type = "error"); return(invisible()) }
    if (!nzchar(usr)) { showNotification("Please enter your name before saving.", type = "error"); return(invisible()) }
    
    spec_val <- input$species %||% ""
    tribe_val <- taxonomy$tribe[taxonomy$species == spec_val][1]; tribe_val <- ifelse(is.na(tribe_val), "", tribe_val)
    ts_iso <- sub('(\\+|\\-)(\\d{2})(\\d{2})$', '\\1\\2:\\3',
                  format(lubridate::with_tz(lubridate::now(), TIME_TZ), '%Y-%m-%dT%H:%M:%S%z'))
    entry <- tibble::tibble(gbifID = gid, user = usr, rating = rate,
                            species = spec_val, tribe = tribe_val, timestamp_iso = ts_iso)
    
    if (is.null(ss)) ss <<- tryCatch(gs4_get(SHEET_ID), error = function(e) NULL)
    tr <- try({
      if (!is.null(ss)) sheet_append(ss, entry, sheet = SHEET_TAB)
      else sheet_append(SHEET_ID, entry, sheet = SHEET_TAB)
    }, silent = TRUE)
    
    if (inherits(tr, "try-error")) {
      msg <- conditionMessage(attr(tr, "condition"))
      message("[GS] sheet_append failed: ", msg)
      showNotification(paste("Could not save rating:", msg), type = "error", duration = 10)
      return(invisible())
    }
    
    showNotification("Rating saved.", type = "message")
    ratings_df(bind_rows(ratings_df(), entry))
  }, ignoreInit = TRUE)
  
  # ---- UI text: species summary (now with rating counts) ----
  output$species_summary <- renderUI({
    req(input$show_map)
    d <- species_data(); req(d)
    sp_value <- selected_species()
    sp_pretty <- pretty_species(sp_value)
    tribe <- taxonomy$tribe[taxonomy$species == sp_value][1]
    total <- nrow(d)
    gbif <- sum(tolower(d$coordinate_source) == "gbif", na.rm = TRUE)
    geoc <- sum(tolower(d$coordinate_source) != "gbif", na.rm = TRUE)
    k <- if ("cluster" %in% names(d)) length(unique(na.omit(d$cluster))) else 0
    
    lr <- latest_ratings_for(unique(as.character(d$gbifID)), ratings_df())
    # Count latest ratings
    latest_vec <- ifelse(is.na(lr$latest), "unrated", lr$latest)
    n_good <- sum(latest_vec == "good", na.rm = TRUE)
    n_neu  <- sum(latest_vec == "neutral", na.rm = TRUE)
    n_bad  <- sum(latest_vec == "bad", na.rm = TRUE)
    n_unr  <- max(total - (n_good + n_neu + n_bad), 0L)
    
    is_tol <- normalize_species(sp_pretty) %in% normalize_species(BRASSITOL_SPECIES)
    is_woody <- normalize_species(sp_pretty) %in% normalize_species(DERIVED_WOODY_SPECIES)
    phylo_html <- if (is_tol) sprintf(" <img src='%s' class='sel-ico' alt='phylo'/>", PHYLO_ICON_URL) else ""
    woody_html <- if (is_woody) sprintf(" <img src='%s' class='sel-ico' alt='woody'/>", WOODY_ICON_URL) else ""
    
    # Simple stacked bar to visualize rating distribution
    w_good <- if (total > 0) round(100 * n_good / total, 1) else 0
    w_neu  <- if (total > 0) round(100 * n_neu  / total, 1) else 0
    w_bad  <- if (total > 0) round(100 * n_bad  / total, 1) else 0
    w_unr  <- max(0, 100 - w_good - w_neu - w_bad)
    
    HTML(glue("
      <div class='kpi'>
        <b>Species overview</b><br>
        <b>Species:</b> {sp_pretty}{phylo_html}{woody_html}<br>
        <b>Tribe:</b> {tribe}<br><br>
        <b>Total records:</b> {total}<br>
        <b>Good / Neutral / Bad / Unrated:</b> {n_good} / {n_neu} / {n_bad} / {n_unr}<br>
        <div class='bar'>
          <div style='width:{w_good}%; background:#4CAF50; float:left;'></div>
          <div style='width:{w_neu}%; background:#ADD8E6; float:left;'></div>
          <div style='width:{w_bad}%; background:#e53935; float:left;'></div>
          <div style='width:{w_unr}%; background:#d3d3d3; float:left;'></div>
        </div>
        <small>&nbsp;Green = good, Blue = neutral, Red = bad, Grey = unrated</small><br><br>
        <b>Records with coordinates from GBIF:</b> {gbif}<br>
        <b>Records with coordinates from geocoding:</b> {geoc}<br>
        <b>Geographic clusters:</b> {k}
      </div>
    "))
  })
  
  # ---- Help modal ----
  observeEvent(input$help_btn, {
    phylo_img <- sprintf("<img src='%s' class='sel-ico' alt='phylo'/>", PHYLO_ICON_URL)
    woody_img <- sprintf("<img src='%s' class='sel-ico' alt='woody'/>", WOODY_ICON_URL)
    
    showModal(modalDialog(
      title = "How this tool works",
      easyClose = TRUE, footer = modalButton("Close"), size = "l",
      HTML(paste0(
        "<p><b>Welcome!</b> This is the <i>Brassicaceae Woodiness</i> project taxonomic curation tool.</p>",
        "<p>We aim to assemble expert-curated occurrence datasets to study the ecological niches of as many Brassicaceae species as possible.</p>",
        
        "<p><b>Background (data sources, preprocessing & sampling)</b></p>",
        "<ul style='margin-top:-8px;'>",
        "<li><b>Sources:</b> We compiled records primarily from GBIF (herbaria, surveys, iNaturalist, etc.).</li>",
        "<li><b>Cleaning:</b> Automated filters were used to remove erroneous records: (1) outside the native range (POWO); (2) with clearly false locations (botanic gardens, city centers, herbaria, oceans).</li>",
        "<li><b>Geocoding:</b> For records without coordinates, coordinates were inferred automatically using reformatting of locality text (using the ChatGPT API) and subsequent geocoding (using the Google Maps API); reported confidence was translated into coordinate uncertainty.</li>",
        "<li><b>Sampling:</b> To keep curation manageable for species with many records, we selected up to 100 accurate, spatially spread records, using ≤10 geographic clusters and round-robin sampling.</li>",
        "</ul>",
        
        "<p><b>What we kindly ask you to do</b></p>",
        "<ul style='margin-top:-8px;'>",
        "<li>Select a tribe and species, then click <b>SHOW MAP</b>.</li>",
        "<li>Click points to rate them <b>good</b>, <b>neutral</b>, or <b>bad</b>. The latest rating wins and you can overwrite earlier ones.</li>",
        "<li>Icons in the species dropdown: ", phylo_img, " = included in BrassiToL; ", woody_img, " = derived woody species. These icons indicate priority species.</li>",
        "<li><b>Bulk rate strategy:</b> a practical workflow is to first bulk-rate <i>all</i> (or all GBIF) points as “good” if the overall pattern looks correct; then adjust by bulk-rating <i>geocoded</i> or specific <i>clusters</i>, and finally change individual outliers to “bad”.</li>",
        "<li>We don’t expect line-by-line auditing. Please balance speed and care: focus on flagging the clearly wrong records; detailed metadata are available in the popups and via GBIF links.</li>",
        "</ul>",
        
        "<p><b>Final notes</b></p>",
        "<ul style='margin-top:-8px;'>",
        "<li><b>Legend:</b> Marker <b>fill</b> = rating (good/neutral/bad/unrated); marker <b>shape</b> = source (circle = GBIF; triangle = geocoded). POWO native ranges and buffered outlines can be toggled in the layer control.</li>",
        "<li><b>Saving & transparency:</b> Ratings are written to a Google Sheet. Each record popup shows the last three ratings with names and timestamps.</li>",
        "<li>Questions or ideas? Reach us at <a href='mailto:", CONTACT_EM, "'>", CONTACT_EM, "</a>.</li>",
        "</ul>"
      ))
    ))
  })
}

# Run the app locally
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))

# ## STEP 5: UPLOAD R SHINY APP ----
# 
# # Connect to shinyapps.io
# rsconnect::setAccountInfo(name='kasperhendriks',
#                           token='81E54EAF5C914A9FFD68A858B0C84884',
#                           secret='w+K+K86UXAc74t6bmIwQ7uTqbBG3j7K9G06jNwVQ')
# 
# # Deploy
# rsconnect::deployApp("./")
# 

