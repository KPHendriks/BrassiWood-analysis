#!/usr/bin/env Rscript

# 2d.RecordsToNiche.R ------------------------------------------------------
# Goal:
#   Build species-level niche summaries for Brassicaceae from curated occurrence
#   records, using one retained occurrence per raster cell and summarising
#   climate across those retained cells.
#
# Overall approach:
#   1. Load the Brassicaceae species list from the frozen local publication CSV.
#   2. Load manual record ratings from a local ratings CSV.
#   3. Load optional expert occurrence records from an external CSV.
#   4. Load environmental rasters (WorldClim monthly precipitation, monthly
#      temperature, elevation, and monthly ET0 from Zomer v3).
#   5. For each species:
#        - read the curated occurrence CSV from the curation app,
#        - join ratings by gbifID,
#        - append expert records,
#        - assign raster cells,
#        - retain one best record per raster cell under the current filtering
#          rules,
#        - extract monthly environmental values at retained cells,
#        - calculate per-point climate metrics,
#        - summarise those metrics per species using median + MAD.
#   6. Write complete and tiered output tables for downstream analyses.
#   7. Produce simple QC plots to inspect data quantity, record quality,
#      geocoding dependence, and stability of the climate summaries.
#
# Study focus and scope:
#   This script now focuses on climate-derived niche summaries only. We do not
#   calculate or retain range-size proxies here anymore (e.g. number of
#   botanical countries, number of Level-3 regions, summed WGSRPD area).
#   This keeps the niche pipeline aligned with the main focus of the study:
#   growth form in relation to climate, especially drought and frost.
#
# Key output:
#   WP2_BrassiNiche/results_final/1_niche_modelling_preparation/
#     - species_niche_summary_tbl.{rds,csv}
#     - species_niche_summary_tbl_all_ok.{rds,csv}
#     - species_niche_summary_tbl_ge5.{rds,csv}
#     - species_niche_summary_tbl_ge10.{rds,csv}
#     - species_niche_summary_ge1.{rds,csv}
#     - species_niche_summary_ge5.{rds,csv}
#     - species_niche_summary_ge10.{rds,csv}
#
# Notes for future use:
#   - The filtering rules for occurrence selection are intentionally conservative.
#   - Unrated GBIF/expert records are treated as neutral.
#   - Geocoded records are excluded unless needed to reach the target number of
#     retained non-bad cells, and then only when sufficiently precise and rated
#     good/neutral.
#   - Climate summaries are based on retained raster cells, not all raw records.


# 1. Prepare paths, output directories, and packages -----------------------

wd <- "~/Google Drive/My Drive/Publications/2026_Hendriks_et_al_BrassiWood/BrassiWood"
setwd(wd)

out_dir <- "WP2_BrassiNiche/results_final/1_niche_modelling_preparation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

occ_dir <- "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data"
cache_dir <- "WP2_BrassiNiche/results_intermediate/niche_cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Optional expert occurrence records, external to the curation app
expert_csv <- "WP2_BrassiNiche/data/expert_occurrence_records/Brassicaceae_expert_occurrence_records.csv"

# Preferred local directory for environmental rasters
env_dir <- "WP2_BrassiNiche/data/environmental_data"

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(terra)
  library(sf)
  library(ggplot2)
  library(scales)
})

# If needed:
# install.packages("geodata")
suppressPackageStartupMessages(library(geodata))


# 2. Load species list from frozen local publication CSV -------------------

# Rationale:
#   The species list defines the study scope. We restrict to accepted species in
#   Brassicaceae and exclude infraspecific taxa here, so the niche summaries are
#   generated at the species level.

species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
stopifnot(file.exists(species_file))

species_details <- readr::read_csv(
  species_file,
  show_col_types = FALSE
)

needed_species_cols <- c(
  "FAMILY",
  "ACCEPTED",
  "SUBSP_VAR",
  "SPECIES_NAME_PRINT"
)

missing_species_cols <- setdiff(needed_species_cols, names(species_details))
if (length(missing_species_cols) > 0) {
  stop(
    "[ERROR] Missing columns in species input file: ",
    paste(missing_species_cols, collapse = ", ")
  )
}

species_brassicaceae <- species_details %>%
  dplyr::filter(
    FAMILY == "Brassicaceae",
    ACCEPTED == "Y",
    is.na(SUBSP_VAR) | stringr::str_squish(as.character(SUBSP_VAR)) == ""
  ) %>%
  dplyr::pull(SPECIES_NAME_PRINT) %>%
  unique()

message("[INFO] Species list loaded from local CSV: n=", length(species_brassicaceae))


# 3. Load record ratings from local CSV ------------------------------------

# Rationale:
#   Ratings represent manual quality control of occurrence records. When multiple
#   ratings exist for the same gbifID, we keep only the most recent one.
#
#   The public/reviewer version expects ratings to be supplied as a local CSV.
#   This avoids connecting to the private Google Sheet used during curation.

ratings_csv <- "WP2_BrassiNiche/data/ratings_publication.csv"

if (file.exists(ratings_csv)) {
  
  ratings_raw <- readr::read_csv(
    ratings_csv,
    show_col_types = FALSE
  )
  
  needed_rating_cols <- c("gbifID", "rating", "timestamp_iso")
  missing_rating_cols <- setdiff(needed_rating_cols, names(ratings_raw))
  
  if (length(missing_rating_cols) > 0) {
    stop(
      "[ERROR] Missing columns in ratings CSV: ",
      paste(missing_rating_cols, collapse = ", ")
    )
  }
  
  ratings_tbl <- ratings_raw %>%
    dplyr::transmute(
      gbifID = as.character(gbifID),
      rating = tolower(as.character(rating)),
      timestamp_iso = as.POSIXct(timestamp_iso, tz = "UTC")
    ) %>%
    dplyr::filter(!is.na(gbifID), gbifID != "") %>%
    dplyr::arrange(gbifID, timestamp_iso) %>%
    dplyr::group_by(gbifID) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gbifID, rating)
  
  message("[INFO] Ratings loaded from local CSV: n=", nrow(ratings_tbl))
  
} else {
  
  warning(
    "[WARN] No local ratings CSV found at: ", ratings_csv,
    ". Proceeding with all GBIF/expert records treated as unrated/neutral unless rated elsewhere."
  )
  
  ratings_tbl <- tibble::tibble(
    gbifID = character(),
    rating = character()
  )
}



# 4. Load optional expert occurrence records -------------------------------

# Rationale:
#   Expert records are appended before thinning so they can participate in the
#   same cell-level filtering and selection procedure as app records.

expert_tbl <- tibble()

if (file.exists(expert_csv)) {
  expert_tbl <- suppressMessages(readr::read_csv(expert_csv, show_col_types = FALSE)) %>%
    mutate(
      longitude = iconv(as.character(longitude), from = "", to = "UTF-8", sub = ""),
      latitude  = iconv(as.character(latitude),  from = "", to = "UTF-8", sub = ""),
      longitude = gsub("[^0-9\\.-]+", "", longitude),
      latitude  = gsub("[^0-9\\.-]+", "", latitude)
    ) %>%
    transmute(
      species = gsub(" ", "_", as.character(species)),
      source  = as.character(source),
      decimalLatitude  = as.numeric(latitude),
      decimalLongitude = as.numeric(longitude),
      date   = as.character(date),
      rating = tolower(as.character(rating))
    ) %>%
    filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
    mutate(
      coordinate_source = "expert",
      gbifID = paste0("expert_", row_number()),
      coordinateUncertaintyInMeters = NA_real_
    ) %>%
    select(
      gbifID, species, decimalLatitude, decimalLongitude,
      coordinate_source, coordinateUncertaintyInMeters, rating,
      source, date
    )
  
  message(
    "[EXPERT] Loaded expert records: n=", nrow(expert_tbl),
    " (species=", length(unique(expert_tbl$species)), ")"
  )
} else {
  message("[EXPERT] No expert CSV found at: ", expert_csv, " (skipping)")
}


# 5. Load or prepare environmental rasters --------------------------------

# Goal:
#   Use elevation plus monthly precipitation, temperature, and ET0 to calculate
#   climate metrics per retained occurrence cell.
#
# Rationale:
#   We do not pre-build large derived global rasters here. Instead, we extract
#   monthly base values at retained points and compute summary climate metrics
#   per point afterwards. This is simpler to maintain and keeps the workflow
#   transparent.

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

# 5.1 Elevation ------------------------------------------------------------

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

names(elev) <- "elev"
template_r <- elev

# 5.2 Monthly precipitation and temperature --------------------------------

et0_v3_dir <- file.path(env_dir, "Global-ET0_v3_monthly")

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

# WorldClim tavg is typically stored as °C * 10
tavg_c <- tavg / 10

# 5.3 Monthly ET0 ----------------------------------------------------------

if (file.exists(et0_cache)) {
  et0 <- terra::rast(et0_cache)
  message("[ENV] ET0: using cached multi-layer tif ", et0_cache)
} else {
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
  
  message("[ENV] ET0 months found: ", paste(months, collapse = ", "))
  
  if (!identical(months, 1:12)) {
    stop(
      "[ENV] ET0 monthly files incomplete.\nFound months: ",
      paste(months, collapse = ", "),
      "\nFiles:\n",
      paste(basename(et0_files), collapse = "\n")
    )
  }
  
  et0 <- terra::rast(et0_files)
  names(et0) <- paste0("ET0_", sprintf("%02d", 1:12))
  message("[ENV] ET0: loaded Zomer v3 monthly ET0 (01–12) from ", et0_v3_dir)
  
  terra::writeRaster(et0, et0_cache, overwrite = TRUE)
}

stopifnot(terra::nlyr(et0) == 12)


# 6. Define helper functions -----------------------------------------------

# Goal:
#   Keep all repeated logic in one place so the per-species workflow remains
#   easy to read and easy to update later.

# 6.1 Robust summary helpers ------------------------------------------------

mad_robust <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  stats::mad(x, center = stats::median(x), constant = 1, na.rm = TRUE)
}

# 6.2 Coordinate uncertainty helpers ---------------------------------------

get_unc_m <- function(x) {
  u <- suppressWarnings(as.numeric(x))
  ifelse(is.na(u), Inf, u)
}

coerce_unc_m <- function(x) {
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_chr[x_chr == ""] <- NA_character_
  x_chr <- gsub("[^0-9eE+\\-\\.]", "", x_chr)
  suppressWarnings(as.numeric(x_chr))
}

# 6.3 Cell-level record selection rules ------------------------------------

# Rationale:
#   Selection is done per raster cell, not per raw record. Explicitly bad records
#   are always excluded. Unrated GBIF/expert records are treated as neutral.
#   Geocoded records are only considered if needed and if sufficiently precise.

pick_best_per_cell3 <- function(df, allow_geocoding = FALSE, geo_unc_max_m = 50000) {
  
  df <- df %>%
    mutate(
      rating_raw = tolower(as.character(rating)),
      rating = case_when(
        rating_raw %in% c("good", "neutral", "bad") ~ rating_raw,
        is.na(rating_raw) & coordinate_source %in% c("GBIF", "expert") ~ "neutral",
        TRUE ~ NA_character_
      ),
      unc_m = get_unc_m(coordinateUncertaintyInMeters),
      rating_ord = case_when(
        rating == "good" ~ 0L,
        rating == "neutral" ~ 1L,
        is.na(rating) ~ 2L,
        TRUE ~ 3L
      ),
      src_ord = case_when(
        coordinate_source == "expert" ~ 0L,
        coordinate_source == "GBIF" ~ 1L,
        coordinate_source == "geocoding" ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    filter(!(rating == "bad"))
  
  if (!allow_geocoding) {
    df <- df %>% filter(coordinate_source != "geocoding")
  } else {
    df <- df %>%
      filter(
        coordinate_source != "geocoding" |
          (unc_m <= geo_unc_max_m & rating %in% c("good", "neutral"))
      )
  }
  
  df %>%
    arrange(cell_id, rating_ord, src_ord, unc_m) %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup()
}

enforce_rules3 <- function(df, min_nonbad_cells = 20, geo_unc_max_m = 50000) {
  best_no_geo <- pick_best_per_cell3(
    df,
    allow_geocoding = FALSE,
    geo_unc_max_m = geo_unc_max_m
  )
  n_nonbad_no_geo <- nrow(best_no_geo)
  
  if (n_nonbad_no_geo >= min_nonbad_cells) {
    return(list(
      kept = best_no_geo,
      rule = "no_geocoding",
      n_cells_nonbad = n_nonbad_no_geo,
      n_cells_good = sum(best_no_geo$rating == "good", na.rm = TRUE),
      n_cells_neutral = sum(best_no_geo$rating == "neutral", na.rm = TRUE),
      n_geocoding_kept = sum(best_no_geo$coordinate_source == "geocoding", na.rm = TRUE)
    ))
  }
  
  best_with_geo <- pick_best_per_cell3(
    df,
    allow_geocoding = TRUE,
    geo_unc_max_m = geo_unc_max_m
  )
  n_nonbad_with_geo <- nrow(best_with_geo)
  
  list(
    kept = best_with_geo,
    rule = "allow_geocoding_if_needed_(rated_good_or_neutral_and_unc<=50km)",
    n_cells_nonbad = n_nonbad_with_geo,
    n_cells_good = sum(best_with_geo$rating == "good", na.rm = TRUE),
    n_cells_neutral = sum(best_with_geo$rating == "neutral", na.rm = TRUE),
    n_geocoding_kept = sum(best_with_geo$coordinate_source == "geocoding", na.rm = TRUE)
  )
}


# 7. Define per-species niche summarisation --------------------------------

# Goal:
#   Encapsulate the full workflow for a single species: load records, apply
#   selection rules, extract environmental values, compute climate metrics, and
#   return one summarised row.

process_one_species <- function(sp_print, ratings_tbl, expert_tbl,
                                template_r, elev, prec, et0, tavg_c,
                                occ_dir,
                                min_good_cells = 20,
                                min_cells_to_compute = 1,
                                low_n_threshold = 10,
                                geo_unc_max_m = 50000) {
  
  sp_file <- gsub(" ", "_", sp_print)
  f <- file.path(occ_dir, paste0(sp_file, ".csv"))
  
  if (!file.exists(f)) {
    return(tibble(species = sp_file, status = "missing_occurrence_csv"))
  }
  
  occ <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  
  need <- c("gbifID", "decimalLatitude", "decimalLongitude", "coordinate_source")
  miss <- setdiff(need, names(occ))
  
  if (length(miss) > 0) {
    return(
      tibble(
        species = sp_file,
        status = paste0("missing_columns: ", paste(miss, collapse = ", "))
      )
    )
  }
  
  if (!("coordinateUncertaintyInMeters" %in% names(occ))) {
    occ$coordinateUncertaintyInMeters <- NA_real_
  }
  
  occ$coordinateUncertaintyInMeters <- coerce_unc_m(occ$coordinateUncertaintyInMeters)
  
  occ <- occ %>%
    mutate(
      gbifID = as.character(gbifID),
      decimalLatitude  = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      coordinate_source = as.character(coordinate_source)
    ) %>%
    filter(is.finite(decimalLatitude), is.finite(decimalLongitude)) %>%
    mutate(species = sp_file)
  
  occ2 <- occ %>%
    left_join(ratings_tbl, by = "gbifID") %>%
    mutate(
      coordinateUncertaintyInMeters = coerce_unc_m(coordinateUncertaintyInMeters)
    )
  
  exp_sp <- expert_tbl %>% filter(species == sp_file)
  
  if (nrow(exp_sp) > 0) {
    exp_sp2 <- exp_sp %>%
      mutate(
        species = sp_file,
        gbifID  = as.character(gbifID),
        decimalLatitude  = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        coordinate_source = "expert",
        coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
        rating = tolower(as.character(rating))
      )
    
    occ2 <- occ2 %>%
      mutate(coordinateUncertaintyInMeters = coerce_unc_m(coordinateUncertaintyInMeters))
    
    common_cols <- union(names(occ2), names(exp_sp2))
    
    occ2 <- dplyr::bind_rows(
      dplyr::select(occ2, dplyr::any_of(common_cols)),
      dplyr::select(exp_sp2, dplyr::any_of(common_cols))
    )
  }
  
  if (nrow(occ2) < 1) {
    return(tibble(species = sp_file, status = "too_few_coords", n_rows = 0L))
  }
  
  xy <- cbind(occ2$decimalLongitude, occ2$decimalLatitude)
  occ2$cell_id <- terra::cellFromXY(template_r, xy)
  occ2 <- occ2 %>% filter(!is.na(cell_id))
  
  if (nrow(occ2) < 1) {
    return(tibble(species = sp_file, status = "too_few_after_cell_assign", n_rows = 0L))
  }
  
  sel <- enforce_rules3(
    occ2,
    min_nonbad_cells = min_good_cells,
    geo_unc_max_m = geo_unc_max_m
  )
  kept <- sel$kept
  
  if (nrow(kept) < min_cells_to_compute) {
    return(tibble(
      species = sp_file,
      status = "too_few_after_thinning",
      n_before = nrow(occ2),
      n_after = nrow(kept),
      rule = sel$rule,
      n_cells_nonbad = sel$n_cells_nonbad,
      n_cells_good = sel$n_cells_good,
      n_cells_neutral = sel$n_cells_neutral,
      n_geocoding_kept = sel$n_geocoding_kept
    ))
  }
  
  pts <- terra::vect(
    kept,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  )
  
  elev_vals <- terra::extract(elev, pts)[, -1, drop = FALSE]
  prec_vals <- terra::extract(prec, pts)[, -1, drop = FALSE]
  et0_vals  <- terra::extract(et0, pts)[, -1, drop = FALSE]
  tavg_vals <- terra::extract(tavg_c, pts)[, -1, drop = FALSE]
  
  P  <- as.matrix(prec_vals)
  E0 <- as.matrix(et0_vals)
  T  <- as.matrix(tavg_vals)
  Z  <- elev_vals[[1]]
  
  # Per-point climate metrics
  Pann  <- rowSums(P,  na.rm = TRUE)
  ETann <- rowSums(E0, na.rm = TRUE)
  AI <- ifelse(is.finite(ETann) & ETann > 0, Pann / ETann, NA_real_)
  
  def <- pmin(P - E0, 0)
  MCWD <- apply(def, 1, function(x) min(cumsum(x)))
  
  BIO15 <- apply(P, 1, function(v) {
    if (any(!is.finite(v))) return(NA_real_)
    m <- mean(v)
    if (!is.finite(m) || m == 0) return(NA_real_)
    100 * stats::sd(v) / m
  })
  
  BIO17 <- apply(P, 1, function(v) {
    if (any(!is.finite(v))) return(NA_real_)
    q <- sapply(1:12, function(i) sum(v[((i:(i + 2) - 1) %% 12) + 1]))
    min(q)
  })
  
  BIO10 <- apply(T, 1, function(v) {
    if (any(!is.finite(v))) return(NA_real_)
    q <- sapply(1:12, function(i) mean(v[((i:(i + 2) - 1) %% 12) + 1]))
    max(q)
  })
  
  BIO6 <- apply(T, 1, function(v) {
    if (any(!is.finite(v))) return(NA_real_)
    min(v)
  })
  
  kept3 <- kept %>%
    mutate(
      elev = Z,
      AI_P_over_ET0 = AI,
      MCWD = MCWD,
      BIO15 = BIO15,
      BIO17 = BIO17,
      BIO10 = BIO10,
      BIO6 = BIO6
    )
  
  low_n_flag <- nrow(kept3) < low_n_threshold
  
  tibble(
    species = sp_file,
    status = "ok",
    rule = sel$rule,
    low_n_flag = low_n_flag,
    n_before = nrow(occ2),
    n_cells_kept = nrow(kept3),
    n_cells_good = sum(kept3$rating == "good", na.rm = TRUE),
    n_cells_neutral = sum(kept3$rating == "neutral", na.rm = TRUE),
    n_cells_unrated = sum(is.na(kept3$rating)),
    n_geocoding_kept = sum(kept3$coordinate_source == "geocoding", na.rm = TRUE),
    n_expert_kept = sum(kept3$coordinate_source == "expert", na.rm = TRUE),
    
    elev_med = median(kept3$elev, na.rm = TRUE),
    elev_mad = mad_robust(kept3$elev),
    
    MCWD_med = median(kept3$MCWD, na.rm = TRUE),
    MCWD_mad = mad_robust(kept3$MCWD),
    
    AI_med = median(kept3$AI_P_over_ET0, na.rm = TRUE),
    AI_mad = mad_robust(kept3$AI_P_over_ET0),
    
    BIO15_med = median(kept3$BIO15, na.rm = TRUE),
    BIO15_mad = mad_robust(kept3$BIO15),
    
    BIO17_med = median(kept3$BIO17, na.rm = TRUE),
    BIO17_mad = mad_robust(kept3$BIO17),
    
    BIO10_med = median(kept3$BIO10, na.rm = TRUE),
    BIO10_mad = mad_robust(kept3$BIO10),
    
    BIO6_med = median(kept3$BIO6, na.rm = TRUE),
    BIO6_mad = mad_robust(kept3$BIO6)
  )
}


# 8. Run niche summarisation across all species -----------------------------

# Goal:
#   Iterate over the species list, summarise niches species by species, and then
#   prepare complete and thresholded analysis tables.

species_vec <- gsub(" ", "_", species_brassicaceae)

summary_tbl <- purrr::map_dfr(
  seq_along(species_vec),
  function(i) {
    if (i %% 10 == 0) {
      message(sprintf(
        "[PROGRESS] %d / %d (%.1f%%)",
        i, length(species_vec), 100 * i / length(species_vec)
      ))
    }
    
    process_one_species(
      sp_print = species_vec[i],
      ratings_tbl = ratings_tbl,
      expert_tbl = expert_tbl,
      template_r = template_r,
      elev = elev,
      prec = prec,
      et0 = et0,
      tavg_c = tavg_c,
      occ_dir = occ_dir,
      min_good_cells = 20,
      min_cells_to_compute = 1,
      low_n_threshold = 10,
      geo_unc_max_m = 50000
    )
  }
)

summary_tbl_all <- summary_tbl %>%
  filter(status == "ok")

summary_tbl_ge5 <- summary_tbl_all %>%
  filter(n_cells_kept >= 5)

summary_tbl_ge10 <- summary_tbl_all %>%
  filter(n_cells_kept >= 10)


# 9. Save niche summary tables ---------------------------------------------

# Rationale:
#   We keep both the full table (including failure states) and filtered
#   analysis-ready tiers. This is useful for both reproducibility and downstream
#   modelling.

# 9.1 Full and intermediate summary tables ---------------------------------

saveRDS(summary_tbl,      file.path(out_dir, "species_niche_summary_tbl.rds"))
saveRDS(summary_tbl_all,  file.path(out_dir, "species_niche_summary_tbl_all_ok.rds"))
saveRDS(summary_tbl_ge5,  file.path(out_dir, "species_niche_summary_tbl_ge5.rds"))
saveRDS(summary_tbl_ge10, file.path(out_dir, "species_niche_summary_tbl_ge10.rds"))

readr::write_csv(summary_tbl,      file.path(out_dir, "species_niche_summary_tbl.csv"))
readr::write_csv(summary_tbl_all,  file.path(out_dir, "species_niche_summary_tbl_all_ok.csv"))
readr::write_csv(summary_tbl_ge5,  file.path(out_dir, "species_niche_summary_tbl_ge5.csv"))
readr::write_csv(summary_tbl_ge10, file.path(out_dir, "species_niche_summary_tbl_ge10.csv"))

# 9.2 Final analysis tiers -------------------------------------------------

summary_ge1 <- summary_tbl %>%
  filter(status == "ok", n_cells_kept >= 1)

summary_ge5 <- summary_tbl %>%
  filter(status == "ok", n_cells_kept >= 5)

summary_ge10 <- summary_tbl %>%
  filter(status == "ok", n_cells_kept >= 10)

saveRDS(summary_ge1,  file.path(out_dir, "species_niche_summary_ge1.rds"))
saveRDS(summary_ge5,  file.path(out_dir, "species_niche_summary_ge5.rds"))
saveRDS(summary_ge10, file.path(out_dir, "species_niche_summary_ge10.rds"))

write_csv(summary_ge1,  file.path(out_dir, "species_niche_summary_ge1.csv"))
write_csv(summary_ge5,  file.path(out_dir, "species_niche_summary_ge5.csv"))
write_csv(summary_ge10, file.path(out_dir, "species_niche_summary_ge10.csv"))


# 10. Print final summaries to console -------------------------------------

message("[DONE] Status counts (all species):")
print(summary_tbl %>% count(status))

message("\n[DONE] Status counts (species with ≥1 raster cell retained):")
print(
  summary_tbl %>%
    filter(status == "ok", n_cells_kept >= 1) %>%
    count(status)
)

message("\n[DONE] Status counts (species with ≥5 raster cells retained):")
print(
  summary_tbl %>%
    filter(status == "ok", n_cells_kept >= 5) %>%
    count(status)
)

message("\n[DONE] Status counts (species with ≥10 raster cells retained):")
print(
  summary_tbl %>%
    filter(status == "ok", n_cells_kept >= 10) %>%
    count(status)
)

message("\n[DONE] Detailed breakdown:")

summary_breakdown <- tibble(
  category = c(
    "All species processed",
    "OK species (≥1 cell)",
    "OK species (≥5 cells)",
    "OK species (≥10 cells)",
    "Low-N OK species (<10 cells)",
    "Too few after thinning",
    "Too few coords",
    "Missing occurrence CSV"
  ),
  n = c(
    nrow(summary_tbl),
    sum(summary_tbl$status == "ok" & summary_tbl$n_cells_kept >= 1, na.rm = TRUE),
    sum(summary_tbl$status == "ok" & summary_tbl$n_cells_kept >= 5,  na.rm = TRUE),
    sum(summary_tbl$status == "ok" & summary_tbl$n_cells_kept >= 10, na.rm = TRUE),
    sum(summary_tbl$status == "ok" & summary_tbl$n_cells_kept < 10,  na.rm = TRUE),
    sum(summary_tbl$status == "too_few_after_thinning"),
    sum(summary_tbl$status == "too_few_coords"),
    sum(summary_tbl$status == "missing_occurrence_csv")
  )
)

print(summary_breakdown)


# 11. Create diagnostic plots for niche data quality -----------------------

# Goal:
#   Provide quick QC figures to inspect whether niche summaries are based on
#   adequate sampling, whether geocoding is disproportionately used in low-N
#   species, and whether the main climate summaries behave sensibly.
#
# Note:
#   Range-size diagnostics have been removed together with all range-size
#   variables from this script.

# 11.1 Load summary tables for plotting ------------------------------------

summary_all  <- read_csv(
  file.path(out_dir, "species_niche_summary_tbl_all_ok.csv"),
  show_col_types = FALSE
)
summary_ge5  <- read_csv(
  file.path(out_dir, "species_niche_summary_tbl_ge5.csv"),
  show_col_types = FALSE
)
summary_ge10 <- read_csv(
  file.path(out_dir, "species_niche_summary_tbl_ge10.csv"),
  show_col_types = FALSE
)

save_plot <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

interpretations <- c(
  "Plot_11_1_sampling_depth_by_tier: Histograms of retained raster cells for ge1+, ge5+, ge10+. Expect right tail; tier differences show how much data are lost by stricter thresholds.",
  "Plot_11_2_geocoding_vs_sampling_by_tier: n_geocoding_kept vs n_cells_kept, coloured by tier. Desired: most points low y; geocoding concentrated at low N (mostly ge1+).",
  "Plot_11_3_fraction_good_vs_sampling_by_tier: frac_good vs N, coloured by tier. Desired: well-sampled species not all near 0; low-N scatter expected.",
  "Plot_11_4_geocoding_fraction_by_Nclass_lt10: Boxplot of frac_geo for <10 vs ≥10. Desired: ≥10 generally low.",
  "Plot_11_4b_geocoding_fraction_by_Nclass_lt5: Boxplot of frac_geo for <5 vs ≥5. Desired: ≥5 lower than <5.",
  "Plot_11_5_AI_vs_Mcwd_by_tier: AI vs MCWD coloured by tier. Desired: coherent dry↔wet trend in all tiers; major tier separation could suggest low-N artefacts.",
  "Plot_11_6_Mcwd_vs_sampling_by_tier: MCWD median vs N with smooth per tier. Desired: stabilization after ~10; ge5 should already look more stable than ge1.",
  "Plot_11_7_McwdMAD_vs_sampling_by_tier: MCWD MAD vs N with smooth per tier. Desired: noisy at low N, stabilizing with more cells."
)

writeLines(
  interpretations,
  con = file.path(out_dir, "Plot_interpretations_11_niche_QC.txt")
)

summary_tiers <- bind_rows(
  summary_all  %>% mutate(tier = "ge1+"),
  summary_ge5  %>% mutate(tier = "ge5+"),
  summary_ge10 %>% mutate(tier = "ge10+")
) %>%
  mutate(
    tier = factor(tier, levels = c("ge1+", "ge5+", "ge10+")),
    frac_good = ifelse(n_cells_kept > 0, n_cells_good / n_cells_kept, NA_real_),
    frac_geo  = ifelse(n_cells_kept > 0, n_geocoding_kept / n_cells_kept, NA_real_)
  )

# 11.2 Plot 1: sampling depth by tier --------------------------------------

p1 <- summary_tiers %>%
  ggplot(aes(n_cells_kept)) +
  geom_histogram(binwidth = 2, color = "white") +
  facet_wrap(~ tier, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 5,  linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  labs(
    x = "Number of raster cells retained",
    y = "Number of species",
    title = "11.1 Sampling depth per species (tiers ge1/ge5/ge10)",
    subtitle = "Dashed lines: 5 and 10 cells. Compare tier-specific sampling distributions."
  ) +
  theme_minimal()

save_plot(p1, "Plot_11_1_sampling_depth_by_tier", w = 10, h = 8)

# 11.3 Plot 2: geocoding usage vs sample size ------------------------------

p2 <- summary_tiers %>%
  ggplot(aes(n_cells_kept, n_geocoding_kept, colour = tier)) +
  geom_point(alpha = 0.35) +
  labs(
    x = "Total cells retained",
    y = "Cells from geocoding",
    title = "11.2 Reliance on geocoded records (by tier)",
    subtitle = "Desired: geocoding concentrated at low N; higher tiers should show reduced geocoding reliance."
  ) +
  theme_minimal()

save_plot(p2, "Plot_11_2_geocoding_vs_sampling_by_tier", w = 10, h = 5)

# 11.4 Plot 3: fraction good-rated vs sample size --------------------------

p3 <- summary_tiers %>%
  ggplot(aes(n_cells_kept, frac_good, colour = tier)) +
  geom_point(alpha = 0.35) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    x = "Cells retained",
    y = "Fraction GOOD-rated",
    title = "11.3 Data quality vs sampling depth (by tier)",
    subtitle = "Desired: well-sampled species not all near 0; low-N variability is expected."
  ) +
  theme_minimal()

save_plot(p3, "Plot_11_3_fraction_good_vs_sampling_by_tier", w = 10, h = 5)

# 11.5 Plot 4: geocoding fraction by sampling class (<10 vs ≥10) ----------

p4 <- summary_all %>%
  mutate(
    N_class  = ifelse(low_n_flag, "<10 cells", "≥10 cells"),
    frac_geo = ifelse(n_cells_kept > 0, n_geocoding_kept / n_cells_kept, NA_real_)
  ) %>%
  ggplot(aes(N_class, frac_geo)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    x = "Sampling class",
    y = "Fraction geocoded",
    title = "11.4 Geocoding dependence by sampling depth (<10 vs ≥10)",
    subtitle = "Desired: <10 often higher is OK; ≥10 should generally be low."
  ) +
  theme_minimal()

save_plot(p4, "Plot_11_4_geocoding_fraction_by_Nclass_lt10", w = 10, h = 5)

# 11.6 Plot 4b: geocoding fraction by sampling class (<5 vs ≥5) -----------

p4b <- summary_all %>%
  mutate(
    N_class  = ifelse(n_cells_kept < 5, "<5 cells", "≥5 cells"),
    frac_geo = ifelse(n_cells_kept > 0, n_geocoding_kept / n_cells_kept, NA_real_)
  ) %>%
  ggplot(aes(N_class, frac_geo)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    x = "Sampling class",
    y = "Fraction geocoded",
    title = "11.4b Geocoding dependence (<5 vs ≥5)",
    subtitle = "Desired: ≥5 should show lower geocoding dependence than <5."
  ) +
  theme_minimal()

save_plot(p4b, "Plot_11_4b_geocoding_fraction_by_Nclass_lt5", w = 10, h = 5)

# 11.7 Plot 5: climate plausibility AI vs MCWD -----------------------------

p5 <- summary_tiers %>%
  ggplot(aes(AI_med, MCWD_med, colour = tier)) +
  geom_point(alpha = 0.35) +
  labs(
    x = "Aridity Index (P / ET0)",
    y = "MCWD (most negative running deficit)",
    title = "11.5 Climate plausibility check: AI vs MCWD (by tier)",
    subtitle = "Desired: coherent dry↔wet trend. Tier separation may indicate low-N artefacts."
  ) +
  theme_minimal()

save_plot(p5, "Plot_11_5_AI_vs_Mcwd_by_tier", w = 10, h = 5)

# 11.8 Plot 6: sensitivity of MCWD median to sample size -------------------

p6 <- summary_tiers %>%
  ggplot(aes(n_cells_kept, MCWD_med, colour = tier)) +
  geom_point(alpha = 0.25) +
  geom_smooth(se = FALSE) +
  labs(
    x = "Cells retained",
    y = "Median MCWD",
    title = "11.6 Sampling sensitivity: MCWD median vs N (by tier)",
    subtitle = "Desired: noise at low N; stabilization with higher N. ge5/ge10 should look progressively more stable."
  ) +
  theme_minimal()

save_plot(p6, "Plot_11_6_Mcwd_vs_sampling_by_tier", w = 10, h = 5)

# 11.9 Plot 7: sensitivity of MCWD MAD to sample size ----------------------

p7 <- summary_tiers %>%
  ggplot(aes(n_cells_kept, MCWD_mad, colour = tier)) +
  geom_point(alpha = 0.25) +
  geom_smooth(se = FALSE) +
  labs(
    x = "Cells retained",
    y = "MCWD MAD (niche breadth proxy)",
    title = "11.7 Niche breadth vs sampling depth (by tier)",
    subtitle = "Desired: breadth is noisy at low N and stabilizes with more cells."
  ) +
  theme_minimal()

save_plot(p7, "Plot_11_7_McwdMAD_vs_sampling_by_tier", w = 10, h = 5)

message("[DONE] Niche summaries and QC plots written to: ", out_dir)









