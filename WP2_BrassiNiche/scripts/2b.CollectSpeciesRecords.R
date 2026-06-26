#In this script, we will collect and clean species occurrence records for as many species of Brassicaceae as possible.
#These data will serve as input to define species niches later.


## STEP 0: PREPARE ----

# Install the following packages (only first time).
# install.packages(c("rWCVP", "sf", "ggplot2", "ggnewscale", "dplyr", "ggforce", "viridis"))
# install.packages("CoordinateCleaner")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")

# Load required packages.
library(ggplot2)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(rgbif)
library(sf)
library(CoordinateCleaner)
library(ggspatial)
library(geosphere)
library(countrycode)
library(httr)
library(jsonlite)
library(googleway)
library(readr)
library(purrr)

# Define the following user settings
buffer_zone <- 100 # set a bufferzone around the native range countries to allow records from just outside the native range
min_records_target <- 100
max_records_to_geocode <- 100 # set a maximum number of records to geocode using the OpenAI API key; using OpenAI for this task is not free!

# Refer to any necessary API keys needed to use ChatGPT and GoogleMaps APIs
# One option is to save these in the .Renviron, so they don't need to be written down in this script
google_api_key <- Sys.getenv("google_api_key_private")
google_api_key <- Sys.getenv("google_api_key_Naturalis")
openai_api_key <- Sys.getenv("openai_api_key")
# opencage_api_key <- Sys.getenv("opencage_api_key")

# Similarly, set the details for using the GBIF API
# An account on GBIF can be created on the website for free
gbif_user <- "khendriks"
gbif_pwd <- Sys.getenv("gbif_pwd")
gbif_email <- "kasper.hendriks@naturalis.nl"



## STEP 1: LOAD ALL SPECIES DATA ----

# Load the species data.

species_file <- "WP2_BrassiNiche/data/brassiwood_species_publication.csv"
stopifnot(file.exists(species_file))

species_details <- readr::read_csv(
  species_file,
  show_col_types = FALSE
)

needed_cols <- c(
  "ACCEPTED",
  "SPECIES_NAME_PRINT",
  "SUBFAMILY",
  "SUPERTRIBE",
  "TRIBE_FULL",
  "GENUS",
  "SUBSP_VAR",
  "WCVP_WGSRPD_LEVEL_3_native",
  "GROWTH_FORM",
  "gbif_matchType",
  "gbif_usageKey",
  "checked_KasperH_if_not_ACCEPTED"
)

missing_cols <- setdiff(needed_cols, names(species_details))
if (length(missing_cols) > 0) {
  stop("[ERROR] Missing columns in species input file: ",
       paste(missing_cols, collapse = ", "))
}

species_details_cleaned <- species_details %>%
  dplyr::select(dplyr::all_of(needed_cols)) %>%
  dplyr::filter(ACCEPTED == "Y") %>%
  dplyr::filter(is.na(SUBSP_VAR) | stringr::str_squish(as.character(SUBSP_VAR)) == "") %>%
  dplyr::mutate(
    GROWTH_FORM_binary = grepl("W", GROWTH_FORM),
    include = gbif_matchType == "EXACT" |
      grepl("Y", checked_KasperH_if_not_ACCEPTED),
    count_records_passed_tests_gbif = NA_integer_,
    count_records_passed_tests_geocoded = NA_integer_,
    count_records_passed_tests_total = NA_integer_
  ) %>%
  dplyr::arrange(SPECIES_NAME_PRINT)

# Check if there are any duplicated gbif_usageKey in the dataset (which could result for example from erroneously including synonyms in the database).
# If we find any duplicates, these need to be checked and corrected in the database prior to continuing this script!
# View(species_details_cleaned[which(duplicated(species_details_cleaned$gbif_usageKey)), ])


# # Create a list object to store downloaded records for this species
# records_list <- list()
#
# # Or, in case the below script has been run before, load any previously saved results first
# # This way, we can prevent to repeat the same requests from the GBIF database and from tools like OpenAI and GoogleMaps
# records_list <- readRDS(file = "WP2_BrassiNiche/results_intermediate/records_list.Rds")



## STEP 2: LOAD WGSRPD LEVEL 3 COUNTRY SHAPES ----

# Please, refer to the following documentation to learn more about the use of the
# World Geographical Scheme for Recording Plant Distributions (WGSRPD):
# https://en.wikipedia.org/wiki/World_Geographical_Scheme_for_Recording_Plant_Distributions
# and https://en.wikipedia.org/wiki/List_of_codes_used_in_the_World_Geographical_Scheme_for_Recording_Plant_Distributions

# Get the shape files following the WGSRPD Level 3 first from https://github.com/tdwg/wgsrpd/archive/refs/heads/master.zip
# (Same as in script "2a.IntroductionOverview.R")
wgsrpd_L3 <- st_read("WP2_BrassiNiche/data/wgsrpd-master/level3/level3.shp")
names(wgsrpd_L3)



## STEP 3: START SPECIES LOOP ----

# Loop through the species_details_cleaned dataframe.


for (row in sample(1:nrow(species_details_cleaned))) {
  # For later niche modelling, we target a total of 100 records per species, with 20 being the mininum.
  # For some species there will be plenty data on GBIF with (accurate) coordinates and we can select those first;
  # for other species, there will be few or no records with coordinates, and we shall try to infer coordinates from any textual description of the location.
  
  # Ta avoid the for loop to break, we will use tryCatch
  tryCatch({
    ### STEP 3.1: SELECT SPECIES ----
    
    # Get row number in the species_details_cleaned table for this usage_key
    # s <- which(species_details_cleaned$gbif_usageKey == usage_key)
    usage_key <- species_details_cleaned$gbif_usageKey[row]
    species <- species_details_cleaned$SPECIES_NAME_PRINT[row]
    tribe <- species_details_cleaned$TRIBE_FULL[row]
    
    # Check if records for this species have been collected before by checking the existence of expected output files
    if (file.exists(
      paste0(
        "WP2_BrassiNiche/results_intermediate/species_occurrence_records_tables/",
        gsub(" ", "_", species),
        ".csv"
      )
    ) &
    file.exists(
      paste0(
        "WP2_BrassiNiche/results_intermediate/species_occurrence_records_maps/",
        gsub(" ", "_", species),
        "_map.pdf"
      )
    )) {
      # This species has been handled before; print this to the screen
      message(
        paste0(
          "\n\n\nRecords for species ",
          species,
          " (tribe ",
          tribe,
          ") have been collected before; moving on to the next species."
        )
      )
      next # skip this species and move on to the next
    }
    
    # Else, print details that we will now collect records for this species.
    cat(
      "\n\n\n",
      "Now starting to collect records for species ",
      species,
      " (tribe ",
      tribe,
      ") from GBIF data.\n",
      sep = ""
    )
    
    
    
    ### STEP 3.2: DEFINE NATIVE RANGE ----
    
    # Collect and set native range information for this species based on the World Checklist of Vascular Plant data that we saved before
    native_range <- species_details_cleaned$WCVP_WGSRPD_LEVEL_3_native[row]
    
    # Filter the botanical countries boundaries object (wgsrpd_L3) based on `native_range`
    native_range_codes <- unlist(strsplit(native_range, ",\\s*"))  # Split native range into individual codes
    filtered_countries <- wgsrpd_L3 %>%
      filter(LEVEL3_COD %in% native_range_codes)  # Keep only countries that are in the native range
    
    # Reproject to a CRS with kilometers (e.g., UTM or a CRS that uses meters or kilometers)
    # Transform to UTM Zone 48S (EPSG:32748)
    filtered_countries_projected_32748 <- st_transform(filtered_countries, crs = 32748)
    
    # We allow record locations just outside the botanical country range and create a buffer of 100 kilometers (in the projected CRS, where units are in meters)
    buffered_countries <- st_buffer(filtered_countries_projected_32748$geometry,
                                    dist = 1000 * buffer_zone)  # 100 km
    
    # Reproject back to WGS 84 (degrees)
    native_range_buffered_wgs84 <- st_transform(buffered_countries, crs = 4326)
    
    # Do a quick visual check to see if the buffered range is what we want
    world <- ne_countries(scale = "medium", returnclass = "sf")
    ggplot() +
      geom_sf(data = world,
              fill = "lightgray",
              color = "white") +
      geom_sf(data = filtered_countries,
              fill = "blue",
              alpha = 0.5) +
      geom_sf(data = native_range_buffered_wgs84,
              fill = "green",
              alpha = 0.3) +
      theme_minimal() +
      labs(
        title = paste0("Native and Buffered Range of ", species, " (tribe ", tribe, ")"),
        subtitle = "Buffered range is 100km outside native range",
        caption = "Source: World Checklist of Vascular Plants (WCVP)"
      ) +
      coord_sf(
        xlim = c(-180, 180),
        ylim = c(-90, 90),
        expand = TRUE
      ) +
      theme(legend.position = "none")
    
    
    
    ## STEP 4: DEFINE CUSTOM FUNCTIONS ----
    
    # Several analyses will be run multiple times in the below script.
    # Therefore, it is easier to first define them as functions, that can later be called.
    
    # Function 1: Native range test
    
    # This function allows to check wether records are from the species native range,
    # thereby helping to remove any outliers and records of invasive specimens,
    # The input are the location's coordinates and the native range.
    
    
    # ##FOIR TESTSING
    # gbif_data_with_coord <- gbif_data[gbif_data$hasCoordinate ==
    #                                     TRUE, ]
    #   native_range <- native_range
    #
    #   native_range_test(gbif_data_with_coord = gbif_data[gbif_data$hasCoordinate ==
    #                                                        TRUE, ],
    #                     native_range = native_range,
    #                     buffer_zone = 100)
    #
    
    native_range_test <- function(gbif_data_with_coord,
                                  native_range,
                                  buffer_zone = 100) {
      # Collect and set native range information based on WCVP data
      native_range_codes <- unlist(strsplit(native_range, ",\\s*"))  # Split native range into individual codes
      filtered_countries <- wgsrpd_L3 %>%
        filter(LEVEL3_COD %in% native_range_codes)  # Keep only countries that are in the native range
      
      # Reproject to a CRS with kilometers (e.g., UTM or a CRS that uses meters or kilometers)
      filtered_countries_projected_32748 <- st_transform(filtered_countries, crs = 32748)
      
      # Create a buffer zone in meters (buffer_zone parameter is in kilometers, so multiply by 1000 to convert to meters)
      if (buffer_zone > 0) {
        message(paste(
          "Attempting to create buffered native range with",
          buffer_zone,
          "km..."
        ))
        buffered_countries <- st_buffer(filtered_countries_projected_32748$geometry,
                                        dist = buffer_zone * 1000)  # Buffer in meters
        
        # Reproject back to WGS 84 (degrees)
        native_range_buffered_wgs84 <- st_transform(buffered_countries, crs = 4326)
        
        # Filter out invalid geometries before proceeding
        native_range_buffered_wgs84 <- native_range_buffered_wgs84[st_is_valid(native_range_buffered_wgs84), ]
        
        # Try running the function with the buffered native range
        tryCatch({
          gbif_data_with_coord_sf <- st_as_sf(
            gbif_data_with_coord,
            coords = c("decimalLongitude", "decimalLatitude"),
            crs = 4326
          )
          
          gbif_data_with_coord_sf <- gbif_data_with_coord_sf %>%
            mutate(
              is_in_native_range = st_within(
                gbif_data_with_coord_sf,
                native_range_buffered_wgs84,
                sparse = FALSE
              )
            ) %>%
            mutate(is_in_native_range = rowSums(is_in_native_range) > 0)  # Convert matrix of logicals to a single TRUE/FALSE
          
          message("Function ran successfully with buffered native range.")
          return(gbif_data_with_coord_sf$is_in_native_range)
          
        }, error = function(e) {
          # If an error occurs with the buffered range, fallback to the unbuffered range
          message("Error occurred with buffered range: ", e$message)
          message("Running the function with the unbuffered native range...")
          
          # Create native range using unbuffered geometry (filtered_countries_projected_32748)
          filtered_countries_projected_32748 <- filtered_countries_projected_32748[st_is_valid(filtered_countries_projected_32748$geometry), ]
          
          gbif_data_with_coord_sf <- st_as_sf(
            gbif_data_with_coord,
            coords = c("decimalLongitude", "decimalLatitude"),
            crs = 4326
          )
          
          gbif_data_with_coord_sf <- gbif_data_with_coord_sf %>%
            mutate(
              is_in_native_range = st_within(
                gbif_data_with_coord_sf,
                filtered_countries_projected_32748,
                sparse = FALSE
              )
            ) %>%
            mutate(is_in_native_range = rowSums(is_in_native_range) > 0)  # Convert matrix of logicals to a single TRUE/FALSE
          
          message("Function ran successfully with unbuffered native range.")
          return(gbif_data_with_coord_sf$is_in_native_range)
        })
      } else {
        # If buffer_zone is 0 or negative, skip buffering and just use the original unbuffered geometry
        message("Skipping buffer zone creation, using original unbuffered native range...")
        
        # Filter out invalid geometries in the unbuffered range
        filtered_countries_projected_32748 <- filtered_countries_projected_32748[st_is_valid(filtered_countries_projected_32748$geometry), ]
        
        # Create native range using unbuffered geometry (filtered_countries_projected_32748)
        gbif_data_with_coord_sf <- st_as_sf(
          gbif_data_with_coord,
          coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326
        )
        
        gbif_data_with_coord_sf <- gbif_data_with_coord_sf %>%
          mutate(
            is_in_native_range = st_within(
              gbif_data_with_coord_sf,
              filtered_countries_projected_32748,
              sparse = FALSE
            )
          ) %>%
          mutate(is_in_native_range = rowSums(is_in_native_range) > 0)  # Convert matrix of logicals to a single TRUE/FALSE
        
        message("Function ran successfully with unbuffered native range.")
        return(gbif_data_with_coord_sf$is_in_native_range)
      }
    }
    
    # Function 2: CoordinateCleaner wrapper
    
    # The R package CoordinateCleaner helps flag any suspicious records (near botanical garden, museum, capital, etc.).
    # Here we define a wrapper function that will simply take as input the dataframe, assuming that we always
    # use the same column names in the dataframe to specify the coordinates and the species name,
    # and returning only the CoordinateCleaner output in new columns.
    # Note that we do not use the outliers test, since we already flagged records outside native range and often have species with few records
    
    coordinate_cleaner_wrapper <- function(gbif_data_with_coord) {
      # Run the clean_coordinates function
      result_coordinate_cleaner <- clean_coordinates(
        x = gbif_data_with_coord,
        lon = "decimalLongitude",
        lat = "decimalLatitude",
        species = "SPECIES_NAME_PRINT",
        tests = c(
          "capitals",
          "centroids",
          "equal",
          "gbif",
          "institutions",
          "seas",
          "zeros"
        )
      )
      
      # Select the relevant columns and rename them with a prefix "coordinate_cleaner_"
      result_coordinate_cleaner %>%
        dplyr::select(.val, .equ, .zer, .cap, .cen, .sea, .gbf, .inst, .summary) %>%
        dplyr::rename_with(~ paste0("coordinate_cleaner_", .))
    }
    
    # Function 3: Direction to bearing lookup
    direction_to_bearing <- function(direction) {
      bearings <- c(
        N = 0,
        NNE = 22.5,
        NE = 45,
        ENE = 67.5,
        E = 90,
        ESE = 112.5,
        SE = 135,
        SSE = 157.5,
        S = 180,
        SSW = 202.5,
        SW = 225,
        WSW = 247.5,
        W = 270,
        WNW = 292.5,
        NW = 315,
        NNW = 337.5
      )
      dir <- toupper(direction)
      if (!dir %in% names(bearings)) {
        warning(paste("Unrecognized direction:", direction))
        return(NA)
      }
      bearings[[dir]]
    }
    
    # Function 4: Heuristic accuracy estimate based on confidence label
    estimate_accuracy_km <- function(confidence) {
      switch(
        confidence,
        "geocode_failed" = NA,
        "country_only" = 250,
        "region_only" = 100,
        "low_specificity_feature" = 50,
        "moderate_geocode" = 15,
        "adjusted_with_direction" = 10,
        "exact_geocode" = 5,
        "uncertain" = 25,
        25
      )
    }
    
    # Function 5: Parse and geocode location_string
    parse_and_geocode_location <- function(location_string,
                                           openai_api_key,
                                           google_api_key,
                                           model = "gpt-3.5-turbo") {
      # Prompt construction
      prompt <- paste0(
        "You are an expert in parsing locality descriptions from biodiversity records.\n",
        "Always identify a known place suitable for geocoding.\n",
        "If a location is described relative to another (e.g. '10 km NW of X'), treat X as the reference_place and preserve the distance and direction.\n",
        "In case the French word 'lieues' is used for distance, interpret it as 4 km\n",
        "If no clearly geocodable place is mentioned, use the best nearby town as fallback_place.\n",
        "The reference_place should be the most specific geocodable location.\n",
        "Handle input strings that may be in a variety of languages, including English, French, German, Spanish, Portuguese, Italian, Latin, Russian, Turkish, Arabic, or Chinese. Translate or interpret local expressions or directional phrases when needed. For example, 'ob' in German often means 'above' or 'north of'.\n",
        "Make sure to translate directions to their abbreviations (e.g., north should be N); if a direction cannot be clearly inferred as one of the standard compass headings (e.g. N, NW, SE), leave it out entirely (do not include a direction key in the JSON at all). This helps avoid geocoding errors., leave it empty (null). This helps avoid geocoding errors.\n",
        "Make sure all backslashes inside string values are escaped as double backslashes (\\\\) to conform with JSON format.\n",
        "Return only this JSON format:\n",
        "{ \"reference_place\": \"...\", \"fallback_place\": \"...\", \"direction\": \"...\", \"distance_km\": ..., \"country\": \"...\" }\n",
        "Now parse:\n\"",
        location_string,
        "\""
      )
      
      openai_response <- POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers(Authorization = paste("Bearer", openai_api_key)),
        content_type_json(),
        body = list(
          model = model,
          messages = list(list(
            role = "user", content = prompt
          )),
          temperature = 0.3
        ),
        encode = "json"
      )
      
      parsed_full <- content(openai_response, as = "parsed", encoding = "UTF-8")
      content_string <- parsed_full$choices[[1]]$message$content
      content_string <- gsub("\\\\", "\\\\\\\\", content_string)  # Escape slashes
      parsed_json <- tryCatch({
        fromJSON(content_string)
      }, error = function(e) {
        warning(paste(
          "JSON parse failed for string:",
          location_string,
          ":",
          e$message
        ))
        return(NULL)
      })
      
      # # If parsing failed, return NULL early
      # if (is.null(parsed_json))
      #   return(NULL)
      
      # Geocode with Google API
      google_url <- "https://maps.googleapis.com/maps/api/geocode/json"
      place_query <- paste(
        parsed_json$reference_place,
        parsed_json$fallback_place,
        parsed_json$country
      )
      
      google_response <- GET(google_url,
                             query = list(address = place_query, key = google_api_key))
      google_content <- content(google_response, as = "parsed", encoding = "UTF-8")
      
      # Handle empty results or ZERO_RESULTS
      if (google_content$status == "ZERO_RESULTS" ||
          length(google_content$results) == 0) {
        stop(paste(
          "Google Geocoding returned ZERO_RESULTS for:",
          place_query
        ))
      }
      
      base_coords <- google_content$results[[1]]$geometry$location
      lat <- base_coords$lat
      lon <- base_coords$lng
      location_type <- google_content$results[[1]]$geometry$location_type
      
      # If directional offset is given, apply it
      if (!is.null(parsed_json$direction) &&
          !is.null(parsed_json$distance_km)) {
        bearing <- direction_to_bearing(parsed_json$direction)
        if (!is.na(bearing)) {
          offset_point <- destPoint(
            p = c(lon, lat),
            b = bearing,
            d = parsed_json$distance_km * 1000
          )
          lat <- offset_point[2]
          lon <- offset_point[1]
          confidence <- "adjusted_with_direction"
        } else {
          confidence <- "direction_uninterpretable"
        }
      } else {
        confidence <- "exact_geocode"
      }
      
      types <- google_content$results[[1]]$types
      if ("country" %in% types || "continent" %in% types) {
        confidence <- "country_only"
      } else if ("administrative_area_level_1" %in% types) {
        confidence <- "region_only"
      } else if (any(types %in% c("natural_feature", "establishment", "political"))) {
        confidence <- "low_specificity_feature"
      } else if (any(types %in% c("locality", "sublocality", "postal_code"))) {
        confidence <- "moderate_geocode"
      } else if (any(
        types %in% c(
          "street_address",
          "route",
          "intersection",
          "premise",
          "point_of_interest"
        )
      )) {
        confidence <- "exact_geocode"
      } else {
        confidence <- "uncertain"
      }
      
      accuracy_km <- estimate_accuracy_km(confidence)
      
      `%||%` <- function(a, b)
        if (length(a) == 0 || is.null(a))
          b
      else
        a
      
      geocoding_result <- data.frame(
        original_string = location_string,
        reference_place = parsed_json$reference_place %||% NA,
        fallback_place = parsed_json$fallback_place %||% NA,
        country = parsed_json$country %||% NA,
        lat = lat,
        lon = lon,
        location_types = paste(unlist(types), collapse = ", "),
        confidence = confidence,
        accuracy_km = accuracy_km,
        location_type = location_type,
        source = paste0("OpenAI API (", model, ") + Google Geocode API"),
        direction = parsed_json$direction %||% NA,
        distance_km = parsed_json$distance_km %||% NA
      )
      
      return(geocoding_result)
    }
    
    
    
    ## STEP 5: DOWNLOAD GBIF RECORDS ----
    
    # Check if the taxonKey is valid
    if (is.na(usage_key)) {
      message("\n\n\nUsage key for this species NA; skip and continue with next species.")
      next # skip this species and move on to the next
    }
    
    # Set a relatively high timeout to allow a slow response from the GBIF server
    httr::set_config(httr::timeout(60)) # 60 seconds
    
    # Download all GBIF records WITH coordinates for this species
    # Submit the download request
    download_key <- occ_download(
      pred("taxonKey", usage_key),
      user = gbif_user,
      pwd = gbif_pwd,
      email = gbif_email
    )
    
    # GBIF will later download a zip file with data; to prevent any mixing with previous files, we will
    # first remove any previous files
    download_gbif_zip_filaname <- paste0(download_key, ".zip")
    if (file.exists(download_gbif_zip_filaname)) {
      message("\n\n\nGBIF download zip file already exists. Removing it to avoid confusion.")
      file.remove(download_gbif_zip_filaname)
    }
    
    # Loop: check status until ready
    repeat {
      status <- occ_download_meta(download_key)$status
      cat("Status: ", status, "\n")
      
      if (status %in% c("SUCCEEDED", "KILLED", "CANCELLED"))
        break
      
      Sys.sleep(10)  # Wait 10 seconds before checking again
    }
    
    # If successful, download the file
    if (status == "SUCCEEDED") {
      occ_file <- occ_download_get(download_key, overwrite = TRUE)
      
      # Import as data.frame
      occ_data <- occ_download_import(occ_file)
      
      cat("Data imported. Number of records: ", nrow(occ_data), "\n")
      
      # occ_data is now your data.frame
    } else {
      cat("Download did not succeed. Status was: ", status, "\n")
    }
    
    # Store only the dataframe and from that only the relevant columns (there can be >150 columns!)
    gbif_data <- occ_data %>%
      select(any_of(
        c(
          "gbifID",
          "acceptedScientificName",
          "scientificName",
          "decimalLatitude",
          "decimalLongitude",
          "hasCoordinate",
          "issues",
          "coordinateUncertaintyInMeters",
          "elevation",
          "elevationAccuracy",
          "basisOfRecord",
          "eventDate",
          "recordedBy",
          "identifiedBy",
          "collectionID",
          "recordNumber",
          "references",
          "datasetName",
          "occurrenceID",
          "continent",
          "countryCode",
          "country",
          "county",
          "locality",
          "locationRemarks",
          "stateProvince",
          "municipality",
          "higherGeography",
          "verbatimLocality",
          "verbatimElevation",
          "verbatimCoordinateSystem",
          "verbatimSRS",
          "locationID"
        )
      ))
    
    # GBIF downloaded a zip file with data; to prevent any mixing with future files, we will remove it
    if (file.exists(download_gbif_zip_filaname)) {
      file.remove(download_gbif_zip_filaname)
    }
    
    # Translate any country codes to full country names for use in geocoding later
    # We use R package countrycode for this
    gbif_data$country_full_name <- countrycode(gbif_data$countryCode,
                                               origin = "iso2c",
                                               destination = "country.name")
    
    # We add taxonomic details that we can use later to filter the data
    gbif_data$SPECIES_NAME_PRINT <- species_details_cleaned$SPECIES_NAME_PRINT[row]
    gbif_data$TRIBE_FULL <- species_details_cleaned$TRIBE_FULL[row]
    gbif_data$gbif_usageKey <- usage_key
    
    # Add details of where the coordinates derive from
    gbif_data$coordinate_source <- ifelse(gbif_data$hasCoordinate, "GBIF", "geocoding")
    
    # Add the new columns to store results from two location tests and initialize them to NA
    gbif_data$native_range_test <- NA
    gbif_data <- gbif_data %>%
      mutate(
        coordinate_cleaner_.val = NA,
        coordinate_cleaner_.equ = NA,
        coordinate_cleaner_.zer = NA,
        coordinate_cleaner_.cap = NA,
        coordinate_cleaner_.cen = NA,
        coordinate_cleaner_.sea = NA,
        coordinate_cleaner_.gbf = NA,
        coordinate_cleaner_.inst = NA,
        coordinate_cleaner_.summary = NA
      )
    gbif_data$passed_all_tests <- NA
    
    # Perform two location tests using the predefined custom functions and add results to the table
    # Check if there are any records with coordinates
    if (sum(gbif_data$hasCoordinate) > 0) {
      gbif_data$native_range_test[which(gbif_data$hasCoordinate == TRUE)] <- native_range_test(
        gbif_data_with_coord = gbif_data[gbif_data$hasCoordinate ==
                                           TRUE, ],
        native_range = native_range,
        buffer_zone = 100
      )
      
      # Run the CoordinateCleaner wrapper and merge the results back into the dataframe
      cleaned_results <- coordinate_cleaner_wrapper(gbif_data_with_coord = gbif_data[gbif_data$hasCoordinate == TRUE, ])
      
      # Merge the results into the original dataframe
      gbif_data[gbif_data$hasCoordinate == TRUE, c(
        "coordinate_cleaner_.val",
        "coordinate_cleaner_.equ",
        "coordinate_cleaner_.zer",
        "coordinate_cleaner_.cap",
        "coordinate_cleaner_.cen",
        "coordinate_cleaner_.sea",
        "coordinate_cleaner_.gbf",
        "coordinate_cleaner_.inst",
        "coordinate_cleaner_.summary"
      )] <- cleaned_results
    }
    
    # Check what records passed both tests so far
    gbif_data$passed_all_tests[gbif_data$hasCoordinate == TRUE] <- gbif_data$native_range_test[gbif_data$hasCoordinate == TRUE] &
      gbif_data$coordinate_cleaner_.summary[gbif_data$hasCoordinate == TRUE]
    
    # Count the number of records that passed both tests; we aim for at least 20 records to safely define the niche later
    number_records_passed_gbif <- gbif_data %>%
      filter(hasCoordinate, passed_all_tests) %>%
      distinct(decimalLatitude, decimalLongitude) %>% # this helps to NOT count records with exact same coordinates, which are most likely duplicates from e.g. herbaria
      nrow()
    
    cat(
      "\n\n\nFor ",
      species,
      ", GBIF holds ",
      number_records_passed_gbif,
      " unique records with coordinates that pass both native range and CoordinateCleaner tests for this species\n",
      sep = ""
    )
    
    
    
    ## STEP 6: CONTINUE WITH GEOCODING IF NEEDED ----
    
    # In case GBIF did not have enough records with coordinates, we can try estimating the coordinates from a location description.
    # We want to use at least 20 records per species to later define the species niche, so in case we have <20 records passing
    # our tests in the previous step, we continue now and  try applying geocoding on the GBIF records that do not
    # have coordinates listed, but only have a textual description of their location.
    
    if (number_records_passed_gbif >= min_records_target) {
      # Print a message if we do not need to continue with geocoding
      message("We do not need to continue trying to geocode any records without coordinates.")
    } else {
      # Print a message if we will continue with geocoding
      message("We will continue trying to geocode any records without coordinates.")
      
      # Create a location string for each record for which no coordinates were reported from GBIF that we can use for geocoding
      
      # Define any GBIF columns that can have descriptive location details that could be useful
      cols_location <- c(
        "continent",
        "country",
        "country_full_name",
        "county",
        "locality",
        "locationRemarks",
        "stateProvince",
        "municipality",
        "higherGeography",
        "verbatimLocality",
        "verbatimElevation",
        "elevation",
        "verbatimCoordinateSystem",
        "verbatimSRS",
        "locationID"
      )
      
      # Keep only columns that actually exist
      cols_location_exist <- intersect(cols_location, colnames(gbif_data))
      
      # Create and add location_string and the count of non-empty columns to the gbif_data dataframe
      # The count of non-empty strings can later be used to prioritize records for geocoding, assuming
      # that records with richer location details are more likely to be successfully and accurately geocoded
      gbif_data <- gbif_data %>%
        mutate(
          location_string = ifelse(
            hasCoordinate == FALSE,
            apply(select(., all_of(
              cols_location_exist
            )), 1, function(row_vals) {
              vals <- as.character(row_vals)
              vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
              if (length(vals) > 0)
                paste(vals, collapse = ", ")
              else
                NA_character_
            }),
            NA_character_
          ),
          non_empty_location_count = ifelse(
            hasCoordinate == FALSE,
            apply(select(., all_of(
              cols_location_exist
            )), 1, function(row_vals) {
              vals <- as.character(row_vals)
              count <- sum(!is.na(vals) & nzchar(trimws(vals)))
              if (count > 0)
                count
              else
                NA_integer_
            }),
            NA_integer_
          )
        )
      
      # To avoid geocoding mutiple records with the exact same location_string, let's check for duplicates
      gbif_data$location_string_duplicated <- duplicated(gbif_data$location_string)
      
      
      # We can now loop through any records and try to geocode these using our custom functions
      # We start with the record with the richest location_string and then continue with records with less details
      
      # Define a vector of records that we can check, starting with the record with the richest location details;
      # ignore diplicated strings, as these will lead to the same geocoded location and thus yield no additional information
      records_to_geocode_ordered <- !is.na(gbif_data$non_empty_location_count) &
        !gbif_data$location_string_duplicated
      records_to_geocode_ordered <- which(records_to_geocode_ordered)[order(gbif_data$non_empty_location_count[records_to_geocode_ordered], decreasing = TRUE)]
      
      # Apply the set max of records to geocode in case there are more records available than this set maximum
      # This way, we can avoid excessive costs from using the OpenAI API key
      records_to_geocode_ordered_subset <- records_to_geocode_ordered[intersect(1:max_records_to_geocode,
                                                                                seq_along(records_to_geocode_ordered))]
      
      # Add a column to the dataframe that allows storing any error messages from our custom geocoding function
      gbif_data$geocoding_errors <- NA
      
      # Then loop using this vector of records
      for (r in records_to_geocode_ordered_subset) {
        message(
          "\n\n\nStarting geocoding for species ",
          species,
          ", gbif record usage key ",
          usage_key,
          ", which is record ",
          which(records_to_geocode_ordered_subset == r),
          " of ",
          length(records_to_geocode_ordered_subset)
        )
        
        # Check if this gbif record has not been geocoded before, in which case we can skip it
        gbif_id <- gbif_data$gbifID[r]
        
        if (!is.na(gbif_data$decimalLatitude[r]) &
            !is.na(gbif_data$decimalLongitude[r])) {
          message("This gbif record has been geocoded before and will be skipped")
          next
        }
        
        # Initialize variables to store the result and error message
        result_geocoding <- NULL
        error_message <- NULL
        coordinate_cleaner_geocoding_result <- NULL
        
        # Use tryCatch to handle potential errors
        result_geocoding <- tryCatch({
          # Try running the custom geocoding function
          parse_and_geocode_location(
            location_string = gbif_data$location_string[r],
            openai_api_key,
            google_api_key,
            model = "gpt-3.5-turbo"
          )
        }, error = function(e) {
          # If an error occurs, store the error message in the error_message variable
          error_message <- paste("Error at record", r, ":", e$message)
          return(NULL)  # Return NULL in case of error
        })
        
        # After the tryCatch, you can check if there was an error
        # Check if error occurred OR result is NULL OR missing key fields
        if (!is.null(error_message) ||
            is.null(result_geocoding) ||
            !is.data.frame(result_geocoding) ||
            !"lat" %in% names(result_geocoding) ||
            length(result_geocoding$lat) == 0 ||
            is.na(result_geocoding$lat)) {
          # Catch silent NULL or malformed returns
          if (is.null(error_message)) {
            error_message <- paste("Invalid or empty geocoding result for this record's location string")
          }
          message(error_message)
          gbif_data$geocoding_errors[r] <- error_message
        } else {
          # Use geocoding result safely
          print(result_geocoding)
          gbif_data$decimalLatitude[r] <- result_geocoding$lat
          gbif_data$decimalLongitude[r] <- result_geocoding$lon
          gbif_data$coordinateUncertaintyInMeters[r] <- 1000 * result_geocoding$accuracy_km
          
          gbif_data$native_range_test[r] <- native_range_test(
            gbif_data_with_coord = gbif_data[r, ],
            native_range = native_range,
            buffer_zone = buffer_zone
          )
          
          coordinate_cleaner_geocoding_result <- coordinate_cleaner_wrapper(gbif_data_with_coord = gbif_data[r, ])
          
          gbif_data[r, c(
            "coordinate_cleaner_.val",
            "coordinate_cleaner_.equ",
            "coordinate_cleaner_.zer",
            "coordinate_cleaner_.cap",
            "coordinate_cleaner_.cen",
            "coordinate_cleaner_.sea",
            "coordinate_cleaner_.gbf",
            "coordinate_cleaner_.inst",
            "coordinate_cleaner_.summary"
          )] <- coordinate_cleaner_geocoding_result
          
          gbif_data$passed_all_tests[r] <- gbif_data$native_range_test[r] &
            gbif_data$coordinate_cleaner_.summary[r]
        }
      }
    }
    
    
    # ### TEMP SOME LINE FOR TESTING SPECIFIC r VALUES FOR THIS SPECIES s
    # t(result_geocoding)
    # message(result_geocoding$lat, ",", result_geocoding$lon)
    # message(result_geocoding$reference_place, ",", result_geocoding$fallback_place, ",", result_geocoding$country)
    
    
    
    ## STEP 7: SAVE & PLOT RESULTS ----
    
    # Summarise results and add to species table
    species_details_cleaned$count_records_passed_tests_gbif[row] <- number_records_passed_gbif
    species_details_cleaned$count_records_passed_tests_geocoded[row] <- sum(gbif_data$coordinate_source ==
                                                                              "geocoding" &
                                                                              gbif_data$passed_all_tests,
                                                                            na.rm = TRUE)
    species_details_cleaned$count_records_passed_tests_total[row] <- species_details_cleaned$count_records_passed_tests_gbif[row] + species_details_cleaned$count_records_passed_tests_geocoded[row]
    
    # Save species table
    write.csv(
      gbif_data,
      file = paste0(
        "WP2_BrassiNiche/results_intermediate/species_occurrence_records_tables/",
        gsub(" ", "_", species),
        ".csv"
      )
    )
    
    # Later we will create an R Shiny app for expert curation of occurrence data
    # If the folder structure already exists, let's also save the same csv file directly to the R Shiny app directory
    if (dir.exists("WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data")) {
      write.csv(
        gbif_data,
        file = paste0(
          "WP2_BrassiNiche/Brassicaceae_taxonomic_curation_app/occurrence_data/",
          gsub(" ", "_", species),
          ".csv"
        )
      )
    }
    
    
    
    saveRDS(object = species_details_cleaned, file = "WP2_BrassiNiche/results_intermediate/species_details_cleaned.Rds")
    
    # # Save the results for this species to our records_list list object
    # records_list$gbif[[species]] <- gbif_data
    
    # # Store locally for further use in niche modelling
    # saveRDS(object = records_list, file = "WP2_BrassiNiche/results_intermediate/records_list_backup2.Rds")
    
    # Create a species plot from these results
    
    # Create an sf object for records that have coordinates
    gbif_points <- gbif_data %>%
      filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
      st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
               crs = 4326)
    
    # Add categories for plotting
    gbif_points$plot_category <- case_when(
      gbif_points$coordinate_source == "GBIF" &
        gbif_points$passed_all_tests ~ "GBIF - passed",
      gbif_points$coordinate_source == "GBIF" &
        !gbif_points$passed_all_tests ~ "GBIF - failed",
      gbif_points$coordinate_source == "geocoding" &
        gbif_points$passed_all_tests ~ "Geocoded - passed",
      gbif_points$coordinate_source == "geocoding" &
        !gbif_points$passed_all_tests ~ "Geocoded - failed"
    )
    
    # Define color palette
    plot_colors <- c(
      "GBIF - passed" = "darkgreen",
      "GBIF - failed" = "grey50",
      "Geocoded - passed" = "orange",
      "Geocoded - failed" = "red"
    )
    
    # Plot and save
    # Plot and save with country labels
    p <- ggplot() +
      # Plot the native range countries
      geom_sf(data = filtered_countries,
              fill = "blue",
              alpha = 0.3) +
      
      # Plot the native range buffer
      geom_sf(data = native_range_buffered_wgs84,
              fill = "lightgreen",
              alpha = 0.2) +
      
      # Plot the records
      geom_sf(
        data = gbif_points,
        aes(color = plot_category),
        size = 2,
        alpha = 0.7
      ) +
      
      # Plot the country borders outside the native range
      geom_sf(
        data = wgsrpd_L3,
        color = "black",
        fill = NA,
        size = 0.5
      ) +
      
      # Add labels for the countries outside the native range
      geom_sf_text(
        data = wgsrpd_L3,
        aes(label = LEVEL3_NAM),
        size = 2,
        color = "black",
        fontface = "italic",
        check_overlap = TRUE
      ) +
      
      # Apply the color palette to the plot categories
      scale_color_manual(values = plot_colors, name = "Record type") +
      
      coord_sf(
        xlim = st_bbox(native_range_buffered_wgs84)[c("xmin", "xmax")],
        ylim = st_bbox(native_range_buffered_wgs84)[c("ymin", "ymax")],
        expand = TRUE
      ) +
      
      labs(
        title = paste0("Occurrence records of ", species, " (tribe ", tribe, ")"),
        subtitle = paste0(
          "Colored by source and test results; ",
          species_details_cleaned$count_records_passed_tests_total[row],
          " records passed the location tests"
        ),
        caption = "Blue = native range; green = native range with buffer",
        x = "Longitude",
        y = "Latitude"
      ) +
      
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Show the plot
    p
    
    # Save the plot
    plot_path <- paste0(
      "WP2_BrassiNiche/results_intermediate/species_occurrence_records_maps/",
      gsub(" ", "_", species),
      "_map.pdf"
    )
    dir.create(dirname(plot_path),
               recursive = TRUE,
               showWarnings = FALSE)
    ggsave(plot_path, p, width = 10, height = 6)
  })
}



## STEP 8: COLLECT AL GBIF DOWNLOAD RECORDS TO REQUEST DOI FOR PUBLICATION ----

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tools)
})

setwd("~/Google Drive/My Drive/Publications/2026_Hendriks_et_al_BrassiWood/BrassiWood")

in_dir <- "WP2_BrassiNiche/results_intermediate/species_occurrence_records_tables"

out_all_csv     <- "WP2_BrassiNiche/results_intermediate/all_gbifIDs_by_file.csv"
out_unique_txt  <- "WP2_BrassiNiche/results_intermediate/all_unique_gbifIDs.txt"
out_summary_csv <- "WP2_BrassiNiche/results_intermediate/all_gbifIDs_summary.csv"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
n_files <- length(files)

message("Found ", n_files, " csv files.")

# Counters
n_files_nonempty <- 0L
n_files_with_gbifID_column <- 0L
n_files_with_rows <- 0L
n_files_with_ids <- 0L

res_list <- vector("list", n_files)

for (i in seq_along(files)) {
  f <- files[i]
  file_name <- basename(f)
  species_name <- file_path_sans_ext(file_name)
  
  # Progress update every 50 files
  if (i %% 50 == 0 || i == 1 || i == n_files) {
    message(
      "[", i, "/", n_files, "] ",
      file_name
    )
  }
  
  # Skip truly empty files
  finfo <- file.info(f)
  if (is.na(finfo$size) || finfo$size == 0) {
    next
  }
  n_files_nonempty <- n_files_nonempty + 1L
  
  # Read only header
  header_names <- tryCatch(
    suppressMessages(
      names(read_csv(
        f,
        n_max = 0,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "minimal"
      ))
    ),
    error = function(e) NULL
  )
  
  # Skip malformed files or files without gbifID
  if (is.null(header_names) || !("gbifID" %in% header_names)) {
    next
  }
  n_files_with_gbifID_column <- n_files_with_gbifID_column + 1L
  
  # Read only gbifID column
  dat <- tryCatch(
    suppressMessages(
      read_csv(
        f,
        col_select = "gbifID",
        show_col_types = FALSE,
        progress = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  n_files_with_rows <- n_files_with_rows + 1L
  
  out <- data.frame(
    file = file_name,
    species = species_name,
    gbifID = as.character(dat$gbifID),
    stringsAsFactors = FALSE
  ) %>%
    mutate(gbifID = trimws(gbifID)) %>%
    filter(!is.na(gbifID), gbifID != "")
  
  if (nrow(out) == 0) {
    next
  }
  n_files_with_ids <- n_files_with_ids + 1L
  
  res_list[[i]] <- out
}

# Combine results
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
res <- bind_rows(res_list)

# Final stats
n_total_ids <- nrow(res)
n_unique_ids <- dplyr::n_distinct(res$gbifID)

summary_tbl <- data.frame(
  metric = c(
    "n_files_scanned",
    "n_files_nonempty",
    "n_files_with_gbifID_column",
    "n_files_with_rows_in_gbifID_column",
    "n_files_with_at_least_one_gbifID",
    "n_total_gbifID_entries",
    "n_unique_gbifIDs"
  ),
  value = c(
    n_files,
    n_files_nonempty,
    n_files_with_gbifID_column,
    n_files_with_rows,
    n_files_with_ids,
    n_total_ids,
    n_unique_ids
  )
)

write_csv(res, out_all_csv)
writeLines(sort(unique(res$gbifID)), out_unique_txt)
write_csv(summary_tbl, out_summary_csv)

message("Done.")
message("Files scanned: ", n_files)
message("Non-empty files: ", n_files_nonempty)
message("Files with gbifID column: ", n_files_with_gbifID_column)
message("Files with rows in gbifID column: ", n_files_with_rows)
message("Files with at least one gbifID: ", n_files_with_ids)
message("Total gbifID entries: ", n_total_ids)
message("Unique gbifIDs: ", n_unique_ids)
