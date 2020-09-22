#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Preliminary assessment of the regional conservation status of snubfin
#'         dolphins (Orcaella heinsohni) in Western Australia
#'         
#'         Bouchet et al. (2020). Frontiers in Marine Science
#' 
#'                        --- 02: DATA ---
#'                  
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

# Note: sf requires GDAL, which can be tricky to install on some OS
# For Mac, see: https://stackoverflow.com/questions/44973639/trouble-installing-sf-due-to-gdal

require(remotes)
remotes::install_github("dkahle/ggmap", ref = "tidyup")
remotes::install_github("jjvanderwal/SDMTools")

pacman::p_load("tidyverse",             # Tidyverse
               "magrittr",              # Ceci n'est pas un pipe
               "lubridate",             # Date handling
               "raster",                # Raster and GIS operations
               "sp",                    # Methods for spatial data
               "rgeos",                 # Topology operations on geometries
               "sf",                    # Simple features
               "spThin",                # Spatial thinning of points
               "bbplot",                # ggplot2 charts in the style used by the BBC News data team
               "ConR",                  # Preliminary Assessment of Conservation Status
               "rgdal",                 # Working with geospatial data
               "spatialEco",            # Spatial Analysis and Modelling Utilities
               "gdistance",             # Geographic distance calculations
               "patchwork",             # The Composer of Plots
               "furrr",                 # Apply Mapping Functions in Parallel using Futures 
               "alphahull",             # Generalization of the Convex Hull
               "ggmap",                 # Plotting maps using ggplot 
               "ggspatial",             # Spatial Data Framework for ggplot2
               "smoothr",               # Smooth and tidy spatial features
               "KernSmooth",            # Kernel Smoothing
               "ks",                    # Kernel density estimation
               "SDMTools",              # Tools for processing data from SDMs
               "scales",                # Scale Functions for Visualization
               "pals")                  # Colour palettes

# Root directory

root <- getwd()

#'------------------------------------
# Set tibble options
#'------------------------------------

options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

set.seed(45) # Set the random seed

Sys.setenv(TZ = "Australia/West") # Set time zone

source(file.path(root, "R/01_Functions.R")) # Load functions

# GIS ----------------------------------------------------------------

#'------------------------------------
# Define coordinate systems
#'------------------------------------

CRSKim <- sp::CRS("+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") # Lambert conformal conic

CRSll <- sp::CRS("+proj=longlat +datum=WGS84")

#'------------------------------------
# IUCN range polygon
#'------------------------------------

iucn.range <- raster::shapefile("gis/iucn_range.shp")

#'------------------------------------
# Kimberley region
#'------------------------------------

kimb <- raster::shapefile(x = file.path(root, "gis", "kimb_coast_highres.shp")) %>%
  sp::spTransform(x = ., CRSobj = CRSll)

kimb_ocean <- raster::shapefile(x = file.path(root, "gis", "kimb_ocean.shp")) %>%
  sp::spTransform(x = ., CRSobj = CRSll)

kim_ocean_sf <- sf::st_as_sf(x = sp::spTransform(x = kimb_ocean, CRSobj = CRSKim))

#'-------------------------------------------------
# Bathymetry and distance to river rasters
#'-------------------------------------------------

depth <- raster::raster(file.path(root, "gis/kimb_depth_30m.tif"))
dfresh <- raster::raster(file.path(root, "gis/kimb_dfresh_30m.tif"))

#'-------------------------------------------------
# River mouths
#'-------------------------------------------------

rivers <- raster::shapefile(file.path("gis/river_mouths.shp"))
rivers <- sp::spTransform(x = rivers, CRSobj = CRSll)

# Survey tracks ----------------------------------------------------------------

tracks <- purrr::map(
  .x = list.files(path = file.path(root, "data/tracks"), pattern = ".shp$", full.names = TRUE),
  .f = ~ raster::shapefile(.x) %>% sp::spTransform(x = ., CRSobj = CRSll)
) %>%
  purrr::set_names(x = ., nm = gsub(pattern = "tracks_", replacement = "", x = gsub(pattern = ".shp", replacement = "", x = list.files(path = file.path(root, "data/tracks"), pattern = ".shp$"))))

all.tracks <- purrr::map(.x = tracks, .f = ~ as(.x, "SpatialLines")) %>%
  do.call(rbind, .)

# Effort (length of surveyed tracklines)

effort.km <- purrr::map(.x = tracks[1:7], .f = ~ .x %>%
  sf::st_as_sf(.) %>%
  sf::st_length(.) %>%
  sum(.) / 1000) %>%
  reshape2::melt(.) %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::rename(track_length = value, dataset_ID = L1)

# Generate 1 point per 10 km of trackline surveyed

effort.km <- effort.km %>% dplyr::mutate(numOfPoints = round(track_length / 10))

effort.points <- purrr::map2(
  .x = tracks[1:7], .y = effort.km$numOfPoints,
  .f = ~ sp::spsample(sp::spTransform(.x, CRSobj = CRSKim), n = .y, type = "regular")
) %>%
  do.call(rbind, .) %>%
  sp::spTransform(., CRSobj = CRSll)

# Convert to spatial object

effort.points <- sp::SpatialPointsDataFrame(effort.points, data.frame(ID = 1:length(effort.points)))

# Snubfin sightings -------------------------------------------------------------

#'-----------------------------------
# Import, compile and standardise sightings data
#'-----------------------------------

# Every dataset is formatted differently, so the first step is to compile the data
# into a standardised format

l.paths <- list.files(path = file.path(root, "data/sightings"), recursive = TRUE, pattern = ".csv")
l.files <- basename(l.paths)

# Import the data

csv.dat <- purrr::map(
  .x = l.paths,
  .f = ~ {
    suppressWarnings(readr::read_csv(file = file.path(root, "data/sightings", .x), col_types = cols())) %>%
      janitor::clean_names(.)
  }
) %>%
  purrr::set_names(., nm = l.paths)

# List dataset names and IDs

dataset.names <- as.list(gsub(pattern = ".csv", replacement = "", x = l.files))

dataset.IDs <- purrr::map(
  .x = l.paths,
  .f = ~ paste0("OH_", substr(
    x = .x,
    start = regexpr(pattern = "sightings_", text = .x)[1] + 10,
    stop = regexpr(pattern = "sightings_", text = .x)[1] + 12
  ))
)

# Extract the relevant columns

lon.names <- c("longitude", "lon", "proj_lon", "lon_proj", "long", "point_x", "gda_long", "longitude_processed")
lat.names <- c("latitude", "lat", "proj_lat", "lat_proj", "point_y", "gda_lat", "latitude_processed")
date.names <- c(
  "date_time_wa", "datetime", "date_time", "date", "date_time_raw", "datetime_local", "survey_date",
  "secondary_sampling_day", "sighting_start", "event_date_parsed"
)
size.names <- c(
  "bestno", "size", "best", "total_no_individuals", "grpsize_best", "group_size", "number_of_dolphins",
  "estimated_number", "estimate_number", "no", "n_tot", "pop_count", "numseen", "number_sighted", "individual_count", "indiv_cnt"
)

csv.filtered <- purrr::map(
  .x = csv.dat,
  .f = ~ .x %>% dplyr::select_if(names(.) %in% c(lon.names, lat.names, date.names, size.names))
)

colcount <- purrr::map_dbl(csv.filtered, ncol)
names(colcount) <- NULL
extracols <- which(colcount > 4)
missingcols <- which(colcount < 4)

for (i in extracols) {
  if (i %in% c(5, 9, 11, 12, 30, 31)) csv.filtered[[i]] <- csv.filtered[[i]] %>% dplyr::select(-date)
  if (i == 8) csv.filtered[[i]] <- csv.filtered[[i]] %>% dplyr::select(-date, -date_time)
  if (i == 10) csv.filtered[[i]] <- csv.filtered[[i]] %>% dplyr::select(-lat, -lon)
  if (i %in% c(32:37)) csv.filtered[[i]] <- csv.filtered[[i]] %>% dplyr::select(-longitude, -latitude)
  if (i == 32) csv.filtered[[i]] <- csv.filtered[[i]] %>% dplyr::select(-size)
}

# Rename and re-order columns

for (j in 1:length(csv.filtered)) {
  for (k in names(csv.filtered[[j]])) {
    if (k %in% lon.names) names(csv.filtered[[j]])[which(names(csv.filtered[[j]]) == k)] <- "longitude"
    if (k %in% lat.names) names(csv.filtered[[j]])[which(names(csv.filtered[[j]]) == k)] <- "latitude"
    if (k %in% date.names) names(csv.filtered[[j]])[which(names(csv.filtered[[j]]) == k)] <- "date_time"
    if (k %in% size.names) names(csv.filtered[[j]])[which(names(csv.filtered[[j]]) == k)] <- "group_size"
  }
  csv.filtered[[j]] <- csv.filtered[[j]] %>% dplyr::select("longitude", "latitude", "date_time", "group_size")
}

# Add dataset IDs and combine all data

snub <- purrr::map2(.x = csv.filtered, .y = dataset.IDs, .f = ~ .x %>% dplyr::mutate(dataset_ID = .y)) %>%
  purrr::map2(.x = ., .y = dataset.names, .f = ~ .x %>% dplyr::mutate(dataset_name = .y)) %>%
  do.call(rbind, .)

# Add sighting classes

class.df <- tibble::tibble(dataset_ID = unique(unlist(dataset.IDs)), sighting_class = c(rep("primary", 12), rep("secondary", 7)))

snub <- snub %>%
  dplyr::left_join(x = ., y = class.df, by = "dataset_ID") %>%
  dplyr::mutate(
    dataset_ID = factor(dataset_ID),
    dataset_name = factor(dataset_name),
    sighting_class = factor(sighting_class)
  )

# Print summary

summary(snub)


#'-----------------------------------
# Extract date & year
#'-----------------------------------

snub <- snub %>% dplyr::mutate(date = as.Date(date_time, tz = "Australia/West"),
                               year = lubridate::year(date_time))

#'---------------------------------------------
# Filter the data geographically and temporally
#'---------------------------------------------

snub <- snub %>% 
  dplyr::filter(!is.na(longitude) | !is.na(latitude)) %>% # Missing lat/lon
  dplyr::filter(!is.na(date)) %>% # Missing date
  dplyr::filter(longitude <= 129, longitude >=121) %>% 
  dplyr::filter(latitude >= -18.25, latitude <= -13) %>% 
  dplyr::filter(year >= 2004) %>% 
  dplyr::mutate(sighting_ID = paste0("obs_", addlzero(string = 1:nrow(.), n = 4)))

#'---------------------------------------------
# Remove duplicates
#'---------------------------------------------

# Based on lon, lat, date/time, and group size

snub.by.dataset <- split(snub, snub$dataset_ID) %>%
  purrr::map(.x = ., .f = ~ dplyr::distinct(.data = .x, latitude, longitude, date_time, group_size, .keep_all = TRUE)) %>%
  do.call(rbind, .)

snub <- dplyr::distinct(.data = snub.by.dataset, latitude, longitude, date_time, group_size, .keep_all = TRUE)
rm(snub.by.dataset)

# The code below is an alternative to using distinct() that counts and identifies
# duplicate records

# snub <- snub %>% dplyr::mutate(group_size = ifelse(is.na(group_size), -9999, group_size))
# 
# snub <- snub %>%
#   dplyr::group_by(latitude, longitude, date_time, group_size) %>%
#   dplyr::mutate(n_dup = n()) %>%
#   dplyr::ungroup()
# 
# snub$dup_obs <- NA
# 
# for(i in 1:nrow(snub)){
#   tmp <- snub %>% filter(longitude == snub$longitude[i],
#                          latitude == snub$latitude[i],
#                          date_time == snub$date_time[i],
#                          group_size == snub$group_size[i])
#   snub$dup_obs[i] <- paste0(tmp$sighting_ID, collapse = ", ")
# }
# 
# snub$dup_obs <- factor(snub$dup_obs)

#'---------------------------------------------
# Project points on land to the nearest location on water
#'---------------------------------------------

# Convert to spatial object first

snublocs <- sp::SpatialPointsDataFrame(coords = cbind(snub$longitude, snub$latitude),
                                       data = snub, 
                                       proj4string = CRSll)

# Points over land - some are really on land, and should be discarded
# Others are actually really on the water, but the available landmass
# shapefile isn't of high enough resolution to discriminate.

land.snub <- snub[is.na(sp::over(snublocs, kimb_ocean)[,1]),] %>%
  dplyr::filter(!is.na(longitude) & !is.na(latitude))

# Prince Regent River sightings must be kept

prr.land <- land.snub %>% dplyr::filter(longitude>125.1616 & longitude<125.3459 & latitude > -15.6463 & latitude < -15.4803) %>% dplyr::pull(sighting_ID)     
land.snub <- land.snub %>% dplyr::filter(!sighting_ID%in%prr.land)
  
# Generate and clip 1 km buffers around each point

land.locs <- sp::SpatialPointsDataFrame(
  coords = cbind(land.snub$longitude, land.snub$latitude),
  data = land.snub,
  proj4string = CRSll
) %>%
  sp::spTransform(x = ., CRSobj = CRSKim) %>%
  sf::st_as_sf(.)

land.buffers <- sf::st_buffer(x = land.locs, dist = 1000) %>%
  sf::st_intersection(x = ., y = kim_ocean_sf)

# Sample a random point within each buffer

land.newpos <- sf::st_sample(x = land.buffers, size = rep(1, nrow(land.buffers)), type = "random") %>%
  sf::st_transform(x = ., crs = CRSll) %>%
  sf::st_coordinates(.)

# Replace the coordinates accordingly

snub.onland <- land.buffers %>%
  tibble::as_tibble(.) %>%
  cbind(., land.newpos) %>%
  dplyr::select(., tidyselect::all_of(c(names(snub)[!names(snub) %in% c("longitude", "latitude")], c("X", "Y")))) %>%
  tibble::as_tibble(.) %>%
  dplyr::rename(longitude = X, latitude = Y) %>%
  dplyr::select(tidyselect::all_of(names(snub)))

# Update the main dataset

snub <- snub %>%
  dplyr::filter(!sighting_ID %in% c(land.snub$sighting_ID)) %>%
  rbind(., snub.onland)

#' ---------------------------------------------
# Thin data within 100 m buffers
#' ---------------------------------------------

# Duplicates are defined as secondary sightings made
# within a 100 m of each other on the same date

# First split the data by date and sighting class

snub <- snub %>% dplyr::mutate(species = "Orcaella heinsohni")
snub.by.date.class <- split(x = snub, f = snub$sighting_class) %>%
  purrr::map(.x = ., .f = ~ split(.x, .x$date))

# Perform spatial thinning

pb <- dplyr::progress_estimated(n = length(snub.by.date.class$secondary))

thinned.locs <- purrr::map(
  .x = snub.by.date.class$secondary,
  .f = ~ {
    pb$tick()$print()
    spThin::thin(
      loc.data = .x,
      long.col = "longitude",
      lat.col = "latitude",
      spec.col = "species",
      thin.par = 0.1,
      reps = 1,
      locs.thinned.list.return = TRUE,
      write.files = FALSE, verbose = FALSE
    )
  }
)

# Put the data back together

snub.thinned <- suppressMessages(purrr::imap(
  .x = thinned.locs,
  .f = ~ .x %>%
    as.data.frame(.) %>%
    dplyr::rename(longitude = Longitude, latitude = Latitude) %>%
    dplyr::left_join(., snub.by.date.class$secondary[[.y]]) %>%
    tibble::as_tibble(.)
)) %>% do.call(rbind, .)

snub <- snub %>%
  dplyr::filter(sighting_class == "primary") %>%
  rbind(., snub.thinned)
