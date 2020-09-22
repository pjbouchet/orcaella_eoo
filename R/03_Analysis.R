#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Preliminary assessment of the regional conservation status of snubfin
#'         dolphins (Orcaella heinsohni) in Western Australia
#'         
#'         Bouchet et al. (2020). Frontiers in Marine Science
#' 
#'                       --- 03: ANALYSIS ---
#'                  
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------

# Extract covariates  ---------------------------------------------------------

# Dolphin sightings

snub <- snub %>%
  dplyr::mutate(
    depth = raster::extract(x = depth, snub[, c("longitude", "latitude")]),
    dfresh = raster::extract(x = dfresh, snub[, c("longitude", "latitude")]) / 1000
  )

# Effort points

effort.pts <- read.csv("data/effort_pts.csv")

# Dolphins/boats can only occur in a minimum depth of 1 m

snub$depth <- ifelse(snub$depth>=-1, -1, snub$depth)
effort.pts$depth <- ifelse(effort.pts$depth>=-1, -1, effort.pts$depth)

# Minor adjustments

snub[snub$sighting_ID%in%c("obs_3085","obs_2927", "obs_3027", "obs_0361", "obs_0355"),]$depth <- -1
snub[snub$sighting_ID=="obs_0772",]$depth <- -15.7
snub[snub$sighting_ID=="obs_2829",]$depth <- -9.33
snub[snub$sighting_ID=="obs_2829",]$dfresh <- 55.113

#'---------------------------------------------
# Code to compute distance to freshwater inputs
#'---------------------------------------------

# This code is used to generate the dfresh raster imported earlier
# It requires the input raster to be up-scaled (x 35) for faster execution.

# Convert depth raster to a uniform surface of value 1

depth.surface <- depth
depth.surface[!is.na(depth.surface)] <- 1

# Upscale the raster by a factor of 35

depth.35 <- raster::aggregate(x = depth.surface, fact = 35)

# Create list objects to store transition matrices and related rasters
# This helps considerably speed up the code

list.rasters <- list()
list.conductance <- list()
lat.extent <- list()

# Extract the coordinates of each raster cell and sort by latitude (for efficiency)

env.df <- raster::as.data.frame(depth.35, xy = TRUE, na.rm = TRUE)
names(env.df) <- c('longitude', 'latitude', 'depth')
env.df <- env.df %>% dplyr::arrange(latitude)

# Create a safe version of the shortest_distance function, which will
# return NULL if any error arises

shortestdist_safe <- purrr::safely(shortest_distance)

# Run calculations for each raster cell

future::plan(multiprocess) # Parallel processing

dist.fresh <- furrr::future_map(
  .x = 1:nrow(env.df),
  .f = ~ {
    dat <- env.df[.x, ] # Extract each data.point

    if (dat$latitude <= -18.2) lcd <- FALSE else lcd <- TRUE

    # Convert to spatial object

    locs <- sp::SpatialPointsDataFrame(
      coords = cbind(dat$longitude, dat$latitude),
      data = dat, proj4string = CRSll
    )

    # Compute geodesic distance

    shortestdist_safe(
      least.cost = lcd,
      in.pts = locs,
      out.pts = rivers,
      closest.pts = 3,
      clip.width = 20000,
      cost.surface = depth.35
    )
  },
  .progress = TRUE
)

# Check whether any errors occurred and if so, identify where

lc.errors <- purrr::map(.x = dist.fresh, .f = "result") %>%
  purrr::map_dbl(., ~ as.numeric(is.null(.)))

errors.ind <- which(lc.errors == 1)

# If errors are found, replace values with straight line distances

straight.d <- furrr::future_map(.x = errors.ind, 
                                .f = ~{
                                  
                                  dat <- env.df[.x, ]
                                  
                                  if(dat$latitude <= -18.2) lcd <- FALSE else lcd <- TRUE
                                  
                                  locs <- sp::SpatialPointsDataFrame(coords = cbind(dat$longitude, dat$latitude),
                                                                     data = dat, proj4string = CRSll)
                                  
                                  shortestdist_safe(least.cost = FALSE,
                                                    in.pts = locs, 
                                                    out.pts = rivers, 
                                                    closest.pts = 3, 
                                                    clip.width = 20000, 
                                                    cost.surface = depth.35)},
                                .progress = TRUE)

straight.d <- purrr::map(straight.d, 'result') %>% do.call(c, .)
dist.fresh.corrected <- dist.fresh
for(i in 1:length(errors.ind)) dist.fresh.corrected[[errors.ind[i]]]$result <- straight.d[i]
distf.r <- purrr::map(dist.fresh.corrected, 'result') %>% do.call(c, .)

# One area near the top end of WA returned inconsistent results
# See the code in 05_Corrections for a fix.


# Bivariate kernel --------------------------------------------------------

# Extract primary sightings

snub.primary <- snub %>% split(x = ., f = .$sighting_class)
snub.primary <- snub.primary$primary

#'---------------------------------------------
# Calculate plug-in bandwidth
#'---------------------------------------------

h.biv <- purrr::map(
  .x = list(snub.primary, effort.pts),
  .f = ~ ks::Hpi(x = cbind(.x$dfresh, .x$depth))) %>%
  purrr::set_names(., nm = c("sightings", "effort")) %>%
  purrr::map(.x = ., .f = ~ {
    tmp <- .x
    colnames(tmp) <- c("dfresh", "depth")
    tmp
  })

#' ---------------------------------------------
# Generate 2D kernels (adjusted bandwidth values)
#' ---------------------------------------------

bivariate.kernel <- purrr::imap(
  .x = list(snub.primary, effort.pts),
  .f = ~ KernSmooth::bkde2D(
    x = cbind(.x$dfresh, .x$depth),
    bandwidth = list(c(1, 1.5), c(15, 25))[[.y]],
    range.x = list(c(0, 135), c(-110, -1)),
    gridsize = c(80, 95)
  )
)

# Convert to raster and crop

bivariate.r <- purrr::map(.x = bivariate.kernel, .f = ~ raster::raster(list(x = .x$x1, y = .x$x2, z = .x$fhat)))
bivariate.r <- purrr::map(
  .x = bivariate.r,
  .f = ~ raster::crop(x = .x, raster::extent(c(extent(.x)[1], 80, -95, extent(.x)[4])))
)

# Standardise by effort

bivariate.r <- bivariate.r[[1]] / bivariate.r[[2]]
bivariate.r <- bivariate.r / sum(bivariate.r[])

# Rescale to 0-1 range

bivariate.r <- rescale_raster(bivariate.r)
plot(bivariate.r, col = pals::parula(100))

#'---------------------------------------------
# Back-transform to geographic space
#'---------------------------------------------

# Load upscaled rasters

depth35 <- raster::raster("gis/kimb_depth_x35.tif")
dfresh35 <- raster::raster("gis/kimb_dfresh_x35.tif")

# Rasters need to have the same extents

dfresh35 <- raster::resample(x = dfresh35, y = depth35)

# Stack rasters and convert to df

kimb.grid <- raster::stack(depth35, dfresh35) %>% 
  raster::as.data.frame(., xy = TRUE, na.rm = TRUE) %>% 
  tibble::as_tibble(.) %>% 
  dplyr::rename(depth = kimb_depth_x35, dfresh = kimb_dfresh_x35) %>% 
  dplyr::mutate(dfresh = dfresh/1000)

# Compute percent volume contours (PVC)

k90 <- pvc(k = 0.9, convert.to.poly = TRUE, smoothing = 25)
k50 <- pvc(k = 0.5, convert.to.poly = TRUE, smoothing = 25)
k25 <- pvc(k = 0.25, convert.to.poly = TRUE, smoothing = 25)

# What percentage of the IUCN range polygon does this represent?

iucn.range.wa <- raster::crop(x = iucn.range, y = raster::extent(k90$poly))
rgeos::gArea(spgeom = sp::spTransform(x = k90$poly, CRSobj = CRSKim))/rgeos::gArea(spgeom = sp::spTransform(x = iucn.range.wa, CRSobj = CRSKim))


#'---------------------------------------------
# Assign inclusion probabilities
#'---------------------------------------------

snub.dat <- snub

snub.dat$incl.prob <- raster::extract(bivariate.r, snub.dat[, c("dfresh", "depth")])
snub.dat$incl.prob <- ifelse(snub.dat$sighting_class == "secondary", snub.dat$incl.prob, 1)

boxplot(incl.prob~dataset_ID, data = droplevels(snub.dat[snub.dat$sighting_class=="secondary",]), ylim = c(0,1))

# Bootstrap --------------------------------------------------------

# Number of iterations

n.iter <- 1000 

#'-------------------------------------------------
# Generate bootstrap resamples 
#'-------------------------------------------------

snub.boot <- split(x = snub.dat, f = snub.dat$sighting_class)
snub.boot <- list(data = snub.boot)

# Primary sightings are resampled with replacement

snub.boot$boot$primary <- purrr::map(
  .x = 1:n.iter,
  .f = ~ dplyr::sample_n(
    tbl = snub.boot$data[[1]],
    size = nrow(snub.boot$data[[1]]),
    replace = TRUE
  )
)

# Secondary sightings are resampled without replacement

select.mat <- matrix(nrow = nrow(snub.boot$data[[2]]), ncol = n.iter)

for (i in 1:n.iter) {
  select.mat[, i] <- purrr::map_int(
    .x = snub.boot$data[[2]]$incl.prob,
    .f = ~ rbinom(n = 1, size = 1, prob = .x)
  )
}

snub.boot$boot$secondary$prob <- purrr::map(.x = 1:n.iter, .f = ~ snub.boot$data[[2]][select.mat[, .x], ])

#' -------------------------------------------------
# Recombine datasets
#' -------------------------------------------------

snub.boot$combined <- purrr::map(
  .x = 1:n.iter,
  .f = ~ rbind(
    snub.boot$boot$primary[[.x]],
    snub.boot$boot$secondary$prob[[.x]]
  )
)


# Extent of occurrence (EOO) --------------------------------------------------------

#' -------------------------------------------------
# EOO as a minimum convex polygon (MCP)
#' -------------------------------------------------

plan(multiprocess)
eoo.mcp <- calc.eoo(dataset.list = snub.boot$combined, convex.hull = TRUE)

bci(eoo.mcp)

median(eoo.mcp)/(rgeos::gArea(spgeom = sp::spTransform(x = iucn.range, CRSobj = CRSKim))/1000000)

#'-------------------------------------------------
# EOO as an alpha-hull
#'-------------------------------------------------

# Define range of alpha values to test 

alpha.values <- seq(0.1, 2, by = 0.1)

# Find lowest alpha corresponding to an alpha-hull with no hollow spaces

pb <- progress_estimated(1000)

alpha.param <- purrr::map_dbl(
  .x = snub.boot$combined,
  .f = ~ smallest_alpha(alphaval = alpha.values, df = .),
  .progress = TRUE
)

# Calculate alpha-EOO

plan(multiprocess)
eoo.alpha <- calc.eoo(dataset.list = snub.boot$combined, alphaval = alpha.param, convex.hull = FALSE)

bci(eoo.alpha)


# Area of occupancy (AOO) -------------------------------------------------

pb <- progress_estimated(1000)

aoo <- purrr:::map_dbl(
  .x = snub.boot$combined,
  .f = ~ calc.aoo(
    input.data = .,
    coordinate.system = CRSKim,
    Cell_size_AOO = 2,
    nbe.rep.rast.AOO = 50
  )
)

# Mean and 95% confidence interval

bci(aoo)

# Conservation status -------------------------------------------------

# Using MCP-EOO

iucn.mcp <- eoo.mcp %>% 
  purrr:::map2_chr(.x = ., .y = aoo, .f = ~classify.threat(EOO = .x, AOO = .y))

# Using alpha-hull EOO

iucn.alpha <- eoo.alpha %>% 
  purrr:::map2_chr(.x = ., .y = aoo, .f = ~classify.threat(EOO = .x, AOO = .y))

table(iucn.mcp) %>% barplot(.)
table(iucn.alpha) %>% barplot(.)
