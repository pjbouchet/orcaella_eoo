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

#'---------------------------------------------
# Distance to freshwater inputs
#'---------------------------------------------

# Convert depth raster to a uniform surface of value 1

depth.surface <- depth
depth.surface[!is.na(depth.surface)] <- 1

# Create list objects to store transition matrices and related rasters
# This helps considerably speed up the code

list.rasters <- list()
list.conductance <- list()

# Create a safe version of the shortest_distance function, which will
# return NULL if any error arises

shortestdist_safe <- purrr::safely(shortest_distance)

pb <- dplyr::progress_estimated(n = nrow(snub)) # Set up progress bar

dist.rivers <- purrr::map(.x = 1:nrow(snub), 
                          .f = ~{
                            
                            dat <- snub[.x, ]
                            
                            if(dat$latitude <= -18.2) lcd <- FALSE else lcd <- TRUE
                            
                            locs <- sp::SpatialPointsDataFrame(coords = cbind(dat$longitude, dat$latitude),
                                                               data = dat, proj4string = CRSll)
                            
                            shortestdist_safe(least.cost = lcd,
                                              in.pts = locs, 
                                              out.pts = rivers, 
                                              closest.pts = 3, 
                                              clip.width = 20000, 
                                              cost.surface = depth.surface, 
                                              show.plot = FALSE, 
                                              verbose = FALSE)},
                          .pb = pb)

# Check if any errors occurred

# lc.errors <- purrr::map(.x = dist.rivers, .f = "error") %>%
#   purrr::map_dbl(., ~as.numeric(is.null(.)))
# errors.ind <- which(lc.errors==0)

snub$dfresh <- purrr::map_dbl(.x = dist.rivers, .f = "result") %>% do.call(c, .)/1000

#'---------------------------------------------
# Bathymetric depth
#'---------------------------------------------

snub <- snub %>% 
  dplyr::mutate(depth = raster::extract(x = depth, snub[, c("longitude", "latitude")]))

# Dolphins can only occur in a minimum depth of 1 m

snub$depth <- ifelse(snub$depth>=-1, -1, snub$depth)

#'---------------------------------------------
# Effort
#'---------------------------------------------

# Can do the same for effort points

effort.pts <- readr::read_csv("data/effort_pts.csv", col_types = cols())
effort.pts$depth <- ifelse(effort.pts$depth>=-1, -1, effort.pts$depth)

rm(depth.surface, list.rasters, list.conductance)

# Bivariate kernel --------------------------------------------------------

# Extract primary sightings

snub.primary <- snub %>% split(x = ., f = .$sighting_class)
snub.primary <- snub.primary$primary

# Calculate plug-in bandwidth

h.biv <- purrr::map(.x = list(snub.primary, effort.pts),
                   .f = ~ks::Hpi(x = cbind(.x$dfresh, .x$depth))) %>% 
  purrr::set_names(., nm = c("sightings", "effort")) %>% 
  purrr::map(.x = ., .f = ~{tmp <- .x
  colnames(tmp) <- c("dfresh", "depth")
  tmp})

# Generate 2D kernels (adjusted bandwidth values)

bivariate.kernel <- purrr::imap(.x = list(snub.primary, effort.pts),
                       .f = ~KernSmooth::bkde2D(x = cbind(.x$dfresh, .x$depth),
                                                bandwidth = list(c(1, 1.5), c(15, 25))[[.y]],
                                                range.x = list(c(0, 135), c(-110, -1)),
                                                gridsize = c(100, 100)))


# Convert to raster and crop

bivariate.r <- purrr::map(.x = bivariate.kernel, .f = ~raster::raster(list(x = .x$x1, y = .x$x2,z = .x$fhat)))
bivariate.r <- purrr::map(.x = bivariate.r, 
                          .f = ~raster::crop(x = .x, raster::extent(c(extent(.x)[1], 80, -95,extent(.x)[4]))))

# Standardise by effort

bivariate.r <- bivariate.r[[1]]/bivariate.r[[2]]

# Rescale to 0-1 range

bivariate.r <- rescale_raster(bivariate.r)

# Compute contours / isopleths

bivariate.contours <- raster::rasterToContour(x = bivariate.r, levels = c(0.5, 0.9))
bivariate.contours <- as(bivariate.contours, "SpatialLines")

bivariate.contours <- sp::SpatialPolygons(
  lapply(1:length(bivariate.contours), 
         function(i) Polygons(lapply(coordinates(bivariate.contours)[[i]], function(y) Polygon(y)), as.character(i))))


# Assign inclusion probabilities

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

snub.boot$boot$primary <- purrr::map(.x = 1:n.iter, 
                          .f = ~dplyr::sample_n(tbl = snub.boot$data[[1]], 
                                                size = nrow(snub.boot$data[[1]]), 
                                                replace = TRUE)) 

# # Only used for comparative purposes
# 
# snub.boot$boot$secondary$noprob <- purrr::map(.x = 1:n.iter, 
#                                               .f = ~dplyr::sample_n(tbl = snub.boot$data[[2]], 
#                                                     size = nrow(snub.boot$data[[2]]), 
#                                                     replace = TRUE)) 

# Secondary sightings are resampled without replacement

select.mat <- matrix(nrow = nrow(snub.boot$data[[2]]), ncol = n.iter)

for(i in 1:n.iter){
select.mat[,i] <- purrr::map_int(.x = snub.boot$data[[2]]$incl.prob, 
                                 .f = ~rbinom(n = 1, size = 1, prob = .x))}

snub.boot$boot$secondary$prob <- purrr::map(.x = 1:n.iter, .f = ~snub.boot$data[[2]][select.mat[,.x],])

#'-------------------------------------------------
# Recombine datasets
#'-------------------------------------------------

snub.boot$combined <- purrr::map(.x = 1:n.iter, 
                                 .f = ~rbind(snub.boot$boot$primary[[.x]], 
                                             snub.boot$boot$secondary$prob[[.x]]))


# Extent of occurrence (EOO) --------------------------------------------------------

#'-------------------------------------------------
# EOO as a minimum convex polygon (MCP)
#'-------------------------------------------------

eoo.mcp <- calc.eoo(dataset.list = snub.boot$combined, convex.hull = TRUE)     
bci(eoo.mcp)

#'-------------------------------------------------
# EOO as an alpha-hull
#'-------------------------------------------------

# Define range of alpha values to test 

alpha.values <- seq(0.1, 5, by = 0.1)

# Find lowest alpha corresponding to an alpha-hull with no hollow spaces

pb <- progress_estimated(n.iter, 0) # Set up progress bar
alpha.param <- purrr::map_dbl(.x = snub.boot$combined, .f = ~smallest_alpha(alphaval = alpha.values, df = .))
                                      
# Calculate alpha-EOO

eoo.alpha <- calc.eoo(dataset.list = snub.boot$combined, alphaval = alpha.param, convex.hull = FALSE)                
bci(eoo.alpha)


# Area of occupancy (AOO) -------------------------------------------------

aoo <- purrr:::map_dbl(.x = snub.boot$combined, 
                       .f = ~calc.aoo(input.data = ., 
                                      coordinate.system = CRSKim,
                                      Cell_size_AOO = 2, 
                                      nbe.rep.rast.AOO = 50))


# Mean and 95% confidence interval

bci(aoo$all)


# Conservation status -------------------------------------------------

# Using MCP-EOO

iucn.mcp <- eoo.mcp %>% 
  purrr:::map2_chr(.x = ., .y = aoo, .f = ~classify.threat(EOO = .x, AOO = .y))

iucn.alpha <- eoo.alpha %>% 
  purrr:::map2_chr(.x = ., .y = aoo, .f = ~classify.threat(EOO = .x, AOO = .y))

table(iucn.alpha) %>% barplot(.)
