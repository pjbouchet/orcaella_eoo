#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#' 
#'  Preliminary assessment of the regional conservation status of snubfin
#'         dolphins (Orcaella heinsohni) in Western Australia
#'         
#'         Bouchet et al. (2020). Frontiers in Marine Science
#' 
#'                  --- 01: FUNCTIONS ---
#'                  
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------


# Data manipulation ----------------------------------------------------------

#'---------------------------------------------
# Adaptation of the movecost function from package {movecost}
# Calculates least cost path from input to target points/features
#'---------------------------------------------

shortest_distance <- function(least.cost = TRUE,
                              in.pts, 
                              out.pts, 
                              closest.pts = 3, 
                              clip.width = 2500,
                              cost.surface,
                              proj.coord.system = CRSKim){
  
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param least.cost logical. If TRUE, returns least-cost distances. If FALSE, returns straight-line distances.
  #' @param in.pts Input points, from which distances are calculated.
  #' @param out.pts Target points, to which distances are calculated.
  #' @param closest.pts Number of target points considered in the calculations.
  #' @param clip.width Width of the buffer used to define crop the cost surface.
  #' @param cost.surface Input cost surface.
  #' @param proj.coord.system Projected coordinate system.
  #'---------------------------------------------

  
  #'-------------------------------------------------
  # Project points
  #'-------------------------------------------------
  
  in.pts <- sp::spTransform(x = in.pts, CRSobj = proj.coord.system)
  out.pts <- sp::spTransform(x = out.pts, CRSobj = proj.coord.system)
  
  #'-------------------------------------------------
  # Calculate Cartesian distances
  #'-------------------------------------------------
  
  dd <- rgeos::gDistance(spgeom1 = in.pts, spgeom2 = out.pts, byid = TRUE)
  
  #'-------------------------------------------------
  # Calculate least cost paths (shortest distances), if required 
  #'-------------------------------------------------
  
  if(!least.cost){
    
    min(dd)  
    
  }else{
    
    # Because the calculation of least cost paths is computer-intensive (and therefore slow)
    # and most sightings are clustered in space, such that the nearest river will be the same
    # for multiple points, we can be more efficient by saving the transition matrices generated
    # in one run and re-use them in subsequent runs, where appropriate.
    
    # Save the X closest river mouths
    
    target.pts <- sp::SpatialPoints(coords = out.pts@coords[which(dd%in%sort(dd)[1:closest.pts]),],
                                    proj4string = proj.coord.system)
    
    # Combine all points
    
    all.pts <- rgeos::gUnion(spgeom1 = in.pts, target.pts)
    LL <- sp::spTransform(all.pts, CRSll)
    
    # If rasters have been previously been saved, check whether the points overlap any of them 
    
    raster.index <- NULL
    
    if(!purrr::is_empty(lat.extent)){
      
      x1 <- purrr::map_int(.x = lat.extent, .f = ~{
        raster::extent(.x)[1] <= raster::extent(LL)[1]
      })
      x2 <- purrr::map_int(.x = lat.extent, .f = ~{
        raster::extent(.x)[2] >= raster::extent(LL)[2]
      })
      y1 <- purrr::map_int(.x = lat.extent, .f = ~{
        raster::extent(.x)[3] <= raster::extent(LL)[3]
      })
      y2 <- purrr::map_int(.x = lat.extent, .f = ~{
        raster::extent(.x)[4] >= raster::extent(LL)[4]
      })
      
      ex.tbl <- cbind(x1, x2, y1, y2)
      
      totals <- apply(X = ex.tbl, MARGIN = 1, FUN = sum)
      
      
      if(any(totals==4)){
        
        raster.index <- max(which(totals==4))
        # is.r <<- TRUE
        Conductance <- list.conductance[[raster.index]]}
      
    } 
   
    # If suitable rasters are identified, use the corresponding transition matrix
    # instead of recalculating it
    
    if(is.null(raster.index)){
    
      #'-------------------------------------------------
      # Create buffer
      #'-------------------------------------------------
      
      r.ext <- as(extent(all.pts), "SpatialPolygons")
      r.ext <- rgeos::gBuffer(spgeom = r.ext, width = clip.width)
      sp::proj4string(r.ext) <- proj.coord.system
      
      #'-------------------------------------------------
      # Clip raster
      #'-------------------------------------------------
      
      r.ext <- sp::spTransform(x = r.ext, CRSobj = CRSll)
      r <- raster::crop(x = cost.surface, y = r.ext)
      r <- raster::mask(x = r, mask = r.ext)
      
      list.rasters <<- append(x = list.rasters, values = r) # Add to raster list
      lat.extent <<- append(x = lat.extent, values = raster::extent(r)) # Add to raster list
      
      #'-------------------------------------------------
      # Compute transition matrix, with geographic correction
      #'-------------------------------------------------
      
      f <- function(x) 1
      rr <- gdistance::transition(r, f, directions = 16)
      Conductance <- gdistance::geoCorrection(rr)
      
      list.conductance <<- append(x = list.conductance, values = Conductance) # Add to 
      
    }
    
    
    #'-------------------------------------------------
    # Re-project points
    #'-------------------------------------------------
    
    in.pts <- sp::spTransform(x = in.pts, CRSobj = CRSll)
    target.pts <- sp::spTransform(x = target.pts, CRSobj = CRSll)
    
    #'-------------------------------------------------
    # Identify shortest path
    #'-------------------------------------------------
    
    safe_path <- purrr::safely(gdistance::shortestPath)
    
    sPath <- purrr::map(1:3, ~safe_path(x = Conductance, 
                                     origin = coordinates(in.pts), 
                                     goal = coordinates(target.pts)[.x,], 
                                     output = "SpatialLines"))
    
    sPath <- purrr::map(sPath, "result") %>% purrr::compact(.) %>% do.call(rbind, .)
    
    # sPath <- suppressWarnings(purrr::map(.x = 1:closest.pts,
    #                                      .f = ~gdistance::shortestPath(x = Conductance, 
    #                                                   origin = coordinates(in.pts), 
    #                                                   goal = coordinates(target.pts)[.x,], 
    #                                                   output = "SpatialLines")))

    # ind <- purrr::map_dbl(.x = sPath@lines, .f = ~.x@Lines[[1]]@coords %>% nrow())
    # sPath <- sPath[which(ind==min(ind[ind>1])),]
    
    #'-------------------------------------------------
    # Calculate length
    #'-------------------------------------------------
    
    sPath.min <- min(rgeos::gLength(sp::spTransform(sPath, CRSKim), byid = TRUE))
    return(sPath.min)
    # rgeos::gLength(sp::spTransform(sPath, CRSKim))
    
  } 
}

#'---------------------------------------------
# Function to compute percent volume contours
#'---------------------------------------------

# http://www.spatialecology.com/htools/pctvolcontour.php
# A percent volume contour represents the boundary of the area that 
# contains x% of the volume of a probability density distribution. 

pvc <- function(k, convert.to.poly = FALSE, smoothing = 20){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param k Desired % of the volume  
  #' @param poly Logical. Converts raster to polygons.
  #' @param smoothing Amount of smoothing for polygons
  #'---------------------------------------------
  
  # Can also be done by numeric optimisation as per the code below
  
  # extract_k <-  function(input.raster, kern = 0.9){
  #   
  #   # https://gis.stackexchange.com/questions/272950/95-contour-from-kernel-density-estimates
  #   
  #   cover <- function(r){
  #     Total <- sum(r[]) * kern
  #     function(t) sum(r[][r[]>t])-Total
  #   }
  #   
  #   # One-dimensional root finding
  #   val <- uniroot(cover(input.raster), c(0, max(input.raster[])))
  #   output.raster <- input.raster
  #   output.raster[output.raster<val$root] <- NA
  #   return(output.raster)
  # }
  # 
  # k90.r <- extract_k(input.raster = bivariate.r, kern = 0.9)
  # k50.r <- extract_k(input.raster = bivariate.r, kern = 0.5)
  
  # Calculate a percent volume on a raster or based on a systematic sample
  
  kr <- spatialEco::raster.vol(x = bivariate.r, p = k)
  
  k.pts <- raster::extract(x = kr, y = kimb.grid[, c("dfresh", "depth")])
  k.pts[is.na(k.pts)] <- 0
  k.pts <- kimb.grid[k.pts==1,]
  k.pts$val <- 1
  
  res <- raster::rasterFromXYZ(xyz = k.pts[, c("x", "y", "val")])
  sp::proj4string(res) <- CRSll
  
  if(convert.to.poly){
    
    HR.poly <- raster::rasterToPolygons(res, n = 16, na.rm = TRUE, digits = 4, dissolve = TRUE)
    
    s <- "polyID"
    HR.polys<-list(0)
    pr <- 0
    
    ##define the area of tiny fragment polygons (4 pixles or smaller)
    min.poly.size <- 4*((res@extent@xmax-res@extent@xmin)/res@ncols)^2
    
    for (p in 1:length(HR.poly@polygons[[1]]@Polygons)) {
      if (HR.poly@polygons[[1]]@Polygons[[p]]@area > min.poly.size) {
        pr<-pr+1
        HR.polys[[pr]]<-HR.poly@polygons[[1]]@Polygons[[p]]
      }
    }
    
    HR.polys <- sp::Polygons(HR.polys, ID = s) # convert to polygons object
    HR.polys <- maptools::checkPolygonsHoles(HR.polys) #assign holes correctly
    HR.polys <- sp::SpatialPolygons(list(HR.polys)) #convert to spatial polygons
    
    HR.polys <- smoothr::smooth(HR.polys, method = "ksmooth", smoothness = smoothing)
    sp::proj4string(HR.polys) <- CRSll
    
    res <- list(raster = res, poly = HR.polys)
  }
  
  return(res)
  
}

# Conservation assessment -------------------------------------------------

#'---------------------------------------------
# Functions to count the number of points outside a given spatial polygon
#'---------------------------------------------

count_points <- function(polyg, pts){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param polyg Spatial polygon within which points are counted.
  #' @param pts Input points.
  #'---------------------------------------------
  
  # Assign projection
  if(is.na(proj4string(polyg))) proj4string(polyg) <- CRSll
  
  # Turn polygon into data.frame object
  polyg <- broom::tidy(polyg)
  
  # Extract lat/lon of points
  if(ncol(pts)>2) pts <- pts %>% dplyr::select(., longitude, latitude)
  
  # Determine the number of points outside the alpha hull
  pts.out <- SDMTools::pnt.in.poly(pts, polyg[, c("long", "lat")])
  pts.out <- pts.out %>% 
    dplyr::filter(pip==0) %>% nrow()
  
  return(pts.out)
}

#'---------------------------------------------
# Functions to find the smallest alpha value
#'---------------------------------------------

smallest_alpha <- function(alphaval, df){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param alphaval Alpha value(s) to test.
  #' @param df Input data.frame, containing occurrence records (incl. coordinates)
  #'---------------------------------------------
  
  if(exists("pb")) pb$tick()$print()
  
  # Calculate Alpha hulls for each value of alpha
  
  alpha_test <- purrr:::map(.x = alphaval,
                            .f = ~dplyr::select(df, latitude, longitude) %>% 
                              ConR::EOO.computing(XY = .,
                                                  Name_Sp = "oh",
                                                  exclude.area = FALSE,
                                                  export_shp = TRUE,
                                                  country_map = kimb_ocean,
                                                  method.range = "alpha.hull",
                                                  alpha = .x,
                                                  buff.alpha = 0.001,
                                                  write_results = FALSE, 
                                                  show_progress = FALSE))
  alpha_test <- alpha_test %>% 
    purrr::map(., ~.x[grep("spatial", names(.x))]) %>% # Extract the polygon shapefiles
    purrr::flatten(.) %>% 
    purrr::map_dbl(.x = .,
                   .f = ~dplyr::select(df, longitude, latitude) %>% 
                     count_points(polyg = .x, pts = .)) # Count points outside polygon
  
  alpha_small <- alpha.values[[which(alpha_test==0)[1]]]
  
  return(alpha_small)
}


#'---------------------------------------------
# Function to calculate EOO for bootstrapped datasets
#'---------------------------------------------

calc.eoo <- function(dataset.list, alphaval = 1, convex.hull = FALSE){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param dataset.list List object containing the replicate datasets.
  #' @param alphaval Alpha value used to compute alpha-hulls (only when mcp = FALSE)
  #' @param convex.hull Logical. If TRUE, the EOO is calculated as a minimum convex polygon.
  #'---------------------------------------------
  
  furrr:::future_map2_dbl(.x = dataset.list, 
                   .y = alphaval,
                   .f = ~{
                     
                     # Extract the lat and lon coordinates
                     xy <- dplyr::select(.x, latitude, longitude)
                     
                     if(convex.hull) met <- "convex.hull" else met <-  "alpha.hull"
                     
                     # Calculate the EOO based on an alpha hull, saving the SHP
                     eoo <- ConR::EOO.computing(XY = xy,
                                                Name_Sp = "Orcaella heinsohni",
                                                exclude.area = TRUE,
                                                export_shp = TRUE,
                                                country_map = kimb_ocean,
                                                method.range = met,
                                                alpha = .y,
                                                buff.alpha = 0.001,
                                                write_results = FALSE, 
                                                show_progress = FALSE)
                     
                     eoo[[1]]}, .progress = TRUE)
  
}

#'---------------------------------------------
# Function to calculate the Area of Occupancy (AOO)
#'---------------------------------------------

# Adapted from ConR functions

calc.aoo <- function (input.data,
                      coordinate.system = NULL,
                      Cell_size_AOO = 2, 
                      nbe.rep.rast.AOO = NULL) 
{
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param input.data Input data.frame.
  #' @param coordinate.system Relevant projected coordinate system.
  #' @param Cell_size_AOO Grid size (in km) used for estimating AOO. Defaults to 2 as per IUCN guidelines.
  #' @param nbe.rep.rast.AOO Number of rasters with random starting position for estimating the AOO. NULL by default, but some minimal translation of the raster are still done.
  #'---------------------------------------------
  
  if(exists("pb")) pb$tick()$print()
  
  if(is.null(coordinate.system)) coordinate.system <- raster::crs("+proj=cea +lon_0=Central Meridian+lat_ts=Standard Parallel+x_0=False Easting+y_0=False Northing +ellps=WGS84")
  
  #'---------------------------------------------
  # Extract coordinates
  #'---------------------------------------------
  
  XY <- input.data[, c("longitude", "latitude")]
  
  #'---------------------------------------------
  # Project points
  #'---------------------------------------------
  
  coordEAC <- SpatialPoints(coords = XY, proj4string = CRSll)
  coordEAC <- sp::spTransform(x = coordEAC, CRSobj = coordinate.system)
  
  coordEAC <- coordinates(coordEAC)
  rownames(coordEAC) <- seq(1, nrow(coordEAC), 1)
  
  #'---------------------------------------------
  # Spatial extent
  #'---------------------------------------------
  
  Corners <- rbind(c(min(coordEAC[, 1]), max(coordEAC[, 1])), 
                   c(min(coordEAC[, 2]), max(coordEAC[, 2])))
  
  #'---------------------------------------------
  # Random translations of raster starting position
  #'---------------------------------------------
  
  if (is.null(nbe.rep.rast.AOO)) {
    
    Occupied_cells <- vector(mode = "numeric", length = nbe.rep.rast.AOO)
    decal <- c(0, 1, 2, 3)
    
    for (h in decal) {
      
      ext = extent(floor(Corners[1, 1]) - h * (Cell_size_AOO * 
                                                 1000/4) - 2 * Cell_size_AOO * 1000, floor(Corners[1, 2]) + h * (Cell_size_AOO * 1000/4) + 2 * Cell_size_AOO * 
                     1000, floor(Corners[2, 1]) - h * (Cell_size_AOO * 1000/4) - 2 * Cell_size_AOO * 1000, floor(Corners[2, 2]) + h * (Cell_size_AOO * 1000/4) + 2 * Cell_size_AOO * 1000)
      
      r = raster::raster(ext, resolution = Cell_size_AOO * 1000, 
                         crs = crs(coordinate.system))
      
      r2_AOO <- raster::rasterize(coordEAC, r)
      
      OCC <- length(which(!is.na(raster::getValues(r2_AOO))))
      Occupied_cells[h] <- OCC
      
      if (OCC == 1) 
        break
    }
    
  } # End if is.null
  
  if (!is.null(nbe.rep.rast.AOO)) {
    
    Occupied_cells <- vector(mode = "numeric", length = nbe.rep.rast.AOO)
    
    for (h in 1:nbe.rep.rast.AOO) {
      
      rd.1 <- runif(1) * Cell_size_AOO * 1000
      rd.2 <- runif(1) * Cell_size_AOO * 1000
      ext = extent(floor(Corners[1, 1]) - rd.1 - 2 * Cell_size_AOO * 
                     1000, floor(Corners[1, 2]) + rd.1 + 2 * Cell_size_AOO * 
                     1000, floor(Corners[2, 1]) - rd.2 - 2 * Cell_size_AOO * 
                     1000, floor(Corners[2, 2]) + rd.2 + 2 * Cell_size_AOO * 
                     1000)
      r = raster::raster(ext, resolution = Cell_size_AOO * 1000, 
                         crs = crs(coordinate.system))
      
      r2_AOO <- raster::rasterize(coordEAC, r)
      OCC <- length(which(!is.na(raster::getValues(r2_AOO))))
      Occupied_cells[h] <- OCC
      # Occupied_cells <- c(Occupied_cells, OCC)
      if (OCC == 1) 
        break
    }
  } # End if is.null
  
  #'---------------------------------------------
  # Calculate AOO
  #'---------------------------------------------
  
  Occupied_cells <- Occupied_cells[Occupied_cells > 0]
  Occupied_cells <- min(Occupied_cells) # Minimum value across random starts
  AOO <- Occupied_cells * Cell_size_AOO * Cell_size_AOO
  return(AOO)
}

#'---------------------------------------------
# Function to determine threat status based on
# EOO + AOO (according to IUCN criterion B)
#'---------------------------------------------

classify.threat <- function (EOO, AOO){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param EOO Extent of occurrence.
  #' @param AOO Area of occupancy.
  #'---------------------------------------------
  
  #'---------------------------------------------
  # Apply IUCN criteria
  #'---------------------------------------------
  
  if (EOO < AOO) EOO <- AOO
  
  #'---------------------------------------------
  # Criterion B1 - Extent of Occurrence
  #'---------------------------------------------
  
  if (EOO < 20000) {
    Rank_EOO <- 3
    if (EOO < 5000) {
      Rank_EOO <- 2
      if (EOO < 100) {
        Rank_EOO <- 1
      }
    }
  } else {
    (Rank_EOO <- 4)}
  
  #'---------------------------------------------
  # Criterion B2 - Area of Occupancy
  #'---------------------------------------------
  
  if (AOO < 2000) {
    Rank_AOO <- 3
    if (AOO < 500) {
      Rank_AOO <- 2
      if (AOO < 10) {
        Rank_AOO <- 1
      }
    }
  } else {
    Rank_AOO <- 4}
  
  #'---------------------------------------------
  # Number of locations
  #'---------------------------------------------
  
  # if (Nbe_locs <= 10) {
  #   Rank_Loc <- 3
  #   if (Nbe_locs <= 5) {
  #     Rank_Loc <- 2
  #     if (Nbe_locs == 1) {
  #       Rank_Loc <- 1
  #     }
  #   }
  # } else {
  #   Rank_Loc <- 4}
  
  #'---------------------------------------------
  # Work out threat level
  #'---------------------------------------------
  
  # Rank_B1a <- max(Rank_EOO, Rank_Loc)
  # Rank_B2a <- max(Rank_AOO, Rank_Loc)
  # Rank_CriteriaB <- min(Rank_B1a, Rank_B2a)
  
  Rank_CriteriaB <- min(Rank_EOO, Rank_AOO)
  
  if (Rank_CriteriaB == 1) Cat <- "CR"
  if (Rank_CriteriaB == 2) Cat <- "EN"
  if (Rank_CriteriaB == 3) Cat <- "VU"
  if (Rank_CriteriaB > 3) Cat <- "LC or NT"
  
  # if (Rank_CriteriaB == 3 && Nbe_locs > 0 && Nbe_locs < 11) Cat <- "VU"
  # if (Rank_CriteriaB > 3 && Nbe_locs >= 0) Cat <- "LC or NT"
  
  # if (Rank_B1a > Rank_B2a) Cat_Code <- paste(Cat, "B2a")
  # if (Rank_B1a < Rank_B2a) Cat_Code <- paste(Cat, "B1a")
  # if (Rank_B1a == Rank_B2a) Cat_Code <- paste(Cat, "B1a+B2a")
  
  if (Rank_EOO > Rank_AOO) Cat_Code <- paste(Cat, "(B2)")
  if (Rank_EOO < Rank_AOO) Cat_Code <- paste(Cat, "(B1)")
  if (Rank_EOO == Rank_AOO) Cat_Code <- paste(Cat, "(B1+B2)")
  
  return(Cat_Code)
  
}


#'---------------------------------------------
# Mapping functions
#'---------------------------------------------

# Plot dolphin sightings in ggmap using the correct base layer

map_sightings <- function(baselayer, 
                          col.primary = "#f1a340", 
                          col.secondary = "#5ab4ac",
                          pts.alpha = 0.75, 
                          pts.size = 3){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param baselayer Input baselayer object, saved from ggmap.
  #' @param col.primary Colour to use for primary sightings.
  #' @param col.secondary Colour to use for secondary sightings.
  #' @param pts.alpha Point transparency (0 = transparent to 1 = fully opaque).
  #' @param pts.size Point size.
  #'---------------------------------------------
  
  # Defines axis limits and labels
  
  if(deparse(substitute(baselayer))=="gmap.roebuck"){
    
    map.x <- seq(122.15,122.4,0.05)
    map.y <- seq(-18.15,-17.95,0.05)
    
    lab.x <- c(paste(map.x, "°E", sep = ""))
    lab.y <- c(paste(rev(-1*map.y), "°S", sep = ""))
    
  }else if(deparse(substitute(baselayer))=="gmap.kimb"){
    
    map.x <- seq(121,130,1)
    map.y <- seq(-19,-13,1)
    
    lab.x <- c(paste(map.x, "°E", sep = ""))
    lab.y <- c(paste(rev(-1*map.y), "°S", sep = ""))
    
  }else if(deparse(substitute(baselayer))=="gmap.cygnet"){
    
    map.x <- seq(122.9,123.7,0.2)
    map.y <- seq(-16.9,-16.1,0.1)
    
    lab.x <- c(paste(map.x, "°E", sep = ""))
    lab.y <- c(paste(rev(-1*map.y), "°S", sep = ""))
    
  }else if(deparse(substitute(baselayer))=="gmap.prince"){
    
    map.x <- seq(124.7,125.5,0.1)
    map.y <- seq(-15.7,-15,0.1)
    
    lab.x <- c(paste(map.x, "°E", sep = ""))
    lab.y <- c(paste(rev(-1*map.y), "°S", sep = ""))
    
  }else if(deparse(substitute(baselayer))=="gmap.cambridge"){
    
    map.x <- seq(127.8,128.6,0.2)
    map.y <- seq(-15.5,-14.7,0.2)
    
    lab.x <- c(paste(map.x, "°E",sep = ""))
    lab.y <- c(paste(rev(-1*map.y), "°S", sep = ""))
    
  }
  
  ggmap(baselayer)+
    
    geom_point(data = snub,
               aes(x = longitude, y = latitude, fill = sighting_class),
               colour = "black",
               pch = 21,
               size = pts.size, 
               alpha = pts.alpha)+
    
    coord_equal() + # Needs to be ### to produce map in right dimension pair
    
    theme_sleek() + # ggsidekick magic happens here
    
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 13),
          axis.text.x = element_text(size = 13),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
          legend.key = element_rect(fill = "transparent"),
          legend.position = c(0.1, 0.1),
          legend.background = element_rect(fill = "transparent", size = 2),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 13, face = "bold", colour = "black"),
          legend.key.size = unit(0.75,"cm")) +
    
    # Colour scale  
    
    scale_fill_manual(values = c(col.secondary, col.primary))+
    
    xlab("")+
    ylab("")+
    
    scale_x_continuous(limits = range(map.x), 
                       breaks = map.x,
                       labels = lab.x, expand = c(0,0))+
    
    scale_y_continuous(limits = range(map.y), 
                       breaks = map.y,
                       labels = rev(lab.y), expand = c(0,0))
  
  
}


# Convenience functions ---------------------------------------------------

#'---------------------------------------------
# Function to add leading zeroes
#'---------------------------------------------

addlzero <- function(string, n = 4){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param string Input string.
  #' @param n Number of leading zeroes.
  #'---------------------------------------------
  
  sapply(X = string, FUN = function(x){
    if(x%%1==0) res <- stringr::str_pad(string = x, width = n, pad = "0")
    if(x%%1>0) {
      a <- do.call(rbind, strsplit(as.character(x),"\\."))
      b <- stringr::str_pad(string = a[,1], width = n, pad = "0")
      res <- paste0(c(b, a[,2]), collapse = "-")}
    return(res)})
}

#'---------------------------------------------
# Function to rescale a raster object to 0-1 range
#'---------------------------------------------

rescale_raster <- function(r){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param r Input raster.
  #'---------------------------------------------
  
  (r-raster::cellStats(r,"min"))/(raster::cellStats(r,"max")-raster::cellStats(r,"min"))}

#'---------------------------------------------
# Function to calculate a percentile CI
#'---------------------------------------------

bci <- function(dat, perc = 95){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param dat Input values.
  #' @param perc Width of percentile confidence interval. Defauls to 0.95 for 95% CI.
  #'---------------------------------------------
  
  dat.mean <- mean(dat, na.rm = TRUE)
  dat.median <- median(dat, na.rm = TRUE)
  
  dat.low <- quantile(sort(dat), (100-perc)/(2*100)) 
  dat.high <- quantile(sort(dat), 1-(100-perc)/(2*100)) 

  res <- c(dat.mean, dat.median, dat.low, dat.high)
  names(res)[1] <- "mean"
  names(res)[2] <- "median"
  return(res)
}

#'---------------------------------------------
# Functions to convert column types
#'---------------------------------------------

chr2fac <- function(y){y %>% dplyr::mutate_if(is.character, as.factor)}
fac2chr <- function(y){y %>% dplyr::mutate_if(is.factor, as.character)}