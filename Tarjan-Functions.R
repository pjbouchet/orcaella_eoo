##code to write the phre function for generalized layers
##Max Tarjan; ltarjan@ucsc.edu
##June 26, 2015

##requires that packages are installed: raster, ks, amap, maptools
##run all code in this file and then apply function to arguments

phre<-function(locs, rast, smoother, percent){ # phre is a function that uses these arguments
  
  library(raster) #required for rast conversion & extract
  library(ks) #required for kde
  library(amap) #required for dist function
  library(maptools) #required for checkpolygonholes function
  
  landscape<-dim(0) #empty matrix for landscape variables
  
  ##define dimensions for area where the kd will be calculated
  x.range<-max(locs[,1])-min(locs[,1])
  y.range<-max(locs[,2])-min(locs[,2])
  
  min.x<-min(locs[,1])-x.range*.5
  max.x<-max(locs[,1])+x.range*.5
  min.y<-min(locs[,2])-y.range*.5
  max.y<-max(locs[,2])+y.range*.5
  
  ##create coordinates for grid array (used to back-translate kde estimate)
  seq.x<-seq(min.x, max.x, (max.x-min.x)/250)
  seq.y<-seq(min.y, max.y, (max.x-min.x)/250)
  
  array<-dim(0)
  
  for (a in 1:length(seq.x)){
    array<-rbind(array, cbind(rep(seq.x[a],length(seq.y)), seq.y))
  }
  
  ##transform resight locations to landscape variables; also extract landscape variables to grid array
  for (r in 1:length(rast)){ ##for each raster layer
    landscape<-cbind(landscape, raster::extract(rast[[r]], locs[,1:2], na.rm=F, method='bilinear'))##calculate the raster values at locs
    array<-cbind(array, raster::extract(rast[[r]], array[,1:2], na.rm=F, method='bilinear'))##add raster values to array
  }
  
  if (length(is.na(landscape[,1]))>length(landscape)/2) {
    print('Warning: Raster layers do not cover the extent of the location data')
  }
  
  ##calculate kernel density values at each point in the grid array; kde function is generated using landscape variables of re-sight locations
  
  if (smoother[1]=='default') {
    kd<-kde(x=na.omit(landscape), eval.points=array[,3:dim(array)[2]])
  } else {
    kd<-kde(x=na.omit(landscape), H=smoother, eval.points=array[,3:dim(array)[2]])
  }
  
  smoother<-kd$H
  density.array<-kd$estimate
  
  ##transform density values to probability values
  z<-density.array/sum(density.array, na.rm=TRUE) ##make array sum to one
  array<-cbind(array,z) 
  
  ##make a matrix of array points that are within the 90% probability kernel
  array.order<-array[order(array[,dim(array)[2]], decreasing=TRUE),]
  if (exists('percent')) {
    percent<-percent/100
  } else {percent<-0.9}
  
  for (i in 1:length(array.order[,dim(array)[2]])) {
    if (sum(array.order[1:i,dim(array)[2]])>=percent) {
      HRpoints<-array.order[1:i,]
      critval<-array.order[i,dim(array)[2]]
      break
    }
  }
  
  ##create a raster layer of the grid array with probability values
  HR.grid<-rasterFromXYZ(cbind(array[,1:2],array[,dim(array)[2]]))
  ##create a raster layer of the grid array points within the 90% kernel
  HR.rast<-rasterFromXYZ(cbind(HRpoints[,1:2],rep(1,dim(HRpoints)[1])), res=abs(seq.y[1]-seq.y[2]))
  ##define the area of tiny fragment polygons (4 pixles or smaller)
  min.poly.size<-4*((HR.rast@extent@xmax-HR.rast@extent@xmin)/HR.rast@ncols)^2
  ##convert HR.rast to polygons, which denote the permissible home range
  HR.poly<-rasterToPolygons(HR.rast, n=16, na.rm=T, digits=4, dissolve=T)
  
  ##remove polygons that are 4 pixels or smaller
  s<-"polyID"
  HR.polys<-list(0)
  pr<-0
  for (p in 1:length(HR.poly@polygons[[1]]@Polygons)) {
    if (HR.poly@polygons[[1]]@Polygons[[p]]@area > min.poly.size) {
      pr<-pr+1
      HR.polys[[pr]]<-HR.poly@polygons[[1]]@Polygons[[p]]
    }
  }
  HR.polys<-Polygons(HR.polys, ID=s) #convert to polygons object
  HR.polys<-checkPolygonsHoles(HR.polys) #assign holes correctly
  HR.polys<-SpatialPolygons(list(HR.polys)) #convert to spatial polygons
  HR.polys<-gSimplify(HR.polys, tol=20)
  
  list(Polygon=HR.polys,locs=locs,HRpoints=HRpoints,array=HR.grid, smoother=smoother) ##outputs; list of four elements: 1) 90% kernel polygon, 2) location data used for hr estimate, 3) array points that fall within the 90% home range, 4)all points in the array with landscape and z values
}