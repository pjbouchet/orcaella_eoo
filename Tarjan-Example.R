##apply phre function to sea otter data
##Max Tarjan; ltarjan@ucsc.edu
##July 29, 2014

##load data; specify file paths to supporting information 9 data
data<-read.csv("~/Dropbox/S9_File/BSR/F4_994_1346_data.xls"); 
data<-cbind(data$Xcoord, data$Ycoord)
grid<-read.csv("~/Dropbox/S9_file/SeaOtter_Grid_Array_100m.csv")

##create list of rasters from the grid array
library(raster)
rasters<-list(0)
rasters[[1]]<-rasterFromXYZ(cbind(grid$TealeX, grid$TealeY, grid$ATOScal))
rasters[[2]]<-rasterFromXYZ(cbind(grid$TealeX, grid$TealeY, log(grid$Distance+0.01)))

##specify the baseline smoothing parameters
smoother<-cbind(c(1,0), c(0,0.3)) #(ATOS, 0) (0, log(dist))

##adapt the ATOS smoothing parameter based on the nearest neighbor distances between re-sight locations
library(amap)
dist<-as.matrix(Dist(data))
nn<-dim(0)
for (c in 1:dim(dist)[2]) {
	min.temp<-min(dist[which(dist[,c]!=0),c])
	nn<-c(nn,min.temp)
}
mean.nn<-mean(log(nn))
smoother[1,1]<-smoother[1,1]*2*(mean.nn/4)^2.5

##apply phre function
HR<-phre(locs=data, rast=rasters, smooth=smoother, percent=90)

##plot of phre list objects
##zoomed in on polygons
plot(HR$Poly); plot(HR$array, add=T); plot(HR$Polygon, col="transparent", border='red', add=T); points(HR$locs, pch=20, cex=.1, col="black")