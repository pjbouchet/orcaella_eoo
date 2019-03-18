## R code by César Capinha(1, 2) and Beatriz Pateiro-Lopez(3)
# (1) Centro de Biologia Ambiental, Faculdade de Ciências da Universidade de Lisboa, Campo Grande 1749-016 Lisboa, Portugal.
#(2) ZoologischesForschungsmuseum Alexander Koenig, Museumsmeile Bonn, Adenauerallee 160, 53113 Bonn, Germany
# (3) Departamento de Estatística e Investigación Operativa, Facultade de Matemáticas, Universidade de Santiago de Compostela, Rúa Lope Gómez de Marzoa s/n, 15782 Santiago de Compostela, Spain.


## DESCRIPTION
# Supplementary R code of the manuscript: "Capinha C., Pateiro-Lopez B. (submitted) Predicting species distributions in new areas or time periods with alpha-shapes" 
# 30th of January 2014. 


library(alphashape3d) 
library(raster) 

## DATA SET UP
# In our case study we use four data inputs: 1) a ".csv" file having the latitude and longitude of each species' occurrence; 2) a ".csv" file having the latitude and longitude of the areas that were surveyed (i.e. the sampling coverage); 3) a ".csv" file with the latitude and longitude of locations in which the species was absent and 4) a set of six climatic predictors in ASCII raster format. 

#Load ".cvs" files
setwd("D://Data_Xlaevis") # Set working directory
occur.data <- read.table("xlaevis_occ.csv", header = TRUE, sep = ",", na.strings = "NA", dec = ".")
cover.data <- read.table("frogmap_cover.csv", header = TRUE, sep = ",", na.strings = "NA", dec = ".")
abs.data <- read.table("xlaevis_abs.csv", header = TRUE, sep = ",", na.strings = "NA", dec = ".")

#Load climatic predictors
grids <- list.files(pattern = "asc", full.names = T) #Load all predictors in working folder
predictors <- stack(grids)

#Reduce the dimensionality of the climatic space to 3D using a Principal Components Analysis (PCA)
pca <- princomp(na.omit(values(predictors)), cor = TRUE)

#Scores of the three first axes of the PCA
PCs <- predict(predictors, pca, index = 1:3)
occur.PCs <- extract(PCs, occur.data)
abs.PCs <- extract(PCs, abs.data)
cover.PCs <- extract(PCs, cover.data)
PCs.points <- rasterToPoints(PCs)

#Visualize occurrences and absences in the 3D space (optional)
plot3d(occur.PCs, col = "blue", type = "s", size = 0.5)
plot3d(abs.PCs, col = "orange", type = "s", size = 0.5, add = TRUE)

##SELECTION OF AN ALPHA-SHAPE
#First method: Alpha-shape achieving highest discrimination between species presences and absences. The vector "values" will store the true skill statistics for each value of alpha.
i <- seq(0, 0.2, by = 0.001)
ashape3d.occ <- ashape3d(unique(as.matrix(occur.PCs)), alpha = i)
result.occurrence <- inashape3d(ashape3d.occ, indexAlpha = "all", points = occur.PCs)
result.absence <- inashape3d(ashape3d.occ, indexAlpha = "all", points = abs.PCs)
sensitivity <- unlist(lapply(result.occurrence, "mean"))
specificity <- 1 - unlist(lapply(result.absence, "mean"))
TSS <- sensitivity + specificity - 1
values <- cbind(i, TSS)

#Second method: the minimum bounding envelope (MBE). Tests the percentage of occurrences that are bounded by the envelope. The vector "values" will store the results of this test. 
i <- seq(0, 0.4, by = 0.001)
ashape3d.occ <- ashape3d(unique(as.matrix(occur.PCs)), alpha = i)
result.occurrence <- inashape3d(ashape3d.occ, indexAlpha = "all", points = as.matrix(occur.PCs))
percentage <- unlist(lapply(result.occurrence, "mean"))
values <- data.frame(i, percentage * 100)

#Visualize the MBE in search of holes or voids.
alpha <- 0.292
ashape3d.occ <- ashape3d(unique(as.matrix(occur.PCs)), alpha = alpha)
plot(ashape3d.occ, bycom = TRUE, tran = 0.2, shininess = 0)
rgl.bbox(color = c("#333377", "white"), emission = "#333377", specular = "#3333FF", shininess = 5, alpha = 0.8)

#Save prediction to geographic space
inside.occur.shape <- inashape3d(ashape3d.occ, points = PCs.points[, 3:5])
occur.pred <- cbind(PCs.points[, 1:2], inside.occur.shape)
raster.result <- rasterize(x = occur.pred[, 1:2], y = predictors[[1]], field = occur.pred[, 3])
writeRaster(raster.result, filename = "Mapped_prediction_a0292.asc", overwrite = TRUE)

##BOOTSTRAP SUPPORT (it may take a while)
R <- 100
boot.supp <- matrix(nrow = nrow(PCs.points), ncol = R)
for (i in 1:R) {
  boot.sample <- occur.PCs[sample(nrow(occur.PCs), nrow(occur.PCs), replace = TRUE), ]
  ashape3d.boot <- ashape3d(unique(as.matrix(boot.sample)), alpha = alpha)
  boot.supp[, i] <- inashape3d(ashape3d.boot, points = PCs.points[, 3:5])
}

#Save bootstrap result to geographic space.
boot.results <- cbind(PCs.points[, 1:2], rowSums(boot.supp))
raster.boot <- rasterize(x = boot.results[, 1:2], y = predictors[[1]], field = boot.results[, 3])
writeRaster(raster.boot, filename = "Bootstrap_a292.asc", overwrite = TRUE)

##IDENTIFICATION OF ANALOG AND NON-ANALOG CONDITIONS
#The calculation is identical to the MBE (see above) but for the full range of climatic conditions that were sampled. 
i <- seq(0.2, 0.3, by = 0.001)
ashape3d.cover <- ashape3d(unique(as.matrix(cover.PCs)), alpha = i)
result.cover <- inashape3d(ashape3d.cover, indexAlpha = "all", points = cover.PCs)
percentage <- unlist(lapply(result.cover, "mean"))
values <- data.frame(i, percentage * 100)

#Visualize the climatic envelope in search of holes or voids.
alpha <- 0.292
ashape3d.coverage <- ashape3d(unique(as.matrix(cover.PCs)), alpha = alpha)
plot(ashape3d.coverage, bycom = TRUE, tran = 0.2, shininess = 0)
rgl.bbox(color = c("#333377", "white"), emission = "#333377", specular = "#3333FF", shininess = 5, alpha = 0.8)

#Identify analog and non-analog conditions in the geographic space
inside.cover.shape <- inashape3d(ashape3d.coverage, points = PCs.points[, 3:5])
cover.pred <- cbind(PCs.points[, 1:2], inside.cover.shape)
raster.result.cover <- rasterize(x = cover.pred[, 1:2], y = predictors[[1]], field = cover.pred[, 3])
writeRaster(raster.result.cover, filename = "Analog_conditions.asc", overwrite = TRUE)
