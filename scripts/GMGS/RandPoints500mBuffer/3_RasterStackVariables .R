#Create a raster stack from which to extract variable values at used (observations) and available (random points) points.  
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Read saved rasters that have been cropped and resampled layers to match extent. 
elevation <- raster("SNV/ResampledRasterFiles/elevation30x30")
tri <- raster("SNV/ResampledRasterFiles/tri30x30")
slope <- raster("SNV/ResampledRasterFiles/slope30x30")
ppt0609 <- raster("SNV/ResampledRasterFiles/ppt060930x30")
ppt1005 <- raster("SNV/ResampledRasterFiles/ppt100530x30")
tmin01 <- raster("SNV/ResampledRasterFiles/tmin0130x30")
tmax07 <- raster("SNV/ResampledRasterFiles/tmax0730x30")
snow <- raster("SNV/ResampledRasterFiles/snow30x30")
ndvimeanmax <-raster("SNV/ResampledRasterFiles/ndvimeanmax30x30")
ndvicv <- raster("SNV/ResampledRasterFiles/ndvicv30x30")
aspen <-raster("SNV/ResampledRasterFiles/aspen30x30")
conifer <-raster("SNV/ResampledRasterFiles/conifer30x30")
meadow <-raster("SNV/ResampledRasterFiles/meadow30x30")
mixed <-raster("SNV/ResampledRasterFiles/mixed30x30")
rock<- raster("SNV/ResampledRasterFiles/rock30x30")
shrub <-raster("SNV/ResampledRasterFiles/shrub30x30")
# Compare the rasters and make sure they have the same structure (extent,
# dimensions, projection, resolution, etc.)

compareRaster(c(elevation,tri,slope,ppt1005,ppt0609,tmin01,tmax07,snow,ndvimeanmax,ndvicv,aspen,conifer,meadow,mixed,rock,shrub))

compareRaster(elevation, shrub)
dim(elevation)
dim(tri)
dim(slope)
dim(ppt1005)
dim(ppt0609)
dim(tmin01)
dim(tmax07)
dim(snow)
dim(ndvimeanmax)
dim(ndvicv)
dim(aspen)
dim(conifer)
dim(meadow)
dim(mixed)
dim(rock)
dim(shrub)

#Raster Stack - can stack all the variable files
#Raster package to extract
#library(raster)

GMGSStack <- stack (elevation, slope, tri, ndvicv, ndvimeanmax, ppt0609, ppt1005, tmin01, tmax07, snow, aspen, conifer, meadow, mixed, rock, shrub)

nlayers(GMGSStack)

names(GMGSStack)

# Function to make a vector of the variable names and capitalize each one
proper <- function(x) {
  s <- strsplit(x," ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

var.names.x <- names(GMGSStack)
var.names <- var.names.x
for(i in 1:length(var.names.x)) {
  var.names.i <- proper(var.names.x[i]) 
  var.names[i] <- var.names.i
}
var.names


names(GMGSStack) <- var.names
GMGSStack
