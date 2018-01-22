#MarmotData

setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Get Marmot Data
#Marfla <- readOGR("Locations/Transects/Marfla", "MarflaCoordsSNV")
Marfla <- shapefile("Locations/Transects/Marfla/MarflaCoordsSNV.shp")
class(Marfla)

#Get Topographical Data 
Elevation <- raster("SNV/Topography/Elevation30SNV.tif") 
Elevation
Slope <- raster("SNV/Topography/Slope30SNV.tif") 
Slope
TRI <- raster("SNV/Topography/TRI30SNV.tif") 
#Productivity Data ("vegetation condition")
#Productivity across growing season
##NDVICV is a measure of the variation of this productivity
NDVICV<- raster("SNV/NDVI/NDVICVMeanMax1989_2015SNV.tif") 
NDVICV
#Max productivity
NDVIMeanMax <- raster("SNV/NDVI/NDVIMeanMax1989_2015SNV.tif") 
NDVIMeanMax

#Get Climate Data
#GrowingSeasonRainfall
Precip0609 <- raster("SNV/Climate/Precip/precip0609SNV.tif")
Precip0609

#Non-GrowingSeasonPrecip
Precip1005 <- raster("SNV/Climate/Precip/precip1005SNV.tif")
Precip1005
#Min Jan Temp
#Note, the January data layers only had Tmin for 07 (July), so used a file from Oct.  Checked the dimensions against Elevation, and they are the same.  
Tmin01<- raster("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_UpdatedNDVIlayerOct2017/AlpineMammals/NicheAnalysis/HabitatData/Climate/Temp/tmin01.tif") 

#Max July Temp
Tmax07 <- raster("SNV/Climate/Tmax/Tmax07SNV.tif") 
Tmax07
dim(Tmax07)
#Snow cover (<15%)
Snow <- raster("SNV/Climate/Snow/SnowFreeDaysSNV.tif")
dim(Snow)
#Habitat
#Note: January data no longer had tif files for this, use shapefiles instead? #Do NOT USE alpmamveg <- readOGR("SNV/Vegetation/Shapefiles", "alpmamveg") - it causes the computer to freeze.  But the rasters were derived from that file.  

Aspen <- raster("SNV/Vegetation/Rasterveg/aspen30m.img") 
dim(Aspen)
dim(Elevation)
Conifer <- raster("SNV/Vegetation/Rasterveg/conifer30m.img") 
Meadow <- raster("SNV/Vegetation/Rasterveg/meadow30m.img")
Mixed <- raster("SNV/Vegetation/Rasterveg/mixed30m.img") 
Rock <- raster("SNV/Vegetation/Rasterveg/rock30m.img") 
Shrub <- raster("SNV/Vegetation/Rasterveg/shrub30m.img") 

# Give all the rasters the same resolution (pixel size) and extent.
#Matching to elevation, because we know that elevation has the pixel size and extent we want for the project. 
# Try to do this with setExtent, with keepres=TRUE.  Likely that it will give phantom rows/columns.  This will result in the rasterStack working, but then Failure at layer IO when I tried to extract.  So tried resample, resammple is working. Resample does distort to fit, so don't use if the dimension difference is too large. 

# Give all rasters the same spatial extent
elevation <- setExtent(Elevation, Elevation, keepres=TRUE)
tri <- setExtent(TRI, Elevation, keepres=TRUE)
slope <- setExtent(Slope, Elevation, keepres=TRUE)
precip1005 <- setExtent(Precip1005, Elevation, keepres=TRUE)
precip0609 <- setExtent(Precip0609, Elevation, keepres=TRUE)
tmax07 <- setExtent(Tmax07, Elevation, keepres=TRUE)
tmin01 <- setExtent(Tmin01, Elevation, keepres=TRUE)
snow <- setExtent(Snow, Elevation, keepres=TRUE)
ndvimeanmax <- resample(NDVIMeanMax, Elevation)
ndvicv <- resample(NDVICV, Elevation)
aspen.0 <- extend(Aspen,Elevation)
aspen <- setExtent(aspen.0, Elevation, keepres=TRUE)
conifer.0 <- extend(Conifer,Elevation)
conifer <- setExtent(conifer.0, Elevation, keepres=TRUE)
meadow.0 <- extend(Meadow,Elevation)
meadow <- setExtent(meadow.0, Elevation, keepres=TRUE)
mixed.0 <- extend(Mixed,Elevation)
mixed <- setExtent(mixed.0, Elevation, keepres=TRUE)
rock.0 <- extend(Rock,Elevation)
rock <- setExtent(rock.0, Elevation, keepres=TRUE)
shrub.0 <- extend(Shrub,Elevation)
shrub <- setExtent(shrub.0, Elevation, keepres=TRUE)
e1d.0 <- extend(E1d,Elevation)
e1d <- setExtent(e1d.0, Elevation, keepres=TRUE)
types.0 <- extend(Types,Elevation)
types <- setExtent(types.0, Elevation, keepres=TRUE)

# Compare the rasters and make sure they have the same structure (extent,
# dimensions, projection, resolution, etc.)

compareRaster(c(elevation,tri,slope,precip1005,precip0609,tmin01,tmax07,snow,ndvimeanmax,ndvicv,aspen,conifer,meadow,mixed,rock,shrub))

#Raster Stack - can stack all the variable files
#Raster package to extract
#library(raster)

MarmotStack <- stack (elevation, slope, tri, ndvicv, ndvimeanmax, precip0609, precip1005, tmin01, tmax07, snow, aspen, conifer, meadow, mixed, rock,shrub)
##Before doing resample, with setExtent, Ran into problem due to different number or columns.  It is the vegetation layers that aren't lining up. 
crs(MarmotStack)
class(MarmotStack)
nlayers(MarmotStack)

names(MarmotStack)
# Function to make a vector of the variable names and capitalize each one

proper <- function(x) {
  s <- strsplit(x," ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

var.names.x <- names(MarmotStack)
var.names <- var.names.x
for(i in 1:length(var.names.x)) {
  var.names.i <- proper(var.names.x[i]) 
  var.names[i] <- var.names.i
}
var.names


names(MarmotStack) <- var.names
MarmotStack

#Get Available habitat for comparison.
#Extract values for MarmotStack layers at Random points
#This is going to take a long time with the 50,000 points
#Used Points generated within a 500mBuffer.  Script is in "Analysis/RandPoints500mBuffer.R"
Available500mBuff <-extract(MarmotStack, RandPtsDF)
str(Available500mBuff)
head(Available500mBuff)

#The extract didn't work, Failure during raster IO.  Possibly due to the phantom rows.  