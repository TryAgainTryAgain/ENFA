#Get random points for the whole project area.  
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI/")

library(rgdal)
library(raster)

#==============================================================================
# Niche factor analysis based on spatial locations of XXX collected along 21
# transects in the Sierra Nevada from 2009 - 2012. In addition to XXX
# locations ("used points") 50,000 random points were generated as "available
# locations" from within the entire project boundary. 
# The available points then need to be filtered for elevations > 2500 m (about 8500
# feet). 
#==============================================================================

library(dismo)
#Read Project Boundary File into R
A <- readOGR("SNV/Boundary", layer = "SNVWorkingBoundary")
plot(A)
class(A)

# select 50000 random points
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)

RandPtsProj <- spsample(A, 50000, type = "random")
plot(RandPtsProj)

str(RandPtsProj)
head(RandPtsProj)

#####
class(RandPtsProj)
RandPtsProj

#write OGR requires a SpatialPointsDataFrame, not just SpatialPoints objet

#Here is the best way to convert a SpatialPoints object to a SpatialPointsDataFrame
spobj <- RandPtsProj #this is your SpatialPoints object
df <- data.frame(id=1:length(spobj)) #creating a very small data.frame 
spdf <- SpatialPointsDataFrame(spobj, data=df) #combining your data.frame and your SpatialPoints object
class(spdf) #check to see if the result really is a SpatialPointsDataFrame
head(spdf)
#Here is how you would save that object as a shapefile
library(rgdal)
writeOGR(spdf, dsn='Analysis/RandomPoints', layer='RandPtsProj',
         driver='ESRI Shapefile')

RandPtsProj <- shapefile("Analysis/RandomPoints/RandPtsProj.shp")
plot(RandPtsProj)
class(RandPtsProj)


#This works, but it saves it as a non-spatial object. Need to add spatial info. 
RandPtsProjDF <- as.data.frame(RandPtsProj)  
class(RandPtsProjDF)
head(RandPtsProjDF)
?coordinates
coordinates(RandPtsProjDF) <- c(1,2)
projection(RandPtsProjDF) <- CRS('+proj=utm +zone=11 +datum=NAD83')

str(RandPtsProjDF)
head(RandPtsProjDF)
class(RandPtsProjDF)
