#Marmot Data
#setwd("/Users/avivarossi/Desktop/AlpineMammalsRScript/My3Mammals/")
#setwd("C:/Users/rcklinger/Documents/Aviva/")
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Add Random Points (unused habitat) generated earlier.  50,000 random points 500m from transects. 


RandPtsDF <- shapefile("Analysis/RandomPoints/RandPtsTransect500mBuffer.shp")

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
#PRR <- raster("NicheAnalysis/HabitatData/Topography/prr.tif") - PotentialRelativeRadiation - this layer wasn't in Jan Layers. 
#Hillshade <- raster("NicheAnalysis/HabitatData/Topography/hillshade.tif")  - this layer wasn't in Jan Layers. 

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
dim(Tmin01)
dim(Elevation)
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
# Tried to do this with setExtent, with keepres=TRUE, but that likely gave phantom rows/columns.  This resulted in the rasterStack working, but then Failure at layer IO when I tried to extract.  But resammple is working. Resample does distort to fit, so don't use if the dimension difference is too large. 

#Note, I put the #in front of some of these because I saved the resampled versions.  Will need to undo this if you don't have access to the resampled versions, and put code in to load them below.  

elevation <- Elevation
dim(elevation)
slope <- Slope
dim(slope)
tri <- TRI
dim(TRI)
#prr <- PRR
#hillshade <- Hillshade
dim(NDVICV)
#ndvicv <- resample(NDVICV, elevation)
dim(NDVIMeanMax)
#ndvimeanmax <- resample(NDVIMeanMax, elevation)
dim(Precip0609)
#ppt0609 <- resample(Precip0609, elevation)
dim(Precip1005)
#ppt1005 <- resample(Precip1005, elevation)
#tmin01 <- resample(Tmin01, elevation)
#tmax07 <- resample(Tmax07, elevation)
snow <- Snow
#aspen <- resample(Aspen, elevation)
#conifer <- resample(Conifer, elevation)
#meadow <- resample(Meadow, elevation)
#mixed <- resample(Mixed, elevation)
#rock <- resample(Rock, elevation)
#shrub <- resample(Shrub, elevation)
dim(Tmin01)
dim(snow)


#Save resampled files so that I don't have to resample each time I run this. 
#writeRaster(ndvicv, "SNV/ResampledRasterFiles/ndvicv30x30")
#writeRaster(ndvimeanmax, "SNV/ResampledRasterFiles/ndvimeanmax30x30")
#writeRaster(ppt0609, "SNV/ResampledRasterFiles/ppt060930x30")
#writeRaster(ppt1005, "SNV/ResampledRasterFiles/ppt100530x30")
#writeRaster(tmin01, "SNV/ResampledRasterFiles/tmin0130x30")
#writeRaster(tmax07, "SNV/ResampledRasterFiles/tmax0730x30")
#writeRaster(aspen, "SNV/ResampledRasterFiles/aspen30x30")
#writeRaster(conifer, "SNV/ResampledRasterFiles/conifer30x30")
#writeRaster(meadow, "SNV/ResampledRasterFiles/meadow30x30")
#writeRaster(mixed, "SNV/ResampledRasterFiles/mixed30x30")
#writeRaster(rock, "SNV/ResampledRasterFiles/rock30x30")
#writeRaster(shrub, "SNV/ResampledRasterFiles/shrub30x30")

#Then this is the code to read those saved resamples back in
ndvicv <- raster("SNV/ResampledRasterFiles/ndvicv30x30")
ndvimeanmax <-raster("SNV/ResampledRasterFiles/ndvimeanmax30x30")
tmin01 <- raster("SNV/ResampledRasterFiles/tmin0130x30")
tmax07 <- raster("SNV/ResampledRasterFiles/tmax0730x30")
ppt0609 <- raster("SNV/ResampledRasterFiles/ppt060930x30")
ppt1005 <- raster("SNV/ResampledRasterFiles/ppt100530x30")
aspen <- raster("SNV/ResampledRasterFiles/aspen30x30")
conifer <-raster("SNV/ResampledRasterFiles/conifer30x30")
meadow <-raster("SNV/ResampledRasterFiles/meadow30x30")
mixed <- ("SNV/ResampledRasterFiles/mixed30x30")
rock<- raster("SNV/ResampledRasterFiles/rock30x30")
shrub <- raster("SNV/ResampledRasterFiles/shrub30x30")

# Compare the rasters and make sure they have the same structure (extent,
# dimensions, projection, resolution, etc.)

compareRaster(c(elevation,tri,slope,ppt1005,ppt0609,tmin01,tmax07,snow,ndvimeanmax,ndvicv,aspen,conifer,meadow,mixed,rock,shrub))

crs(meadow)
dim(elevation)
dim(Aspen)
compareRaster(elevation, Aspen)
compareRaster(elevation, rock)


#Raster Stack - can stack all the variable files
#Raster package to extract
#library(raster)

MarmotStack <- stack (elevation, slope, tri, ndvicv, ndvimeanmax, ppt0609, ppt1005, tmin01, tmax07, snow, aspen, conifer, meadow, mixed, rock,shrub)
dim(elevation)
dim(slope)
dim(tri)
dim(ndvicv)
dim(ndvimeanmax)
dim(ppt0609)
dim(ppt1005)
dim(tmin01)
dim(tmax07)
dim(snow)
dim(aspen)
dim(conifer)
dim(meadow)
dim(mixed)

##Before doing resample, Ran into problem due to different number or columns.  It is the vegetation layers that aren't lining up. 
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

#Extract values for MarmotStack layers at Random points
#This is going to take a long time with the 50,000 points
#Used Points generated within a 500mBuffer
Available500mBuff <-extract(MarmotStack, RandPtsDF)
str(Available500mBuff)
head(Available500mBuff)
class(Available500mBuff)

#This takes 10 minutes, so use the csv you made last time.
#write.csv(Available500mBuff, "SNV/Available500mBuff", row.names = FALSE)
Available500mBuff <- read.csv("SNV/Available500mBuff") 
read.csv
 
#Extract values for MarmotStack layers at Marfla points (this takes 10min, skip to code below to read in csv you made last time)
HabitatVariablesMarfla <- extract(MarmotStack, Marfla)

class(HabitatVariablesMarfla)

#write.csv(HabitatVariablesMarfla, "SNV/HabitatVariablesMarfla", row.names = FALSE)
HabitatVariablesMarfla <- read.csv("SNV/HabitatVariablesMarfla")

HabitatVariablesMarfla.dfm <- data.frame(HabitatVariablesMarfla)
str(HabitatVariablesMarfla.dfm)
colnames(HabitatVariablesMarfla.dfm) <- var.names
str(HabitatVariablesMarfla.dfm)

marfla.coords <- data.frame(coordinates(Marfla))
colnames(marfla.coords) <- c("UTMEW", "UTMNS")
str(marfla.coords)
nrow(unique(marfla.coords))

habvarsMarfla <- cbind(marfla.coords, HabitatVariablesMarfla.dfm)
#str(habvarsMarfla)

#Get them into a CSV? 
#write.csv(habvarsMarfla, file = "HabitatVariablesMarfla.csv", row.names = FALSE)
habvarsMarfla <- read.csv("HabitatVariablesMarfla.csv")

str(habvarsMarfla)
#Niche Factor Analysis
library(adehabitatHS)
library(adehabitatMA)
library(vegan)
library(raster)
library(sp)
library(maptools)

#==============================================================================
# Niche factor analysis based on spatial locations of Marmot collected along 21
# transects in the Sierra Nevada from 2009 - 2012. In addition to 1703
# locations ("used points") 50,000 random points were generated as "available
# locations". Values of thirteen habitat variables were extracted for each point.
# The available points were then filtered for elevations > 2500 m (about 8500
# feet). 
#==============================================================================


Marflahabvars.temp <- read.csv("HabitatVariablesMarfla.csv", header=TRUE)
str (Marflahabvars.temp)
summary (Marflahabvars.temp)

#only use the below is unused already in dataframe


##HELP not sure what is happening here.  Looks like overwriting what we did above.  
Marflahabvars.temp <- habvarsMarfla
str (Marflahabvars.temp)
summary (Marflahabvars.temp)

Marflahabvars.temp2 <- Marflahabvars.temp
str(Marflahabvars.temp2)

# Available habitat
names(Available500mBuff)
names(Marflahabvars.temp2)
class(Available500mBuff)
head(Available500mBuff)
class(Marflahabvars.temp2)

#There is a problem combining the available (Available500mBuff) and the used(Marflahabvars).  I think it is because the Available is a matrix, the used is a data.frame, and I haven't added the coordinates like I did for the used - resulting in different column numbers.  So attempting to use the code to add UTM to the available.

Available500mBuff.dfm <- data.frame(Available500mBuff)
str(Available500mBuff.dfm)
colnames(Available500mBuff.dfm) <- var.names
str(Available500mBuff.dfm)

Avail.coords <- data.frame(coordinates(RandPtsDF))
colnames(Avail.coords) <- c("UTMEW", "UTMNS")
str(Avail.coords)
nrow(unique(Avail.coords))

habvarsAvail500mBuff <- cbind(Avail.coords, Available500mBuff.dfm)

# Combine used and available habitat 
Marflahabvars.temp3 <- rbind(Marflahabvars.temp2, habvarsAvail500mBuff)
str(Marflahabvars.temp3)

#The below Gets rid of observations below 2500.  This pulls out random points below there.  All transects are above this.  (around 15k points got pulled out)
elevfilter <- Marflahabvars.temp3$Elevation30SNV > 2500
Marflahabvars.temp3elevfilt <- Marflahabvars.temp3[elevfilter, ]
str(Marflahabvars.temp3elevfilt)

Marflahabvars.temp3 <- Marflahabvars.temp3elevfilt

# You need to generate a vector for "present" (i.e. observed; = 1) and random
# (available; = 0) and it needs to be the first field (column) of the dataframe
# Since you are going to remove the UTM coordinates you might as well just use
# one of those columns after renaming it and put the vector of 0's and 1's in
# it

Marflahabvars <- Marflahabvars.temp3
str(Marflahabvars)
names(Marflahabvars)[1] <- "Marfla"
Marflahabvars$Marfla <- 0
Marflahabvars$Marfla[1:nrow(habvarsMarfla)] <- 1
Marflahabvars$UTMNS <- NULL
str(Marflahabvars)

marfla.pa <- Marflahabvars$Marfla
marfla.pa
Marflahabvars$Marfla <- NULL
str(Marflahabvars)

#This gets rid of the information we no longer need (marmot siting ID, because variables already extracted)
str(Marflahabvars)

# Pool Aspen and Mixed land cover types with Shrub)
Marflahabvars$Shrub30m <- Marflahabvars$Shrub30m + Marflahabvars$Aspen30m + Marflahabvars$Mixed30m
#Get rid of columns for Aspen and Mixed after pooling into Shrub
Marflahabvars$Aspen30m <- NULL
Marflahabvars$Mixed30m <- NULL


cor(Marflahabvars)#checking how correlated variables are
str(Marflahabvars)

# Used habitat
# If using a CSV that has both used and unused, Go back and get used and unused from before we pulled that column out. 
# weights is something that AdeHabitatHS needs, possibly even be named this.  
weights <- marfla.pa
# weights is equivalent to pr in histniche example, it is the number of detections at that location

str(weights)

# Comparison of used (dark gray bars) and available (white bars) frequency for
# each habitat variable

#histniche is a function in AdeHabitatHS that lets you compare used vs available for each of the variables.  (Note to self - Manly operates off of categories, this is using continuous variables)

histniche(Marflahabvars, weights)

length(weights[weights == 0])

# Next step requires PCA first
# Initial PCA of the available habitat variables
#dudi is from ADE package, for doing PCA
habavail.pca <- dudi.pca(Marflahabvars, scannf=FALSE)

# Generalized niche factor analysis (used habitat is the Focus group)
#Of the three versions of Niche Factor Analysis, this may be the least useful, but valuable for me to see all three versions.  
#Doesn't separate marginality and specialization, synthesizes them. 
#uses PCA values we just did
#Focus is actual marmot locations, the used
#centering, center all variables so they are on the same scale. 

marmot.gn <- gnesfa(habavail.pca, Focus=weights, centering="single", nfFirst=1, nfLast=1, scannf=FALSE)

#nfFirst is marginality, 1 means include it.  nfLast is specialization, 1 means include it. 
#Look at eigenvalues.  Measure of the variation accounted for by that gradient.  If the second one is much lower than the first (like 37 to 2), only use the first.  Always use the last.  
marmot.gn
marmot.gn$cor

#Note, if you save the previous plots, make sure to close those plots, otherwise future plots won't work. 

scatterniche(marmot.gn$li, weights, pts=TRUE)
#can see Axis 1 is marginality, and Axis 2 is specialization

#These combine marginality and specialization
s.arrow(marmot.gn$co)
#raw values, doesn't account for correlations
s.arrow(marmot.gn$cor)

#the above accounts for correlation among variables, and is the one of these graphs that matters the most. 

# For the marmot, this is showing that they are associated with temperature (negatively with high T max, and positively with high T min).  The next important ones have to do with roughness and slope.  

# Mahalanobis analysis (available habitat is the Focus group)
# Way of measuring in multivariate space how different groups (used vs available are).  The focus right now is on the available (as opposed to the used in the last one)
marmot.mad <- madifa(habavail.pca, weights, scan=FALSE)
##Sometimes getting an error message - Error in eigen(W) : infinite or missing values in 'x'.  Is this because the Marfla is in the analysis.  Ask Rob if Marfla is supposed to be there, and if not, how to run the analysis without including Marfla as a factor.  
marmot.mad
?madifa

# li is like coordinates to plot in multivariate space
scatterniche(marmot.mad$li, weights, pts=TRUE)
s.arrow(marmot.mad$co)
s.arrow(marmot.mad$cor)

#Gray is available.  Black is used.  X axis is marginality, y axis is specialization.  

# Niche factor analysis (ENFA; used habitat is the Focus group)
#This is the third NFA, and the one Rob finds most useful.
#This one explicitly discriminates between marginality and specialization
marmot.enfa <- enfa(habavail.pca, weights, scan=FALSE)
marmot.enfa
# co shows how important that variable is loading for the 1st and 2nd.  The li is the actual coordinate of the variable in multivariate space.  
marmot.enfa$co
hist(marmot.enfa, type="h")
scatter(marmot.enfa)
scatter(marmot.enfa, nc=T, percent=95, clabel=0, Acol="green", Ucol="blue",
        side="none", Aborder="black", Uborder="black")
scatterniche(marmot.enfa$li, weights)
s.arrow(marmot.enfa$co)
marmot

#==============================================================================
# Graph of used vs. availability locations from the PCA
#==============================================================================

str(habavail.pca)
pcacols <- rep("green", nrow(habavail.pca$li))
pcacols[weights > 0] <- "blue"
ua.fac.avail <- rep("Available", length(weights))
ua.fac.used <- rep("Available", length(weights))
ua.fac.used[weights > 0] <- "Used"
plot(habavail.pca$li, pch=21, bg=pcacols, cex=0.5, xlim=c(-10,10), ylim=c(-9,11),
     xlab="Axis 1", ylab="Axis 2", cex.main=0.9,
     main="Marfla - Sierra Nevada Core Landscape\nPCA")
ordihull(habavail.pca$li, groups=ua.fac.avail, draw="lines", col="green",
         lwd=2)
ordihull(habavail.pca$li, groups=ua.fac.used, show.groups="Used", col="blue",
         draw="lines", lwd=3)
legend(7,0, legend=c("Available","Used"), pch=21, pt.bg=c("green","blue"),
       cex=0.9, box.lty=0)

