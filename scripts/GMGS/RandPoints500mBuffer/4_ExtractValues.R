#GMGS Data
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Add Random Points (unused habitat) generated earlier.  50,000 random points 500m from transects. 
RandPtsDF <- shapefile("Analysis/RandomPoints/RandPtsTransect500mBuffer.shp")
plot(RandPtsDF)

#Get GMGS Data
Callat <- shapefile("Locations/Transects/Spelat/SpelatCoordsSNV.shp")
class(Callat)

#Extract values for GMGSStack layers at Random points
#This is going to take a long time with the 50,000 points
#Used Points generated within a 500mBuffer
#Available500mBuff <-extract(GMGSStack, RandPtsDF)
str(Available500mBuff)
head(Available500mBuff)
class(Available500mBuff)

#This takes 10 minutes, so use the csv you made last time.
#write.csv(Available500mBuff, "SNV/Available500mBuff", row.names = FALSE)
Available500mBuff <- read.csv("SNV/Available500mBuff") 

#There is a problem combining the available (Available500mBuff) and the used(Callathabvars).  I think it is because the Available is a matrix, the used is a data.frame, and I haven't added the coordinates like I did for the used - resulting in different column numbers.  So attempting to use the code to add UTM to the available.

Available500mBuff.dfm <- data.frame(Available500mBuff)
str(Available500mBuff.dfm)
colnames(Available500mBuff.dfm) <- var.names
str(Available500mBuff.dfm)

Avail.coords <- data.frame(coordinates(RandPtsDF))
colnames(Avail.coords) <- c("UTMEW", "UTMNS")
str(Avail.coords)
nrow(unique(Avail.coords))

habvarsAvail500mBuff <- cbind(Avail.coords, Available500mBuff.dfm)
str(habvarsAvail500mBuff)

#Extract values for GMGSStack layers at Callat points (this takes 10min, skip to code below to read in csv you made last time)
#HabitatVariablesCallat <- extract(GMGSStack, Callat)
#class(HabitatVariablesCallat)

#write.csv(HabitatVariablesCallat, "SNV/HabitatVariablesCallat", row.names = FALSE)
HabitatVariablesCallat <- read.csv("SNV/HabitatVariablesCallat")

#Need it in a data frame format
HabitatVariablesCallat.dfm <- data.frame(HabitatVariablesCallat)
str(HabitatVariablesCallat.dfm)
colnames(HabitatVariablesCallat.dfm) <- var.names
str(HabitatVariablesCallat.dfm)

#Making it a spatial dataframe?
callat.coords <- data.frame(coordinates(Callat))
colnames(callat.coords) <- c("UTMEW", "UTMNS")
str(callat.coords)
nrow(unique(callat.coords))

#Combine Used Points and UTM into one data frame. 
#habvarsCallat <- cbind(callat.coords, HabitatVariablesCallat.dfm)
#str(habvarsCallat)

#Get them into a CSV? 
#write.csv(habvarsCallat, file = "HabitatVariablesCallat.csv", row.names = FALSE)
habvarsCallat <- read.csv("HabitatVariablesCallat.csv")

str(habvarsCallat)
class(habvarsCallat)

##HELP - The following 10 lines of script was inherited, and I'm not sure what the relevance is of the temporary name assignments.   
Callathabvars.temp <- read.csv("HabitatVariablesCallat.csv", header=TRUE)
str (Callathabvars.temp)

Callathabvars.temp <- habvarsCallat
str (Callathabvars.temp)

Callathabvars.temp2 <- Callathabvars.temp
str(Callathabvars.temp2)
###END HELP

# Combine used and available habitat 
Callathabvars.temp3 <- rbind(Callathabvars.temp2, habvarsAvail500mBuff)
str(Callathabvars.temp3)

#The below Gets rid of observations below 2500.  This pulls out random points below there.  All transects are above this.  (around 15k points got pulled out)
elevfilter <- Callathabvars.temp3$Elevation30SNV > 2500
Callathabvars.temp3elevfilt <- Callathabvars.temp3[elevfilter, ]
str(Callathabvars.temp3elevfilt)

Callathabvars.temp3 <- Callathabvars.temp3elevfilt

write.csv(Callathabvars.temp3, file = "Used_And_Avail_500mBuff_HabitatVariablesCallat.csv", row.names = FALSE)
