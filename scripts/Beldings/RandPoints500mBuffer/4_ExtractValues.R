#Beldings Data
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Add Random Points (unused habitat) generated earlier.  50,000 random points 500m from transects. 
RandPtsDF <- shapefile("Analysis/RandomPoints/RandPtsTransect500mBuffer.shp")
plot(RandPtsDF)

#Get Beldings Data
Urobel <- shapefile("Locations/Transects/Spebel/SpebelCoordsSNV.shp")
class(Urobel)

#Extract values for BeldingsStack layers at Random points
#This is going to take a long time with the 50,000 points
#Used Points generated within a 500mBuffer
#Available500mBuff <-extract(BeldingsStack, RandPtsDF)
str(Available500mBuff)
head(Available500mBuff)
class(Available500mBuff)

#This takes 10 minutes, so use the csv you made last time.
#write.csv(Available500mBuff, "SNV/Available500mBuff", row.names = FALSE)
Available500mBuff <- read.csv("SNV/Available500mBuff") 

#There is a problem combining the available (Available500mBuff) and the used(Urobelhabvars).  I think it is because the Available is a matrix, the used is a data.frame, and I haven't added the coordinates like I did for the used - resulting in different column numbers.  So attempting to use the code to add UTM to the available.

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

#Extract values for BeldingsStack layers at Urobel points (this takes 10min, skip to code below to read in csv you made last time)
#HabitatVariablesUrobel <- extract(BeldingsStack, Urobel)
#class(HabitatVariablesUrobel)

#write.csv(HabitatVariablesUrobel, "SNV/HabitatVariablesUrobel", row.names = FALSE)
HabitatVariablesUrobel <- read.csv("SNV/HabitatVariablesUrobel")

#Need it in a data frame format
HabitatVariablesUrobel.dfm <- data.frame(HabitatVariablesUrobel)
str(HabitatVariablesUrobel.dfm)
colnames(HabitatVariablesUrobel.dfm) <- var.names
str(HabitatVariablesUrobel.dfm)

#Making it a spatial dataframe?
urobel.coords <- data.frame(coordinates(Urobel))
colnames(urobel.coords) <- c("UTMEW", "UTMNS")
str(urobel.coords)
nrow(unique(urobel.coords))

#Combine Used Points and UTM into one data frame. 
#habvarsUrobel <- cbind(urobel.coords, HabitatVariablesUrobel.dfm)
#str(habvarsUrobel)

#Get them into a CSV? 
#write.csv(habvarsUrobel, file = "HabitatVariablesUrobel.csv", row.names = FALSE)
habvarsUrobel <- read.csv("HabitatVariablesUrobel.csv")

str(habvarsUrobel)
class(habvarsUrobel)

##HELP - The following 10 lines of script was inherited, and I'm not sure what the relevance is of the temporary name assignments.   
Urobelhabvars.temp <- read.csv("HabitatVariablesUrobel.csv", header=TRUE)
str (Urobelhabvars.temp)

Urobelhabvars.temp <- habvarsUrobel
str (Urobelhabvars.temp)

Urobelhabvars.temp2 <- Urobelhabvars.temp
str(Urobelhabvars.temp2)
###END HELP

# Combine used and available habitat 
Urobelhabvars.temp3 <- rbind(Urobelhabvars.temp2, habvarsAvail500mBuff)
str(Urobelhabvars.temp3)

#The below Gets rid of observations below 2500.  This pulls out random points below there.  All transects are above this.  (around 15k points got pulled out)
elevfilter <- Urobelhabvars.temp3$Elevation30SNV > 2500
Urobelhabvars.temp3elevfilt <- Urobelhabvars.temp3[elevfilter, ]
str(Urobelhabvars.temp3elevfilt)

Urobelhabvars.temp3 <- Urobelhabvars.temp3elevfilt

write.csv(Urobelhabvars.temp3, file = "Used_And_Avail_500mBuff_HabitatVariablesUrobel.csv", row.names = FALSE)
