#Marmot Data
setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

library(rgdal)
library(raster)

#Add Random Points (unused habitat) generated earlier.  50,000 random points 500m from transects. 
RandPtsDF <- shapefile("Analysis/RandomPoints/RandPtsTransect500mBuffer.shp")
plot(RandPtsDF)

#Get Marmot Data
Marfla <- shapefile("Locations/Transects/Marfla/MarflaCoordsSNV.shp")
class(Marfla)

#Extract values for MarmotStack layers at Random points
#This is going to take a long time with the 50,000 points
#Used Points generated within a 500mBuffer
#Available500mBuff <-extract(MarmotStack, RandPtsDF)
str(Available500mBuff)
head(Available500mBuff)
class(Available500mBuff)

#This takes 10 minutes, so use the csv you made last time.
#write.csv(Available500mBuff, "SNV/Available500mBuff", row.names = FALSE)
Available500mBuff <- read.csv("SNV/Available500mBuff") 

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
str(habvarsAvail500mBuff)


#The below Gets rid of observations below 2500.  This pulls out random points below there.  All transects are above this.  (around 15k points got pulled out)
elevfilter <- Marflahabvars.temp3$Elevation30SNV > 2500
Marflahabvars.temp3elevfilt <- Marflahabvars.temp3[elevfilter, ]
str(Marflahabvars.temp3elevfilt)

Marflahabvars.temp3 <- Marflahabvars.temp3elevfilt

#Extract values for MarmotStack layers at Marfla points (this takes 10min, skip to code below to read in csv you made last time)
#HabitatVariablesMarfla <- extract(MarmotStack, Marfla)
#class(HabitatVariablesMarfla)

#write.csv(HabitatVariablesMarfla, "SNV/HabitatVariablesMarfla", row.names = FALSE)
HabitatVariablesMarfla <- read.csv("SNV/HabitatVariablesMarfla")

#Need it in a data frame format
HabitatVariablesMarfla.dfm <- data.frame(HabitatVariablesMarfla)
str(HabitatVariablesMarfla.dfm)
colnames(HabitatVariablesMarfla.dfm) <- var.names
str(HabitatVariablesMarfla.dfm)

#Making it a spatial dataframe?
marfla.coords <- data.frame(coordinates(Marfla))
colnames(marfla.coords) <- c("UTMEW", "UTMNS")
str(marfla.coords)
nrow(unique(marfla.coords))

#Combine Used Points and UTM into one data frame. 
#habvarsMarfla <- cbind(marfla.coords, HabitatVariablesMarfla.dfm)
#str(habvarsMarfla)

#Get them into a CSV? 
#write.csv(habvarsMarfla, file = "HabitatVariablesMarfla.csv", row.names = FALSE)
habvarsMarfla <- read.csv("HabitatVariablesMarfla.csv")

str(habvarsMarfla)
class(habvarsMarfla)

##HELP - The following 10 lines of script was inherited, and I'm not sure what the relevance is of the temporary name assignments.   
Marflahabvars.temp <- read.csv("HabitatVariablesMarfla.csv", header=TRUE)
str (Marflahabvars.temp)

Marflahabvars.temp <- habvarsMarfla
str (Marflahabvars.temp)

Marflahabvars.temp2 <- Marflahabvars.temp
str(Marflahabvars.temp2)
###END HELP

# Combine used and available habitat 
Marflahabvars.temp3 <- rbind(Marflahabvars.temp2, habvarsAvail500mBuff)
str(Marflahabvars.temp3)

write.csv(Marflahabvars.temp3, file = "Used_And_Avail_500mBuff_HabitatVariablesMarfla.csv", row.names = FALSE)
