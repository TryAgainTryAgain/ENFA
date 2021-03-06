setwd("/Volumes/AvivasDissertation_2018_4158472461/Dissertation/GISData_Jan2018UpdatedNDVI") 

#Niche Factor Analysis
#GMGS Data
library(adehabitatHS)
library(adehabitatMA)
library(vegan)
library(raster)
library(sp)
library(maptools)

#==============================================================================
# Niche factor analysis based on spatial locations of GMGS collected along 21
# transects in the Sierra Nevada from 2009 - 2012. In addition to 1703
# locations ("used points") 50,000 random points were generated as "available
# locations". Values of thirteen habitat variables were extracted for each point.
# The available points were then filtered for elevations > 2500 m (about 8500
# feet). 
#==============================================================================


#Read in CSV with variables from extracted Raster Stack for Used (500m Buffer) and Available
Callathabvars.temp3 <- read.csv("Used_And_Avail_500mBuff_HabitatVariablesCallat.csv", header=TRUE)
str(Callathabvars.temp3)
head(Callathabvars.temp3)

# You need to generate a vector for "present" (i.e. observed; = 1) and random
# (available; = 0) and it needs to be the first field (column) of the dataframe
# Since you are going to remove the UTM coordinates you might as well just use
# one of those columns after renaming it and put the vector of 0's and 1's in
# it

Callathabvars <- Callathabvars.temp3
str(Callathabvars)
names(Callathabvars)[1] <- "Callat"
Callathabvars$Callat <- 0
Callathabvars$Callat[1:nrow(habvarsCallat)] <- 1
Callathabvars$UTMNS <- NULL
str(Callathabvars)
head(Callathabvars)

#This gets rid of the information we no longer need (gmgs siting ID, because variables already extracted)
callat.pa <- Callathabvars$Callat
callat.pa
Callathabvars$Callat <- NULL
str(Callathabvars)

# Pool Aspen and Mixed land cover types with Shrub, because there is not enough Aspen and Mixed is effectively shrub
Callathabvars$Shrub30m <- Callathabvars$Shrub30m + Callathabvars$Aspen30m + Callathabvars$Mixed30m
#Get rid of columns for Aspen and Mixed after pooling into Shrub
Callathabvars$Aspen30m <- NULL
Callathabvars$Mixed30m <- NULL


cor(Callathabvars)#checking how correlated variables are
str(Callathabvars)
str(callat.pa)
head(callat.pa)
tail(callat.pa)


# Used habitat
# If using a CSV that has both used and unused, Go back and get used and unused from before we pulled that column out. 
# weights is something that AdeHabitatHS needs, possibly even be named this.  
weights <- callat.pa
# weights is equivalent to pr in histniche example, it is the number of detections at that location

str(weights)

# Comparison of used (dark gray bars) and available (white bars) frequency for
# each habitat variable

#histniche is a function in AdeHabitatHS that lets you compare used vs available for each of the variables.  (Note to self - Manly operates off of categories, this is using continuous variables)

histniche(Callathabvars, weights)

#Checking how many of the weights is 0, should be 5,000 (the number of random points)
length(weights[weights == 0])

# Next step requires PCA first
# Initial PCA of the available habitat variables
#dudi is from ADE package, for doing PCA
habavail.pca <- dudi.pca(Callathabvars, scannf=FALSE)

# Generalized niche factor analysis (used habitat is the Focus group)
#https://www.researchgate.net/publication/5454530_A_general_framework_for_the_statistical_exploration_of_the_ecological_niche 
#Of the three versions of Niche Factor Analysis, this may be the least useful, but valuable for me to see all three versions.  
#Doesn't separate marginality and specialization, synthesizes them. 
#uses PCA values we just did
#Focus is actual gmgs locations, the used
#centering, center all variables so they are on the same scale. 
#nfFirst is marginality, 1 means include it.  nfLast is specialization, 1 means include it.
gmgs.gn <- gnesfa(habavail.pca, Focus=weights, centering="single", nfFirst=1, nfLast=1, scannf=FALSE)


#Look at eigenvalues.  Measure of the variation accounted for by that gradient.  If the second one is much lower than the first (like 37 to 2), only use the first.  Always use the last.  
gmgs.gn
gmgs.gn$cor

scatterniche(gmgs.gn$li, weights, pts=TRUE)
#can see Axis 1 is marginality, and Axis 2 is specialization

#These s.arrow plots combine marginality and specialization
#raw values, doesn't account for correlations
s.arrow(gmgs.gn$co)

#the below accounts for correlation among variables, and is the one of these graphs that matters the most. 
s.arrow(gmgs.gn$cor)


# For the gmgs, this is showing that they are associated with temperature (negatively with high T max, and positively with high T min).  The next important ones have to do with roughness and slope.  

# Mahalanobis analysis (available habitat is the Focus group)
# Way of measuring in multivariate space how different groups (used vs available are).  The focus right now is on the available (as opposed to the used in the last one)
gmgs.mad <- madifa(habavail.pca, weights, scan=FALSE)
  
gmgs.mad
?madifa

# li is like coordinates to plot in multivariate space
scatterniche(gmgs.mad$li, weights, pts=TRUE)
s.arrow(gmgs.mad$co)
s.arrow(gmgs.mad$cor)

#Gray is available.  Black is used.  X axis is marginality, y axis is specialization.  

# Niche factor analysis (ENFA; used habitat is the Focus group)
#This is the third NFA, and the one Rob finds most useful.
#This one explicitly discriminates between marginality and specialization
gmgs.enfa <- enfa(habavail.pca, weights, scan=FALSE)
gmgs.enfa
str(gmgs.enfa)
# co shows how important that variable is loading for the 1st and 2nd.  The li is the actual coordinate of the variable in multivariate space.  
gmgs.enfa$co

hist(gmgs.enfa, type="h")

scatter(gmgs.enfa)
scatter(gmgs.enfa, nc=T, percent=95, clabel=0, Acol="green", Ucol="blue",
        side="none", Aborder="black", Uborder="black")
scatterniche(gmgs.enfa$li, weights)
s.arrow(gmgs.enfa$co)
gmgs.enfa

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
     main="Callat - Sierra Nevada Core Landscape\nPCA")
ordihull(habavail.pca$li, groups=ua.fac.avail, draw="lines", col="green",
         lwd=2)
ordihull(habavail.pca$li, groups=ua.fac.used, show.groups="Used", col="blue",
         draw="lines", lwd=3)
legend(7,0, legend=c("Available","Used"), pch=21, pt.bg=c("green","blue"),
       cex=0.9, box.lty=0)



