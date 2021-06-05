################################################################################
## 2_Spatial_Habitat_Filtering.R
##
## This code groups species into habitat-based guilds and conducts species-
##      specific filtering based on habitat, breeding/wintering dates, and
##      distributional ranges. 
##
## R code accompanies the manuscript:
## Michel, N.L., S.P. Saunders, T.D. Meehan, C.B. Wilsey. 2021. Effects of
##      stewardship on protected area effectiveness for coastal birds. 
##      Conservation Biology, https://doi.org/10.1111/cobi.13698. 
################################################################################

Sys.setenv(TZ='GMT')

library(raster)
library(rgeos)
library(sf)
library(rgdal)
library(sp)
library(data.table)
options(scipen=999999)

# define species
#spec <- "AMOY"

# define all species
coasts.spec <- c("AMOY","BLSK","BRPE","CLRA","LBCU","LETE", "PIPL","REKN","REEG","SNPL","WESA","MAGO")

# group into habitat-based guilds to determine areas from which to draw observations and zeroes
beach.spec <- c("AMOY","BLSK","BRPE","LETE", "PIPL","REKN","REEG","SNPL","WESA","MAGO")  
marsh.spec <- c("CLRA")
both.spec <- c("LBCU")

# load coast buffer shapefile for cropping
# Load coastline buffered by 2 km. We used the GSHHG coastline from NOAA
#   available at: https://www.ngdc.noaa.gov/mgg/shorelines/
coast_buff <- readOGR(dsn="Coasts_StudyArea", layer="coastline_2km_GulfAtlantic")

# Load estuarine wetlands buffered by 1 km. We extracted estuarine wetlands 
#    (including estuarine forested, scrub/shrub, and emergent wetlands) from
#    NOAA's 2016 C-CAP Regional Land Cover dataset: 
#    https://coast.noaa.gov/digitalcoast/data/ccapregional.html
marsh_buff <- readOGR(dsn="Coasts_StudyArea", layer="wetlands_GulfAtlantic")

for (spec in coasts.spec){
  #read in zero-filled species data (spp-specific seasons)
  ebd.summ1 <- fread(file=paste0("",spec,"_zf.csv"))  
  ebd.summ <- as.data.frame(ebd.summ1) #convert to df
  
  #extract only counts within study area 
  ebd_spdf <- SpatialPointsDataFrame(coords=ebd.summ[,c("longitude","latitude")], data=ebd.summ , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  if (spec %in% beach.spec){
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(coast_buff))
    ebd_coast <- ebd_spdf.nad[coast_buff,] 
  } else if (spec %in% marsh.spec){
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(marsh_buff))
    ebd_coast <- ebd_spdf.nad[marsh_buff,] 
  } else {
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(coast_buff))
    ebd_coast <- ebd_spdf.nad[coast_buff,] 
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(marsh_buff))
    ebd_marsh <- ebd_spdf.nad[marsh_buff,]
    ebd_coast.df <- as.data.frame(ebd_coast)
    ebd_marsh.df <- as.data.frame(ebd_marsh)
    ebd_full.df <- rbind(ebd_coast.df, ebd_marsh.df)
    
  }
  ebd_coast.df <- as.data.frame(ebd_coast)
  
  if (spec %in% both.spec){
    fwrite(ebd_full.df, file=paste0(spec,"_zf_BreedWinter_GulfAtlantic.csv")) #write.csv
  } else {
    fwrite(ebd_coast.df, file=paste0(spec,"_zf_BreedWinter_GulfAtlantic.csv")) #write.csv
  }
}
rm(ebd_coast)
rm(ebd.summ)


## filter by seasonal dates on a species-specific basis
spec <- "MAGO"

# read in species files cropped to study area
ebd.gulfat <- fread(file=paste0(spec,"_zf_BreedWinter_GulfAtlantic.csv"))  
ebd.gulfatdf <- as.data.frame(ebd.gulfat)

# separate into breeding and wintering datasets using species-specific dates
# seasonal date ranges from the eBird Status and Trends products (Fink et al., 2020)
if (spec=="AMOY"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 129 &ebd.gulfatdf$Jdate < 209),]
} else if (spec=="BLSK") {
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 150 &ebd.gulfatdf$Jdate < 216),]
} else if (spec=="BRPE"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 129 &ebd.gulfatdf$Jdate < 180),]
} else if (spec=="CLRA"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 143 &ebd.gulfatdf$Jdate < 258),]
} else if (spec=="LETE"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 150 &ebd.gulfatdf$Jdate < 202),]
} else if (spec=="PIPL"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 150 &ebd.gulfatdf$Jdate < 173),]
} else if (spec=="REEG"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 59 &ebd.gulfatdf$Jdate < 228),]
} else if (spec=="SNPL"){
  ebd.breed <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 143 &ebd.gulfatdf$Jdate < 180),]
}

# winter
if (spec=="AMOY"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 326 | ebd.gulfatdf$Jdate < 54),]
} else if (spec=="BLSK") {
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 333 | ebd.gulfatdf$Jdate < 82),]
} else if (spec=="BRPE"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 347 | ebd.gulfatdf$Jdate < 61),]
} else if (spec=="CLRA"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 340 | ebd.gulfatdf$Jdate < 89),]
} else if (spec=="LBCU"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 312 | ebd.gulfatdf$Jdate < 54),]
} else if (spec=="MAGO"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 326 | ebd.gulfatdf$Jdate < 82),]
} else if (spec=="PIPL"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 312 | ebd.gulfatdf$Jdate < 54),]
} else if (spec=="REEG"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 334 | ebd.gulfatdf$Jdate < 47),]
} else if (spec=="REKN"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 319 | ebd.gulfatdf$Jdate < 103),]
} else if (spec=="SNPL"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 298 | ebd.gulfatdf$Jdate < 75),]
} else if (spec=="WESA"){
  ebd.wint <- ebd.gulfatdf[which(ebd.gulfatdf$Jdate > 326 | ebd.gulfatdf$Jdate < 61),]
}

# only keep seasonal points within appropriate seasonal distributional ranges
# Seasonal ranges derived by intersecting eBird Status and Trends range maps:
#     (https://ebird.org/science/status-and-trends/range-maps) with Bird Life 
#     International range maps (http://datazone.birdlife.org/species/requestdis)
seas_buff <- readOGR(dsn="RangeMaps", layer="MAGO_Breeding") 

ebd_spdf <- SpatialPointsDataFrame(coords=ebd.breed[,c("longitude","latitude")], data=ebd.breed , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(seas_buff))
ebd_breed <- ebd_spdf.nad[seas_buff,] 
ebd_breeddf <- as.data.frame(ebd_breed)

#winter
seas_buff <- readOGR(dsn="RangeMaps", layer="MAGO_Winter") 

ebd_spdf <- SpatialPointsDataFrame(coords=ebd.wint[,c("longitude","latitude")], data=ebd.wint , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(seas_buff))
ebd_winter <- ebd_spdf.nad[seas_buff,] 
ebd_winterdf <- as.data.frame(ebd_winter)

#merge breeding and wintering back together or save brdg/winter separately if only one season for that species
ebd_full.df <- rbind(ebd_breeddf, ebd_winterdf)

fwrite(ebd_full.df, file=paste0(spec,"_zf_BreedWinter_GulfAtlantic_FINAL.csv")) #write.csv
fwrite(ebd_winterdf, file=paste0(spec,"_zf_Winter_GulfAtlantic_FINAL.csv")) #write.csv
fwrite(ebd_breeddf, file=paste0(spec,"_zf_Breed_GulfAtlantic_FINAL.csv")) #write.csv

