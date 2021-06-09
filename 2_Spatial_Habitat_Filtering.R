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

# scroll through species and conduct spatial filtering
for (spec in coasts.spec){
  #read in zero-filled species data (spp-specific seasons)
  ebd.summ1 <- fread(file=paste0("",spec,"_zf.csv"))  
  ebd.summ <- as.data.frame(ebd.summ1) #convert to df
  
  ######################################
  # filter counts spatially using coast/wetland study area boundaries
  ebd_spdf <- SpatialPointsDataFrame(coords=ebd.summ[,c("longitude","latitude")], data=ebd.summ , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  if (spec %in% beach.spec){
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(coast_buff))
    ebd_hab.df <- as.data.frame(ebd_spdf.nad[coast_buff,])
  } else if (spec %in% marsh.spec){
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(marsh_buff))
    ebd_hab.df <- as.data.frame(ebd_spdf.nad[marsh_buff,] )
  } else { # species in both habitats
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(coast_buff))
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(marsh_buff))
    ebd_hab.df <- rbind(as.data.frame(ebd_spdf.nad[coast_buff,]), as.data.frame(ebd_spdf.nad[marsh_buff,]))
  }
  
  # remove old seasonal files
  rm(ebd.breed)
  rm(ebd.wint)
  ######################################
  # separate into breeding and wintering datasets using species-specific dates
  # seasonal date ranges from the eBird Status and Trends products (Fink et al., 2020)
  if (spec=="AMOY"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 129 &ebd_hab.df$Jdate < 209),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 326 | ebd_hab.df$Jdate < 54),]
  } else if (spec=="BLSK") {
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 150 &ebd_hab.df$Jdate < 216),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 333 | ebd_hab.df$Jdate < 82),]
  } else if (spec=="BRPE"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 129 &ebd_hab.df$Jdate < 180),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 347 | ebd_hab.df$Jdate < 61),]
  } else if (spec=="CLRA"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 143 &ebd_hab.df$Jdate < 258),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 340 | ebd_hab.df$Jdate < 89),]
  } else if (spec=="LBCU"){
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 312 | ebd_hab.df$Jdate < 54),]
  } else if (spec=="LETE"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 150 &ebd_hab.df$Jdate < 202),]
  } else if (spec=="MAGO"){
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 326 | ebd_hab.df$Jdate < 82),]
  } else if (spec=="PIPL"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 150 &ebd_hab.df$Jdate < 173),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 312 | ebd_hab.df$Jdate < 54),]
  } else if (spec=="REEG"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 59 &ebd_hab.df$Jdate < 228),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 334 | ebd_hab.df$Jdate < 47),]
  } else if (spec=="REKN"){
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 319 | ebd_hab.df$Jdate < 103),]
  } else if (spec=="SNPL"){
    ebd.breed <- ebd_hab.df[which(ebd_hab.df$Jdate > 143 &ebd_hab.df$Jdate < 180),]
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 298 | ebd_hab.df$Jdate < 75),]
  } else if (spec=="WESA"){
    ebd.wint <- ebd_hab.df[which(ebd_hab.df$Jdate > 326 | ebd_hab.df$Jdate < 61),]
  }
  
  ######################################
  # only keep seasonal points within appropriate seasonal distributional ranges
  # Seasonal ranges derived by intersecting eBird Status and Trends range maps:
  #     (https://ebird.org/science/status-and-trends/range-maps) with Bird Life 
  #     International range maps (http://datazone.birdlife.org/species/requestdis)
  # breeding
  if (exists("ebd.breed")){
    seas_buff <- readOGR(dsn="RangeMaps", layer=paste0(spec,"_Breeding")) 
    ebd_spdf <- SpatialPointsDataFrame(coords=ebd.breed[,c("longitude","latitude")], data=ebd.breed , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(seas_buff))
    ebd_breed.spdf <- ebd_spdf.nad[seas_buff,] 
    ebd.breed <- as.data.frame(ebd_breed.spdf)
  }
  #winter
  if (exists("ebd.wint")){
    seas_buff <- readOGR(dsn="RangeMaps", layer=paste0(spec,"_Winter"))
    ebd_spdf <- SpatialPointsDataFrame(coords=ebd.wint[,c("longitude","latitude")], data=ebd.wint , proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    ebd_spdf.nad <- spTransform(ebd_spdf, proj4string(seas_buff))
    ebd_wint.spdf <- ebd_spdf.nad[seas_buff,] 
    ebd.wint<- as.data.frame(ebd_wint.spdf)
  }
  
  
  #######################################
  ## identify stewardship site checklists
  
  # recombine breeding and wintering datasets, identified by season
  ebd.breed$Season <- "Breed"
  ebd.wint$Season <- "Winter"
  dat1 <- rbind(ebd.breed, ebd.wint)
  
  # Load in shapefile with stewardship site locations, included in GitHub repository
  sites_buff <- readOGR(dsn="Spatial", layer="StewardshipSitesBuffers")
  
  # convert to SPDF
  dat1 <- dat1[,-c(grep("latitude.",names(dat1)), grep("longitude.",names(dat1)))]
  dat1_spdf <- SpatialPointsDataFrame(coords=dat1[,c("longitude","latitude")], data=dat1, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  dat1_spdf.nad <- spTransform(dat1_spdf, crs(sites_buff))
  dat1_stew <- dat1_spdf.nad[sites_buff,] 
  dat1_stew.df <- as.data.frame(dat1_stew)
  dat1_stew.df$Site <- 1
  dat1$Site <- 0
  dat1_stew.df <- dat1_stew.df[,-c(grep("latitude.",names(dat1_stew.df)), grep("longitude.",names(dat1_stew.df)))] 
  dat1 <- rbind(dat1_stew.df, dat1[which(!(dat1$sampling_event_identifier %in% dat1_stew.df$sampling_event_identifier)),])
  
  #######################################
  ## identify protected area checklists
  
  # Load in protected area shapefiles. We derived habitat-specific protected area
  #   shapefiles by intersecting the habitat shapefiles used in 2_Spatial_Habitat_Filtering.R
  #   and protected areas from the Protected Areas Database 2.0 with GAP status 
  #   codes of 1-3: https://www.sciencebase.gov/catalog/item/5b030c7ae4b0da30c1c1d6de
  if (spec %in% beach.spec){
    pad <- readOGR(dsn="Spatial", layer="PAD_coastline")
  } else if (spec %in% marsh.spec){
    pad <- readOGR(dsn="Spatial", layer="PAD_wetlands")
  } else if (spec %in% both.spec){
    pad <- readOGR(dsn="Spatial", layer="PAD_coastline")
    pad.w <- readOGR(dsn="Spatial", layer="PAD_wetlands")
  }
  
  # convert to SPDF
  dat1_spdf <- SpatialPointsDataFrame(coords=dat1[,c("longitude","latitude")], data=dat1, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  dat1_spdf.nad <- spTransform(dat1_spdf, crs(pad))
  dat1_stew <- dat1_spdf.nad[pad,] 
  dat1_stew.df <- as.data.frame(dat1_stew)
  dat1_stew.df$SitePA <- 1
  if (spec %in% both.spec){ # add in wetland PAs for species in both.spec
    dat1_stew2 <- dat1_spdf.nad[pad.w,]
    dat1_stew.df2 <- as.data.frame(dat1_stew2)
    dat1_stew.df2$SitePA <- 1
    dat1_stew.df2 <- dat1_stew.df2[,-c(grep("latitude.",names(dat1_stew.df2)), grep("longitude.",names(dat1_stew.df2)))]
    dat1_stew.df <- rbind(dat1_stew.df, dat1_stew.df2[which(!(dat1_stew.df2$sampling_event_identifier %in% dat1_stew.df$sampling_event_identifier)),])
  }
  dat1$SitePA <- 0
  dat1_stew.df <- dat1_stew.df[,-c(grep("latitude.",names(dat1_stew.df)), grep("longitude.",names(dat1_stew.df)))] 
  dat1 <- rbind(dat1_stew.df, dat1[which(!(dat1$sampling_event_identifier %in% dat1_stew.df$sampling_event_identifier)),])
  
  ################################
  ## save filtered files
  if (spec %in% c("LBCU","MAGO","WESA","REKN")){ # WINTER ONLY
    fwrite(dat1, file=paste0(spec,"_zf_Winter_GulfAtlantic_FINAL.csv"))
  } else if (spec =="LETE") { # BREEDING ONLY
    fwrite(dat1, file=paste0(spec,"_zf_Breed_GulfAtlantic_FINAL.csv"))
  } else { # BREEDING AND WINTER
    fwrite(dat1, file=paste0(spec,"_zf_BreedWinter_GulfAtlantic_FINAL.csv")) 
  }
  
} # end loop through species


