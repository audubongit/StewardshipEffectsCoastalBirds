################################################################################
## 3_Trend_Modeling.R
##
## This code filters by stewardship and protected area status, conducts 
##      spatial case control sampling, then conducts seasonal population 
##      trend modeling in INLA. Before starting, manually download and
##      install the r-inla repository (not available on CRAN).
##      See instructions here: https://www.r-inla.org/download-install
##
## R code accompanies the manuscript:
## Michel, N.L., S.P. Saunders, T.D. Meehan, C.B. Wilsey. 2021. Effects of
##      stewardship on protected area effectiveness for coastal birds. 
##      Conservation Biology, https://doi.org/10.1111/cobi.13698. 
################################################################################


# set up ----
# options
setwd("~/StewardshipTrends/")
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

# libs
library(scales)
library(ggplot2)
library(stringr)
library(dplyr)
library(sf)
library(INLA)
library(brinla)

# inla posterior sampler function
inla.posterior.sampler <- function(result=inla1,
                                   full.ss=1000,
                                   part.ss=500,
                                   exclude_random="overdisp"){
  # set up
  out <- c()
  n.loops <- full.ss / part.ss
  all.par.names <- row.names(as.data.frame(
    inla.posterior.sample(1, result, num.threads = 3)[[1]]$latent))
  fix_samp_names <-row.names(result$summary.fixed)
  fitted_samp_names <- gsub(pattern="fitted.", replacement="",
                            row.names(result$summary.fitted.values))
  fitted_samp_names <- gsub("(?<![0-9])0+", "", fitted_samp_names, perl = TRUE)
  fitted_samp_names <- gsub(pattern=".", replacement=":", fitted_samp_names,
                            fixed=T)
  random_samp_names <- names(result$summary.random)[names(
    result$summary.random)!=exclude_random]
  # loop through samples
  for(i in 1:n.loops){
    post.1 <- inla.posterior.sample(part.ss, result, num.threads = 1)
    post.2 <- as.matrix(sapply(post.1, function(x) x$latent))
    block.1 <- t(post.2)
    colnames(block.1) <- all.par.names
    out <- rbind(out, block.1)
    rm(post.1)
    rm(post.2)
    gc()
    print(paste("Subsample", i, "of", n.loops))
  }
  # extract parts
  out_fix <- out[, grep(paste(fix_samp_names, collapse="|"), colnames(out))]
  colnames(out_fix) <- fix_samp_names
  out_fit <-  matrix(NA, 1, 1)
  #out[ , grep(paste(fitted_samp_names, collapse="|"), colnames(out))]
  random_samples <- list()
  for(i in 1:length(random_samp_names)){
    random_samples[[i]] <- out[, grep(random_samp_names[i], colnames(out))]
  }
  names(random_samples) <- random_samp_names
  # piece together output
  out_lst1 <- list(fixed_samples=out_fix,
                   random_samples=random_samples,
                   fitted_samples=out_fit)
  return(out_lst1)
}

# scatter plot style
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# pretty axes functions
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
nnz <- function(x) sum(as.numeric(x>0))
pnz <- function(x) sum(as.numeric(x>0)) / length(x)
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}
# ----

#################################
## MODEL BREEDING SEASON TRENDS
#################################

beach.spec <- c("AMOY","BLSK","BRPE","LETE", "PIPL","REKN","REEG","SNPL","WESA","MAGO")  
marsh.spec <- c("CLRA")
both.spec <- c("LBCU")

spec <- "PIPL"

#get data ---
if (spec %in% c("LBCU","MAGO","WESA","REKN")){
  dat1 <- read.csv(file=paste0(spec,"_zf_Winter_GulfAtlantic_FINAL.csv"))
} else if (spec %in% c("LETE","SOSH")) {
  dat1 <- read.csv(file=paste0(spec,"_zf_Breed_GulfAtlantic_FINAL.csv"))
} else {
  dat1 <- read.csv(file=paste0(spec,"_zf_BreedWinter_GulfAtlantic_FINAL.csv"))
}

# remove duplicates, if needed (shouldn't be necessary)
dat1 <- dat1[!(duplicated(dat1$sampling_event_identifier)),]

# Load in shapefile with stewardship site locations, included in GitHub repository
library(sp)
library(raster)
library(rgdal)
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
detach("package:raster", unload=TRUE)

# separate winter and breeding data depending on species
dat1.breed <- dat1[which(dat1$Jdate > 134 &dat1$Jdate < 182),]  

if (spec=="AMOY"){
  dat1.breed <- dat1[which(dat1$Jdate > 129 &dat1$Jdate < 209),]
} else if (spec=="BLSK") {
  dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 216),]
} else if (spec=="BRPE"){
  dat1.breed <- dat1[which(dat1$Jdate > 129 &dat1$Jdate < 180),]
} else if (spec=="CLRA"){
  dat1.breed <- dat1[which(dat1$Jdate > 143 &dat1$Jdate < 258),]
} else if (spec=="LETE"){
  dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 202),]
} else if (spec=="PIPL"){
  dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 173),]
} else if (spec=="REEG"){
  dat1.breed <- dat1[which(dat1$Jdate > 59 &dat1$Jdate < 228),]
} else if (spec=="SNPL"){
  dat1.breed <- dat1[which(dat1$Jdate > 143 &dat1$Jdate < 180),]
} 


#####################################################
# BREEDING TREND MODELING
#####################################################

dat1 <- dat1.breed 
# clean up data for modeling 
library(lubridate)
# clean up variables
dat1 <- dat1 %>% 
  mutate(
    # convert X observations to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# handle data and make derived vars
dat1$Strata <- paste(dat1$state_code, dat1$bcr_code, sep=".")
dat1 <- dat1 %>%
  select(count=observation_count, site=Site, sitePA=SitePA,year, duration=duration_minutes,
         day=Jdate, #time=time_dec,
         stratum=Strata, lon=longitude, lat=latitude,
         distance=effort_distance_km) %>%
  mutate(state=str_split(stratum, "[.]", simplify=T)[,1],
         bcr=str_split(stratum, "[.]", simplify=T)[,2],
         site=site-0.5,
         sitePA=sitePA-0.5,
         cyear=year-min(year)+1,
         lduration=log(duration+0.5),
         ldistance=log(distance+0.001),
         site_by_cyear=site*cyear,
         sitePA_by_cyear=sitePA*cyear
         ) %>%
  select(count, site, sitePA, year, cyear, site_by_cyear, sitePA_by_cyear, ldistance, lduration,
         state, bcr, lon, lat, stratum) %>%
  arrange(cyear, state, lon, lat) %>%
  as.data.frame()

#dat1 <- dat1 %>%
  #mutate(sitePA2=sitePA+-0.5,
         #sitePA2=sitePA2)

dat2 <- dat1

## remove NA counts (those checklists that only marked X instead of a number)
dat2 <- dat2 %>% filter(!is.na(count))

# make histograms to show distribution of counts
png(filename=paste0("Figs/",spec,"_CountDist.png"), width=3.4, height=3.4, units="in", res=300)
hist(dat2$count,breaks=50)
dev.off()

#-------------SPATIO TEMPORAL CASE CONTROL SAMPLING--------------------------------------------------------------------------------------
## Spatial case control sampling 
## From Fink et al. 2020 supp 2

# Draw one observation randomly from each 10 km grid cell, separately for site type & presences/absences within site types
# read in raster. We used a 10 km fishnet grid for North America
library(raster)
grid <- raster("grid_10km_land.tif")

# create spatial points object from lat and long fields
dat2_pts <- SpatialPoints(cbind(dat2$lon, dat2$lat), proj4string=CRS("+proj=longlat"))

# extract grid cell IDs at points
dat2_pts <- extract(grid, dat2_pts, method="simple", buffer=NULL, small=FALSE, cellnumbers=TRUE, fun=NULL, na.rm=TRUE, df=TRUE)

# add grid cell IDs back to data
dat2 <- cbind(dat2, Grid10km=dat2_pts$cells)

detach("package:raster", unload=TRUE)
library(dplyr)

## Randomly sample one observation from each grid cell for stewardship (-0.5) vs. non-stewardship (0.5) sites & PA vs. non-PA & absences vs. presences
# First group by year for subsequent subsetting
dat2a <- dat2 %>%
  group_by(year)

#then separate out by site type
dat2other <- subset(dat2a,site==-0.5)
dat2stew <- subset(dat2a,site==0.5)
dat2otherpa <- subset(dat2a,sitePA==-0.5)
dat2pa <- subset(dat2a,sitePA==0.5)

#then separate out presences and absences per site type
dat2other.p <- subset(dat2other,count>0)
dat2other.a <- subset(dat2other,count==0)
dat2stew.p <- subset(dat2stew,count>0)
dat2stew.a <- subset(dat2stew,count==0)

dat2otherpa.p <- subset(dat2otherpa,count>0)
dat2otherpa.a <- subset(dat2otherpa,count==0)
dat2pa.p <- subset(dat2pa,count>0)
dat2pa.a <- subset(dat2pa,count==0)

#group subsetted dataframes into list for function below
df.list <- list(dat2other.p,dat2other.a,dat2stew.p,dat2stew.a,dat2otherpa.p,dat2otherpa.a,dat2pa.p,dat2pa.a)

################################################
### ****START HERE WITH FUNCTION*****
parallel.f <- function(){
  #generate indices of random rows to keep within each of the subsets based on year & grid cell; repeat for each subsetted dataframe
  func <- function(df){
    df1 <- df[sample(1:dim(df)[1],replace=FALSE),] #randomly permutes data frame
    subset.new <- df1[!duplicated(df1[,c("year","Grid10km")]),] #removes duplicate observations per year/cell so only one remains
  }
  
  #perform above function for each dataframe and rbind the reduced subsets back together and call it dat3a
  dat3a <- do.call(rbind,lapply(df.list,func))
  
  ## Temporal component: balance per year sample size-------------------------------------------
  ssizes <- dat3a %>% 
    group_by(year) %>%
    summarise(no_rows = length(year))
  
  avg.ss <- mean(ssizes$no_rows)
  
  #years with less than average obs are sampled with replacement; years with more than average obs are sampled w/o replacement
  with.r <- which(ssizes$no_rows<avg.ss)
  without.r <- which(ssizes$no_rows>avg.ss)
  yrs.with <- with.r + 2006
  yrs.without <- without.r + 2006
  
  #sample each of those yrs.with with replacement to the average
  func2 <- function(yr){
    df <- dat3a %>% filter(year==yr)
    df2 <- df[sample(nrow(df),size=round(avg.ss),replace=TRUE),]
  }
  tt <- do.call(rbind,lapply(yrs.with,func2))
  
  #sample each of yrs.without without replacement to the average
  func3 <- function(yr){
    df <- dat3a %>% filter(year==yr)
    df2 <- df[sample(nrow(df),size=round(avg.ss),replace=FALSE),]
  }
  tt2 <- do.call(rbind,lapply(yrs.without,func3))
  
  dat3 <- rbind(tt,tt2)
  
  ##------------END OF SPATIO TEMPORAL CC SAMPLING-----------------------------------------------------------
  
  # add in random effect of state
  dat3 <- dat3 %>%
    mutate(state_int=as.integer(factor(state)),
           state_slp_st=as.integer(factor(state)),
           state_slp_yr=as.integer(factor(state)),
           state_slp_styr=as.integer(factor(state)),
           stratum_int=as.integer(factor(stratum))) %>%
    as.data.frame()
  
  # re-scale year to be centered so year effect is based on mid-point of time series
  dat3 <- dat3 %>%
    mutate(cyear_cent=scale(cyear,scale=FALSE),                       
           site_by_cyearcent=site*cyear_cent,
           sitePA_by_cyearcent=sitePA*cyear_cent) %>%
    select(count,site,sitePA,cyear_cent,site_by_cyearcent,sitePA_by_cyearcent,ldistance,lduration,state_int,lon,lat)
  
  # modeling ---------------------------------------------------------------------------------------------
  # make model formula
  form1 <- count ~ 1 + site + sitePA + cyear_cent + site_by_cyearcent + sitePA_by_cyearcent +
    lduration + ldistance +
    f(state_int, model="iid")
  
  famname <- "zeroinflatednbinomial1"  #use ZINFNBIN1 for all spp-seasons
  # run model  
  inla1 <- inla(form1, data=dat3, family=famname,
                verbose=F,
                control.compute=list(config=T), #cpo = T
                #control.inla=list(int.strategy="auto", strategy="adaptive"),
                num.threads = 3)
  
  # posterior sampling and plots -------------------------------------------------------
  # get posterior samples
  posterior_ss <- 3000
  post_samp <- inla.posterior.sampler(result=inla1,
                                      full.ss=posterior_ss,
                                      part.ss=1000)
  
  # get fixed intercept posterior
  fixed_samp <- as.data.frame(post_samp$fixed_samples)
  fix_int <- fixed_samp[,"(Intercept)"]
  
  # get fixed site posterior
  fix_site <- fixed_samp$site
  
  # get fixed site PA posterior
  fix_sitePA <- fixed_samp$sitePA
  
  # get fixed year posterior
  fix_year <- fixed_samp$cyear_cent
  
  # get fixed interaction posterior
  fix_cross <- fixed_samp$site_by_cyearcent
  
  # get fixed interaction posterior
  fix_crossPA <- fixed_samp$sitePA_by_cyearcent
  
  # fixed distance
  fix_dist <- fixed_samp$ldistance
  
  # fixed duration
  fix_dur <- fixed_samp$lduration
  
  # random state effect
  ran_state <- post_samp$random_samples$state_int
  
  # hyperparameter for zinf component
  hyperp <- rep(inla1$summary.hyperpar$mean[2],3000)
  
  out1 <-cbind(as.matrix(data.frame(fix_int=fix_int,fix_site=fix_site,fix_sitePA=fix_sitePA,
                                    fix_year=fix_year,fix_cross=fix_cross,fix_crossPA=fix_crossPA,fix_dist=fix_dist,fix_dur=fix_dur,hyperp=hyperp)),ran_state)
  colnames(out1)[10:ncol(out1)] <- sort(unique(dat2$state))
  rm(post_samp, inla1)
  gc()
  return(out1)
  
}
## END FUNCTION HERE #####
#test <- parallel.f() #make sure test is 3000 by ~20 (variable by number of state REs)

ptm <- proc.time()
sample.out <- replicate(100, parallel.f()) #should be 3000 by 20 by 100
proc.time() - ptm

runs <- 100 #specify how many runs were done above
post_all <- c()
for(i in 1:runs){
  dat = sample.out[,,i]
  post_all = rbind(post_all,dat)
} #should be 300k by ~20 (variable by number of REs)

#------------------------------------------------------------------------------------------------------
## Post-processing to determine significance; plotting of trends
#-----------------------------------------------------------------------------------------------------

#Subsample the full 300k posterior to 30k
subsamp <- sample_n(as.data.frame(post_all),0.1*nrow(post_all),replace=FALSE)

#Saving posteriors of parameter estimates (subsample of full 300k posterior)-----------------------------
#full summary (20 by 7)
post_summary <- data.frame(t(apply(subsamp, 2, function(x)
  quantile(x, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925), na.rm=T))))
names(post_summary) <- c("median", "lcl95", "ucl95","lcl50" ,"ucl50","lcl85","ucl85")

write.csv(post_summary, file=paste0(spec,"_Breed_ParamSummary.csv"))

## abundance estimates ---------------------------
#stewardship avg abundance (non protected)
StewAvgAb <- quantile((1-mean(subsamp$hyperp))*exp(subsamp$fix_int + subsamp$fix_site*0.5 + subsamp$fix_sitePA*-0.5), probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

#protected area avg abundance (non stewardship)
PAAvgAb <- quantile((1-mean(subsamp$hyperp))*exp(subsamp$fix_int + subsamp$fix_site*-0.5 + subsamp$fix_sitePA*0.5), probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

#non stewardship non protected abundance
NeitherAvgAb <- quantile((1-mean(subsamp$hyperp))*exp(subsamp$fix_int + subsamp$fix_site*-0.5 + subsamp$fix_sitePA*-0.5), probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

# stewardship and protected abundance
BothAvgAb <- quantile((1-mean(subsamp$hyperp))*exp(subsamp$fix_int + subsamp$fix_site*0.5 + subsamp$fix_sitePA*0.5), probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

abun_est <- rbind(StewAvgAb,PAAvgAb,NeitherAvgAb,BothAvgAb)
write.csv(abun_est, file=paste0(spec,"_Breed_AbundanceSummary.csv"))

## trend estimates (% per year) -------------------------------
# trend stewardship (not protected)
StewTrend <- quantile((1-mean(subsamp$hyperp))*(exp(subsamp$fix_year + subsamp$fix_cross*0.5 + subsamp$fix_crossPA*-0.5)-1)*100, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

# protected trend (not stewardship)
PATrend <- quantile((1-mean(subsamp$hyperp))*(exp(subsamp$fix_year + subsamp$fix_cross*-0.5 + subsamp$fix_crossPA*0.5)-1)*100, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

# non stewardship non protected trend
NeitherTrend <- quantile((1-mean(subsamp$hyperp))*(exp(subsamp$fix_year + subsamp$fix_cross*-0.5 + subsamp$fix_crossPA*-0.5)-1)*100, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

# stewardship and protected trend
BothTrend <- quantile((1-mean(subsamp$hyperp))*(exp(subsamp$fix_year + subsamp$fix_cross*0.5 + subsamp$fix_crossPA*0.5)-1)*100, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))

trend_est <- rbind(StewTrend,PATrend,NeitherTrend,BothTrend)
write.csv(trend_est, file=paste0(spec,"_Breed_TrendSummary.csv"))

# make prediction matrices for plotting-------------------------
posterior_ss=nrow(subsamp)
vyear <- seq(-5.5,5.5,by=1) #centered year
myear <- matrix(vyear, ncol=length(vyear), nrow=posterior_ss, byrow=T)
mint <- matrix(subsamp$fix_int, ncol=length(vyear), nrow=posterior_ss)
msteff <- matrix(subsamp$fix_site, ncol=length(vyear), nrow=posterior_ss)
msteffPA <- matrix(subsamp$fix_sitePA, ncol=length(vyear), nrow=posterior_ss)
myreff <- matrix(subsamp$fix_year, ncol=length(vyear), nrow=posterior_ss)
mcrosseff <- matrix(subsamp$fix_cross, ncol=length(vyear), nrow=posterior_ss)
mcrosseffPA <- matrix(subsamp$fix_crossPA, ncol=length(vyear), nrow=posterior_ss)
mst0 <- matrix(-0.5, ncol=length(vyear), nrow=posterior_ss, byrow=T)
mst1 <- matrix(0.5, ncol=length(vyear), nrow=posterior_ss, byrow=T)

# calculate predictions and make df
pred00 <- (1-mean(subsamp$hyperp))*exp(mint + msteff*mst0 + msteffPA*mst0 + myreff*myear + mcrosseff*mst0*myear + mcrosseffPA*mst0*myear) #non non
pred01 <- (1-mean(subsamp$hyperp))*exp(mint + msteff*mst0 + msteffPA*mst1 + myreff*myear + mcrosseff*mst0*myear + mcrosseffPA*mst1*myear) #protected
pred10 <- (1-mean(subsamp$hyperp))*exp(mint + msteff*mst1 + msteffPA*mst0 + myreff*myear + mcrosseff*mst1*myear + mcrosseffPA*mst0*myear) #stewardship
pred11 <- (1-mean(subsamp$hyperp))*exp(mint + msteff*mst1 + msteffPA*mst1 + myreff*myear + mcrosseff*mst1*myear + mcrosseffPA*mst1*myear) #both

neither <- data.frame(year=rep(1:12,1), t(apply(pred00, 2, function(x)
  quantile(x, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925), na.rm=T))), Cross=rep("Neither",12))

patrend <- data.frame(year=rep(1:12,1), t(apply(pred01, 2, function(x)
  quantile(x, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925), na.rm=T))), Cross=rep("PA",12))

stewtrend <- data.frame(year=rep(1:12,1), t(apply(pred10, 2, function(x)
  quantile(x, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925), na.rm=T))), Cross=rep("Stew",12))

both <- data.frame(year=rep(1:12,1), t(apply(pred11, 2, function(x)
  quantile(x, probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925), na.rm=T))), Cross=rep("Both",12))

pred_all <- rbind(neither,patrend,stewtrend,both)
pred_all$year <- pred_all$year + 2006

write.csv(pred_all, file=paste0(spec,"_Breed_PredictionSummary.csv"))


## figs -------------------------------------------------------------------------------
# I. full plot of all four trends with 50% CIs
plot.new()
p1 <- ggplot(data=pred_all) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), #85 CIs
  alpha=0.1) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("Neither", "PA", "Stew","Both"), 
                     labels=c("Neither", "Protected", "Stewardship", "Both"), values=c("#66A61E","#E6AB02","#7570B3","#b370af")) + #colors need to be in alph order
  scale_fill_manual(name="", breaks=c("Neither", "PA", "Stew","Both"), 
                    labels=c("Neither", "Protected", "Stewardship", "Both"), values=c("#66A61E","#E6AB02","#7570B3","#b370af")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(),
                     labels = prettyNum) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1

#save final fig w/ 50 and 85% CIs
png(filename=paste0("Figs/",spec,"_4Trends_50CIs.png"), width=3.6, height=3.4, units="in", res=600)
p1
dev.off()

# II. focal plot of Stew vs. Protected trends with 50% and 85% CIs
plot.new()
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X25., ymax=X75., fill=Cross), #50 CIs
              alpha=0.25) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), #85 CIs
              alpha=0.15) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#66A61E","#b370af")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#66A61E","#b370af")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(),
                     labels = prettyNum) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2

#save final fig w/ 50 and 85% CIs
png(filename=paste0("Figs/",spec,"_StewvsPA_50_85CIs.png"), width=3.6, height=3.4, units="in", res=600)
p2
dev.off()

###############################################
## Goodness of fit evaluation
################################################
#First generate a new dataset for testing

new.data <- function (){
  func <- function(df){
    df1 <- df[sample(1:dim(df)[1],replace=FALSE),] #randomly permutes data frame
    subset.new <- df1[!duplicated(df1[,c("year","Grid10km")]),] #removes duplicate observations per year/cell so only one remains
  }
  
  #perform above function for each dataframe and rbind the reduced subsets back together and call it dat3a
  dat3a <- do.call(rbind,lapply(df.list,func))
  
  ## Temporal component: balance per year sample size-------------------------------------------
  ssizes <- dat3a %>% 
    group_by(year) %>%
    summarise(no_rows = length(year))
  
  avg.ss <- mean(ssizes$no_rows)
  
  #years with less than average obs are sampled with replacement; years with more than average obs are sampled w/o replacement
  with.r <- which(ssizes$no_rows<avg.ss)
  without.r <- which(ssizes$no_rows>avg.ss)
  yrs.with <- with.r + 2006
  yrs.without <- without.r + 2006
  
  #sample each of those yrs.with with replacement to the average
  func2 <- function(yr){
    df <- dat3a %>% filter(year==yr)
    df2 <- df[sample(nrow(df),size=round(avg.ss),replace=TRUE),]
  }
  tt <- do.call(rbind,lapply(yrs.with,func2))
  
  #sample each of yrs.without without replacement to the average
  func3 <- function(yr){
    df <- dat3a %>% filter(year==yr)
    df2 <- df[sample(nrow(df),size=round(avg.ss),replace=FALSE),]
  }
  tt2 <- do.call(rbind,lapply(yrs.without,func3))
  
  dat3 <- rbind(tt,tt2)
  
  dat3 <- dat3 %>%
    mutate(state_int=as.integer(factor(state)),
           state_slp_st=as.integer(factor(state)),
           state_slp_yr=as.integer(factor(state)),
           state_slp_styr=as.integer(factor(state)),
           stratum_int=as.integer(factor(stratum))) %>%
    as.data.frame()
  
  # re-scale year to be centered so year effect is based on mid-point of time series
  dat3 <- dat3 %>%
    mutate(cyear_cent=scale(cyear,scale=FALSE),                       
           site_by_cyearcent=site*cyear_cent,
           sitePA_by_cyearcent=sitePA*cyear_cent) %>%
    select(count,site,sitePA,cyear_cent,site_by_cyearcent,sitePA_by_cyearcent,ldistance,lduration,state_int,lon,lat)
  
  return(dat3)
}

#dat4 <- new.data()

#next use param estimates (from ensemble run of 100) to fit new data & compare obs vs. fitted
library(fastDummies)

fit <- c()
for (i in 1:500){ #do this 500 times to yield a posterior of goodness of fit stats
  dat4 <- new.data()
  dat4 <- dummy_cols(dat4, select_columns = "state_int")
  par_sub <- sample_n(subsamp,1)
  #par_sub <- as.data.frame(matrix(apply(subsamp,2,median),nrow=1))
  #colnames(par_sub) <- c("fix_int","fix_site","fix_sitePA","fix_year","fix_cross","fix_crossPA","fix_dist","fix_dur")
  fitted.new <- (1-mean(subsamp$hyperp))* exp(par_sub$fix_int + par_sub$fix_site*dat4$site + par_sub$fix_sitePA*dat4$sitePA + par_sub$fix_year*dat4$cyear_cent + 
                                                par_sub$fix_cross*dat4$site_by_cyearcent + par_sub$fix_crossPA*dat4$sitePA_by_cyearcent +
                                                par_sub$fix_dist*dat4$ldistance + par_sub$fix_dur*dat4$lduration +
                                                par_sub[,10]*dat4$state_int_1 + par_sub[,11]*dat4$state_int_2 + par_sub[,12]*dat4$state_int_3 + par_sub[,13]*dat4$state_int_4 + 
                                                par_sub[,14]*dat4$state_int_5 + par_sub[,15]*dat4$state_int_6 + par_sub[,16]*dat4$state_int_7 + par_sub[,17]*dat4$state_int_8 + 
                                                par_sub[,18]*dat4$state_int_9 + par_sub[,19]*dat4$state_int_10 + par_sub[,20]*dat4$state_int_11 #+ 
                                              #par_sub[,21]*dat4$state_int_12 + par_sub[,22]*dat4$state_int_13 + par_sub[,23]*dat4$state_int_14 + par_sub[,24]*dat4$state_int_15
  ) #NOTE: need to manually change this equation depending on number of state REs for species/season
  
  fit.q <- quantile(fitted.new,probs=c(0,1))
  obs.q <- quantile(dat4$count,probs=c(0,1))
  fit.dif <- fit.q[2]-fit.q[1]
  obs.dif <- obs.q[2]-obs.q[1]
  test.q <- fit.q[2]-obs.q[2] #difference in spread of fitted vs. obs values
  mean.dif <- mean(fitted.new)-mean(dat4$count)
  disp.stat <- sum(((dat4$count)-fitted.new)^2)/(nrow(dat4)-ncol(par_sub)) #dispersion stat
  sse <- sum((dat4$count - fitted.new)^2) #sum of squares error
  spcor <- cor(fitted.new,dat4$count,method="spearman") #spearman corr
  chisq <- sum((dat4$count-fitted.new)^2/fitted.new) #chi-squared
  freeTuke <- sum((sqrt(dat4$count) - sqrt(fitted.new))^2) #freemans Tukey
  df = data.frame(spcor,chisq,sse,freeTuke,test.q,fit.dif,obs.dif,mean.dif,disp.stat)
  fit = rbind(fit,df)
} #run time approx 5 mins

#dim(fit)

#hist(fit$spcor)
#hist(fit$chisq)
#hist(fit$sse)
#hist(fit$freeTuke)
#hist(fit$test.q)

sc <- quantile(fit$spcor,probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
chis <- quantile(fit$chisq,probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
ssq <- quantile(fit$sse,probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
ftuk <- quantile(fit$freeTuke,probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
testq <- quantile(fit$test.q,probs=c(0.5, 0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
fdifs <- quantile(fit$fit.dif,probs=c(0.5,0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
obsdifs <- quantile(fit$obs.dif,probs=c(0.5,0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
meandifs <- quantile(fit$mean.dif,probs=c(0.5,0.025, 0.975, 0.25, 0.75, 0.075, 0.925))
dispstat <- quantile(fit$disp.stat,probs=c(0.5,0.025, 0.975, 0.25, 0.75, 0.075, 0.925)) #should be relatively close to 1
#pseudo-chat to measure overdispersion
chat<-fit$chisq/mean(fit$chisq)
chats <- quantile(chat,probs=c(0.5,0.025, 0.975, 0.25, 0.75, 0.075, 0.925)) #should be relatively close to 1

fitstats <- rbind(sc,chis,ssq,ftuk,testq,chats,fdifs,obsdifs,meandifs,dispstat) #most informative metrics: sc, chats, fdifs vs. obsdifs, meandifs
write.csv(fitstats, file=paste0(spec,"_Breed_GOFstats.csv"))


## End ##-------------------------


