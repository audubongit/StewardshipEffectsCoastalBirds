################################################################################
## 4_Plotting_Results.R
##
## This code summarizes the results of abundance and trend models
##      and produces tables and figures included in publication and supplement. 
##
## R code accompanies the manuscript:
## Michel, N.L., S.P. Saunders, T.D. Meehan, C.B. Wilsey. 2021. Effects of
##      stewardship on protected area effectiveness for coastal birds. 
##      Conservation Biology, https://doi.org/10.1111/cobi.13698. 
######################################################################################################################################

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(cowplot)
library(scales)
library(rgdal)
library(sp)
library(raster)

setwd("~/StewardshipTrends/")
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

pathtofigs <- "Figs/"


##########################################################################
##  TREND - SUMMARIZE AND PLOT FIGURE 2
##########################################################################

## summarize trend and trend differences
setwd(pathtooutput)
trenddiff <- data.frame()
trendall <- data.frame()
lf <- list.files(pattern = "*_TrendSummary.csv")
for (i in 1:length(lf)){
  temp <- read.csv(lf[i])
  trenddiff <- rbind(trenddiff, data.frame(Species=substr(lf[i],1,4), Season=substr(lf[i],6,11), StewDiff=(temp$X50.[which(temp$X=="StewTrend")]-temp$X50.[which(temp$X=="NeitherTrend")])/temp$X50.[which(temp$X=="NeitherTrend")],
                                           PADiff=(temp$X50.[which(temp$X=="PATrend")]-temp$X50.[which(temp$X=="NeitherTrend")])/temp$X50.[which(temp$X=="NeitherTrend")]))
  temp$Species <- substr(lf[i],1,4)
  temp$Season <- substr(lf[i],6,11)
  trendall <- rbind(trendall, temp)
}
trenddiff$Season <- as.character(trenddiff$Season)
trenddiff$Season[which(trenddiff$Season=="Breed_")] <- "Breed"
write.csv(trenddiff, file=paste0(pathtooutput, "AllSpecies_TrendDifferences.csv"))
trendall$Season <- as.character(trendall$Season)
trendall$Season[which(trendall$Season=="Breed_")] <- "Breed"
write.csv(trendall, file=paste0(pathtooutput, "AllSpecies_AllTrends.csv"))

mean(trenddiff$StewDiff) # 40.98492
sd(trenddiff$StewDiff)/sqrt(nrow(trenddiff)) # se: 43.84515
# 85% ci: multiply se by 1.44  (-22.1521, 104.1219)

mean(trenddiff$PADiff) # 1.472768
sd(trenddiff$PADiff)/sqrt(nrow(trenddiff)) # se: 3.928569
# 85% ci: multiply se by 1.44  (-4.184371, 7.129907)


######################
### plot breeding trend
library(ggplot2)
trendbreed <- trendall[which(trendall$X != "BothTrend" & trendall$Season=="Breed"),]
names(trendbreed)[which(names(trendbreed)=="X")] <- "Site"
trendbreed <- trendbreed[order(trendbreed$Species, decreasing=FALSE),]
trendbreed$Species <- as.factor(trendbreed$Species)
trendbreed$siteord <- 1
trendbreed$siteord[which(trendbreed$Site=="NeitherTrend")] <- 2
trendbreed$siteord[which(trendbreed$Site=="PATrend")] <- 3
trendbreed <- trendbreed[order(trendbreed$Species, trendbreed$siteord),]
pbrcols <- c("#1b9e77", "#d95f02", "#7570b3")
# make vertical axis (x) numeric
trendbreed$SpecNum <- rev(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3)))
png(file=paste0(pathtofigs,"Trend_All_Breed.png"), width=6.4, height=4.74, units="in" , res=600)
p <- ggplot(trendbreed12, aes(fill=rev(Site), y=X50., x=SpecNum)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  geom_vline(xintercept= 1.5, col="#CAD4D1")+
  geom_vline(xintercept= 2.5, col="#CAD4D1")+
  geom_vline(xintercept= 3.5, col="#CAD4D1")+
  geom_vline(xintercept= 4.5, col="#CAD4D1")+
  geom_vline(xintercept= 5.5, col="#CAD4D1")+
  geom_vline(xintercept= 6.5, col="#CAD4D1")+
  geom_vline(xintercept= 7.5, col="#CAD4D1")+
  scale_color_manual(name="", breaks=c("StewTrend", "NeitherTrend", "PATrend"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="", breaks=c("StewTrend", "NeitherTrend", "PATrend"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_x_continuous(limits=c(0.7,8.3), breaks=c(1:8),
                     labels=rev(c("American Oystercatcher","Black Skimmer","Brown Pelican", "Clapper Rail", "Least Tern", "Piping Plover", "Reddish Egret", "Snowy Plover")))+
  coord_flip()+
  labs(x="", y="Trend (median + 85% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = c(0.85,0.95), legend.background=element_blank(), 
        legend.text=element_text(size=10), legend.key.size = unit(0.45, 'cm'))
print(p)
dev.off()

#### produce donut plots
# Create new dataset for each donut - manually enter the number of species meeting each criteria
# shown for STW>UP
data <- data.frame( category=c("YES", "NO"),count=c(4,4))
# Compute percentages
data$fraction <- data$count / sum(data$count)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$category, "\n (n = ", data$count, ")")

# Make the plot - create each plot manually and name brd1 - brd4
# green (stw>up): "#1b9e77"; purple (pr>up): "#7570b3"; blue (stw>pr): "#629FD7"; gray (null): "#CAD4D1" 
brd1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=0.9, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer) #blue: #629FD7
  scale_fill_manual(values=c("#CAD4D1","#1b9e77")) + # first color is grey, change second color value
  scale_color_manual(values=c("#000000","#000000")) + # text color = black
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

# save strip of donut figures
png(file=paste0(pathtofigs,"Donut_Trend_Breed.png"), width=11, height=3.74, units="in" , res=600)
grid.arrange(brd1,brd2,brd3,brd4, ncol=4)
dev.off()

#make labels - print once for all
par(mar = c(0,0,0,0))
png(file=paste0(pathtofigs,"Donut_Label4.png"), width=4, height=4, units="in" , res=300)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("PR > STW"), 
     cex = 1.1, col = "black")
dev.off()
#restore settings
par(mar = c(5, 4, 4, 2) + 0.1)


#######################
#### plot winter trend 
trendwint <- trendall[which(trendall$X != "BothTrend" & trendall$Season=="Winter"),]
names(trendwint)[which(names(trendwint)=="X")] <- "Site"
trendwint$Species <- as.character(trendwint$Species)
trendwint$Species[which(trendwint$Species=="REKN")] <- "REAN"
trendwint <- trendwint[order(trendwint$Species, decreasing=FALSE),]
trendwint$Species <- as.factor(trendwint$Species)
trendwint$siteord <- 1
trendwint$siteord[which(trendwint$Site=="NeitherTrend")] <- 2
trendwint$siteord[which(trendwint$Site=="PATrend")] <- 3
trendwint <- trendwint[order(trendwint$Species, trendwint$siteord),]
pbrcols <- c("#1b9e77", "#d95f02", "#7570b3")
trendwint$SpecNum <- rev(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3),
                           rep(9,3), rep(10,3), rep(11,3)))
png(file=paste0(pathtofigs,"Trend_All_Winter.png"), width=6.4, height=4.74, units="in" , res=600)
p <- ggplot(trendwint, aes(fill=rev(Site), y=X50., x=SpecNum)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  geom_vline(xintercept= 1.5, col="#CAD4D1")+
  geom_vline(xintercept= 2.5, col="#CAD4D1")+
  geom_vline(xintercept= 3.5, col="#CAD4D1")+
  geom_vline(xintercept= 4.5, col="#CAD4D1")+
  geom_vline(xintercept= 5.5, col="#CAD4D1")+
  geom_vline(xintercept= 6.5, col="#CAD4D1")+
  geom_vline(xintercept= 7.5, col="#CAD4D1")+
  geom_vline(xintercept= 8.5, col="#CAD4D1")+
  geom_vline(xintercept= 9.5, col="#CAD4D1")+
  geom_vline(xintercept= 10.5, col="#CAD4D1")+
  scale_color_manual(name="", breaks=c("StewTrend", "NeitherTrend", "PATrend"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="", breaks=c("StewTrend", "NeitherTrend", "PATrend"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_x_continuous(limits = c(0.7,11.3),breaks=c(1:11),
                     labels=rev(c("American Oystercatcher","Black Skimmer","Brown Pelican", "Clapper Rail", "Long-billed Curlew", "Marbled Godwit", "Piping Plover", "Red Knot", "Reddish Egret", "Snowy Plover", "Western Sandpiper")))+
  coord_flip()+
  labs(x="", y="Trend (median + 85% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = c(0.85,0.95), legend.background=element_blank(), 
        legend.text=element_text(size=10), legend.key.size = unit(0.45, 'cm'))
print(p)
dev.off()

#### produce donut plots
# Create new dataset for each donut - manually enter the number of species meeting each criteria
# shown for STW>UP
data <- data.frame( category=c("YES", "NO"),count=c(2,9))
# Compute percentages
data$fraction <- data$count / sum(data$count)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$category, "\n (n = ", data$count, ")")

# Make the plot - create each plot manually and name brd1 - brd4
# green (stw>up): "#1b9e77"; purple (pr>up): "#7570b3"; blue (stw>pr): "#629FD7"; gray (null): "#CAD4D1" 
win1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=0.9, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer) #blue: #629FD7
  scale_fill_manual(values=c("#CAD4D1","#1b9e77")) + # first color is grey, change second color value
  scale_color_manual(values=c("#000000","#000000")) + # text color = black
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

# save strip of donut figures
png(file=paste0(pathtofigs,"Donut_Trend_Winter.png"), width=11, height=3.74, units="in" , res=600)
grid.arrange(win1,win2,win3,win4, ncol=4)
dev.off()



##########################################################################
##  PRODUCE TREND PLOTS FOR SPECIES WHERE STEWARDSHIP > PROTECTED AREA
##    FOR FIGURE 3
##########################################################################

### PIPL Winter
pred_all <- read.csv(file="PIPL_Winter_PredictionSummary.csv")
png(file=paste0(pathtofigs,"PBR_PIPLWinTrend.png"), width=3.6, height=3.4, units="in", res=600)
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), 
              alpha=0.3) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=c(0.001, 0.005, 0.050, 0.500, 1.000),
                     labels = c("0.001", "0.005", "0.050", "0.500", "1.000")) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

### PIPL Summer
pred_all <- read.csv(file="PIPL_Breed_PredictionSummary.csv")
png(file=paste0(pathtofigs,"PBR_PIPLBreedTrend.png"), width=3.6, height=3.4, units="in", res=600)
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross),
              alpha=0.3) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=c(0.001, 0.005, 0.050, 0.500, 1.000, 5.000),
                     labels = c("0.001", "0.005", "0.050", "0.500", "1.000", "5.000")) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

### BLSK Summer
pred_all <- read.csv(file="E:/PriorityBirds_Report/Output_AudubonProAreas/BLSK_Breed_PredictionSummary.csv")
png(file=paste0(pathtofigs,"PBR_BLSKBreedTrend.png"), width=3.6, height=3.4, units="in", res=600)
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), 
              alpha=0.3) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=c(0.005, 0.020, 0.050, 0.200, 0.500, 1.000),
                     labels = c("0.005", "0.020", "0.050", "0.200", "0.500", "1.000")) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

### BRPE Summer
pred_all <- read.csv(file="BRPE_Breed_PredictionSummary.csv")
png(file=paste0(pathtofigs,"PBR_BRPEBreedTrend.png"), width=3.6, height=3.4, units="in", res=600)
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), #85 CIs
              alpha=0.3) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=c(0.5, 1.0, 2.0, 5.0, 10.0),
                     labels = c("0.500", "1.000", "2.000", "5.000", "10.000")) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

### LETE Summer
pred_all <- read.csv(file="LETE_Breed_PredictionSummary.csv")
png(file=paste0(pathtofigs,"PBR_LETEBreedTrend.png"), width=3.6, height=3.4, units="in", res=600)
p2 <- ggplot(data=pred_all[13:36,]) +
  geom_ribbon(aes(x=year, ymin=X7.5., ymax=X92.5., fill=Cross), #85 CIs
              alpha=0.3) +
  geom_line(aes(x=year, y=X50., group=Cross, col=Cross),lwd=1.5) +
  scale_color_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  scale_fill_manual(name="", breaks=c("PA","Stew"), labels=c("Protected", "Stewardship"), values=c("#7570B3","#1B9E77")) +
  xlab("Year") + ylab("Abundance index") +
  scale_y_continuous(trans=log_trans(), breaks=c(0.05, 0.1, 0.2, 0.5, 1.0),
                     labels = c("0.050", "0.100", "0.200", "0.500", "1.000")) +
  scale_x_continuous(breaks=c(2008, 2010, 2012, 2014, 2016, 2018), labels = c(2008, 2010, 2012, 2014, 2016, 2018)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 14),
        axis.title.x=element_text(size = 14))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()



##########################################################################
##  SUMMARIZE MODEL PARAMETERS FOR TABLE 2 AND APPENDIX S4
##########################################################################

speclist <- c("AMOY","BLSK","BRPE","LETE", "PIPL","REKN","REEG","SNPL","WESA","MAGO", "CLRA", "LBCU")  

alldat <- data.frame()
for (spec in speclist){
  #get data ---
  if (spec %in% c("LBCU","MAGO","WESA","REKN")){
    dat1 <- read.csv(file=paste0(spec,"_zf_Winter_GulfAtlantic_FINAL.csv"))
    dat1$Spec <- spec
  } else if (spec %in% c("LETE","SOSH")) {
    dat1 <- read.csv(file=paste0(spec,"_zf_Breed_GulfAtlantic_FINAL.csv"))
    dat1$Spec <- spec
  } else {
    dat1 <- read.csv(file=paste0(spec,"_zf_BreedWinter_GulfAtlantic_FINAL.csv"))
    dat1$Spec <- spec
  }
  
  # remove duplicates, if needed (shouldn't be necessary)
  dat1 <- dat1[!(duplicated(dat1$sampling_event_identifier)),]
  alldat <- rbind(alldat, dat1)
}


dat1 <- alldat

all.breed <- data.frame()
all.wint <- data.frame()
for (spec in speclist){
  if (spec=="AMOY"){
    dat1.breed <- dat1[which(dat1$Jdate > 129 &dat1$Jdate < 209),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 129 &dat1$Jdate < 209)),]
  } else if (spec=="BLSK") {
    dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 216),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 150 &dat1$Jdate < 216)),]
  } else if (spec=="BRPE"){
    dat1.breed <- dat1[which(dat1$Jdate > 129 &dat1$Jdate < 180),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 129 &dat1$Jdate < 180)),]
  } else if (spec=="CLRA"){
    dat1.breed <- dat1[which(dat1$Jdate > 143 &dat1$Jdate < 258),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 143 &dat1$Jdate < 258)),]
  } else if (spec=="LETE"){
    dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 202),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 150 &dat1$Jdate < 202)),]
  } else if (spec=="PIPL"){
    dat1.breed <- dat1[which(dat1$Jdate > 150 &dat1$Jdate < 173),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 150 &dat1$Jdate < 173)),]
  } else if (spec=="REEG"){
    dat1.breed <- dat1[which(dat1$Jdate > 59 &dat1$Jdate < 228),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 59 &dat1$Jdate < 228)),]
  } else if (spec=="SNPL"){
    dat1.breed <- dat1[which(dat1$Jdate > 143 &dat1$Jdate < 180),]
    dat1.wint <- dat1[which(!(dat1$Jdate > 143 &dat1$Jdate < 180)),]
  }
  
  all.breed <- rbind(all.breed, dat1.breed)
  all.wint <- rbind(all.wint, dat1.wint)
}

# remove duplicates and produce summaries

spec <- "AMOY"
specwint <- all.wint[which(all.wint$Spec==spec),]
specwint <- specwint[!(duplicated(specwint$sampling_event_identifier)),]
specwint <- specwint[which(specwint$observation_count != "X"),]
nrow(specwint[which(specwint$observation_count!="0"),])
nrow(specwint)

specbr <- all.breed[which(all.breed$Spec==spec),]
specbr <- specbr[!(duplicated(specbr$sampling_event_identifier)),]
nrow(specbr[which(specbr$observation_count!="0"),])
nrow(specbr)

# check against old table
if (spec %in% c("LBCU","MAGO","WESA","REKN")){
  nrow(specwint) 
} else if (spec == "LETE") {
  nrow(specbr)
} else {
  nrow(specwint) + nrow(specbr)
}

# compile parameter summaries
param.breed <- data.frame()
param.wint <- data.frame()
for (spec in speclist){
  #get data ---
  if (spec %in% c("LBCU","MAGO","WESA","REKN")){
    pw <- read.csv(file=paste0(spec,"_Winter_ParamSummary.csv"))
    pw$Spec <- spec
    param.wint <- rbind(param.wint, pw[,c("Spec","X","median","lcl85","ucl85")])
  } else if (spec == "LETE") {
    pb <- read.csv(file=paste0(spec,"_Breed_ParamSummary.csv"))
    pb$Spec <- spec
    param.breed <- rbind(param.breed, pb[,c("Spec","X","median","lcl85","ucl85")])
  } else {
    pw <- read.csv(file=paste0(spec,"_Winter_ParamSummary.csv"))
    pw$Spec <- spec
    param.wint <- rbind(param.wint, pw[,c("Spec","X","median","lcl85","ucl85")])
    pb <- read.csv(file=paste0(spec,"_Breed_ParamSummary.csv"))
    pb$Spec <- spec
    param.breed <- rbind(param.breed, pb[,c("Spec","X","median","lcl85","ucl85")])
  }
}
# Breeding
param.breed$BreedParam <- paste0(round(param.breed$median, 2), " (", round(param.breed$lcl85, 2), ", ", round(param.breed$ucl85,2), ")")
param.breed$ParamName <- "Intercept"
param.breed$ParamName[which(param.breed$X == "fix_site")] <- "Stewardship site type"
param.breed$ParamName[which(param.breed$X == "fix_sitePA")] <- "Protected site type"
param.breed$ParamName[which(param.breed$X == "fix_year")] <- "Year"
param.breed$ParamName[which(param.breed$X == "fix_cross")] <- "Stewardship site x year"
param.breed$ParamName[which(param.breed$X == "fix_crossPA")] <- "Protected site x year"
param.breed$ParamName[which(param.breed$X == "fix_dist")] <- "Survey distance"
param.breed$ParamName[which(param.breed$X == "fix_dur")] <- "Survey duration"
param.breed$ParamName[which(param.breed$X == "hyperp")] <- "Zero-inflation hyperparameter"
param.breed$ParamName[which(param.breed$X == "US-AL")] <- "Alabama"
param.breed$ParamName[which(param.breed$X == "US-CT")] <- "Connecticut"
param.breed$ParamName[which(param.breed$X == "US-DE")] <- "Delaware"
param.breed$ParamName[which(param.breed$X == "US-FL")] <- "Florida"
param.breed$ParamName[which(param.breed$X == "US-GA")] <- "Georgia"
param.breed$ParamName[which(param.breed$X == "US-LA")] <- "Louisiana"
param.breed$ParamName[which(param.breed$X == "US-MA")] <- "Massachusetts"
param.breed$ParamName[which(param.breed$X == "US-MD")] <- "Maryland"
param.breed$ParamName[which(param.breed$X == "US-ME")] <- "Maine"
param.breed$ParamName[which(param.breed$X == "US-MS")] <- "Mississippi"
param.breed$ParamName[which(param.breed$X == "US-NC")] <- "North Carolina"
param.breed$ParamName[which(param.breed$X == "US-NH")] <- "New Hampshire"
param.breed$ParamName[which(param.breed$X == "US-NJ")] <- "New Jersey"
param.breed$ParamName[which(param.breed$X == "US-RI")] <- "Rhode Island"
param.breed$ParamName[which(param.breed$X == "US-SC")] <- "South Carolina"
param.breed$ParamName[which(param.breed$X == "US-TX")] <- "Texas"
param.breed$ParamName[which(param.breed$X == "US-VA")] <- "Virginia"
write.csv(param.breed, file="All_Breed_ParamSummary.csv")

# Winter
param.wint$WintParam <- paste0(round(param.wint$median, 2), " (", round(param.wint$lcl85, 2), ", ", round(param.wint$ucl85,2), ")")
param.wint$ParamName <- "Intercept"
param.wint$ParamName[which(param.wint$X == "fix_site")] <- "Stewardship site type"
param.wint$ParamName[which(param.wint$X == "fix_sitePA")] <- "Protected site type"
param.wint$ParamName[which(param.wint$X == "fix_year")] <- "Year"
param.wint$ParamName[which(param.wint$X == "fix_cross")] <- "Stewardship site x year"
param.wint$ParamName[which(param.wint$X == "fix_crossPA")] <- "Protected site x year"
param.wint$ParamName[which(param.wint$X == "fix_dist")] <- "Survey distance"
param.wint$ParamName[which(param.wint$X == "fix_dur")] <- "Survey duration"
param.wint$ParamName[which(param.wint$X == "hyperp")] <- "Zero-inflation hyperparameter"
param.wint$ParamName[which(param.wint$X == "US-AL")] <- "Alabama"
param.wint$ParamName[which(param.wint$X == "US-CT")] <- "Connecticut"
param.wint$ParamName[which(param.wint$X == "US-DE")] <- "Delaware"
param.wint$ParamName[which(param.wint$X == "US-FL")] <- "Florida"
param.wint$ParamName[which(param.wint$X == "US-GA")] <- "Georgia"
param.wint$ParamName[which(param.wint$X == "US-LA")] <- "Louisiana"
param.wint$ParamName[which(param.wint$X == "US-MA")] <- "Massachusetts"
param.wint$ParamName[which(param.wint$X == "US-MD")] <- "Maryland"
param.wint$ParamName[which(param.wint$X == "US-ME")] <- "Maine"
param.wint$ParamName[which(param.wint$X == "US-MS")] <- "Mississippi"
param.wint$ParamName[which(param.wint$X == "US-NC")] <- "North Carolina"
param.wint$ParamName[which(param.wint$X == "US-NH")] <- "New Hampshire"
param.wint$ParamName[which(param.wint$X == "US-NJ")] <- "New Jersey"
param.wint$ParamName[which(param.wint$X == "US-RI")] <- "Rhode Island"
param.wint$ParamName[which(param.wint$X == "US-SC")] <- "South Carolina"
param.wint$ParamName[which(param.wint$X == "US-TX")] <- "Texas"
param.wint$ParamName[which(param.wint$X == "US-VA")] <- "Virginia"
write.csv(param.wint, file="All_Winter_ParamSummary.csv")




##########################################################################
## SUMMARIZE PROPORTION DIFFERENCE IN TRENDS AND ABUNDANCE AMONG 
##    MANAGEMENT TYPES FOR EACH SPECIES AND SEASON FOR REPORTING
##########################################################################

abun.all <- read.csv(file="AllSpecies_AllAbundance.csv")

# summer
abun.summ <- abun.all[which(abun.all$Season == "Breed"),]

summ.abundiff <- data.frame()
for (spec in unique(abun.summ$Species)){
  temp <- abun.summ[which(abun.summ$Species==spec),]
  StewPDiff=((temp$X50.[which(temp$X=="StewAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  PAPDiff=((temp$X50.[which(temp$X=="PAAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  BothPDiff=((temp$X50.[which(temp$X=="BothAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  summ.abundiff <- rbind(summ.abundiff, data.frame(Species=spec, Season="Summer", StewPDiff=StewPDiff, PAPDiff=PAPDiff, BothPDiff=BothPDiff))
}
summ.abundiff

# winter
abun.wint <- abun.all[which(abun.all$Season == "Winter"),]

wint.abundiff <- data.frame()
for (spec in unique(abun.wint$Species)){
  temp <- abun.wint[which(abun.wint$Species==spec),]
  StewPDiff=((temp$X50.[which(temp$X=="StewAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  PAPDiff=((temp$X50.[which(temp$X=="PAAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  BothPDiff=((temp$X50.[which(temp$X=="BothAvgAb")] - temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")])
  wint.abundiff <- rbind(wint.abundiff, data.frame(Species=spec, Season="Winter", StewPDiff=StewPDiff, PAPDiff=PAPDiff, BothPDiff=BothPDiff))
}
wint.abundiff

##########
## TRENDS

trend.all <- read.csv(file="AllSpecies_AllTrends.csv")

# summer
trend.summ <- trend.all[which(trend.all$Season == "Breed"),]

summ.trenddiff <- data.frame()
for (spec in unique(trend.summ$Species)){
  temp <- trend.summ[which(trend.summ$Species==spec),]
  StewPDiff=((temp$X50.[which(temp$X=="StewTrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  PAPDiff=((temp$X50.[which(temp$X=="PATrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  BothPDiff=((temp$X50.[which(temp$X=="BothTrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  summ.trenddiff <- rbind(summ.trenddiff, data.frame(Species=spec, Season="Summer", StewPDiff=StewPDiff, PAPDiff=PAPDiff, BothPDiff=BothPDiff))
}
summ.trenddiff

# winter
trend.wint <- trend.all[which(trend.all$Season == "Winter"),]

wint.trenddiff <- data.frame()
for (spec in unique(trend.wint$Species)){
  temp <- trend.wint[which(trend.wint$Species==spec),]
  StewPDiff=((temp$X50.[which(temp$X=="StewTrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  PAPDiff=((temp$X50.[which(temp$X=="PATrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  BothPDiff=((temp$X50.[which(temp$X=="BothTrend")] - temp$X50.[which(temp$X=="NeitherTrend")])/abs(temp$X50.[which(temp$X=="NeitherTrend")]))
  wint.trenddiff <- rbind(wint.trenddiff, data.frame(Species=spec, Season="Winter", StewPDiff=StewPDiff, PAPDiff=PAPDiff, BothPDiff=BothPDiff))
}
wint.trenddiff



##########################################################################
##  ABUNDANCE - SUMMARIZE AND PLOT FIGURE S5
##########################################################################

## summarize abundance and abundance differences
setwd(pathtooutput)
abundiff <- data.frame()
abunall <- data.frame()
lf <- list.files(pattern = "*_AbundanceSummary.csv")
for (i in 1:length(lf)){
  temp <- read.csv(lf[i])
  spec <- substr(lf[i],1,4)
  abundiff <- rbind(abundiff, data.frame(Species=substr(lf[i],1,4), Season=substr(lf[i],6,11), StewDiff=(temp$X50.[which(temp$X=="StewAvgAb")]-temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")],
                                         PADiff=(temp$X50.[which(temp$X=="PAAvgAb")]-temp$X50.[which(temp$X=="NeitherAvgAb")])/temp$X50.[which(temp$X=="NeitherAvgAb")]))
  temp$Species <- substr(lf[i],1,4)
  temp$Season <- substr(lf[i],6,11)
  abunall <- rbind(abunall, temp)
  
}
abundiff$Season <- as.character(abundiff$Season)
abundiff$Season[which(abundiff$Season=="Breed_")] <- "Breed"
write.csv(abundiff, file=paste0(pathtooutput, "AllSpecies_AbundanceDifferences.csv"))
abunall$Season <- as.character(abunall$Season)
abunall$Season[which(abunall$Season=="Breed_")] <- "Breed"
write.csv(abunall, file=paste0(pathtooutput, "AllSpecies_AllAbundance.csv"))

##########################
# plot breeding abundance
abunbreed <- abunall[which(abunall$X != "BothAvgAb" & abunall$Season=="Breed"),]
names(abunbreed)[which(names(abunbreed)=="X")] <- "Site"
abunbreed <- abunbreed[order(abunbreed$Species, decreasing=FALSE),]
abunbreed$Species <- factor(abunbreed$Species)
abunbreed$siteord <- 1
abunbreed$siteord[which(abunbreed$Site=="NeitherAvgAb")] <- 2
abunbreed$siteord[which(abunbreed$Site=="PAAvgAb")] <- 3
abunbreed <- abunbreed[order(abunbreed$Species, abunbreed$siteord),]
pbrcols <- c("#1b9e77", "#d95f02", "#7570b3")
#### plot main figure
abunbreed$SpecNum <- rev(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3)))
png(file=paste0(pathtofigs,"Abundance_All_Breed.png"), width=6.4, height=4.74, units="in" , res=300)
p <- ggplot(abunbreed, aes(fill=rev(Site), y=X50., x=SpecNum)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  geom_vline(xintercept= 1.5, col="#CAD4D1")+
  geom_vline(xintercept= 2.5, col="#CAD4D1")+
  geom_vline(xintercept= 3.5, col="#CAD4D1")+
  geom_vline(xintercept= 4.5, col="#CAD4D1")+
  geom_vline(xintercept= 5.5, col="#CAD4D1")+
  geom_vline(xintercept= 6.5, col="#CAD4D1")+
  geom_vline(xintercept= 7.5, col="#CAD4D1")+
  scale_color_manual(name="", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship", "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_x_continuous(limits = c(0.7,8.3),breaks=c(1:8),
                     labels=c("American Oystercatcher","Black Skimmer","Brown Pelican", "Clapper Rail", "Least Tern", "Piping Plover", "Reddish Egret", "Snowy Plover"))+
  coord_flip()+
  labs(x="", y="Relative Abundance (median + 85% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
print(p)
dev.off()

# manually add legend 
legend <- cowplot::get_legend(p)
grid.newpage()
png(file=paste0(pathtofigs,"StewardshipTypeLegend.png"), width=3, height=2.5, units="in", res=300)
grid.draw(legend)
dev.off()

#### plot inset - omit BRPE
splevel <- c("SNPL","REEG","PIPL","LETE","CLRA","BLSK","AMOY")
abunbreed.inset <- abunbreed %>% dplyr::filter(Species!=c("BRPE"))
png(file=paste0(pathtofigs,"Abundance_All_Breed_Inset.png"), width=3.4, height=2.74, units="in" , res=300)
pinset <- ggplot(abunbreed.inset, aes(fill=rev(Site), y=X50., x=Species)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  scale_color_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship", "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_x_discrete(limits = splevel,breaks=c("AMOY","BLSK","CLRA","LETE", "PIPL", "REEG", "SNPL"),
                   labels=c("AMOY","BLSK", "CLRA", "LETE", "PIPL", "REEG", "SNPL"))+
  coord_flip()+
  labs(x="", y="") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y=element_text(size = 12), 
        axis.title.x=element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
pinset
dev.off()

#### produce donut plots
# Create new dataset for each donut - manually enter the number of species meeting each criteria
# shown for STW>UP
data <- data.frame( category=c("YES", "NO"),count=c(6,2))
# Compute percentages
data$fraction <- data$count / sum(data$count)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$category, "\n (n = ", data$count, ")")

# Make the plot - create each plot manually and name brd1 - brd4
# green (stw>up): "#1b9e77"; purple (pr>up): "#7570b3"; blue (stw>pr): "#629FD7"; gray (null): "#CAD4D1" 
brd1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=0.9, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer) #blue: #629FD7
  scale_fill_manual(values=c("#CAD4D1","#1b9e77")) + # first color is grey, change second color value
  scale_color_manual(values=c("#000000","#000000")) + # text color = black
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
# save strip of donut figures
png(file=paste0(pathtofigs,"Donut_Abundance_Breed.png"), width=11, height=3.74, units="in" , res=600)
grid.arrange(brd1,brd2,brd3,brd4, ncol=4)
dev.off()

##########################
# plot winter abundance
abunwint <- abunall[which(abunall$X != "BothAvgAb" & abunall$Season=="Winter"),]
names(abunwint)[which(names(abunwint)=="X")] <- "Site"
abunwint$Species <- as.character(abunwint$Species)
abunwint$Species[which(abunwint$Species=="REKN")] <- "REAN" # rename REKN to order correctly
abunwint <- abunwint[order(abunwint$Species, decreasing=FALSE),]
abunwint$Species <- factor(abunwint$Species)
abunwint$siteord <- 1
abunwint$siteord[which(abunwint$Site=="NeitherAvgAb")] <- 2
abunwint$siteord[which(abunwint$Site=="PAAvgAb")] <- 3
abunwint <- abunwint[order(abunwint$Species, abunwint$siteord),]
pbrcols <- c("#1b9e77", "#d95f02", "#7570b3")

# plot main figure
abunwint$SpecNum <- rev(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3), rep(9,3), rep(10,3), rep(11,3)))
png(file=paste0(pathtofigs,"Abundance_All_Winter.png"), width=6.4, height=4.74, units="in" , res=300)
p <- ggplot(abunwint, aes(fill=rev(Site), y=X50., x=Species)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  geom_vline(xintercept= 1.5, col="#CAD4D1")+
  geom_vline(xintercept= 2.5, col="#CAD4D1")+
  geom_vline(xintercept= 3.5, col="#CAD4D1")+
  geom_vline(xintercept= 4.5, col="#CAD4D1")+
  geom_vline(xintercept= 5.5, col="#CAD4D1")+
  geom_vline(xintercept= 6.5, col="#CAD4D1")+
  geom_vline(xintercept= 7.5, col="#CAD4D1")+
  geom_vline(xintercept= 8.5, col="#CAD4D1")+
  geom_vline(xintercept= 9.5, col="#CAD4D1")+
  geom_vline(xintercept= 10.5, col="#CAD4D1")+
  scale_color_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship","Unprotected","Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship","Unprotected","Protected"), values=c(rev(pbrcols))) +
  scale_x_discrete(limits = rev(levels(abunwint$Species)),breaks=c("AMOY","BLSK","BRPE", "CLRA", "LBCU", "MAGO", "PIPL", "REAN", "REEG", "SNPL", "WESA"),
                   labels=c("American Oystercatcher","Black Skimmer","Brown Pelican", "Clapper Rail", "Long-billed Curlew", "Marbled Godwit", "Piping Plover", "Red Knot", "Reddish Egret", "Snowy Plover", "Western Sandpiper"))+
  coord_flip()+
  labs(x="", y="Relative Abundance (median + 85% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
print(p)
dev.off()

# Plot inset, omit BLSK BRPE WESA 
splevel <- c("SNPL","REEG","REAN","PIPL","MAGO","LBCU","CLRA","AMOY")
abunwint1 <- abunwint %>% dplyr::filter(Species!=c("BRPE")) %>% dplyr::filter(Species!=c("BLSK")) %>% dplyr::filter(Species!=c("WESA"))
png(file=paste0(pathtofigs,"Abundance_All_Winter_Inset.png"), width=3.7, height=3.04, units="in" , res=300)
pinset <- ggplot(abunwint1, aes(fill=rev(Site), y=X50., x=Species)) + 
  geom_point(aes(color=rev(Site)), size=2.5, shape=21, stroke=1.4, col="black", alpha=0.8, position=position_dodge(width=0.9), stat="identity")+
  geom_errorbar(aes(ymin=X7.5., ymax=X92.5., color=rev(Site)), position=position_dodge(width=0.9), width=.1, size=1) +
  scale_color_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship", "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_fill_manual(name="Protected\n Area Type", breaks=c("StewAvgAb", "NeitherAvgAb", "PAAvgAb"), labels=c("Stewardship",  "Unprotected", "Protected"), values=c(rev(pbrcols))) +
  scale_x_discrete(limits = splevel,breaks=c("AMOY","CLRA", "LBCU", "MAGO", "PIPL", "REAN", "REEG", "SNPL"),
                   labels=c("AMOY","CLRA", "LBCU", "MAGO", "PIPL", "REKN", "REEG", "SNPL"))+
  coord_flip()+
  labs(x="", y="Relative Abundance") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y=element_text(size = 12), 
        axis.title.x=element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
pinset
dev.off()
#### produce donut plots
# Create new dataset for each donut - manually enter the number of species meeting each criteria
# shown for STW>UP
data <- data.frame( category=c("YES", "NO"),count=c(5,6))
# Compute percentages
data$fraction <- data$count / sum(data$count)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)
# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# Compute a good label
data$label <- paste0(data$category, "\n (n = ", data$count, ")")

# Make the plot - create each plot manually and name brd1 - brd4
# green (stw>up): "#1b9e77"; purple (pr>up): "#7570b3"; blue (stw>pr): "#629FD7"; gray (null): "#CAD4D1" 
win1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text(x=0.9, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer) #blue: #629FD7
  scale_fill_manual(values=c("#CAD4D1","#1b9e77")) + # first color is grey, change second color value
  scale_color_manual(values=c("#000000","#000000")) + # text color = black
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

# save strip of donut figures
png(file=paste0(pathtofigs,"Donut_Abundance_Winter.png"), width=11, height=3.74, units="in" , res=600)
grid.arrange(win1,win2,win3,win4, ncol=4)
dev.off()


##########################################################################
##  SUMMARIZE GOODNESS OF FIT STATISTICS FOR TABLE S6
##########################################################################

## summarize GOF stats
gofall <- data.frame()
lf <- list.files(pattern = "*_GOFstats.csv")
for (i in 1:length(lf)){
  temp <- read.csv(lf[i])
  names(temp) <- c("X", "X50.","X2.5.", "X97.5.", "X25.", "X75." ,   "X7.5." ,  "X92.5.")
  temp$Species <- substr(lf[i],1,4)
  temp$Season <- substr(lf[i],6,11)
  gofall <- rbind(gofall, temp)
}
gofall$Season <- as.character(gofall$Season)
gofall$Season[which(gofall$Season=="Breed_")] <- "Breed"
write.csv(gofall, file=paste0(pathtooutput, "AllSpecies_AllGOFStats.csv"))

summary(gofall[which(gofall$X=="spcor"),c("X50.")]) 
summary(gofall[which(gofall$X=="chats"),c("X50.")]) 
summary(gofall[which(gofall$X=="fdifs"),c("X50.")]) 
summary(gofall[which(gofall$X=="obsdifs"),c("X50.")])
summary(gofall[which(gofall$X=="meandifs"),c("X50.")])



##########################################################################
## PRODUCE RANDOM EFFECTS PLOTS FOR APPENDIX S7
##########################################################################

setwd("Spatial")
# load a shapefile containing US states, e.g., https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip
stprov <- readOGR(".","US_States")
stprov2 <- stprov
stprov2@data$id <- rownames(stprov2@data)
# create a data.frame from our spatial object
stprovs <- fortify(stprov2, region = "id")
# merge the "fortified" data with the data from our spatial object
stprovDF <- merge(stprovs, stprov2@data, by = "id")

# For strata plotting
stprov3<-stprov2
stprov.fort <- fortify(stprov3, region = "ObjectID")
idList <- stprov3@data$ObjectID
centroids.df <- as.data.frame(coordinates(stprov3))
names(centroids.df) <- c("Longitude", "Latitude")  

## load random effect parameter summaries ##--------------------------
state.lookup <- stprovDF[!(duplicated(stprovDF$ObjectID)),c("ObjectID","GMI_ADMIN")]
re.breed <- read.csv(file="All_Breed_ParamSummary.csv")
re.breed <- re.breed[which(re.breed$X %in% c("US-AL", "US-CT", "US-DE", "US-FL", "US-GA", "US-LA", 
                                             "US-MA", "US-MD", "US-ME", "US-MS", "US-NC", "US-NH", 
                                             "US-NJ", "US-NY", "US-RI", "US-SC", "US-TX", "US-VA")),]
re.breed$ObjectID <- 0
re.breed$ObjectID[which(re.breed$X=="US-AL")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-ALB")]
re.breed$ObjectID[which(re.breed$X=="US-CT")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-CNN")]
re.breed$ObjectID[which(re.breed$X=="US-DE")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-DEL")]
re.breed$ObjectID[which(re.breed$X=="US-FL")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-FLA")]
re.breed$ObjectID[which(re.breed$X=="US-GA")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-GEO")]
re.breed$ObjectID[which(re.breed$X=="US-LA")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-LOU")]
re.breed$ObjectID[which(re.breed$X=="US-MA")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-MSS")]
re.breed$ObjectID[which(re.breed$X=="US-MD")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-MRY")]
re.breed$ObjectID[which(re.breed$X=="US-ME")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-MAI")]
re.breed$ObjectID[which(re.breed$X=="US-MS")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-MSP")]
re.breed$ObjectID[which(re.breed$X=="US-NC")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-NCR")]
re.breed$ObjectID[which(re.breed$X=="US-NH")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-NHM")]
re.breed$ObjectID[which(re.breed$X=="US-NJ")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-NJR")]
re.breed$ObjectID[which(re.breed$X=="US-NY")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-NYO")]
re.breed$ObjectID[which(re.breed$X=="US-RI")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-RHI")]
re.breed$ObjectID[which(re.breed$X=="US-SC")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-SCR")]
re.breed$ObjectID[which(re.breed$X=="US-TX")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-TEX")]
re.breed$ObjectID[which(re.breed$X=="US-VA")] <- state.lookup$ObjectID[which(state.lookup$GMI_ADMIN=="USA-VRG")]
# identify significant REs
re.breed$Sig <- 0
re.breed$Sig[which(re.breed$lcl85<0 & re.breed$ucl85<0)] <- 1
re.breed$Sig[which(re.breed$lcl85>0 & re.breed$ucl85>0)] <- 2

re.wint <- read.csv(file="All_Winter_ParamSummary.csv")
re.wint <- re.wint[which(re.wint$X %in% c("US-AL", "US-CT", "US-DE", "US-FL", "US-GA", "US-LA", 
                                          "US-MA", "US-MD", "US-ME", "US-MS", "US-NC", "US-NH", 
                                          "US-NJ", "US-NY", "US-RI", "US-SC", "US-TX", "US-VA")),]
re.wint <- merge(re.wint, re.breed[!(duplicated(re.breed$ObjectID)),c("ObjectID","X")], by="X")
# identify significant REs
re.wint$Sig <- 0
re.wint$Sig[which(re.wint$lcl85<0 & re.wint$ucl85<0)] <- 1
re.wint$Sig[which(re.wint$lcl85>0 & re.wint$ucl85>0)] <- 2

speclist <- c(  "AMOY","AMOY","BLSK","BLSK","BRPE","BRPE","CLRA","CLRA","LBCU","LETE","MAGO","PIPL","PIPL","REKN","REEG","REEG","SNPL","SNPL","WESA")  
seasonlist <- c("Summer","Winter","Summer","Winter","Summer","Winter","Summer","Winter","Winter","Summer","Winter","Summer","Winter","Winter","Summer","Winter","Summer","Winter","Winter")
speclookup <- data.frame(Spec=unique(speclist), SpecNames=c("American Oystercatcher", "Black Skimmer", "Brown Pelican", "Clapper Rail", "Long-billed Curlew",
                                                            "Least Tern", "Marbled Godwit", "Piping Plover", "Red Knot", "Reddish Egret", "Snowy Plover", "Western Sandpiper"))

### COMBINE plot and map for each SPECIES-SEASON
for (i in 1:length(speclist)){
  spec <- speclist[i]
  season <- seasonlist[i]
  
  if (season=="Summer"){
    tempspec <- re.breed[which(re.breed$Spec==spec),]
  } else {
    tempspec <- re.wint[which(re.wint$Spec==spec),]
  }
  tempspec <- tempspec[order(tempspec$ParamName),]
  
  re.means <- tempspec$median
  re.lc <- tempspec$lcl85
  re.uc <- tempspec$ucl85
  
  stratas <- tempspec$ObjectID
  sig.stratas <- tempspec$Sig
  statenames <- tempspec$ParamName
  ObjectID <- tempspec$ObjectID
  
  #add season and species columns for faceting
  Season <- c(rep(season,length(stratas)))
  Species <- c(rep(spec, length(stratas)))
  
  re.df <- data.frame(re.means,re.lc,re.uc,stratas,sig.stratas,Season,Species, statenames, ObjectID)
  re.df$stateorder <- seq(1:nrow(re.df))
  
  replot <- ggplot(data = re.df, aes(x = reorder(statenames, rev(stateorder)), col = factor(sig.stratas))) +
    geom_errorbar(aes(ymin=re.lc, ymax=re.uc),size=0.9,width=0)+
    geom_point(aes(y=re.means),size=3, shape=21, stroke=1.4,fill="White")+
    scale_color_manual(values = c("grey72", "#68023F", "#008169"))+ 
    coord_flip()+
    labs(y="State-specific Intercept",x="")+
    theme_bw() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y=element_text(size = 12), 
          axis.title.x=element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black"),
          panel.background = element_rect(fill = 'white', colour = 'white'))+
    theme(legend.position="none")+
    theme(plot.margin=unit(c(1,-1,0,0),"cm"))

  #Add in map-----------------------------------------------------------
  
  #create subsets of stprovDF for each significance level
  highocc <- stprovDF[which(stprovDF$ObjectID %in% re.df$ObjectID[which(re.df$sig.stratas=="2")]),] # sig higher abundance
  lowocc <- stprovDF[which(stprovDF$ObjectID %in% re.df$ObjectID[which(re.df$sig.stratas=="1")]),] # sig lower abundance
  highoccn <- stprovDF[which(stprovDF$ObjectID %in% re.df$ObjectID[which(re.df$sig.stratas=="0")]),] # not sig
  
  highocc$Season <- rep(season, nrow(highocc))
  lowocc$Season <- rep(season, nrow(lowocc))
  highoccn$Season <- rep(season, nrow(highoccn))
  
  highocc$Species <- rep(spec, nrow(highocc))
  lowocc$Species <- rep(spec, nrow(lowocc))
  highoccn$Species <- rep(spec, nrow(highoccn))
  
  # plot map
  remap <- ggplot(data = stprovDF, aes(x=long, y=lat)) +
    geom_polygon(aes(group=group),fill="white")  +
    geom_polygon(aes(group=group),data = highocc, fill = "#008169") + 
    geom_polygon(aes(group=group),data = lowocc, fill ="#68023F") + 
    geom_polygon(aes(group=group),data = highoccn, fill ="grey72") +
    geom_path(aes(group=group),color = "darkgray",size=0.2) +
    coord_equal() +
    theme_bw()+
    theme(plot.margin=unit(c(1,0.9,0,-2),"cm"))+
    theme(legend.position = "none", title = element_blank(),
          axis.text = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())

  
  
  fig<- ggarrange(replot,remap, ncol=2, nrow=1, widths=c(1, 2), align="v")
  annotate_figure(fig,
                  top = grobTree(rectGrob(gp = gpar(fill = "gray86"),height=2.7,width=0.75,vjust=1.2), text_grob(season, color = "black", size = 12, face = "bold", family = "sans", vjust=1.9)),
                  right = grobTree(rectGrob(gp = gpar(fill = "gray86"),width=2.7,height=0.75,hjust=1.2), text_grob(speclookup$SpecNames[which(speclookup$Spec==spec)], color = "black", family="sans",size = 12,rot=-90, vjust=1.9))
  )
  ggsave(paste0("Figs/", spec, "_", season, "_RE_MapFig.png"), width = 11, height = 5, units = "in",dpi=600)
  
  
} # end speclist loop


