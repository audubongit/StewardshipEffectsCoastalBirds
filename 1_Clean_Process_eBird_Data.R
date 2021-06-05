################################################################################
## 1_CLEAN_PROCESS_EBIRD_DATA.R
##
## This code extracts and processes eBird data for trend analysis. Before
##      running, download and unzip the latest eBird Basic Dataset (EBD) and
##      sampling dataset, and install gawk (https://www.gnu.org/software/gawk/).
##
## R code accompanies the manuscript:
## Michel, N.L., S.P. Saunders, T.D. Meehan, C.B. Wilsey. 2021. Effects of
##      stewardship on protected area effectiveness for coastal birds. 
##      Conservation Biology, https://doi.org/10.1111/cobi.13698. 
################################################################################

# load tidyverse and auk
library(tidyverse)
library(auk)


# set working directory and see what's inside
ebd_dir <- "~/eBird/Sept2019Raw"
setwd(ebd_dir)
dir()

# set ebird ebd count file
f_ebd <- file.path(ebd_dir, "ebd_relSep-2019.txt")

# set ebird sampling file
f_sampling <- file.path(ebd_dir, "ebd_sampling_relSep-2019.txt")

# define the stuff you want to keep
target_bbox <- c()
target_bcr <- c()
target_country <- c("US")
target_date <- c("2007-01-01", "2018-12-31")
target_distance <- c(0, 8.1)
target_duration <- c(0, 180)
target_extent <- c()
target_protocol <- c("Traveling","Stationary","Area", "Random","International Shorebird Survey (ISS)",
                      "Coastal Shorebird Survey","Audubon Coastal Bird Survey","Great Texas Birding Classic",
                      "California Brown Pelican Survey","Greater Gulf Refuge Waterbird Count",
                     "TNC California Waterbird Count","Waterbird Count") #"eBird Pelagic Protocol" to include for SOSH

target_species <- c("Long-billed Curlew", "Marbled Godwit")

#target_species <- c("Brown Pelican", "Reddish Egret", "Clapper Rail", 
#                    "Piping Plover", "Snowy Plover", "American Oystercatcher", 
#                    "Hudsonian Godwit", "Red Knot", "Western Sandpiper", 
#                    "Semipalmated Sandpiper", "Least Tern", "Black Skimmer", "Saltmarsh Sparrow")
target_state <- c()
target_time <- c()

# create your filters that will be filter the ebd and sampling files
ebd_filters <- auk_ebd(f_ebd, file_sampling=f_sampling) %>% 
  auk_country(country=target_country) %>% 
  auk_date(date=target_date) %>% 
  auk_distance(distance=target_distance) %>%
  auk_duration(duration=target_duration) %>%
  auk_protocol(protocol = target_protocol) %>%
  auk_species(species=target_species) %>%
  #auk_state(state=target_state) %>%
  auk_complete()

# define the filtered output files and run the filtering step
f_out_ebd <- "./script_output/ebd_lbcu_mago.txt"
f_out_sampling <- "./script_output/ebd_sampling_lbcu_mago.txt"
ebd_sed_filtered <- auk_filter(ebd_filters, file=f_out_ebd, 
                           file_sampling=f_out_sampling,
                           overwrite=T)

# you now have a cleaned ebd and ebd_sampling text file in the script output 
# folder in the working directory. if you just want data for where the birds 
# occur, you can open the ebd file with auk_read() and it will remove group 
# duplicates and do the taxonomic rollup step. if you want to zero fill the
# ebd file using the ebd_sampling file, then continue. in the case, below, the
# ebd file is filtered for summer tanager, then read and cleaned like when 
# using the auk_read() function, then zero filled, and saved as a csv.

f_out_ebd <- read_ebd("./script_output/ebd_lbcu_mago.txt") #if crashes, need to read in saved ebd files again
f_out_sampling <- read_sampling("./script_output/ebd_sampling_lbcu_mago.txt") #if crashes, need to read in saved ebd files again

# species specific zero-filling (fills in 0s for complete checklists where
# species was not observed)
mago_zf_df <- auk_zerofill(x=f_out_ebd, sampling_events=f_out_sampling,
                           species= "Marbled Godwit",
                           collapse=T, unique=T, rollup=T, drop_higher=T,
                           complete=T)

# format dates
library(lubridate)
mago_zf_df$date_start <- as.Date(mago_zf_df$observation_date, format="%Y-%m-%d")
mago_zf_df$Jdate <- yday(mago_zf_df$date_start)

# separate into breeding and winter season datasets
ebd.breed <- mago_zf_df[which(mago_zf_df$Jdate > 151 & mago_zf_df$Jdate < 244),] #standard BBS: 134 - 182
ebd.wint <- mago_zf_df[which(mago_zf_df$Jdate > 326 | mago_zf_df$Jdate < 82),]  #standard CBC: > 334 | < 32
ebd.summ <- rbind(ebd.breed, ebd.wint)

# save cleaned, zero-filled file as CSV
library(data.table)
fwrite(ebd.wint, "./script_output/MAGO_zf.csv", row.names=F)

#then remove species data frames before moving on to next species to clear memory
rm(mago_zf_df,ebd.breed,ebd.wint,ebd.summ)
rm(f_out_ebd, f_out_sampling)




