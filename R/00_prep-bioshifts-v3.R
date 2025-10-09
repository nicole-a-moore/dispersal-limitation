## add climate velocities to bioshifts v3
library(tidyverse)

## read in bioshifts v3
v3 = read.csv("data-raw/BIOSHIFTSv3/BIOSHIFTS_v3.csv")
v3$scientificName_checked = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$scientificName = str_replace_all(v3$sp_name_database, "\\_", " ")
v3$scientificName

## fix ecosystem column 
## everything != marine is ter
v3$Eco = ifelse(!v3$Eco %in% c("Mar", "Ter"), "Ter", v3$Eco)

## except: A164_P1 = "Mar"
v3$Eco[which(v3$ID %in% c("A164_P1", "A183_P1"))] = "Mar"

## save 
write.csv(v3, "data-processed/v3_shifts.csv", row.names = FALSE)

####################################################
####   prep species-level climate velocity data ####
####################################################
## prep species-specific climate velocities 
## read in file with preliminary climate velocities 
cvs_spp = read.csv("data-raw/BIOSHIFTSv3/clim_vel/vel_SA_sps_all.csv")

## get rid of velocities for just elevation studies
key = v3 %>%
  select(ID, Type) %>%
  distinct() %>% 
  filter(Type != "ELE")

cvs_spp = cvs_spp[cvs_spp$ID %in% key$ID,]

## reformat and rename species scienfitic name column so it matches v3
cvs_spp = cvs_spp %>%
  rename("scientificName_checked" = sps) %>%
  mutate(scientificName_checked = str_replace_all(scientificName_checked, "\\_", " "))

colnames(cvs_spp)

cvs_spp = cvs_spp[,which(!str_detect(colnames(cvs_spp), "Ele"))] ## get rid of elevation velocities
cvs_spp = cvs_spp[,which(!str_detect(colnames(cvs_spp), "map"))] ## get rid of precipitation velocities
cvs_spp = cvs_spp[,which(!str_detect(colnames(cvs_spp), "baseline"))] ## get rid of baseline temps
cvs_spp = cvs_spp[,which(!str_detect(colnames(cvs_spp), "trend"))] ## get rid of trends 
cvs_spp = cvs_spp[,which(!str_detect(colnames(cvs_spp), "Vel_"))]## get rid of velocity that contains longitudinal component

colnames(cvs_spp)

## reformat data 
cvs_spp <- cvs_spp %>%
  gather(key = cv_type, value = cv_spp, colnames(cvs_spp)[4:ncol(cvs_spp)]) %>%
  filter(!is.na(cv_spp)) %>%
  mutate(spp_id = paste(scientificName_checked, ID), 
         Type = "LAT") %>%
  mutate(Eco = ifelse(str_detect(cv_type, "sst"), "Mar", 
                      ifelse(str_detect(cv_type, "mat"), "Ter", NA))) %>%
  mutate(cv_type = str_split_fixed(cv_type, "\\_", 3)[,2]) %>%
  rename("cv_res" = res) %>%
  select(ID, scientificName_checked, Eco, Type, cv_res, cv_type, cv_spp) %>%
  mutate(cv_type = paste0(cv_type, "_cv_spp")) %>%
  spread(key = "cv_type", value = "cv_spp") %>%
  rename("q1_cv_spp" = `1Q_cv_spp`, "q3_cv_spp" = `3Q_cv_spp`) %>%
  distinct()

## get rid of marine studies
cvs_spp = cvs_spp[which(cvs_spp$Eco == "Ter"), ]

## save 
write.csv(cvs_spp, "data-processed/v3_lat-spp-specific-cvs.csv", row.names = FALSE)

## get range of SD of temps across study area 
cvs_study = read.csv("data-raw/BIOSHIFTSv3/clim_vel/vel_SA_all.csv")

min(cvs_study$VelLat_sd_mat, na.rm = T)
max(cvs_study$VelLat_sd_mat, na.rm = T)

# GARBAGE

####################################################
####   prep species-level climate velocity data ####
####################################################
## prep species-specific climate velocities 
## read in file with preliminary climate velocities 
cvs_spp <- read.csv("~/Documents/bioshifts-traits/data-processed/new-cvs_preliminary.csv")

cvs_spp <- cvs_spp %>%
  select(-range_source) %>%
  mutate(scientificName_checked = str_split_fixed(species_studyid, "\\_", 2)[,1],
         ID = str_split_fixed(species_studyid, "\\_", 2)[,2]) %>%
  filter(!str_detect(cv_type, "map")) %>% ## get rid of precipitation
  filter(!str_detect(cv_type, 'Ele')) %>% ## get rid of elevation
  mutate(cv_res =  ifelse(str_detect(cv_type, "25km"), "25km", 
                          ifelse(str_detect(cv_type, "50km"), "50km",
                                 ifelse(str_detect(cv_type, "1km"), "1km", 
                                        ifelse(str_detect(cv_type, "110km"), "110km", 
                                               NA)))),
         Eco = ifelse(str_detect(cv_res, "sst"), "Mar", "Ter")) %>%
  select(scientificName_checked, ID, Type, Eco, mean_cv_studylevel, sd_cv_studylevel, 
         mean_cv_sppspecific, sd_cv_sppspecific, cv_res) %>%
  distinct()


## save 
write.csv(cvs_spp, "data-processed/v3_lat-spp-specific-cvs.csv", row.names = FALSE)

## get names of all csv files in the climate velocity folder
files <- list.files(path = "~/Documents/bioshifts-traits/data-raw/bioshiftsv3/Velocity_SA", pattern="csv$", recursive=TRUE, full.names=TRUE)

## filter out elevation studies
sas <- list.files(path = "~/Documents/bioshifts-traits/data-raw/bioshiftsv3/Velocity_SA", pattern="csv$", recursive=TRUE)
sas <- paste0(str_split_fixed(sas, "\\_", 3)[,1], "_", str_split_fixed(sas, "\\_", 3)[,2])
files = files[which(sas %in% v3$ID[which(v3$Type == "LAT")])]

## pick out columns to keep for terrestrial/marine study areas 
terr_cols <- c("ID", "v.lat.mean.mat", # mean latitudinal velocity of air temperature 
               "v.ele.mean.mat", # mean elevational velocity of air temperature
               "v.lat.sd.mat", # sd of both 
               "v.ele.sd.mat")
mar_cols <- c("ID", "v.lat.mean.sst",
              "v.lat.sd.sst")

cvs <- c()
for(f in 1:length(files)) {
  
  ## read in a file 
  file <- read.csv(files[f])
  
  ## get scale
  scale = ifelse(str_detect(files[f], "25km"), "25km", 
                 ifelse(str_detect(files[f], "50km"), "50km",
                        ifelse(str_detect(files[f], "1km"), "1km", 
                               ifelse(str_detect(files[f], "110km"), "110km", 
                                      NA))))
  
  ## if study is terrestrial, extract air temperature columns for ele and lat 
  if(any(str_detect(colnames(file), "mat"))) {
    
    file <- select(file, c("ID", "v.lat.mean.mat", "v.lat.sd.mat"))
    
    ## rearrange
    file <- gather(file, key = "cv_type", value = "val", 
                   c(v.lat.mean.mat, v.lat.sd.mat)) %>%
      mutate(Type = "LAT",
             Eco = "Ter") %>%
      mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                              "sd_cv_studylevel")) %>%
      select(-cv_type) %>%
      spread(key = measure, value = val) %>%
      mutate(cv_res = scale)
    
    ## bind 
    cvs <- rbind(cvs, file)
  } 
  else {
    file <- select(file, mar_cols)
    
    ## rearrange
    file <- file %>%
      mutate(Type = "LAT",
             Eco = "Mar") %>%
      rename("mean_cv_studylevel" = v.lat.mean.sst,
             "sd_cv_studylevel" = v.lat.sd.sst) %>%
      mutate(cv_res = scale)
    
    ## bind 
    cvs <- rbind(cvs, file)
  }
  
  print(paste0("On file number ", f))
}

## get rid of 1km cvs
cvs <- filter(cvs, cv_res != "1km")

## save 
write.csv(cvs, "data-processed/v3_lat-study-level-cvs.csv", row.names = FALSE)

## make sure all studies are there 
unique(v3$ID[which(v3$Type == "LAT")])[which(!unique(v3$ID[which(v3$Type == "LAT")]) %in% cvs$ID)]
## missing A66_P1
