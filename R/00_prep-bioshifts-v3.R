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
write.csv(v3, "data-processed/intermediate_files/bioshifts/v3_shifts.csv", row.names = FALSE)

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
write.csv(cvs_spp, "data-processed/intermediate_files/bioshifts/v3_lat-spp-specific-cvs.csv", row.names = FALSE)

## get range of SD of temps across study area 
cvs_study = read.csv("data-raw/BIOSHIFTSv3/clim_vel/vel_SA_all.csv")

min(cvs_study$VelLat_sd_mat, na.rm = T)
max(cvs_study$VelLat_sd_mat, na.rm = T)
