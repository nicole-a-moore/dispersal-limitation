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
## download species-level climate velocity data from BioShiftR and reformat to match data 

## read in species shifts and cv 
cvs_spp = read.csv("data-raw/BIOSHIFTSv3/cvs/sp_shifts_and_cvs_for_nikki.csv")
## note: this is just birds and plants 

## ID, sps, res, type
## get rid of unnecessary columns 
cvs_spp = cvs_spp %>%
  select(old_id, id, eco, sp_name_checked, type, colnames(.)[33:55]) %>%
  rename("ID" = old_id, "NewID" = id) 

## reformat data to create resolution column 
cvs_spp = cvs_spp %>%
  gather(key = "cv_type", value = "cv", 
         c("cv_temp_mean_res1km", "cv_temp_sd_res1km", "cv_temp_q25_res1km", "cv_temp_median_res1km", "cv_temp_q75_res1km",     
           "cv_temp_mean_res25km","cv_temp_sd_res25km","cv_temp_q25_res25km","cv_temp_median_res25km" ,
           "cv_temp_q75_res25km","cv_temp_mean_res50km","cv_temp_sd_res50km","cv_temp_q25_res50km"   , 
           "cv_temp_median_res50km","cv_temp_q75_res50km","cv_temp_mean_res110km","cv_temp_sd_res110km"  ,  
           "cv_temp_q25_res110km","cv_temp_median_res110km","cv_temp_q75_res110km")) %>%
  mutate(cv_res = str_split_fixed(cv_type, "_res", 2)[,2],
         cv_type = str_split_fixed(cv_type, "_res", 2)[,1]) %>%
  distinct()

## spread back 
cvs_spp = cvs_spp %>%
  spread(key = cv_type, value = cv)

## rename to match old names
cvs_spp = cvs_spp %>%
  rename("mean_cv_spp" = cv_temp_mean,
         "median_cv_spp" = cv_temp_median,
          "q1_cv_spp" = cv_temp_q25,
          "q3_cv_spp" = cv_temp_q75,
          "sd_cv_spp" = cv_temp_sd)

## get rid of marine studies
cvs_spp = cvs_spp[which(cvs_spp$eco == "Ter"), ]

## save 
write.csv(cvs_spp, "data-processed/v3_with-cv.csv", row.names = FALSE)

## get range of SD of temps across study areas 
#devtools::install_github("Bioshifts/BioShiftR", force=T)
library(BioShiftR)

shifts = get_shifts() %>%
  add_cv(type = "SA",
         stat = "sd")

min(shifts$cv_temp_sd, na.rm = T) # 2.805362e-07
max(shifts$cv_temp_sd, na.rm = T) # 21.51143
