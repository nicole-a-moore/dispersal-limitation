## join bioshfitsv3 with dispersal distance to latitudinal climate velocities 
## select the most appropriate spatial resolution of climate velocity based on dispersal distance 
library(tidyverse)


#############################
####   read in v3 data   ####
#############################
## read in bioshfifts with dispersal distance estimates 
v3 <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitudinal shifts 
v3 <- filter(v3, Type == "LAT")
length(unique(v3$scientificName_checked)) ## 399 species 
v3 %>% 
  select(group, scientificName_checked) %>%
  distinct() %>%
  group_by(group) %>%
  tally()

######################################
####   prep climate velocity data ####
######################################
## read in file with preliminary climate velocities across different scales
## includes mean across study area and across each species range within study area 
cvs <- read.csv("~/Documents/bioshifts-traits/data-processed/v3_new-cvs_preliminary.csv")

cvs <- cvs %>%
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
         mean_cv_sppspecific, sd_cv_sppspecific, p90_cv_sppspecific, p10_cv_sppspecific, cv_res) %>%
  distinct()

## save 
write.csv(cvs, "data-processed/v3_lat-spp-specific-cvs.csv", row.names = FALSE)



##############################################
####   join v3 with climate velocity data ####
##############################################
## first, make column that represents climate velocity calculated at the spatial scale that is most relevant 
## if dispersal distance < 25km, then 25 km
## if 25 km < dispersal distance < 50km, then 50 km
## if dispersal distance > 50km, then 110 km

## first, get rid of mean study level columns 
cvs <- select(cvs, -mean_cv_studylevel, -sd_cv_studylevel)

v3 <- v3 %>% 
  mutate(cv_res = ifelse(MaxDispersalDistanceKm <= (25/2+35.355), "25km",
                            ifelse(MaxDispersalDistanceKm <= (50/2+70.71), "50km",
                                   ifelse(MaxDispersalDistanceKm > (50/2+70.71), "110km",
                                          NA)))) %>%
  select(MaxDispersalDistanceKm, cv_res, everything()) %>%
  left_join(., cvs) %>%
  rename("ClimVeloKmY_RelScale" = mean_cv_sppspecific, 
         "sdClimVeloKmY_RelScale" = sd_cv_sppspecific, 
         "p90ClimVeloKmY_RelScale" = p90_cv_sppspecific, 
         "p10ClimVeloKmY_RelScale" = p10_cv_sppspecific,
         "RelScale" = cv_res) %>%
  select(ClimVeloKmY_RelScale, sdClimVeloKmY_RelScale, p90ClimVeloKmY_RelScale, p10ClimVeloKmY_RelScale, RelScale, everything())

## now add all other scales 
## pivot wider: 
cvs <- cvs %>%
  rename("ClimVeloKmY" = mean_cv_sppspecific, "sdClimVeloKmY" = sd_cv_sppspecific, 
         "p90ClimVeloKmY" = p90_cv_sppspecific, "p10ClimVeloKmY" = p10_cv_sppspecific) %>%
  pivot_wider(names_from = cv_res, 
              values_from = c(ClimVeloKmY, sdClimVeloKmY, p90ClimVeloKmY, p10ClimVeloKmY),
              names_glue = "{.value}_{cv_res}") %>%
  select(scientificName_checked, ID, Type, Eco,
         ClimVeloKmY_25km, sdClimVeloKmY_25km, p90ClimVeloKmY_25km, p10ClimVeloKmY_25km,
         ClimVeloKmY_50km, sdClimVeloKmY_50km, p90ClimVeloKmY_50km, p10ClimVeloKmY_50km,
         ClimVeloKmY_110km, sdClimVeloKmY_110km, p90ClimVeloKmY_110km, p10ClimVeloKmY_110km)

## make sure each unique species x study area combination has only 1 row 
length(unique(paste(cvs$scientificName_checked, cvs$ID))) == nrow(cvs) # TRUE
## :-)

## now merge to v3
v3_merged <- left_join(v3, cvs)

## make sure no observations were duplicated
nrow(v3) == nrow(v3_merged) # TRUE

## rearrange columns 
v3_merged <- v3_merged %>%
  select(scientificName_checked, ID, Type, Eco, Param, group, 
         ShiftKmY, DispersalPotentialKmY, MedianDispersalPotentialKmY,
         ClimVeloKmY_RelScale, sdClimVeloKmY_RelScale, p90ClimVeloKmY_RelScale, p10ClimVeloKmY_RelScale, RelScale, 
         ClimVeloKmY_25km, sdClimVeloKmY_25km, p90ClimVeloKmY_25km, p10ClimVeloKmY_25km,
         ClimVeloKmY_50km, sdClimVeloKmY_50km, p90ClimVeloKmY_50km, p10ClimVeloKmY_50km,
         ClimVeloKmY_110km, sdClimVeloKmY_110km, p90ClimVeloKmY_110km, p10ClimVeloKmY_110km, everything())
## save 
write.csv(v3_merged, "data-processed/v3_with-cv.csv", row.names = FALSE)

