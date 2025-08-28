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
length(unique(v3$scientificName_checked)) ## 398 species 
v3 %>% 
  select(group, scientificName_checked) %>%
  distinct() %>%
  group_by(group) %>%
  tally()

######################################
####   prep climate velocity data ####
######################################
## read in file with species-specific climate velocities across different scales
## includes mean across study area and across each species range within study area 
cvs <- read.csv("data-processed/v3_lat-spp-specific-cvs.csv")

##############################################
####   join v3 with climate velocity data ####
##############################################
## first, make column that represents climate velocity calculated at the spatial scale that is most relevant 
## if dispersal distance < 25km, then 25 km
## if 25 km < dispersal distance < 50km, then 50 km
## if dispersal distance > 50km, then 110 km
v3 <- v3 %>% 
  mutate(cv_res = ifelse(MaxDispersalDistanceKm <= (25/2+35.355), "25km",
                            ifelse(MaxDispersalDistanceKm <= (50/2+70.71), "50km",
                                   ifelse(MaxDispersalDistanceKm > (50/2+70.71), "110km",
                                          NA)))) %>%
  select(MaxDispersalDistanceKm, cv_res, everything()) %>%
  left_join(., cvs) %>%
  select(-median_cv_spp, -min_cv_spp, -max_cv_spp) %>% ## get rid of cv metrics we don't want
  rename("ClimVeloKmY_RelScale" = mean_cv_spp, 
         "sdClimVeloKmY_RelScale" = sd_cv_spp, 
         "q1ClimVeloKmY_RelScale" = q1_cv_spp, 
         "q3ClimVeloKmY_RelScale" = q3_cv_spp,
         "RelScale" = cv_res) %>%
  select(ClimVeloKmY_RelScale, sdClimVeloKmY_RelScale, q1ClimVeloKmY_RelScale, q3ClimVeloKmY_RelScale, 
         RelScale, everything())

## now add all other scales 
## pivot wider: 
cvs <- cvs %>%
  select(-median_cv_spp, -min_cv_spp, -max_cv_spp) %>%
  rename("ClimVeloKmY" = mean_cv_spp, "sdClimVeloKmY" = sd_cv_spp, 
         "q1ClimVeloKmY" = q1_cv_spp, "q3ClimVeloKmY" = q3_cv_spp) %>%
  pivot_wider(names_from = cv_res, 
              values_from = c(ClimVeloKmY, sdClimVeloKmY, q1ClimVeloKmY, q3ClimVeloKmY),
              names_glue = "{.value}_{cv_res}") %>%
  select(scientificName_checked, ID, Type, Eco,
         ClimVeloKmY_25km, sdClimVeloKmY_25km, q1ClimVeloKmY_25km, q3ClimVeloKmY_25km,
         ClimVeloKmY_50km, sdClimVeloKmY_50km, q1ClimVeloKmY_50km, q3ClimVeloKmY_50km,
         ClimVeloKmY_110km, sdClimVeloKmY_110km, q1ClimVeloKmY_110km, q3ClimVeloKmY_110km)

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
         ClimVeloKmY_RelScale, sdClimVeloKmY_RelScale, q1ClimVeloKmY_RelScale, q3ClimVeloKmY_RelScale, RelScale, 
         ClimVeloKmY_25km, sdClimVeloKmY_25km, q1ClimVeloKmY_25km, q3ClimVeloKmY_25km,
         ClimVeloKmY_50km, sdClimVeloKmY_50km, q1ClimVeloKmY_50km, q3ClimVeloKmY_50km,
         ClimVeloKmY_110km, sdClimVeloKmY_110km, q1ClimVeloKmY_110km, q3ClimVeloKmY_110km, everything())
## save 
write.csv(v3_merged, "data-processed/v3_with-cv.csv", row.names = FALSE)

## see how many have neg clim velo at leading edge 
v3_merged = read.csv("data-processed/v3_with-cv.csv")

