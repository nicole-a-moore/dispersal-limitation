## dig in to the data 
library(tidyverse)
library(sf)

## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to latitudinal 
dd <- filter(dd, Type == "LAT")

## filter 
dd <- dd %>%
  ## filter to cases where we expect dispersal to matter (range expansion): 
  ## leading edge shifts where climate velocity is positive
  ## trailing edge shifts where climate velocity is negative 
  ## centroid shifts
  filter(Param == "TE" & mean_cv_studylevel <= 0 |
           Param == "LE" & mean_cv_studylevel >= 0 |
           Param == "O") 

## split into plants and birds
plants <- filter(dd, group == "Plants")
birds <- filter(dd, group == "Birds")

## what types of plants are we dealing with? 
plants %>% 
  ggplot(aes(x = class)) + geom_bar()
## mostly Magnoliopsida (dicotyledons - flowering plants), small number of liliopsida (monocotyledons) 
## and pinopsida (gymnosperms)

## annual/biennial/perennial? 
plants %>% 
  ggplot(aes(x = YearOfMaturity, fill = class)) + geom_histogram()
## most mature after 1 year  

## which ones mature really late?
plants %>%
  filter(YearOfMaturity >= 10) %>%
  select(scientificName, YearOfMaturity, class.x) %>%
  View()
## the trees 

## what's the main mode of dispersal? 
tamme <- read.csv("data-processed/Tamme_harmonized.csv")

tamme_sub <- tamme %>% select(scientificName, Growth_form, Dispersal_syndrome, Seed_release_height_.m.)

plants_traits <- plants %>% 
  left_join(., tamme_sub, by = c("scientificName")) 

## GROWTH FORM
plants_traits %>% 
  ggplot(aes(x = Growth_form, fill = class)) + geom_bar()
## most are herbaceous plants or trees

## DISPERSAL SYNDROME
plants_traits %>% 
  ggplot(aes(x = Dispersal_syndrome, fill = class)) + geom_bar()
## most are wind dispersed with special seed adaptations or with no special adaptation
## then animal, then wind
## very few ant
plants_traits %>% 
  ggplot(aes(x = Dispersal_syndrome, y = ShiftKmY, colour = Growth_form)) + 
  geom_point() +
  facet_wrap(~Growth_form) 

## HEIGHT
plants_traits %>% 
  ggplot(aes(x = Seed_release_height_.m., fill = class)) + geom_histogram()

plants_traits %>%
  filter(Seed_release_height_.m. < 10) %>%
  ggplot(aes(x = Seed_release_height_.m., fill = class)) + geom_histogram()

## most plants are pretty short (release seeds from < 5m)

## are the flowering plants capable of selfing?
## get from TRY - trait 204, describes whether flowers are xenogamous or autogamous 
try = read_delim("data-raw/dispersal-data/TRY-query/31695_12022024015607/31695.txt") %>%
  filter(!is.na(TraitName))
unique(try$TraitName)

## subset to bioshifts species with dispersal distance
try_bs <- filter(try, SpeciesName %in% dd$scientificName)

try_bs %>%
  ggplot(aes(x = OrigValueStr)) +
  geom_bar()
## a lot are allogamous 

try_bs %>%
  filter(SpeciesName %in% try_bs$SpeciesName[duplicated(try_bs$SpeciesName)]) %>%
  View
## some have multiple reported mating types 



plants %>% 
  filter(class == "Magnoliopsida") %>% 
  ggplot(aes(x = order)) + geom_bar() +
  coord_flip()
## lots of asterales 
## caryophyllales

plants %>% 
  filter(class == "Pinopsida") %>% 
  ggplot(aes(x = family)) + geom_bar() +
  coord_flip()
## pine trees

plants %>% 
  filter(class == "Liliopsida") %>% 
  ggplot(aes(x = order)) + geom_bar() +
  coord_flip()
## most are poales (grasses, bromeliads, rushes, sedges)


plants %>%
 ggplot(aes(x = DispersalPotentialKmY)) + geom_histogram()

plants %>%
  ggplot(aes(x = ShiftKmY)) + geom_histogram()

plants %>%
  ggplot(aes(y = ShiftKmY, x = DispersalPotentialKmY, colour = ClimVeloTKmY_study)) + geom_point() +
  scale_x_log10() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 8))

plants %>%
  filter(Code != "MeanDispersalDistance") %>%
  filter(Param == "LE") %>%
  ggplot(aes(y = ShiftKmY*1000, x = DispersalPotentialKmY*1000, colour = Grain_size)) + geom_point() +
  scale_x_log10() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 4000)) +
  labs(x = "Potential dispersal rate (m/y)", y = "Range shift rate (m/y)")


## plants %>%
plants %>%
  mutate(diff = ShiftKmY - DispersalPotentialKmY) %>% 
  arrange(-diff) %>%
  select(scientificName, class, ShiftKmY, DispersalPotentialKmY, diff, YearOfMaturity,AgeAtMaturityDays,
         DispersalDistance, DispersalUnit,  
         DispersalPotentialKmY, Database, Code) %>%
  View


plants_traits %>%
  mutate(diff = ShiftKmY - DispersalPotentialKmY) %>% 
  ggplot(aes(x = diff*1000, fill = Dispersal_syndrome)) + geom_histogram() +
  facet_wrap(~Dispersal_syndrome)
## wind dispersed seeds with no special adaptations seem largely underestimated 
## weird 


birds %>%
  filter(Code != "MeanDispersalDistance") %>%
  filter(Param == "LE") %>%
  ggplot(aes(y = ShiftKmY, x = DispersalPotentialKmY, colour = Grain_size)) + geom_point() +
  scale_x_log10() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 8)) 




## what types of birds are we dealing with?
## site fidelity limitations? habitat limitations? 
birds %>% 
  ggplot(aes(x = order)) + geom_bar() +
  coord_flip()
## mostly passeriformes

## what types of dispersal studies are we dealing with?
plants %>%
  ggplot(aes(x = Code)) +
  geom_bar()
## mostly mean or max observations 
plants %>%
  ggplot(aes(x = Database)) +
  geom_bar()
## most observations are from Tamme 



## what types of range shift studies are we dealing with?

## plot area 
plants %>%
  ggplot(aes(x = Areakm2, fill = class)) + geom_histogram() 
## plot area versus shift distance 
plants %>%
  ggplot(aes(x = Areakm2, y = ShiftKmY*Duration, colour = class)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

## plot latitudinal extent versus shift distance 
plants %>%
  ggplot(aes(x = LatExtentk, y = ShiftKmY*Duration, colour = class)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

## plot latitudinal extent versus potential shift based on dispersal
plants %>%
  ggplot(aes(x = LatExtentk, y = DispersalPotentialKmY*Duration, colour = class)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)


plants %>%
  ggplot(aes(x = Duration)) + geom_histogram() 

plants %>%
  ggplot(aes(x = ShiftKmY)) + geom_histogram() 
## all are less than km/y
plants %>%
  ggplot(aes(x = ShiftKmY*1000)) + geom_histogram() # in m/y

min(plants$ShiftKmY*1000)
max(plants$ShiftKmY*1000)


#######
## garbage 
#######

## read in study area shapefiles 
layers <- st_layers("data-raw/BIOSHIFTSv1/bioshifts-download/Bioshifts/Study_Areas.gdb")

sf_use_s2(FALSE)

polys <- c()
names = c()
i=1
while(i <= length(layers$name)) {
  
  curr = st_read("data-raw/BIOSHIFTSv1/bioshifts-download/Bioshifts/Study_Areas.gdb", layer = layers$name[i]) %>%
    select(Shape)
  
  names = append(names, rep(layers$name[i], nrow(curr)))
  
  polys <- rbind(polys, curr)
  
  i=i+1
}

polys$ID = names

## filter to only plants we have data for
polys <- filter(polys, ID %in% dd$Source)

polys_plants <- filter(polys, ID %in% plants$Source)

## plot a study area
library(rnaturalearth)
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')

## crop map to species range
## plot study area on top
curr = polys %>%
  filter(ID == "A116_P1")

bbox = st_bbox(curr$Shape)

cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10,
                                   bbox$xmax + 10, bbox$ymax + 10))

curr %>%
  ggplot(data = ., aes(geometry = Shape)) +
  geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) +
  geom_sf(fill = "red") +
  theme_bw() +
  theme(panel.grid = element_blank())







tamme <- read_csv("data-raw/dispersal-data/Tamme/Tamme_DispersalDistanceData.csv")
colnames(tamme) <- str_replace_all(colnames(tamme), "\\ ", "_")
length(unique(tamme$Species)) #576 spp

tamme %>%
  ggplot(aes(x = `Maximum_recorded_dispersal_distance_(m)`, fill = Data_type)) +
  geom_histogram() +
  scale_x_log10()
## field data tend to be greater than modelled data



