## calculate potential dispersal rates by dividing dispersal distance by reproductive frequency (age at maturity / generation time)
## combine with range expansion data 
library(tidyverse)
source("R/taxonomic-harmonization/clean_taxa_functions.R")
## read function to harmonize taxonomy
source("R/taxonomic-harmonization/harmonize.R")

#########################################
####   prep dispersal distance data  ####
#########################################
#----------------------
## read in dispersal distance data 
dscale = read.csv("data-processed/intermediate_files/dispersal-distance/dispersal-distance-collated.csv") 

## make sure it's clean
unique(dscale$Unit)
unique(dscale$DispersalDistance)
length(which(is.na(dscale$DispersalDistance)))

## plot distribution of distances 
dscale %>%
  ggplot(aes(x = DispersalDistanceKm, fill = class)) + geom_histogram() +
  theme_bw() +
  scale_x_log10()

dscale %>%
  ggplot(aes(x = Code, fill = class)) + geom_bar() +
  theme_bw() + coord_flip()

## how variable are dispersal potential estimates within species?
sd <- dscale %>%
  group_by(scientificName_checked) %>%
  summarise(sd_dp = sd(DispersalDistanceKm), class = unique(class)) %>%
  filter(!is.na(sd_dp))

sd %>%
  ggplot(aes(x = sd_dp, fill = class)) + geom_histogram() + 
  scale_x_log10() 
## birds have highest standard deviation in dispersal potential within species 
max(sd$sd_dp)


## calculate the max and median dispersal distance per species 
## keep rows that correspond to original data for the max value 
dscale_temp <- dscale %>%
  group_by(scientificName_checked) %>%
  mutate(MaxDispersalDistanceKm = max(DispersalDistanceKm),
         MedianDispersalDistanceKm = median(DispersalDistanceKm),
         NObs_MedianDispersalDistanceKm = n()) %>% 
  ungroup() 

dscale <- dscale_temp %>%
  filter(DispersalDistanceKm == MaxDispersalDistanceKm)

## look at duplicates
dups = dscale %>%
  group_by(scientificName_checked) %>%
  filter(length(scientificName_checked) >= 2) %>% 
  mutate(Field = ifelse(Source == "Chu 2021", "geom/arimean", Field),
         ObservationTypeSpecific = ifelse(Source == "Chu 2021", "geometric/arithmetic mean natal dispersal distance",
                                          ObservationTypeSpecific)) %>%
  filter(Source != "unpub.") %>%
  filter(!(length(scientificName_checked) >= 2 & Database == "Vittoz & Engler 2007")) %>% ## get rid of vittoz duplicates 
  ungroup() %>%
  distinct()

leftover = dups %>% group_by(scientificName_checked) %>%
  filter(length(scientificName_checked) >=2) %>% distinct() 

leftover <- filter(leftover, !str_detect(reported_name, "subsp"))
leftover = leftover[1:2,]

dups <- filter(dups, !scientificName_checked %in% leftover$scientificName_checked)

dscale <- dscale %>%
  group_by(scientificName_checked) %>%
  filter(length(scientificName_checked) < 2) %>%
  rbind(., dups) %>%
  rbind(., leftover) %>%
  select(-DispersalDistanceKm)

## check that each species has 1 obs
length(unique(dscale$scientificName_checked)) == nrow(dscale)

###########################################
####   prep range expansion rate data  ####
###########################################
## read clean version of bioshifts v3 
v3 = read.csv("data-processed/intermediate_files/bioshifts/v3_shifts.csv")

## subset v3 to only species with dispersal distance  
v3 <- filter(v3, scientificName_checked %in% dscale$scientificName_checked)
length(unique(v3$scientificName_checked)) #601 species 


#########################################################
####   join range shift and dispersal distance data  ####
#########################################################
#----------------------
## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "scientificName", "db", "db_code", "kingdom",
                            "phylum", "class", "family")) %>%
  unique()

## rename some columns that have the same names 
dscale <- rename(dscale, "DispersalSource"= Source, "DispersalUnit"= Unit)

## join
v3 = left_join(v3, dscale, by = "scientificName_checked") 

## check on merge
length(unique(v3$scientificName_checked[which(!is.na(v3$MaxDispersalDistanceKm))])) # 601 species have max dispersal distance
length(which(is.na(v3$MaxDispersalDistanceKm))) # 0 missing max dispersal distance


###################################################
####   calculating annual dispersal potential  ####
###################################################
### join age at maturity data with dispersal data 
am <- read.csv("data-processed/intermediate_files/dispersal-frequency/age-at-maturity.csv")

## plot
am %>% 
  ggplot(aes(x = AgeAtMaturity, fill = class)) + geom_histogram() +
  theme_bw() +
  scale_x_log10()

## if multiple estimates of age at maturity per species, keep the lowest 
am_join <- am %>%
  mutate(AgeAtMaturity = as.numeric(as.character(AgeAtMaturity))) %>%
  mutate(AgeAtMaturityDays = ifelse(Unit == "yrs", 
                                    AgeAtMaturity*365,
                                    ifelse(Unit == "weeks",
                                           AgeAtMaturity*7,
                                           AgeAtMaturity))) %>% # convert all to days 
  group_by(scientificName_checked) %>%
  filter(AgeAtMaturityDays == min(AgeAtMaturityDays)) %>% # select minimum per species 
  slice(1) %>% ## keep first row 
  ungroup() %>%
  select(scientificName_checked, AgeAtMaturityDays, AgeAtMaturity, Database, Unit, Code) %>%
  unique() %>%
  rename("AgeAtMaturyCode" = Code, "AgeAtMaturityDatabase" = Database) %>%
  mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 
# 
# am_join = am_join_temp %>%
#   select(scientificName_checked, AgeAtMaturityDays) %>%
#   unique() %>%
#   rename("AgeAtMaturyCode" = Code, "AgeAtMaturitySource" = Source) %>%
#   mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 

## join to dispersal data:
v3 <- left_join(v3, am_join, by = "scientificName_checked")

length(unique(v3$scientificName_checked)) ## still have all the species
length(unique(v3$scientificName_checked[which(is.na(v3$AgeAtMaturityDays))])) 
## 138 / 601 species with dispersal estimates do not have age at maturity data 

unique(v3$scientificName[which(is.na(v3$AgeAtMaturityDays))])

## filter to only species with age at maturity 
v3 <- filter(v3, !is.na(AgeAtMaturityDays))

## calculate dispersal potential for species with age at maturity 
v3 = v3 %>%
  mutate(DispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                        MaxDispersalDistanceKm/YearOfMaturity,
                                        NA)) %>%
  mutate(MedianDispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                        MedianDispersalDistanceKm/YearOfMaturity,
                                        NA)) 


## convert units of latitudinal and elevation shifts to km/y:
v3$ShiftKmY <- ifelse(v3$Type == "ELE", v3$Rate / 1000, 
                               v3$Rate)

## make database that has species, dispersal distance measurements + source, longevity measurement + source, Max and Median dispersal measurement

db <- v3 %>%
  select(sp_name_checked, Database, 
         MaxDispersalDistanceKm, MedianDispersalDistanceKm, NObs_MedianDispersalDistanceKm,
         AgeAtMaturity, Unit, AgeAtMaturyCode, AgeAtMaturityDatabase, YearOfMaturity, DispersalPotentialKmY) %>% 
  distinct() %>%
  rename("MaxDispersalDistance_Database" = Database, "DispersalFrequency_MeasurementType" = AgeAtMaturyCode,
         "DispersalFrequency_Database" = AgeAtMaturityDatabase, "DispersalFrequency" = AgeAtMaturity, 
         "DispersalFrequency_Unit" = Unit, "DispersalFrequencyY" = YearOfMaturity,
         "PotentialDispersalRateKmY" = DispersalPotentialKmY) %>%
  mutate(sp_name_checked = str_replace_all(sp_name_checked, "_", " ")) 

## add in raw dispersal estimates (before mean and median were calculated)
db = dscale_temp %>%
  select(scientificName_checked, DispersalDistance, Unit, Code, Database, ObservationTypeGeneral) %>%
  rename("sp_name_checked" = scientificName_checked, "DispersalDistance_ObservationType" = ObservationTypeGeneral,
         "DispersalDistance_Database" = Database,
         "DispersalDistance_MeasurementType" = Code, "DispersalDistance_Unit" = Unit) %>%
  left_join(., db) %>%
  filter(!is.na(DispersalFrequency)) %>%
  select(sp_name_checked, DispersalDistance, DispersalDistance_Unit, DispersalDistance_ObservationType,
         DispersalDistance_MeasurementType, DispersalDistance_Database,
         DispersalFrequency, DispersalFrequency_Unit, DispersalFrequency_MeasurementType,
         DispersalFrequency_Database, MaxDispersalDistanceKm, MaxDispersalDistance_Database, 
         MedianDispersalDistanceKm, NObs_MedianDispersalDistanceKm, 
         DispersalFrequencyY, PotentialDispersalRateKmY) %>%
  arrange(sp_name_checked) %>%
  rename("GenusSpecies" = sp_name_checked)

length(which(db$DispersalFrequency_MeasurementType %in% c("MaturityFromLifespan", "MaturityFromGrowthForm")))
## 24 species 

## save the database as a csv:
write.csv(db, "figures/databaseS1_dispersal-data.csv", row.names = F)


##########################
####   plot the data  ####
##########################
## colour pal
mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

nrow(v3) # 3952 range shifts
length(unique(v3$scientificName_checked)) # 465 species 
unique(v3$Param)

v3 %>%
  group_by(Type) %>% tally()

groups <- v3 %>%
  ggplot(aes(x = group)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") 

groups

ggsave(groups, path = "figures/sotm", filename = "barplot-groups.png", 
       device = "png", height = 2, width = 4)

v3 %>%
  group_by(group) %>% tally()

## save 
write.csv(v3, "data-processed/v3_potential-dispersal-rate.csv", row.names = FALSE)
  