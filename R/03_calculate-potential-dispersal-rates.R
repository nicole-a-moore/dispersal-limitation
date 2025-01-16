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
dscale = read.csv("data-processed/dispersal-distance-collated.csv") 

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
dscale <- dscale %>%
  group_by(scientificName_checked) %>%
  mutate(MaxDispersalDistanceKm = max(DispersalDistanceKm),
         MedianDispersalDistanceKm = median(DispersalDistanceKm),
         NObs_MedianDispersalDistanceKm = n()) %>% 
  ungroup() %>%
  filter(DispersalDistanceKm == MaxDispersalDistanceKm)

## look at duplicates
dups = dscale %>%
  group_by(scientificName_checked) %>%
  filter(length(scientificName_checked) >= 2) %>% 
  mutate(Field = ifelse(Source == "Chu 2021", "geom/arimean", Field),
         ObservationTypeSpecific = ifelse(Source == "Chu 2021", "geometric/arithmetic mean natal dispersal distance",
                                          ObservationTypeSpecific)) %>%
  filter(Source != "unpub.") %>%
  filter(Database != "Vittoz & Engler 2007") %>% ## get rid of vittoz duplicates 
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
v3 = read.csv("data-processed/v3_shifts.csv")

## subset v3 to only species with dispersal distance  
v3 <- filter(v3, scientificName_checked %in% dscale$scientificName_checked)
length(unique(v3$scientificName_checked)) #602 species 


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
length(unique(v3$scientificName_checked[which(!is.na(v3$MaxDispersalDistanceKm))])) # 602 species have max dispersal distance
length(which(is.na(v3$MaxDispersalDistanceKm))) # 0 missing max dispersal distance


###################################################
####   calculating annual dispersal potential  ####
###################################################
### join age at maturity data with dispersal data 
am <- read.csv("data-processed/age-at-maturity.csv")

## plot
am %>% 
  ggplot(aes(x = AgeAtMaturity, fill = class)) + geom_histogram() +
  theme_bw() +
  scale_x_log10()

## if multiple estimates of age at maturity per species, keep the lowest 
am_join <- am %>%
  group_by(scientificName_checked) %>%
  mutate(AgeAtMaturity = as.numeric(as.character(AgeAtMaturity))) %>%
  mutate(AgeAtMaturityDays = ifelse(Unit == "yrs", 
                                    AgeAtMaturity*365,
                                    ifelse(Unit == "weeks",
                                           AgeAtMaturity*7,
                                           AgeAtMaturity))) %>% # convert all to days 
  mutate(AgeAtMaturityDays = min(AgeAtMaturityDays)) %>% # select minimum per species 
  ungroup() %>%
  select(scientificName_checked, AgeAtMaturityDays) %>%
  unique() %>%
  mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 

## join to dispersal data:
v3 <- left_join(v3, am_join, by = "scientificName_checked")

length(unique(v3$scientificName_checked)) ## still have all the species
length(unique(v3$scientificName_checked[which(is.na(v3$AgeAtMaturityDays))])) 
## 136 / 602 species with dispersal estimates do not have age at maturity data 

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

##########################
####   plot the data  ####
##########################
## colour pal
mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

nrow(v3) # 3957 range shifts
length(unique(v3$scientificName_checked)) # 466 species 
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

v3 <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter out elevational studies 
lat <- v3  %>%
  filter(Type != "ELE")
  
nrow(lat) # 2764 shifts
length(unique(lat$scientificName_checked)) # 401 species 
  
## how many are centroid vs. leading edge 
lat %>% 
  group_by(Param) %>%
  tally()

## plot distribution of dispersal scale for proposal document:
lat %>%
  select(DispersalPotentialKmY, scientificName_checked, group) %>%
  unique() %>%
  ggplot(aes(x = DispersalPotentialKmY, fill = group)) + geom_histogram() +
  theme_bw() +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.0001", "0.001", "0.01", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(fill = "", x = "Potential dispersal rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank())

ggsave(path = "figures/proposal", filename = "disp-pot-distribution.png",
       height = 3, width = 6)

lat %>%
  select(DispersalPotentialKmY, scientificName_checked, group) %>%
  unique() %>%
  ggplot(aes(x = DispersalPotentialKmY)) + geom_histogram() +
  theme_bw() +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(fill = "", x = "Potential dispersal rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank()) 

lat %>%
  ggplot(aes(x = DispersalPotentialKmY, fill = group)) + geom_histogram() +
  theme_bw() +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.0001", "0.001", "0.01", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(fill = "", x = "Potential dispersal rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Param) +
  theme(legend.position = "none")

ggsave(path = "figures/proposal", filename = "disp-pot-distribution-by-position.png",
       height = 3, width = 6)





