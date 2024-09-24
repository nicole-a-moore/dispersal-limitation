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
  group_by(scientificName) %>%
  summarise(sd_dp = sd(DispersalDistanceKm), class = unique(class)) %>%
  filter(!is.na(sd_dp))

sd %>%
  ggplot(aes(x = sd_dp, fill = class)) + geom_histogram() + 
  scale_x_log10() 
## birds have highest standard deviation in dispersal potential within species 
max(sd$sd_dp)


## choose one dispersal estimate per species
## if multiple estimates of dispersal, keep the highest 
## test for sensitivity later 
dscale <- dscale %>%
  group_by(scientificName) %>%
  mutate(MaxDispersalDistanceKm = max(DispersalDistanceKm)) %>% 
  ungroup() %>%
  filter(DispersalDistanceKm == MaxDispersalDistanceKm)

## look at duplicates
dups = dscale %>%
  group_by(scientificName) %>%
  filter(length(scientificName) >=2) %>% 
  mutate(Field = ifelse(Source == "Chu 2021", "geom/arimean", Field),
         ObservationTypeSpecific = ifelse(Source == "Chu 2021", "geometric/arithmetic mean natal dispersal distance",
                                          ObservationTypeSpecific)) %>%
  ungroup() %>%
  filter(Source != "unpub.") %>%
  distinct()  

leftover = dups %>% group_by(scientificName) %>%
  filter(length(scientificName) >=2) %>% distinct() 
leftover = leftover[1,]

dups <- filter(dups, !scientificName %in% leftover$scientificName)

dscale <- dscale %>%
  group_by(scientificName) %>%
  filter(length(scientificName) <2) %>%
  rbind(., dups) %>%
  rbind(., leftover)


###########################################
####   prep range expansion rate data  ####
###########################################
## read clean version of bioshifts v3 with climate velocities
v3 = read.csv("data-processed/v3_with-cv.csv")

## subset v3 to only species with dispersal distance  
v3 <- filter(v3, scientificName %in% dscale$scientificName)
length(unique(v3$scientificName)) #581 species 


#########################################################
####   join range shift and dispersal distance data  ####
#########################################################
#----------------------
## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "db", "db_code", "kingdom",
                            "phylum", "class", "family")) %>%
  unique()

## rename some columns that have the same names 
dscale <- rename(dscale, "DispersalSource"= Source, "DispersalUnit"= Unit)

## join
v3 = left_join(v3, dscale) 

## check on merge
length(unique(v3$scientificName[which(!is.na(v3$DispersalDistanceKm))])) # 581 species have dispersal distance
length(which(is.na(v3$DispersalDistanceKm))) # 0 missing dispersal distance


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
  group_by(scientificName) %>%
  mutate(AgeAtMaturity = as.numeric(as.character(AgeAtMaturity))) %>%
  mutate(AgeAtMaturityDays = ifelse(Unit == "yrs", 
                                    AgeAtMaturity*365,
                                    ifelse(Unit == "weeks",
                                           AgeAtMaturity*7,
                                           AgeAtMaturity))) %>% # convert all to days 
  mutate(AgeAtMaturityDays = min(AgeAtMaturityDays)) %>% # select minimum per species 
  ungroup() %>%
  select(scientificName, AgeAtMaturityDays) %>%
  unique() %>%
  mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 

## join to dispersal data:
v3 <- left_join(v3, am_join, by = "scientificName")

length(unique(v3$scientificName)) ## still have all the species
length(unique(v3$scientificName[which(is.na(v3$AgeAtMaturityDays))])) 
## 132 / 581 species with dispersal estimates do not have age at maturity data 

unique(v3$scientificName[which(is.na(v3$AgeAtMaturityDays))])

## filter to only species with age at maturity 
v3 <- filter(v3, !is.na(AgeAtMaturityDays))

## calculate dispersal potential for species with age at maturity 
v3 = v3 %>%
  mutate(DispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                        DispersalDistanceKm/YearOfMaturity,
                                        NA)) %>%
  mutate(DispersalPotentialmY = ifelse(!is.na(YearOfMaturity), 
                                       (DispersalDistanceKm*1000)/YearOfMaturity,
                                       NA)) %>%
  mutate(DispersalDistancem = DispersalDistanceKm*1000)

## make new vars 
v3 <- v3 %>%
  ## make variable for whether shift is in line with climate velocity 
  mutate(tracking_climate = sign(Rate) == sign(mean_cv_studylevel)) %>%
  ## make variable for whether shift is a contraction
  mutate(is_contraction = ifelse(Param == "LE" & Rate >= 0, "NO",
                                 ifelse(Param == "TE" & Rate <= 0, "NO",
                                        ifelse(Param == "O", "UNKNOWN", "YES")))) %>%
  ## make one column for annual dispersal potential
  mutate(annual_dispersal_pot = ifelse(Type == "ELE",
                                       DispersalPotentialmY,
                                       ifelse(Type == "LAT",
                                              DispersalPotentialKmY,
                                              NA))) 

## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
v3$ShiftKmY <- ifelse(v3$Type == "ELE", v3$Rate / 1000, 
                               v3$Rate)
v3$ClimVeloTKmY_study <- ifelse(v3$Type == "ELE", v3$mean_cv_studylevel / 1000,
                                   v3$mean_cv_studylevel)
v3$ClimVeloTKmY_spp <- ifelse(v3$Type == "ELE", v3$mean_cv_sppspecific / 1000,
                                       v3$mean_cv_sppspecific)
v3$AnnualDispPotKmY <- ifelse(v3$Type == "ELE", v3$annual_dispersal_pot / 1000,
                                       v3$annual_dispersal_pot)


##########################
####   plot the data  ####
##########################
## colour pal
mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

nrow(v3) # 3823 range shifts
length(unique(v3$scientificName)) # 449 species 
unique(v3$Param)

grad <- v3 %>%
  ggplot(aes(x = Type, fill = Type)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") +
  scale_x_discrete(labels = c("Elevation studies", "Latitudinal studies")) +
  scale_fill_manual(values = c("#28587B", "#9FB798")) 

grad

ggsave(grad, path = "figures/sotm", filename = "barplot-gradients.png", 
       device = "png", height = 2, width = 4)

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

v3 %>%
  ggplot(aes(x = mean_cv_studylevel, y = mean_cv_sppspecific)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

## how often are the signs of study versus species-level climate velocity different?
v3 %>%
  filter(!is.na(mean_cv_sppspecific)) %>%
  mutate(same_sign = ifelse(sign(mean_cv_studylevel) == sign(mean_cv_sppspecific), "same", "different")) %>%
  ggplot(aes(same_sign)) +
  geom_bar()
## not often at all

## 
v3 %>%
  ggplot(aes(x = ClimVeloTKmY_study, y = ShiftKmY, colour = is_contraction)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(Type~Param)

## save 
write.csv(v3, "data-processed/v3_potential-dispersal-rate.csv", row.names = FALSE)

v3 <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter out elevational studies 
lat <- v3  %>%
  filter(Type != "ELE")
  
nrow(lat) # 2722 shifts
length(unique(lat$scientificName)) # 391 species 
  
## how many are centroid vs. leading edge 
lat %>% 
  group_by(Param) %>%
  tally()

## plot distribution of dispersal scale for proposal document:
lat %>%
  select(DispersalPotentialKmY, scientificName, group) %>%
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
  select(DispersalPotentialKmY, scientificName, group) %>%
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





