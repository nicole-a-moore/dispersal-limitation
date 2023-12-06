## calculate potential dispersal rates by dividing dispersal distance by reproductive frequency (age at maturity / generation time)
## combine with range expansion data 
library(tidyverse)
source("R/taxonomic-harmonization/clean_taxa_functions.R")

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
## read clean version of bioshifts 
v1 = read.csv("data-processed/BIOSHIFTSv1_harmonized.csv")

## get start and end / duration data from other version of bioshifts 
other_rs_data <- read.delim('data-raw/BIOSHIFTSv1/Shifts2018_checkedtaxo.txt')

other_rs_data <- select(other_rs_data, c(New_name, ID, START, END, DUR)) %>%
  distinct()

## if there are multiple durations per study, choose shortest 
other_rs_data = other_rs_data %>%
  group_by(ID, New_name) %>%
  mutate(DUR = min(DUR), 
         START = min(START),
         END = min(END)) %>%
  ungroup() %>%
  distinct() %>%
  rename("Source" = ID) %>%
  mutate(New_name = str_replace_all(New_name, "\\_", " "))

## clean names in the orginal database 
other_rs_data$reported_name <- other_rs_data$New_name
other_rs_data$New_name = Clean_Names(other_rs_data$New_name, return_gen_sps = F)

##harmonize the new name
bs_harm <- harmonize(other_rs_data$New_name)
notfound <- filter(bs_harm, is.na(db_code))

## add back species that weren't found but are resolved to species level
notfound <- filter(notfound, !str_detect(notfound$species, "sp."))
notfound$scientificName = notfound$species
bs_harm <- rbind(bs_harm, notfound) %>%
  distinct()

other_rs_data <- other_rs_data %>%
  rename("reported_name_fixed" = New_name)

other_rs_data_corrected <- left_join(other_rs_data, bs_harm, by = c("reported_name_fixed" = "species"),
                          relationship = "many-to-many") %>%
  unique() %>%
  filter(!is.na(scientificName))

v1$scientificName[which(!v1$scientificName %in% other_rs_data_corrected$scientificName)]
## all are there 

## left join with original database 
v1 = left_join(v1,  other_rs_data_corrected)

## subset v1 to only species with dispersal distance  
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$scientificName)) #586 species 


#########################################################
####   join range shift and dispersal distance data  ####
#########################################################
#----------------------
## add dispersal distance 
## get rid of old taxonomy columns from v1 (they aren't right)
v1 <- select(v1, -c("Kingdom", "Phylum", "Class", "Order", "Family"))
v1_saved = v1

## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "db", "db_code")) %>%
  unique()

## rename some columns that have the same names 
dscale <- rename(dscale, "DispersalSource"= Source, "DispersalUnit"= Unit)

## join
v1 = left_join(v1, dscale, by = "scientificName") 

## check on merge
length(unique(v1$scientificName[which(!is.na(v1$DispersalDistanceKm))])) # 586 species have dispersal scale
length(which(is.na(v1$DispersalDistanceKm))) # 0 missing dispersal scale



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
v1 <- left_join(v1, am_join, by = "scientificName")

length(unique(v1$scientificName)) ## still have all the species
length(unique(v1$scientificName[which(is.na(v1$AgeAtMaturityDays))])) 
## 137 / 586 species with dispersal estimates do not have age at maturity data 

unique(v1$scientificName[which(is.na(v1$AgeAtMaturityDays))])

## filter to only species with age at maturity 
v1 <- filter(v1, !is.na(AgeAtMaturityDays))

## calculate dispersal potential for species with age at maturity 
v1 = v1 %>%
  mutate(DispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                        DispersalDistanceKm/YearOfMaturity,
                                        NA)) %>%
  mutate(DispersalPotentialmY = ifelse(!is.na(YearOfMaturity), 
                                       (DispersalDistanceKm*1000)/YearOfMaturity,
                                       NA)) %>%
  mutate(DispersalDistancem = DispersalDistanceKm*1000)


## filter data 
v1_filtered <- v1 %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>% 
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make one column for annual dispersal potential, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       DispersalPotentialmY,
                                       ifelse(Gradient == "Latitudinal",
                                              DispersalPotentialKmY,
                                              NA))) %>%
  mutate(climate_velocity = ifelse(Gradient == "Elevation",
                                   EleVeloT,
                                   ifelse(Gradient == "Latitudinal",
                                          LatVeloT,
                                          NA)))

## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
v1_filtered$ShiftKmY <- ifelse(v1_filtered$Gradient == "Elevation", v1_filtered$ShiftR / 1000, 
                               v1_filtered$ShiftR)
v1_filtered$ClimVeloTKmY <- ifelse(v1_filtered$Gradient == "Elevation", v1_filtered$EleVeloT / 1000,
                                   v1_filtered$LatVeloT)
v1_filtered$AnnualDispPotKmY <- ifelse(v1_filtered$Gradient == "Elevation", v1_filtered$annual_dispersal_pot / 1000,
                                       v1_filtered$annual_dispersal_pot)


##########################
####   plot the data  ####
##########################
## colour pal
mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

nrow(v1_filtered) # 1938 range shifts
length(unique(v1_filtered$scientificName)) # 411 species 
unique(v1_filtered$Position)

grad <- v1_filtered %>%
  ggplot(aes(x = Gradient, fill = Gradient)) +
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

v1_filtered %>%
  group_by(Gradient) %>% tally()

groups <- v1_filtered %>%
  ggplot(aes(x = group)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") 

groups

ggsave(groups, path = "figures/sotm", filename = "barplot-groups.png", 
       device = "png", height = 2, width = 4)

v1_filtered %>%
  group_by(group) %>% tally()

## save 
write.csv(v1_filtered, "data-processed/v1_potential-dispersal-rate.csv", row.names = FALSE)

