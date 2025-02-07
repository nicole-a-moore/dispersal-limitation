## brings together published data on empirical dispersal distance that are representative of natal/seed dispersal
library(readr)
library(tidyverse)
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

## read function to harmonize taxonomy
source("R/taxonomic-harmonization/harmonize.R")


## get list of all species in bioshifts v3 
#################################################
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

sp <- select(v3, scientificName, scientificName_checked) %>%
  distinct()
length(unique(sp$scientificName)) ## 12617 spp

#----------------------------
# Collate dispersal data 
#----------------------------
## data from:
## Paradis et al 2002 - mean bird natal/breeding dispersal distances 
## Sutherland - mean/median natal dispersal distances
## Whitmee & Orme 2013 - mammal natal breeding dispersal distances 
## Jenkins - bird and plant dispersal 
## Flores - plant seed dispersal
## Tamme - plant seed dispersal 
## TRY - seed or whole plant dispersal
## Chu - bird natal dispersal distances 
## Vittoz and Engler 2007 - seed dispersal distances 

## write vector of columns we want to be in final database
cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")

#---------------------
# Paradis et al 2002
#---------------------
## read in data
par = read_csv("data-raw/dispersal-data/Paradis/Paradis_et_al_2002.csv") 
colnames(par) <- str_replace_all(colnames(par), "\\ ", "_")
length(unique(par$Species)) #75 spp

## prep data so function doesn't throw error
sp_names = par$Species
tofind <- data.frame(matrix(nrow = length(sp_names), ncol = 8))
names(tofind) = c("scientificName", "kingdom", "phylum", "class", "order", "family", "db", "db_code")

tofind <- data.frame(species = sp_names, tofind)

tofind[,1:ncol(tofind)] = lapply(tofind[,1:ncol(tofind)], as.character) 

togo <- tofind[which(is.na(tofind$scientificName)),]

## harmonize reported names
par_harm <- harmonize(par$Species)
notfound <- filter(par_harm, is.na(db_code)) ## 0 species not found

## rename columns 
par <- par %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

par <- left_join(par, par_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## reorganize the data 
par_dd <- par %>%
  select(all_of(cols_to_keep), 
         AM_breeding, AM_natal, GM_breeding, GM_natal) %>%
  gather(key = "Field", value = "DispersalDistance", c(AM_breeding, AM_natal, GM_breeding, GM_natal)) %>%
  mutate(Field = ifelse(Field == "AM_breeding", "ArithmeticMeanBreedingDispersal",
                        ifelse(Field == "AM_natal", 
                               "ArithmeticMeanNatalDispersal",
                               ifelse(Field == "GM_natal", 
                                      "GeometricMeanNatalDispersal",
                                      ifelse(Field == "GM_breeding", 
                                             "GeometricMeanBreedingDispersal", 
                                             NA))))) %>%
  mutate(Code = "MeanDispersalDistance") %>%
  mutate(ObservationTypeSpecific = ifelse(Field %in% c("ArithmeticMeanNatalDispersal", 
                                                       "GeometricMeanNatalDispersal"), "natal dispersal", 
                                          "breeding dispersal")) %>%
  filter(!is.na(DispersalDistance)) %>%
  mutate(Sex = NA, Source = NA, Unit = "km",
         Database = "Paradis et al. 2002") 

## check how many species in bioshifts 
length(which(unique(par_dd$scientificName) %in% unique(sp$scientificName))) ## 74 :-)
par_sp <- unique(par_dd$scientificName)[which(unique(par_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# Sutherland 2000
#---------------------
## birds and mammals natal and breeding dispersal
suth_mamm <- read_csv("data-raw/dispersal-data/Sutherland/Sutherland_2000_mammals.csv")
colnames(suth_mamm) <- str_replace_all(colnames(suth_mamm), "\\ ", "_")
suth_mamm$Species = ifelse(suth_mamm$Species == "", NA, suth_mamm$Species)

## fill
suth_mamm = fill(suth_mamm, Species, .direction = "down")

sm_harm <- harmonize(suth_mamm$Species)

notfound <- filter(sm_harm, is.na(db_code))

## rename columns 
suth_mamm <- suth_mamm %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

suth_mamm <- left_join(suth_mamm, sm_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Sutherland mammals:
suth_mamm_dd <- suth_mamm %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type, Source) %>%
  mutate(ObservationTypeSpecific = "natal dispersal") %>%
  gather(key = "Field", value = "DispersalDistance", c(`Natal_dispersal_median_distance_(km)`,
                                                       `Natal_dispersal_maximum_distance_(km)`)) %>%
  mutate(Code = ifelse(Field == "Natal_dispersal_median_distance_(km)",
                       "MedianDispersalDistance", 
                       ifelse(Field == "Natal_dispersal_maximum_distance_(km)",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Database = "Sutherland (mammals)", Unit = "km") %>%
  filter(!is.na(DispersalDistance), DispersalDistance != "...") %>%
  mutate(Sex = str_split_fixed(DispersalDistance, " ", 2)[,2], 
         DispersalDistance = str_split_fixed(DispersalDistance, " ", 2)[,1]) %>%
  select(-Obs_type) %>%
  mutate()

## check how many species in bioshifts 
length(which(unique(suth_mamm_dd$scientificName) %in% unique(sp$scientificName))) ## 11
suth_mamm_sp <- unique(suth_mamm_dd$scientificName)[which(unique(suth_mamm_dd$scientificName) %in% unique(sp$scientificName))]


suth_bird <- read_csv("data-raw/dispersal-data/Sutherland/Sutherland_2000_birds.csv")
colnames(suth_bird) <- str_replace_all(colnames(suth_bird), "\\ ", "_")
length(unique(suth_bird$Species)) #78 spp
suth_bird$Species = ifelse(suth_bird$Species == "", NA, suth_bird$Species)

## fill so each row has a species name
suth_bird = fill(suth_bird, Species, .direction = "down")

## fix taxonomy manually
suth_bird$reported_name_fixed = suth_bird$Species
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Accipiter cooperi")] <- "Accipiter cooperii"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Bonasa bonasia")] <- "Tetrastes bonasia"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Dendragapus canadensis")] <- "Falcipennis canadensis"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Actitis macularia")] <- "Actitis macularius"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Speotyto cunicularia")] <- "Athene cunicularia"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Delichon urbica")] <- "Delichon urbicum"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Parus atricapillus")] <- "Poecile atricapillus" 
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Parus palustris")] <- "Poecile palustris"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Sitta europea")] <- "Sitta europaea"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Carduelis chloris")] <- "Chloris chloris"

## harmonize 
sb_harm <- harmonize(suth_bird$reported_name_fixed)
notfound <- filter(sb_harm, is.na(db_code)) ## 0 species not found

## rename columns 
suth_bird <- suth_bird %>%
  rename("reported_name" = Species) 

suth_bird <- left_join(suth_bird, sb_harm, by = c("reported_name_fixed" = "species"),
                       relationship = "many-to-many") %>%
  unique()

## Sutherland birds:
suth_bird_dd <- suth_bird %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type,
         Source) %>%
  filter(Source != "Paradis et al. (1998)") %>% ## get rid of estimates that are from Paradis 
  mutate(ObservationTypeSpecific = "natal dispersal") %>%
  gather(key = "Field", value = "DispersalDistance", c(`Natal_dispersal_median_distance_(km)`,
                                                       `Natal_dispersal_maximum_distance_(km)`)) %>%
  mutate(Code = ifelse(Field == "Natal_dispersal_median_distance_(km)",
                       "MedianDispersalDistance", 
                       ifelse(Field == "Natal_dispersal_maximum_distance_(km)",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Database = "Sutherland (birds)", Unit = "km") %>%
  filter(!is.na(DispersalDistance), DispersalDistance != "...") %>%
  mutate(Sex = str_split_fixed(DispersalDistance, " ", 2)[,2], 
         DispersalDistance = str_split_fixed(DispersalDistance, " ", 2)[,1]) %>%
  select(-Obs_type)

## check how many species in bioshifts 
length(which(unique(suth_bird_dd$scientificName) %in% unique(sp$scientificName))) ## 68
suth_bird_sp <- unique(suth_bird_dd$scientificName)[which(unique(suth_bird_dd$scientificName) %in% unique(sp$scientificName))]

#---------------------
# Whitmee & Orme 2013
#---------------------
## mammal natal dispersal distances
wo <- read_csv("data-raw/dispersal-data/Whitmee&Orme/Whitmee_and_Orme_2013.csv")
colnames(wo) <- str_replace_all(colnames(wo), "\\ ", "_")

length(unique(wo$Species)) #104 spp

wo_harm <- harmonize(wo$Species)

## fix some taxonomy
notfound <- filter(wo_harm, is.na(db_code))
# none

## bind:
wo_harm <-filter(wo_harm, !is.na(db_code)) 

## rename columns 
wo <- wo %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

wo <- left_join(wo, wo_harm, by = c("reported_name" = "species")) %>%
  unique()

## Whitmee & Orme
wo_dd <- wo %>%
  select(all_of(cols_to_keep), Value, Units, Sex, Measure, Ref_no) %>%
  mutate(Code = ifelse(Measure == "Mean",
                       "MeanDispersalDistance", 
                       ifelse(Measure == "Maximum",
                              "MaxDispersalDistance", 
                              ifelse(Measure == "Median",
                                     "MedianDispersalDistance", 
                                     NA)))) %>%
  mutate(Unit = ifelse(str_detect(.$Units, "km"), "km",
                       ifelse(str_detect(.$Units, "Metres"), "m", 
                              ifelse(str_detect(.$Units, "Miles"), "miles", NA)))) %>%
  mutate(Database = "Whitmee & Orme 2013", ObservationTypeSpecific = "individual movement distance", 
         Field = "Value") %>%
  rename("DispersalDistance" = Value, "Source" = Ref_no) %>%
  filter(!is.na(DispersalDistance)) %>%
  select(-Measure, -Units) 

## get rid of species in Sutherland to avoid duplicating observations
length(which(unique(wo_dd$scientificName) %in% unique(suth_mamm_dd$scientificName))) ## 48
length(which(unique(wo_dd$scientificName) %in% unique(suth_bird_dd$scientificName))) ## 0

wo_dd <- filter(wo_dd, !scientificName %in% c(suth_mamm_dd$scientificName, suth_bird_dd$scientificName))

## check how many species in bioshifts 
length(which(unique(wo_dd$scientificName) %in% unique(sp$scientificName))) ## 2
wo_sp <- unique(wo_dd$scientificName)[which(unique(wo_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# Jenkins et al. 2007
#---------------------
## for birds, either maximum of either natal or breeding dispersal distance
## for plants, maximum seed dispersal distance
jenkins <- read_csv("data-raw/dispersal-data/Jenkins/Jenkins_et_al_2007.csv")
colnames(jenkins) <- str_replace_all(colnames(jenkins), "\\ ", "_")
length(unique(jenkins$Scientific_Name)) #792 spp

## fix some taxonomy
## get rid of subspecies 
jenkins$genus = str_split_fixed(jenkins$Scientific_Name, " ", 2)[,1]
jenkins$species = str_split_fixed(jenkins$Scientific_Name, " ", 3)[,2]
jenkins$subspecies = str_split_fixed(jenkins$Scientific_Name, " ", 3)[,3]

jenkins = jenkins %>%
  mutate(genus = ifelse(genus == "", NA, genus)) %>%
  mutate(species = ifelse(species == "", NA, species)) %>%
  mutate(subspecies = ifelse(subspecies == "", NA, subspecies)) %>%
  mutate(genus_species = ifelse(is.na(species), genus,
                                paste(genus, species, sep = " ")))

## fix taxonomy manually 
jenkins$reported_name_fixed = jenkins$Scientific_Name
jenkins$reported_name_fixed[which(jenkins$reported_name_fixed == "Aedes taenerhencys")] <- "Aedes taeniorhynchus"
jenkins$reported_name_fixed[which(jenkins$reported_name_fixed == "Artimesia frigida")] <- "Artemisia frigida"
jenkins$reported_name_fixed[which(jenkins$reported_name_fixed == "Artimesia herba-alba")] <- "Artemisia herba-alba"
jenkins$reported_name_fixed[which(jenkins$reported_name_fixed == "Paruline hochequeue")] <- "Parkesia motacilla"
jenkins$reported_name_fixed[which(jenkins$reported_name_fixed == "Solidago alpestris")] <- "Solidago virgaurea"

## harmonize taxonomy in each 
jenkins_harm <- harmonize(jenkins$reported_name_fixed)
jenkins_harm$scientificName[which(jenkins_harm$scientificName == "NA ")] <- jenkins_harm$species[which(jenkins_harm$scientificName == "NA ")]
notfound <- filter(jenkins_harm, is.na(db_code)) ## quite a few not found - many not classified to sp level

## rename columns 
jenkins <- jenkins %>%
  rename("reported_name" = Scientific_Name)

jenkins <- left_join(jenkins, jenkins_harm, by = c("reported_name_fixed" = "species"),
                     relationship = "many-to-many") %>%
  unique()

## Jenkins
jenkins_dd <- jenkins %>%
  select(all_of(cols_to_keep), 
         `Max._Indiv._Dispersal_Distance_(m)`) %>%
  rename("DispersalDistance" = `Max._Indiv._Dispersal_Distance_(m)`) %>%
  mutate(ObservationTypeSpecific = ifelse(kingdom == "Animalia", "individual movement distance",
                                          "seed/plant dispersal (unknown)")) %>%
  mutate(Sex = NA, Unit = "m", 
         Field = "Max._Indiv._Dispersal_Distance_(m)",
         Code = "MaxDispersalDistance", Database = "Jenkins et al. 2007",
         Source = NA) %>%
  filter(!is.na(DispersalDistance))%>%
  filter(!is.na(scientificName)) %>%
  unique()

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 299
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# Flores et al. 2013
#---------------------
## seed traps, tracked individual seeds, marked and recaptured seeds and estimated dispersal distances based on tracking vectors and calculating gut or fur retention times
flores <- read.delim("data-raw/dispersal-data/Flores/Flores_et_al_2013.txt") %>%
  select(-family, -order)
colnames(flores) <- str_replace_all(colnames(flores), "\\ ", "_")
length(unique(flores$Species)) #56 spp

## fix taxonomy manually 
flores$reported_name_fixed = flores$Species
flores$reported_name_fixed[which(flores$reported_name_fixed == "Erodyum cicutarium")] <- "Erodium cicutarium"

## harmonize taxonomy
flores_harm <- harmonize(flores$reported_name_fixed)
notfound <- filter(flores_harm, is.na(db_code)) #one not found

## rename columns 
flores <- flores %>%
  rename("reported_name" = Species) 

flores <- left_join(flores, flores_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique()

## Flores
flores_dd <- flores %>%
  select(all_of(cols_to_keep), Mean.dispersal.distance.m., Maximum.dispersal.distance.m.,
         Reference) %>%
  gather(key = "Field", value = "DispersalDistance", 
         c(Mean.dispersal.distance.m., Maximum.dispersal.distance.m.)) %>%
  mutate(Code = ifelse(Field == "Mean.dispersal.distance.m.", 
                       "MeanDispersalDistance", 
                       ifelse(Field == "Maximum.dispersal.distance.m.",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Sex = NA, ObservationTypeSpecific = "seed/plant dispersal (field)", Unit = "m", 
         Database = "Flores et al. 2013") %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("Source" = Reference) %>%
  mutate(Source = ifelse(str_detect(Source, "Bossard CC. The Role of Habitat Disturbance"),
                         "Bossard CC. The Role of Habitat Disturbance, Seed Predation and Ant Dispersal on Establishment of the Exotic Shrub Cytisus scoparius in California. American Midland Naturalist. 1991;126(1):1",
                         Source)) %>%
  distinct() %>% 
  filter(!is.na(scientificName))

## some duplicates from Flores in Jenkins
## exclude observations from Jenkins with same Code and DispersalDistance as those from Flores since has more info
key = flores_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(flores_dd$scientificName) %in% unique(sp$scientificName))) ## 34
flores_sp <- unique(flores_dd$scientificName)[which(unique(flores_dd$scientificName) %in% unique(sp$scientificName))]

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 297
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# Tamme et al 2014
#---------------------
tamme <- read_csv("data-raw/dispersal-data/Tamme/Tamme_DispersalDistanceData.csv")
colnames(tamme) <- str_replace_all(colnames(tamme), "\\ ", "_")
length(unique(tamme$Species)) #576 spp

## get rid of subspecies 
tamme$genus = str_split_fixed(tamme$Species, " ", 3)[,1]
tamme$species = str_split_fixed(tamme$Species, " ", 3)[,2]
tamme$subspecies = str_split_fixed(tamme$Species, " ", 3)[,3]

tamme = tamme %>%
  mutate(genus = ifelse(genus == "", NA, genus)) %>%
  mutate(species = ifelse(species == "", NA, species)) %>%
  mutate(subspecies = ifelse(subspecies == "", NA, subspecies)) %>%
  mutate(genus_species = ifelse(is.na(species), genus,
                                paste(genus, species, sep = " ")))

## fix taxonomy 
tamme$reported_name_fixed = tamme$genus_species
tamme$reported_name_fixed[which(tamme$reported_name_fixed == "Viola stricta")] <- "Viola striata"
tamme$reported_name_fixed[which(tamme$reported_name_fixed == "Pachira quinata")] <- "Bombax ceiba"
tamme$reported_name_fixed[which(tamme$reported_name_fixed == "Inula conyzae")] <- "Pentanema squarrosum"

## harmonize taxonomy in each 
tamme_harm <- harmonize(tamme$reported_name_fixed)
tamme_harm$scientificName = ifelse(tamme_harm$scientificName == "NA ", 
                                   tamme_harm$species, 
                                   tamme_harm$scientificName)
tamme_harm$scientificName = ifelse(is.na(tamme_harm$scientificName), 
                                   tamme_harm$species, 
                                   tamme_harm$scientificName)
notfound <- filter(tamme_harm, is.na(db_code))

## rename columns 
tamme <- tamme %>%
  rename("reported_name" = Species)

tamme <- left_join(tamme, tamme_harm, by = c("reported_name_fixed" = "species"),
                     relationship = "many-to-many") %>%
  unique()

## write
write.csv(tamme, "data-processed/Tamme_harmonized.csv", row.names = FALSE)


## Tamme
tamme_dd <- tamme %>%
  select(all_of(cols_to_keep), 
         `Maximum_recorded_dispersal_distance_(m)`,
         `Mean_dispersal_distance_(m)`,
         `Median_dispersal_distance_(m)`, 
         `Mode_dispersal_distance_(m)`,
         `90th_percentile_dispersal_distance_(m)`,
         `99th_percentile_dispersal_distance_(m)`,
         Data_type, Reference) %>%
  gather(key = "Field", value = "DispersalDistance", c(`Maximum_recorded_dispersal_distance_(m)`,
                                                       `Mean_dispersal_distance_(m)`,
                                                       `Median_dispersal_distance_(m)`, 
                                                       `Mode_dispersal_distance_(m)`,
                                                       `90th_percentile_dispersal_distance_(m)`,
                                                       `99th_percentile_dispersal_distance_(m)`)) %>%
  mutate(Code = ifelse(str_detect(.$Field, "Mean"),
                       "MeanDispersalDistance", 
                       ifelse(str_detect(.$Field, "Maximum"),
                              "MaxDispersalDistance", 
                              ifelse(str_detect(.$Field, "Median"),
                                     "MedianDispersalDistance", 
                                     ifelse(str_detect(.$Field, "90"),
                                            "90thPercentileDispersalDistance", 
                                            ifelse(str_detect(.$Field, "99"),
                                                   "99thPercentileDispersalDistance",
                                                   ifelse(str_detect(.$Field, "Mode"),
                                                          "ModeDispersalDistance",
                                                          NA))))))) %>%
  mutate(Unit = "m", Database = "Tamme et al. 2014", Sex = NA) %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("ObservationTypeSpecific" = Data_type, "Source" = Reference) %>%
  mutate(ObservationTypeSpecific = ifelse(ObservationTypeSpecific == "field", 
                                          "seed dispersal (field)", 
                                          ifelse(ObservationTypeSpecific == "model", 
                                                 "seed dispersal (model)",
                                                 NA)))

## Tamme took from Jenkins
## exclude observations from Jenkins with same Code and DispersalDistance as Tamme since has more info
key = tamme_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## also exclude ones with different Code from source Augspurger 1986
key = tamme_dd %>%
  filter(Source == "Augspurger 1986") %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 217
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


## looks like Tamme and Flores have some overlapping data but rounded differently 
## if duplicated after rounding and source is from Stamp 1989 or Yumoto 1999, get rid of it in Tamme
key = flores_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, round(DispersalDistance, 2), sep = "_")) %>%
  select(data_id)

tamme_dd <- tamme_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!(paste(scientificName, Code, round(DispersalDistance, 2), sep = "_") %in% 
             key$data_id & Source %in% c("Stamp 1989", "Yumoto 1999")), 
         Source != "Smith & Kok 1984")

## check how many species in bioshifts 
length(which(unique(tamme_dd$scientificName) %in% unique(sp$scientificName))) ## 334
tamme_sp <- unique(tamme_dd$scientificName)[which(unique(tamme_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# TRY
#---------------------
## read in TRY query
try = read_delim("data-raw/dispersal-data/TRY-query/23454.txt")
unique(try$TraitName)

## get only dispersal data 
try_dd <- try %>%
  filter(TraitName == "Dispersal distance")
length(unique(try_dd$SpeciesName)) #152 spp

## harmonize
try_dd_harm <- harmonize(try_dd$SpeciesName)
notfound <- filter(try_dd_harm, is.na(db_code)) ## 0 species not found

## rename columns 
try_dd <- try_dd %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

try_dd <- left_join(try_dd, try_dd_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique()

## TRY
try_dd_dd <- try_dd %>%
  select(all_of(cols_to_keep), 
         OrigValueStr,
         OriglName, OrigUnitStr, Reference) %>%
  rename("Source" = Reference, "Field" = OriglName, "DispersalDistance" = OrigValueStr, "Unit" = OrigUnitStr) %>%
  mutate(Code = ifelse(Field %in% c("MeanDispersalDistanceMean", "MeanDispersalDistanceMax", 
                                    "MeanDispersalDistanceMin"), 
                       "MeanDispersalDistance", 
                       ifelse(Field %in%c("MaxDispersalDistanceMax", "Max Seed Dispersal Distance",
                                          "MaxDispersalDistanceMin", "MaxDispersalDistanceMean",
                                          "Seed_dispersal_distance_95",
                                          "Seed_dispersal_distance_5", "Max Seed Dispersal Distance"),
                              "MaxDispersalDistance", 
                              ifelse(Field %in% c("Effective Seed Dispersal Distance", "Seed Dispersal distance"),
                                     "DispersalDistance",
                                     NA)))) %>%
  mutate(Sex = NA, ObservationTypeSpecific = "seed/plant dispersal (unknown)",  
         Database = "TRY database") %>% 
  filter(!is.na(Unit)) %>%
  filter(!is.na(scientificName)) %>%
  filter(Source != "unpub.") ## get rid of unpub. observations

## check how many species in bioshifts 
length(which(unique(try_dd_dd$scientificName) %in% unique(sp$scientificName))) ## 38
try_dd_sp <- unique(try_dd_dd$scientificName)[which(unique(try_dd_dd$scientificName) %in% unique(sp$scientificName))]


#---------------------
# Chu 2021
#---------------------
## read in data
chu = read_csv("data-raw/dispersal-data/Chu/Chu_Natal_Dispersal_Estimates_Nikki.csv")
colnames(chu)[1] <- "common name"
length(unique(chu$sciname)) #104 spp

## harmonize taxonomy
chu_harm <- harmonize(chu$sciname)
notfound <- filter(chu_harm, is.na(db_code)) ## no species not found 

## rename columns 
chu <- chu %>%
  rename("reported_name" = sciname) %>%
  mutate(reported_name_fixed = reported_name)

chu <- left_join(chu, chu_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## if n=1, label as single observation of dispersal
chu$Code = ifelse(chu$n == 1, "SingleObservation", "MeanDispersalDistance")

## reorganize
chu_dd <- chu %>%
  select(all_of(cols_to_keep), geom, arimean, Code) %>%
  gather(key = "Field", value = "DispersalDistance", c(geom, arimean)) %>%
  mutate(Source = "Chu 2021",
         Sex = NA, Unit = "km",
         ObservationTypeSpecific = "natal dispersal distance",
         Database = "Chu 2021") 

## check how many species in bioshifts 
length(which(unique(chu_dd$scientificName) %in% unique(sp$scientificName))) ## 81
chu_sp <- unique(chu_dd$scientificName)[which(unique(chu_dd$scientificName) %in% unique(sp$scientificName))]

#---------------------
# Vittoz and Engler 2007
#---------------------
## read in data
vitt = read_csv("data-raw/dispersal-data/Vittoz&Engler/Vittoz_and_Engler_2007.csv")[,1:6]
length(unique(vitt$Species)) #175 spp

## fix species name  
vitt$reported_name_fixed <- vitt$Species
vitt$reported_name_fixed[which(vitt$Species == "Viola stricta")] <- "Viola striata" ## Viola stricta <- viola striata 

## harmonize taxonomy
vitt_harm <- harmonize(vitt$reported_name_fixed)
notfound <- filter(vitt_harm, is.na(db_code)) ## many species not found, but a lot are only identified to genus 
notfound <- filter(notfound, !str_detect(notfound$species, "sp."))
notfound

## rename columns 
vitt <- vitt %>%
  rename("reported_name" = Species) 

vitt <- left_join(vitt, vitt_harm, by = c("reported_name_fixed" = "species")) %>%
  unique() %>%
  filter(!is.na(scientificName), scientificName != "NA ") 

## clean the data 
vitt <- vitt %>%
  mutate(Type_of_measurement_general = ifelse(str_detect(Type_of_measurement, "Maxim"), "MaxDispersalDistance",
                                              ifelse(str_detect(Type_of_measurement, "Mean"), "MeanDispersalDistance",
                                                     ifelse(str_detect(Type_of_measurement, "Median"), "MedianDispersalDistance",
                                                            ifelse(str_detect(Type_of_measurement, "median"), "MedianDispersalDistance",
                                                            ifelse(str_detect(Type_of_measurement, "Minimum"), "MinDispersalDistance",
                                                     ifelse(str_detect(Type_of_measurement, "99"), "99thPercentileDispersalDistance",
                                                            ifelse(str_detect(Type_of_measurement, "90"), "90thPercentileDispersalDistance",
                                                                   ifelse(str_detect(Type_of_measurement, "80"), "80thPercentileDispersalDistance",
                                                                          ifelse(str_detect(Type_of_measurement, "70"), "70thPercentileDispersalDistance",
                                                            NA)))))))))) %>%
  filter(!is.na(Type_of_measurement_general))

## values with asterisks represent long distance dispersal events 
## get rid of < and > 
vitt$Distance_m <- str_replace_all(vitt$Distance_m, "\\<", "")
vitt$Distance_m <- str_replace_all(vitt$Distance_m, "\\>", "")
vitt$Distance_m <- str_replace_all(vitt$Distance_m, "\\ ", "")
vitt <- vitt %>% mutate(is_LDD = ifelse(str_detect(Distance_m, "\\*"), 1, 0))
vitt$Distance_m <- str_replace_all(vitt$Distance_m, "\\*", "")
vitt$Distance_m <- str_replace_all(vitt$Distance_m, "\\'", "")
vitt$Distance_m[which(str_detect(vitt$Distance_m, "\\-"))] <- str_split_fixed(vitt$Distance_m[which(str_detect(vitt$Distance_m, "\\-"))], "\\-", 2)[,2]
vitt$Distance_m = as.numeric(vitt$Distance_m)

vitt %>%
  ggplot(aes(x = Distance_m)) + geom_histogram() + facet_wrap(~is_LDD)

## reorganize
vitt_dd <- vitt %>%
  select(all_of(cols_to_keep), Type_of_measurement, Type_of_measurement_general, Distance_m, Reference) %>%
  rename("Field" = Type_of_measurement, "DispersalDistance" = Distance_m, "Source" = Reference,
         "Code" = Type_of_measurement_general) %>%
  mutate(Sex = NA, Unit = "m",
         ObservationTypeSpecific = "seed dispersal (unknown)",
         Database = "Vittoz & Engler 2007") 

## check how many species in bioshifts 
length(which(unique(vitt_dd$scientificName) %in% unique(sp$scientificName))) ## 133
vitt_sp <- unique(vitt_dd$scientificName)[which(unique(vitt_dd$scientificName) %in% unique(sp$scientificName))]

## check how many unique species we have data for:
species_with_dd <- append(tamme_sp, jenkins_sp) %>%
  append(., suth_bird_sp) %>%
  append(., suth_mamm_sp) %>%
  append(., wo_sp) %>%
  append(., flores_sp) %>%
  append(., try_dd_sp) %>%
  append(., par_sp) %>%
  append(., chu_sp) %>%
  append(., vitt_sp)

length(unique(species_with_dd)) # 615 species 

## now: make subsets of each database with only bioshifts species 
tamme_sub = filter(tamme_dd, scientificName %in% tamme_sp)
suth_bird_sub = filter(suth_bird_dd, scientificName %in% suth_bird_sp)
suth_mamm_sub = filter(suth_mamm_dd, scientificName %in% suth_mamm_sp)
wo_sub = filter(wo_dd, scientificName %in% wo_sp)
jenkins_sub = filter(jenkins_dd, scientificName %in% jenkins_sp)
flores_sub = filter(flores_dd, scientificName %in% flores_sp)
try_dd_sub = filter(try_dd_dd, scientificName %in% try_dd_sp)
par_sub = filter(par_dd, scientificName %in% par_sp)
chu_sub = filter(chu_dd, scientificName %in% chu_sp)
vitt_sub = filter(vitt_dd, scientificName %in% vitt_sp)

## collate all bioshifts data
dd_collated <- rbind(tamme_sub, suth_bird_sub) %>%
  rbind(., suth_mamm_sub) %>%
  rbind(., wo_sub) %>%
  rbind(., jenkins_sub) %>%
  rbind(., flores_sub) %>%
  rbind(., try_dd_sub) %>%
  rbind(., par_sub) %>%
  rbind(., chu_sub) %>%
  rbind(., vitt_sub) %>%
  unique()

length(unique(dd_collated$scientificName)) # 615

length(which(!dd_collated$scientificName %in% sp$scientificName))

## write intermediate file 
write.csv(dd_collated, "data-processed/dispersal-distance-collated_intermediate.csv", row.names = FALSE)
dd_collated <- read.csv("data-processed/dispersal-distance-collated_intermediate.csv")

## fix classes that are missing
v3 <- v3 %>%
  select(scientificName, class, kingdom, phylum, family) %>%
  distinct()

missing_class <- dd_collated %>%
  filter(is.na(class)) 

missing_class_sp <- filter(v3, scientificName %in% missing_class$scientificName) 

missing_class <- select(missing_class, -class) %>%
  left_join(., missing_class_sp)

dd_collated <- filter(dd_collated, !scientificName %in% missing_class$scientificName) %>%
  rbind(., missing_class)

## fix Chloris chloris 
dd_collated[which(dd_collated$scientificName == "Chloris chloris"), c(4,5,6,7,8)] <-
  rep(c("Animalia", "Chordata", "Aves", "Passeriformes", "Fringillidae"), each = 7)

## get rid of sea anemone, shark, fish
dd_collated <- filter(dd_collated, class != "Ascidiacea", class != "Elasmobranchii", class != "Actinopteri")

## make column for group and for type of dispersal estimate
dd_collated$group <- ifelse(dd_collated$class %in% c("Magnoliopsida", "Pinopsida",
                                                     "Liliopsida"), "Plants", dd_collated$class)

dd_collated <- filter(dd_collated, !(group == "Aves" & str_detect(ObservationTypeSpecific, "seed")))

dd_collated <- dd_collated %>%
  mutate(type = ifelse(Code %in% c("MaxDispersalDistance", "90thPercentileDispersalDistance",
                                                   "99thPercentileDispersalDistance", 
                                   "70thPercentileDispersalDistance"), "Max",
                           ifelse(Code == "MeanDispersalDistance", "Mean", 
                                  ifelse(Code == "MedianDispersalDistance", "Median",
                                         ifelse(Code == "ModeDispersalDistance", "Mode",
                                                ifelse(Code == "MinDispersalDistance", "Min",
                                                "Unknown")))))) %>%
  mutate(ObservationTypeSpecific = ifelse(Field %in% c("arimean", "ArithmeticMeanNatalDispersal"), 
                                          "arithmetic mean natal dispersal distance",
                                          ifelse(Field %in% c("geom", "GeometricMeanNatalDispersal"), 
                                                              "geometric mean natal dispersal distance",
                                                 ObservationTypeSpecific))) %>%
  filter(type != "Min")


length(unique(dd_collated$scientificName)) #604 species

## make general observation type column 
## general types of studies:
## 1. natal dispersal
## 2. breeding dispersal
## 3. seed dispersal

unique(dd_collated$ObservationTypeSpecific)
dd_collated <- dd_collated %>%
  mutate(ObservationTypeGeneral = ifelse(str_detect(ObservationTypeSpecific, "seed"),
                                         "seed dispersal", ObservationTypeSpecific)) %>%
  mutate(ObservationTypeGeneral = ifelse(str_detect(ObservationTypeSpecific, "natal"),
                                         "natal dispersal", ObservationTypeGeneral)) %>%
  mutate(ObservationTypeGeneral = ifelse(str_detect(ObservationTypeSpecific, "breeding"),
                                         "breeding dispersal", ObservationTypeGeneral)) %>%
  mutate(ObservationTypeGeneral = ifelse(str_detect(ObservationTypeSpecific, "individual"),
                                         "natal/breeding dispersal", ObservationTypeGeneral))

unique(dd_collated$ObservationTypeGeneral)
  
## convert to km 
unique(dd_collated$Unit)

dd_collated <- dd_collated %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(DispersalDistanceKm = ifelse(Unit == "m",
                                      DispersalDistance/1000, 
                                      ifelse(Unit == "miles",
                                             DispersalDistance*1.60934,
                                             DispersalDistance))) %>%
  filter(!is.na(DispersalDistanceKm)) 

## add back checked scientific name 
dd_collated <- left_join(dd_collated, sp) %>%
  select(reported_name, reported_name_fixed, scientificName, scientificName_checked, everything())

write.csv(dd_collated, "data-processed/dispersal-distance-collated.csv", row.names = FALSE)
dd_collated <- read.csv("data-processed/dispersal-distance-collated.csv")


#----------------------------
# Make some plots 
#----------------------------
dd_collated %>%
  select(scientificName, group, class, DispersalDistanceKm) %>%
  unique() %>% 
  ggplot(aes(x = log(DispersalDistanceKm), fill = class)) + geom_histogram() +
  theme_classic() + coord_flip() + facet_wrap(~group)

dd_collated %>%
  select(scientificName, group, Code, DispersalDistanceKm) %>%
  unique() %>% 
  ggplot(aes(x = log(DispersalDistanceKm), fill = Code)) + geom_histogram() +
  theme_classic() + coord_flip() + facet_wrap(~group)

dd_collated %>%
  select(scientificName, group, Code, DispersalDistanceKm) %>%
  unique() %>% 
  ggplot(aes(x = log(DispersalDistanceKm), fill = group)) + geom_histogram() +
  theme_classic() + coord_flip() + facet_wrap(~Code)

## compare mean versus max versus median versus mode
dd_collated %>%
  select(scientificName, group, Code, type, DispersalDistanceKm) %>%
  unique() %>% 
  ggplot(aes(x = log(DispersalDistanceKm), fill = group)) + geom_histogram() +
  theme_classic() + coord_flip() + facet_wrap(~type)

## how different are geometric versus arithmetic mean distances? 
dd_collated %>%
  filter(ObservationTypeGeneral == "natal dispersal",
         type == "Mean") %>%
  select(scientificName, ObservationTypeSpecific, Field, DispersalDistanceKm) %>%
  unique() %>% 
  ggplot(aes(x = log(DispersalDistanceKm), fill = ObservationTypeSpecific)) + geom_histogram() +
  theme_classic() + coord_flip() + facet_wrap(~ObservationTypeSpecific)
## not that different 


