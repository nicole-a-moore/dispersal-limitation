## extract proxy traits for dispersal
library(tidyverse)
library(readr)
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

## read function to harmonize taxonomy
source("R/taxonomic-harmonization/harmonize.R")

## write vector of columns we want to be in final database
cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")


## read in dispersal data
dd <- read.csv("data-processed/dispersal-distance-collated.csv")

key <- select(dd, scientificName, scientificName_checked) %>% unique()

#---------------------
# PLANTS 
#---------------------
## collate body size data for plants
####################################
#----------------------------
# TRY: Vegetative height 
#----------------------------
try = read.csv("data-raw/dispersal-data/TRY-query/TRY_VegetativeHeight.csv", sep = ";")

## harmonize 
try_harm <- harmonize(try$SpeciesChecked)
notfound <- filter(try_harm, is.na(db_code)) ## 3 species not found

## rename columns 
try <- try %>%
  rename("reported_name" = SpeciesChecked) %>%
  mutate(reported_name_fixed = reported_name)

try_bs <- left_join(try, try_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique()

write.csv(try_bs, "data-processed/TRY_body-size_harmonized.csv", row.names = FALSE)

## subset to our species
try <- try_bs %>%
  select(all_of(cols_to_keep),
         m) %>%
  rename("BodySize"= m) %>%
  mutate(Database = "TRY", 
         Code = "MeanBodySize",
         Unit = "mm",
         Field = "VegetativeHeight") %>%
  filter(scientificName %in% dd$scientificName_checked)

length(unique(try$scientificName)) #344
nrow(try)

## for later
#----------------------------
# Tamme et al.: other traits
#----------------------------
tamme <- read.csv("data-processed/Tamme_harmonized.csv")
colnames(tamme)

## subset to our species
tamme <- tamme %>%
  select(all_of(cols_to_keep),
         Growth_form, Dispersal_syndrome, Seed_weight_.mg., Seed_release_height_.m.) %>%
  rename("GrowthForm" = Growth_form, 
         "DispersalSyndrome" = Dispersal_syndrome, 
         "SeedWeightmg" = Seed_weight_.mg., 
         "SeedReleaseHeightm" = Seed_release_height_.m.) %>%
  mutate(Database = "Tamme et al.") %>%
  filter(scientificName %in% dd$scientificName_checked)

write.csv(tamme, "data-processed/dispersal-proxy-trait-compilation_tamme.csv", row.names = FALSE)


#---------------------
# BIRDS 
#---------------------
## collate body size data for birds
####################################
#------------------------
# AVONET: Tarsus.Length
#------------------------
avonet <- read.csv("data-raw/dispersal-data/AVONET/AVONET_raw.csv")
colnames(avonet)

## harmonize 
avo_harm <- harmonize(avonet$Species1)
notfound <- filter(avo_harm, is.na(db_code)) ## 0 species not found

## rename columns 
avonet <- avonet %>%
  rename("reported_name" = Species1) %>%
  mutate(reported_name_fixed = reported_name)

avonet <- left_join(avonet, avo_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique()

write.csv(avonet, "data-processed/AVONET_body-size_harmonized.csv", row.names = FALSE)

## subset to our species
avo_bs <- avonet %>%
  select(all_of(cols_to_keep),
         Tarsus.Length) %>%
  rename("BodySize"= Tarsus.Length) %>%
  mutate(Database = "AVONET", 
         Code = "BodySize",
         Unit = "mm",
         Field = "TarsusLength") %>%
  filter(scientificName %in% dd$scientificName_checked)

length(unique(avo_bs$scientificName)) #170
nrow(avo_bs)
## some have multiple measures
## take mean for now

avo_bs <- avo_bs %>%
  group_by(scientificName) %>%
  mutate(BodySize = mean(BodySize)) %>%
  ungroup() %>%
  distinct()

## other traits:
avonet <- avonet %>%
  select(all_of(cols_to_keep),
         Hand.Wing.Index, Kipps.Distance) %>%
  rename("HandWingIndex" = Hand.Wing.Index,
         "KippsDistance" = Kipps.Distance) %>%
  mutate(Database = "AVONET",
         Unit = "mm") %>%
  filter(scientificName %in% dd$scientificName_checked)


write.csv(avonet, "data-processed/dispersal-proxy-trait-compilation_avonet.csv", row.names = FALSE)

## good for now 
## combine plants and birds
bodysize <- rbind(avo_bs, try) %>%
  select(-reported_name, -reported_name_fixed, -db_code) %>%
  distinct()

bodysize <- left_join(bodysize, key, by = c("scientificName" = "scientificName_checked")) %>%
  rename("scientificName_checked" = scientificName)

## write
write.csv(bodysize, "data-processed/dispersal-proxy-trait-compilation.csv", row.names = FALSE)

