## brings together published data on the age at maturity / generation time of birds and plants 
## this will serve as a measure of the frequency of dispersal events
library(readr)
library(tidyverse)
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

## read function to harmonize taxonomy
source("R/taxonomic-harmonization/harmonize.R")

## read in dispersal data
dd_collated <- read.csv("data-processed/dispersal-distance-collated.csv")

cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")

#---------------------
# PLANTS 
#---------------------
## collate age at maturity data
###############################
## for plants 
#---------------------
# TRY
#---------------------
try_am = read_delim("data-raw/dispersal-data/TRY-query/23659.txt")
unique(try_am$TraitName)

## age at maturity 
try_am <- try_am %>%
  filter(TraitName %in% c("Plant ontogeny: age of maturity (first flowering)", "Dispersal frequency")) 

## harmonize names 
try_am_harm <- harmonize(try_am$SpeciesName)
notfound <- filter(try_am_harm, is.na(db_code))

## fuzzy match to make sure we don't accidentally miss a species 
match_dd <- select(dd_collated, scientificName) %>%
  unique() %>%
  mutate(key1 = 1:nrow(.))

match_notfound <- select(notfound, species) %>%
  unique() %>%
  mutate(key2 = 1:nrow(.))

fuzz = fedmatch::merge_plus(data1 = match_dd, 
           data2 = match_notfound,
           by.x = "scientificName",
           by.y = "species", 
           match_type = "fuzzy", 
           unique_key_1 = "key1",
           unique_key_2 = "key2",
           fuzzy_settings = fedmatch::build_fuzzy_settings(maxDist = .1))

fuzz = fuzz[[1]]
## none are the same species

## rename columns 
try_am <- try_am %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

try_am <- left_join(try_am, try_am_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique()

## write 
write.csv(try_am, "data-processed/TRY_age-at-maturity_harmonized.csv", row.names = FALSE)

## subset to bioshifts species with dispersal distance
try_am_bs <- filter(try_am, scientificName %in% dd_collated$scientificName)

length(unique(try_am_bs$scientificName)) # 268 species 

length(unique(dd_collated$scientificName[which(dd_collated$group == "Plants")])) # out of the 399 plants - nice!!! 

## clean the data 
unique(try_am_bs$OriglName)
# keep:
## WholePlant SMINP = age to first seed production
## primary juvenile period = Plant age at first flowering
## Annual seed dispersal, maturity age, floqwering time, first flowering time

# get rid of:
## WholePlant SMINP = frequency of good crop
# flowering = indicates whether flowered during first year in mesocosm

try_am_bs = try_am_bs %>%
  ## split units from value
  filter(!OrigValueStr %in% c("mast", "regular", "irregular", "alternate")) %>% 
  mutate(Unit = str_split_fixed(OrigValueStr, "\\ ", 2)[,2],
         AgeAtMaturity = str_split_fixed(OrigValueStr, "\\ ", 2)[,1]) %>% 
  mutate(Unit = ifelse(Unit == "", OrigUnitStr, Unit)) %>%
  filter(!OriglName %in% c("flowering", "Secondary juvenile period", "WholePlant SMINP")) %>%
  mutate(Unit = ifelse(str_detect(Unit, "yrs") | Unit %in% c("years?", "yr", "years"), 
                       "yrs",
                       ifelse(str_detect(Unit, "day"), 
                              "days", 
                              Unit))) %>%
  mutate(Unit = ifelse(is.na(Unit) & 
                         DataName %in% c("Frequency of good seeding years", "Age of maturity of saplings"), 
                              "yrs", 
                       ifelse(is.na(Unit) & Dataset %in% c("Fazlioglu et al. 2017", "Fazlioglu et al. 2016"), 
                       "days", 
                       ifelse(is.na(Unit), "yrs", Unit))))

## make sure no NA units 
length(which(is.na(try_am_bs$Unit)))

unique(try_am_bs$Unit)
## fix ones that are ranges 
## want minimum age at maturity to give a maximum dispersal estimate 

try_am_bs = try_am_bs %>%
  mutate(AgeAtMaturity = ifelse(Unit == "5 years", 5, 
                                ifelse(Unit %in% c("1 and 5 years", "1 year"), 1,
                                       AgeAtMaturity))) %>%
  mutate(Unit = ifelse(Unit %in% c("5 years", "1 and 5 years", "1 year"), "yrs", Unit)) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity,"\\<", ""),
         AgeAtMaturity = str_replace_all(AgeAtMaturity,"\\>", ""),
         AgeAtMaturity = str_split_fixed(AgeAtMaturity, "\\-", 2)[,1]) 

unique(try_am_bs$AgeAtMaturity)
unique(try_am_bs$Unit)

try_am_bs <- try_am_bs %>%
  select(all_of(cols_to_keep), 
         Unit, AgeAtMaturity, OriglName, 
         Reference) %>%
  mutate(Database = "TRY", Code = "AgeAtMaturity") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  rename("Source" = Reference, "Field" = OriglName)

## save 
write.csv(try_am_bs, "data-processed/age-at-maturity-TRY.csv", row.names = FALSE)

## fill in the gaps:
missing_am = dd_collated %>%
  filter(kingdom == "Plantae" & !scientificName %in% try_am_bs$scientificName) %>%
  select(scientificName, family, order, class, phylum, kingdom) %>%
  unique()


## use longevity/growth form information to infer age at maturity for annual plants
####################################################################################
## for species that live 1 year, age at maturity = 1 year 
## longevity: 
lifespan <- read.csv("data-raw/CompilationLifeSpan_02242022.csv")
unique(lifespan$Database)

## filter to species with dispersal distance
lifespan <- filter(lifespan, SpeciesChecked %in% dd_collated$scientificName)
length(unique(lifespan$SpeciesChecked)) ## 444 of our species 

## how many are plants?
lifespan = filter(lifespan, phylum %in% c("Streptophyta", "Tracheophyta"))
length(unique(lifespan$SpeciesChecked)) ## 286 of our species 

## harmonize names 
lifespan_harm <- harmonize(lifespan$SpeciesChecked)
notfound <- filter(lifespan_harm, is.na(db_code)) # none

## rename columns 
lifespan <- lifespan %>%
  rename("reported_name" = SpeciesChecked) %>%
  mutate(reported_name_fixed = reported_name) %>%
  select(-family, -genus, -phylum, -class, -order)

lifespan <- left_join(lifespan, lifespan_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## clean lifespan
unique(lifespan$LifeSpan)
unique(lifespan$LifeSpan)[str_detect(unique(lifespan$LifeSpan), "1\\-")]
unique(lifespan$LifeSpan)[str_detect(unique(lifespan$LifeSpan), "1\\,")]

lifespan_am <- lifespan %>%
  mutate(Unit = ifelse(LifeSpan %in% c("annuals", "summer annuals",
                                       "winter annuals"),
                       "years",
                       Unit),
         Code = "MaturityFromLifespan") %>%
  mutate(AgeAtMaturity = ifelse(LifeSpan %in% c("annuals", "summer annuals",
                                                "winter annuals"),
                                1,
                                ifelse(str_detect(LifeSpan, "1\\-"), 
                                       1,
                                       ifelse(str_detect(LifeSpan, "1\\,"), 
                                              1,
                                              ifelse(str_detect(LifeSpan, "\\<1"), 
                                                     1,
                                                     ifelse(LifeSpan == 1,
                                                            1, 
                                                            NA)))))) %>%
  filter(AgeAtMaturity == 1) %>% 
  select(all_of(cols_to_keep), -X, -Taxonomic.Group, -Ecosystem.Type, -Group, -ecotype,
         AgeAtMaturity, Unit, Database, Field, Code) %>%
  unique() %>%
  filter(scientificName %in% missing_am$scientificName) %>%
  mutate(Source = "Lifespan compilation")

length(unique(lifespan_am$scientificName)) #8 more species 

#---------------------
# BROT
#---------------------
brot = read.csv("data-raw/dispersal-data/BROT/BROT2_dat.csv")
unique(brot$Trait)

## GrowthForm, DispMode
brot = brot %>%
  filter(Trait %in% c("GrowthForm"))

## harmonize names 
brot_harm <- harmonize(brot$Taxon)
notfound <- filter(brot_harm, is.na(db_code)) # lots

## fuzzy match to make sure we don't accidentally miss a species 
match_dd <- select(dd_collated, scientificName) %>%
  unique() %>%
  mutate(key1 = 1:nrow(.))

match_notfound <- select(notfound, species) %>%
  unique() %>%
  mutate(key2 = 1:nrow(.))

fuzz = fedmatch::merge_plus(data1 = match_dd, 
                            data2 = match_notfound,
                            by.x = "scientificName",
                            by.y = "species", 
                            match_type = "fuzzy", 
                            unique_key_1 = "key1",
                            unique_key_2 = "key2",
                            fuzzy_settings = fedmatch::build_fuzzy_settings(maxDist = .1))

fuzz = fuzz[[1]]
## none are the same species

## rename columns 
brot <- brot %>%
  rename("reported_name" = Taxon) %>%
  mutate(reported_name_fixed = reported_name)

brot <- left_join(brot, brot_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts plant species missing data on age at maturity 
brot_sub <- filter(brot, scientificName %in% missing_am$scientificName) %>%
  filter(!scientificName %in% lifespan_am$scientificName)

## use information on whether each species is annual or perennial 
unique(brot_sub$Data)

brot_am <- brot_sub %>%
  select(all_of(cols_to_keep), Data, SourceID) %>%
  mutate(AgeAtMaturity = ifelse(str_detect(Data, "annual"),
                                "1", 
                                NA)) %>%
  rename("Source" = SourceID) %>%
  mutate(Field = "Growth form",
         Unit = "yrs", Database = "TRY", Code = "MaturityFromGrowthForm") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  select(-Data)

length(unique(brot_am$scientificName))# 1 more

#-----------------
# BIRDS
#------------------
## age at maturity
#-------
# AnAge
#-------
anage <- read.delim("data-raw/dispersal-data/AnAge/anage_data.txt")
anage$genus_species <- paste(anage$Genus, anage$Species, sep = " ")

## harmonize names 
anage_harm <- harmonize(as.character(anage$genus_species))
notfound <- filter(anage_harm, is.na(db_code)) #none

## rename columns 
anage <- anage %>%
  rename("reported_name" = genus_species) %>%
  mutate(reported_name_fixed = reported_name)

anage <- left_join(anage, anage_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts species with dispersal distance
anage_bs <- filter(anage, scientificName %in% dd_collated$scientificName)
length(unique(anage_bs$scientificName)) # 179 species 

anage_am = anage_bs %>%
  select(all_of(cols_to_keep), Female.maturity..days., Male.maturity..days.,References) %>%
  gather(key = "Field", value = "AgeAtMaturity", 
         c(Male.maturity..days., Female.maturity..days.)) %>%
  rename("Source" = References) %>%
  mutate(Database = "AnAge", 
         Unit = "days", 
         Code = "AgeAtMaturity")%>%
  filter(!is.na(AgeAtMaturity))

length(unique(anage_am$scientificName)) # 169 spp


#-----------------
# Amniota
#------------------
amniota <- read.csv("data-raw/dispersal-data/amniota/ECOL_96_269/Data_Files/Amniote_Database_Aug_2015.csv") %>%
  select(-class, -order, -family)
amniota$genus_species <- paste(amniota$genus, amniota$species, sep = " ")

## harmonize names 
amni_harm <- harmonize(as.character(amniota$genus_species))
notfound <- filter(amni_harm, is.na(db_code)) # lots not found

## fuzzy match to make sure we don't accidentally miss a species 
match_dd <- select(dd_collated, scientificName) %>%
  unique() %>%
  mutate(key1 = 1:nrow(.))

match_notfound <- select(notfound, species) %>%
  unique() %>%
  mutate(key2 = 1:nrow(.))

fuzz = fedmatch::merge_plus(data1 = match_dd, 
                            data2 = match_notfound,
                            by.x = "scientificName",
                            by.y = "species", 
                            match_type = "fuzzy", 
                            unique_key_1 = "key1",
                            unique_key_2 = "key2",
                            fuzzy_settings = fedmatch::build_fuzzy_settings(maxDist = .1))
## none are the same species

## rename columns 
amniota <- amniota %>%
  rename("reported_name" = genus_species) %>%
  mutate(reported_name_fixed = reported_name)

amniota <- left_join(amniota, amni_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts species with dispersal distance
amniota_bs <- filter(amniota, scientificName %in% dd_collated$scientificName)
length(unique(amniota_bs$scientificName)) # 180 spp

amniota_am = amniota_bs %>%
  select(all_of(cols_to_keep), female_maturity_d, male_maturity_d) %>%
  gather(key = "Field", value = "AgeAtMaturity", 
         c(female_maturity_d, male_maturity_d)) %>%
  mutate(Database = "Amniota", 
         Unit = "days", 
         Code = "AgeAtMaturity",
         Source = NA,
         Sex = ifelse(Field == "female_maturity_d", "f", 
                      ifelse(Field == "male_maturity_d", "m", NA)))%>%
  filter(!is.na(AgeAtMaturity)) %>%
  filter(AgeAtMaturity != -999)

length(unique(amniota_am$scientificName)) # 175 spp


## combine
all_am <- rbind(try_am_bs, lifespan_am) %>%
  rbind(., brot_am) %>%
  rbind(., anage_am)
all_am$Sex = NA

all_am <- rbind(all_am, amniota_am)

all_am$GenerationLength = NA

length(unique(all_am$scientificName)) ## 449 spp

## write 
write.csv(all_am, "data-processed/age-at-maturity.csv", row.names = FALSE)

