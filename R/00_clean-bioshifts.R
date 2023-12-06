## hamronize names in bioshifts v1 
library(readr)
library(tidyverse)
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

## read function to harmonize taxonomy
source("R/taxonomic-harmonization/harmonize.R")

#----------------------
## read in bioshifts v1 
## version downloaded from original paper 
v1 = read.table("data-raw/BIOSHIFTSv1/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",h=T,dec=".",
                stringsAsFactors = FALSE) 
v1$reported_name = v1$Species

## clean names in the orginal database 
v1$Species = Clean_Names(v1$Species, return_gen_sps = F)

## manual cleaning 
v1$Species[which(v1$Species == "Quercus x")] = "Quercus" 
v1$Species[which(v1$Species == "Mentha x")] = "Mentha" 
v1$Species[which(v1$Species == "Circaea x intermedia")] = "Circaea intermedia" 

v1$Species = Clean_Names(v1$Species, return_gen_sps = F)

## now, harmonize names
bs_harm <- harmonize(v1$Species)
notfound <- filter(bs_harm, is.na(db_code))

## add back species that weren't found but are resolved to species level
notfound <- filter(notfound, !str_detect(notfound$species, "sp."))
notfound$scientificName = notfound$species
bs_harm <- rbind(bs_harm, notfound) %>%
  distinct()

## rename columns 
v1 <- v1 %>%
  rename("reported_name_fixed" = Species)

v1_corrected <- left_join(v1, bs_harm, by = c("reported_name_fixed" = "species"),
                    relationship = "many-to-many") %>%
  unique() %>%
  filter(!is.na(scientificName))

length(unique(v1_corrected$scientificName)) # 12189 species 

## save version of bioshifts with fixed taxonomy 
#################################################
# 'reported_name' is the name originally in the database 
# 'reported_name_fixed' is the name after cleaning
# 'scientificName' is the name after cleaning and harmonizing the reported name
write.csv(v1_corrected, "data-processed/BIOSHIFTSv1_harmonized.csv", row.names = FALSE)
