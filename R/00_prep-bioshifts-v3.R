## add climate velocities to bioshifts v3


## read in bioshifts v3
v3 = read.csv("data-raw/BIOSHIFTSv3/BIOSHIFTS_v3.csv")
v3$scientificName = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$scientificName

## fix ecosystem column 
## everything != marine is ter
v3$Eco = ifelse(!v3$Eco %in% c("Mar", "Ter"), "Ter", v3$Eco)

## except: A164_P1 = "Mar"
v3$Eco[which(v3$ID %in% c("A164_P1", "A183_P1"))] = "Mar"


##################################################
####   prep study-level climate velocity data ####
##################################################
## combine all the csv files from Brunno's calculations
## get names of all csv files in the climate velocity folder
files <- list.files(path = "~/Documents/bioshifts-traits/data-raw/bioshiftsv3/Velocity_SA", pattern="csv$", recursive=TRUE, full.names=TRUE)

## pick out columns to keep for terrestrial/marine study areas 
terr_cols <- c("ID", "v.lat.mean.mat", # mean latitudinal velocity of air temperature 
               "v.ele.mean.mat", # mean elevational velocity of air temperature
               "v.lat.sd.mat", # sd of both 
               "v.ele.sd.mat")
mar_cols <- c("ID", "v.lat.mean.sst",
              "v.lat.sd.sst")

cvs <- c()
for(f in 1:length(files)) {
  
  ## read in a file 
  file <- read.csv(files[f])
  
  ## if study is terrestrial, extract air temperature columns for ele and lat 
  if(any(str_detect(colnames(file), "mat"))) {
  
    ## if calcualted at 25 km 
    if(str_detect(files[f], "25km")) {
      file <- select(file, c("ID", "v.lat.mean.mat", "v.lat.sd.mat"))
      
      ## rearrange
      file <- gather(file, key = "cv_type", value = "val", 
                     c(v.lat.mean.mat, v.lat.sd.mat)) %>%
        mutate(Type = ifelse(str_detect(cv_type, "ele"), "ELE", 
                             "LAT"),
               Eco = "Ter") %>%
        mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                                "sd_cv_studylevel")) %>%
        select(-cv_type) %>%
        spread(key = measure, value = val) %>%
        mutate(cv_res = "25km")
    } 
    else {
      file <- select(file, terr_cols)
      ## rearrange
      file <- gather(file, key = "cv_type", value = "val", 
                     c(v.lat.mean.mat, v.ele.mean.mat, v.lat.sd.mat, v.ele.sd.mat)) %>%
        mutate(Type = ifelse(str_detect(cv_type, "ele"), "ELE", 
                             "LAT"),
               Eco = "Ter") %>%
        mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                                "sd_cv_studylevel")) %>%
        select(-cv_type) %>%
        spread(key = measure, value = val) %>%
        mutate(cv_res = "1km")
    }
    
    ## bind 
    cvs <- rbind(cvs, file)
    
  }
  else {
    file <- select(file, mar_cols)
    
    ## rearrange
    file <- file %>%
      mutate(Type = "LAT",
             Eco = "Mar") %>%
      rename("mean_cv_studylevel" = v.lat.mean.sst,
             "sd_cv_studylevel" = v.lat.sd.sst) %>%
      mutate(cv_res = "25km")
    
    ## bind 
    cvs <- rbind(cvs, file)
  }
  
  print(paste0("On file number ", f))
}

## if a latitudinal study has both a 25km and a 1km resolution latitudinal climate velocity, keep the 25km one
sub <- cvs %>%
  group_by(Type, ID) %>%
  mutate(n = n()) %>%
  filter(n == 2, cv_res == "25km") 

cvs <- cvs %>%
  group_by(Type, ID) %>%
  mutate(n = n()) %>%
  filter(n != 2) %>% 
  rbind(sub, .) %>%
  select(-n)

## save 
write.csv(cvs, "data-processed/v3_study-level-cvs.csv", row.names = FALSE)

  

## make sure all studies are there 
unique(v3$ID)[which(!unique(v3$ID) %in% cvs$ID)]
## missing some - Brunno said some too small 
## check back later 

## join to v3
join <- left_join(v3, cvs)
nrow(v3) == nrow(join) 
length(which(is.na(join$mean_cv_studylevel))) ## 49 do not have climate velocity 
## make sure these are the ones without csvs 
unique(join$ID[which(is.na(join$mean_cv_studylevel))]) %in% unique(v3$ID)[which(!unique(v3$ID) %in% cvs$ID)]
## one is not: B602_P1
## v3 says it is a marine study, but the climate velocities I have are terrestrial 
## ignore for now 


####################################################
####   prep species-level climate velocity data ####
####################################################
## add species-specific climate velocities 
## read in file with preliminary climate velocities 
cvs_spp <- read.csv("~/Documents/bioshifts-traits/data-processed/v3_new-cvs_preliminary.csv")
length(unique(cvs_spp$species_name))

cvs_spp <- cvs_spp %>%
  select(species_studyid, cv_type, range_source, 
         mean_cv_sppspecific, sd_cv_sppspecific, Type, ID) %>%
  mutate(scientificName = str_split_fixed(species_studyid, "\\_", 2)[,1]) %>%
  filter(!str_detect(cv_type, "map")) %>% ## get rid of precipitation
  distinct()

## do small test join
join <- join %>%
  select(ID, mean_cv_studylevel, sd_cv_studylevel, everything()) 

join <- left_join(join, cvs_spp, 
                  relationship = "many-to-many")

length(which(is.na(join$mean_cv_sppspecific))) ## 3724 shifts missing species-specific cv estimate

## check nothing got duplicated 
nrow(join) == nrow(v3)

## save 
write.csv(join, "data-processed/v3_with-cv.csv", row.names = FALSE)
