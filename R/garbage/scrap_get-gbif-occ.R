## get GBIF occurrences during study period 
library(tidyverse)
library(rgbif)
library(parallel)
library(data.table)

# set the settings
computer = "matrics"
work_dir <- getwd()
source("~/Documents/bioshifts-traits/R/calculating-climate-velocities/settings.R")
source("~/Documents/bioshifts-traits/R/calculating-climate-velocities/decompress_file.R")
range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))


####################################
#    GET GBIF CODES FOR SPECIES    #
####################################
## get list of species
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

## filter to leading edge shifts with positive climate velocity 
dd <- filter(dd, Param == "LE" & ClimVeloTKmY_spp >= 0)

missing <- select(dd, scientificName, class)

## get their gbif codes
codes <- read.csv("data-raw/bioshiftsv3/splist_v3.csv") %>%
  select(species, kingdom, phylum, class, order, family, db, db_code) %>%
  unique()

missing <- left_join(missing, codes, by = c("scientificName"= "species", "class" = "class")) %>%
  filter(!is.na(db_code)) %>%
  distinct()

missing$db_code <- gsub("GBIF:","", missing$db_code)

missing <- filter(missing, !str_detect(missing$db_code, "\\:"))


##############################
#    SUBMIT GBIF REQUESTS    #
##############################
## request occurrence downloads for each species 
requests <- data.frame()
inforequest <- list()
for(i in 96:nrow(missing)){ 
  
  print(paste0("Requesting species ", i))
  
  inforequest[[i]] <- occ_download(
    #pred("taxonKey", 8417931),
    pred_in("speciesKey", missing$db_code[i]),
    pred_in("basisOfRecord",
            basisOfRecord),
    pred_gte("year", range_env_data[1]),
    pred_lte("year", range_env_data[2]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("occurrenceStatus", "PRESENT"),
    format = "SIMPLE_CSV",
    user = Sys.getenv('GBIF_USER'), 
    pwd = Sys.getenv('GBIF_PWD'),
    email = "nicole.moore@mail.mcgill.ca")
  
  tmp <- inforequest[[i]]
  
  # save information about request
  tmp.save <- attributes(tmp)
  tmp.save <- data.frame(download_key = tmp[1],
                         created = tmp.save$created,
                         download_link = tmp.save$downloadLink,
                         doi = tmp.save$doi,
                         citation = tmp.save$citation,
                         format = tmp.save$format,
                         user = tmp.save$user,
                         email = tmp.save$email)
  requests <- rbind(requests,tmp.save)
  
  occ_download_wait(inforequest[[i]])
}

## check  
library(purrr)
df <-
  inforequest %>%
  map_df(
    ~ data.frame(
      download_key = .[1],
      created = attr(., "created"),
      citation = attr(., "citation"),
      format = attr(., "format"),
      user = attr(., "user"),
      email = attr(., "email"),
      doi    =  attr(., "doi"),
      download_link    =  attr(., "downloadLink")
    )
  )

requests <- requests[-which(!requests$doi %in% df$doi),]

requests <- rbind(requests, df[which(!df$doi %in% requests$doi),])

write.csv(requests, "data-processed/gbif_requests.csv", row.names = F)

############################################
#    DOWNLOAD OCCURRENCES FROM REQUESTS    #
############################################
# Download requests at the HPC ----
requests <- read.csv("data-processed/gbif_requests.csv")
head(requests)

GBIF_zip_dir <- "/Volumes/NIKKI/Dispersal_GBIF/GBIF_data_NM"
if(!dir.exists(GBIF_zip_dir)){
  dir.create(GBIF_zip_dir,recursive = T)
}
ncores = parallelly::availableCores()

check <- "Error"
attempt <- 0
errors <- TRUE

keystogo <- requests$download_key
while( any(errors) ) {
  attempt <- attempt + 1
  
  cat("Attempt", attempt, "\n")
  
  check <- mclapply(1:length(keystogo), function(i){
    
    keystogo_i = as.character(keystogo[i])
    
    test = occ_download_meta(keystogo_i)
    test = test$status
    
    while(test == "RUNNING"){
      
      test = occ_download_meta(keystogo_i)
      test = test$status
      
      Sys.sleep(60)
      
    }
    
    check[[i]] <- try (
      {
        rgbif::occ_download_get(key = keystogo_i,
                                path = GBIF_zip_dir,
                                overwrite = TRUE)
        return(paste("OK", i))
      },
      silent = TRUE
    )
    
  }
  , mc.cores = ifelse(ncores>length(keystogo),length(keystogo),ncores))
  
  errors <- sapply(check, function(x) class(x)=="try-error")
  
  if(any(errors)){
    errors <- requests$download_key[which(errors)]
    keystogo <- errors
  }
  cat("Error on keys:", errors)
} 


#################################################
#    TEST IF EVERYTHING DOWNLOADED CORRECTLY    #
#################################################
# Test if everything downloaded correctly
# Compare file sizes of remote and local files
download_size <- function(url) as.numeric(httr::HEAD(url)$headers$`content-length`)

for(i in 1:length(keystogo)){ cat("\rChecking file", i, "from", length(keystogo))
  # local file size
  size_local <- file.size(here::here(GBIF_zip_dir, paste0(keystogo[i],".zip")))
  # remote file size
  size_remote <- download_size(requests$download_link[i])
  # download again if files dont have the same size
  if(!size_local == size_remote){
    cat("\rFiles don't have the same size\nDownloading file", i)
    rgbif::occ_download_get(key = keystogo[i],
                            path = GBIF_zip_dir,
                            overwrite = TRUE)
  }
}

