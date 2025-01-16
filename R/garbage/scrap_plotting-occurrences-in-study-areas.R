## delineating historical range edges from GBIF data 
library(tidyverse)
library(rgbif)
library(parallel)
library(data.table)
library(sf)
library(cowplot)

# set the settings
computer = "matrics"
work_dir <- getwd()
source("~/Documents/bioshifts-traits/R/calculating-climate-velocities/settings.R")
source("~/Documents/bioshifts-traits/R/calculating-climate-velocities/decompress_file.R")
range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))


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

#############################
## filter the occurrences  ##
#############################
## add ecosystem 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$scientificName = str_replace_all(v3$sp_name_checked, "\\_", " ")

v3 <- v3 %>%
  mutate(Eco = ifelse(Eco %in% c("Mar"), "Mar", "Ter")) %>%
  select(scientificName, Eco) %>%
  distinct()

missing <- left_join(missing, v3)

GBIF_zip_dir <- "/Volumes/NIKKI/Dispersal_GBIF/GBIF_data_NM"

myselection <- c("basisOfRecord", "speciesKey", "decimalLongitude", "decimalLatitude", "year", "month","species")
basisOfRecord_selected <- basisOfRecord

terrestrials <- missing$species_name[which(missing$Eco == "Ter")]
marines <- missing$species_name[which(missing$Eco == "Mar")]

# my terrestrial raster
ter.ras <- terra::rast("~/Documents/bioshifts-traits/data-raw/model_raster_ter_1km.tif")
# my marine raster
mar.ras <- terra::rast("~/Documents/bioshifts-traits/data-raw/model_raster_mar.tif")

# test
any(terrestrials %in% marines)

# zipped occ files
occs <- list.files(here::here(GBIF_zip_dir), pattern = ".zip")

# create dir to save sps occurrences
occ_dir = "/Volumes/NIKKI/Dispersal_GBIF/GBIF_occurrences_dispersal"
if(!dir.exists(occ_dir)){
  dir.create(occ_dir, recursive = T)
}

# create temp dir to decompress zip file
tmp.dir <- here::here(occ_dir,"tmp")
if(!dir.exists(tmp.dir)){
  dir.create(tmp.dir)
}


ncores = parallelly::availableCores()


i = 1
while(i <= length(occs)) {
  
  # decompress zip file
  zippedfile <- here::here(GBIF_zip_dir, occs[i])
  decompress_file(directory = tmp.dir, file = zippedfile)
  unzippedfile <- gsub(".zip",".csv",occs[i])
  unzippedfile <- gsub(".zip",".csv",here::here(tmp.dir,unzippedfile))
  
  # read in
  tmp <- fread(unzippedfile, select = myselection, nThread = 10, quote="")
  tmp <- na.omit(tmp) 
  
  # filter basisOfRecord
  tmp <- tmp %>% dplyr::filter(basisOfRecord %in% basisOfRecord_selected)
  
  # get cells
  cells. <-  terra::extract(ter.ras,
                            tmp[,c("decimalLongitude", "decimalLatitude")],
                            cells=TRUE)
  tmp$cell <- cells.$cell
  tmp$layer <- cells.[,2]
  # remove cells falling in the ocean
  tmp <- tmp[which(tmp$layer==1),] # keep land / remove ocean
  tmp = tmp[,-"layer"]
  
  # Remove duplicates: same species in the same location and date
  rm <- duplicated(tmp[,c("species", "cell", "year", "month")])
  if(any(rm)){
    tmp <- tmp[-which(rm),]
  }
  
  # Remove other potential issues with the package CoordinateCleaner
  tmp = CoordinateCleaner::clean_coordinates(
    x = tmp,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    tests = c("capitals", "centroids", "equal",
              "gbif", "institutions", "zeros"),
    value = "clean")

  ## filter dispersal data to studies with this species
  dd_sp <- filter(dd, scientificName %in% tmp$species) %>%
    select(scientificName, ID, START_firstperiod, END_firstperiod, START_secondperiod, END_secondperiod, 
           ShiftKmY) %>%
    distinct()
  
  ## loop through studies 
  for(l in 1:nrow(dd_sp)) {
    year1 <- dd$START_firstperiod[l]
    year2 <- dd$END_firstperiod[l]
    
    year3 <- dd$START_secondperiod[l]
    year4 <- dd$END_secondperiod[l]
    
    ## filter occurrences to study beginning and end of study period 
    occ_firstper <- filter(tmp,  year >= year1 & year <= year2)
    occ_secondper <- filter(tmp,  year >= year3 & year <= year4)
    
    ## plot with study area shapefile 
    file = paste0("~/Documents/bioshifts-traits/data-raw/bioshiftsv3/ShapefilesBioShiftsv3/", dd_sp$ID[l], ".shp")
    study_area = st_read(file)
    
    st_firstper <- st_multipoint(as.matrix(occ_firstper[,3:4])) %>%
      st_sfc()
    st_crs(st_firstper) = st_crs(study_area)
    
    st_secondper <- st_multipoint(as.matrix(occ_secondper[,3:4])) %>%
      st_sfc()
    st_crs(st_secondper) = st_crs(study_area)
    
    ## crop to occurrences within study area
    st_firstper_crop <- st_intersection(st_firstper, study_area)
    st_secondper_crop <- st_intersection(st_secondper, study_area)
    
    # ## first period
    # study_area %>%
    #   ggplot() +
    #   geom_sf(fill = "transparent") +
    #   geom_sf(data = st_firstper, size = 1)
    # 
    # ## second period 
    # study_area %>%
    #   ggplot() +
    #   geom_sf(fill = "transparent") +
    #   geom_sf(data = st_secondper, size = 1)
    
    ## first period
    p1 <- study_area %>%
      ggplot() +
      geom_sf(fill = "transparent") +
      geom_sf(data = st_firstper_crop, size = 1) +
      labs(title = dd_sp$scientificName[l], subtitle = paste(year1, "-", year2), 
           caption =  paste(dd_sp$ShiftKmY[l], "km/y"))
    
    ## second period 
    p2 <- study_area %>%
      ggplot() +
      geom_sf(fill = "transparent") +
      geom_sf(data = st_secondper_crop, size = 1) +
      labs(subtitle = paste(year3, "-", year4)) 
    
    ## arrange in grid 
    p <- plot_grid(p1, p2, align = "h")
    
    ## save plot
    ggsave(p, path = "figures/study-area_occurrences", 
           filename = paste0(dd_sp$ID[l], "_", dd_sp$scientificName[l], ".png"),
           width = 8, height = 4)
  }
  
  i = i + 1
}


# delete tmp dir used to decompress zipfiles
unlink(tmp.dir, recursive = TRUE)
unlink(GBIF_zip_dir, recursive = TRUE)
