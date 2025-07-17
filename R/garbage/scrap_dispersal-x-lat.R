## look at dispersal speed versus latitude 
library(tidyverse)
library(sf)

## read in dispersal dataset
v3 <- read.csv("data-processed/v3_potential-dispersal-rate.csv")
v3$sp_name_checked = str_replace_all(v3$sp_name_checked, "\\_", " ")

## read in range polygons 
ranges = readRDS("~/Documents/bioshifts-traits/data-processed/large-data/collated-ranges.rds")
colnames(ranges)

## make sure one range per species
range_sources <- c("IUCN", "BOTW", "GARD", "BIEN", "Fishmap", "Aquamaps", "Butterfly Maps", "GIFT", "GBIF occurrence")

ranges <- ranges %>%
  mutate(range_source = factor(.$range_source, ordered = TRUE, 
                               levels = range_sources)) %>%
  group_by(binomial) %>% ## group by species_studyid
  arrange(range_source, .by_group = TRUE) %>% ## sort from high to low priority within group
  filter(row_number() == 1) ## keep first row of each group (highest priority)

length(unique(ranges$binomial)) == nrow(ranges)

## get latitude of the centroid of each range 
ranges <- ranges %>%
  mutate(lat_centroid = map_dbl(geometry, ~st_point_on_surface(.x)[[2]])) 

ranges$lat_min = as.vector(lapply(1:nrow(ranges), FUN = function(x) {return(st_bbox(ranges[x,])$ymin)}))
ranges$lat_min = as.numeric(as.character(unlist(ranges$lat_min)))
ranges$lat_max = as.vector(lapply(1:nrow(ranges), FUN = function(x) {return(st_bbox(ranges[x,])$ymax)}))
ranges$lat_max = as.numeric(as.character(unlist(ranges$lat_max)))


## plot to make sure:
ranges %>%
  filter(binomial == ranges$binomial[10]) %>%
  ggplot() +
  geom_sf() +
  geom_hline(aes(yintercept = lat_centroid)) +
  geom_hline(aes(yintercept = lat_max), colour = "blue") +
  geom_hline(aes(yintercept = lat_min), colour = "red") 

## join to dispersal data
v3_sub <- inner_join(ranges, v3, by = c("binomial" = "sp_name_checked"))

## plot across latitude 
v3_sub %>%
  ggplot(aes(x = abs(lat_max - lat_min), y = DispersalPotentialKmY, colour = group)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10()


