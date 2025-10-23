## trying out Jonathan's methodological figure idea
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
library(ggplot2)
library(stringr)
theme_set(theme_bw())

## start by reading in data
data = read.csv("data-processed/model-data_main.csv")

## now read in species-specific study polygons 
polys = st_read("data-raw/BIOSHIFTSv3/ShapefilesBioShiftsv3_species.gpkg")

## and full study polygons 
studies = st_read("data-raw/BIOSHIFTSv3/ShapefilesBioShiftsv3_studyareas.gpkg")
                
## filter to specices x study area combinations in the dispersal data 
polys = polys %>%
  filter(ID %in% paste(data$sp_name_checked, data$ID, sep = "_"))

nrow(polys) ## 406 
length(unique(paste(data$sp_name_checked, data$ID, sep = "_"))) ## 406 - all are there 

studies = studies %>%
  filter(Name %in% data$ID)
nrow(studies) ## 24 

## plot the polygons on a map:
## first, rasterize: 
## create template to define the extent and resolution of the raster 
r = rast(vect(polys), res=0.5)

## rasterize 
polys_raster <- terra::rasterize(vect(polys), r, field = 1, fun = "sum")
crs(polys_raster) <- "EPSG:4326" # if missing
polys_raster

my_breaks <- c(1, 10, 100, 3000)

## add group
group_key = data %>%
  select(group, sp_name_checked) %>%
  unique()

polys = polys %>%
  mutate(sp_name_checked = str_replace_all(Species_sdt, " ", "_")) %>%
  left_join(., group_key) %>%
  mutate(studyID = str_split_fixed(ID, "\\_", 3)[,3]) 

## crop countries to northern hemisphere 
world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
northern_hemisphere_map <- world_map[st_coordinates(st_centroid(world_map))[, "Y"] > 0, ]

## isolate the china study polygon
china = filter(polys, str_detect(ID, "A145_P1"))

all_polys = ggplot() + 
  geom_spatraster(data = polys_raster) + 
  geom_sf(data = world_map, fill = "grey90", colour = "gray70", linewidth = 0.1) + 
  geom_spatraster(data = polys_raster) +
  coord_sf() +
  scale_fill_viridis_c(option = "plasma", 
                       na.value = "transparent", 
                       alpha = 0.5,
                       name = "No. species") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-155, 136)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,80),
                     breaks = c(0, 30, 60, 90)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "gray98"), 
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "bottom") +
  ## add outlines of clipped study polys 
  geom_sf(data = polys, fill = "transparent", colour = "black", linewidth = 0.1, linetype = "longdash",
          inherit.aes = F) +
  ## add study area polygon outlines
  geom_sf(data = studies, fill = "transparent", colour = "black", linewidth = 0.1, inherit.aes = F) +
  ## add china studies with a base colour
  geom_sf(data = china, fill = "#7976b6", colour = "#7976b6", linewidth = 0.1, linetype = "longdash",
          inherit.aes = F) +
  ## add china studies with a base colour
  geom_sf(data = china, fill = "transparent", colour = "black", linewidth = 0.15, linetype = "longdash",
          inherit.aes = F) 

## save:
ggsave(all_polys, path = "figures/methodological-figure", filename = "map_all-polys.png", width = 7.8, height = 4)

## make legend 
legend_plot = ggplot(world_map) + 
  geom_sf(fill = "grey90", linewidth = 0.1) + 
  coord_sf() +
  geom_spatraster(data = polys_raster) + 
  scale_x_continuous(expand = c(0,0),
                     limits = c(-155, 136)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,80),
                     breaks = c(0, 30, 60, 90)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "gray98"), 
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "bottom") +
  ## add study area polygon outlines
  geom_sf(data = filter(studies, Name %in% c("A129_P1", "A156_P1")), aes(linetype = Name), 
          fill = "transparent", colour = "black", linewidth = 0.5, inherit.aes = F) +
  scale_linetype_manual(values = c("solid", "longdash"), labels = c("Study area", "Species range\nwithin study area")) +
  guides(fill = "none") +
  labs(linetype = "")

ggsave(legend_plot, path = "figures/methodological-figure", filename = "map_legend.png", 
       width = 4, height = 5)

## choose one study with a few species to zoom in on 
polys %>% 
  as.data.frame() %>%
  mutate(studyID = str_split_fixed(ID, "\\_", 3)[,3]) %>% 
  group_by(studyID) %>%
  count(studyID)

data %>% 
  group_by(ID, group) %>%
  count(ID) %>% view()

## zoom in on a plant and a bird:
candidates = polys %>% 
  filter(studyID %in% c("A129_P1", "A129_P3", "A32_P1", "A181_P1", "A85_P1"))

candidates = polys %>% 
  filter(studyID %in% c("A156_P1"))

## Bird in Britain and plant in Quebec 
#########################################
## Carduelis_carduelis - common goldfinch
## Acer rubrum - red maple
## save two separate insets 

## european goldfinch
one = ggplot(world_map) + 
  geom_sf(fill = "grey90", colour = "gray70", linewidth = 0.1) + 
  coord_sf() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-8, 2)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(49.8,59)) +
  theme(panel.grid = element_blank()) +
  ## add outlines of clipped study polys 
  geom_sf(data = filter(polys, ID %in% c("Carduelis_carduelis_A129_P1")), colour = "black",
          aes(colour = Species_sdt, fill = Species_sdt), linewidth = 0.2, 
          alpha = 0.5, inherit.aes = F,  linetype = "longdash") +
  theme(panel.background = element_rect(fill = "gray98"), legend.position = "none") +
  labs(colour = "", fill = "") +
  scale_fill_manual(values = c("#e49e00")) + 
  scale_colour_manual(values = c("#e49e00")) +
  ## add study area polygon outlines
  geom_sf(data = filter(studies, Name %in% c("A129_P1")), fill = "transparent", colour = "black", 
          linewidth = 0.2, inherit.aes = F) 

two = ggplot(world_map) + 
  geom_sf(fill = "grey90", colour = "gray70", linewidth = 0.1) + 
  coord_sf() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-80.2, -59.1)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(43.5,55.7)) +
  theme(panel.grid = element_blank()) +
  ## add outlines of clipped study polys 
  geom_sf(data = filter(polys, ID %in% c("Acer_rubrum_A85_P2")), colour = "black",
          aes(colour = Species_sdt, fill = Species_sdt),  linetype = "longdash",
          linewidth = 0.2, alpha = 0.5, inherit.aes = F) +
  theme(panel.background = element_rect(fill = "gray98"), legend.position = "none") +
  labs(colour = "", fill = "") +
  scale_fill_manual(values = c("#A2B06D")) + 
  scale_colour_manual(values = c("#A2B06D")) +
  ## add study area polygon outlines
  geom_sf(data = filter(studies, Name %in% c("A85_P1")), fill = "transparent", colour = "black", 
          linewidth = 0.2, inherit.aes = F) 

ggsave(one, path = "figures/methodological-figure", filename = "map_goldfinch.png", width = 3, height = 3)
ggsave(two, path = "figures/methodological-figure", filename = "map_redmaple.png", width = 4, height = 3)

## get the data
stats = data %>%
  filter(sp_name_checked %in% c("Carduelis_carduelis", "Acer_rubrum"), 
         ID %in% c("A129_P1", "A85_P1")) %>%
  select(sp_name_checked, ShiftKmY, DispersalPotentialKmY, MedianDispersalPotentialKmY, 
         ClimVeloKmY_RelScale, q3ClimVeloKmY_RelScale)

stats[,2:6] = round(stats[,2:6], digits = 1)  
stats

## write
write.csv(stats, "figures/methodological-figure/stats-for-methodological-figure.csv", row.names = F)
