## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

# dd$ClimVeloTKmY_study = dd$ClimVeloTKmY_spp
# dd = filter(dd, !is.na(ClimVeloTKmY_study))

## filter to leading edge with + climate velocity
dd <- filter(dd, Param == "LE" & ClimVeloTKmY_spp >= 0)

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_spp)) %>%
  mutate(ClimVeloTKmY_spp = abs(ClimVeloTKmY_spp)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_spp, NA)) %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

## plot
dd %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) 
+
  facet_wrap(~is_contraction) 

dd %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Climate velocity (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~is_contraction) 

dd %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Climate velocity (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~is_contraction) 

dd %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = N_periodes)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_hline(yintercept = 0) 



## get rid of contractions
dd = filter(dd, is_contraction %in% c("UNKNOWN", "NO"))

## filter to expansions in same direction as cv 
dd <- filter(dd, tracking_climate == TRUE)

## standardize so CV > 0 means away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ShiftKmY = abs(dd$ShiftKmY)




dd %>%
  filter(is_contraction == "YES") %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Minimum of rate of climate change and\npotential dispersal rate (km/y)",
       y = "Range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(aes(group = is_contraction), method = "lm") +
  facet_wrap(~N_periodes) 




## plot climate velocity raster for  A10_P1
library(terra)
rast <- rast("~/Documents/bioshifts-traits/data-raw/bioshiftsv3/Velocity_SA/A10_P1_map_25km_gVelLat.tif")

rast

rast2 <- rast$Vel > 0

ggplot() +
  tidyterra::geom_spatraster(data = rast) 


dd %>%
   arrange(Rate) %>%
   select(Rate, ID, scientificName, everything()) %>%
   View




join %>%
  filter(ID == "A10_P1") %>%
  filter(Param == "LE") %>%
  arrange(-Rate) %>%
  View()
  
  
  
join %>%
  filter(ID == "A10_P1") %>%
  filter(Param == "LE") %>%
  ggplot(aes(x = Rate)) +
  geom_histogram() +
  geom_vline(xintercept = 0)




## get 95th and 5th perctile (most extreme shifts and contractions)
quant <- quantile(dd$Rate, p  = c(0.05,0.95))

dd_sub <- filter(dd, Rate <= quant[1] | Rate >= quant[2])

dd = dd %>%
  mutate(quant = ifelse(Rate <= quant[1], "5th", 
                        ifelse(Rate >= quant[2], "95th", 
                               "Not in the quantile")) )



## plot with 5th and 95th quantile coloured differently 
dd %>%
  filter(quant!= "5th") %>%
  filter(is_contraction == "NO") %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Limiting rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  #scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "lm")

dd %>%
  filter(ID == "A138_P1") %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-26, 40), expand = c(0,0)) +
  labs(x = "Limiting rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  #scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  geom_smooth(method = "lm")


dd %>%
  ggplot(aes(x = Rate)) +
  geom_histogram()


View(dd[which(sign(dd$mean_cv_sppspecific) != sign(dd$mean_cv_studylevel)),])

