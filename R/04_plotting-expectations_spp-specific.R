## plotting predictions and data
## script and analysis for species on the move presentation (May 2023)
library(cowplot)
library(lme4)
library(MuMIn)

## read in data
v3 <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter 
v3 <- v3 %>%
  ## filter to cases where we expect dispersal to matter (range expansion): 
  ## leading edge shifts where climate velocity is positive
  ## trailing edge shifts where climate velocity is negative 
  ## centroid shifts
  filter(Param == "TE" & mean_cv_studylevel <= 0 |
           Param == "LE" & mean_cv_studylevel >= 0 |
           Param == "O") 

## get rid of non-expansions 
v3 = filter(v3, is_contraction %in% c("UNKNOWN", "NO"))

## get rid of centroid shifts in opposite direction to cv 
v3 = filter(v3, !(Param == "0" & !sign(ClimVeloTKmY_spp) == sign(ShiftKmY)))

## standardize so CV > 0 meanings away from range centre 
v3$ClimVeloTKmY_spp = abs(v3$ClimVeloTKmY_spp)
v3$ShiftKmY = abs(v3$ShiftKmY)

## add gradient 
v3 <- v3 %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

## get rid of shifts that don't have a spp-specific climate velocity
v3 <- filter(v3, !is.na(ClimVeloTKmY_spp))

hist(v3$ClimVeloTKmY_study)
hist(v3$ClimVeloTKmY_spp)

##################################
##       plot expectations      ##
##################################
hist(v3$ClimVeloTKmY_spp)

## rug data
rug <- v3 %>%
  mutate(x = 0, y = 0)

## climate velocity 
rug %>%
  ggplot(aes(x = x, y = ClimVeloTKmY_spp, colour = Type)) +
  geom_point(shape = "_", size = 8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 7), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)")

## make histograms
v3 %>%
  ggplot(aes(x = ClimVeloTKmY_spp, fill = ..x..)) +
  geom_histogram() + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-1, 15), expand = c(0,0)) +
  labs(x = "Mean climate velocity across study area (km/y)",
       y = "Number of range expansion observations") +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_fill_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 4.5) 

ggsave(path = "figures/sotm", filename = "climate-velo-hists_spp.png", 
       device = "png", height = 4, width = 8)


## make new dataframe for expectation plot
clim_velos <- c(v3$ClimVeloTKmY_spp)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 5, by = 0.01)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 5, by = 0.01)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 5, by = 0.01))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.01)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(v3), ClimVeloTKmY_spp, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY_spp"))

one <- clim_data %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

one_legend <- get_legend(one)

one <- one +
  theme(legend.position = "none")

ggsave(one, path = "figures/sotm", filename = "expectation-climate-1_spp.png", 
       device = "png", height = 4, width = 5)
ggsave(one_legend, path = "figures/sotm", filename = "expectation-climate-1-legend_spp.png", 
       device = "png", height = 2, width = 1)



## split by latitude versus elevation
two = clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

two_legend = get_legend(two)

two <- two +
  theme(legend.position = "none")

ggsave(two, path = "figures/sotm", filename = "expectation-climate-2_spp.png", 
       device = "png", height = 4, width = 8)
ggsave(two_legend, path = "figures/sotm", filename = "expectation-climate-2-legend_spp.png", 
       device = "png", height = 2.4, width = 1.5)

## zoom in
filtered_lags <- filter(v3, ClimVeloTKmY_spp < 0.5)
clim_velos <- filtered_lags$ClimVeloTKmY_spp

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1, by = 0.0001)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1, by = 0.0001)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1, by = 0.0001))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.0001)) 
  return(rep(x, reps))
}))

clim_data_zoom <- data.frame(y_values = append(y_values_11, y_values), 
                             x_values = append(x_values_11, x_values), 
                             velos = append(velos, y_values))

gradients <- select(ungroup(filtered_lags), ClimVeloTKmY_spp, Gradient) %>%
  distinct()

clim_data_zoom <- left_join(clim_data_zoom, gradients, by = c("velos" = "ClimVeloTKmY_spp"))

three = clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  labs(x = "",
       y = "", 
       colour = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

three_legend = get_legend(three)

three <- three +
  theme(legend.position = "none")

ggsave(three, path = "figures/sotm", filename = "expectation-climate-3_spp.png", 
       device = "png",  height = 2.8, width = 3.1)
ggsave(three, path = "figures/sotm", filename = "expectation-climate-3-large_spp.png", 
       device = "png", height = 4, width = 4.2)
ggsave(three_legend, path = "figures/sotm", filename = "expectation-climate-3-legend_spp.png", 
       device = "png", height = 1, width = 1.5)


## now split by latitude and elevation
ele_split <- clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none")

ggsave(ele_split, path = "figures/sotm", filename = "expectation-elevation_spp.png", 
       device = "png", height = 4, width = 4.2)

lat_split = clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

ggsave(lat_split, path = "figures/sotm", filename = "expectation-latitude_spp.png", 
       device = "png", height = 4, width = 8)


## add dispersal potential rug plot 
rug_lat <- filter(rug, Gradient == "Latitudinal") %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE))
rug_ele <- filter(rug, Gradient == "Elevation") %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE))

rugplot_lat <- clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

ggsave(rugplot_lat, path = "figures/sotm", filename = "expectation-latitude-rug-zoom_spp.png", 
       device = "png", height = 4, width = 8)

rugplot_ele_zoom_squish <- clim_data %>%
  #filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                             add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_zoom_squish, path = "figures/sotm", filename = "expectation-elevation-rug_spp.png", 
       device = "png",  height = 2.8, width = 3.1)


rugplot_ele <- clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele, path = "figures/sotm", filename = "expectation-elevation-rug-zoom_spp.png", 
       device = "png",  height = 2.8, width = 3.1)

rugplot_lat_inset <- clim_data %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0,5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0))+
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_lat, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_lat_inset, path = "figures/sotm", filename = "expectation-latitude-rug-zoom_spp.png", 
       device = "png",  height = 2.8, width = 3)

## zoom out
## make new dataframe for expectation plot
clim_velos <- c(v3$ClimVeloTKmY_spp)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1400, by = 1)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 5, by = 0.01)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1400, by = 1))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.01)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(v3), ClimVeloTKmY_spp, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY_spp"))


rugplot_lat_unzoom <- clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                             add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

ggsave(rugplot_lat_unzoom, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom_spp.png", 
       device = "png", height = 4, width = 8)

rugplot_ele_unzoom <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                             add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_unzoom, path = "figures/sotm", filename = "expectation-elevation-rug_spp.png", 
       device = "png", height = 4, width = 4.2)

clim_data_test <- clim_data %>%
  arrange(velos, x_values, y_values)

## make log version
rugplot_lat_unzoom_log <- clim_data %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1, linejoin = "mitre") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), 
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") 

ggsave(rugplot_lat_unzoom_log, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom-log_spp.png", 
       device = "png", height = 4, width = 8)

rugplot_lat_unzoom_log_line <- clim_data %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1, linejoin = "mitre") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), 
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 5), expand = c(0,0)) +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  geom_vline(aes(xintercept = max(velos)))

ggsave(rugplot_lat_unzoom_log_line, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom-log-line_spp.png", 
       device = "png", height = 4, width = 8)


clim_velos <- c(v3$ClimVeloTKmY_spp)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1400, by = 1)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 5, by = 0.0001)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1400, by = 1))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.0001)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(v3), ClimVeloTKmY_spp, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY_spp"))


rugplot_ele_unzoom_log <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  ggplot(aes(x = x_values, y = y_values, colour = velos, group = velos)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), expand = c(0,0),
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_unzoom_log, path = "figures/sotm", filename = "expectation-elevation-rug-log_spp.png", 
       device = "png", height = 2.8, width = 3.1)

rugplot_ele_unzoom_log_line <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  ggplot(aes(x = x_values, y = y_values, colour = velos, group = velos)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), expand = c(0,0),
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_vline(aes(xintercept = max(velos)))

ggsave(rugplot_ele_unzoom_log_line, path = "figures/sotm", filename = "expectation-elevation-rug-log-line_spp.png", 
       device = "png", height = 4, width = 4.2)


## compare distributions directly 
dens_grad <- v3 %>%
  gather(key = "Measure", value = "Measurement", c(AnnualDispPotKmY, ClimVeloTKmY_spp)) %>%
  ggplot(aes(x = Measurement, fill = Measure)) +
  geom_density(alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Gradient) +
  scale_x_log10() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.8)) +
  scale_fill_discrete(labels = c("Annual dispersal potential",
                                 "Climate velocity")) +
  labs(x = "", y = "Density", fill = "")


ggsave(dens_grad, path = "figures/sotm", filename = "density_gradient_spp.png", 
       device = "png", height = 3.5, width = 8.5)


dens_taxa <- v3 %>%
  gather(key = "Measure", value = "Measurement", c(AnnualDispPotKmY, ClimVeloTKmY_spp)) %>%
  ggplot(aes(x = Measurement, fill = Measure)) +
  geom_density(alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(group~Gradient) +
  scale_x_log10() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2)) +
  scale_fill_discrete(labels = c("Annual dispersal potential",
                                 "Climate velocity")) +
  labs(x = "", y = "Density", fill = "")

ggsave(dens_taxa, path = "figures/sotm", filename = "density_taxa_spp.png", 
       device = "png", height = 3.5, width = 8.5)

cv_vs_dp <- v3 %>%
  ggplot(aes(y = ClimVeloTKmY_spp, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp, path = "figures/sotm", filename = "cv_vs_dp_spp.png", 
       device = "png", height = 3.5, width = 4.5)


cv_vs_dp_zoom <- v3 %>%
  ggplot(aes(y = ClimVeloTKmY_spp, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  scale_x_continuous(limits = c(0, 7)) +  
  scale_y_continuous(limits = c(0, 7)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp_zoom, path = "figures/sotm", filename = "cv_vs_dp_zoom_spp.png", 
       device = "png", height = 3.5, width = 4.5)

cv_vs_dp_zoomzoom <- v3 %>%
  ggplot(aes(y = ClimVeloTKmY_spp, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  scale_x_continuous(limits = c(0, 0.015)) +  
  scale_y_continuous(limits = c(0, 0.015)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp_zoomzoom, path = "figures/sotm", filename = "cv_vs_dp_zoomzoom_spp.png", 
       device = "png", height = 3.5, width = 4.5)


## make a bar plot
bar_group <- v3 %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY_spp <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  mutate(expect_limitation = factor(expect_limitation, levels = c("Yes", "No"), 
                                    ordered = TRUE)) %>%
  ggplot(aes(x = expect_limitation, fill = Gradient)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(group~Gradient) +
  scale_fill_manual(values = c("#9FB798", "#28587B")) +
  labs(x = "Expect dispersal limitation?", y = "Number of range expansion observations") +
  theme(legend.position = "none")

ggsave(bar_group, path = "figures/sotm", filename = "bar_group_spp.png", 
       device = "png", height = 3.5, width = 4.5)

## get numbers 
v3 %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY_spp <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  mutate(expect_limitation = factor(expect_limitation, levels = c("Yes", "No"), 
                                    ordered = TRUE)) %>%
  group_by(group, Gradient, expect_limitation) %>%
  tally()

bar_gradient <- v3 %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY_spp <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  ggplot(aes(x = expect_limitation, fill = Gradient)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(~Gradient) +
  scale_fill_manual(values = c("#28587B", "#9FB798")) +
  labs(x = "Expect dispersal limitation?", y = "Number of observed expansions") +
  theme(legend.position = "none")

ggsave(bar_gradient, path = "figures/sotm", filename = "bar_gradient_spp.png", 
       device = "png", height = 3.5, width = 4.5)


##################################
##         plot the data        ##
##################################
## now plot the real points onto these plots 
lat <- filter(v3, Gradient == "Latitudinal") %>%
  mutate(Group = factor(group, ordered = F))

lat$ClimVeloTKmY_og <- lat$ClimVeloTKmY_spp
hist(lat$ClimVeloTKmY_spp)
q = quantile(lat$ClimVeloTKmY_spp, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
lat$ClimVeloTKmY_spp = cut(lat$ClimVeloTKmY_spp,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY_spp)

lat$ClimVeloTKmY_spp <- str_replace_all(lat$ClimVeloTKmY_spp, "\\[", "(") 
lat$ClimVeloTKmY_spp <- str_replace_all(lat$ClimVeloTKmY_spp, "\\]", ")") 

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)


#### hypothesis testing 
data_unlog <- v3 %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                        add = c(20, 0)),
                     limits = c(0, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none")+
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

ggsave(data_unlog, path = "figures/sotm", filename = "data-unlog_spp.png", 
       device = "png", height = 4, width = 8)

data_log <- v3 %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

ggsave(data_log, path = "figures/sotm", filename = "data-log_spp.png", 
       device = "png", height = 4, width = 8)


## make elevation one with diff scale 
data_ele_log <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5)  + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) 

ggsave(data_ele_log, path = "figures/sotm", filename = "data-ele-log_spp.png", 
       device = "png", height = 2.8, width = 3.1)


## now, plot line representing max climate velocity across elev and lat:
data_log <- data_log +
  geom_vline(aes(xintercept = max(ClimVeloTKmY_spp)))

ggsave(data_log, path = "figures/sotm", filename = "data-log-line_spp.png", 
       device = "png", height = 4, width = 8)

data_ele_log <- data_ele_log +
  geom_vline(aes(xintercept = max(ClimVeloTKmY_spp)))

ggsave(data_ele_log, path = "figures/sotm", filename = "data-ele-log-line_spp.png", 
       device = "png", height = 4, width = 4.2)


v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY_spp) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  theme(panel.grid = element_blank()) +
  # scale_x_continuous(expand = expansion(mult = c(0, 0), 
  #                                       add = c(20, 0)),
  #                    limits = c(0, 1400)) +
  scale_x_log10(limits = c(0.001, 1400)) +
  scale_y_log10() +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0)) +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#9FB798")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_vline(aes(xintercept = max(ClimVeloTKmY_spp)))

ggsave(path = "figures/sotm", filename = "data-lat-log-gt_spp.png", 
       device = "png", height = 3.5, width = 4.5)


## model data 
notlimited_lat <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY_spp) %>%
  mutate(ClimVeloTKmY_spp = factor(ClimVeloTKmY_spp, levels = unique(.$ClimVeloTKmY_spp),
                               ordered = TRUE))

limited_lat <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY_spp)  %>%
  mutate(ClimVeloTKmY_spp = factor(ClimVeloTKmY_spp, levels = unique(.$ClimVeloTKmY_spp),
                               ordered = TRUE))


library(nlme)
mod_nl_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY_spp, 
                  data = notlimited_lat)
summary(mod_nl_lat)

R2m=r.squaredGLMM(mod_nl_lat)[[1]] 
R2c=r.squaredGLMM(mod_nl_lat)[[2]]
R2m
R2c

mod_l_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY_spp, 
                 data = limited_lat)
summary(mod_l_lat)

R2m=r.squaredGLMM(mod_l_lat)[[1]] 
R2c=r.squaredGLMM(mod_l_lat)[[2]]
R2m
R2c

mod_l_lat_lm <- lm(ShiftKmY ~ AnnualDispPotKmY, 
                   data = limited_lat)
summary(mod_l_lat_lm)

R2m=r.squaredGLMM(mod_l_lat_lm)[[1]] 
R2m

AIC(mod_l_lat_lm, 
    mod_l_lat) ## this is better fit 

## elevation
notlimited_ele <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY_spp) %>%
  mutate(ClimVeloTKmY_spp = factor(ClimVeloTKmY_spp, levels = unique(.$ClimVeloTKmY_spp),
                               ordered = TRUE))

limited_ele <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY_spp)  %>%
  mutate(ClimVeloTKmY_spp = factor(ClimVeloTKmY_spp, levels = unique(.$ClimVeloTKmY_spp),
                               ordered = TRUE))


library(nlme) 
mod_nl_ele <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY_spp,
                  data = notlimited_ele)
summary(mod_nl_ele)

R2m=r.squaredGLMM(mod_nl_ele)[[1]]
R2c=r.squaredGLMM(mod_nl_ele)[[2]]
R2m
R2c

mod_l_ele <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY_spp,
                 data = limited_ele)
summary(mod_l_ele)

R2m=r.squaredGLMM(mod_l_ele)[[1]]
R2c=r.squaredGLMM(mod_l_ele)[[2]]
R2m
R2c

mod_l_ele_lm <- lm(ShiftKmY ~ AnnualDispPotKmY,
                   data = limited_ele)
summary(mod_l_ele_lm)

R2m=r.squaredGLMM(mod_l_ele_lm)[[1]]
R2m

AIC(mod_l_ele_lm, ## this is better fit
    mod_l_ele)

# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## plot predictions  
new_data <- expand.grid(AnnualDispPotKmY = seq(min(notlimited_lat$AnnualDispPotKmY),
                                               max(notlimited_lat$AnnualDispPotKmY), 
                                               by = 0.01))

pred_nl_lat <- predict(mod_nl_lat, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_nl_lat <- new_data %>%
  mutate(pred_expansion = pred_nl_lat$fit) %>%
  mutate(pred_expansion_SE = pred_nl_lat$se.fit)


nl_lat_nopred <- v3 %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY_spp) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

nl_lat_pred <- nl_lat_nopred +
  geom_ribbon(data = fitted_pred_nl_lat, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_nl_lat, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(nl_lat_nopred, path = "figures/sotm", filename = "data-notlimited-lat-nopredictions_spp.png", 
       device = "png", height = 4, width = 8)
ggsave(nl_lat_pred, path = "figures/sotm", filename = "data-notlimited-lat-predictions_spp.png", 
       device = "png", height = 4, width = 8)


## now ele
new_data <- expand.grid(AnnualDispPotKmY = seq(min(notlimited_ele$AnnualDispPotKmY),
                                               max(notlimited_ele$AnnualDispPotKmY), 
                                               by = 0.01))

pred_nl_ele <- predict(mod_nl_ele, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_nl_ele <- new_data %>%
  mutate(pred_expansion = pred_nl_ele$fit) %>%
  mutate(pred_expansion_SE = pred_nl_ele$se.fit)

nl_ele_nopred <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY_spp) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none")  +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

nl_ele_pred <- nl_ele_nopred +
  geom_ribbon(data = fitted_pred_nl_ele, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_nl_ele, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(nl_ele_nopred, path = "figures/sotm", filename = "data-notlimited-ele-nopredictions_spp.png", 
       device = "png", height = 2.8, width = 3.1)
ggsave(nl_ele_pred, path = "figures/sotm", filename = "data-notlimited-ele-predictions_spp.png", 
       device = "png", height = 2.8, width = 3.1)


new_data <- expand.grid(AnnualDispPotKmY = seq(min(limited_lat$AnnualDispPotKmY),
                                               max(limited_lat$AnnualDispPotKmY), 
                                               by = 0.0001))

pred_l_lat <- predict(mod_l_lat, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_l_lat <- new_data %>%
  mutate(pred_expansion = pred_l_lat$fit) %>%
  mutate(pred_expansion_SE = pred_l_lat$se.fit)


l_lat_nopred <- v3 %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY_spp) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

l_lat_pred <- l_lat_nopred +
  geom_ribbon(data = fitted_pred_l_lat, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_l_lat, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(l_lat_nopred, path = "figures/sotm", filename = "data-limited-lat-nopredictions_spp.png", 
       device = "png", height = 4, width = 8)
ggsave(l_lat_pred, path = "figures/sotm", filename = "data-limited-lat-predictions_spp.png", 
       device = "png", height = 4, width = 8)


## now ele
new_data <- expand.grid(AnnualDispPotKmY = seq(min(limited_ele$AnnualDispPotKmY),
                                               max(limited_ele$AnnualDispPotKmY), 
                                               by = 0.0001))

pred_l_ele <- predict(mod_l_ele, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_l_ele <- new_data %>%
  mutate(pred_expansion = pred_l_ele$fit) %>%
  mutate(pred_expansion_SE = pred_l_ele$se.fit)

l_ele_nopred <- v3 %>%
  filter(Rate > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY_spp) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 1.5) 

l_ele_pred <- l_ele_nopred +
  geom_ribbon(data = fitted_pred_l_ele, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_l_ele, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(l_ele_nopred, path = "figures/sotm", filename = "data-limited-ele-nopredictions_spp.png", 
       device = "png", height = 2.8, width = 3.1)
ggsave(l_ele_pred, path = "figures/sotm", filename = "data-limited-ele-predictions_spp.png", 
       device = "png", height = 2.8, width = 3.1)





