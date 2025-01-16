## mock bayesian model for committee meeting 
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)

theme_set(theme_minimal())


########################
## read in the data
########################
## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## filter 
dd <- dd %>%
  ## filter to cases where we expect dispersal to matter (range expansion): 
  ## leading edge shifts where climate velocity is positive
  ## trailing edge shifts where climate velocity is negative 
  ## centroid shifts
  filter(Param == "TE" & mean_cv_studylevel <= 0 |
           Param == "LE" & mean_cv_studylevel >= 0 |
           Param == "O") 

## get rid of non-expansions 
dd = filter(dd, is_contraction %in% c("UNKNOWN", "NO"))

## get rid of centroid shifts in opposite direction to cv 
dd = filter(dd, !(Param == "O" & !sign(ClimVeloTKmY_spp) == sign(ShiftKmY)))

## standardize so CV > 0 meanings away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ShiftKmY = abs(dd$ShiftKmY)

## add gradient 
dd <- dd %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

hist(dd$ClimVeloTKmY_study)
hist(dd$ClimVeloTKmY_spp)

dd = dd %>%
  mutate(limiting_rate = ifelse(DispersalPotentialKmY < ClimVeloTKmY_study, DispersalPotentialKmY, ClimVeloTKmY_study)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == limiting_rate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_study, NA))

###########################
## make pretty plots
###########################
dd %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_study)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
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
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

ggsave(filename = "minimum-model_obs.png", path = "figures/proposal",
       width = 4.1, height = 3.25)

## predictions plot:
clim_velos <- c(dd$ClimVeloTKmY_study)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1400, by = 1)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 40, by = 0.01)) 
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

gradients <- select(ungroup(dd), ClimVeloTKmY_study, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY_study"))

pred_plot <- clim_data %>%
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
  scale_x_log10(limits = c(0.001, 1400), 
                breaks = c(0.001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("0", "10", "20", "30", "40")) +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range shift rate (km/y)", 
       colour = "") 

ggsave(pred_plot, filename = "minimum-model_pred.png", path = "figures/proposal",
       width = 3.5, height = 3.25)


dd %>%
  ggplot(aes(x = limiting_rate, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Minimum of rate of climate change and\npotential dispersal rate (km/y)",
       y = "Range shift rate (km/y)", 
       colour = 'Limiting\nfactor', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") 

ggsave(filename = "minimum-model_limiting-rate.png", path = "figures/proposal",
       width = 4.25, height = 3.5)

dd %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_point(data = filter(dd, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Range shift rate (km/y)", 
       colour = 'Limiting\nfactor', 
       shape = "") +
  scale_y_continuous(limits = c(0,40)) 

ggsave(filename = "minimum-model_dispersal-rate.png", path = "figures/proposal",
       width = 4.35, height = 3.5)


## some notes: 
## largest shift in wild turkey - habitat restoration and reintroduction since 1980
## second is bald eagle 

