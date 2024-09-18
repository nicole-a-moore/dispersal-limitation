## mock bayesian model for committee meeting 
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)

theme_set(theme_minimal())

########################
## read in model
########################
## in this model, range expansion rate = the minimum of clim velo and dispersal + sampling error
## read in stan model
## this model has gamma error distributions + two slopes 
## before breakpoint = 1, after breakpoint = 0
breakpoint_minimum <- cmdstan_model("stan/breakpoint_minimum.stan")
breakpoint_minimum

########################
## feed the model data! 
########################
## read in data 
dd <- read.csv("data-processed/v1_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Gradient == "Latitudinal")

## choose a single study
## pick the one with the most observations 
dd %>% 
  group_by(Source) %>%
  tally() %>%
  arrange(-n)
## study A32_P1

study <- dd %>%
  filter(Source == "A32_P1")

## get vector of observed shifts 
shift_obs <- dd$ShiftKmY

## feed the data to the model
breakpoint_minimum_posterior <- breakpoint_minimum$sample(parallel_chains = 4, 
                                                          data = list(disp = dd$DispersalPotentialKmY,
                                                                      n = nrow(dd), 
                                                                      climVelo = dd$LatVeloT,
                                                                      shift = dd$ShiftKmY)) 

breakpoint_minimum_posterior
## sigma = 1.83


shinystan::launch_shinystan(breakpoint_minimum_posterior)


###########################
## make pretty plot
###########################
prep <- dd %>%
  select(AnnualDispPotKmY, ShiftKmY, ClimVeloTKmY) %>%
  mutate(type = "Observed sample") 

## get 11 samples 
pretty_data <- breakpoint_minimum_posterior %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 1)%>% 
  mutate(AnnualDispPotKmY = dd$DispersalPotentialKmY[i],
         ClimVeloTKmY = dd$ClimVeloTKmY[i]) %>%
  rename("ShiftKmY" = .value) %>%
  mutate(type = "Simulated sample", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(AnnualDispPotKmY, ClimVeloTKmY, type, ShiftKmY) %>%
  rbind(., prep)

pretty_data %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
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
       y = "Range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  facet_wrap(~type)

ggsave(filename = "minimum-model_obs-vs-simulated.png", path = "figures/proposal",
       width = 7, height = 3.25)


dd %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
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
clim_velos <- c(dd$ClimVeloTKmY)

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

gradients <- select(ungroup(dd), ClimVeloTKmY, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY"))

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



dd = dd %>%
  mutate(limiting_rate = ifelse(DispersalPotentialKmY < ClimVeloTKmY, DispersalPotentialKmY, ClimVeloTKmY)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == limiting_rate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY, NA)) 

dd %>%
  ggplot(aes(x = limiting_rate, y = ShiftKmY, colour = ClimVeloTKmY, shape = group)) +
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
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY, shape = group)) +
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






###########################
## making my own PPC plots
###########################
## plot shift_obs and 11 samples of shift_pred
breakpoint_minimum_posterior

## get 11 samples 
samples <- breakpoint_minimum_posterior %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 11) %>% 
  mutate(disp = dd$DispersalPotentialKmY[i]) %>%
  rename("shift" = .value, "DispersalPotentialKmY" = disp) %>%
  mutate(type = "Simulated sample", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(DispersalPotentialKmY, type, shift, .draw)

## add real data to the dataframe 
samples <- dd %>%
  rename("shift" = ShiftKmY) %>%
  mutate(type = "Observed sample", .draw = "shift_obs") %>%
  select(shift, type, DispersalPotentialKmY, .draw) %>%
  rbind(., samples)

## plot scatter plot
samples %>% 
  ggplot(aes(x = DispersalPotentialKmY, y = shift, colour = type)) + 
  geom_point() + 
  facet_wrap(~.draw) +
  labs(colour = "",
       x = "Potential dispersal rate (km/y)", 
       y = "Range shift rate (km/y)") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_color_manual(values = c("Red", "Black"))


## plot histograms & density plots 
samples %>% 
  ggplot(aes(x = shift, fill = type)) + 
  geom_histogram() + 
  facet_wrap(~.draw)

samples %>% 
  ggplot(aes(x = shift, fill = type)) + 
  geom_density() + 
  facet_wrap(~.draw)

## overlay 
samples %>% 
  mutate(type = factor(.$type, ordered = TRUE, levels = c("Observed sample", "Simulated sample"))) %>%
  ggplot(aes(x = shift, 
             fill = type, alpha = type, 
             group = .draw, 
             colour = type)) + 
  scale_fill_manual(values = c("red", "transparent")) +
  scale_colour_manual(values = c("red", "black")) +
  scale_alpha_manual(values = c(0.4, 0)) +
  geom_density() +
  labs(colour = "", alpha = "", fill = "", 
       x = "Potential dispersal rate (km/y)", 
       y = "Density") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

ggsave(filename = "minimum-model_density.png", path = 'figures/proposal',
       width = 6, height = 4)

all_samples = breakpoint_minimum_posterior %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 4000) %>% 
  mutate(disp = dd$DispersalPotentialKmY[i]) %>%
  rename("shift" = .value, "DispersalPotentialKmY" = disp) %>%
  mutate(type = "Simulated sample", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(DispersalPotentialKmY, type, shift, .draw)

samples_summary <- all_samples %>%
  group_by(.draw) %>%
  mutate(mean_shift = mean(shift), max_shift = max(shift), min_shift = min(shift), sd_shift = sd(shift)) %>%
  ungroup() %>%
  select(-.draw, -shift, -DispersalPotentialKmY) %>%
  distinct()

## add real data to the dataframe 
real_mean <- mean(dd$ShiftKmY)
real_max <- max(dd$ShiftKmY)
real_min <- min(dd$ShiftKmY)
real_sd <- sd(dd$ShiftKmY)

## MEAN
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = mean_shift)) + 
  geom_histogram(colour = "black", fill = "black", alpha = 0.05) +
  geom_vline(xintercept = real_mean, colour = "red") +
  labs(x = "Mean range shift rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1400))

ggsave(filename = "minimum-model_hist-mean-shift.png", path = 'figures/proposal',
       width = 2.5, height = 2)

samples_summary %>%
  ggplot(aes(x = mean_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_mean) +
  theme(legend.position = "none") 


## MAX
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = max_shift)) + 
  geom_histogram(colour = "black", fill = "black", alpha = 0.05) +
  geom_vline(xintercept = real_max, colour = "red") +
  labs(x = "Max range shift rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3500))

ggsave(filename = "minimum-model_hist-max-shift.png", path = 'figures/proposal',
       width = 2.5, height = 2)

samples %>%
  ggplot(aes(x = max_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_max) +
  theme(legend.position = "none") 

## MIN
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = min_shift)) + 
  geom_histogram(colour = "black", fill = "black", alpha = 0.05) +
  geom_vline(xintercept = real_min, colour = "red") +
  labs(x = "Min range shift rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4100))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1,
                                           decimal.mark = ','), 
                     limits = c(-0.05, 1))

ggsave(filename = "minimum-model_hist-min-shift.png", path = 'figures/proposal',
       width = 2.5, height = 2)

samples_summary %>%
  ggplot(aes(x = min_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_min) +
  theme(legend.position = "none") 

## SD
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = sd_shift)) + 
  geom_histogram(colour = "black", fill = "black", alpha = 0.05) +
  geom_vline(xintercept = real_sd, colour = "red") +
  labs(x = "SD range shift rate (km/y)", y = "Frequency") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4000))

ggsave(filename = "minimum-model_hist-sd-shift.png", path = 'figures/proposal',
       width = 2.5, height = 2)

samples_summary %>%
  ggplot(aes(x = sd_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_sd) +
  theme(legend.position = "none") 


## make plot of average shift_pred versus shift_obs 
avg <- all_samples %>%
  group_by(DispersalPotentialKmY) %>%
  mutate(avg_shift_pred = mean(shift)) %>%
  ungroup() %>%
  select(-.draw, -shift, -type) %>%
  distinct()

dd %>%
  select(DispersalPotentialKmY, ShiftKmY) %>%
  rename("shift_obs" = ShiftKmY) %>%
  left_join(., avg) %>% 
  ggplot(aes(x = shift_obs, y = avg_shift_pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed range shift rate (km/y)", y = "Mean simulated\nrange shift rate (km/y)") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) 
  
ggsave(filename = "minimum-model_scatter-obs-pred.png", path = 'figures/proposal',
       width = 2.8, height = 2)






