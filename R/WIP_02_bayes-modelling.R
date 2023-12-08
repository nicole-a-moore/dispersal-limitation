## writing my own models 
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)

theme_set(theme_minimal())


########################
## read in my model
########################
## read in stan model
## this model has gamma error distributions
breakpoint_B_gamma_Nikki <- cmdstan_model("stan/breakpoint_B_gamma_Nikki.stan")
breakpoint_B_gamma_Nikki


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
shift_obs <- study$ShiftKmY

## feed the data to the model
breakpoint_B_gamma_Nikki_posterior <- breakpoint_B_gamma_Nikki$sample(parallel_chains = 4, 
                                                                      data = list(disp = study$DispersalPotentialKmY,
                                                                                  n = nrow(study), 
                                                                                  climVelo = unique(study$LatVeloT),
                                                                                  shift = study$ShiftKmY)) 

breakpoint_B_gamma_Nikki_posterior
## b2 = 3.9
## sigma = 0.17

shinystan::launch_shinystan(breakpoint_B_gamma_Nikki_posterior)

## rhat = like an F ratio, want it to be close to 1 

## posterior predictive check: make sure shift_obs looks like shift_pred


## PPC test statistic
## interpretation of tests based on the posterior predictive distribution is straightforward. It represents the probability that a future observation will exceed the existing data, given the model. An extreme p value, therefore, implies poor model fit
## mean fake shift is less than mean real shift 
## sd fake shift is less than sd real shift 
## min fake shift is less than min real shift 
## max fake shift is less than max real shift 
## average y rep > observed y


## try a different study 
study2 <- dd %>%
  filter(Source == "A10_P1")

## get vector of observed shifts 
shift_obs2 <- study2$ShiftKmY

## feed the data to the model
breakpoint_B_gamma_Nikki_posterior2 <- breakpoint_B_gamma_Nikki$sample(parallel_chains = 4, 
                                                                      data = list(disp = study2$DispersalPotentialKmY,
                                                                                  n = nrow(study2), 
                                                                                  climVelo = unique(study2$LatVeloT),
                                                                                  shift = study2$ShiftKmY)) 

breakpoint_B_gamma_Nikki_posterior2
## b2 = 4.01
## sigma = 4.46

shinystan::launch_shinystan(breakpoint_B_gamma_Nikki_posterior2)

## better fit, but still not great
hist(shift_obs2)




###########################
## making my own PPC plots
###########################

## plot shift_obs and 11 samples of shift_pred
breakpoint_B_gamma_Nikki_posterior2

## get 11 samples 
samples <- breakpoint_B_gamma_Nikki_posterior2 %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 11) %>% 
  mutate(disp = study2$DispersalPotentialKmY[i]) %>%
  rename("shift" = .value, "DispersalPotentialKmY" = disp) %>%
  mutate(type = "predicted", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(DispersalPotentialKmY, type, shift, .draw)

## add real data to the dataframe 
samples <- study2 %>%
  rename("shift" = ShiftKmY) %>%
  mutate(type = "observed", .draw = "shift_obs") %>%
  select(shift, type, DispersalPotentialKmY, .draw) %>%
  rbind(., samples)

## plot scatter plot
samples %>% 
  ggplot(aes(x = DispersalPotentialKmY, y = shift, colour = type)) + 
  geom_point() + 
  facet_wrap(~.draw)

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
  mutate(type = factor(.$type, ordered = TRUE, levels = c("observed", "predicted"))) %>%
  ggplot(aes(x = shift, 
             fill = type, alpha = type, 
             group = .draw)) + 
  scale_fill_manual(values = c("red", "transparent")) +
  scale_alpha_manual(values = c(0.4, 0)) +
  geom_density() 


all_samples = breakpoint_B_gamma_Nikki_posterior2 %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 4000) %>% 
  mutate(disp = study2$DispersalPotentialKmY[i]) %>%
  rename("shift" = .value, "DispersalPotentialKmY" = disp) %>%
  mutate(type = "predicted", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(DispersalPotentialKmY, type, shift, .draw)

samples_summary <- all_samples %>%
  group_by(.draw) %>%
  mutate(mean_shift = mean(shift), max_shift = max(shift), min_shift = min(shift), sd_shift = sd(shift)) %>%
  ungroup() %>%
  select(-.draw, -shift, -DispersalPotentialKmY) %>%
  distinct()

## add real data to the dataframe 
real_mean <- mean(study2$ShiftKmY)
real_max <- max(study2$ShiftKmY)
real_min <- min(study2$ShiftKmY)
real_sd <- sd(study2$ShiftKmY)

## MEAN
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = mean_shift, fill = type)) + 
  geom_histogram() +
  geom_vline(xintercept = real_mean) +
  theme(legend.position = "none")

samples_summary %>%
  ggplot(aes(x = mean_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_mean) +
  theme(legend.position = "none") 


## MAX
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = max_shift, fill = type)) + 
  geom_histogram() +
  geom_vline(xintercept = real_max) +
  theme(legend.position = "none")

samples %>%
  ggplot(aes(x = max_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_max) +
  theme(legend.position = "none") 

## MIN
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = min_shift, fill = type)) + 
  geom_histogram() +
  geom_vline(xintercept = real_min) +
  theme(legend.position = "none")

samples_summary %>%
  ggplot(aes(x = min_shift, fill = type, colour = type)) + 
  geom_density() +
  geom_vline(xintercept = real_min) +
  theme(legend.position = "none") 

## SD
## plot mean of predicted versus mean of observed 
samples_summary %>%
  ggplot(aes(x = sd_shift, fill = type)) + 
  geom_histogram() +
  geom_vline(xintercept = real_sd) +
  theme(legend.position = "none")

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

study2 %>%
  select(DispersalPotentialKmY, ShiftKmY) %>%
  rename("shift_obs" = ShiftKmY) %>%
  left_join(., avg) %>% 
  ggplot(aes(x = shift_obs, y = avg_shift_pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)



## next steps
## regression line drawing :-)
## add more parameters 
## write new one that can work on multiple studies 



#####################################################
## write model that accepts data from across studies 
#####################################################
## read in stan model
## this model has gamma error distributions + two slopes 
## before breakpoint = 1, after breakpoint = 0
breakpoint_all_studies <- cmdstan_model("stan/breakpoint_all_studies.stan")
breakpoint_all_studies

shift_obs = dd$ShiftKmY

## feed the data to the model
breakpoint_all_studies_posterior <- breakpoint_all_studies$sample(parallel_chains = 4, 
                                                            data = list(disp = dd$DispersalPotentialKmY,
                                                                        n = nrow(dd), 
                                                                        climVelo = dd$LatVeloT,
                                                                        shift = dd$ShiftKmY)) 

breakpoint_all_studies_posterior
## b2 = 1.04 
## sigma = 1.89

shinystan::launch_shinystan(breakpoint_all_studies_posterior)

## simulated data still doesn't look perfectly like observed 
## mostly, observed shifts are still often faster than predicted shifts 
## I could've guessed this would happen
## I wonder - what if we try and control for methods / use corrected shift?

## get 11 samples 
samples <- breakpoint_all_studies_posterior %>%
  tidybayes::gather_draws(shift_pred[i], ndraws = 11) %>% 
  mutate(disp = dd$DispersalPotentialKmY[i]) %>%
  rename("shift" = .value, "DispersalPotentialKmY" = disp) %>%
  mutate(type = "predicted", .draw = paste0("shift_pred_", .draw)) %>%
  ungroup() %>%
  select(DispersalPotentialKmY, type, shift, .draw)

## add real data to the dataframe 
samples <- dd %>%
  rename("shift" = ShiftKmY) %>%
  mutate(type = "observed", .draw = "shift_obs") %>%
  select(shift, type, DispersalPotentialKmY, .draw) %>%
  rbind(., samples)

## plot scatter plot
samples %>% 
  ggplot(aes(x = DispersalPotentialKmY, y = shift, colour = type)) + 
  geom_point() + 
  facet_wrap(~.draw)



