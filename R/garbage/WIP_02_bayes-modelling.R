## writing my own models 
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(loo)

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

## filter to only leading edge observations 
#dd = filter(dd, Position == 'Leading edge')

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
  facet_wrap(~.draw) +
  scale_x_log10()




#####################################################
## write model without slope parameter 
#####################################################
## in this model, range expansion rate = the minimum of clim velo and dispersal + sampling error
## read in stan model
## this model has gamma error distributions + two slopes 
## before breakpoint = 1, after breakpoint = 0
breakpoint_minimum <- cmdstan_model("stan/breakpoint_minimum.stan")
breakpoint_minimum


shift_obs = dd$ShiftKmY

## feed the data to the model
breakpoint_minimum_posterior <- breakpoint_minimum$sample(parallel_chains = 4, 
                                                          data = list(disp = dd$DispersalPotentialKmY,
                                                                      n = nrow(dd), 
                                                                      climVelo = dd$LatVeloT,
                                                                      shift = dd$ShiftKmY)) 

breakpoint_minimum_posterior
## sigma = 1.83


shinystan::launch_shinystan(breakpoint_minimum_posterior)

## get 11 samples 
samples <- breakpoint_minimum_posterior %>%
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
  facet_wrap(~.draw) +
  scale_x_log10() 

samples %>% 
  ggplot(aes(x = DispersalPotentialKmY, y = shift, colour = type)) + 
  geom_point() + 
  facet_wrap(~.draw) +
  scale_x_log10() + 
  scale_y_log10() 
## model allows for really small shifts because climate velocity + dispersal potential is allowed to be really small

sort(samples$shift)[which(sort(samples$shift) != 0)]
sort(dd$ShiftKmY)[which(sort(dd$ShiftKmY) != 0)]
## but we know that studies can't detect shifts of say 6.284515e-321 km/y
## could we somehow implement a minimum detectable shift? 
## like shifts that are less a certain amount are either 0 shifts or are the minimum detected shift in the data?


min(dd$ShiftKmY) ## minimum detected shift = 2.18e-17
hist(dd$ShiftKmY)

exp(1)^quantile(log(dd$ShiftKmY), 0.05)

samples %>% 
  filter(shift >= 0.02579858,
         DispersalPotentialKmY >= 0.02579858) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = shift, colour = type)) + 
  geom_point() + 
  facet_wrap(~.draw) +
  scale_x_log10() + 
  scale_y_log10() 






## falsifying a model:

### could compare to a model with just an intercept 
### could fit and see if it fails 
### could start with simplest model and add complexity, then compete models 






## read in stan model
## this model has gamma error distributions + two slopes 
## before breakpoint = 1, after breakpoint = 0
breakpoint_all_studies <- cmdstan_model("stan/breakpoint_all_studies_noneg.stan", pedantic = TRUE)
breakpoint_all_studies

## feed the data to the model
breakpoint_all_studies_posterior <- breakpoint_all_studies$sample(parallel_chains = 4, 
                                                                  data = list(disp = dd$DispersalPotentialKmY,
                                                                              n = nrow(dd), 
                                                                              climVelo = dd$LatVeloT,
                                                                              shift = dd$ShiftKmY)) 
## throws errors about non-finite parameters, but only during warmup - so it's okay 

breakpoint_all_studies_posterior
## b2 = 1.04 
## sigma = 1.89


### plotting 
## choose one study
## choose one with the highest climate veloicity 

study <- dd %>%
  filter(ClimVeloTKmY == max(ClimVeloTKmY))

study %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point() 


## generate quantities 
n = 100
disp_vec <- modelr::seq_range(study$DispersalPotentialKmY, n)
disp_vec <- seq(0.00001, max(study$DispersalPotentialKmY), length.out = n)

breakpoint_generate_quantities <- cmdstan_model("stan/breakpoint_generate_quantities.stan", pedantic = TRUE)

breakpoint_generate_quantities <- breakpoint_generate_quantities$generate_quantities(fitted_params = breakpoint_all_studies_posterior,
                                                                                     parallel_chains = 4, 
                                                                                     data = list(disp = disp_vec,
                                                                                                 n = n, 
                                                                                                 climVelo = rep(unique(study$ClimVeloTKmY), n))) 

breakpoint_generate_quantities

library(tidybayes)

# breakpoint_generate_quantities %>%
#   gather_rvars(shift_pred[i]) %>% head(2) %>% 
#   #mutate(disp = disp_vec[i]) %>%
#   ggplot(aes(y = i, dist = .value)) +
#   stat_histinterval()

## this is like broom
## .variable = model parameter 
## .value is value of parameter 

breakpoint_generate_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = DispersalPotentialKmY, y = ShiftKmY), inherit.aes = FALSE) 

## zoom in
breakpoint_generate_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = DispersalPotentialKmY, y = ShiftKmY), inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0,10), ylim = c(0, 10))


## plot all the predictions:
n = 100
cvs = unique(dd$ClimVeloTKmY)
disp_vec <- c()
for(i in cvs) {
  study = filter(dd, ClimVeloTKmY == i)
  disp_vec = append(disp_vec, seq(0.00001, max(study$DispersalPotentialKmY), length.out = n))
}

breakpoint_generate_quantities <- cmdstan_model("stan/breakpoint_generate_quantities.stan", pedantic = TRUE)

breakpoint_generate_quantities <- breakpoint_generate_quantities$generate_quantities(fitted_params = breakpoint_all_studies_posterior,
                                                                                     parallel_chains = 4, 
                                                                                     data = list(disp = disp_vec,
                                                                                                 n = n*length(unique(dd$ClimVeloTKmY)), 
                                                                                                 climVelo = rep(unique(dd$ClimVeloTKmY), each = n))) 

breakpoint_generate_quantities

## plot one 
study_1 = filter(dd, ClimVeloTKmY == first(dd$ClimVeloTKmY))
breakpoint_generate_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i],
         climVelo = rep(unique(dd$ClimVeloTKmY), each = n)[i]) %>%
  filter(climVelo == first(dd$ClimVeloTKmY)) %>% 
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study_1, aes(x = DispersalPotentialKmY, y = ShiftKmY), inherit.aes = FALSE) 


## plot all
preds <- breakpoint_generate_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i],
         climVelo = rep(unique(dd$ClimVeloTKmY), each = n)[i]) 

preds %>%
  ggplot(aes(x = disp, dist = .value, group = climVelo)) +
  stat_lineribbon() +
  geom_point(data = dd, aes(x = DispersalPotentialKmY, y = ShiftKmY), inherit.aes = FALSE) 

## zoom in
preds %>%
  ggplot(aes(x = disp, dist = .value, group = climVelo)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = DispersalPotentialKmY, y = ShiftKmY), inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0,10), ylim = c(0, 10))

## facet by study 
dd$climVelo = dd$ClimVeloTKmY

preds %>%
  ggplot(aes(x = disp, dist = .value, group = climVelo)) +
  stat_lineribbon() +
  geom_point(data = dd, aes(x = DispersalPotentialKmY, y = ShiftKmY, dist = NA)) + 
  facet_wrap(~climVelo)

## zoom in 
preds %>%
  ggplot(aes(x = disp, dist = .value, group = climVelo)) +
  stat_lineribbon() +
  geom_point(data = dd, aes(x = DispersalPotentialKmY, y = ShiftKmY, dist = NA)) + 
  facet_wrap(~climVelo) +
  coord_cartesian(xlim = c(0,10), ylim = c(0, 10))
## data do not really follow mean - hmm 

## other thoughts:
## low signal to noise ratio when climate velocity is low - could lead to bias in shift detection 




#### null model
#############################
null_model <- cmdstan_model("stan/null-model.stan", pedantic = TRUE)
null_model

## feed the data to the model
null_posterior <- null_model$sample(parallel_chains = 4, 
                                                  data = list(disp = dd$DispersalPotentialKmY,
                                                                              n = nrow(dd), 
                                                                              climVelo = dd$LatVeloT,
                                                                              shift = dd$ShiftKmY)) 

shinystan::launch_shinystan(null_posterior)

null_loo <- null_posterior$loo()
null_loo$pointwise
plot(null_loo)

#### model with b2 
#############################
breakpoint_all_studies <- cmdstan_model("stan/breakpoint_all_studies.stan")
breakpoint_all_studies


breakpoint_all_studies_posterior <- breakpoint_all_studies$sample(parallel_chains = 4, 
                                                                  data = list(disp = dd$DispersalPotentialKmY,
                                                                              n = nrow(dd), 
                                                                              climVelo = dd$LatVeloT,
                                                                              shift = dd$ShiftKmY)) 
                                                                                   
bp_loo <- breakpoint_all_studies_posterior$loo()
bp_loo$pointwise
plot(bp_loo)

# compare the models
comparison <- loo_compare(null_loo, bp_loo, intercept_only_loo, intercept_slope_loo)

null_posterior
breakpoint_all_studies_posterior


#### null model with random effect of study ID
###############################################
## add random effects
## need to make sure when you add/subtract, value doesn't go below zero
## one way to do this: use a log scale 
intercept_only_model <- cmdstan_model("stan/intercept-only-model.stan", pedantic = TRUE)
intercept_only_model

## feed the data to the model
intercept_only_model_posterior <- intercept_only_model$sample(parallel_chains = 4, 
                                    data = list(n = nrow(dd), 
                                                num_study = length(unique(dd$Source)),
                                                disp = dd$DispersalPotentialKmY,
                                                climVelo = dd$LatVeloT,
                                                shift = dd$ShiftKmY,
                                                study_id = as.numeric(factor(dd$Source)))) 

shinystan::launch_shinystan(intercept_only_model_posterior)

intercept_only_loo <- intercept_only_model_posterior$loo()
intercept_only_loo$pointwise
plot(intercept_only_loo)

# compare the models
comparison <- loo_compare(null_loo, intercept_only_loo)

#### slope model with random effect of study ID
###############################################
intercept_slope_model <- cmdstan_model("stan/intercept-slope-model.stan", pedantic = TRUE)
intercept_slope_model

## feed the data to the model
intercept_slope_model_posterior <- intercept_slope_model$sample(parallel_chains = 4, 
                                                                data = list(n = nrow(dd), 
                                                                            num_study = length(unique(dd$Source)),
                                                                            disp = dd$DispersalPotentialKmY,
                                                                            climVelo = dd$LatVeloT,
                                                                            shift = dd$ShiftKmY,
                                                                            study_id = as.numeric(factor(dd$Sourc )))) 

shinystan::launch_shinystan(intercept_slope_model_posterior)

intercept_slope_loo <- intercept_slope_model_posterior$loo()
intercept_slope_loo$pointwise
plot(intercept_slope_loo)

# compare the models
comparison <- loo_compare(null_loo, intercept_only_loo, intercept_slope_loo)

intercept_only_model_posterior$summary("sigma_log_study_diffs")
intercept_slope_model_posterior$summary("sigma_log_study_diffs")


#### null model with random effect of study ID + species 
###############################################
## add random effects
## need to make sure when you add/subtract, value doesn't go below zero
## one way to do this: use a log scale 
intercept_only_model_sp_reff <- cmdstan_model("stan/intercept-only-model-species-reff.stan", pedantic = TRUE)
intercept_only_model_sp_reff

## feed the data to the model
intercept_only_model_sp_reff_posterior <- intercept_only_model_sp_reff$sample(parallel_chains = 4, 
                                                              data = list(n = nrow(dd), 
                                                                          num_study = length(unique(dd$Source)),
                                                                          num_sp = length(unique(dd$scientificName)),
                                                                          disp = dd$DispersalPotentialKmY,
                                                                          climVelo = dd$LatVeloT,
                                                                          shift = dd$ShiftKmY,
                                                                          study_id = as.numeric(factor(dd$Source)),
                                                                          sp_id = as.numeric(factor(dd$scientificName)))) 

shinystan::launch_shinystan(intercept_only_model_sp_reff_posterior)

intercept_only_model_sp_reff_posterior$summary(c("sigma_log_study_diffs", "sigma_log_sp_diffs"))

#### slope model with random effect of study ID + species 
###############################################
## add random effects
## need to make sure when you add/subtract, value doesn't go below zero
## one way to do this: use a log scale 
intercept_slope_model_sp_reff <- cmdstan_model("stan/intercept-slope-model-species-reff.stan", pedantic = TRUE)
intercept_slope_model_sp_reff

## feed the data to the model
intercept_slope_model_sp_reff_posterior <- intercept_slope_model_sp_reff$sample(parallel_chains = 4, 
                                                                              data = list(n = nrow(dd), 
                                                                                          num_study = length(unique(dd$Source)),
                                                                                          num_sp = length(unique(dd$scientificName)),
                                                                                          disp = dd$DispersalPotentialKmY,
                                                                                          climVelo = dd$LatVeloT,
                                                                                          shift = dd$ShiftKmY,
                                                                                          study_id = as.numeric(factor(dd$Source)),
                                                                                          sp_id = as.numeric(factor(dd$scientificName)))) 

shinystan::launch_shinystan(intercept_slope_model_sp_reff_posterior)

intercept_slope_model_sp_reff_posterior$summary(c("sigma_log_study_diffs", "sigma_log_sp_diffs"))


#### compare models
###############################################



