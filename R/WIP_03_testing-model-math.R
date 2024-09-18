## investigating why the model is making weird prediction intervals 
## generate simulated data that follows the expected relationship 
## then see what prediction interval looks like 
library(tidybayes)

## simulate data 
########################################
## randomly generate 100 climate velocities 
climVelo = runif(100, 0.01, 7)
climVelo = rep(climVelo, each = 50)

## randomly generate 50 dispersal rates between 0-15 for each clim velo
dispRs = runif(length(climVelo), 0.01, 15)


## generate shifts 
########################################
## figure out which is the min of clim velo and dispersal rate 
min_rate = ifelse(climVelo < dispRs, climVelo, dispRs)

## set sigma and b2
sigma = 1.89
b2 = 1.04

## IMPRORTANT:
###################################
## increasing sigma causes the number of 0 shifts to increase - why? could this be the problem?

## it also prevents the model from fitting:
##  Error evaluating the log probability at the initial value.
##  Random variable is 0, but must be positive finite!
##  Inverse scale parameter is -96.8433, but must be positive finite!

## seems like when sd is high, model passes zero values to the gamma function (which is a problem)
## get around this by adding a constant 0.1 to shifts 

## calculate shifts 
shift = rgamma(length(min_rate), (min_rate*b2)^2 / sigma^2, (min_rate*b2) / sigma^2)
shift = shift + 0.01
hist(shift)

## fit model to shifts 
########################################
fake_data <- data.frame(dispRs, climVelo, shift)

## plot it 
fake_data %>%
  ggplot(aes(x = dispRs, y = shift, colour = climVelo)) +
  geom_point()

## plot one
fake_data %>%
  filter(climVelo == max(climVelo)) %>%
  ggplot(aes(x = dispRs, y = shift, colour = climVelo)) +
  geom_point()

breakpoint_all_studies <- cmdstan_model("stan/breakpoint_all_studies_noneg.stan")
breakpoint_all_studies

## feed the data to the model
breakpoint_all_studies_posterior <- breakpoint_all_studies$sample(parallel_chains = 4, 
                                                                  data = list(disp = fake_data$dispRs,
                                                                              n = nrow(fake_data), 
                                                                              climVelo = fake_data$climVelo,
                                                                              shift = fake_data$shift)) 
breakpoint_all_studies_posterior

#shinystan::launch_shinystan(breakpoint_all_studies_posterior)
## looks good 

## the moment of truth: plot the predictions
#############################################
## get data from one "study" - one with fastest climate velocity 
study = filter(fake_data, climVelo == max(climVelo))
n = 50
disp_vec <- modelr::seq_range(study$dispRs, n)

breakpoint_generate_quantities <- cmdstan_model("stan/breakpoint_generate_quantities.stan", pedantic = TRUE)

breakpoint_generate_quantities_quantities <- breakpoint_generate_quantities$generate_quantities(fitted_params = breakpoint_all_studies_posterior,
                                                                                     parallel_chains = 4, 
                                                                                     data = list(disp = disp_vec,
                                                                                                 n = n, 
                                                                                                 climVelo = rep(unique(study$climVelo), n))) 

breakpoint_generate_quantities_quantities

breakpoint_generate_quantities_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = dispRs, y = shift), inherit.aes = FALSE) 

## now from one with slowest climate velocity 
study = filter(fake_data, climVelo == min(climVelo))
n = 50
disp_vec <- modelr::seq_range(study$dispRs, n)

breakpoint_generate_quantities_quantities <- breakpoint_generate_quantities$generate_quantities(fitted_params = breakpoint_all_studies_posterior,
                                                                                     parallel_chains = 4, 
                                                                                     data = list(disp = disp_vec,
                                                                                                 n = n, 
                                                                                                 climVelo = rep(unique(study$climVelo), n))) 

breakpoint_generate_quantities_quantities %>% 
  gather_rvars(shift_pred[i]) 

breakpoint_generate_quantities_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = dispRs, y = shift), inherit.aes = FALSE) 
  

## smallest clim velo looks weird!!
## i think maybe the problem happens when clim velocity is less than 1?


## filter to one with slowest clim velo that's greater than 1 
gt_1 <- climVelo[which(climVelo >= 1)]

study = filter(fake_data, climVelo == min(gt_1))
n = 50
disp_vec <- modelr::seq_range(study$dispRs, n)

breakpoint_generate_quantities_quantities <- breakpoint_generate_quantities$generate_quantities(fitted_params = breakpoint_all_studies_posterior,
                                                                                                parallel_chains = 4, 
                                                                                                data = list(disp = disp_vec,
                                                                                                            n = n, 
                                                                                                            climVelo = rep(unique(study$climVelo), n))) 

breakpoint_generate_quantities_quantities %>% 
  gather_rvars(shift_pred[i]) 

breakpoint_generate_quantities_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = dispRs, y = shift), inherit.aes = FALSE) 

## zoom in 
breakpoint_generate_quantities_quantities %>%
  gather_rvars(shift_pred[i]) %>%
  mutate(disp = disp_vec[i]) %>%
  ggplot(aes(x = disp, dist = .value)) +
  stat_lineribbon() +
  geom_point(data = study, aes(x = dispRs, y = shift), inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))




