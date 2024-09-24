## running Andrew's code 
# https://aammd.github.io/posts/2023-07-05-breakpoint/
# https://github.com/aammd/aammd.github.io/tree/main/posts/2023-07-05-breakpoint
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)

############################
## Write that in mathematics
############################

## some things:
############################
## range shift speed (y variable) has normally distributed error 
## can we change this? climate velocity has a strong right skew (most places are experiencing slow change) and so would expect range shift speed to follow a similar distribution since if species are responding to climate velocity they should shift slower or the same speed 

## our range shift values (and climate velocity values) are 0-bounded since I kept only cases when we would expect range expansions and where range expansions occurred


###########################################
## Simple Stan code with prior predictions
###########################################
## read in stan model
breakpoint_B <- cmdstan_model("stan/breakpoint_B.stan")
breakpoint_B

## n = number of sites
## B = potential dispersal rate
## x[i] = climate velocity at site i 
## b2 = slope of the line when dispersal limitation is acting (should be 1)

## questions:
############################
## what is sigma?

## "model { b2 ~ normal(1, 1);
##          sigma ~ exponential(1);
##          }
## ^ does this bit set the prior distributions for the breakpoint and error around it? think so

## another thing:
############################
## in this model, B has only one value but we have one dispersal rate per species 
## in other words, this model is showing what data from one species might look like when exposed to different climate velocities

## sample x variable (climate velocity) from uniform distribution
xvar <- runif(42, min = 1, max=55)
hist(xvar)

## feed the model our fake data
breakpoint_B_prior <- breakpoint_B$sample(chains=1,
                                          data = list(x = xvar,
                                                      n = 42, ## 42 observations 
                                                      B = 25)) ## potential dispersal rate = 25


## look at the priors 
breakpoint_B_prior
## yay! making sense 
## b2 ~ 1 
## sigma ~ 1

## except:
## what is lp? log posterior 

## if we were happy with our model, what's the next step? how do we feed our data in? this generates y values, but we already know our y values ?


## plot samples from the prior 
prior_draws <- breakpoint_B_prior |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])

prior_draws |> 
  ggplot(aes(x = x, y = .value)) + geom_point() + 
  facet_wrap(~.draw)
## is B = 25 for each of these examples



###############################
## Rewriting this with step()
###############################
breakpoint_step_prior <- cmdstan_model("stan/breakpoint_step_prior.stan")
breakpoint_step_prior

## not sure exactly what the step function does, but let's assume it does what I suspect
## and takes on a value of 0 if argument is negative, 1 if positive

## if so, this model says:
## if dispersal rate is slower than climate velocity at site i (B>x), let x2 = 0
## if dispersal rate is faster than climate velocity at site i (B<x), let x2 = 1

## then, let mean of range shift rate by centred around dipsersal rate + (climate velocity - dispersal rate)*x2
## so if slower is slower than climate velocity, x2 = 0 and mean range shift rate = dispersal rate
## but if species is faster than climate velocity, x2 = 1 and mean range shift rate = dispersal rate + (climate velocity - dispersal rate)*1

## let's say climate velocity = 3
## dispersal rate = 5
## then x2 = 1
## and y = 5 + (3 - 5)*1 = 3

## and the opposite
## climate velocity = 5
## dispersal rate = 3
## then x2 = 0
## and y = 3 + (3 - 5)*0 = 3

## run it!
set.seed(4812)
## sample climate velocities 
xvar <- runif(42, min = 1, max=55)

breakpoint_step_prior_samples <- breakpoint_step_prior$sample(chains=1,
                                                              data = list(x = xvar,
                                                                          n = 42, 
                                                                          B = 10), ## set dispersal to 10
                                                              refresh = 0L)
breakpoint_step_prior_samples
## this model doesn't estimate the slope b2 

## plot draws 
prior_draws <- breakpoint_step_prior_samples |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])

prior_draws |> 
  ggplot(aes(x = x, y = .value)) + geom_point() + 
  facet_wrap(~.draw)


####################################
## thinking about how to adapt model
#####################################
theme_set(theme_minimal())

## read in data 
dd <- read.csv("data-processed/v1_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Gradient == "Latitudinal")

## each species = 1 potential dispersal rate
## each study area = 1 climate velocity 
## each species x study area = 1 range shift rate

## uneven sampling of climate velocities within species; sometimes, all areas species has been studied in have climate velocities slower than dispersal rate or the opposite

## how many species are sampled across multiple sites?
tally <-dd %>%
  select(scientificName, Source, LatVeloT) %>%
  distinct() %>%
  group_by(scientificName) %>%
  tally() %>%
  arrange(-n) %>%
  mutate(scientificName = factor(scientificName, levels = unique(.$scientificName), ordered = TRUE))

tally %>%
  ggplot(aes(x = scientificName, y = n)) + geom_col() +
  coord_flip()

dd %>%
  filter(scientificName %in% tally$scientificName[1:12]) %>%
  ggplot(aes(x = LatVeloT, y = ShiftKmY, colour = AnnualDispPotKmY)) +
  geom_point() +
  facet_wrap(~scientificName) + ## facet by species
  theme(legend.position = "none") +
  geom_abline(slope = 1, yintercept = 0) ## plot 1:1 line

dd %>%
  filter(scientificName %in% tally$scientificName[1:12]) %>%
  ggplot(aes(x = LatVeloT, y = ShiftKmY, colour = AnnualDispPotKmY)) +
  geom_point() +
  facet_wrap(~scientificName) + ## facet by species
  theme(legend.position = "none") +
  geom_abline(slope = 1, yintercept = 0) + ## plot 1:1 line
  geom_vline(aes(xintercept = AnnualDispPotKmY)) ## plot 1:1 line
  
## see, for these species all observations are in places where climate is slower than dispersal potential
temp = dd %>%
  filter(scientificName %in% tally$scientificName[1:12]) %>%
  mutate(diff = AnnualDispPotKmY - LatVeloT)
length(which(temp$diff < 0))

## how many species are sampled in at least 2 places where we expect dispersal limitation AND in at least 2 places where we don't?
dd %>%
  select(scientificName, Source, LatVeloT, AnnualDispPotKmY) %>%
  mutate(diff = AnnualDispPotKmY - LatVeloT) %>%
  distinct() %>%
  group_by(scientificName) %>%
  mutate(n_limited = length(which(diff < 0)),
         n_notlimited = length(which(diff > 0)),
         even_sample = n_limited >= 2 & n_notlimited >=2) %>%
  ungroup() %>%
  select(scientificName, even_sample) %>%
  distinct() %>%
  group_by(even_sample) %>% 
  tally()
## only 2 species - Pinus sylvestris, Anthriscus sylvestris

## look at them
dd %>%
  filter(scientificName %in% c("Pinus sylvestris", "Anthriscus sylvestris")) %>%
  ggplot(aes(x = LatVeloT, y = ShiftKmY, colour = AnnualDispPotKmY)) +
  geom_point() +
  facet_wrap(~scientificName) + ## facet by species
  theme(legend.position = "none") +
  geom_abline(slope = 1, yintercept = 0) + ## plot 1:1 line
  geom_vline(aes(xintercept = AnnualDispPotKmY)) ## plot 1:1 line
  

## could we just change the sample size for each species at each site to reflect the actual sample sizes we have?
## for example:
## for a species with dispersal potential = 25 sampled at 3 different sites
xvar <- runif(3, min = 1, max=55)
breakpoint_B_prior <- breakpoint_B$sample(chains=1,
                                          data = list(x = xvar,
                                                      n = 3, ## 3 observations 
                                                      B = 25)) ## potential dispersal rate = 25


## look at the priors 
breakpoint_B_prior
## still does a good job
## slope b2 is still ~1 with sigma = ~1

prior_draws <- breakpoint_B_prior |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])

prior_draws |> 
  ggplot(aes(x = x, y = .value)) + geom_point() + 
  facet_wrap(~.draw)


## what happens if we only have data sampled on one side of the theoretical breakpoint?
## for example:
## for a species with dispersal potential = 25 sampled at 3 different sites
## and we measured climate velocity at sites with more than 25 km/y climate velocity 
xvar <- runif(3, min = 30, max=55)
breakpoint_B_prior <- breakpoint_B$sample(chains=1,
                                          data = list(x = xvar,
                                                      n = 3, ## 3 observations 
                                                      B = 25)) ## potential dispersal rate = 25


## look at the priors 
breakpoint_B_prior
## somehow still works? why?

prior_draws <- breakpoint_B_prior |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])

prior_draws |> 
  ggplot(aes(x = x, y = .value)) + geom_point() + 
  facet_wrap(~.draw)

## try with real climate velocities

data = dd %>%
  filter(scientificName %in% c("Anthriscus sylvestris")) 

unique(data$DispersalPotentialKmY)
  
xvar <- data$LatVeloT
breakpoint_B_prior <- breakpoint_B$sample(chains=1,
                                          data = list(x = xvar,
                                                      n = 7, ## 7 observations 
                                                      B = 0.12)) ## potential dispersal rate = 0.12


## look at the priors 
breakpoint_B_prior
## somehow still works!!!!! what is this sorcery 

prior_draws <- breakpoint_B_prior |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])


## plot the real data on top:
prior_draws |> 
  ggplot(aes(x = x, y = .value)) + geom_point() + 
  facet_wrap(~.draw) +
  geom_point(data = data, aes(x = LatVeloT, y = ShiftKmY), inherit.aes = FALSE, colour = "red")
## red are the real data 


## overall questions:
## how to feed in our data to the model once we decide on one 
## and how to interpret the model results :-)

## what do we want to estimate about our data?
## slope before theoretical breakpoint - is it 1? 
## slope after theoretical breakpoint - is it 0?
## how close is theoretical to actual breakpoint? or, does actual breakpoint increase with climate velocity?

## can we even do this with our data? probably not enough data across sites + within species to estimate 2 slopes 


############################
## switch for gamma errors
############################
## dispersal rate can never be negative, so use a gamma distribution around the mean instead of a normal distribution

breakpoint_step_gamma_prior <- cmdstan_model("stan/breakpoint_step_gamma_prior.stan")
breakpoint_step_gamma_prior

set.seed(4812)
xvar <- runif(42, min = 1, max = 55)

breakpoint_step_gamma_prior_samples <- breakpoint_step_gamma_prior$sample(
  chains=1,
  data = list(x = xvar,
              n = 42, B = 33), refresh = 1000L)

prior_draws <- breakpoint_step_gamma_prior_samples |> 
  tidybayes::gather_draws(y[i], ndraws = 12) |> 
  mutate(x = xvar[i])

prior_draws |> 
  ggplot(aes(x = x, y = .value)) + 
  geom_point() + 
  facet_wrap(~.draw)


vec = rgamma(400, 1^2/0.1^2, 1/0.1^2)
sd(vec)
mean(vec)

hist(vec)
min(vec)



