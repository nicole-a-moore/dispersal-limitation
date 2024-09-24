## trying nls modelling 
library(tidyverse)
library(AICcmodavg)

########################
## read in the data
########################
dd <- read.csv("data-processed/v1_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Gradient == "Latitudinal")

## read in proxy trait data 
proxyt <- read.csv("data-processed/dispersal-proxy-trait-compilation.csv") %>%
  select(scientificName, kingdom, BodySize)

## read in other  traits 
tamme = read.csv("data-processed/dispersal-proxy-trait-compilation_tamme.csv") %>%
  select(scientificName, kingdom, GrowthForm, DispersalSyndrome, 
         SeedWeightmg, SeedReleaseHeightm)
avonet = read.csv("data-processed/dispersal-proxy-trait-compilation_avonet.csv") %>%
  select(scientificName, kingdom, HandWingIndex, KippsDistance)

## join
dd <- left_join(dd, proxyt, relationship = "many-to-many") %>%
  left_join(., avonet, relationship = "many-to-many") %>%
  left_join(., tamme,  relationship = "many-to-many")

dd$min = ifelse(dd$DispersalPotentialKmY < dd$ClimVeloTKmY, dd$DispersalPotentialKmY, dd$ClimVeloTKmY)


#########################################
## compete a minimum-is-the-maximum model
#########################################
## subset to data with body size 
data <- dd[which(!is.na(dd$BodySize)),]

## model where range shift = climate velocity 
fit_clim_velo <- nls(formula = ShiftKmY ~ ClimVeloTKmY + error,
                     data = data, 
                     start = list(error = 0),
                     algorithm = "port")

summary(fit_clim_velo)


## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit_logic_kmy <- nls(formula = ShiftKmY ~ min + int,
           data = data, 
           start = list(int = 0),
           algorithm = "port")

summary(fit_logic_kmy)

# fit_logic_kmy_slope <- nls(formula =  ShiftKmY ~ 
#                              ifelse(DispersalPotentialKmY < ClimVeloTKmY, 
#                                     DispersalPotentialKmY*m + int, 
#                                     ClimVeloTKmY*m + int),
#                      data = data, 
#                      start = list(m = 1, int = 0),
#                      algorithm = "port")
# 
# summary(fit_logic_kmy_slope)

fit_logic_km <- nls(formula =  ShiftKmY ~ ifelse(DispersalDistanceKm < ClimVeloTKmY, 
                                              DispersalDistanceKm + int,
                                              ClimVeloTKmY + int),
           data = data, 
           start = list(int = 0),
           algorithm = "port")


summary(fit_logic_kmy)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit_proxy_disp_dist <- nls(formula =  ShiftKmY ~ m*DispersalDistanceKm + int,
            data = data, 
            start = list(m = 1, int = 0),
            algorithm = "port")

summary(fit_proxy_disp_dist)
## should be the same as this:
# lm_disp_dist <- lm(ShiftKmY ~ DispersalDistanceKm, data)
# summary(lm_disp_dist)

## write model where range shift rate = estimated slope*body size + error
## (this is a linear model where range shift ~ body size)
fit_proxy_bodysize <- nls(formula =  ShiftKmY ~ m*log(BodySize) + int,
            data = data, 
            start = list(m = 1, int = 0),
            algorithm = "port")

summary(fit_proxy_bodysize)
## should be the same as this:
# lm_bodysize <- lm(ShiftKmY ~ log(BodySize), data)
# summary(lm_bodysize)


## compete the models 
models <- list(fit_clim_velo, fit_logic_kmy, fit_logic_km, fit_proxy_disp_dist, fit_proxy_bodysize)
model.names <- c('clim_velo', 'logic_kmy', 'logic_km', 'proxy_disp_dist', 'proxy_body_size')

aictab(cand.set = models, modnames = model.names)
## YAY!!!!
## logical model performs better than just climate velocity, dispersal as a distance, and linear models where dispersal distance and body size explain range shift rate 


## plot the diagnostics 
## are relationships between proxy traits and range shifts in the expected direction? 
plot(fit_clim_velo)
hist(residuals(fit_clim_velo))
plot(x = data$ClimVeloTKmY, y = data$ShiftKmY)

plot(fit_logic_kmy)
hist(residuals(fit_logic_kmy))
plot(x = data$min, y = data$ShiftKmY)

plot(fit_logic_km)
hist(residuals(fit_logic_km))
plot(x = ifelse(data$DispersalDistanceKm < data$ClimVeloTKmY, data$DispersalDistanceKm, data$ClimVeloTKmY),
     y = data$ShiftKm)

plot(fit_proxy_bodysize)
hist(residuals(fit_proxy_bodysize))
plot(x = log(data$BodySize), y = data$ShiftKmY)

plot(fit_proxy_disp_dist)
hist(residuals(fit_proxy_disp_dist))
plot(x = data$DispersalDistanceKm, y = data$ShiftKmY)


## calculate r2
## r2 = 1 - ss residuals/ ss total
ss_res = sum(residuals(fit_logic_kmy)^2)
ss_tot = sum((data$ShiftKmY - mean(data$ShiftKmY))^2)
r2 = 1 - ss_res/ss_tot
r2
## very low


##########################################################
## check whether residual variation is greater among birds  
##########################################################
data$resid = residuals(fit_logic_kmy)

data %>%
  ggplot(aes(x = group, y = resid)) + geom_boxplot() +
  theme_bw()

data %>%
  group_by(group) %>%
  mutate(n = length(unique(scientificName))) %>%
  mutate(avg_res = sum(resid^2)/n) %>%
  select(avg_res, group) %>%
  distinct()




########################
## compete bird models
########################
## subset to plant data
birds <- dd[which(group == "Birds"),]
## with body size
birds <- birds[which(!is.na(birds$BodySize)),]

## model where range shift = climate velocity 
fit_birds_clim_velo <- nls(formula = ShiftKmY ~ ClimVeloTKmY + error,
                            data = birds, 
                            start = list(error = 0),
                            algorithm = "port")

summary(fit_birds_clim_velo)

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit_birds_logic_kmy <- nls(formula = ShiftKmY ~ min + int,
                            data = birds, 
                            start = list(int = 0),
                            algorithm = "port")

summary(fit_birds_logic_kmy)

fit_birds_logic_km <- nls(formula =  ShiftKmY ~ ifelse(DispersalDistanceKm < ClimVeloTKmY, 
                                                        DispersalDistanceKm + int,
                                                        ClimVeloTKmY + int),
                           data = birds, 
                           start = list(int = 0),
                           algorithm = "port")


summary(fit_birds_logic_km)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit_birds_proxy_disp_dist <- nls(formula =  ShiftKmY ~ m*DispersalDistanceKm + int,
                                  data = birds, 
                                  start = list(m = 1, int = 0),
                                  algorithm = "port")

summary(fit_birds_proxy_disp_dist)

## write model where range shift rate = estimated slope*body size + error
## (this is a linear model where range shift ~ body size)
fit_birds_proxy_bodysize <- nls(formula =  ShiftKmY ~ m*log(BodySize) + int,
                                 data = birds, 
                                 start = list(m = 1, int = 0),
                                 algorithm = "port")

summary(fit_birds_proxy_bodysize)


## compete the models 
models <- list(fit_birds_clim_velo, fit_birds_logic_kmy, fit_birds_logic_km, 
               fit_birds_proxy_disp_dist, fit_birds_proxy_bodysize)
model.names <- c('clim_velo', 'logic_kmy', 'logic_km', 'proxy_disp_dist', 'proxy_body_size')

aictab(cand.set = models, modnames = model.names)
## proxy traits are best fit but slope of relationship is 0

## plot the diagnostics 
## are relationships between proxy traits and range shifts in the expected direction? 
plot(fit_birds_clim_velo)
hist(residuals(fit_birds_clim_velo))
plot(x = birds$ClimVeloTKmY, y = birds$ShiftKmY)

plot(fit_birds_logic_kmy)
hist(residuals(fit_birds_logic_kmy))
plot(x = birds$min, y = birds$ShiftKmY)

plot(fit_birds_logic_km)
hist(residuals(fit_birds_logic_km))
plot(x = ifelse(birds$DispersalDistanceKm < birds$ClimVeloTKmY, birds$DispersalDistanceKm, birds$ClimVeloTKmY),
     y = birds$ShiftKm)

plot(fit_birds_proxy_bodysize)
hist(residuals(fit_birds_proxy_bodysize))
plot(x = log(birds$BodySize), y = birds$ShiftKmY)

plot(fit_birds_proxy_disp_dist)
hist(residuals(fit_birds_proxy_disp_dist))
plot(x = birds$DispersalDistanceKm, y = birds$ShiftKmY)

########################
## compete plant models
########################
## subset to plant data
plants <- dd[which(group == "Plants"),]
## with body size
plants <- plants[which(!is.na(plants$BodySize)),]

## model where range shift = climate velocity 
fit_plants_clim_velo <- nls(formula = ShiftKmY ~ ClimVeloTKmY + error,
                     data = plants, 
                     start = list(error = 0),
                     algorithm = "port")

summary(fit_plants_clim_velo)

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit_plants_logic_kmy <- nls(formula = ShiftKmY ~ min + int,
                     data = plants, 
                     start = list(int = 0),
                     algorithm = "port")

summary(fit_plants_logic_kmy)

fit_plants_logic_km <- nls(formula =  ShiftKmY ~ ifelse(DispersalDistanceKm < ClimVeloTKmY, 
                                                 DispersalDistanceKm + int,
                                                 ClimVeloTKmY + int),
                    data = plants, 
                    start = list(int = 0),
                    algorithm = "port")


summary(fit_plants_logic_km)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit_plants_proxy_disp_dist <- nls(formula =  ShiftKmY ~ m*DispersalDistanceKm + int,
                           data = plants, 
                           start = list(m = 1, int = 0),
                           algorithm = "port")

summary(fit_plants_proxy_disp_dist)

## write model where range shift rate = estimated slope*body size + error
## (this is a linear model where range shift ~ body size)
fit_plants_proxy_bodysize <- nls(formula =  ShiftKmY ~ m*log(BodySize) + int,
                          data = plants, 
                          start = list(m = 1, int = 0),
                          algorithm = "port")

summary(fit_plants_proxy_bodysize)


## compete the models 
models <- list(fit_plants_clim_velo, fit_plants_logic_kmy, fit_plants_logic_km, 
               fit_plants_proxy_disp_dist, fit_plants_proxy_bodysize)
model.names <- c('clim_velo', 'logic_kmy', 'logic_km', 'proxy_disp_dist', 'proxy_body_size')

aictab(cand.set = models, modnames = model.names)

## plot the diagnostics 
## are relationships between proxy traits and range shifts in the expected direction? 
plot(fit_plants_clim_velo)
hist(residuals(fit_plants_clim_velo))
plot(x = plants$ClimVeloTKmY, y = plants$ShiftKmY)

plot(fit_plants_logic_kmy)
hist(residuals(fit_plants_logic_kmy))
plot(x = plants$min, y = plants$ShiftKmY)

plot(fit_plants_logic_km)
hist(residuals(fit_plants_logic_km))
plot(x = ifelse(plants$DispersalDistanceKm < plants$ClimVeloTKmY, plants$DispersalDistanceKm, plants$ClimVeloTKmY),
     y = plants$ShiftKm)

plot(fit_plants_proxy_bodysize)
hist(residuals(fit_plants_proxy_bodysize))
plot(x = log(plants$BodySize), y = plants$ShiftKmY)

plot(fit_plants_proxy_disp_dist)
hist(residuals(fit_plants_proxy_disp_dist))
plot(x = plants$DispersalDistanceKm, y = plants$ShiftKmY)


## plot 
coefs = coefficients(fit_logic_kmy)

ggdata %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x + coefs[1]}) + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
 # scale_x_log10() +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 

ggdata %>%
  ggplot(aes(x = clim_velo, y = shift)) + geom_point() +
  theme_bw() +
  # scale_x_log10() +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 

ggdata %>%
  ggplot(aes(x = disp_rate, y = shift)) + geom_point() +
  theme_bw() +
  # scale_x_log10() +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 

coefs_birds = coefficients(fit_birds_logic_kmy)
coefs_plants = coefficients(fit_plants_logic_kmy)

ggdata %>%
  filter(group == "Birds") %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x + coefs_birds[1]}) + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
 # scale_x_log10() +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 

ggdata %>%
  filter(group == "Plants") %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x + coefs_plants[1]}) + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  #scale_x_log10() +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 



##########################################################
## check whether plant dispersal rates > range shift rates  
##########################################################
co = coefficients(fit_plants_logic_kmy)[[1]]

plants %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, 
             colour = ShiftKmY > DispersalPotentialKmY)) + geom_point() +
  theme_bw() +
  scale_x_log10() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") +
  stat_function(colour = "black", fun = function(x){x + co[1]}) + 
  scale_y_continuous(limits = c(0, 20)) +
  labs(colour = "Shifting faster\nthan dispersal rate:") 

summary(fit_plants_logic_kmy)
plot(fit_plants_logic_kmy)
hist(residuals(fit_plants_logic_kmy))
plot(x = plants$min, y = plants$ShiftKmY)



##########################################################
## subset to only leading edge shifts   
##########################################################
## filter to leading edge only 
le <- data[which(data$Position == "Leading edge"),]

## model where range shift = climate velocity 
fit_clim_velo <- nls(formula = ShiftKmY ~ ClimVeloTKmY + error,
                     data = le, 
                     start = list(error = 0),
                     algorithm = "port")

summary(fit_clim_velo)


## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit_logic_kmy <- nls(formula = ShiftKmY ~ min + int,
                     data = le, 
                     start = list(int = 0),
                     algorithm = "port")

summary(fit_logic_kmy)

# fit_logic_kmy_slope <- nls(formula =  ShiftKmY ~ 
#                              ifelse(DispersalPotentialKmY < ClimVeloTKmY, 
#                                     DispersalPotentialKmY*m + int, 
#                                     ClimVeloTKmY*m + int),
#                      data = data, 
#                      start = list(m = 1, int = 0),
#                      algorithm = "port")
# 
# summary(fit_logic_kmy_slope)

fit_logic_km <- nls(formula =  ShiftKmY ~ ifelse(DispersalDistanceKm < ClimVeloTKmY, 
                                                 DispersalDistanceKm + int,
                                                 ClimVeloTKmY + int),
                    data = le, 
                    start = list(int = 0),
                    algorithm = "port")


summary(fit_logic_kmy)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit_proxy_disp_dist <- nls(formula =  ShiftKmY ~ m*DispersalDistanceKm + int,
                           data = le, 
                           start = list(m = 1, int = 0),
                           algorithm = "port")

summary(fit_proxy_disp_dist)
## should be the same as this:
# lm_disp_dist <- lm(ShiftKmY ~ DispersalDistanceKm, data)
# summary(lm_disp_dist)

## write model where range shift rate = estimated slope*body size + error
## (this is a linear model where range shift ~ body size)
fit_proxy_bodysize <- nls(formula =  ShiftKmY ~ m*log(BodySize) + int,
                          data = le, 
                          start = list(m = 1, int = 0),
                          algorithm = "port")

summary(fit_proxy_bodysize)
## should be the same as this:
# lm_bodysize <- lm(ShiftKmY ~ log(BodySize), data)
# summary(lm_bodysize)


## compete the models 
models <- list(fit_clim_velo, fit_logic_kmy, fit_logic_km, fit_proxy_disp_dist, fit_proxy_bodysize)
model.names <- c('clim_velo', 'logic_kmy', 'logic_km', 'proxy_disp_dist', 'proxy_body_size')

aictab(cand.set = models, modnames = model.names)
## YAY!!!!
## logical model performs better than just climate velocity, dispersal as a distance, and linear models where dispersal distance and body size explain range shift rate 


## plot the diagnostics 
## are relationships between proxy traits and range shifts in the expected direction? 
plot(fit_clim_velo)
hist(residuals(fit_clim_velo))
plot(x = le$ClimVeloTKmY, y = le$ShiftKmY)

plot(fit_logic_kmy)
hist(residuals(fit_logic_kmy))
plot(x = le$min, y = le$ShiftKmY)

plot(fit_logic_km)
hist(residuals(fit_logic_km))
plot(x = ifelse(le$DispersalDistanceKm < le$ClimVeloTKmY, le$DispersalDistanceKm, le$ClimVeloTKmY),
     y = le$ShiftKm)

plot(fit_proxy_bodysize)
hist(residuals(fit_proxy_bodysize))
plot(x = log(le$BodySize), y = le$ShiftKmY)

plot(fit_proxy_disp_dist)
hist(residuals(fit_proxy_disp_dist))
plot(x = le$DispersalDistanceKm, y = le$ShiftKmY)

## dispersal distance model is better at explaining variation, but slope is 0

## things for later:
## what to do about models with non-normal residuals? 

## calculate r2
## r2 = 1 - ss residuals/ ss total
ss_res = sum(residuals(fit_logic_kmy)^2)
ss_tot = sum((le$ShiftKmY - mean(le$ShiftKmY))^2)
r2 = 1 - ss_res/ss_tot
r2
## very low - but higher than with centroid 



#### next:
#### see if linear model can be fit without variables slope 






##########################################################
## garbage
##########################################################


## try manually fitting model without parameters
## predicted values:
pred = ifelse(data$disp_rate < data$clim_velo, data$disp_rate, data$clim_velo)

## calculate residuals:
res = data$shift - pred
hist(res)

## calculate residual sum of squares 
rss = sum(res^2)

## calculate AIC:
n = length(data$shift)
AIC_nopar = 2*1 + n*(log(rss/n))
AIC = 2*2 + n*(log(sum(residuals(fit)^2)/n))
AIC2 = 2*2 + n*(log(sum(residuals(fit2)^2)/n))
AIC3 = 2*2 + n*(log(sum(residuals(fit3)^2)/n))

## the no parameter model doesn't do very well 



## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*growth_form_num + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ growth_form_num, data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'growth_form')

aictab(cand.set = models, modnames = model.names)

plot(x = data_sub$growth_form_num, y = data_sub$shift, )
levels <- levels(as.factor(data_sub$growth_form))

## subset to data with growth form
data_sub <- data[which(!is.na(seed_release_height)),]

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*log(seed_release_height) + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ log(seed_release_height), data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'seed_release_height')

aictab(cand.set = models, modnames = model.names)

plot(x = log(data_sub$seed_release_height), y = data_sub$shift)


## subset to data with growth form
data_sub <- data[which(!is.na(seed_weight)),]

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*log(seed_weight) + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ log(seed_weight), data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'seed_weight')

aictab(cand.set = models, modnames = model.names)

plot(x = log(data_sub$seed_weight), y = data_sub$shift)



## subset to data with dispersal syndrome
data_sub <- data[which(!is.na(dispersal_syndrome)),]

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*dispersal_syndrome_num + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ dispersal_syndrome_num, data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'dispersal_syndrome')

aictab(cand.set = models, modnames = model.names)

plot(x = log(data_sub$dispersal_syndrome_num), y = data_sub$shift)
levels <- levels(as.factor(data_sub$dispersal_syndrome))




## now for birds


disp_rate <- dd$DispersalPotentialKmY
clim_velo <- dd$ClimVeloTKmY
shift <- dd$ShiftKmY
group <- dd$group # as.numeric(as.factor(dd$group)) -1
bodysize <- dd$BodySize
disp_dist <- dd$DispersalDistanceKm
min = ifelse(disp_rate < clim_velo, disp_rate, clim_velo)
which_min = as.numeric(disp_rate == min)
growth_form = dd$GrowthForm
growth_form_num = as.numeric(as.factor(dd$GrowthForm)) - 1
dispersal_syndrome = dd$DispersalSyndrome
dispersal_syndrome_num = as.numeric(as.factor(dd$DispersalSyndrome)) - 1
seed_release_height = dd$SeedReleaseHeightm
seed_weight = dd$SeedWeightmg
hwi = dd$HandWingIndex
kipps_dist = dd$KippsDistance

## make dataframe
data <- data.frame(disp_rate, clim_velo, shift, group, bodysize, disp_dist, min, which_min, 
                   growth_form, growth_form_num, seed_release_height, seed_weight, dispersal_syndrome,
                   dispersal_syndrome_num, hwi, kipps_dist)

## subset to data with hwi
data_sub <- data[which(!is.na(hwi)),]

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*hwi + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ hwi, data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'hwi')

aictab(cand.set = models, modnames = model.names)

## subset to data with kipps
data_sub <- data[which(!is.na(kipps_dist)),]

## write model where range shift rate = minimum of climate velocity and potential dipsersal rate + error
## error is the only estimated parameter
fit <- nls(formula =  shift ~ ifelse(disp_rate < clim_velo, ye = disp_rate + error, no = clim_velo + error),
           data = data_sub, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## write model where range shift rate = estimated slope*dispersal distance + error
## (this is a linear model where range shift ~ dispersal distance )
fit2 <- nls(formula =  shift ~ m*kipps_dist + error,
            data = data_sub, 
            start = list(m = 1, error = 0),
            algorithm = "port")

summary(fit2)
## should be the same as this:
lm <- lm(shift ~ hwi, data_sub)
summary(lm)

## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'kipps')

aictab(cand.set = models, modnames = model.names)




plot(x = data_sub$kipps_dist, y = data_sub$shift)
plot(x = data_sub$hwi, y = data_sub$shift)


plot(x = data_sub$min, y = data_sub$shift)

## seems like linear models explain more error, but slope ~ 0 (no expected positive relationship - except dispersal syndrome)






## mess around
fit <- nls(formula =  shift ~ min + error,
           data = data, 
           start = list(error = 0),
           algorithm = "port")

summary(fit)

## plot 
ggdata <- data.frame(shift = data$shift, 
                     logic_pred = data$min, disp_rate = data$disp_rate,
                     clim_velo = data$clim_velo, 
                     model_pred = data$min + 0.8452, group = data$group, 
                     what_limits = ifelse(disp_rate < clim_velo, "dispersal rate", 
                                          "climate velocity"))

## theoretically, shift = min
## our model allows shift = min + error 
ggdata %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point() +
  theme_bw() +
  geom_smooth(aes(x = logic_pred, y = model_pred), colour = "red", method = "lm") +
  geom_abline(intercept = 0, slope = 1)
## red line = model predictions 
## black line = 1:1 line
## distance from red line to points = residual variation


ggdata %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "red", fun = function(x){x + 0.8452}) + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  scale_x_log10(limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 


ggdata %>%
  ggplot(aes(x = disp_rate, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  scale_x_log10(limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40)) +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 


## try allowing slope to vary
fit2 <- nls(formula =  shift ~ min*slope + int,
            data = data, 
            start = list(int = 0, slope = 1),
            algorithm = "port")

summary(fit2)

ggdata <- data.frame(shift = data$shift, 
                     logic_pred = data$min, 
                     disp_rate = data$disp_rate,
                     clim_velo = data$clim_velo, 
                     model_pred = data$min*0.57352 + 1.36832, 
                     group = data$group, 
                     what_limits = ifelse(disp_rate < clim_velo, "dispersal rate", 
                                          "climate velocity"),
                     residuals = residuals(fit))

ggdata %>%
  ggplot(aes(x = logic_pred, y = shift)) + geom_point(aes(colour = what_limits)) +
  theme_bw() +
  stat_function(colour = "red", fun = function(x){x*0.57352 + 0.8452}) + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 40)) +
  labs(x = "Maximum theoretical shift rate (km/y)", "Observed shift rate (km/y)", 
       colour = "Theoretical limiting factor") 


## compete the models 
models <- list(fit, fit2)
model.names <- c('logic', 'logic_slope')

aictab(cand.set = models, modnames = model.names)


ggdata %>%
  ggplot(aes(x = abs(residuals), fill = group)) + geom_histogram()




