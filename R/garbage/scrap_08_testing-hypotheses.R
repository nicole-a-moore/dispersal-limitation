## TESTING EXPANSIONS
## - fit models to range expansions and compare 
##    rs ~ dd
##    rs ~ cv
##    rs ~ min(dd, cv)

##    rs ~ slope*dd + int
##    rs ~ slope*cv + int 
##    rs ~ slope*min(dd, cv) + int 

##    rs ~ slope*proxy trait + int

## macroecological proxy traits: body size, range size


##    expect:
##    - under hyp 1 (proxy traits poorly estimate dispersal): dispersal rate model has slope of 1 and is better fit than proxy trait models 
##    - under hyp 2 (some spp not dispersal limited): limiting rate model has slope of 1 & is better fit than climate velocity or dispersal rate alone

## - test whether results are robust:

## TESTING CONTRACTIONS
##    -  modelling range contractions
##       expect:
##       - hyp 1 (traits poorly estimate dispersal): dispersal rate model is not a good fit
##       - hyp 2 (some spp not dispersal limited): limiting rate model does not have a slope of 1 (although.. it might because of climate velocity) and is not better fit than climate velocity alone

## RANDOMIZATION
##    - randomization
##      - sample range shift rate randomly without replacement 1000 times 
##      - fit all models and compare the distribution of their coefficients to the real ones

## ERROR
##    - in climate velocity
##    - in dispersal

##    - study versus species-specific climate velocity 

##    - methods (e.g., grain size/study area - bird study has bigger shifts and bigger grain size)

##    - including centroid shifts 
##      - remove them and see what happens 
##      - or just forget about them??


## Qs for Jenn:
## treatment of expansions and contractions 
##  - bounding of data by 0
## models and hypothesis testing framework
## structure in residual variation

###################################
##       TESTING EXPANSIONS      ##
###################################
## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

## get rid of contractions
dd = filter(dd, is_contraction %in% c("NO"))

## filter to leading edge shifts with positive climate velocity 
dd <- filter(dd, Param == "LE" & ClimVeloTKmY_spp >= 0)

## standardize so CV > 0 means away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ClimVeloTKmY_spp = abs(dd$ClimVeloTKmY_spp)
dd$ShiftKmY = abs(dd$ShiftKmY)

## add gradient 
dd <- dd %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

hist(dd$ClimVeloTKmY_study)
hist(dd$ClimVeloTKmY_spp)

lat <- filter(dd, Type == "LAT") %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_study,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_study)) %>%
  mutate(ClimVeloTKmY_study = abs(ClimVeloTKmY_study)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_study, NA)) %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

## look at distribution of response variable 
hist(lat$ShiftKmY)
## bounded by 0 
## use a gamma distribution? 

mod <- glm(ShiftKmY ~ LimitingRate, 
           family = Gamma(link="log"),
           data = lat)
summary(mod, dispersion = 1)

## relationship with DispersalPotentialKmY definitely not well-represented by a linear model - what to do?

mod_disp <- lm(ShiftKmY ~ DispersalPotentialKmY, data = lat)
mod_climvelo <- lm(ShiftKmY ~ ClimVeloTKmY_study, data = lat)
mod_limrate <- lm(ShiftKmY ~ LimitingRate, data = lat)

## model summary
summary(mod_disp)
summary(mod_climvelo)
summary(mod_limrate)

## r squared
summary(mod_disp)$r.squared
summary(mod_climvelo)$r.squared
summary(mod_limrate)$r.squared

## adjusted r squared
summary(mod_disp)$adj.r.squared
summary(mod_climvelo)$adj.r.squared
summary(mod_limrate)$adj.r.squared

AIC(mod_disp, mod_climvelo, mod_limrate)

## plot residuals 
plot(mod_disp)
plot(mod_climvelo)
plot(mod_limrate)

hist(residuals(mod_disp))
hist(residuals(mod_climvelo))
hist(residuals(mod_limrate))
## uh oh non-normal

## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(lat$DispersalPotentialKmY),
                                                  by = 0.001))

df_disp$pred <- predict(mod_disp, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_study = seq(min(lat$ClimVeloTKmY_study), max(lat$ClimVeloTKmY_study),
                                                  by = 0.001))

df_climvelo$pred <- predict(mod_climvelo, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(lat$LimitingRate), max(lat$LimitingRate),
                                                   by = 0.001))

df_limrate$pred <- predict(mod_limrate, se.fit = FALSE, newdata = df_limrate)

## plot
lat %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_study)) +
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
  scale_y_continuous(limits = c(0, 40), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~is_contraction) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred), inherit.aes = FALSE)

lat %>%
  ggplot(aes(x = ClimVeloTKmY_study, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_study, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Mean climate velocity (km/y)",
       y = "Range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_study, y = pred), inherit.aes = FALSE)

lat %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
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
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred), inherit.aes = FALSE)

## methodological effects
lat %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ID, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
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
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred), inherit.aes = FALSE)

## study ID as a random effect?
library(nlme)
lme_disp <- lme(ShiftKmY ~ DispersalPotentialKmY,
                random = ~1|ID, 
                data = lat)
lme_climvelo <- lme(ShiftKmY ~ ClimVeloTKmY_study,
                    random = ~1|ID, 
                    data = lat)
lme_limrate <- lme(ShiftKmY ~ LimitingRate,
                random = ~1|ID, 
                data = lat)

summary(lme_disp)
summary(lme_climvelo)
summary(lme_limrate)

AIC(lme_disp, lme_climvelo, lme_limrate)

r.squaredGLMM(lme_disp)
r.squaredGLMM(lme_climvelo)
r.squaredGLMM(lme_limrate)

#####################################
##       TESTING CONTRACTIONS      ##
#####################################
## model range contractions that are in line with climate velocity the same way
## should not explain range contraction as well as it explains range expansion

## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

## get rid of expansions
dd = filter(dd, is_contraction %in% c("UNKNOWN", "YES"))

## filter to contractions in same direction as cv 
dd <- filter(dd, tracking_climate == TRUE)

## standardize so CV > 0 means away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ShiftKmY = abs(dd$ShiftKmY)

lat <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_study,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_study)) %>%
  mutate(ClimVeloTKmY_study = abs(ClimVeloTKmY_study)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_study, NA)) %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

mod_disp <- lm(ShiftKmY ~ DispersalPotentialKmY, data = lat)
mod_climvelo <- lm(ShiftKmY ~ ClimVeloTKmY_study, data = lat)
mod_limrate <- lm(ShiftKmY ~ LimitingRate, data = lat)

## model summary
summary(mod_disp)
summary(mod_climvelo)
summary(mod_limrate)

## r squared
summary(mod_disp)$r.squared
summary(mod_climvelo)$r.squared
summary(mod_limrate)$r.squared

## adjusted r squared
summary(mod_disp)$adj.r.squared
summary(mod_climvelo)$adj.r.squared
summary(mod_limrate)$adj.r.squared

AIC(mod_disp, mod_climvelo, mod_limrate)


## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(lat$DispersalPotentialKmY),
                                                  by = 0.001))

df_disp$pred <- predict(mod_disp, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_study = seq(min(lat$ClimVeloTKmY_study), max(lat$ClimVeloTKmY_study),
                                                   by = 0.001))

df_climvelo$pred <- predict(mod_climvelo, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(lat$LimitingRate), max(lat$LimitingRate),
                                            by = 0.001))

df_limrate$pred <- predict(mod_limrate, se.fit = FALSE, newdata = df_limrate)

## plot
lat %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_study)) +
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
  scale_y_continuous(limits = c(0, 40), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~is_contraction) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred), inherit.aes = FALSE)

lat %>%
  ggplot(aes(x = ClimVeloTKmY_study, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(lat, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_study, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Mean climate velocity (km/y)",
       y = "Range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_study, y = pred), inherit.aes = FALSE)

lat %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_study, shape = group)) +
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
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred), inherit.aes = FALSE)

library(nlme)
lme_disp <- lme(ShiftKmY ~ DispersalPotentialKmY,
                random = ~1|ID, 
                data = lat)
lme_climvelo <- lme(ShiftKmY ~ ClimVeloTKmY_study,
                    random = ~1|ID, 
                    data = lat)
lme_limrate <- lme(ShiftKmY ~ LimitingRate,
                   random = ~1|ID, 
                   data = lat)

summary(lme_disp)
summary(lme_climvelo)
summary(lme_limrate)

AIC(lme_disp, lme_climvelo, lme_limrate)

r.squaredGLMM(lme_disp)
r.squaredGLMM(lme_climvelo)
r.squaredGLMM(lme_limrate)


###################################
##         RANDOMIZATION         ##
###################################
## randomization test
#########################
## randomize the y variable 1000 times
## fit the limiting slope model and extract model coefficients 
## compare to observed coefficients 

## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

## get rid of contractions
dd = filter(dd, is_contraction %in% c("UNKNOWN", "NO"))

## filter to expansions in same direction as cv 
dd <- filter(dd, tracking_climate == TRUE)

## standardize so CV > 0 means away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ShiftKmY = abs(dd$ShiftKmY)

## add gradient 
dd <- dd %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

hist(dd$ClimVeloTKmY_study)
hist(dd$ClimVeloTKmY_spp)

lat <- filter(dd, Type == "LAT") %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_study,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_study)) %>%
  mutate(ClimVeloTKmY_study = abs(ClimVeloTKmY_study)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_study, NA)) %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 


## randomize range shift rate 
coeffs <- c()
for(i in 1:1000) {
  
  ## randomize shifts 
  rando_df <- data.frame(LimitingRate = lat$LimitingRate,
                         ClimVeloTKmY_study = lat$ClimVeloTKmY_study, 
                         DispersalPotentialKmY = lat$DispersalPotentialKmY, 
                         rando_shifts = sample(lat$ShiftKmY, size = nrow(lat), replace = FALSE))
  
  ## fit linear model
  lm1 <- lm(rando_shifts ~ LimitingRate,
            data = rando_df)
  lm2 <- lm(rando_shifts ~ ClimVeloTKmY_study,
            data = rando_df)
  lm3 <- lm(rando_shifts ~ DispersalPotentialKmY,
            data = rando_df)
  
  mod_summary1 = summary(lm1)
  mod_summary2 = summary(lm2)
  mod_summary3 = summary(lm3)
  
  ## extract info
  slopes = c(mod_summary1$coefficients[2,1], mod_summary2$coefficients[2,1], mod_summary3$coefficients[2,1])
  intercepts = c(mod_summary1$coefficients[1,1], mod_summary2$coefficients[1,1], mod_summary3$coefficients[1,1])
  r_squ = c(mod_summary1$r.squared, mod_summary2$r.squared, mod_summary3$r.squared)
  adj_r_squ = c(mod_summary1$adj.r.squared, mod_summary2$adj.r.squared, mod_summary3$adj.r.squared)
  
  ## get and save coeffs
  coeffs <- rbind(coeffs, 
                  data.frame(mod_num = rep(i, 3),
                             mod_type = c("LimitingRate", "ClimVeloTKmY_study", "DispersalPotentialKmY"),
                             slope = slopes,
                             intercept = intercepts,
                             r_squared = r_squ,
                             adj_r_squared = adj_r_squ))
  print(i)
}

## fit models to real data
mod_limrate <- lm(ShiftKmY ~ LimitingRate, data = lat)
sum_limrate <- summary(mod_limrate)

mod_climvelo <- lm(ShiftKmY ~ ClimVeloTKmY_study, data = lat)
sum_climvelo <- summary(mod_climvelo)

mod_disp <- lm(ShiftKmY ~ DispersalPotentialKmY, data = lat)
sum_disp <- summary(mod_disp)

## extract coefficients and add to database
coeffs <- coeffs %>%
  mutate(real_slope = ifelse(mod_type == "LimitingRate", sum_limrate$coefficients[2,1],
                             ifelse(mod_type == "ClimVeloTKmY_study", sum_climvelo$coefficients[2,1],
                                    sum_disp$coefficients[2,1])),
         real_intercept = ifelse(mod_type == "LimitingRate", sum_limrate$coefficients[1,1],
                                 ifelse(mod_type == "ClimVeloTKmY_study", sum_climvelo$coefficients[1,1],
                                        sum_disp$coefficients[1,1])),
         real_r_squ = ifelse(mod_type == "LimitingRate", sum_limrate$r.squared,
                             ifelse(mod_type == "ClimVeloTKmY_study", sum_climvelo$r.squared,
                                    sum_disp$r.squared)),
         real_adj_r_squ = ifelse(mod_type == "LimitingRate", sum_limrate$adj.r.squared,
                             ifelse(mod_type == "ClimVeloTKmY_study", sum_climvelo$adj.r.squared,
                                    sum_disp$adj.r.squared)))

## slope
coeffs %>%
  ggplot(aes(x = slope)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_slope), colour = "red") +
  facet_wrap(~mod_type)

## intercept
coeffs %>%
  ggplot(aes(x = intercept)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_intercept), colour = "red") +
  facet_wrap(~mod_type)

## r squared
coeffs %>%
  ggplot(aes(x = r_squared)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_r_squ), colour = "red") +
  facet_wrap(~mod_type)

## adj r squared
coeffs %>%
  ggplot(aes(x = adj_r_squared)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_adj_r_squ), colour = "red") +
  facet_wrap(~mod_type)



#### body size model
bodysize <- read.csv("data-processed/dispersal-proxy-trait-compilation.csv")

## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## get rid of centroid shifts
dd = filter(dd, Param != "O")

## get rid of contractions
dd = filter(dd, is_contraction %in% c("UNKNOWN", "NO"))

## filter to expansions in same direction as cv 
dd <- filter(dd, tracking_climate == TRUE)

## standardize so CV > 0 means away from range centre 
dd$ClimVeloTKmY_study = abs(dd$ClimVeloTKmY_study)
dd$ShiftKmY = abs(dd$ShiftKmY)

## add gradient 
dd <- dd %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

hist(dd$ClimVeloTKmY_study)
hist(dd$ClimVeloTKmY_spp)

lat <- filter(dd, Type == "LAT") %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_study,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_study)) %>%
  mutate(ClimVeloTKmY_study = abs(ClimVeloTKmY_study)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_study, NA)) %>%
  mutate(Gradient = ifelse(Type == "ELE", "Elevation", "Latitudinal")) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) 

bs <- bodysize %>%
  select(BodySize, scientificName) %>%
  left_join(lat, .) %>%
  select(BodySize, everything()) %>%
  filter(!is.na(BodySize))

mod_bs <- lm(ShiftKmY ~ log(BodySize), data = bs)

summary(mod_bs)
summary(mod_bs)$r.squared
summary(mod_bs)$adj.r.squared

## predict using each model
df_bs <- data.frame(BodySize = log(seq(min(bs$BodySize), max(bs$BodySize),
                                                  by = 0.001)))

df_bs$pred <- predict(mod_bs, se.fit = FALSE, newdata = df_bs)

## plot
bs %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, colour = ClimVeloTKmY_study)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(bs, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = BodySize, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0)) +
  labs(x = "Body size (m)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~is_contraction) +
  geom_line(data = df_bs, aes(x = exp(BodySize), y = pred), inherit.aes = FALSE)


dd %>%
  filter(ID == "A181_P1") %>%
  ggplot(aes(x = Rate)) + 
  geom_histogram()



