library(AICcmodavg)
library(tidyverse)
library(ggplot2)
library(dotwhisker)
theme_set(theme_bw())

## TESTING EXPANSIONS
## filter to leading edge range shifts with positive climate velocity (expect expansion)
## - fit models to range expansions and compare 
##    rs ~ dd
##    rs ~ cv
##    rs ~ min(dd, cv)
# (...are these slopeless models really necessary? could just comment on whether estimated slope = 1)

##    rs ~ slope*dd + int
##    rs ~ slope*cv + int 
##    rs ~ slope*min(dd, cv) + int 

##    rs ~ slope*proxy trait + int

## macroecological proxy traits: body size, range size

##    expect:
##    - under hyp 1 (proxy traits poorly estimate dispersal): dispersal rate model has slope of 1 and is better fit than proxy trait models 
##    - under hyp 2 (some spp not dispersal limited): limiting rate model has slope of 1 & is better fit than climate velocity or dispersal rate alone

## ACCOUNTING FOR ERROR 
## see whether including slow, negative contractions (which might be improperly measured expansions) changes the result 

## ROBUSTNESS TO RANDOMIZATION
##    - randomization
##      - sample range shift expansions randomly without replacement 1000 times 
##      - fit all models and compare the distribution of their coefficients to the real ones

## TESTING CONTRACTIONS
##  - fit models to range contractions and compare 
##       expect:
##       - hyp 1 (traits poorly estimate dispersal): dispersal rate model is not a good fit
##       - hyp 2 (some spp not dispersal limited): limiting rate model does not have a slope of 1 (although.. it might because of climate velocity) and is not better fit than climate velocity alone

## OTHER THINGS
##    - error in climate velocity
##    - methods (e.g., grain size/study area - bird study has bigger shifts and bigger grain size)


## Qs for Jenn:
## bounding of data by 0 in the linear models 
## models and hypothesis testing framework
## structure in residual variation

#############################
##       PREPARE DATA      ##
#############################
## read in data 
dd <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## filter to leading edge shifts with positive climate velocity 
dd <- filter(dd, !Param %in% c("O", "TE") & ClimVeloTKmY_spp >= 0)

hist(dd$ClimVeloTKmY_spp)
hist(dd$ClimVeloTKmY_spp)

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_spp)) %>%
  mutate(ClimVeloTKmY_spp = abs(ClimVeloTKmY_spp)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_spp, NA)) 

## join to body size and range size data 
bodysize <- read.csv("data-processed/dispersal-proxy-trait-compilation.csv")
rangesize <- read.csv("data-raw/BIOSHIFTSv3/polygon_areas.csv")

## filter to species in dispersal dataset 
rangesize <- filter(rangesize, rangesize$binomial %in% dd$scientificName) %>%
  select(ID, binomial, range_source, Area_m2_ID, Area_m2_range) %>%
  distinct()

## select range size from one source
## create a list of priority for range sources:
range_sources <- c("IUCN", "BOTW", "GARD", "BIEN", "Fishmap", "Aquamaps", "Butterfly Maps", "GIFT", "GBIF occurrence")

rangesize <- rangesize %>%
  group_by(binomial, ID) %>%
  arrange(match(range_source, range_sources),.by_group = TRUE) %>%
  filter(row_number()==1) %>%
  ungroup()

## join traits to dd 
dd <- dd %>%
  select(-range_source) %>%
  left_join(., rangesize, by = c("ID", "scientificName" = "binomial")) %>%
  left_join(., select(bodysize, scientificName, class, kingdom, BodySize))


###################################
##       TESTING EXPANSIONS      ##
###################################
## get rid of range contractions
data <- filter(dd, is_contraction == "NO")

## save data
write.csv(data, "data-processed/model-data_expansions.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)
## bounded by 0 
## use a gamma distribution? 

# mod <- glm(ShiftKmY ~ LimitingRate, 
#            family = Gamma(link="log"),
#            data = data)
# summary(mod, dispersion = 1)

## redataionship with DispersalPotentialKmY definitely not well-represented by a linear model - what to do?

## fit models without slope parameter 
nls_disp <- nls(formula = ShiftKmY ~ DispersalPotentialKmY + int,
                data = data, 
                start = list(int = 0),
                algorithm = "port")

nls_climvelo <- nls(formula = ShiftKmY ~ ClimVeloTKmY_spp + int,
                data = data, 
                start = list(int = 0),
                algorithm = "port")

nls_limrate <- nls(formula = ShiftKmY ~ LimitingRate + int,
                    data = data, 
                    start = list(int = 0),
                    algorithm = "port")

## fit models with estimated slope parameter 
lm_disp <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                data = data, 
                start = list(int = 0, m = 1),
                algorithm = "port")

lm_climvelo <- nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                    data = data, 
                    start = list(int = 0, m = 1),
                    algorithm = "port")

lm_limrate <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                   data = data, 
                   start = list(int = 0, m = 1),
                   algorithm = "port")

## trait models
lm_dr <- nls(formula = ShiftKmY ~ m*log(DispersalPotentialKmY) + int,
             data = data, 
             start = list(int = 0, m = 1),
             algorithm = "port")
lm_bs <- nls(formula = ShiftKmY ~ m*log(BodySize) + int,
                  data = data, 
                  start = list(int = 0, m = 1),
                  algorithm = "port")
lm_rs <- nls(formula = ShiftKmY ~ m*log(Area_m2_range) + int,
             data = data, 
             start = list(int = 0, m = 0),
             algorithm = "port")


## save models 
expansion_models <- list(nls_disp, nls_climvelo, nls_limrate,
                        lm_disp, lm_climvelo, lm_limrate,
                        lm_dr, lm_bs, lm_rs)
names(expansion_models) <- c("nls_disp", "nls_climvelo", "nls_limrate",
                             "lm_disp", "lm_climvelo", "lm_limrate",
                             "lm_dr", "lm_bs", "lm_rs")

saveRDS(expansion_models, "data-processed/modelfits_expansions.rds")

## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                  by = 0.001))

bs = nlraa::boot_nls(nls_disp, R = 1e2, data = data)
confint(bs, type = "perc")

df_disp$pred_nls <- nlraa::predict2_nls(nls_disp, interval = "confidence", level = 0.95,
                            newdata = df_disp)[,1]


df_disp$pred_nls <- predict(nls_disp, se.fit = T, interval = "confidence", level = 0.95,
                                newdata = df_disp)[,1]
df_disp$pred_lm <- predict(lm_disp, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                  by = 0.001))

df_climvelo$pred_nls <- predict(nls_climvelo, se.fit = FALSE, newdata = df_climvelo)
df_climvelo$pred_lm <- predict(lm_climvelo, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                                   by = 0.001))

df_limrate$pred_nls <- predict(nls_limrate, se.fit = FALSE, newdata = df_limrate)
df_limrate$pred_lm <- predict(lm_limrate, se.fit = FALSE, newdata = df_limrate)

df_dr <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                by = 0.001))
df_bs <- data.frame(BodySize = seq(min(data$BodySize, na.rm = TRUE), max(data$BodySize, na.rm = TRUE),
                                   by = 0.001))
df_rs <- data.frame(Area_m2_range = seq(min(data$Area_m2_range, na.rm = TRUE), 
                                        max(data$Area_m2_range, na.rm = TRUE),
                                            by = 100000000))

df_dr$pred <- predict(lm_dr, se.fit = FALSE, newdata = df_dr)
df_bs$pred <- predict(lm_bs, se.fit = FALSE, newdata = df_bs)
df_rs$pred <- predict(lm_rs, se.fit = FALSE, newdata = df_rs)


## plot
## DISPERSAL
disp_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE) 
disp_plot_2 <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE) 

## CLIMATE
clim_velo_plot <- data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black",
             fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_lm), colour = "black", inherit.aes = FALSE)

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_nls), colour = "grey", 
            inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", 
            inherit.aes = FALSE)

plot_grid(disp_plot, clim_velo_plot, limrate_plot, 
         ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_expansions.png", 
       width = 12.5, height = 4)

## proxy 
## dispersal rate
dr_plot = data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_dr, aes(x = DispersalPotentialKmY, y = pred), colour = "black", inherit.aes = FALSE)

## body size
bs_plot = data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = BodySize, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  labs(x = "Body size (mm)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_bs, aes(x = BodySize, y = pred), colour = "black", inherit.aes = FALSE)

## range size
rs_plot = data %>%
  ggplot(aes(x = Area_m2_range, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = Area_m2_range, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 26), expand = c(0,0.5)) +
  labs(x = "Realized ange size (m2)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_rs, aes(x = Area_m2_range, y = pred), colour = "black", inherit.aes = FALSE)

plot_grid(bs_plot, rs_plot, dr_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_proxy-trait_expansions.png", 
       width = 12.5, height = 4)


#####################################
##       ACCOUNTING FOR ERROR      ##
#####################################
## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

## get rid of contractions faster than 5km/y
data <- filter(dd, Rate >= -5)

## save data
write.csv(data, "data-processed/model-data_error.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

## fit models without slope parameter 
nls_disp_error <- nls(formula = ShiftKmY ~ DispersalPotentialKmY + int,
                data = data, 
                start = list(int = 0),
                algorithm = "port")

nls_climvelo_error <- nls(formula = ShiftKmY ~ ClimVeloTKmY_spp + int,
                    data = data, 
                    start = list(int = 0),
                    algorithm = "port")

nls_limrate_error <- nls(formula = ShiftKmY ~ LimitingRate + int,
                   data = data, 
                   start = list(int = 0),
                   algorithm = "port")

## fit models with estimated slope parameter 
lm_disp_error <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
               data = data, 
               start = list(int = 0, m = 1),
               algorithm = "port")

lm_climvelo_error <- nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                   data = data, 
                   start = list(int = 0, m = 1),
                   algorithm = "port")

lm_limrate_error <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                  data = data, 
                  start = list(int = 0, m = 1),
                  algorithm = "port")

## trait models
lm_dr_error <- nls(formula = ShiftKmY ~ m*log(DispersalPotentialKmY) + int,
             data = data, 
             start = list(int = 0, m = 1),
             algorithm = "port")
lm_bs_error <- nls(formula = ShiftKmY ~ m*log(BodySize) + int,
             data = data, 
             start = list(int = 0, m = 1),
             algorithm = "port")
lm_rs_error <- nls(formula = ShiftKmY ~ m*log(Area_m2_range) + int,
             data = data, 
             start = list(int = 0, m = 0),
             algorithm = "port")

## save models 
error_models <- list(nls_disp_error, nls_climvelo_error, nls_limrate_error,
                         lm_disp_error, lm_climvelo_error, lm_limrate_error,
                         lm_dr_error, lm_bs_error, lm_rs_error)
names(error_models) <- c("nls_disp_error", "nls_climvelo_error", "nls_limrate_error",
                             "lm_disp_error", "lm_climvelo_error", "lm_limrate_error",
                             "lm_dr_error", "lm_bs_error", "lm_rs_error")

saveRDS(error_models, "data-processed/modelfits_error.rds")

## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                  by = 0.001))

df_disp$pred_nls <- predict(nls_disp_error, se.fit = FALSE, newdata = df_disp)
df_disp$pred_lm <- predict(lm_disp_error, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                 by = 0.001))

df_climvelo$pred_nls <- predict(nls_climvelo_error, se.fit = FALSE, newdata = df_climvelo)
df_climvelo$pred_lm <- predict(lm_climvelo_error, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))

df_limrate$pred_nls <- predict(nls_limrate_error, se.fit = FALSE, newdata = df_limrate)
df_limrate$pred_lm <- predict(lm_limrate_error, se.fit = FALSE, newdata = df_limrate)

df_dr <- data.frame(DispersalPotentialKmY = seq(0.0001, max(data$DispersalPotentialKmY),
                                                by = 0.001))
df_bs <- data.frame(BodySize = seq(min(data$BodySize, na.rm = TRUE), max(data$BodySize, na.rm = TRUE),
                                   by = 0.001))
df_rs <- data.frame(Area_m2_range = seq(min(data$Area_m2_range, na.rm = TRUE), 
                                        max(data$Area_m2_range, na.rm = TRUE),
                                        by = 100000000))

df_dr$pred <- predict(lm_dr_error, se.fit = FALSE, newdata = df_dr)
df_bs$pred <- predict(lm_bs_error, se.fit = FALSE, newdata = df_bs)
df_rs$pred <- predict(lm_rs_error, se.fit = FALSE, newdata = df_rs)

## plot
## DISPERSAL
disp_plot = data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE) 

data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_hline(yintercept = 0) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "red", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "purple", inherit.aes = FALSE) 

## CLIMATE
clim_velo_plot = data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_lm), colour = "black", inherit.aes = FALSE)

## LIMITING RATE
limrate_plot =data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE)

plot_grid(disp_plot, clim_velo_plot, limrate_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_expansions-with error.png", 
       width = 12.5, height = 4)

## proxy 
## dispersal rate
dr_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_dr, aes(x = DispersalPotentialKmY, y = pred), colour = "black", inherit.aes = FALSE)

## body size
bs_plot <- data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = BodySize, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Body size (mm)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_bs, aes(x = BodySize, y = pred), colour = "black", inherit.aes = FALSE)

## range size
rs_plot <- data %>%
  ggplot(aes(x = Area_m2_range, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = Area_m2_range, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Realized range size (m2)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_rs, aes(x = Area_m2_range, y = pred), colour = "black", inherit.aes = FALSE)

plot_grid(bs_plot, rs_plot, dr_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_proxy-trait_expansions-with-error.png", 
       width = 12.5, height = 4)

#####################################
##       TESTING CONTRACTIONS      ##
#####################################
## model leading edge and trailing edge range contractions that are in line with climate velocity the same way
## the limiting rate should not explain range contraction as well as it explains range expansion

## read in data 
cont <- read.csv("data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
cont <- filter(cont, Type == "LAT")

## filter to leading and trailing edge shifts that are contractions
cont <- filter(cont, !Param %in% c("O") & Direction == "Contraction")

## filter to shifts in same direction as climate velocity 
cont <- filter(cont, Param == "LE" & ClimVeloTKmY_spp <= 0 | Param == "TE" & ClimVeloTKmY_spp >= 0)

hist(cont$ClimVeloTKmY_spp)
hist(cont$ClimVeloTKmY_spp)

## take absolute value of clim velo and shift rate
cont$ClimVeloTKmY_spp = abs(cont$ClimVeloTKmY_spp)
cont$ShiftKmY = abs(cont$ShiftKmY)

cont <- cont %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_spp)) %>%
  mutate(ClimVeloTKmY_spp = abs(ClimVeloTKmY_spp)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_spp, NA))

## join traits to dd 
cont <- cont %>%
  select(-range_source) %>%
  left_join(., rangesize, by = c("ID", "scientificName" = "binomial")) %>%
  left_join(., select(bodysize, scientificName, class, kingdom, BodySize))

data <- cont

## save data 
write.csv(cont, "data-processed/model-data_contractions.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

## fit models without slope parameter 
nls_disp_cont <- nls(formula = ShiftKmY ~ DispersalPotentialKmY + int,
                      data = data, 
                      start = list(int = 0),
                      algorithm = "port")

nls_climvelo_cont <- nls(formula = ShiftKmY ~ ClimVeloTKmY_spp + int,
                          data = data, 
                          start = list(int = 0),
                          algorithm = "port")

nls_limrate_cont <- nls(formula = ShiftKmY ~ LimitingRate + int,
                         data = data, 
                         start = list(int = 0),
                         algorithm = "port")

## fit models with estimated slope parameter 
lm_disp_cont <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                     data = data, 
                     start = list(int = 0, m = 1),
                     algorithm = "port")

lm_climvelo_cont <- nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                         data = data, 
                         start = list(int = 0, m = 1),
                         algorithm = "port")

lm_limrate_cont <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                        data = data, 
                        start = list(int = 0, m = 1),
                        algorithm = "port")

lm_dr_cont <- nls(formula = ShiftKmY ~ m*log(DispersalPotentialKmY) + int,
                    data = data, 
                    start = list(int = 0, m = 1),
                    algorithm = "port")
lm_bs_cont <- nls(formula = ShiftKmY ~ m*log(BodySize) + int,
                  data = data, 
                  start = list(int = 0, m = 1),
                  algorithm = "port")
lm_rs_cont <- nls(formula = ShiftKmY ~ m*log(Area_m2_range) + int,
                  data = data, 
                  start = list(int = 0, m = 1),
                  algorithm = "port")

## save models
cont_models <- list(nls_disp_cont, nls_climvelo_cont, nls_limrate_cont,
                     lm_disp_cont, lm_climvelo_cont, lm_limrate_cont,
                     lm_dr_cont, lm_bs_cont, lm_rs_cont)
names(cont_models) <- c("nls_disp_cont", "nls_climvelo_cont", "nls_limrate_cont",
                         "lm_disp_cont", "lm_climvelo_cont", "lm_limrate_cont",
                         "lm_dr_cont", "lm_bs_cont", "lm_rs_cont")

saveRDS(cont_models, "data-processed/modelfits_contractions.rds")


## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                  by = 0.001))

df_disp$pred_nls <- predict(nls_disp_cont, se.fit = FALSE, newdata = df_disp)
df_disp$pred_lm <- predict(lm_disp_cont, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                 by = 0.001))

df_climvelo$pred_nls <- predict(nls_climvelo_cont, se.fit = FALSE, newdata = df_climvelo)
df_climvelo$pred_lm <- predict(lm_climvelo_cont, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))

df_limrate$pred_nls <- predict(nls_limrate_cont, se.fit = FALSE, newdata = df_limrate)
df_limrate$pred_lm <- predict(lm_limrate_cont, se.fit = FALSE, newdata = df_limrate)

df_dr <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                by = 0.001))
df_bs <- data.frame(BodySize = seq(min(data$BodySize, na.rm = TRUE), max(data$BodySize, na.rm = TRUE),
                                   by = 0.001))
df_rs <- data.frame(Area_m2_range = seq(min(data$Area_m2_range, na.rm = TRUE), 
                                        max(data$Area_m2_range, na.rm = TRUE),
                                        by = 100000000))

df_dr$pred <- predict(lm_dr_cont, se.fit = FALSE, newdata = df_dr)
df_bs$pred <- predict(lm_bs_cont, se.fit = FALSE, newdata = df_bs)
df_rs$pred <- predict(lm_rs_cont, se.fit = FALSE, newdata = df_rs)

## plot
## DISPERSAL
disp_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE)

## CLIMATE
clim_velo_plot <- data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", fill = "transparent", 
             pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_lm), colour = "black", inherit.aes = FALSE)

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  labs(x = "Minimum of rate of climate change and\npotential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE)


plot_grid(disp_plot, clim_velo_plot, limrate_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_contractions.png", 
       width = 12.5, height = 4)


## proxy 
## dispersal rate
dr_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_dr, aes(x = DispersalPotentialKmY, y = pred), colour = "black", inherit.aes = FALSE)

## body size
bs_plot <- data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = BodySize, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  labs(x = "Body size (mm)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_bs, aes(x = BodySize, y = pred), colour = "black", inherit.aes = FALSE)

## range size
rs_plot <- data %>%
  ggplot(aes(x = Area_m2_range, y = ShiftKmY)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
             fill = "transparent", pch = 1,
             aes(x = Area_m2_range, y = ShiftKmY, shape = group)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0.5)) +
  labs(x = "Realized range size (m2)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_line(data = df_rs, aes(x = Area_m2_range, y = pred), colour = "black", inherit.aes = FALSE)

plot_grid(bs_plot, rs_plot, dr_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_proxy-trait_contractions.png", 
       width = 12.5, height = 4)

###################################
##         RANDOMIZATION         ##
###################################
## randomization test
#########################
## randomize the dependent variable (range shift rate) 1000 times
## fit the model and extract model coefficients 
## compare to observed coefficients 
data <- filter(dd, is_contraction == "NO")

## randomize range shift rate 
coeffs <- c()
for(i in 1:1000) {
  
  ## randomize shifts 
  rando_df <- data.frame(LimitingRate = data$LimitingRate,
                         ClimVeloTKmY_spp = data$ClimVeloTKmY_spp, 
                         DispersalPotentialKmY = data$DispersalPotentialKmY, 
                         rando_shifts = sample(data$ShiftKmY, size = nrow(data), replace = FALSE))
  
  ## fit linear model
  lm1 <- lm(rando_shifts ~ DispersalPotentialKmY,
            data = rando_df)
  lm2 <- lm(rando_shifts ~ ClimVeloTKmY_spp,
            data = rando_df)
  lm3 <- lm(rando_shifts ~ LimitingRate,
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
                             mod_type = c("DispersalPotentialKmY", "ClimVeloTKmY_spp", "LimitingRate"),
                             slope = slopes,
                             intercept = intercepts,
                             r_squared = r_squ,
                             adj_r_squared = adj_r_squ))
  print(i)
}


## extract coefficients from real models and add to database
coeffs <- coeffs %>%
  mutate(real_slope = ifelse(mod_type == "LimitingRate", coefficients(lm_limrate)[2],
                             ifelse(mod_type == "ClimVeloTKmY_spp", coefficients(lm_climvelo)[2],
                                    coefficients(lm_disp)[2])),
         real_intercept = ifelse(mod_type == "LimitingRate", coefficients(lm_limrate)[1],
                                 ifelse(mod_type == "ClimVeloTKmY_spp", coefficients(lm_climvelo)[1],
                                        coefficients(lm_disp)[1])),
         real_r_squ = ifelse(mod_type == "LimitingRate", summary(lm(ShiftKmY ~ predict(lm_limrate), 
                                                                    data = data))$r.squared,
                             ifelse(mod_type == "ClimVeloTKmY_spp", summary(lm(ShiftKmY ~ predict(lm_climvelo),
                                                                                 data = data))$r.squared,
                                    summary(lm(ShiftKmY ~ predict(lm_disp), data = data))$r.squared)))

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



##################################
##       PLOT MODEL COEFFS      ##
##################################
## make dot whisker plots 

## dispersal limitation models
################################
## expansions 
disp_nls = dwplot(nls_disp) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(-44,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

climvelo_nls = dwplot(nls_climvelo) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept"))  +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

limrate_nls = dwplot(nls_limrate) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

disp = dwplot(lm_disp) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)")

climvelo = dwplot(lm_climvelo) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept"))  +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)")

limrate = dwplot(lm_limrate) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)")

plot_grid(disp_nls, climvelo_nls, limrate_nls, 
          disp, climvelo, limrate, ncol = 1, align = "v")

ggsave(path = "figures/model_results", filename = "dwplot_expansions.png", 
       width = 6, height = 7)

dwplot(list(lm_disp, lm_climvelo, lm_limrate)) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  scale_color_discrete(labels = c("Maximum potential dispersal rate (km/y)", 
                                  "Mean rate of climate change (km/y)", 
                                  "Minimum of potential dispersal rate and climate velocity (km/y)")) +
  labs(colour = "Model") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

dwplot(list(nls_disp, nls_climvelo, nls_limrate)) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(-44,3)) +
  scale_color_discrete(labels = c("Maximum potential dispersal rate (km/y)", 
                                  "Mean rate of climate change (km/y)", 
                                  "Minimum of potential dispersal rate and climate velocity (km/y)")) +
  labs(colour = "Model") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

## accounting for error
disp_nls = dwplot(nls_disp_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(-44,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

climvelo_nls = dwplot(nls_climvelo_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept"))  +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

limrate_nls = dwplot(nls_limrate_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

disp = dwplot(lm_disp_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)")

climvelo = dwplot(lm_climvelo_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept"))  +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)")

limrate = dwplot(lm_limrate_error) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_x_continuous(limits = c(0,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)")

plot_grid(disp_nls, climvelo_nls, limrate_nls, 
          disp, climvelo, limrate, ncol = 1, align = "v")

ggsave(path = "figures/model_results", filename = "dwplot_expansions-with-error.png", 
       width = 6, height = 7)

## contractions
disp_nls = dwplot(nls_disp_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(-46,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

climvelo_nls = dwplot(nls_climvelo_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept"))  +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

limrate_nls = dwplot(nls_limrate_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Intercept")) +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)") +
  scale_colour_manual(values = c("skyblue1"))

disp = dwplot(lm_disp_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept")) +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Maximum potential dispersal rate (km/y)")

climvelo = dwplot(lm_climvelo_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", "Intercept"))  +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Mean rate of climate change (km/y)")

limrate = dwplot(lm_limrate_cont) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_x_continuous(limits = c(-1,3)) +
  labs(subtitle = "Minimum of potential dispersal rate and climate velocity (km/y)")

plot_grid(disp_nls, climvelo_nls, limrate_nls, 
          disp, climvelo, limrate, ncol = 1, align = "v")

ggsave(path = "figures/model_results", filename = "dwplot_contractions.png", 
       width = 6, height = 7)


## proxy traits models
################################
## expansions
dwplot(list(lm_dr, lm_bs, lm_rs)) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_color_discrete(labels = c("Log maximum potential dispersal rate (km/y)", 
                                  "Log body size (mm)", 
                                  "Log realized range size (m2)")) +
  labs(colour = "Model") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

ggsave(path = "figures/model_results", filename = "dwplot_proxy-trait_expansions.png", 
       width = 6, height = 3)

## contractions 
dwplot(list(lm_dr_cont, lm_bs_cont, lm_rs_cont)) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("Slope", 
                              "Intercept")) +
  scale_color_discrete(labels = c("Log maximum potential dispersal rate (km/y)", 
                                  "Log body size (mm)", 
                                  "Log realized range size (m2)")) +
  labs(colour = "Model") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

ggsave(path = "figures/model_results", filename = "dwplot_proxy-trait_contractions.png", 
       width = 6, height = 3)





