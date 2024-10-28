library(AICcmodavg)
library(tidyverse)
library(ggplot2)
library(dotwhisker)
library(nlraa)
library(cowplot)
theme_set(theme_bw())


## functions to extract predictions from nls models 
fun_disp = function(x, mod) {
  predict2_nls(mod, interval = "conf", 
               newdata = data.frame(ClimVeloTKmY_spp = median(data$ClimVeloTKmY_spp), 
                                    DispersalPotentialKmY = x,
                                    m = coefficients(mod)[2], 
                                    m_cv = coefficients(mod)[3], 
                                    m_int = coefficients(mod)[4],
                                    int = coefficients(mod)[1]))
}
fun_bs = function(x, mod) {
  predict2_nls(mod, interval = "conf", 
               newdata = data.frame(ClimVeloTKmY_spp = median(data$ClimVeloTKmY_spp), 
                                    BodySize = x,
                                    m = coefficients(mod)[2], 
                                    m_cv = coefficients(mod)[3], 
                                    m_int = coefficients(mod)[4],
                                    int = coefficients(mod)[1]))
}
fun_rs = function(x, mod) {
  predict2_nls(mod, interval = "conf", 
               newdata = data.frame(ClimVeloTKmY_spp = median(data$ClimVeloTKmY_spp), 
                                    Area_km2_range = x,
                                    m = coefficients(mod)[2], 
                                    m_cv = coefficients(mod)[3], 
                                    m_int = coefficients(mod)[4],
                                    int = coefficients(mod)[1]))
}

## TESTING EXPANSIONS
## filter to leading edge range shifts with positive climate velocity (expect expansion)
## filter to only expansions + contractions that are < 1sd from the mean

## - fit models to range expansions and compare 
##    rs ~ dd
##    rs ~ cv
##    rs ~ min(dd, cv)
# (...are these slopeless models really necessary? could just comment on whether estimated slope = 1)

##    rs ~ slope*dd + int
##    rs ~ slope*cv + int 
##    rs ~ slope*min(dd, cv) + int 

##    rs ~ slope*dd + slope*cv + slope*dd*cv + int

##    rs ~ slope*proxy trait + int

## macroecological proxy traits: body size, range size

##    expect:
##    - under hyp 1 (proxy traits poorly estimate dispersal): dispersal rate model has slope of 1 and is better fit than proxy trait models 
##    - under hyp 2 (some spp not dispersal limited): limiting rate model has slope of 1 & is better fit than climate velocity or dispersal rate alone

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
hist(dd$ShiftKmY)

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
  ungroup() %>%
  mutate(Area_km2_range = Area_m2_range / 1000000) ## convert range size to km2

## join traits to dd 
dd <- dd %>%
  select(-range_source) %>%
  left_join(., rangesize, by = c("ID", "scientificName" = "binomial")) %>%
  left_join(., select(bodysize, scientificName, class, kingdom, BodySize))

###################################
##       TESTING EXPANSIONS      ##
###################################
## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

## get rid of range contractions that are father than 1 sd from the mean shift 
data <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

## save data
write.csv(data, "data-processed/model-data_expansions.csv", row.names = FALSE)

## make version where all observations have traits 
data_traits <- filter(data, !is.na(BodySize), !is.na(Area_km2_range))

## look at distribution of response variable 
hist(data$ShiftKmY)

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

## fit model where dispersal interacts with climate velocity 
lm_disp_int <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
                     m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
                   data = data, 
                   start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                   algorithm = "port")

## trait models
lm_dr <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
               m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
             data = data_traits, 
             start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
             algorithm = "port")
lm_bs <- nls(formula = ShiftKmY ~ m*BodySize + m_cv*ClimVeloTKmY_spp +
               m_int*ClimVeloTKmY_spp*BodySize + int,
             data = data_traits, 
             start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
             algorithm = "port")
lm_rs <- nls(formula = ShiftKmY ~ m*Area_km2_range + m_cv*ClimVeloTKmY_spp +
               m_int*ClimVeloTKmY_spp*Area_km2_range + int,
             data = data_traits, 
             start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
             algorithm = "port")

## save models 
exp_models <- list(nls_disp, nls_climvelo, nls_limrate,
                         lm_disp, lm_climvelo, lm_limrate,
                     lm_disp_int,
                         lm_dr, lm_bs, lm_rs)
names(exp_models) <- c("nls_disp", "nls_climvelo", "nls_limrate",
                             "lm_disp", "lm_climvelo", "lm_limrate",
                         "lm_disp_int",
                             "lm_dr", "lm_bs", "lm_rs")

saveRDS(exp_models, "data-processed/modelfits_expansions.rds")

## predict using each model
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                  by = 0.001))

df_disp$pred_nls <- predict(nls_disp, se.fit = T, newdata = df_disp)
df_disp$pred_lm <- predict(lm_disp, se.fit = FALSE, newdata = df_disp)

df_climvelo <- data.frame(ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                 by = 0.001))

df_climvelo$pred_nls <- predict(nls_climvelo, se.fit = FALSE, newdata = df_climvelo)
df_climvelo$pred_lm <- predict(lm_climvelo, se.fit = FALSE, newdata = df_climvelo)

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))

df_limrate$pred_nls <- predict(nls_limrate, se.fit = FALSE, newdata = df_limrate)
df_limrate$pred_lm <- predict(lm_limrate, se.fit = FALSE, newdata = df_limrate)

## for models with interactions, predict across all climate velocites:
df_disp_int <- expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                       by = 0.1),
                           ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                  by = 0.1))

df_disp_int$pred_lm <- predict(lm_disp_int, se.fit = FALSE, newdata = df_disp_int)


df_dr <- expand.grid(DispersalPotentialKmY = seq(0.0001, max(data$DispersalPotentialKmY),
                                                 by = 0.1),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))
df_bs <- expand.grid(BodySize = seq(min(data$BodySize, na.rm = TRUE), max(data$BodySize, na.rm = TRUE),
                                    by = 0.001),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))

df_rs <- expand.grid(Area_km2_range = seq(min(data$Area_km2_range, na.rm = TRUE), 
                                          max(data$Area_km2_range, na.rm = TRUE),
                                          by = 100000),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))

df_dr$pred <- predict(lm_dr, se.fit = FALSE, newdata = df_dr)
df_bs$pred <- predict(lm_bs, se.fit = FALSE, newdata = df_bs)
df_rs$pred <- predict(lm_rs, se.fit = FALSE, newdata = df_rs)

## bootstrap 95% confidence intervals for linear models
ints_lm_disp <- predict2_nls(lm_disp, interval = "conf")
ints_lm_climvelo <- predict2_nls(lm_climvelo, interval = "conf")
ints_lm_limrate <- predict2_nls(lm_limrate, interval = "conf")

ints_lm_dr <- predict2_nls(lm_dr, interval = "conf")
ints_lm_bs <- predict2_nls(lm_bs, interval = "conf")
ints_lm_rs <- predict2_nls(lm_rs, interval = "conf")

## ones with continuous interactions 
ints_lm_disp_int <- do.call(rbind, lapply(data$DispersalPotentialKmY, fun_disp, mod = lm_disp_int))
ints_lm_dr <- do.call(rbind, lapply(data_traits$DispersalPotentialKmY, fun_disp, mod = lm_dr))
ints_lm_bs <- do.call(rbind, lapply(data_traits$BodySize, fun_bs, mod = lm_bs))
ints_lm_rs <- do.call(rbind, lapply(data_traits$Area_km2_range, fun_rs, mod = lm_rs))


## plot
## DISPERSAL
# disp_plot = data %>%
#   ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
#   geom_point(alpha = 0.7) +
#   geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"), 
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Maximum potential dispersal rate (km/y)",
#        y = "Observed range shift rate (km/y)", 
#        colour = 'Mean\nclimate\nvelocity\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = cbind(data, ints_lm_disp), aes(x = DispersalPotentialKmY, y = fitted(lm_disp),
#                                                     ymin = Q2.5, ymax = Q97.5), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", 
#             inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE) 

# data %>%
#   ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
#   geom_point(alpha = 0.7) +
#   geom_point(data = filter(data, is.na(colour)), inherit.aes = FALSE, colour = "black", 
#              fill = "transparent", pch = 1,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "black", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines")) +
#   scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
#                 labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
#                 limits = c(0.0001, 1400)) +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Potential dispersal rate (km/y)",
#        y = "Observed range shift rate (km/y)", 
#        colour = 'Mean\nclimate\nvelocity\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_hline(yintercept = 0) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "red", inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "purple", inherit.aes = FALSE) 

## DISPERSAL x CLIMATE
disp_int_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp), 
            inherit.aes = FALSE, alpha = 0.1) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
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
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

disp_int_plot_med <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
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
  geom_ribbon(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                               y = Estimate,
                                               ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                                        y = Estimate), 
              colour = "black", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey",
            inherit.aes = FALSE)

## CLIMATE
clim_velo_plot = data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
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
  geom_ribbon(data = cbind(data, ints_lm_climvelo), aes(x = ClimVeloTKmY_spp, y = fitted(lm_climvelo),
                                                    ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_lm), colour = "black", inherit.aes = FALSE)

## LIMITING RATE
limrate_plot = data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
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
  geom_ribbon(data = cbind(data, ints_lm_limrate), aes(x = LimitingRate, y = fitted(lm_limrate),
                                                        ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE)

plot_grid(disp_int_plot_med, clim_velo_plot, limrate_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_expansions.png", 
       width = 12.5, height = 4)

## proxy 
## dispersal rate
dr_plot = data_traits %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_dr), aes(x = DispersalPotentialKmY, y = Estimate,
                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_dr), aes(x = DispersalPotentialKmY, y = Estimate),
            colour = "black", inherit.aes = FALSE)

## body size
bs_plot = data_traits %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Body size (mm)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_bs), aes(x = BodySize, y = Estimate,
                                                                          ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_bs), aes(x = BodySize, y = Estimate),
              colour = "black", inherit.aes = FALSE)

## range size
rs_plot = data_traits %>%
  filter(Area_km2_range != min(.$Area_km2_range)) %>% ## get rid of smallest range size data point for plotting purposes 
  ggplot(aes(x = Area_km2_range, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_x_log10(limits = c(500000, 1e8)) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Realized range size (km2)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_rs), aes(x = Area_km2_range, y = Estimate,
                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_rs), aes(x = Area_km2_range, y = Estimate),
            colour = "black", inherit.aes = FALSE)

plot_grid(bs_plot, rs_plot, dr_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_proxy-trait_expansions.png", 
       width = 12.5, height = 4)


###############################################################################
####      analyzing influence of outliers / points with high leverage      ####
###############################################################################
## Cook's distance can be used to calculate influence on all parameters at once 
exps <- data

## get predictions from full model
preds_full <- predict(lm_limrate, newdata = data.frame(LimitingRate = exps$LimitingRate))

## calculate cook's distance
MSE = mean(residuals(lm_limrate)^2) ## mean squared error
p = 2 ## number of estimated parameters in model

## loop through points, removing one by one, fitting model, and getting prediction for point that is removed

## note: must assign new data frame to old 'lat' object, or else weird error occurs 
## so must refresh data each time 
cooks <- c() 
x = 1
while(x <= nrow(exps)) {
  data = exps
  
  ## remove point x
  data <- data[-x,]
  
  ##  fit model to data without the point:
  mod_sub <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                 data = data, 
                 start = list(int = 0, m = 1),
                 algorithm = "port") 
  
  ## get model predictions
  preds <- predict(mod_sub, newdata = data.frame(LimitingRate = exps$LimitingRate))
  
  ## calculate cook' distance 
  cook <- sum((preds - preds_full)^2)/p*MSE
  
  ## save 
  cooks <- append(cooks, cook)
  
  print(paste0("On obs. number: ", x))
  
  x = x + 1
}

data <- exps

hist(cooks)

## add cook's distance + coefficients to dataset 
exps$cooks_dist = cooks

## plot
exps_influential <- exps %>%
  mutate(is_bigger = cooks_dist >= 1) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = is_bigger, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  # stat_function(colour = "grey", fun = function(x){x},
  #               linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-6, 25), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = "Cook's distance > 1", 
       shape = "") 

## run model without high influence points 
data_sub <- exps[exps$cooks_dist < 1, ]

mod_lowlev <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                  data = data_sub, 
                  start = list(int = 0, m = 1),
                  algorithm = "port") 
summary(mod_lowlev)

## make predictions from model
df_limrate_lowlev <- data.frame(LimitingRate = seq(min(data_sub$LimitingRate), max(data_sub$LimitingRate),
                                                   by = 0.001))

df_limrate_lowlev$pred_lm <- predict(mod_lowlev, se.fit = FALSE, newdata = df_limrate_lowlev)
df_limrate_lowlev$facet = "Without influential points"

ints_lm_limrate_lowlev <- predict2_nls(mod_lowlev, interval = "conf")

df_limrate$facet = "With influential points"
ints_lm_limrate_lowlev$facet = "Without influential points"
ints_lm_limrate$facet = "With influential points"

## plot predictions from model with and without the high leverage points 
facet <- data_sub %>%
  select(-cooks_dist) %>%
  mutate(facet = "Without influential points") %>%
  rbind(., mutate(data, facet = "With influential points")) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 25), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_ribbon(data = cbind(data, ints_lm_limrate), aes(x = LimitingRate, y = fitted(lm_limrate),
                                                       ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_ribbon(data = cbind(data_sub, ints_lm_limrate_lowlev), aes(x = LimitingRate, y = fitted(mod_lowlev),
                                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE) +
  geom_line(data = df_limrate_lowlev, aes(x = LimitingRate, y = pred_lm), 
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~facet)

ggsave(facet, path = "figures/model_results", filename = "model-predictions_influence_expansions.png", 
       width = 8, height = 4)


################################################
##         RANDOMIZATION - EXPANSIONS         ##
################################################
## randomization test
#########################
## randomize the dependent variable (range shift rate) 1000 times
## fit the model and extract model coefficients 
## compare to observed coefficients 
data <- exps

## randomize range shift rate 
coeffs <- c()
for(i in 1:10000) {
  
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
         real_slope_influential = ifelse(mod_type == "LimitingRate", coefficients(mod_lowlev)[2],NA),
         real_intercept = ifelse(mod_type == "LimitingRate", coefficients(lm_limrate)[1],
                                 ifelse(mod_type == "ClimVeloTKmY_spp", coefficients(lm_climvelo)[1],
                                        coefficients(lm_disp)[1])),
         real_intercept_influential = ifelse(mod_type == "LimitingRate", coefficients(mod_lowlev)[1],NA),
         real_r_squ = ifelse(mod_type == "LimitingRate", summary(lm(ShiftKmY ~ predict(lm_limrate), 
                                                                    data = data))$r.squared,
                             ifelse(mod_type == "ClimVeloTKmY_spp", 
                                    summary(lm(ShiftKmY ~ predict(lm_climvelo),
                                                                               data = data))$r.squared,
                                    summary(lm(ShiftKmY ~ predict(lm_disp), data = data))$r.squared)),
         real_r_squ_influential = ifelse(mod_type == "LimitingRate", summary(lm(ShiftKmY ~ predict(mod_lowlev), 
                                                                                data = data_sub))$r.squared,
                                         NA),)

## slope
rando_slope <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = slope)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_slope), colour = "red") +
  geom_vline(aes(xintercept = real_slope_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Estimated slope", y = "Count") 

## intercept
rando_int <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = intercept)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_intercept), colour = "red") +
  geom_vline(aes(xintercept = real_intercept_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Estimated intercept", y = "Count")

## r squared
rando_r2 <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = r_squared)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_r_squ), colour = "red") +
  geom_vline(aes(xintercept = real_r_squ_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Regression coefficient", y = "Count")

plot_grid(rando_slope, rando_int, rando_r2, 
          nrow = 3, align = "v")

ggsave(path = "figures/model_results", filename = "randomization-histograms_expansions.png", 
       width = 8, height = 6)


#####################################
##       TESTING CONTRACTIONS      ##
#####################################
## model leading edge range contractions that are in line with climate velocity in the same way
## the limiting rate should not explain range contraction at the leading edge as well as it explains range expansion 

## read in data 
cont <- read.csv("data-processed/v3_potential-dispersal-rate.csv") %>%
  filter(!is.na(ShiftKmY))

## filter to only latitude 
cont <- filter(cont, Type == "LAT")

## filter to only leading and trailing edge shifts
cont <- filter(cont, Param != "O")

## filter to shifts where climate velocity predicts contraction should occur (i.e., negative climate velocity at leading edge, positive climate velocity at trailing edge)
cont <- filter(cont, Param == "LE" & ClimVeloTKmY_spp <= 0 | Param == "TE" & ClimVeloTKmY_spp >= 0)

## flip sign of trailing edge shifts
cont$ShiftKmY = ifelse(cont$Param == "TE", -1*cont$ShiftKmY, cont$ShiftKmY) 
## and climate velocity 
cont$ClimVeloTKmY_spp = ifelse(cont$Param == "LE", -1*cont$ClimVeloTKmY_spp, cont$ClimVeloTKmY_spp)

## filter to only contractions + expansions within 1sd of the mean 
cont <- filter(cont, ShiftKmY <= (mean(cont$ShiftKmY) + sd(cont$ShiftKmY)))

hist(cont$ClimVeloTKmY_spp)
hist(cont$ShiftKmY)

cont <- cont %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_spp)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_spp, NA)) %>%
  filter(!is.na(ClimVeloTKmY_spp))

## join traits to dd 
cont <- cont %>%
  select(-range_source) %>%
  left_join(., rangesize, by = c("ID", "scientificName" = "binomial")) %>%
  left_join(., select(bodysize, scientificName, class, kingdom, BodySize))

data <- cont

## save data 
write.csv(cont, "data-processed/model-data_contractions.csv", row.names = FALSE)

## make version with no missing traits
data_traits <- filter(data, !is.na(BodySize), !is.na(Area_km2_range))

## look at distribution of response variable 
hist(data$ShiftKmY)

## fit models without slope parameter 
nls_disp_cont <- nls(formula = ShiftKmY ~ -1*DispersalPotentialKmY + int,
                      data = data, 
                      start = list(int = 0),
                      algorithm = "port")

nls_climvelo_cont <- nls(formula = ShiftKmY ~ -1*ClimVeloTKmY_spp + int,
                          data = data, 
                          start = list(int = 0),
                          algorithm = "port")

nls_limrate_cont <- nls(formula = ShiftKmY ~ -1*LimitingRate + int,
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

## fit model where dispersal interacts with climate velocity 
lm_disp_int_cont <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
                           m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
                         data = data, 
                         start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                         algorithm = "port")

## trait models
lm_dr_cont <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
                     m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
                   data = data_traits, 
                   start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                   algorithm = "port")
lm_bs_cont <- nls(formula = ShiftKmY ~ m*BodySize + m_cv*ClimVeloTKmY_spp +
                     m_int*ClimVeloTKmY_spp*BodySize + int,
                   data = data_traits, 
                   start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                   algorithm = "port")
lm_rs_cont <- nls(formula = ShiftKmY ~ m*Area_km2_range + m_cv*ClimVeloTKmY_spp +
                     m_int*ClimVeloTKmY_spp*Area_km2_range + int,
                   data = data_traits, 
                   start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                   algorithm = "port")
## save models
cont_models <- list(nls_disp_cont, nls_climvelo_cont, nls_limrate_cont,
                     lm_disp_cont, lm_climvelo_cont, lm_limrate_cont,
                    lm_disp_int_cont,
                     lm_dr_cont, lm_bs_cont, lm_rs_cont)
names(cont_models) <- c("nls_disp_cont", "nls_climvelo_cont", "nls_limrate_cont",
                         "lm_disp_cont", "lm_climvelo_cont", "lm_limrate_cont",
                        "lm_disp_int_cont",
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

## for models with interactions, predict across all climate velocites:
df_disp_int <- expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                       by = 0.1),
                           ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                  by = 0.1))

df_disp_int$pred_lm <- predict(lm_disp_int_cont, se.fit = FALSE, newdata = df_disp_int)


df_dr <- expand.grid(DispersalPotentialKmY = seq(0.0001, max(data$DispersalPotentialKmY),
                                                 by = 0.1),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))
df_bs <- expand.grid(BodySize = seq(min(data$BodySize, na.rm = TRUE), max(data$BodySize, na.rm = TRUE),
                                    by = 0.001),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))

df_rs <- expand.grid(Area_km2_range = seq(min(data$Area_km2_range, na.rm = TRUE), 
                                          max(data$Area_km2_range, na.rm = TRUE),
                                          by = 100000),
                     ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                            by = 0.1))

df_dr$pred <- predict(lm_dr_cont, se.fit = FALSE, newdata = df_dr)
df_bs$pred <- predict(lm_bs_cont, se.fit = FALSE, newdata = df_bs)
df_rs$pred <- predict(lm_rs_cont, se.fit = FALSE, newdata = df_rs)

## bootstrap 95% confidence intervals for linear models
ints_lm_disp <- predict2_nls(lm_disp_cont, interval = "conf")
ints_lm_climvelo <- predict2_nls(lm_climvelo_cont, interval = "conf")
ints_lm_limrate <- predict2_nls(lm_limrate_cont, interval = "conf")

ints_lm_dr <- predict2_nls(lm_dr_cont, interval = "conf")
ints_lm_bs <- predict2_nls(lm_bs_cont, interval = "conf")
ints_lm_rs <- predict2_nls(lm_rs_cont, interval = "conf")

## ones with continuous interactions 
ints_lm_disp_int <- do.call(rbind, lapply(data$DispersalPotentialKmY, fun_disp, mod = lm_disp_int_cont))
ints_lm_dr <- do.call(rbind, lapply(data_traits$DispersalPotentialKmY, fun_disp, mod = lm_dr_cont))
ints_lm_bs <- do.call(rbind, lapply(data_traits$BodySize, fun_bs, mod = lm_bs_cont))
ints_lm_rs <- do.call(rbind, lapply(data_traits$Area_km2_range, fun_rs, mod = lm_rs_cont))

## plot
## DISPERSAL
# disp_plot <- data %>%
#   ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
#   geom_point(alpha = 0.7) +
#   geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"), 
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-2, 20), expand = c(0,0.5)) +
#   labs(x = "Maximum potential dispersal rate (km/y)",
#        y = "Observed range shift rate (km/y)", 
#        colour = 'Mean\nclimate\nvelocity\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = cbind(data, ints_lm_disp), aes(x = DispersalPotentialKmY, y = fitted(lm_disp_cont),
#                                                   ymin = Q2.5, ymax = Q97.5), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm), colour = "black", inherit.aes = FALSE)

## DISPERSAL x CLIMATE
disp_int_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp), 
            inherit.aes = FALSE, alpha = 0.1) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

disp_int_plot_med <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                                        y = Estimate,
                                                        ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                                      y = Estimate), 
            colour = "black", inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_nls), colour = "grey",
            inherit.aes = FALSE)

## CLIMATE
clim_velo_plot <- data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_ribbon(data = cbind(data, ints_lm_climvelo), aes(x = ClimVeloTKmY_spp, y = fitted(lm_climvelo_cont),
                                                    ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_climvelo, aes(x = ClimVeloTKmY_spp, y = pred_lm), colour = "black", inherit.aes = FALSE)

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_ribbon(data = cbind(data, ints_lm_limrate), aes(x = LimitingRate, y = fitted(lm_limrate_cont),
                                                        ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_nls), colour = "grey", inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE)


plot_grid(disp_int_plot_med, clim_velo_plot, limrate_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_contractions.png", 
       width = 12.5, height = 4)


## proxy 
## dispersal rate
dr_plot = data_traits %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_dr), aes(x = DispersalPotentialKmY, y = Estimate,
                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_dr), aes(x = DispersalPotentialKmY, y = Estimate),
            colour = "black", inherit.aes = FALSE)

## body size
bs_plot = data_traits %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Body size (mm)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_bs), aes(x = BodySize, y = Estimate,
                                                                          ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_bs), aes(x = BodySize, y = Estimate),
            colour = "black", inherit.aes = FALSE)

## range size
rs_plot = data_traits %>%
  ggplot(aes(x = Area_km2_range, y = ShiftKmY, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Realized range size (km2)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  geom_ribbon(data = cbind(data_traits, ints_lm_rs), 
              aes(x = Area_km2_range, y = Estimate,
                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data_traits, ints_lm_rs), aes(x = Area_km2_range, y = Estimate),
            colour = "black", inherit.aes = FALSE)

plot_grid(bs_plot, rs_plot, dr_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "model-predictions_proxy-trait_contractions.png", 
       width = 12.5, height = 4)

###############################################################################
####      analyzing influence of outliers / points with high leverage      ####
###############################################################################
## Cook's distance can be used to calculate influence on all parameters at once 

## get predictions from full model
preds_full <- predict(lm_limrate_cont, newdata = data.frame(LimitingRate = cont$LimitingRate))

## calculate cook's distance
MSE = mean(residuals(lm_limrate_cont)^2)
p = 2

## loop through points, removing one by one, fitting model, and getting prediction for point that is removed

## note: must assign new data frame to old 'lat' object, or else weird error occurs 
## so must refresh data each time 
cooks <- c()
x = 1
while(x <= nrow(cont)) {
  data = cont

  ## remove point x
  data <- data[-x,]
  
  ##  fit model to data without the point:
  mod_sub <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                         data = data, 
                         start = list(int = 0, m = 1),
                         algorithm = "port") 
  
  ## get model predictions
  preds <- predict(mod_sub, newdata = data.frame(LimitingRate = cont$LimitingRate))
  
  ## calculate cook' distance 
  cook <- sum((preds - preds_full)^2)/p*MSE
    
  ## save 
  cooks <- append(cooks, cook)
  
  print(paste0("On obs. number: ", x))
  
  x = x + 1
}
data <- cont

hist(cooks)

## add cook's distance to dataset 
cont$cooks_dist = cooks

## Cook's distance of > 1 = high influence 
## plot
cont_influential <- cont %>%
  mutate(is_bigger = cooks_dist >= 1) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = is_bigger, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = "Cook's distance > 1", 
       shape = "") 


## run model without high leverage points 
data_sub <- cont[cont$cooks_dist < 1,]

mod_lowlev <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
               data = data_sub, 
               start = list(int = 0, m = 1),
               algorithm = "port") 
summary(mod_lowlev)

## make predictions from model
df_limrate_lowlev <- data.frame(LimitingRate = seq(min(data_sub$LimitingRate), max(data_sub$LimitingRate),
                                            by = 0.001))

df_limrate_lowlev$pred_lm <- predict(mod_lowlev, se.fit = FALSE, newdata = df_limrate_lowlev)
df_limrate_lowlev$facet = "Without influential points"

ints_lm_limrate_lowlev <- predict2_nls(mod_lowlev, interval = "conf")

df_limrate$facet = "With influential points"
ints_lm_limrate_lowlev$facet = "Without influential points"
ints_lm_limrate$facet = "With influential points"

## plot predictions from model with and without the high leverage points 
facet <- data_sub %>%
  select(-cooks_dist) %>%
  mutate(facet = "Without influential points") %>%
  rbind(., mutate(data, facet = "With influential points")) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){-x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-18, 6), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_ribbon(data = cbind(data, ints_lm_limrate), aes(x = LimitingRate, y = fitted(lm_limrate_cont),
                                                       ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_ribbon(data = cbind(data_sub, ints_lm_limrate_lowlev), aes(x = LimitingRate, y = fitted(mod_lowlev),
                                                                  ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE) +
  geom_line(data = df_limrate_lowlev, aes(x = LimitingRate, y = pred_lm), 
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~facet)

## save 
ggsave(facet, path = "figures/model_results", filename = "model-predictions_influence_contractions.png", 
       width = 8, height = 4)


##################################################
##         RANDOMIZATION - CONTRACTIONS         ##
##################################################
## randomization test
#########################
## randomize the dependent variable (range shift rate) 1000 times
## fit the model and extract model coefficients 
## compare to observed coefficients


## randomize range shift rate 
coeffs <- c()
for(i in 1:10000) {
  
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
  mutate(real_slope = ifelse(mod_type == "LimitingRate", coefficients(lm_limrate_cont)[2],
                             ifelse(mod_type == "ClimVeloTKmY_spp", coefficients(lm_climvelo_cont)[2],
                                    coefficients(lm_disp_cont)[2])),
         real_slope_influential =  ifelse(mod_type == "LimitingRate", coefficients(mod_lowlev)[2], NA),
         real_intercept = ifelse(mod_type == "LimitingRate", coefficients(lm_limrate_cont)[1],
                                 ifelse(mod_type == "ClimVeloTKmY_spp", coefficients(lm_climvelo_cont)[1],
                                        coefficients(lm_disp_cont)[1])),
         real_intercept_influential =  ifelse(mod_type == "LimitingRate", coefficients(mod_lowlev)[1], NA),
         real_r_squ = ifelse(mod_type == "LimitingRate", summary(lm(ShiftKmY ~ predict(lm_limrate_cont), 
                                                                    data = data))$r.squared,
                             ifelse(mod_type == "ClimVeloTKmY_spp",
                                    summary(lm(ShiftKmY ~ predict(lm_climvelo_cont),
                                                                               data = data))$r.squared,
                                    summary(lm(ShiftKmY ~ predict(lm_disp_cont), data = data))$r.squared)),
         real_r_squ_influential = ifelse(mod_type == "LimitingRate", 
                                         summary(lm(ShiftKmY ~ predict(mod_lowlev), 
                                                    data = data_sub))$r.squared,
                                         NA))

## slope
rando_slope <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = slope)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_slope), colour = "red") +
  geom_vline(aes(xintercept = real_slope_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Estimated slope", y = "Count")

## intercept
rando_int <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = intercept)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_intercept), colour = "red") +
  geom_vline(aes(xintercept = real_intercept_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Estimated intercept", y = "Count")

## r squared
rando_r2 <- coeffs %>%
  mutate(mod_type = ifelse(mod_type == "ClimVeloTKmY_spp", 
                           "Mean rate of climate change (km/y)", 
                           ifelse(mod_type == "DispersalPotentialKmY",
                                  "Maximum potential dispersal rate (km/y)",
                                  "Minimum of potential dispersal rate\nand climate velocity (km/y)"))) %>%
  ggplot(aes(x = r_squared)) +
  geom_histogram() +
  geom_vline(aes(xintercept = real_r_squ), colour = "red") +
  geom_vline(aes(xintercept = real_r_squ_influential), colour = "red", linetype = "longdash") +
  facet_wrap(~mod_type) +
  labs(x = "Regression coefficient", y = "Count")

plot_grid(rando_slope, rando_int, rando_r2, 
          nrow = 3, align = "v")

ggsave(path = "figures/model_results", filename = "randomization-histograms_contractions.png", 
       width = 8, height = 6)


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
  scale_x_continuous(limits = c(-1,3)) +
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
                                  "Log realized range size (km2)")) +
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
                                  "Log realized range size (km2)")) +
  labs(colour = "Model") +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

ggsave(path = "figures/model_results", filename = "dwplot_proxy-trait_contractions.png", 
       width = 6, height = 3)
