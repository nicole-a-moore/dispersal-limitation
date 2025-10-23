## fit main models with maximum dispersal distance
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringr)
library(gt)
library(AICcmodavg)
library(nlme)
library(ggallin)
theme_set(theme_bw())

get_formula <- function(lm) {
  return(as.character(lm$call)[2])
}


#############################
##       PREPARE DATA      ##
#############################
## read in data 
dd <- read.csv("data-processed/v3_with-cv.csv")

## see how many leading edge shifts with negative clim velocity
dd %>%
  filter(Param == "LE") %>%
  filter(ClimVeloKmY_RelScale < 0)

## filter to leading edge shifts 
dd = filter(dd, !Param %in% c("O", "TE"))

## get info on % of observations per group
n_sp = length(unique(dd$sp_name_checked)) #279
round(length(unique(dd$sp_name_checked[which(dd$group == "Aves")]))/n_sp*100, digits = 2) #45.52%
round(length(unique(dd$sp_name_checked[which(dd$group == "Plants")]))/n_sp*100, digits = 2) #51.25%
round(length(unique(dd$sp_name_checked[which(dd$group == "Squamata")]))/n_sp*100, digits = 2) #1.43%
round(length(unique(dd$sp_name_checked[which(dd$group == "Mammalia")]))/n_sp*100, digits = 2) # 0.36%
round(length(unique(dd$sp_name_checked[which(dd$group == "Amphibia")]))/n_sp*100, digits = 2) # 1.43%

## get number of shift estimates
nrow(dd)

## get stats about the mean, median and max number of shifts per species:
tally = dd %>%
  group_by(sp_name_checked) %>%
  tally()
max(tally$n) # 5
mean(tally$n) # 1.71
median(tally$n) # 1

## filter to shifts with positive mean climate velocity 
dd <- filter(dd, ClimVeloKmY_RelScale >= 0)

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloKmY_RelScale,
                               DispersalPotentialKmY,
                               ClimVeloKmY_RelScale)) %>%
  mutate(LimitingRate_median = ifelse(MedianDispersalPotentialKmY <= ClimVeloKmY_RelScale,
                               MedianDispersalPotentialKmY,
                               ClimVeloKmY_RelScale)) %>%
  mutate(LimitingRate_q3 = ifelse(DispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                                  DispersalPotentialKmY,
                                  q3ClimVeloKmY_RelScale)) %>%
  mutate(LimitingRate_q3_median = ifelse(MedianDispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                                  MedianDispersalPotentialKmY,
                                  q3ClimVeloKmY_RelScale)) %>%
  mutate(ClimVeloKmY_RelScale = abs(ClimVeloKmY_RelScale)) %>%
  mutate(q3ClimVeloKmY_RelScale = abs(q3ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate"),
         what_is_limiting_q3 = ifelse(DispersalPotentialKmY == LimitingRate_q3, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY_RelScale, NA),
         colour_q3 = ifelse(what_is_limiting_q3 == "Climate", q3ClimVeloKmY_RelScale, NA)) %>%
  mutate(group = ifelse(group == "Aves", "Bird", 
                        ifelse(group == "Mammalia", "Mammal", 
                               ifelse(group == "Plants", "Plant", 
                                      ifelse(group == "Squamata", "Squamate",
                                             NA)))))

## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

## get # of contractions 
length(which(dd$ShiftKmY < 0)) # 177

## get rid of range contractions that are father than 1 sd from the mean shift 
contractors <- filter(dd, Rate < (mean(dd$ShiftKmY) - sd(dd$ShiftKmY))) %>%
  arrange(ShiftKmY) %>%
  select(ShiftKmY, scientificName_checked) %>%
  distinct()

nrow(contractors) ## 29 extreme range contractions

## write out:
write.csv(contractors, "data-processed/intermediate-files/extreme-contractors.csv", row.names = F)

data <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

## get rid of non-birds and non-plants
data = filter(data, group %in% c("Bird", "Plant"))

## get stats about the mean, median and max number of shifts per species:
tally = data %>%
  group_by(sp_name_checked) %>%
  tally()
max(tally$n) # 5
median(tally$n) # 1

## save data
write.csv(data, "data-processed/model_data/model-data_main.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

## get % of species with dispersal > clim velo
round(length(which(data$DispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 62 % max disp. mean cv.
round(length(which(data$DispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 55 % max disp. q3 cv
## birds:
birds = filter(data, group == "Bird")
round(length(which(birds$DispersalPotentialKmY > birds$ClimVeloKmY_RelScale))/nrow(birds)*100, digits = 0) ## 97 % max disp. mean cv.
## plants
plants = filter(data, group == "Plant")
round(length(which(plants$DispersalPotentialKmY > plants$ClimVeloKmY_RelScale))/nrow(plants)*100, digits = 0) ## 29% max disp. mean cv.
## median
round(length(which(data$MedianDispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 51% med disp. mean cv.
round(length(which(data$MedianDispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 44% med disp. q3 cv
## plants
round(length(which(data$MedianDispersalPotentialKmY > data$ClimVeloKmY_RelScale & data$group == "Plant"))/nrow(data)*100, digits = 0) ## 10% med disp.
round(length(which(data$MedianDispersalPotentialKmY > data$q3ClimVeloKmY_RelScale & data$group == "Plant"))/nrow(data)*100, digits = 0) ## 6% med disp.

## get min and max 
min(data$DispersalPotentialKmY)
max(data$DispersalPotentialKmY)

## function to sqrt transform with negative numbers:
ssqrt_trans <- scales::trans_new(name      = 'signed square root',
                                 transform = function(x) sign(x) * sqrt(abs(x)),
                                 inverse   = function(y) sign(y) * y^2,
                                 domain    = c(-Inf,Inf))
## density plot of all 3 rates
density = data %>%
  gather(key = "type", value = "measure", c(ClimVeloKmY_RelScale, DispersalPotentialKmY, ShiftKmY)) %>%
  mutate(type = factor(type, levels = c("DispersalPotentialKmY", "ClimVeloKmY_RelScale", "ShiftKmY"), 
                       ordered = T)) %>%
  ggplot(aes(x = type, y = measure, colour = group)) + 
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("Potential\ndispersal rate", 
                              "Local velocity of\nisotherm shifts", 
                              "Local range\nexpansion rate")) +
  geom_point(position = position_jitterdodge(seed = 120), 
             size = 0.1, alpha = 0.4) +
  geom_violin(fill = "transparent", scale = "width") +
  labs(x = "", y = "Rate (km/yr)", colour = "") +
  scale_y_continuous(trans = ssqrt_trans, breaks = c(0, 200, 400, 600, 800, 1000, 1200)) +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#e49e00", "#A2B06D"))

## save as new figure 3
ggsave(density, path = "figures",  filename = "figure3_distribution-of-rates.png",
       width = 7, height = 3)

## make a pie chart
library(ggpattern)
library(ggrepel)
pie = data %>%
  mutate(is_faster = ifelse(ClimVeloKmY_RelScale > DispersalPotentialKmY, 
                            "Potential dispersal rate *slower* than<br />velocity of isotherm shifts", 
                            "Potential dispersal rate *faster* than<br />velocity of isotherm shifts")) %>%
  mutate(is_faster = factor(is_faster, ordered = TRUE, 
                            levels = c("Potential dispersal rate *slower* than<br />velocity of isotherm shifts", 
                                       "Potential dispersal rate *faster* than<br />velocity of isotherm shifts"))) %>%
  count(group, is_faster) %>% 
  mutate(prop = n / sum(n)) %>%
  mutate(group = factor(group)) %>%
  arrange(desc(prop)) %>%
  ggplot(aes(x = "", y = prop, fill = group), colour = "black") +
  labs(y = "% of total range shift estimates", fill = "", pattern = "") +
  geom_bar(stat = "identity",  alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_col_pattern(aes(fill = group, pattern = is_faster), colour = "black", size = 0.5,
                   pattern_colour = "black",
                   pattern_fill = "black",
                   pattern_size = 0.1, 
                   pattern_density = 0.05,  alpha = 0.3) +
  theme_void() +
  scale_pattern_discrete(choices = c("none", "stripe")) +
  theme(legend.text = ggtext::element_markdown()) +
  guides(fill = "none", 
         pattern = guide_legend(override.aes = list(fill = "transparent", colour = "black"))) +
  geom_label_repel(aes(label = scales::percent(prop, accuracy = 1)), y = c(0.75, 0.3, 0.09, 0.99),
                   colour = "black", nudge_x = 0.5, min.segment.length = 10) +
  scale_fill_manual(values = c("white", "white"))

ggsave(pie, path = "figures",  filename = "figure3_pie-chart.png",
       width = 4, height = 3)

###################################
##          FIT MODELS           ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY, 
               random = ~ 1|sp_name_checked,
               data = data, 
               method = "ML")
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY, 
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale, 
             random = ~ 1|sp_name_checked,
             data = data,
             method = "ML")
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~ 1|sp_name_checked,
                 data = data,
                 method = "ML")

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "ML")
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY + ClimVeloKmY_RelScale,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "ML")
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "ML")

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale,
                    random = ~ 1|sp_name_checked,
                    data = data,
                    method = "ML")
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                       random = ~ 1|sp_name_checked,
                       data = data,
                       method = "ML")
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale,
                           random = ~ 1|sp_name_checked,
                           data = data,
                           method = "ML")
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                              random = ~ 1|sp_name_checked,
                              data = data,
                              method = "ML")

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate <- lme(ShiftKmY ~ LimitingRate,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "ML")
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate_median,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "ML")
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3_median,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "ML")

## save models 
main_models_ML <- list(lme_disp, lme_cv, lme_cv_q3, lme_disp_cv, lme_disp_cv_q3, 
                       lme_disp_int, lme_disp_int_q3, lme_limrate, lme_limrate_q3,
                       lme_disp_median, lme_disp_cv_median, lme_disp_cv_q3_median, 
                       lme_disp_int_median, lme_disp_int_q3_median, lme_limrate_median, lme_limrate_q3_median)
names(main_models_ML) <- c("lme_disp","lme_cv", "lme_cv_q3", "lme_disp_cv", "lme_disp_cv_q3", 
                           "lme_disp_int", "lme_disp_int_q3", "lme_limrate", "lme_limrate_q3",
                           "lme_disp_median", "lme_disp_cv_median", 
                           "lme_disp_cv_q3_median", "lme_disp_int_median", "lme_disp_int_q3_median", 
                           "lme_limrate_median", "lme_limrate_q3_median")

saveRDS(main_models_ML, "data-processed/model_fits/modelfits_main_lme_ML.rds")

## get aic and rank
aic_main <- aictab(cand.set = main_models_ML, modnames = names(main_models_ML)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames) %>%
  arrange(AICc) %>%
  mutate(Rank = 1:nrow(.))

## refit using REML 
## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY, 
               random = ~ 1|sp_name_checked,
               data = data, 
               method = "REML")
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY, 
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale, 
             random = ~ 1|sp_name_checked,
             data = data,
             method = "REML")
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~ 1|sp_name_checked,
                 data = data,
                 method = "REML")

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "REML")
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY + ClimVeloKmY_RelScale,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "REML")
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "REML")

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale,
                    random = ~ 1|sp_name_checked,
                    data = data,
                    method = "REML")
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                       random = ~ 1|sp_name_checked,
                       data = data,
                       method = "REML")
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale,
                           random = ~ 1|sp_name_checked,
                           data = data,
                           method = "REML")
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                              random = ~ 1|sp_name_checked,
                              data = data,
                              method = "REML")

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate <- lme(ShiftKmY ~ LimitingRate,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "REML")
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate_median,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "REML")
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3_median,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "REML")

## plot residuals 
plot(lme_disp) 
plot(lme_cv)
plot(lme_cv_q3)
plot(lme_disp_cv)
plot(lme_disp_cv_q3)
plot(lme_disp_int)
plot(lme_disp_int_q3)
plot(lme_limrate)
plot(lme_limrate_q3)
plot(lme_disp_median) 
plot(lme_disp_cv_median)
plot(lme_disp_cv_q3_median)
plot(lme_disp_int_median)
plot(lme_disp_int_q3_median)
plot(lme_limrate_median)
plot(lme_limrate_q3_median)

hist(residuals(lme_disp))
hist(residuals(lme_cv))
hist(residuals(lme_cv_q3))
hist(residuals(lme_disp_cv))
hist(residuals(lme_disp_cv_q3))
hist(residuals(lme_disp_int))
hist(residuals(lme_disp_int_q3))
hist(residuals(lme_limrate))
hist(residuals(lme_limrate_q3))
hist(residuals(lme_disp_median))
hist(residuals(lme_disp_cv_median))
hist(residuals(lme_disp_cv_q3_median))
hist(residuals(lme_disp_int_median))
hist(residuals(lme_disp_int_q3_median))
hist(residuals(lme_limrate_median))
hist(residuals(lme_limrate_q3_median))

## model summary 
summary(lme_disp) 
summary(lme_cv)
summary(lme_cv_q3)
summary(lme_disp_cv)
summary(lme_disp_cv_q3)
summary(lme_disp_int)
summary(lme_disp_int_q3)
summary(lme_limrate)
summary(lme_limrate_q3)
summary(lme_disp_median) 
summary(lme_disp_cv_median)
summary(lme_disp_cv_q3_median)
summary(lme_disp_int_median)
summary(lme_disp_int_q3_median)
summary(lme_limrate_median)
summary(lme_limrate_q3_median)

## save models - REML
main_models <- list(lme_disp, lme_cv, lme_cv_q3, lme_disp_cv, lme_disp_cv_q3, 
                       lme_disp_int, lme_disp_int_q3, lme_limrate, lme_limrate_q3,
                       lme_disp_median, lme_disp_cv_median, lme_disp_cv_q3_median, 
                       lme_disp_int_median, lme_disp_int_q3_median, lme_limrate_median, lme_limrate_q3_median)
names(main_models) <- c("lme_disp","lme_cv", "lme_cv_q3", "lme_disp_cv", "lme_disp_cv_q3", 
                           "lme_disp_int", "lme_disp_int_q3", "lme_limrate", "lme_limrate_q3",
                           "lme_disp_median", "lme_disp_cv_median", 
                           "lme_disp_cv_q3_median", "lme_disp_int_median", "lme_disp_int_q3_median", 
                           "lme_limrate_median", "lme_limrate_q3_median")

saveRDS(main_models, "data-processed/model_fits/modelfits_main_lme_REML.rds")



#############################
##       MAKE A TABLE      ##
#############################
## get model equations
formulas <- unlist(lapply(main_models, get_formula))

## get types
disp_type = ifelse(names(main_models) %in% c("lme_disp", "lme_disp_cv", "lme_disp_cv_q3", 
                                             "lme_disp_int", "lme_disp_int_q3",
                                             "lme_limrate", "lme_limrate_q3"), 
                   "Maximum",
                   ifelse(names(main_models) %in% c("lme_disp_median", "lme_disp_cv_median", "lme_disp_cv_q3_median",
                                                    "lme_disp_int_median", "lme_disp_int_q3_median", "lme_cv_q3_median",
                                                    "lme_limrate_q3_median", "lme_limrate_median", "lme_cv_median"),
                          "Median",
                          "-"))
cv_type = ifelse(names(main_models) %in% c("lme_cv", "lme_disp_cv",  "lme_disp_int", "lme_cv_median",
                                           "lme_disp_cv_median", "lme_disp_int_median",
                                           "lme_limrate", "lme_limrate_median"), 
                 "Mean",
                 ifelse(names(main_models) %in% c("lme_cv_q3", "lme_disp_cv_q3",  "lme_disp_int_q3", "lme_cv_q3_median",
                                                  "lme_disp_cv_q3_median", "lme_disp_int_q3_median",
                                                  "lme_limrate_q3", "lme_limrate_q3_median"),
                        "3rd quartile",
                        "-"))
model_type = ifelse(names(main_models) %in% c("lme_cv", "lme_cv_q3"), 
                    "Velocity of isotherm shift",
                    ifelse(names(main_models) %in% c("lme_limrate_q3_median", "lme_limrate_q3", 
                                                     "lme_limrate", "lme_limrate_median"),
                           "Minimum rate",
                           ifelse(names(main_models) %in% c("lme_disp", "lme_disp_q3", "lme_disp_q3_median", "lme_disp_median"),
                                  "Potential dispersal rate", 
                                  ifelse(names(main_models) %in% c("lme_disp_cv", "lme_disp_cv_q3", 
                                                                   "lme_disp_cv_q3_median", "lme_disp_cv_median"),
                                         "Potential dispersal rate and velocity of isotherm shift (additive)",
                                         ifelse(names(main_models) %in% c("lme_disp_int", "lme_disp_int_q3", 
                                                                          "lme_disp_int_median", "lme_disp_int_q3_median"),
                                                "Potential dispersal rate and velocity of isotherm shift (interactive)",
                                                "-")))))


## get cond and marg r squared 
cond_rsq <- unlist(lapply(main_models, FUN = function(lme) {
 cond_r2 = performance::r2_nakagawa(lme)[1]
return(cond_r2)
}))

marg_rsq <- unlist(lapply(main_models, FUN = function(lme) {
  marg_r2 = performance::r2_nakagawa(lme)[2]
  return(marg_r2)
}))
            
## get n 
n <- c(rep(nrow(data), length(main_models)))

## get aic and rank
aic_main <- aictab(cand.set = main_models_ML, modnames = names(main_models_ML)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames) %>%
  arrange(AICc) %>%
  mutate(Rank = 1:nrow(.))

## get coefs, join to r squared + se and aic table
coefs <- lapply(main_models, FUN = function(x) {summary(x)$tTable})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(main_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared_marg = marg_rsq[i]
  cur$r_squared_cond = cond_rsq[i]
  cur$n = n[i]
  cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
  cur$cv_type = cv_type[i]
  cur$disp_type = disp_type[i]
  cur$model_type = model_type[i]
  coefs[[i]] <- cur
}
coefs <- coefs %>%
  bind_rows() %>%
  left_join(., aic_main) 

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## reorder/rename columns 
coefs <- coefs %>%
  rename("Conditional R2" = r_squared_cond, "Marginal R2" = r_squared_marg, "ΔAICc" = Delta_AICc,
         "AIC weight" = AICcWt, "Estimate" = Value) %>%
  select(-Model) %>%
  rename("Model" = model_type) %>%
  select(Model, cv_type, disp_type, Formula, Parameter, Estimate, `Std.Error`,
         `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, LL, AICc, everything())  %>%
  arrange(`ΔAICc`) %>% 
  mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", 
                            ifelse(str_detect(Parameter, "\\:"), 
                                   "Potential dispersal rate:velocity of isotherm shift", 
                                   ifelse(str_detect(Parameter, "ClimVeloKmY"), "Velocity of isotherm shift",
                                          ifelse(str_detect(Parameter, "DispersalPotential"), 
                                                 "Potential dispersal rate", 
                                                 "Minimum rate"))))) %>%
  select(-Formula, -Cum.Wt)

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  rename("Local velocity of isotherm shifts" = cv_type,
         "Potential dispersal rate" = disp_type)

## make gt table
table_main <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - all observations"
  ) 

gtsave(table_main, path = "figures", filename = "tableS1_all-models_lme.png")
gtsave(table_main, path = "figures", filename = "tableS1_all-models_lme.docx")


## make table summarizing results by model 'type' 
## calculate cumulative weight and mean rank per model 'type'
table1 = coefs %>%
  select(Model, `Potential dispersal rate`, `Local velocity of isotherm shifts`, `AIC weight`) %>%
  distinct() %>%
  mutate(ModelRank = 1:nrow(.)) %>%
  group_by(Model) %>%
  mutate(LowestRank = min(ModelRank),
         MeanRank = mean(ModelRank), 
         CumulativeWeight = sum(as.numeric(`AIC weight`))) %>%
  select(Model, LowestRank, MeanRank, CumulativeWeight) %>%
  distinct() %>%
  ungroup() %>%
  rename("Lowest model rank" = LowestRank, "Mean model rank" = MeanRank, "Cumulative weight" = CumulativeWeight,
         "Model type" = Model)

table1 = table1 %>% gt()

## save table 
gtsave(table1, path = "figures", filename = "table1_summary.docx")


#################################
##       PLOT PREDICTIONS      ##
#################################
## plot predictions 
df_limrate_q3 <- data.frame(LimitingRate_q3 = seq(min(data$LimitingRate_q3), max(data$LimitingRate_q3),
                                                  by = 0.001))
pred =  predictSE(lme_limrate_q3, se.fit = TRUE, newdata = df_limrate_q3, level = 0)
df_limrate_q3$pred_lme <- pred$fit
df_limrate_q3$pred_lme_se <-  pred$se.fit

## new figure 4:
fig4a = data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Local velocity of\nisotherm\nshifts (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, unit = "cm"))

fig4b <- data %>%
  ggplot(aes(x = LimitingRate_q3, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand local velocity of isotherm shifts (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Local velocity of\nisotherm\nshifts (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme, 
                                        ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, unit = "cm"))

fig4 = plot_grid(fig4a, fig4b,
          ncol = 2, align = "h")

ggsave(fig4, path = "figures", filename = "figure4_dispersal-mod-vs-lim-rate-mod.png", 
       width = 6.8, height = 3.2)

##save legend
temp = fig4a +
  theme(legend.position = "right", legend.justification = "center")

legend = ggpubr::get_legend(temp) 

ggsave(legend, path = "figures", filename = "figure4_legend.png",
       width = 2, height = 4)

## calculate ratio and stats about shift versus expansion
round(length(which(data$DispersalPotentialKmY > abs(data$ShiftKmY)))/nrow(data)*100, digits = 0) ## 57
round(length(which(data$MedianDispersalPotentialKmY > abs(data$ShiftKmY)))/nrow(data)*100, digits = 0) ## 47

mean(data$DispersalPotentialKmY[data$ShiftKmY != 0]/abs(data$ShiftKmY[data$ShiftKmY != 0]), na.rm = T)
## 97.85
median(data$DispersalPotentialKmY[data$ShiftKmY != 0]/abs(data$ShiftKmY[data$ShiftKmY != 0]), na.rm = T)
## 3.75


##########################################################
##        FIT MODELS (DISPERSAL < ISOTHERM VELOCITY)    ##
##########################################################
## fit and compete the same linear models uisng best fit climate velocity and dispersal metric to only observations of range expansion where dispersal < climate velocity
## mean+sd climate velocity was best fit
di_data <- filter(data, what_is_limiting_q3 == "Dispersal")

lme_disp_di <- lme(ShiftKmY ~ DispersalPotentialKmY,
                   random = ~1|sp_name_checked,
                   data = di_data,
                   method = "ML")

lme_cv_di <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~1|sp_name_checked,
                 data = di_data,
                 method = "ML")

## save models 
di_models_ML <- list(lme_disp_di, lme_cv_di)
names(di_models_ML) <- c("lme_disp_di", "lme_cv_di")

## get aic and rank
aic_di <- aictab(cand.set = di_models_ML, modnames = names(di_models_ML)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames) %>%
  arrange(AICc) %>%
  mutate(Rank = 1:nrow(.))

## refit using REML
lme_disp_di <- lme(ShiftKmY ~ DispersalPotentialKmY,
                   random = ~1|sp_name_checked,
                   data = di_data,
                   method = "REML")

lme_cv_di <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~1|sp_name_checked,
                 data = di_data,
                 method = "REML")

## plot residuals 
plot(lme_disp_di) 
plot(lme_cv_di)

hist(residuals(lme_disp_di))
hist(residuals(lme_cv_di))

## model summary 
summary(lme_disp_di)
summary(lme_cv_di)

## save models 
di_models <- list(lme_disp_di, lme_cv_di)
names(di_models) <- c("lme_disp_di", "lme_cv_di")

saveRDS(di_models, "data-processed/model_fits/modelfits_main_diobs_lme.rds")


#############################
##       MAKE A TABLE      ##
#############################
## get model equations
formulas <- unlist(lapply(di_models, get_formula))

model_type = ifelse(names(di_models) %in% c("lme_cv_di"), 
                    "Velocity of isotherm shift",
                    ifelse(names(di_models) %in% c("lme_disp_di"),
                           "Potential dispersal rate", 
                           "-"))

## get cond and marg r squared 
di_cond_rsq <- unlist(lapply(di_models, FUN = function(lme) {
  cond_r2 = performance::r2_nakagawa(lme)[1]
  return(cond_r2)
}))

di_marg_rsq <- unlist(lapply(di_models, FUN = function(lme) {
  marg_r2 = performance::r2_nakagawa(lme)[2]
  return(marg_r2)
}))

## get n 
n <- c(rep(nrow(di_data), 2))

## get coefs, join to r squared + se and aic table
coefs <- lapply(di_models, FUN = function(x) {summary(x)$tTable})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(di_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared_marg = marg_rsq[i]
  cur$r_squared_cond = cond_rsq[i]
  cur$n = n[i]
  cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
  coefs[[i]] <- cur
}
coefs <- coefs %>%
  bind_rows() %>%
  left_join(., aic_di) 

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## reorder/rename columns 
coefs <- coefs %>%
  rename("Conditional R2" = r_squared_cond, "Marginal R2" = r_squared_marg, "ΔAICc" = Delta_AICc,
         "AIC weight" = AICcWt, "Estimate" = Value) %>%
  select(Model, Formula, Parameter, Estimate, `Std.Error`,
         `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, LL, AICc, everything())  %>%
  select(-Cum.Wt, -Formula) %>%
  arrange(`ΔAICc`)

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, "")))

## make gt table
table_di <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - dispersal-insufficient observations"
  ) 

#################################
##       PLOT PREDICTIONS      ##
#################################
## plot model predictions
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(di_data$DispersalPotentialKmY),
                                                  by = 0.001))
pred =  predictSE(lme_disp_di, se.fit = TRUE, newdata = df_disp, level = 0)
df_disp$pred_lme <- pred$fit
df_disp$pred_lme_se <- pred$se.fit

df_cv <- data.frame(q3ClimVeloKmY_RelScale = seq(min(di_data$q3ClimVeloKmY_RelScale), max(di_data$q3ClimVeloKmY_RelScale),
                                                 by = 0.001))
pred = predictSE(lme_cv_di, se.fit = TRUE, newdata = df_cv, level = 0)
df_cv$pred_lme <- pred$fit
df_cv$pred_lme_se <- pred$se.fit

## DISPERSAL 
disp_plot_di <- di_data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme, 
                                  ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12)) 

## CLIMATE
cv_plot_di <- di_data %>%
  ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Velocity of isotherm shift (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme, 
                                ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12)) 


#########################################
##       PHYLOGENETIC CORRELATION      ##
#########################################
### check for phylogenetic correlation in residuals of the best model(s)
library(ape)
library(phylobase)
library(phylosignal)

## write out species list to use to get tree from TimeTree
write.csv(str_replace_all(unique(data$sp_name_checked), "\\_", " "), "data-processed/intermediate-files/sp_list.csv", row.names = F)

## read in TimeTree:
tree = read.tree("data-raw/other/timetree.nwk")

## transform into phylo4 object
tree <- as(tree, "phylo4")

## make dataframe of best model residuals 
## best model:
lme_limrate_q3

resid = data.frame(resid = residuals(lme_limrate_q3))

## for species with multiple residuals, take mean 
resid = resid %>%
  mutate(sp_name_checked = names(residuals(lme_limrate_q3))) %>%
  group_by(sp_name_checked) %>%
  mutate(resid = mean(resid)) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()

## get rid of species with insufficient data 
# Larus californicus (insufficient data in TimeTree to place this taxon)
# Strix varia (insufficient data in TimeTree to place this taxon)
# Larus fuscus (insufficient data in TimeTree to place this taxon)
# Carduus crispus (insufficient data in TimeTree to place this taxon)
# Sialia mexicana (insufficient data in TimeTree to place this taxon)
# Accipiter cooperii (insufficient data in TimeTree to place this taxon)
# Quercus petraea (insufficient data in TimeTree to place this taxon)
# Primula elatior (insufficient data in TimeTree to place this taxon)
# Larus delawarensis (insufficient data in TimeTree to place this taxon)
not_in_tree = c("Larus_californicus", "Strix_varia", "Larus_fuscus", "Carduus_crispus", 
                "Sialia_mexicana", "Accipiter_cooperii", "Quercus_petraea", "Primula_elatior",
                "Larus_delawarensis")

resid = resid[which(!resid$sp_name_checked %in% not_in_tree),] 

## change species names that time tree changed
resid$sp_name_checked[which(resid$sp_name_checked == "Grus_canadensis")] <- "Antigone_canadensis"

rownames(resid) = resid$sp_name_checked

resid = select(resid, resid)

## bind tip data to tree
p4 <- phylo4d(tree, tip.data = resid)

## use Pagel's lambda to test the strength of phylogenetic autocorrelation:
physig = phyloSignal(p4, methods = "Lambda")
physig$stat #Lambda = 0.02
physig$pvalue < 0.05 #FALSE

## run same test for residuals of dispersal model for dispersal-insufficient observations
## read in TimeTree:
tree = read.tree("data-raw/other/timetree.nwk")

## get sp list 
sp = tree$tip.label[which(tree$tip.label %in% di_data$sp_name_checked)]

## prune the tree 
tree <- drop.tip(tree, tree$tip.label[-match(sp, tree$tip.label)])

## transform into phylo4 object
tree <- as(tree, "phylo4")

## make dataframe of best model residuals 
## best model:
lme_disp_di

resid = data.frame(resid = residuals(lme_disp_di))

## for species with multiple residuals, take mean 
resid = resid %>%
  mutate(sp_name_checked = names(residuals(lme_disp_di))) %>%
  group_by(sp_name_checked) %>%
  mutate(resid = mean(resid)) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()

## get rid of species with insufficient data 
# Larus californicus (insufficient data in TimeTree to place this taxon)
# Strix varia (insufficient data in TimeTree to place this taxon)
# Larus fuscus (insufficient data in TimeTree to place this taxon)
# Carduus crispus (insufficient data in TimeTree to place this taxon)
# Sialia mexicana (insufficient data in TimeTree to place this taxon)
# Accipiter cooperii (insufficient data in TimeTree to place this taxon)
# Quercus petraea (insufficient data in TimeTree to place this taxon)
# Primula elatior (insufficient data in TimeTree to place this taxon)
# Larus delawarensis (insufficient data in TimeTree to place this taxon)
not_in_tree = c("Larus_californicus", "Strix_varia", "Larus_fuscus", "Carduus_crispus", 
                "Sialia_mexicana", "Accipiter_cooperii", "Quercus_petraea", "Primula_elatior",
                "Larus_delawarensis")

resid = resid[which(!resid$sp_name_checked %in% not_in_tree),] 

## change species names that time tree changed
resid$sp_name_checked[which(resid$sp_name_checked == "Grus_canadensis")] <- "Antigone_canadensis"

rownames(resid) = resid$sp_name_checked

resid = select(resid, resid)

## bind tip data to tree
p4 <- phylo4d(tree, tip.data = resid)

## use Pagel's lambda to test the strength of phylogenetic autocorrelation:
physig = phyloSignal(p4, methods = "Lambda")
physig$stat #Lambda = 0.32
physig$pvalue < 0.05 #FALSE



