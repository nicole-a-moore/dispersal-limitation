## fit models with median instead of maximum dispersal distance with an interaction with group (plant or bird) 
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringr)
library(gt)
library(AICcmodavg)
library(nlme)
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
  mutate(LimitingRate = ifelse(MedianDispersalPotentialKmY <= ClimVeloKmY_RelScale,
                               MedianDispersalPotentialKmY,
                               ClimVeloKmY_RelScale)) %>%
  mutate(LimitingRate_q3 = ifelse(MedianDispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                                  MedianDispersalPotentialKmY,
                                  q3ClimVeloKmY_RelScale)) %>%
  mutate(ClimVeloKmY_RelScale = abs(ClimVeloKmY_RelScale)) %>%
  mutate(q3ClimVeloKmY_RelScale = abs(q3ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(MedianDispersalPotentialKmY == LimitingRate, "Dispersal", "Climate"),
         what_is_limiting_q3 = ifelse(MedianDispersalPotentialKmY == LimitingRate_q3, "Dispersal", "Climate")) %>%
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
write.csv(contractors, "figures/extreme-contractors.csv", row.names = F)

data <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

## get rid of non-birds and non-plants
data = filter(data, group %in% c("Bird", "Plant"))

## get stats about the mean, median and max number of shifts per species:
tally = data %>%
  group_by(sp_name_checked) %>%
  tally()
max(tally$n) # 5
median(tally$n) # 1

## read in data
data = read.csv("data-processed/model-data_main_median-disp.csv")

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY*group, 
               random = ~ 1|sp_name_checked,
               data = data)

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv_median = lme(ShiftKmY ~ ClimVeloKmY_RelScale*group, 
             random = ~ 1|sp_name_checked,
             data = data)
lme_cv_q3_median <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale*group,
                 random = ~ 1|sp_name_checked,
                 data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                    random = ~ 1|sp_name_checked,
                    data = data)
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                       random = ~ 1|sp_name_checked,
                       data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate*group,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3*group,
                      random = ~ 1|sp_name_checked,
                      data = data)

## plot residuals 
plot(lme_disp_median) 
plot(lme_cv_median)
plot(lme_cv_q3_median)
plot(lme_disp_cv_median)
plot(lme_disp_cv_q3_median)
plot(lme_disp_int_median)
plot(lme_disp_int_q3_median)
plot(lme_limrate_median)
plot(lme_limrate_q3_median)

hist(residuals(lme_disp_median))
hist(residuals(lme_cv_median))
hist(residuals(lme_cv_q3_median))
hist(residuals(lme_disp_cv_median))
hist(residuals(lme_disp_cv_q3_median))
hist(residuals(lme_disp_int_median))
hist(residuals(lme_disp_int_q3_median))
hist(residuals(lme_limrate_median))
hist(residuals(lme_limrate_q3_median))

## model summary 
summary(lme_disp_median) 
summary(lme_cv_median)
summary(lme_cv_q3_median)
summary(lme_disp_cv_median)
summary(lme_disp_cv_q3_median)
summary(lme_disp_int_median)
summary(lme_disp_int_q3_median)
summary(lme_limrate_median)
summary(lme_limrate_q3_median)

## save models 
main_models <- list(lme_disp_median, lme_cv_median, lme_cv_q3_median, lme_disp_cv_median, lme_disp_cv_q3_median, 
                    lme_disp_int_median, lme_disp_int_q3_median, lme_limrate_median, lme_limrate_q3_median)
names(main_models) <- c("lme_disp_median","lme_cv_median", "lme_cv_q3_median", "lme_disp_cv_median", 
                        "lme_disp_cv_q3_median", "lme_disp_int_median", "lme_disp_int_q3_median", 
                        "lme_limrate_median", "lme_limrate_q3_median")

saveRDS(main_models, "data-processed/modelfits_main_allobs_median-disp_lme_group.rds")

## get model equations
formulas <- unlist(lapply(main_models, get_formula))

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

## get aic
aic_main <- aictab(cand.set = main_models, modnames = names(main_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, Res.LL) %>%
  rename("Model" = Modnames)

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
  select(Model, Formula, Parameter, Estimate, `Std.Error`,
        `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, Res.LL, AICc, everything())  %>%
  select(-Cum.Wt, -Formula) %>%
  arrange(`ΔAICc`)

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, "")))

## make gt table
table_main <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - all observations - median dispersal (with group interaction)"
  ) 

#gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_main-models_all-observations_lme.png")


## make table that competes all models with group interactions 
max_models <- readRDS("data-processed/modelfits_main_allobs_lme_group.rds")
median_models <- readRDS("data-processed/modelfits_main_allobs_median-disp_lme_group.rds")

## get rid of duplicated model
median_models = median_models[which(!names(median_models) %in% c("lme_cv_q3_median", "lme_cv_median"))]

## combine
main_models = append(max_models, median_models)

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

## get aic
aic_main <- aictab(cand.set = main_models, modnames = names(main_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, Res.LL) %>%
  rename("Model" = Modnames)

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
         `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, Res.LL, AICc, everything())  %>%
  arrange(`ΔAICc`) %>%
  select(-Formula, -Cum.Wt)

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  rename("Velocity of isotherm shift" = cv_type,
         "Dispersal rate" = disp_type)

## change parameter labels:
coefs$Parameter = str_replace_all(coefs$Parameter, "_q3", "")
coefs$Parameter = str_replace_all(coefs$Parameter, "q3", "")
coefs$Parameter = str_replace_all(coefs$Parameter, "DispersalPotentialKmY", "Potential dispersal rate")
coefs$Parameter = str_replace_all(coefs$Parameter, "MedianPotential", "Potential")
coefs$Parameter = str_replace_all(coefs$Parameter, ":ClimVeloKmY_RelScale", ":velocity of isotherm shift")
coefs$Parameter = str_replace_all(coefs$Parameter, "ClimVeloKmY_RelScale", "Velocity of isotherm shift")
coefs$Parameter = str_replace_all(coefs$Parameter, ":groupPlant", ":group (plants)")
coefs$Parameter = str_replace_all(coefs$Parameter, "groupPlant", "Group: plants")
coefs$Parameter = str_replace_all(coefs$Parameter, "\\(Intercept\\)", "Intercept")
coefs$Parameter = str_replace_all(coefs$Parameter, "LimitingRate", "Minimum rate")
 
## make gt table
table_main <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - all observations (with group interactions)"
  ) 

gtsave(table_main, path = "figures/model_results/all-observations", filename = "tableS2_all-models_lme_by-group.png")
gtsave(table_main, path = "figures/model_results/all-observations", filename = "tableS2_all-models_lme_by-group.docx")

## make table summarizing results by model 'type' 
## calculate cumulative weight and mean rank per model 'type'
tableS3 = coefs %>%
  select(Model, `Dispersal rate`, `Velocity of isotherm shift`, `AIC weight`) %>%
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

tableS3 = tableS3 %>% gt()

## save table 
gtsave(tableS3, path = "figures/model_results/all-observations", filename = "tableS3_by-group.docx")


