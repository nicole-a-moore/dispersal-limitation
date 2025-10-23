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

## filter to leading edge shifts 
dd = filter(dd, !Param %in% c("O", "TE"))

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

## get rid of non-birds and non-plants
data = filter(dd, group %in% c("Bird", "Plant"))

## save data
write.csv(data, "data-processed/model_data/model-data_contractions.csv", row.names = FALSE)

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

saveRDS(main_models_ML, "data-processed/model_fits/modelfits_w-contractions_lme_ML.rds")

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

saveRDS(main_models, "data-processed/model_fits/modelfits_w-contractions_lme_REML.rds")


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
    title = "Main model set - all observations (with contractions)"
  ) 

gtsave(table_main, path = "figures", filename = "tableS4_all-models_lme_w-contractions.png")
gtsave(table_main, path = "figures", filename = "tableS4_all-models_lme_w-contractions.docx")

## make table summarizing results by model 'type' 
## calculate cumulative weight and mean rank per model 'type'
tableS5 = coefs %>%
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

tableS5 = tableS5 %>% gt()

## save table 
gtsave(tableS5, path = "figures", filename = "tableS5_summary-w-contractions.docx")
