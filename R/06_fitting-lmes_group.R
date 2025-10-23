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
data <- read.csv("data-processed/model_data/model-data_main.csv")

###################################
##          FIT MODELS           ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY*group, 
               random = ~ 1|sp_name_checked,
               data = data, 
               method = "ML")
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY*group, 
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale*group, 
             random = ~ 1|sp_name_checked,
             data = data,
             method = "ML")
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale*group,
                 random = ~ 1|sp_name_checked,
                 data = data,
                 method = "ML")

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "ML")
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*group  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "ML")
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group  + q3ClimVeloKmY_RelScale,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "ML")

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                    random = ~ 1|sp_name_checked,
                    data = data,
                    method = "ML")
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                       random = ~ 1|sp_name_checked,
                       data = data,
                       method = "ML")
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                           random = ~ 1|sp_name_checked,
                           data = data,
                           method = "ML")
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                              random = ~ 1|sp_name_checked,
                              data = data,
                              method = "ML")

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate <- lme(ShiftKmY ~ LimitingRate*group,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "ML")
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3*group,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "ML")
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate_median*group,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "ML")
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3_median*group,
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

saveRDS(main_models_ML, "data-processed/model_fits/modelfits_main_lme_group_ML.rds")

## get aic and rank
aic_main <- aictab(cand.set = main_models_ML, modnames = names(main_models_ML)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames) %>%
  arrange(AICc) %>%
  mutate(Rank = 1:nrow(.))

## refit using REML 
## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY*group, 
               random = ~ 1|sp_name_checked,
               data = data, 
               method = "REML")
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY*group, 
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale*group, 
             random = ~ 1|sp_name_checked,
             data = data,
             method = "REML")
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale*group,
                 random = ~ 1|sp_name_checked,
                 data = data,
                 method = "REML")

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "REML")
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*group  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "REML")
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*group  + q3ClimVeloKmY_RelScale,
                             random = ~ 1|sp_name_checked,
                             data = data,
                             method = "REML")

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                    random = ~ 1|sp_name_checked,
                    data = data,
                    method = "REML")
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                       random = ~ 1|sp_name_checked,
                       data = data,
                       method = "REML")
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                           random = ~ 1|sp_name_checked,
                           data = data,
                           method = "REML")
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                              random = ~ 1|sp_name_checked,
                              data = data,
                              method = "REML")

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate <- lme(ShiftKmY ~ LimitingRate*group,
                   random = ~ 1|sp_name_checked,
                   data = data,
                   method = "REML")
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3*group,
                      random = ~ 1|sp_name_checked,
                      data = data,
                      method = "REML")
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate_median*group,
                          random = ~ 1|sp_name_checked,
                          data = data,
                          method = "REML")
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3_median*group,
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

saveRDS(main_models, "data-processed/model_fits/modelfits_main_lme_group_REML.rds")


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

gtsave(table_main, path = "figures", filename = "tableS2_all-models_lme_by-group.png")
gtsave(table_main, path = "figures", filename = "tableS2_all-models_lme_by-group.docx")


## make table summarizing results by model 'type' 
## calculate cumulative weight and mean rank per model 'type'
tableS3 = coefs %>%
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

tableS3 = table3 %>% gt()

## save table 
gtsave(tableS3, path = "figures", filename = "tableS3_by-group.docx")
