## fit models with median instead of maximum dispersal distance to observations INCLUDING EXTREME CONTRACTIONS
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringr)
library(gt)
library(AICcmodavg)
theme_set(theme_bw())

get_formula <- function(lm) {
  return(as.character(lm$call)[2])
}


#############################
##       PREPARE DATA      ##
#############################
## read in data 
dd <- read.csv("data-processed/v3_with-cv.csv")

## filter to leading edge shifts with positive climate velocity 
dd <- filter(dd, !Param %in% c("O", "TE") & ClimVeloKmY_RelScale >= 0)

dd <- dd %>%
  mutate(LimitingRate = ifelse(MedianDispersalPotentialKmY <= ClimVeloKmY_RelScale,
                                     MedianDispersalPotentialKmY,
                               ClimVeloKmY_RelScale),
         LimitingRate_p90 = ifelse(MedianDispersalPotentialKmY <= p90ClimVeloKmY_RelScale,
                                     MedianDispersalPotentialKmY,
                               p90ClimVeloKmY_RelScale)) %>%
  mutate(p90ClimVeloKmY_RelScale = abs(p90ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(MedianDispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_p90 = ifelse(MedianDispersalPotentialKmY == LimitingRate_p90, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY_RelScale, NA)) %>%
  mutate(colour_p90 = ifelse(what_is_limiting_p90 == "Climate", p90ClimVeloKmY_RelScale, NA)) 

## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

data = dd

## save data
write.csv(data, "data-processed/model-data_main_median-disp.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear models to all observations of range expansion using median potential dispersal rate
## fit each using mean and 90th percentile of climate velocity 

## 1. range expansion rate ~ potential dispersal rate
lm_disp_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY,
              data = data)

## 2. range expansion rate ~ velocity of climate change (at ecologically-relevant spatial scale)
lm_cv_median <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
            data = data)
lm_cv_p90_median <- lm(ShiftKmY ~ p90ClimVeloKmY_RelScale,
                data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of climate change
lm_disp_cv_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY + ClimVeloKmY_RelScale,
                 data = data)
lm_disp_cv_p90_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY + p90ClimVeloKmY_RelScale,
                     data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of climate change 
lm_disp_int_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale,
                  data = data)
lm_disp_int_p90_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY*p90ClimVeloKmY_RelScale,
                      data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of climate change 
lm_limrate_median <- lm(ShiftKmY ~ LimitingRate,
                 data = data)
lm_limrate_p90_median <- lm(ShiftKmY ~ LimitingRate_p90,
                     data = data)

## plot residuals 
plot(lm_disp_median) 
plot(lm_cv_median)
plot(lm_cv_p90_median)
plot(lm_disp_cv_median)
plot(lm_disp_cv_p90_median)
plot(lm_disp_int_median)
plot(lm_disp_int_p90_median)
plot(lm_limrate_median)
plot(lm_limrate_p90_median)

hist(residuals(lm_disp_median))
hist(residuals(lm_cv_median))
hist(residuals(lm_cv_p90_median))
hist(residuals(lm_disp_cv_median))
hist(residuals(lm_disp_cv_p90_median))
hist(residuals(lm_disp_int_median))
hist(residuals(lm_disp_int_p90_median))
hist(residuals(lm_limrate_median))
hist(residuals(lm_limrate_p90_median))

## model summary 
summary(lm_disp_median) 
summary(lm_cv_median)
summary(lm_cv_p90_median)
summary(lm_disp_cv_median)
summary(lm_disp_cv_p90_median)
summary(lm_disp_int_median)
summary(lm_disp_int_p90_median)
summary(lm_limrate_median)
summary(lm_limrate_p90_median)

## save models 
main_models <- list(lm_disp_median, lm_cv_median, lm_cv_p90_median, lm_disp_cv_median, lm_disp_cv_p90_median, 
                    lm_disp_int_median, lm_disp_int_p90_median, lm_limrate_median, lm_limrate_p90_median)
names(main_models) <- c("lm_disp_median","lm_cv_median", "lm_cv_p90_median", "lm_disp_cv_median", 
                        "lm_disp_cv_p90_median", "lm_disp_int_median", "lm_disp_int_p90_median", 
                        "lm_limrate_median", "lm_limrate_p90_median")

saveRDS(main_models, "data-processed/modelfits_main_allobs_median-disp_w-contractions.rds")

## get model equations
formulas <- unlist(lapply(main_models, get_formula))

## get r squared 
main_rsq <- unlist(lapply(main_models, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(data), length(main_models)))

## get aic
aic_main <- aictab(cand.set = main_models, modnames = names(main_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(main_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(main_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = main_rsq[i]
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
  rename("R2" = r_squared, "p-value" = `Pr(>|t|)`, "t-value" = `t value`, "ΔAICc" = Delta_AICc,
         "AIC weight" = AICcWt) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, `t-value`, 
         `p-value`, `R2`, n, K, LL, AICc, everything())  %>%
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
    title = "Main model set - all observations - median dispersal"
  ) 

#gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_main-models_all-observations_median-disp.png")


## make table that competes all models with all observations
max_models <- readRDS("data-processed/modelfits_main_allobs_w-contractions.rds")
median_models <- readRDS("data-processed/modelfits_main_allobs_median-disp_w-contractions.rds")

## get rid of duplicated model
median_models = median_models[which(!names(median_models) %in% c("lm_cv_p90_median", "lm_cv_median"))]

## combine
main_models = append(max_models, median_models)

## get model equations
formulas <- unlist(lapply(main_models, get_formula))

## get types
disp_type = ifelse(names(main_models) %in% c("lm_disp", "lm_disp_cv", "lm_disp_cv_p90", 
                                             "lm_disp_int", "lm_disp_int_p90",
                                             "lm_limrate", "lm_limrate_p90"), 
                   "Maximum",
                   ifelse(names(main_models) %in% c("lm_disp_median", "lm_disp_cv_median", "lm_disp_cv_p90_median",
                                                    "lm_disp_int_median", "lm_disp_int_p90_median",
                                                    "lm_limrate_p90_median", "lm_limrate_median"),
                          "Median",
                          "-"))
cv_type = ifelse(names(main_models) %in% c("lm_cv", "lm_disp_cv",  "lm_disp_int",
                                            "lm_disp_cv_median", "lm_disp_int_median",
                                           "lm_limrate", "lm_limrate_median"), 
                 "Mean",
                 ifelse(names(main_models) %in% c("lm_cv_p90", "lm_disp_cv_p90",  "lm_disp_int_p90",
                                                  "lm_disp_cv_p90_median", "lm_disp_int_p90_median",
                                                  "lm_limrate_p90", "lm_limrate_p90_median"),
                        "90th percentile",
                        "-"))
model_type = ifelse(names(main_models) %in% c("lm_cv", "lm_cv_p90"), 
                 "Velocity of climate change",
                 ifelse(names(main_models) %in% c("lm_limrate_p90_median", "lm_limrate_p90", 
                                                  "lm_limrate", "lm_limrate_median"),
                        "Minimum rate",
                        ifelse(names(main_models) %in% c("lm_disp", "lm_disp_p90", "lm_disp_p90_median", "lm_disp_median"),
                               "Potential dispersal rate", 
                               ifelse(names(main_models) %in% c("lm_disp_cv", "lm_disp_cv_p90", 
                                                                "lm_disp_cv_p90_median", "lm_disp_cv_median"),
                                      "Potential dispersal rate and velocity of climate change (additive)",
                                      ifelse(names(main_models) %in% c("lm_disp_int", "lm_disp_int_p90", 
                                                                       "lm_disp_int_median", "lm_disp_int_p90_median"),
                                             "Potential dispersal rate and velocity of climate change (interactive)",
                                             "-")))))

## get r squared 
main_rsq <- unlist(lapply(main_models, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(data), length(main_models)))

## get aic
aic_main <- aictab(cand.set = main_models, modnames = names(main_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(main_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(main_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = main_rsq[i]
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
  rename("R2" = r_squared, "p-value" = `Pr(>|t|)`, "t-value" = `t value`, "ΔAICc" = Delta_AICc,
         "AIC weight" = AICcWt) %>%
  select(-Model) %>%
  rename("Model" = model_type) %>%
  select(Model, cv_type, disp_type, Formula, Parameter, Estimate, `Std. Error`, `t-value`, 
         `p-value`, `R2`, n, K, LL, AICc, everything())  %>%
  arrange(`ΔAICc`) %>%
  mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", 
                            ifelse(str_detect(Parameter, "\\:"), 
                                          "Potential dispersal rate:velocity of climate change", 
                                   ifelse(str_detect(Parameter, "ClimVeloKmY"), "Velocity of climate change",
                                          ifelse(str_detect(Parameter, "DispersalPotential"), 
                                                 "Potential dispersal rate", 
                                                 "Minimum rate"))))) %>%
  select(-Cum.Wt, -Formula) 

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  rename("Velocity of climate change" = cv_type,
         "Dispersal rate" = disp_type)

## make gt table
table_main <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - all observations"
  ) 

gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_all-models_w-contractions.png")
gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_all-models_w-contractions.docx")


###################################################
##        DISPERSAL-INSUFFICIENT OBSERVATIONS    ##
###################################################
## fit and compete the same linear models uisng best fit climate velocity to only observations of range expansion where dispersal < climate velocity
## p90 climate velocity was best fit
di_data <- filter(data, what_is_limiting_p90 == "Dispersal")

lm_disp_di <- lm(ShiftKmY ~ MedianDispersalPotentialKmY,
                 data = di_data)
lm_cv_di <- lm(ShiftKmY ~ p90ClimVeloKmY_RelScale,
               data = di_data)

## plot residuals 
plot(lm_disp_di)
plot(lm_cv_di)

hist(residuals(lm_disp_di))
hist(residuals(lm_cv_di))

## model summary 
summary(lm_disp_di)
summary(lm_cv_di)

## save models 
di_models <- list(lm_disp_di, lm_cv_di)
names(di_models) <- c("lm_disp_di", "lm_cv_di")

saveRDS(di_models, "data-processed/modelfits_main_diobs_median-disp_w-contractions.rds")

## get model equations
formulas <- unlist(lapply(di_models, get_formula))

## get r squared 
di_rsq <- unlist(lapply(di_models, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(di_data), 2))

## get aic
aic_di <- aictab(cand.set = di_models, modnames = names(di_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(di_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(di_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = di_rsq[i]
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
  rename("R2" = r_squared, "p-value" = `Pr(>|t|)`, "t-value" = `t value`, "ΔAICc" = Delta_AICc,
         "AIC weight" = AICcWt) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, `t-value`, 
         `p-value`, `R2`, n, K, LL, AICc, everything())  %>%
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
    title = "Main model set - dispersal-insufficient observations - median dispersal"
  ) 

gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q90_median-disp_w-contractions.png")
