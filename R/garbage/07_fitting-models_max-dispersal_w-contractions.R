## fit main models with maximum dispersal distance to observations INCLUDING EXTREME CONTRACTIONS
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

## filter to leading edge shifts with positive mean climate velocity 
dd <- filter(dd, !Param %in% c("O", "TE") & ClimVeloKmY_RelScale >= 0)

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloKmY_RelScale,
                               DispersalPotentialKmY,
                               ClimVeloKmY_RelScale)) %>%
  mutate(LimitingRate_q3 = ifelse(DispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                               DispersalPotentialKmY,
                               q3ClimVeloKmY_RelScale)) %>%
  mutate(ClimVeloKmY_RelScale = abs(ClimVeloKmY_RelScale)) %>%
  mutate(q3ClimVeloKmY_RelScale = abs(q3ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate"),
         what_is_limiting_q3 = ifelse(DispersalPotentialKmY == LimitingRate_q3, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY_RelScale, NA),
         colour_q3 = ifelse(what_is_limiting_q3 == "Climate", q3ClimVeloKmY_RelScale, NA)) 

## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

data = dd 

## get rid of non-birds and non-plants
data = filter(data, group %in% c("Bird", "Plant"))

## save data
write.csv(data, "data-processed/model-data_with-contractions.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity 

## 1. range expansion rate ~ potential dispersal rate
lm_disp <- lm(ShiftKmY ~ DispersalPotentialKmY,
              data = data)

## 2. range expansion rate ~ velocity of climate change (at ecologically-relevant spatial scale)
lm_cv <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
            data = data)
lm_cv_q3 <- lm(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of climate change
lm_disp_cv <- lm(ShiftKmY ~ DispersalPotentialKmY + ClimVeloKmY_RelScale,
                 data = data)
lm_disp_cv_q3 <- lm(ShiftKmY ~ DispersalPotentialKmY + q3ClimVeloKmY_RelScale,
                     data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of climate change 
lm_disp_int <- lm(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale,
                  data = data)
lm_disp_int_q3 <- lm(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                      data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of climate change 
lm_limrate <- lm(ShiftKmY ~ LimitingRate,
                 data = data)
lm_limrate_q3 <- lm(ShiftKmY ~ LimitingRate_q3,
                     data = data)


## plot residuals 
plot(lm_disp) 
plot(lm_cv)
plot(lm_cv_q3)
plot(lm_disp_cv)
plot(lm_disp_cv_q3)
plot(lm_disp_int)
plot(lm_disp_int_q3)
plot(lm_limrate)
plot(lm_limrate_q3)

hist(residuals(lm_disp))
hist(residuals(lm_cv))
hist(residuals(lm_cv_q3))
hist(residuals(lm_disp_cv))
hist(residuals(lm_disp_cv_q3))
hist(residuals(lm_disp_int))
hist(residuals(lm_disp_int_q3))
hist(residuals(lm_limrate))
hist(residuals(lm_limrate_q3))

## model summary 
summary(lm_disp) 
summary(lm_cv)
summary(lm_cv_q3)
summary(lm_disp_cv)
summary(lm_disp_cv_q3)
summary(lm_disp_int)
summary(lm_disp_int_q3)
summary(lm_limrate)
summary(lm_limrate_q3)

## save models 
main_models <- list(lm_disp, lm_cv, lm_cv_q3, lm_disp_cv, lm_disp_cv_q3, 
                    lm_disp_int, lm_disp_int_q3, lm_limrate, lm_limrate_q3)
names(main_models) <- c("lm_disp","lm_cv", "lm_cv_q3", "lm_disp_cv", "lm_disp_cv_q3", 
                        "lm_disp_int", "lm_disp_int_q3", "lm_limrate", "lm_limrate_q3")

saveRDS(main_models, "data-processed/modelfits_main_allobs_w-contractions.rds")

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
    title = "Main model set - all observations"
  ) 

#gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_main-models_all-observations.png")


###################################################
##        DISPERSAL-INSUFFICIENT OBSERVATIONS    ##
###################################################
## fit and compete the same linear models uisng best fit climate velocity to only observations of range expansion where dispersal < climate velocity
## mean+sd climate velocity was best fit
di_data <- filter(data, what_is_limiting_q3 == "Dispersal")

lm_disp_di <- lm(ShiftKmY ~ DispersalPotentialKmY,
              data = di_data)
lm_cv_di <- lm(ShiftKmY ~ q3ClimVeloKmY_RelScale,
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

saveRDS(di_models, "data-processed/modelfits_main_diobs_w-contractions.rds")

## get model equations
formulas <- unlist(lapply(di_models, get_formula))

model_type = ifelse(names(di_models) %in% c("lm_cv_di"), 
                    "Velocity of climate change",
                           ifelse(names(di_models) %in% c("lm_disp_di"),
                                  "Potential dispersal rate", 
                                  "-"))

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
  cur$model_type = model_type[i]
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
  select(-Model) %>%
  rename("Model" = model_type) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, `t-value`, 
         `p-value`, `R2`, n, K, LL, AICc, everything())  %>%
  mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", 
                            ifelse(str_detect(Parameter, "q3ClimVeloKmY_RelScale"), "Velocity of climate change",
                                   ifelse(str_detect(Parameter, "DispersalPotentialKmY"), "Potential dispersal rate", 
                                          "Minimum rate")))) %>%
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

gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3_w-contractions.png")
gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3_w-contractions.docx")


## repeat with mean cv
di_data_mean <- filter(data, what_is_limiting == "Dispersal")

lm_disp_di_mean <- lm(ShiftKmY ~ DispersalPotentialKmY,
                 data = di_data_mean)
lm_cv_di_mean <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
               data = di_data_mean)

## plot residuals 
plot(lm_disp_di_mean) 
plot(lm_cv_di_mean)

hist(residuals(lm_disp_di_mean))
hist(residuals(lm_cv_di_mean))

## model summary 
summary(lm_disp_di_mean)
summary(lm_cv_di_mean)

## save models 
di_models_mean <- list(lm_disp_di_mean, lm_cv_di_mean)
names(di_models_mean) <- c("lm_disp_di_mean", "lm_cv_di_mean")

saveRDS(di_models_mean, "data-processed/modelfits_main_diobs-mean_w-contractions.rds")

## get model equations
formulas <- unlist(lapply(di_models_mean, get_formula))

## get r squared 
di_rsq <- unlist(lapply(di_models_mean, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(di_data_mean), 2))

## get aic
aic_di <- aictab(cand.set = di_models_mean, modnames = names(di_models_mean)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(di_models_mean, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(di_models_mean)[i]
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
table_di_mean <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - dispersal-insufficient observations"
  ) 

gtsave(table_di_mean, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_mean_w-contractions.png")

