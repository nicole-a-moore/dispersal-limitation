## fit main models with maximum dispersal distance, mean rate of climate change 
library(tidyverse)
library(ggplot2)
library(dotwhisker)
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

hist(dd$ClimVeloKmY_RelScale)
hist(dd$ShiftKmY)

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloKmY_RelScale,
                               DispersalPotentialKmY,
                               ClimVeloKmY_RelScale)) %>%
  mutate(ClimVeloKmY_RelScale = abs(ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY_RelScale, NA)) 

## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

## get rid of range contractions that are father than 1 sd from the mean shift 
data <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

## save data
write.csv(data, "data-processed/model-data_main_mean-cv.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear models to all observations of range expansion:
## 1. range expansion rate ~ maximum potential dispersal rate
## 2. range expansion rate ~ rate of climate change (at ecologically-relevant spatial scale)
## 2. range expansion rate ~ minimum of potential dispersal rate and rate of climate change 

lm_disp <- lm(ShiftKmY ~ DispersalPotentialKmY,
              data = data)
lm_cv <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
            data = data)
lm_limrate <- lm(ShiftKmY ~ LimitingRate,
                 data = data)

## plot residuals 
plot(lm_disp) ## not good
plot(lm_cv)
plot(lm_limrate)

hist(residuals(lm_disp))
hist(residuals(lm_cv))
hist(residuals(lm_limrate))

## model summary 
summary(lm_disp)
summary(lm_cv)
summary(lm_limrate)

## save models 
main_models <- list(lm_disp, lm_cv, lm_limrate)
names(main_models) <- c("lm_disp", "lm_cv", "lm_limrate")

saveRDS(main_models, "data-processed/modelfits_main_allobs_mean-cv.rds")

## get model equations
formulas <- unlist(lapply(main_models, get_formula))

## get r squared 
main_rsq <- unlist(lapply(main_models, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(data), 3))

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
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, everything())  %>%
  arrange(Delta_AICc)

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

gtsave(table_main, path = "figures/model_results", filename = "table_main-models_all-observations_mean-cv.png")


## plot model predictions
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY),
                                                  by = 10))
pred =  predict(lm_disp, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(data$ClimVeloKmY_RelScale), max(data$ClimVeloKmY_RelScale),
                                           by = 0.001))
pred = predict(lm_cv, se.fit = TRUE, newdata = df_cv)
df_cv$pred_lm <- pred$fit
df_cv$pred_lm_se <- pred$se.fit

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))
pred =  predict(lm_limrate, se.fit = TRUE, newdata = df_limrate)
df_limrate$pred_lm <- pred$fit
df_limrate$pred_lm_se <-  pred$se.fit


## DISPERSAL 
disp_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
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
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18))

## CLIMATE
cv_plot <- data %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18))

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate, aes(x = LimitingRate, y = pred_lm, 
                                     ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18))


plot_grid(disp_plot, cv_plot, limrate_plot, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results", filename = "main-model-predictions_all-observation_mean-cvs.png", 
       width = 12.5, height = 4)


###################################################
##        DISPERSAL-INSUFFICIENT OBSERVATIONS    ##
###################################################
## fit and compete the same linear models to only observations of range expansion where dispersal < climate velocity
di_data <- filter(data, what_is_limiting == "Dispersal")

lm_disp_di <- lm(ShiftKmY ~ DispersalPotentialKmY,
              data = di_data)
lm_cv_di <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
            data = di_data)

## plot residuals 
plot(lm_disp_di) ## not good
plot(lm_cv_di)

hist(residuals(lm_disp_di))
hist(residuals(lm_cv_di))

## model summary 
summary(lm_disp_di)
summary(lm_cv_di)

## save models 
di_models <- list(lm_disp_di, lm_cv_di)
names(di_models) <- c("lm_disp_di", "lm_cv_di")

saveRDS(di_models, "data-processed/modelfits_main_diobs_mean-cv.rds")

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
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, everything())  %>%
  arrange(Delta_AICc)

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

gtsave(table_di, path = "figures/model_results", filename = "table_main-models_di-observations_mean-cv.png")


## plot model predictions
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(di_data$DispersalPotentialKmY),
                                                  by = 0.001))
pred =  predict(lm_disp_di, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(di_data$ClimVeloKmY_RelScale), max(di_data$ClimVeloKmY_RelScale),
                                           by = 0.001))
pred = predict(lm_cv_di, se.fit = TRUE, newdata = df_cv)
df_cv$pred_lm <- pred$fit
df_cv$pred_lm_se <- pred$se.fit


## DISPERSAL 
disp_plot_di <- di_data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
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
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18))

## CLIMATE
cv_plot_di <- di_data %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18))

plot_grid(disp_plot_di, cv_plot_di, 
          ncol = 2, align = "h")

ggsave(path = "figures/model_results", filename = "main-model-predictions_di-observations_mean-cv.png", 
       width = 8.5, height = 4)

###################################################
##        DISPERSAL-SUFFICIENT OBSERVATIONS      ##
###################################################
## fit and compete the same linear models to only observations of range expansion where dispersal >climate velocity
ds_data <- filter(data, what_is_limiting == "Climate")

lm_disp_ds <- lm(ShiftKmY ~ DispersalPotentialKmY,
                 data = ds_data)
lm_cv_ds <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
               data = ds_data)

## plot residuals 
plot(lm_disp_ds)
plot(lm_cv_ds)

hist(residuals(lm_disp_ds))
hist(residuals(lm_cv_ds))

## model summary 
summary(lm_disp_ds)
summary(lm_cv_ds)

## save models 
ds_models <- list(lm_disp_ds, lm_cv_ds)
names(ds_models) <- c("lm_disp_ds", "lm_cv_ds")

saveRDS(ds_models, "data-processed/modelfits_main_dsobs_mean-cv.rds")

## get model equations
formulas <- unlist(lapply(ds_models, get_formula))

## get r squared 
ds_rsq <- unlist(lapply(ds_models, FUN = function(lm) {summary(lm)$r.squared}))

## get n 
n <- c(rep(nrow(ds_data), 2))

## get aic
aic_ds <- aictab(cand.set = ds_models, modnames = names(ds_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(ds_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(ds_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = ds_rsq[i]
  cur$n = n[i]
  cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
  coefs[[i]] <- cur
}
coefs <- coefs %>%
  bind_rows() %>%
  left_join(., aic_ds) 

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, everything())  %>%
  arrange(Delta_AICc)

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, "")))

## make gt table
table_ds <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - dispersal-sufficient observations"
  ) 

gtsave(table_ds, path = "figures/model_results", filename = "table_main-models_ds-observations_mean-cv.png")


## plot model predictions
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(ds_data$DispersalPotentialKmY),
                                                  by = 10))
pred =  predict(lm_disp_ds, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(ds_data$ClimVeloKmY_RelScale), max(ds_data$ClimVeloKmY_RelScale),
                                           by = 0.001))
pred = predict(lm_cv_ds, se.fit = TRUE, newdata = df_cv)
df_cv$pred_lm <- pred$fit
df_cv$pred_lm_se <- pred$se.fit


## DISPERSAL 
disp_plot_ds <- ds_data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
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
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18))

## CLIMATE
cv_plot_ds <- ds_data %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Aves"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Squamata"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(ds_data, is.na(colour), group == "Mammalia"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18))

plot_grid(disp_plot_ds, cv_plot_ds, 
          ncol = 2, align = "h")

ggsave(path = "figures/model_results", filename = "main-model-predictions_ds-observations_mean-cv.png", 
       width = 8.5, height = 4)



