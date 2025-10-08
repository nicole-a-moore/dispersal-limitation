## fit main models with maximum dispersal distance with an interaction with group (plant or bird) 
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
data <- read.csv("data-processed/model-data_main.csv")

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence
## include interaction term with group 
library(nlme)

## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY*group,
               random = ~ 1|sp_name_checked,
               data = data)

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale*group,
               random = ~ 1|sp_name_checked,
               data = data)
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale*group,
                random = ~ 1|sp_name_checked,
                data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY*group + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*group + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale*group,
                    random = ~ 1|sp_name_checked,
                    data = data)
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale*group,
                       random = ~ 1|sp_name_checked,
                       data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift
lme_limrate <- lme(ShiftKmY ~ LimitingRate*group,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3*group,
                      random = ~ 1|sp_name_checked,
                      data = data)

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

hist(residuals(lme_disp))
hist(residuals(lme_cv))
hist(residuals(lme_cv_q3))
hist(residuals(lme_disp_cv))
hist(residuals(lme_disp_cv_q3))
hist(residuals(lme_disp_int))
hist(residuals(lme_disp_int_q3))
hist(residuals(lme_limrate))
hist(residuals(lme_limrate_q3))

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

## save models 
main_models <- list(lme_disp, lme_cv, lme_cv_q3, lme_disp_cv, lme_disp_cv_q3, 
                    lme_disp_int, lme_disp_int_q3, lme_limrate, lme_limrate_q3)
names(main_models) <- c("lme_disp","lme_cv", "lme_cv_q3", "lme_disp_cv", "lme_disp_cv_q3", 
                        "lme_disp_int", "lme_disp_int_q3", "lme_limrate", "lme_limrate_q3")

saveRDS(main_models, "data-processed/modelfits_main_allobs_lme_group.rds")

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
    title = "Main model set - all observations (with group interaction)"
  ) 

#gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_main-models_all-observations_lme.png")

## plot predictions for best model
df_limrate_birds <- data.frame(expand.grid(LimitingRate = seq(min(data$LimitingRate[which(data$group == "Bird")]),
                                                              max(data$LimitingRate[which(data$group == "Bird")]),
                                            by = 0.001), group = unique(data$group)[1]))
df_limrate_plants <- data.frame(expand.grid(LimitingRate = seq(min(data$LimitingRate[which(data$group == "Plant")]),
                                                              max(data$LimitingRate[which(data$group == "Plant")]),
                                                              by = 0.001), group = unique(data$group)[2]))
df_limrate = rbind(df_limrate_birds, df_limrate_plants)
pred =  predictSE(lme_limrate, se.fit = TRUE, newdata = df_limrate, level = 0)
df_limrate$pred_lme <- pred$fit
df_limrate$pred_lme_se <-  pred$se.fit

data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  # stat_function(colour = "grey", fun = function(x){x},
  #               linetype = "dashed") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand velocity of isotherm shift (km/yr)",
       y = "Estimated range expansion rate (km/yr)",
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate, aes(x = LimitingRate, y = pred_lme,
                                     ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lme, group = group),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  # annotate("text", x = 7, y = 23, hjust = 0,
  #          label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_limrate")[1]],
  #                         "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_limrate")[1]]), colour = "black",
  #          size = 3.5) +
  facet_wrap(~group, scales = "free_x")

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_cv_q3_lme_by-group.png", 
       width = 9.2, height = 3.2)

