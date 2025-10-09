## fit main models with maximum dispersal distance INCLUDING EXTREME CONTRACTIONS
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
  mutate(LimitingRate_q3 = ifelse(DispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                                  DispersalPotentialKmY,
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
write.csv(data, "data-processed/model-data_contractions.csv", row.names = FALSE)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY, 
               random = ~ 1|sp_name_checked,
               data = data)

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv = lme(ShiftKmY ~ ClimVeloKmY_RelScale, 
             random = ~ 1|sp_name_checked,
             data = data)
lme_cv_q3 <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~ 1|sp_name_checked,
                 data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv <- lme(ShiftKmY ~ DispersalPotentialKmY + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_disp_cv_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale,
                    random = ~ 1|sp_name_checked,
                    data = data)
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                       random = ~ 1|sp_name_checked,
                       data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate <- lme(ShiftKmY ~ LimitingRate,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_limrate_q3 <- lme(ShiftKmY ~ LimitingRate_q3,
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

saveRDS(main_models, "data-processed/modelfits_main_allobs_lme_w-contractions.rds")

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
    title = "Main model set - all observations (with contractions)"
  ) 

## plot model predictions
df_disp <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                  ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred = predictSE(lme_disp, se.fit = T, newdata = df_disp, level = 0)
df_disp$pred_lme <- pred$fit
df_disp$pred_lme_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(data$ClimVeloKmY_RelScale), max(data$ClimVeloKmY_RelScale),
                                               by = 0.001))
pred = predictSE(lme_cv, se.fit = TRUE, newdata = df_cv, level = 0)
df_cv$pred_lme <- pred$fit
df_cv$pred_lme_se <- pred$se.fit

df_cv_q3 <- data.frame(q3ClimVeloKmY_RelScale = seq(min(data$q3ClimVeloKmY_RelScale), 
                                                    max(data$q3ClimVeloKmY_RelScale),
                                                    by = 0.001))
pred = predictSE(lme_cv_q3, se.fit = TRUE, newdata = df_cv_q3, level = 0)
df_cv_q3$pred_lme <- pred$fit
df_cv_q3$pred_lme_se <- pred$se.fit

df_disp_cv <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                     ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_cv, se.fit = TRUE, newdata = df_disp_cv, level = 0)
df_disp_cv$pred_lme <- pred$fit
df_disp_cv$pred_lme_se <- pred$se.fit

df_disp_cv_q3 <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                        q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_cv_q3, se.fit = TRUE, newdata = df_disp_cv_q3, level = 0)
df_disp_cv_q3$pred_lme <- pred$fit
df_disp_cv_q3$pred_lme_se <- pred$se.fit


df_disp_int <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                      ClimVeloKmY_RelScale = c(min(data$ClimVeloKmY_RelScale),
                                                               max(data$ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int, se.fit = TRUE, newdata = df_disp_int, level = 0)
df_disp_int$pred_lme <- pred$fit
df_disp_int$pred_lme_se <- pred$se.fit

df_disp_int_q3 <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                         q3ClimVeloKmY_RelScale = c(min(data$q3ClimVeloKmY_RelScale),
                                                                    max(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int_q3, se.fit = TRUE, newdata = df_disp_int_q3, level = 0)
df_disp_int_q3$pred_lme <- pred$fit
df_disp_int_q3$pred_lme_se <- pred$se.fit

df_disp_int_q3_mean <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                              q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int_q3, se.fit = TRUE, newdata = df_disp_int_q3_mean, level = 0)
df_disp_int_q3_mean$pred_lme <- pred$fit
df_disp_int_q3_mean$pred_lme_se <- pred$se.fit

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))
pred =  predictSE(lme_limrate, se.fit = TRUE, newdata = df_limrate, level = 0)
df_limrate$pred_lme <- pred$fit
df_limrate$pred_lme_se <-  pred$se.fit

df_limrate_q3 <- data.frame(LimitingRate_q3 = seq(min(data$LimitingRate_q3), max(data$LimitingRate_q3),
                                                  by = 0.001))
pred =  predictSE(lme_limrate_q3, se.fit = TRUE, newdata = df_limrate_q3, level = 0)
df_limrate_q3$pred_lme <- pred$fit
df_limrate_q3$pred_lme_se <-  pred$se.fit

data$cv_lab_mean = "Mean climate velocity"
data$cv_lab_q3 = "3rd quartile climate velocity"


## DISPERSAL 
disp_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  # geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                 ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp")[1]]), colour = "black",
           size = 3.5)

## CLIMATE
cv_plot <- data %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
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
  geom_ribbon(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lme, 
                                ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_cv")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_cv")[1]]), colour = "black",
           size = 3.5)

cv_plot_q3 <- data %>%
  ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  geom_ribbon(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme, 
                                   ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_cv_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_cv_q3")[1]]), colour = "black",
           size = 3.5)


## DISPERSAL AND CLIMATE - ADDITIVE
disp_cv_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  # geom_ribbon(data = df_disp_cv, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                    ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_cv, aes(x = DispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_cv")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_cv")[1]]), colour = "black",
           size = 3.5)

disp_cv_q3_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  # geom_ribbon(data = df_disp_cv_q3, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                       ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_cv_q3, aes(x = DispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 700, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_cv_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_cv_q3")[1]]), colour = "black",
           size = 3.5)

## DISPERSAL AND CLIMATE - INTERACTIVE
disp_int_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  # geom_ribbon(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                     ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se,
  #                                     group = ClimVeloKmY_RelScale), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_smooth(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                     colour = ClimVeloKmY_RelScale, group = ClimVeloKmY_RelScale),
  #             inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_int")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_int")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
  # geom_ribbon(data = df_disp_int_q3, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                        ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se,
  #                                        group = q3ClimVeloKmY_RelScale), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_smooth(data = df_disp_int_q3, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                        colour = q3ClimVeloKmY_RelScale, group = q3ClimVeloKmY_RelScale),
  #            inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_int_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_int_q3")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_mean_plot <- data %>%
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
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                             ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 650, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_int_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_int_q3")[1]]), colour = "black",
           size = 3.5)

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand velocity of isotherm shift (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate, aes(x = LimitingRate, y = pred_lme, 
                                     ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_limrate")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_limrate")[1]]), colour = "black",
           size = 3.5)

limrate_plot_q3 <- data %>%
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
  labs(x = "Minimum of potential dispersal rate\nand velocity of isotherm shift (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme, 
                                        ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_limrate_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_limrate_q3")[1]]), colour = "black",
           size = 3.5) +
  scale_x_continuous(limits = c(0, 12.1)) 

plot_grid(disp_int_q3_mean_plot, cv_plot_q3, limrate_plot_q3, 
          ncol = 3, align = "h")

###################################################
##        DISPERSAL-INSUFFICIENT OBSERVATIONS    ##
###################################################
## fit and compete the same linear models uisng best fit climate velocity to only observations of range expansion where dispersal < climate velocity
## mean+sd climate velocity was best fit
di_data <- filter(data, what_is_limiting_q3 == "Dispersal")

lme_disp_di <- lme(ShiftKmY ~ DispersalPotentialKmY,
                   random = ~1|sp_name_checked,
                   data = di_data)
lme_cv_di <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~1|sp_name_checked,
                 data = di_data)

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

saveRDS(di_models, "data-processed/modelfits_main_diobs_lme_w-contractions.rds")

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

## get aic
aic_di <- aictab(cand.set = di_models, modnames = names(di_models)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, Res.LL) %>%
  rename("Model" = Modnames)

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
         `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, Res.LL, AICc, everything())  %>%
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
    title = "Main model set - dispersal-insufficient observations (with contractions)"
  ) 

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

coefs_saved = coefs

## repeat with mean cv
di_data_mean <- filter(data, what_is_limiting == "Dispersal")

lme_disp_di_mean <- lme(ShiftKmY ~ DispersalPotentialKmY,
                        random = ~1|sp_name_checked,
                        data = di_data_mean)
lme_cv_di_mean <- lme(ShiftKmY ~ ClimVeloKmY_RelScale,
                      random = ~1|sp_name_checked,
                      data = di_data_mean)

## plot residuals 
plot(lme_disp_di_mean) 
plot(lme_cv_di_mean)

hist(residuals(lme_disp_di_mean))
hist(residuals(lme_cv_di_mean))

## model summary 
summary(lme_disp_di_mean)
summary(lme_cv_di_mean)

## save models 
di_models_mean <- list(lme_disp_di_mean, lme_cv_di_mean)
names(di_models_mean) <- c("lme_disp_di_mean", "lme_cv_di_mean")

saveRDS(di_models_mean, "data-processed/modelfits_main_diobs-mean_lme_w-contractions.rds")

## get model equations
formulas <- unlist(lapply(di_models_mean, get_formula))

## get cond and marg r squared 
di_cond_rsq <- unlist(lapply(di_models_mean, FUN = function(lme) {
  cond_r2 = performance::r2_nakagawa(lme)[1]
  return(cond_r2)
}))

di_marg_rsq <- unlist(lapply(di_models_mean, FUN = function(lme) {
  marg_r2 = performance::r2_nakagawa(lme)[2]
  return(marg_r2)
}))

## get n 
n <- c(rep(nrow(di_data), 2))

## get aic
aic_di <- aictab(cand.set = di_models_mean, modnames = names(di_models_mean)) %>%
  data.frame() %>%
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, Res.LL) %>%
  rename("Model" = Modnames)

## get coefs, join to r squared + se and aic table
coefs <- lapply(di_models_mean, FUN = function(x) {summary(x)$tTable})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(di_models_mean)[i]
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
         `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, Res.LL, AICc, everything())  %>%
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
    title = "Main model set - dispersal-insufficient observations (with contractions)"
  ) 

## plot model predictions
df_disp_mean <- data.frame(DispersalPotentialKmY = seq(0, max(di_data_mean$DispersalPotentialKmY),
                                                       by = 0.001))
pred =  predictSE(lme_disp_di_mean, se.fit = TRUE, newdata = df_disp_mean, level = 0)
df_disp_mean$pred_lme <- pred$fit
df_disp_mean$pred_lme_se <- pred$se.fit

df_cv_mean <- data.frame(ClimVeloKmY_RelScale = seq(min(di_data_mean$ClimVeloKmY_RelScale), 
                                                    max(di_data_mean$ClimVeloKmY_RelScale),
                                                    by = 0.001))
pred = predictSE(lme_cv_di_mean, se.fit = TRUE, newdata = df_cv_mean, level = 0)
df_cv_mean$pred_lme <- pred$fit
df_cv_mean$pred_lme_se <- pred$se.fit


## DISPERSAL 
disp_plot_di_mean <- di_data_mean %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
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
cv_plot_di_mean <- di_data_mean %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
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
  geom_ribbon(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lme, 
                                     ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18))  +
  scale_x_continuous(limits = c(0, 12)) 

disp_plot_di <- disp_plot_di +
  facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent")) 
cv_plot_di <- cv_plot_di +
  facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
disp_plot_di_mean<- disp_plot_di_mean +
  facet_grid(cols = vars(cv_lab_mean), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
cv_plot_di_mean <- cv_plot_di_mean +
  facet_grid(cols = vars(cv_lab_mean), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent"))

plot_grid(cv_plot_di, disp_plot_di, 
          cv_plot_di_mean, disp_plot_di_mean, 
          ncol = 2, align = "h")
