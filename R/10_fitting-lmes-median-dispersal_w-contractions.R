## fit models with median instead of maximum dispersal distance INCLUDING EXTREME CONTRACTIONS
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

## filter to leading edge shifts 
dd = filter(dd, !Param %in% c("O", "TE"))

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

## get rid of non-birds and non-plants
data = filter(dd, group %in% c("Bird", "Plant"))

## save data
write.csv(data, "data-processed/model-data_contractions_median-disp.csv", row.names = FALSE)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp_median = lme(ShiftKmY ~ MedianDispersalPotentialKmY, 
               random = ~ 1|sp_name_checked,
               data = data)

## 2. range expansion rate ~ velocity of isotherm shift (at ecologically-relevant spatial scale)
lme_cv_median = lme(ShiftKmY ~ ClimVeloKmY_RelScale, 
             random = ~ 1|sp_name_checked,
             data = data)
lme_cv_q3_median <- lme(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                 random = ~ 1|sp_name_checked,
                 data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of isotherm shift
lme_disp_cv_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY + ClimVeloKmY_RelScale,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_disp_cv_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY  + q3ClimVeloKmY_RelScale,
                      random = ~ 1|sp_name_checked,
                      data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale,
                    random = ~ 1|sp_name_checked,
                    data = data)
lme_disp_int_q3_median <- lme(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                       random = ~ 1|sp_name_checked,
                       data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of isotherm shift 
lme_limrate_median <- lme(ShiftKmY ~ LimitingRate,
                   random = ~ 1|sp_name_checked,
                   data = data)
lme_limrate_q3_median <- lme(ShiftKmY ~ LimitingRate_q3,
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

saveRDS(main_models, "data-processed/modelfits_main_allobs_median-disp_lme_w-contractions.rds")

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
    title = "Main model set - all observations - median dispersal (with contractions)"
  ) 


## plot model predictions
df_disp <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                  ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred = predictSE(lme_disp_median, se.fit = T, newdata = df_disp, level = 0)
df_disp$pred_lme <- pred$fit
df_disp$pred_lme_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(data$ClimVeloKmY_RelScale), max(data$ClimVeloKmY_RelScale),
                                               by = 0.001))
pred = predictSE(lme_cv_median, se.fit = TRUE, newdata = df_cv, level = 0)
df_cv$pred_lme <- pred$fit
df_cv$pred_lme_se <- pred$se.fit

df_cv_q3 <- data.frame(q3ClimVeloKmY_RelScale = seq(min(data$q3ClimVeloKmY_RelScale), 
                                                    max(data$q3ClimVeloKmY_RelScale),
                                                    by = 0.001))
pred = predictSE(lme_cv_q3_median, se.fit = TRUE, newdata = df_cv_q3, level = 0)
df_cv_q3$pred_lme <- pred$fit
df_cv_q3$pred_lme_se <- pred$se.fit

df_disp_cv <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                     ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_cv_median, se.fit = TRUE, newdata = df_disp_cv, level = 0)
df_disp_cv$pred_lme <- pred$fit
df_disp_cv$pred_lme_se <- pred$se.fit

df_disp_cv_q3 <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                        q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_cv_q3_median, se.fit = TRUE, newdata = df_disp_cv_q3, level = 0)
df_disp_cv_q3$pred_lme <- pred$fit
df_disp_cv_q3$pred_lme_se <- pred$se.fit


df_disp_int <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                      ClimVeloKmY_RelScale = c(min(data$ClimVeloKmY_RelScale),
                                                               max(data$ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int_median, se.fit = TRUE, newdata = df_disp_int, level = 0)
df_disp_int$pred_lme <- pred$fit
df_disp_int$pred_lme_se <- pred$se.fit

df_disp_int_q3 <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                         q3ClimVeloKmY_RelScale = c(min(data$q3ClimVeloKmY_RelScale),
                                                                    max(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int_q3_median, se.fit = TRUE, newdata = df_disp_int_q3, level = 0)
df_disp_int_q3$pred_lme <- pred$fit
df_disp_int_q3$pred_lme_se <- pred$se.fit

df_disp_int_q3_mean <- data.frame(expand.grid(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                              q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predictSE(lme_disp_int_q3_median, se.fit = TRUE, newdata = df_disp_int_q3_mean, level = 0)
df_disp_int_q3_mean$pred_lme <- pred$fit
df_disp_int_q3_mean$pred_lme_se <- pred$se.fit

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))
pred =  predictSE(lme_limrate_median, se.fit = TRUE, newdata = df_limrate, level = 0)
df_limrate$pred_lme <- pred$fit
df_limrate$pred_lme_se <-  pred$se.fit

df_limrate_q3 <- data.frame(LimitingRate_q3 = seq(min(data$LimitingRate_q3), max(data$LimitingRate_q3),
                                                  by = 0.001))
pred =  predictSE(lme_limrate_q3_median, se.fit = TRUE, newdata = df_limrate_q3, level = 0)
df_limrate_q3$pred_lme <- pred$fit
df_limrate_q3$pred_lme_se <-  pred$se.fit

data$cv_lab_mean = "Mean climate velocity"
data$cv_lab_q3 = "3rd quartile climate velocity"


## DISPERSAL 
disp_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                 ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lme),
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
  labs(x = "velocity of isotherm shift (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
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
  labs(x = "velocity of isotherm shift (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
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
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_cv, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                    ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_cv, aes(x = MedianDispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_cv")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_cv")[1]]), colour = "black",
           size = 3.5)

disp_cv_q3_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_cv_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                       ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_cv_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 700, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_cv_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_cv_q3")[1]]), colour = "black",
           size = 3.5)

## DISPERSAL AND CLIMATE - INTERACTIVE
disp_int_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_int, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                     ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se,
  #                                     group = ClimVeloKmY_RelScale), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_smooth(data = df_disp_int, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                     colour = ClimVeloKmY_RelScale, group = ClimVeloKmY_RelScale),
  #             inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_int")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_int")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_int_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                        ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se,
  #                                        group = q3ClimVeloKmY_RelScale), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_smooth(data = df_disp_int_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                        colour = q3ClimVeloKmY_RelScale, group = q3ClimVeloKmY_RelScale),
  #            inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lme_disp_int_q3")[1]],
                          "\nR2 = ", coefs$`Conditional R2`[which(coefs$Model == "lme_disp_int_q3")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_mean_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  labs(x = "Potential dispersal rate (km/yr)",
       y = "Estimated range expansion rate (km/yr)", 
       colour = 'Velocity of\nisotherm\nshift (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_int_q3_mean, aes(x = MedianDispersalPotentialKmY, y = pred_lme, 
  #                                             ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_int_q3_mean, aes(x = MedianDispersalPotentialKmY, y = pred_lme),
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
       y = "Estimated range expansion rate (km/yr)", 
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
       y = "Estimated range expansion rate (km/yr)", 
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

cv_plot_q3 <- cv_plot_q3 +
  facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent")) 
cv_plot <- cv_plot +
  facet_grid(cols = vars(cv_lab_mean), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
limrate_plot_q3<- limrate_plot_q3 +
  facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
limrate_plot <- limrate_plot +
  facet_grid(cols = vars(cv_lab_mean), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent")) 

plot_grid(cv_plot_q3, cv_plot, limrate_plot_q3, limrate_plot,
          ncol = 2, align = "h")

## make table that competes all models with all observations INCLUDING CONTRACTIONS
max_models <- readRDS("data-processed/modelfits_main_allobs_lme_w-contractions.rds")
median_models <- readRDS("data-processed/modelfits_main_allobs_median-disp_lme_w-contractions.rds")

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

## get r squared 
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
  rename("Velocity of isotherm shift" = cv_type,
         "Dispersal rate" = disp_type)

## make gt table
table_main <- coefs %>% 
  gt() %>%
  tab_header(
    title = "Main model set - all observations (with contractions)"
  ) 

gtsave(table_main, path = "figures/model_results/all-observations", filename = "tableS4_all-models_lme_w-contractions.png")
gtsave(table_main, path = "figures/model_results/all-observations", filename = "tableS4_all-models_lme_w-contractions.docx")

## make table summarizing results by model 'type' 
## calculate cumulative weight and mean rank per model 'type'
tableS5 = coefs %>%
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

tableS5 = tableS5 %>% gt()

## save table 
gtsave(tableS5, path = "figures/model_results/all-observations", filename = "tableS5_summary-w-contractions.docx")

