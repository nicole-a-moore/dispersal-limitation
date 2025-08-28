## fit models with median instead of maximum dispersal distance 
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
         LimitingRate_q3 = ifelse(MedianDispersalPotentialKmY <= q3ClimVeloKmY_RelScale,
                                     MedianDispersalPotentialKmY,
                               q3ClimVeloKmY_RelScale)) %>%
  mutate(q3ClimVeloKmY_RelScale = abs(q3ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(MedianDispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_q3 = ifelse(MedianDispersalPotentialKmY == LimitingRate_q3, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY_RelScale, NA)) %>%
  mutate(colour_q3 = ifelse(what_is_limiting_q3 == "Climate", q3ClimVeloKmY_RelScale, NA)) %>%
  mutate(group = ifelse(group == "Aves", "Bird", 
                        ifelse(group == "Mammalia", "Mammal", 
                               ifelse(group == "Plants", "Plant", 
                                      ifelse(group == "Squamata", "Squamate",
                                             NA)))))

## range contractions at the leading edge likely happen for ecological reasons other than climate (e.g., habitat loss)
## but small contractions might be expansions that are measured as contractions due to measurement error 
## to be sure that the bounding of the response variable (range shift rate) at 0 doesn't have an effect on the results, refit models keeping range contractions that are slower than 10km/y

## get rid of range contractions that are father than 1 sd from the mean shift 
data <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

## get rid of non-birds and non-plants
data = filter(data, group %in% c("Bird", "Plant"))

## save data
write.csv(data, "data-processed/model-data_main_median-disp.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear models to all observations of range expansion using median potential dispersal rate
## fit each using mean and 3rd quartile of climate velocity 

## 1. range expansion rate ~ potential dispersal rate
lm_disp_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY,
              data = data)

## 2. range expansion rate ~ velocity of climate change (at ecologically-relevant spatial scale)
lm_cv_median <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
            data = data)
lm_cv_q3_median <- lm(ShiftKmY ~ q3ClimVeloKmY_RelScale,
                data = data)

## 3. range expansion rate ~ potential dispersal rate + velocity of climate change
lm_disp_cv_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY + ClimVeloKmY_RelScale,
                 data = data)
lm_disp_cv_q3_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY + q3ClimVeloKmY_RelScale,
                     data = data)

## 4. range expansion rate ~ potential dispersal rate*velocity of climate change 
lm_disp_int_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY*ClimVeloKmY_RelScale,
                  data = data)
lm_disp_int_q3_median <- lm(ShiftKmY ~ MedianDispersalPotentialKmY*q3ClimVeloKmY_RelScale,
                      data = data)

## 5. range expansion rate ~ minimum of potential dispersal rate and velocity of climate change 
lm_limrate_median <- lm(ShiftKmY ~ LimitingRate,
                 data = data)
lm_limrate_q3_median <- lm(ShiftKmY ~ LimitingRate_q3,
                     data = data)

## plot residuals 
plot(lm_disp_median) 
plot(lm_cv_median)
plot(lm_cv_q3_median)
plot(lm_disp_cv_median)
plot(lm_disp_cv_q3_median)
plot(lm_disp_int_median)
plot(lm_disp_int_q3_median)
plot(lm_limrate_median)
plot(lm_limrate_q3_median)

hist(residuals(lm_disp_median))
hist(residuals(lm_cv_median))
hist(residuals(lm_cv_q3_median))
hist(residuals(lm_disp_cv_median))
hist(residuals(lm_disp_cv_q3_median))
hist(residuals(lm_disp_int_median))
hist(residuals(lm_disp_int_q3_median))
hist(residuals(lm_limrate_median))
hist(residuals(lm_limrate_q3_median))

## model summary 
summary(lm_disp_median) 
summary(lm_cv_median)
summary(lm_cv_q3_median)
summary(lm_disp_cv_median)
summary(lm_disp_cv_q3_median)
summary(lm_disp_int_median)
summary(lm_disp_int_q3_median)
summary(lm_limrate_median)
summary(lm_limrate_q3_median)

## save models 
main_models <- list(lm_disp_median, lm_cv_median, lm_cv_q3_median, lm_disp_cv_median, lm_disp_cv_q3_median, 
                    lm_disp_int_median, lm_disp_int_q3_median, lm_limrate_median, lm_limrate_q3_median)
names(main_models) <- c("lm_disp_median","lm_cv_median", "lm_cv_q3_median", "lm_disp_cv_median", 
                        "lm_disp_cv_q3_median", "lm_disp_int_median", "lm_disp_int_q3_median", 
                        "lm_limrate_median", "lm_limrate_q3_median")

saveRDS(main_models, "data-processed/modelfits_main_allobs_median-disp.rds")

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


## plot model predictions
df_disp <- data.frame(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY),
                                                        by = 10))
pred =  predict(lm_disp_median, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv_q3 <- data.frame(q3ClimVeloKmY_RelScale = seq(min(data$q3ClimVeloKmY_RelScale), 
                                                      max(data$q3ClimVeloKmY_RelScale),
                                                      by = 0.001))
pred = predict(lm_cv_q3_median, se.fit = TRUE, newdata = df_cv_q3)
df_cv_q3$pred_lm <- pred$fit
df_cv_q3$pred_lm_se <- pred$se.fit

df_disp_cv <- data.frame(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY),
                                                  by = 10),
                      ClimVeloKmY_RelScale = mean(data$ClimVeloKmY_RelScale))
pred =  predict(lm_disp_cv_median, se.fit = TRUE, newdata = df_disp_cv)
df_disp_cv$pred_lm <- pred$fit
df_disp_cv$pred_lm_se <- pred$se.fit

df_disp_cv_q3 <- data.frame(MedianDispersalPotentialKmY = seq(0, max(data$MedianDispersalPotentialKmY),
                                                           by = 10),
                         q3ClimVeloKmY_RelScale = mean(data$q3ClimVeloKmY_RelScale))
pred =  predict(lm_disp_cv_q3_median, se.fit = TRUE, newdata = df_disp_cv_q3)
df_disp_cv_q3$pred_lm <- pred$fit
df_disp_cv_q3$pred_lm_se <- pred$se.fit

df_disp_int_q3 <- data.frame(expand.grid(MedianDispersalPotentialKmY = 
                                                 seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                               q3ClimVeloKmY_RelScale = c(min(data$q3ClimVeloKmY_RelScale),
                                                                           max(data$q3ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_int_q3_median, se.fit = TRUE, newdata = df_disp_int_q3)
df_disp_int_q3$pred_lm <- pred$fit
df_disp_int_q3$pred_lm_se <- pred$se.fit

df_disp_int_q3_mean <- data.frame(expand.grid(MedianDispersalPotentialKmY = 
                                                 seq(0, max(data$MedianDispersalPotentialKmY), by = 10), 
                                               q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_int_q3_median, se.fit = TRUE, newdata = df_disp_int_q3_mean)
df_disp_int_q3_mean$pred_lm <- pred$fit
df_disp_int_q3_mean$pred_lm_se <- pred$se.fit

df_limrate_q3 <- data.frame(LimitingRate_q3 = seq(min(data$LimitingRate_q3), max(data$LimitingRate_q3),
                                                        by = 0.001))
pred =  predict(lm_limrate_q3_median, se.fit = TRUE, newdata = df_limrate_q3)
df_limrate_q3$pred_lm <- pred$fit
df_limrate_q3$pred_lm_se <-  pred$se.fit

## DISPERSAL 
disp_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Median potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 400, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_median")[1]]), colour = "black",
           size = 3.5) +
  scale_x_continuous(limits = c(0, 650))

## CLIMATE 
cv_plot_q3 <- data %>%
  ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean + sd\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
                                    ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_cv_q3_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_cv_q3_median")[1]]), colour = "black",
           size = 3.5)

## DISPERSAL AND CLIMATE - ADDITIVE
disp_cv_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Median potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_cv_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp_cv_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 400, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_cv_q3_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_cv_q3_median")[1]]), colour = "black",
           size = 3.5) +
  scale_x_continuous(limits = c(0, 650))

## DISPERSAL AND CLIMATE - INTERACTON
disp_int_q3_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_int_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                          ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se,
                                          group = q3ClimVeloKmY_RelScale), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_smooth(data = df_disp_int_q3, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                          colour = q3ClimVeloKmY_RelScale, group = q3ClimVeloKmY_RelScale),
              inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 375, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_int_q3_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_int_q3_median")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_mean_plot <- data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_int_q3_mean, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                               ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp_int_q3_mean, aes(x = MedianDispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 375, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_int_q3_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_int_q3_median")[1]]), colour = "black",
           size = 3.5)


## LIMITING RATE
limrate_plot_q3 <- data %>%
  ggplot(aes(x = LimitingRate_q3, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_q3, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean + sd\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lm, 
                                           ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_limrate_q3_median")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_limrate_q3_median")[1]]), colour = "black",
           size = 3.5) 

plot_grid(disp_int_q3_mean_plot, cv_plot_q3, limrate_plot_q3, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_median-disp.png", 
       width = 9, height = 3)

## make table that competes all models with all observations
max_models <- readRDS("data-processed/modelfits_main_allobs.rds")
median_models <- readRDS("data-processed/modelfits_main_allobs_median-disp.rds")

## get rid of duplicated model
median_models = median_models[which(!names(median_models) %in% c("lm_cv_q3_median", "lm_cv_median"))]

## combine
main_models = append(max_models, median_models)

## get model equations
formulas <- unlist(lapply(main_models, get_formula))

## get types
disp_type = ifelse(names(main_models) %in% c("lm_disp", "lm_disp_cv", "lm_disp_cv_q3", 
                                             "lm_disp_int", "lm_disp_int_q3",
                                             "lm_limrate", "lm_limrate_q3"), 
                   "Maximum",
                   ifelse(names(main_models) %in% c("lm_disp_median", "lm_disp_cv_median", "lm_disp_cv_q3_median",
                                                    "lm_disp_int_median", "lm_disp_int_q3_median",
                                                    "lm_limrate_q3_median", "lm_limrate_median"),
                          "Median",
                          "-"))
cv_type = ifelse(names(main_models) %in% c("lm_cv", "lm_disp_cv",  "lm_disp_int",
                                            "lm_disp_cv_median", "lm_disp_int_median",
                                           "lm_limrate", "lm_limrate_median"), 
                 "Mean",
                 ifelse(names(main_models) %in% c("lm_cv_q3", "lm_disp_cv_q3",  "lm_disp_int_q3",
                                                  "lm_disp_cv_q3_median", "lm_disp_int_q3_median",
                                                  "lm_limrate_q3", "lm_limrate_q3_median"),
                        "3rd quartile",
                        "-"))
model_type = ifelse(names(main_models) %in% c("lm_cv", "lm_cv_q3"), 
                 "Velocity of climate change",
                 ifelse(names(main_models) %in% c("lm_limrate_q3_median", "lm_limrate_q3", 
                                                  "lm_limrate", "lm_limrate_median"),
                        "Minimum rate",
                        ifelse(names(main_models) %in% c("lm_disp", "lm_disp_q3", "lm_disp_q3_median", "lm_disp_median"),
                               "Potential dispersal rate", 
                               ifelse(names(main_models) %in% c("lm_disp_cv", "lm_disp_cv_q3", 
                                                                "lm_disp_cv_q3_median", "lm_disp_cv_median"),
                                      "Potential dispersal rate and velocity of climate change (additive)",
                                      ifelse(names(main_models) %in% c("lm_disp_int", "lm_disp_int_q3", 
                                                                       "lm_disp_int_median", "lm_disp_int_q3_median"),
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

gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_all-models.png")
gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_all-models.docx")


###################################################
##        DISPERSAL-INSUFFICIENT OBSERVATIONS    ##
###################################################
## fit and compete the same linear models uisng best fit climate velocity to only observations of range expansion where dispersal < climate velocity
## q3 climate velocity was best fit
di_data <- filter(data, what_is_limiting_q3 == "Dispersal")

lm_disp_di <- lm(ShiftKmY ~ MedianDispersalPotentialKmY,
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

saveRDS(di_models, "data-processed/modelfits_main_diobs_median-disp.rds")

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

gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q90_median-disp.png")


## plot model predictions
df_disp <- data.frame(MedianDispersalPotentialKmY = seq(0, max(di_data$MedianDispersalPotentialKmY),
                                                  by = 0.001))
pred =  predict(lm_disp_di, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv <- data.frame(q3ClimVeloKmY_RelScale = seq(min(di_data$q3ClimVeloKmY_RelScale), max(di_data$q3ClimVeloKmY_RelScale),
                                                     by = 0.001))
pred = predict(lm_cv_di, se.fit = TRUE, newdata = df_cv)
df_cv$pred_lm <- pred$fit
df_cv$pred_lm_se <- pred$se.fit


## DISPERSAL 
disp_plot_di <- di_data %>%
  ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Median potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12))  +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_di")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_di")[1]]), colour = "black",
           size = 3.5) 

## CLIMATE
cv_plot_di <- di_data %>%
  ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12)) +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_cv_di")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_cv_di")[1]]), colour = "black",
           size = 3.5) 

plot_grid(disp_plot_di, cv_plot_di,
          ncol = 2, align = "h")

ggsave(path = "figures/model_results/dispersal-insufficient", filename = "model-predictions_di_cv_q3_median-disp.png", 
       width = 6, height = 2.8)


# ###################################################
# ##        DISPERSAL-SUFFICIENT OBSERVATIONS      ##
# ###################################################
# ## fit and compete the same linear models to only observations of range expansion where dispersal >climate velocity
# ds_data <- filter(data, what_is_limiting_q3 == "Climate")
# 
# lm_disp_ds <- lm(ShiftKmY ~ MedianDispersalPotentialKmY,
#                  data = ds_data)
# lm_cv_ds <- lm(ShiftKmY ~ q3ClimVeloKmY_RelScale,
#                data = ds_data)
# 
# ## plot residuals 
# plot(lm_disp_ds) 
# plot(lm_cv_ds)
# 
# hist(residuals(lm_disp_ds))
# hist(residuals(lm_cv_ds))
# 
# ## model summary 
# summary(lm_disp_ds)
# summary(lm_cv_ds)
# 
# ## save models 
# ds_models <- list(lm_disp_ds, lm_cv_ds)
# names(ds_models) <- c("lm_disp_ds", "lm_cv_ds")
# 
# saveRDS(ds_models, "data-processed/modelfits_main_dsobs_median-disp.rds")
# 
# ## get model equations
# formulas <- unlist(lapply(ds_models, get_formula))
# 
# ## get r squared 
# ds_rsq <- unlist(lapply(ds_models, FUN = function(lm) {summary(lm)$r.squared}))
# 
# ## get n 
# n <- c(rep(nrow(ds_data), 2))
# 
# ## get aic
# aic_ds <- aictab(cand.set = ds_models, modnames = names(ds_models)) %>%
#   data.frame() %>%
#   select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
#   rename("Model" = Modnames)
# 
# ## get coefs, join to r squared + se and aic table
# coefs <- lapply(ds_models, FUN = function(x) {summary(x)$coefficients})
# for (i in 1:length(coefs)) {
#   cur <- as.data.frame(coefs[[i]])
#   cur$Model = names(ds_models)[i]
#   cur$Parameter = rownames(cur)
#   cur$r_squared = ds_rsq[i]
#   cur$n = n[i]
#   cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
#   coefs[[i]] <- cur
# }
# coefs <- coefs %>%
#   bind_rows() %>%
#   left_join(., aic_ds) 
# 
# ## round everything to 2 decimal places 
# nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
# coefs[,nums] <- round(coefs[,nums], 2)
# 
# ## reorder/rename columns 
# coefs <- coefs %>%
#   rename("R2" = r_squared, "p-value" = `Pr(>|t|)`, "t-value" = `t value`, "ΔAICc" = Delta_AICc,
#          "AIC weight" = AICcWt) %>%
#   select(Model, Formula, Parameter, Estimate, `Std. Error`, `t-value`, 
#          `p-value`, `R2`, n, K, LL, AICc, everything())  %>%
#   select(-Cum.Wt, -Formula) %>%
#   arrange(`ΔAICc`)
# 
# ## replace NA with blank
# coefs <- coefs %>% 
#   mutate(across(everything(), as.character)) %>%
#   mutate(across(everything(), ~replace_na(.x, "")))
# 
# ## make gt table
# table_ds <- coefs %>% 
#   gt() %>%
#   tab_header(
#     title = "Main model set - dispersal-sufficient observations - median dispersal"
#   ) 
# 
# gtsave(table_ds, path = "figures/model_results/dispersal-sufficient", filename = "table_ds_median-disp.png")
# 
# 
# ## plot model predictions
# df_disp <- data.frame(MedianDispersalPotentialKmY = seq(0, max(ds_data$MedianDispersalPotentialKmY),
#                                                   by = 0.001))
# pred =  predict(lm_disp_ds, se.fit = TRUE, newdata = df_disp)
# df_disp$pred_lm <- pred$fit
# df_disp$pred_lm_se <- pred$se.fit
# 
# df_cv <- data.frame(q3ClimVeloKmY_RelScale = seq(min(ds_data$q3ClimVeloKmY_RelScale), max(ds_data$q3ClimVeloKmY_RelScale),
#                                                      by = 0.001))
# pred = predict(lm_cv_ds, se.fit = TRUE, newdata = df_cv)
# df_cv$pred_lm <- pred$fit
# df_cv$pred_lm_se <- pred$se.fit
# 
# 
# ## DISPERSAL 
# disp_plot_di <- ds_data %>%
#   ggplot(aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
#   geom_point(alpha = 0.7, aes(shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 0,
#              aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 5,
#              aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = MedianDispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"),
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Median potential dispersal rate (km/y)",
#        y = "Observed range expansion rate (km/y)", 
#        colour = 'Mean\nclimate\nvelocity\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm, 
#                                   ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = MedianDispersalPotentialKmY, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5)  +
#   scale_shape_manual(values = c(19,17,15,18)) +
#   theme(strip.text.x = element_blank(),
#         strip.background =  element_blank()) +
#   annotate("text", x = 400, y = 23, hjust = 0,
#            label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_ds")[1]],
#                           "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_ds")[1]]), colour = "black",
#            size = 3.5) +
#   scale_x_continuous(limits = c(0, 650))
# 
# ## CLIMATE
# cv_plot_di <- ds_data %>%
#   ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
#   geom_point(alpha = 0.7, aes(shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 0,
#              aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 5,
#              aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"),
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Velocity of climate change (km/y)",
#        y = "Observed range expansion rate (km/y)", 
#        colour = 'Mean\nclimate\nvelocity\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
#                                 ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5) +
#   scale_shape_manual(values = c(19,17,15,18)) +
#   scale_x_continuous(limits = c(0, 12.1)) +
#   theme(strip.text.x = element_blank(),
#         strip.background =  element_blank()) +
#   annotate("text", x = 7, y = 23, hjust = 0,
#            label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_cv_ds")[1]],
#                           "\nR2 = ", coefs$R2[which(coefs$Model == "lm_cv_ds")[1]]), colour = "black",
#            size = 3.5) 
# 
# plot_grid(cv_plot_di, disp_plot_di, 
#           ncol = 2, align = "h")
# 
# ggsave(path = "figures/model_results/dispersal-sufficient", filename = "model-predictions_ds_median-disp.png", 
#        width = 8.5, height = 4)

