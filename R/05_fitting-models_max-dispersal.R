## fit main models with maximum dispersal distance
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
mean(tally$n) # 1.66
median(tally$n) # 1

## save data
write.csv(data, "data-processed/model-data_main.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

## make a plot 
data %>%
  mutate(is_faster = ifelse(ClimVeloKmY_RelScale > DispersalPotentialKmY, 
                            "Dispersal-limited", "Not dispersal-limited")) %>%
  ggplot(aes(x = is_faster)) +
  geom_bar()

barplot = data %>%
  mutate(is_faster = ifelse(q3ClimVeloKmY_RelScale > DispersalPotentialKmY, 
                            "Dispersal-limited", "Not dispersal-limited")) %>%
  count(group, is_faster) %>% 
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = is_faster, y = prop, fill = group)) +
  geom_col(position = position_stack(), colour = "black", alpha = 0.6) +
  #geom_text(aes(label = round(100 * prop)), vjust = -0.2) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "% of total observations", fill = "", x = "") +
  scale_x_discrete(expand = c(0,0)) 
## 1,5,45,43,6,0.002

ggsave(barplot, path = "figures",  filename = "dispersal-limitation_barplot.png",
       width = 4, height = 3)

## get % of species with dispersal > clim velo
round(length(which(data$DispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 2) ## 61.99 % max disp. mean cv.
round(length(which(data$DispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 2) ## 55.11 % max disp. q3 cv
## birds:
birds = filter(data, group == "Bird")
round(length(which(birds$DispersalPotentialKmY > birds$ClimVeloKmY_RelScale))/nrow(birds)*100, digits = 2) ## 96.62 % max disp. mean cv.
## plants
plants = filter(data, group == "Plant")
round(length(which(plants$DispersalPotentialKmY > plants$ClimVeloKmY_RelScale))/nrow(plants)*100, digits = 2) ## 28.50 % max disp. mean cv.
## median
round(length(which(data$MedianDispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 2) ## 50.83 % med disp. mean cv.
round(length(which(data$MedianDispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 2) ## 43.71 % med disp. q3 cv

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

saveRDS(main_models, "data-processed/modelfits_main_allobs.rds")

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


## plot model predictions
df_disp <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                      ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred =  predict(lm_disp, se.fit = TRUE, newdata = df_disp)
df_disp$pred_lm <- pred$fit
df_disp$pred_lm_se <- pred$se.fit

df_cv <- data.frame(ClimVeloKmY_RelScale = seq(min(data$ClimVeloKmY_RelScale), max(data$ClimVeloKmY_RelScale),
                                               by = 0.001))
pred = predict(lm_cv, se.fit = TRUE, newdata = df_cv)
df_cv$pred_lm <- pred$fit
df_cv$pred_lm_se <- pred$se.fit

df_cv_q3 <- data.frame(q3ClimVeloKmY_RelScale = seq(min(data$q3ClimVeloKmY_RelScale), 
                                                      max(data$q3ClimVeloKmY_RelScale),
                                                      by = 0.001))
pred = predict(lm_cv_q3, se.fit = TRUE, newdata = df_cv_q3)
df_cv_q3$pred_lm <- pred$fit
df_cv_q3$pred_lm_se <- pred$se.fit

df_disp_cv <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                      ClimVeloKmY_RelScale = c(mean(data$ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_cv, se.fit = TRUE, newdata = df_disp_cv)
df_disp_cv$pred_lm <- pred$fit
df_disp_cv$pred_lm_se <- pred$se.fit

df_disp_cv_q3 <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                        q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_cv_q3, se.fit = TRUE, newdata = df_disp_cv_q3)
df_disp_cv_q3$pred_lm <- pred$fit
df_disp_cv_q3$pred_lm_se <- pred$se.fit


df_disp_int <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                     ClimVeloKmY_RelScale = c(min(data$ClimVeloKmY_RelScale),
                                                              max(data$ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_int, se.fit = TRUE, newdata = df_disp_int)
df_disp_int$pred_lm <- pred$fit
df_disp_int$pred_lm_se <- pred$se.fit

df_disp_int_q3 <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                         q3ClimVeloKmY_RelScale = c(min(data$q3ClimVeloKmY_RelScale),
                                                                  max(data$q3ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_int_q3, se.fit = TRUE, newdata = df_disp_int_q3)
df_disp_int_q3$pred_lm <- pred$fit
df_disp_int_q3$pred_lm_se <- pred$se.fit

df_disp_int_q3_mean <- data.frame(expand.grid(DispersalPotentialKmY = seq(0, max(data$DispersalPotentialKmY), by = 10), 
                                          q3ClimVeloKmY_RelScale = c(mean(data$q3ClimVeloKmY_RelScale))))
pred =  predict(lm_disp_int_q3, se.fit = TRUE, newdata = df_disp_int_q3_mean)
df_disp_int_q3_mean$pred_lm <- pred$fit
df_disp_int_q3_mean$pred_lm_se <- pred$se.fit

df_limrate <- data.frame(LimitingRate = seq(min(data$LimitingRate), max(data$LimitingRate),
                                            by = 0.001))
pred =  predict(lm_limrate, se.fit = TRUE, newdata = df_limrate)
df_limrate$pred_lm <- pred$fit
df_limrate$pred_lm_se <-  pred$se.fit

df_limrate_q3 <- data.frame(LimitingRate_q3 = seq(min(data$LimitingRate_q3), max(data$LimitingRate_q3),
                                            by = 0.001))
pred =  predict(lm_limrate_q3, se.fit = TRUE, newdata = df_limrate_q3)
df_limrate_q3$pred_lm <- pred$fit
df_limrate_q3$pred_lm_se <-  pred$se.fit

data$cv_lab_mean = "Mean climate velocity"
data$cv_lab_q3 = "3rd quartile climate velocity"


## DISPERSAL 
disp_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp")[1]]), colour = "black",
           size = 3.5)

## CLIMATE
cv_plot <- data %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_cv")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_cv")[1]]), colour = "black",
           size = 3.5)

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
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
                                    ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv_q3, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_cv_q3")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_cv_q3")[1]]), colour = "black",
           size = 3.5)


## DISPERSAL AND CLIMATE - ADDITIVE
disp_cv_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_cv, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp_cv, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_cv")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_cv")[1]]), colour = "black",
           size = 3.5)

disp_cv_q3_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_cv_q3, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp_cv_q3, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 700, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_cv_q3")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_cv_q3")[1]]), colour = "black",
           size = 3.5)

## DISPERSAL AND CLIMATE - INTERACTIVE
disp_int_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                     ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se,
                                     group = ClimVeloKmY_RelScale), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_smooth(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                      colour = ClimVeloKmY_RelScale, group = ClimVeloKmY_RelScale),
            inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_int")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_int")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_int_q3, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                      ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se,
                                      group = q3ClimVeloKmY_RelScale), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_smooth(data = df_disp_int_q3, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                      colour = q3ClimVeloKmY_RelScale, group = q3ClimVeloKmY_RelScale),
              inherit.aes = FALSE, alpha = 0.5, method = "lm")  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 750, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_int_q3")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_int_q3")[1]]), colour = "black",
           size = 3.5)

disp_int_q3_mean_plot <- data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                          ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lm),
              inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 650, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_disp_int_q3")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_disp_int_q3")[1]]), colour = "black",
           size = 3.5)

## LIMITING RATE
limrate_plot <- data %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Minimum of potential dispersal rate\nand velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate, aes(x = LimitingRate, y = pred_lm, 
                                     ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_limrate")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_limrate")[1]]), colour = "black",
           size = 3.5) 

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
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lm, 
                                     ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs$AICc[which(coefs$Model == "lm_limrate_q3")[1]],
                          "\nR2 = ", coefs$R2[which(coefs$Model == "lm_limrate_q3")[1]]), colour = "black",
           size = 3.5) +
  scale_x_continuous(limits = c(0, 12.1)) 


plot_grid(disp_int_q3_mean_plot, cv_plot_q3, limrate_plot_q3, 
          ncol = 3, align = "h")

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_cv_q3.png", 
       width = 9.2, height = 3.2)

##save legend
temp = cv_plot_q3 +
  theme(legend.position = "right")

legend = ggpubr::get_legend(temp) 

ggsave(legend, path = "figures/model_results/all-observations", filename = "model-predictions_cv_q3_legend.png",
       width = 1, height = 4)

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

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_cv_comparison.png", 
       height = 7.5, width = 8)


  

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

saveRDS(di_models, "data-processed/modelfits_main_diobs.rds")

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

gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3.png")
gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3.docx")



## plot model predictions
df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(di_data$DispersalPotentialKmY),
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
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12)) 

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
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12)) 

coefs_saved = coefs

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

saveRDS(di_models_mean, "data-processed/modelfits_main_diobs-mean.rds")

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

gtsave(table_di_mean, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_mean.png")


## plot model predictions
df_disp_mean <- data.frame(DispersalPotentialKmY = seq(0, max(di_data_mean$DispersalPotentialKmY),
                                                  by = 0.001))
pred =  predict(lm_disp_di_mean, se.fit = TRUE, newdata = df_disp_mean)
df_disp_mean$pred_lm <- pred$fit
df_disp_mean$pred_lm_se <- pred$se.fit

df_cv_mean <- data.frame(ClimVeloKmY_RelScale = seq(min(di_data_mean$ClimVeloKmY_RelScale), 
                                                          max(di_data_mean$ClimVeloKmY_RelScale),
                                                     by = 0.001))
pred = predict(lm_cv_di_mean, se.fit = TRUE, newdata = df_cv_mean)
df_cv_mean$pred_lm <- pred$fit
df_cv_mean$pred_lm_se <- pred$se.fit


## DISPERSAL 
disp_plot_di_mean <- di_data_mean %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
                                  ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12))

## CLIMATE
cv_plot_di_mean <- di_data_mean %>%
  ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 0,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 5,
             aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data_mean, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
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
  labs(x = "Velocity of climate change (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = 'Velocity of\nclimate\nchange\n(km/y)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
                                ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
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

ggsave(path = "figures/model_results/dispersal-insufficient", filename = "model-predictions_di_cv_comparison.png", 
       width = 8, height = 7.5)

disp_plot_di <- disp_plot_di +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "Potential dispersal rate")[1]],
                          "\nR2 = ", coefs_saved$R2[which(coefs_saved$Model == "Potential dispersal rate")[1]]), colour = "black",
           size = 3.5) 
cv_plot_di <- cv_plot_di +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 7, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "Velocity of climate change")[1]],
                          "\nR2 = ", coefs_saved$R2[which(coefs_saved$Model == "Velocity of climate change")[1]]), colour = "black",
           size = 3.5) 

plot_grid(disp_plot_di, cv_plot_di,
          ncol = 2, align = "h")

ggsave(path = "figures/model_results/dispersal-insufficient", filename = "model-predictions_di_cv_q3.png", 
       width = 6, height = 2.8)

disp_plot_di + theme(legend.position = "right")


## calculate some means
mean(di_data$DispersalPotentialKmY)
mean(di_data$ShiftKmY)


# ###################################################
# ##        DISPERSAL-SUFFICIENT OBSERVATIONS      ##
# ###################################################
# ## fit and compete the same linear models to only observations of range expansion where dispersal >climate velocity
# ds_data <- filter(data, what_is_limiting_q3 == "Climate")
# 
# lm_disp_ds <- lm(ShiftKmY ~ DispersalPotentialKmY,
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
# saveRDS(ds_models, "data-processed/modelfits_main_dsobs.rds")
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
#     title = "Main model set - dispersal-sufficient observations"
#   ) 
# 
# gtsave(table_ds, path = "figures/model_results/dispersal-sufficient", filename = "table_ds_cv_q3.png")
# 
# 
# ## plot model predictions
# df_disp <- data.frame(DispersalPotentialKmY = seq(0, max(ds_data$DispersalPotentialKmY),
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
# disp_plot_ds<- ds_data %>%
#   ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
#   geom_point(alpha = 0.7, aes(shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 0,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Squamate"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 5,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data, is.na(colour_q3), group == "Mammal"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"),
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Potential dispersal rate (km/y)",
#        y = "Observed range expansion rate (km/y)", 
#        colour = 'Velocity of\nclimate\nchange\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
#                                   ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5)  +
#   scale_shape_manual(values = c(19,17,15,18)) 
# 
# ## CLIMATE
# cv_plot_ds<- ds_data %>%
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
#        colour = 'Velocity of\nclimate\nchange\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm, 
#                                 ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5) +
#   scale_shape_manual(values = c(19,17,15,18)) +
#   scale_x_continuous(limits = c(0, 12.5)) 
# 
# coefs_saved = coefs
# 
# ## repeat with mean cv
# ds_data_mean <- filter(data, what_is_limiting == "Climate")
# 
# lm_disp_ds_mean <- lm(ShiftKmY ~ DispersalPotentialKmY,
#                       data = ds_data_mean)
# lm_cv_ds_mean <- lm(ShiftKmY ~ ClimVeloKmY_RelScale,
#                     data = ds_data_mean)
# 
# ## plot residuals 
# plot(lm_disp_ds_mean)
# plot(lm_cv_ds_mean)
# 
# hist(residuals(lm_disp_ds_mean))
# hist(residuals(lm_cv_ds_mean))
# 
# ## model summary 
# summary(lm_disp_ds_mean)
# summary(lm_cv_ds_mean)
# 
# ## save models 
# ds_models_mean <- list(lm_disp_ds_mean, lm_cv_ds_mean)
# names(ds_models_mean) <- c("lm_disp_ds_mean", "lm_cv_ds_mean")
# 
# saveRDS(ds_models_mean, "data-processed/modelfits_main_dsobs-mean.rds")
# 
# ## get model equations
# formulas <- unlist(lapply(ds_models_mean, get_formula))
# 
# ## get r squared 
# ds_rsq <- unlist(lapply(ds_models_mean, FUN = function(lm) {summary(lm)$r.squared}))
# 
# ## get n 
# n <- c(rep(nrow(ds_data_mean), 2))
# 
# ## get aic
# aic_ds <- aictab(cand.set = ds_models_mean, modnames = names(ds_models_mean)) %>%
#   data.frame() %>%
#   select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
#   rename("Model" = Modnames)
# 
# ## get coefs, join to r squared + se and aic table
# coefs <- lapply(ds_models_mean, FUN = function(x) {summary(x)$coefficients})
# for (i in 1:length(coefs)) {
#   cur <- as.data.frame(coefs[[i]])
#   cur$Model = names(ds_models_mean)[i]
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
# table_ds_mean <- coefs %>% 
#   gt() %>%
#   tab_header(
#     title = "Main model set - dispersal-sufficient observations"
#   ) 
# 
# gtsave(table_ds_mean, path = "figures/model_results/dispersal-sufficient", filename = "table_ds_cv_mean.png")
# 
# 
# ## plot model predictions
# df_disp_mean <- data.frame(DispersalPotentialKmY = seq(0, max(ds_data_mean$DispersalPotentialKmY),
#                                                        by = 0.001))
# pred =  predict(lm_disp_ds_mean, se.fit = TRUE, newdata = df_disp_mean)
# df_disp_mean$pred_lm <- pred$fit
# df_disp_mean$pred_lm_se <- pred$se.fit
# 
# df_cv_mean <- data.frame(ClimVeloKmY_RelScale = seq(min(ds_data_mean$ClimVeloKmY_RelScale), 
#                                                     max(ds_data_mean$ClimVeloKmY_RelScale),
#                                                     by = 0.001))
# pred = predict(lm_cv_ds_mean, se.fit = TRUE, newdata = df_cv_mean)
# df_cv_mean$pred_lm <- pred$fit
# df_cv_mean$pred_lm_se <- pred$se.fit
# 
# 
# ## DISPERSAL 
# disp_plot_ds_mean <- ds_data_mean %>%
#   ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
#   geom_point(alpha = 0.7, aes(shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 0,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 5,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
#   theme_bw() +
#   stat_function(colour = "grey", fun = function(x){x},
#                 linetype = "dashed") + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(2 , "lines"),
#         legend.position = "none") +
#   scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
#   labs(x = "Potential dispersal rate (km/y)",
#        y = "Observed range expansion rate (km/y)", 
#        colour = 'Velocity of\nclimate\nchange\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm, 
#                                   ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5)  +
#   scale_shape_manual(values = c(19,17,15,18))
# 
# ## CLIMATE
# cv_plot_ds_mean <- ds_data_mean %>%
#   ggplot(aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, colour = ClimVeloKmY_RelScale)) +
#   geom_point(alpha = 0.7, aes(shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Plant"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 0,
#              aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Bird"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 1,
#              aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Squamate"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 5,
#              aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
#   geom_point(data = filter(ds_data_mean, is.na(colour), group == "Mammal"), colour = "black", inherit.aes = FALSE,
#              fill = "transparent", pch = 2,
#              aes(x = ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
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
#        colour = 'Velocity of\nclimate\nchange\n(km/y)') +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_ribbon(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lm, 
#                                      ymin = pred_lm - pred_lm_se, ymax = pred_lm + pred_lm_se), 
#               fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
#   geom_line(data = df_cv_mean, aes(x = ClimVeloKmY_RelScale, y = pred_lm),
#             inherit.aes = FALSE, alpha = 0.5) +
#   scale_shape_manual(values = c(19,17,15,18))  +
#   scale_x_continuous(limits = c(0, 12.5)) 
# 
# disp_plot_ds<- disp_plot_ds +
#   facet_grid(cols = vars(cv_lab_q3), scales = "free") +
#   theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
#         strip.background = element_rect(fill = "lightgrey", colour = "transparent")) 
# cv_plot_ds <- cv_plot_ds +
#   facet_grid(cols = vars(cv_lab_q3), scales = "free") +
#   theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
#         strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
# disp_plot_ds_mean <- disp_plot_ds_mean +
#   facet_grid(cols = vars(cv_lab_mean), scales = "free") +
#   theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
#         strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
# cv_plot_ds_mean <- cv_plot_ds_mean +
#   facet_grid(cols = vars(cv_lab_mean), scales = "free") +
#   theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
#         strip.background = element_rect(fill = "lightgrey", colour = "transparent"))
# 
# plot_grid(cv_plot_ds, disp_plot_ds, 
#           cv_plot_ds_mean, disp_plot_ds_mean, 
#           ncol = 2, align = "h")
# 
# ggsave(path = "figures/model_results/dispersal-sufficient", filename = "model-predictions_ds_cv_comparison.png", 
#        width = 8, height = 7.5)
# 
# disp_plot_ds<- disp_plot_ds+
#   theme(strip.text.x = element_blank(),
#         strip.background =  element_blank()) +
#   annotate("text", x = 800, y = 23, hjust = 0,
#            label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "lm_disp_ds")[1]],
#                           "\nR2 = ", coefs_saved$R2[which(coefs_saved$Model == "lm_disp_ds")[1]]), colour = "black",
#            size = 3.5) 
# cv_plot_ds<- cv_plot_ds+
#   theme(strip.text.x = element_blank(),
#         strip.background =  element_blank()) +
#   annotate("text", x = 7, y = 23, hjust = 0,
#            label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "lm_cv_ds")[1]],
#                           "\nR2 = ", coefs_saved$R2[which(coefs_saved$Model == "lm_cv_ds")[1]]), colour = "black",
#            size = 3.5) 
# 
# plot_grid(cv_plot_ds, disp_plot_ds, 
#           ncol = 2, align = "h")
# 
# ggsave(path = "figures/model_results/dispersal-sufficient", filename = "model-predictions_ds_cv_q3.png", 
#        width = 8.5, height = 4)
# 
# 
