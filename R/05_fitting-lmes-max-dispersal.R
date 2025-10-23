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
median(tally$n) # 1

## save data
write.csv(data, "data-processed/model-data_main.csv", row.names = FALSE)

## look at distribution of response variable 
hist(data$ShiftKmY)

## get % of species with dispersal > clim velo
round(length(which(data$DispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 62 % max disp. mean cv.
round(length(which(data$DispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 55 % max disp. q3 cv
## birds:
birds = filter(data, group == "Bird")
round(length(which(birds$DispersalPotentialKmY > birds$ClimVeloKmY_RelScale))/nrow(birds)*100, digits = 0) ## 97 % max disp. mean cv.
## plants
plants = filter(data, group == "Plant")
round(length(which(plants$DispersalPotentialKmY > plants$ClimVeloKmY_RelScale))/nrow(plants)*100, digits = 0) ## 29% max disp. mean cv.
## median
round(length(which(data$MedianDispersalPotentialKmY > data$ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 51% med disp. mean cv.
round(length(which(data$MedianDispersalPotentialKmY > data$q3ClimVeloKmY_RelScale))/nrow(data)*100, digits = 0) ## 44% med disp. q3 cv
## plants
round(length(which(data$MedianDispersalPotentialKmY > data$ClimVeloKmY_RelScale & data$group == "Plant"))/nrow(data)*100, digits = 0) ## 10% med disp.
round(length(which(data$MedianDispersalPotentialKmY > data$q3ClimVeloKmY_RelScale & data$group == "Plant"))/nrow(data)*100, digits = 0) ## 6% med disp.

## get min and max 
min(data$DispersalPotentialKmY)
max(data$DispersalPotentialKmY)


## function to sqrt transform with negative numbers:
ssqrt_trans <- scales::trans_new(name      = 'signed square root',
                                 transform = function(x) sign(x) * sqrt(abs(x)),
                                 inverse   = function(y) sign(y) * y^2,
                                 domain    = c(-Inf,Inf))
## density plot of all 3 rates
density = data %>%
  gather(key = "type", value = "measure", c(ClimVeloKmY_RelScale, DispersalPotentialKmY, ShiftKmY)) %>%
  mutate(type = factor(type, levels = c("DispersalPotentialKmY", "ClimVeloKmY_RelScale", "ShiftKmY"), 
                       ordered = T)) %>%
  ggplot(aes(x = type, y = measure, colour = group)) + 
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("Potential\ndispersal rate", 
                              "Local velocity of\nisotherm shifts", 
                              "Local range\nexpansion rate")) +
  geom_point(position = position_jitterdodge(seed = 120), 
             size = 0.1, alpha = 0.4) +
  geom_violin(fill = "transparent", scale = "width") +
  labs(x = "", y = "Rate (km/yr)", colour = "") +
  scale_y_continuous(trans = ssqrt_trans, breaks = c(0, 200, 400, 600, 800, 1000, 1200)) +
  coord_flip() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#e49e00", "#A2B06D"))

## save as new figure 3
ggsave(density, path = "figures",  filename = "figure3_distribution-of-rates.png",
       width = 7, height = 3)


## try plotting the difference between dispersal and climate velocity with broken axis
## transform y positions > 25
trans <- function(x){ pmin(x, 50) + 0.05*pmax(x - 50, 0) + ifelse(x > 50, 5, 0)}

## new y labs
y_labs = c(0, 25, 50, 50, 500, 1000)

data %>%
  mutate(diff = DispersalPotentialKmY - ClimVeloKmY_RelScale) %>%
  arrange(diff) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(diff = trans(diff)) %>%
  mutate(ymin = pmin(diff, 0), ymax = pmax(diff, 0)) %>%
  ggplot(aes(x = id, colour = ClimVeloKmY_RelScale)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax)) +
  geom_rect(ymin=50, ymax=55, xmin=-30, xmax=450, fill="white", inherit.aes = F) +
  scale_y_continuous(breaks = c(0, 25, 50, 55, 77.5, 102.5), labels = y_labs) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 50, linetype="dashed", colour = "gray") +
  geom_hline(yintercept = 55, linetype="dashed", colour = "gray") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  labs(x = "", y = "Difference between potential dispersal rate\nand velocity of isotherm shift (km/yr)",
       colour = "Local velocity\nof isotherm\nshift (km/yr)") +
  geom_hline(yintercept = 0, colour = "black") 

## make a plot 
barplot = data %>%
  mutate(is_faster = ifelse(ClimVeloKmY_RelScale > DispersalPotentialKmY, 
                            "Potential dispersal rate *faster* than<br />velocity of isotherm shift", 
                            "Potential dispersal rate *slower* than<br />velocity of isotherm shift")) %>%
  count(group, is_faster) %>% 
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = is_faster, y = prop, fill = group)) +
  geom_col(position = position_stack(), colour = "black", alpha = 0.6) +
  #geom_text(aes(label = round(100 * prop)), vjust = -0.2) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "% of total range shift estimates", fill = "", x = "") +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = ggtext::element_markdown())
# ggsave(barplot, path = "figures",  filename = "dispersal-limitation_barplot.png",
#        width = 4, height = 3)

## make a pie chart
library(ggpattern)
library(ggrepel)
pie = data %>%
  mutate(is_faster = ifelse(ClimVeloKmY_RelScale > DispersalPotentialKmY, 
                            "Potential dispersal rate *slower* than<br />velocity of isotherm shifts", 
                            "Potential dispersal rate *faster* than<br />velocity of isotherm shifts")) %>%
  mutate(is_faster = factor(is_faster, ordered = TRUE, 
                            levels = c("Potential dispersal rate *slower* than<br />velocity of isotherm shifts", 
                                       "Potential dispersal rate *faster* than<br />velocity of isotherm shifts"))) %>%
  count(group, is_faster) %>% 
  mutate(prop = n / sum(n)) %>%
  mutate(group = factor(group)) %>%
  arrange(desc(prop)) %>%
  ggplot(aes(x = "", y = prop, fill = group), colour = "black") +
  labs(y = "% of total range shift estimates", fill = "", pattern = "") +
  geom_bar(stat = "identity",  alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_col_pattern(aes(fill = group, pattern = is_faster), colour = "black", size = 0.5,
                   pattern_colour = "black",
                   pattern_fill = "black",
                   pattern_size = 0.1, 
                   pattern_density = 0.05,  alpha = 0.3) +
  theme_void() +
  scale_pattern_discrete(choices = c("none", "stripe")) +
  theme(legend.text = ggtext::element_markdown()) +
  guides(fill = "none", 
         pattern = guide_legend(override.aes = list(fill = "transparent", colour = "black"))) +
  geom_label_repel(aes(label = scales::percent(prop, accuracy = 1)), y = c(0.75, 0.3, 0.09, 0.99),
                   colour = "black", nudge_x = 0.5, min.segment.length = 10) +
  scale_fill_manual(values = c("white", "white"))

ggsave(pie, path = "figures",  filename = "figure3_pie-chart.png",
       width = 4, height = 3)

###################################
##        ALL OBSERVATIONS       ##
###################################
## fit and compete the following linear mixed effect models to all observations of range expansion
## fit each using mean and 3rd quartile of climate velocity
## add species as a random effect to account for non-independence

## 1. range expansion rate ~ potential dispersal rate
lme_disp = lme(ShiftKmY ~ DispersalPotentialKmY, 
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

## 4. range expansion rate ~ potential dispersal rate*velocity of isotherm shift 
lme_disp_int <- lme(ShiftKmY ~ DispersalPotentialKmY*ClimVeloKmY_RelScale,
                    random = ~ 1|sp_name_checked,
                    data = data,
                    method = "ML")
lme_disp_int_q3 <- lme(ShiftKmY ~ DispersalPotentialKmY*q3ClimVeloKmY_RelScale,
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

saveRDS(main_models, "data-processed/modelfits_main_allobs_lme.rds")

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
  select(Modnames, K, AICc, Delta_AICc, AICcWt, Cum.Wt, LL) %>%
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
        `Conditional R2`, `Conditional R2`, `Marginal R2`, DF, `t-value`, `p-value`, n, K, LL, AICc, everything())  %>%
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

#gtsave(table_main, path = "figures/model_results/all-observations", filename = "table_main-models_all-observations_lme.png")


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
  labs(x = "Minimum of potential dispersal rate\nand local velocity of isotherm shifts (km/yr)",
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
  labs(x = "Minimum of potential dispersal rate\nand local velocity of isotherm shifts (km/yr)",
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

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_cv_q3_lme.png", 
       width = 9.2, height = 3.2)

##save legend
temp = cv_plot_q3 +
  theme(legend.position = "right", legend.justification = "center")

legend = ggpubr::get_legend(temp) 

ggsave(legend, path = "figures/model_results/all-observations", filename = "model-predictions_cv_q3_legend.png",
       width = 1, height = 4)

## new figure 4:
fig4a = data %>%
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
       colour = 'Local velocity of\nisotherm\nshifts (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  # geom_ribbon(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lme, 
  #                                             ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
  #             fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = df_disp_int_q3_mean, aes(x = DispersalPotentialKmY, y = pred_lme),
  #           inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, unit = "cm"))

fig4b <- data %>%
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
  labs(x = "Minimum of potential dispersal rate\nand local velocity of isotherm shifts (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Local velocity of\nisotherm\nshifts (km/yr)') +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme, 
                                        ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate_q3, aes(x = LimitingRate_q3, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 12.1)) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, unit = "cm"))

fig4 = plot_grid(fig4a, fig4b,
          ncol = 2, align = "h")

ggsave(fig4, path = "figures", filename = "figure4_dispersal-mod-vs-lim-rate-mod.png", 
       width = 6.8, height = 3.2)

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

ggsave(path = "figures/model_results/all-observations", filename = "model-predictions_cv_comparison_lme.png", 
       height = 7.5, width = 8)

## calculate ratio and stats about shift versus expansion
round(length(which(data$DispersalPotentialKmY > abs(data$ShiftKmY)))/nrow(data)*100, digits = 0) ## 57
round(length(which(data$MedianDispersalPotentialKmY > abs(data$ShiftKmY)))/nrow(data)*100, digits = 0) ## 47

mean(data$DispersalPotentialKmY[data$ShiftKmY != 0]/abs(data$ShiftKmY[data$ShiftKmY != 0]), na.rm = T)
## 97.85
median(data$DispersalPotentialKmY[data$ShiftKmY != 0]/abs(data$ShiftKmY[data$ShiftKmY != 0]), na.rm = T)
## 3.75

### new figure 2:
plot_grid(cv_plot_q3, cv_plot, limrate_plot_q3, limrate_plot,
          ncol = 2, align = "h")

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

saveRDS(di_models, "data-processed/modelfits_main_diobs_lme.rds")

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
    title = "Main model set - dispersal-insufficient observations"
  ) 

gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3_lme.png")
gtsave(table_di, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_q3_lme.docx")



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

saveRDS(di_models_mean, "data-processed/modelfits_main_diobs-mean_lme.rds")

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
    title = "Main model set - dispersal-insufficient observations"
  ) 

gtsave(table_di_mean, path = "figures/model_results/dispersal-insufficient", filename = "table_di_cv_mean_lme.png")


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

ggsave(path = "figures/model_results/dispersal-insufficient", filename = "model-predictions_di_cv_comparison_lme.png", 
       width = 8, height = 7.5)

disp_plot_di <- di_data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2, alpha = 0.5, 
             aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1, alpha = 0.5, 
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
       colour = 'Local velocity\n of isotherm\nshifts (km/yr)') +
   facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme, 
                                  ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_disp, aes(x = DispersalPotentialKmY, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5)  +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 8.5)) +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 5, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "lme_disp_di")[1]],
                          "\nR2 = ", format(round(as.numeric(coefs_saved$`Conditional R2`[which(coefs_saved$Model == "lme_disp_di")[1]]), 
                                           digits = 2), nsmall = 2)), colour = "black",
           size = 3.5) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, t = 0.1, unit = "cm")) 

cv_plot_di <-  di_data %>%
  ggplot(aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, colour = q3ClimVeloKmY_RelScale)) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Plant"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2, alpha = 0.5, 
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(di_data, is.na(colour_q3), group == "Bird"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1, alpha = 0.5, 
             aes(x = q3ClimVeloKmY_RelScale, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Local velocity of isotherm shifts (km/yr)",
       y = "Local range expansion rate (km/yr)", 
       colour = 'Local velocity of\nisotherm\nshifts (km/yr)') +
  facet_grid(cols = vars(cv_lab_q3), scales = "free") +
  theme(strip.text.x = element_text(size = 11, angle = 0, colour = "black"),
        strip.background = element_rect(fill = "lightgrey", colour = "transparent")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_ribbon(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme, 
                                ymin = pred_lme - pred_lme_se, ymax = pred_lme + pred_lme_se), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_cv, aes(x = q3ClimVeloKmY_RelScale, y = pred_lme),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_shape_manual(values = c(19,17,15,18)) +
  scale_x_continuous(limits = c(0, 8.5)) +
  theme(strip.text.x = element_blank(),
        strip.background =  element_blank()) +
  annotate("text", x = 5, y = 23, hjust = 0,
           label = paste0("AICc = ", coefs_saved$AICc[which(coefs_saved$Model == "lme_cv_di")[1]],
                          "\nR2 = ", coefs_saved$`Conditional R2`[which(coefs_saved$Model == "lme_cv_di")[1]]), colour = "black",
           size = 3.5) +
  theme(plot.margin = margin(r = 0.5, l = 0.2, t = 0.5, unit = "cm"))

plot_grid(disp_plot_di, cv_plot_di,
          ncol = 2, align = "h")

ggsave(path = "figures", filename = "figure5_dispersal-slower-than-cv.png", 
       width = 6.2, height = 3)

disp_plot_di = disp_plot_di + theme(legend.position = "right")

## save legend
legend = ggpubr::get_legend(disp_plot_di) 

ggsave(legend, path = "figures", filename = "figure5_legend.png", 
       width = 1, height = 4)

## calculate some means
mean(di_data$DispersalPotentialKmY)
mean(di_data$ShiftKmY)

#########################################
##       PHYLOGENETIC CORRELATION      ##
#########################################
### check for phylogenetic correlation in residuals of the best model(s)
library(ape)
library(phylobase)
library(phylosignal)

## write out species list to use to get tree from TimeTree
write.csv(str_replace_all(unique(data$sp_name_checked), "\\_", " "), "data-processed/sp_list.csv", row.names = F)

## read in TimeTree:
tree = read.tree("data-raw/other/timetree.nwk")

## transform into phylo4 object
tree <- as(tree, "phylo4")

## make dataframe of best model residuals 
## best model:
lme_limrate_q3

resid = data.frame(resid = residuals(lme_limrate_q3))

## for species with multiple residuals, take mean 
resid = resid %>%
  mutate(sp_name_checked = names(residuals(lme_limrate_q3))) %>%
  group_by(sp_name_checked) %>%
  mutate(resid = mean(resid)) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()

## get rid of species with insufficient data 
# Larus californicus (insufficient data in TimeTree to place this taxon)
# Strix varia (insufficient data in TimeTree to place this taxon)
# Larus fuscus (insufficient data in TimeTree to place this taxon)
# Carduus crispus (insufficient data in TimeTree to place this taxon)
# Sialia mexicana (insufficient data in TimeTree to place this taxon)
# Accipiter cooperii (insufficient data in TimeTree to place this taxon)
# Quercus petraea (insufficient data in TimeTree to place this taxon)
# Primula elatior (insufficient data in TimeTree to place this taxon)
# Larus delawarensis (insufficient data in TimeTree to place this taxon)
not_in_tree = c("Larus_californicus", "Strix_varia", "Larus_fuscus", "Carduus_crispus", 
                "Sialia_mexicana", "Accipiter_cooperii", "Quercus_petraea", "Primula_elatior",
                "Larus_delawarensis")

resid = resid[which(!resid$sp_name_checked %in% not_in_tree),] 

## change species names that time tree changed
resid$sp_name_checked[which(resid$sp_name_checked == "Grus_canadensis")] <- "Antigone_canadensis"

rownames(resid) = resid$sp_name_checked

resid = select(resid, resid)

## bind tip data to tree
p4 <- phylo4d(tree, tip.data = resid)

## use Pagel's lambda to test the strength of phylogenetic autocorrelation:
physig = phyloSignal(p4, methods = "Lambda")
physig$stat #Lambda = 0.02
physig$pvalue < 0.05 #FALSE

## run same test for residuals of dispersal model for dispersal-insufficient observations
## read in TimeTree:
tree = read.tree("data-raw/other/timetree.nwk")

## get sp list 
sp = tree$tip.label[which(tree$tip.label %in% di_data$sp_name_checked)]

## prune the tree 
tree <- drop.tip(tree, tree$tip.label[-match(sp, tree$tip.label)])

## transform into phylo4 object
tree <- as(tree, "phylo4")

## make dataframe of best model residuals 
## best model:
lme_disp_di

resid = data.frame(resid = residuals(lme_disp_di))

## for species with multiple residuals, take mean 
resid = resid %>%
  mutate(sp_name_checked = names(residuals(lme_disp_di))) %>%
  group_by(sp_name_checked) %>%
  mutate(resid = mean(resid)) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()

## get rid of species with insufficient data 
# Larus californicus (insufficient data in TimeTree to place this taxon)
# Strix varia (insufficient data in TimeTree to place this taxon)
# Larus fuscus (insufficient data in TimeTree to place this taxon)
# Carduus crispus (insufficient data in TimeTree to place this taxon)
# Sialia mexicana (insufficient data in TimeTree to place this taxon)
# Accipiter cooperii (insufficient data in TimeTree to place this taxon)
# Quercus petraea (insufficient data in TimeTree to place this taxon)
# Primula elatior (insufficient data in TimeTree to place this taxon)
# Larus delawarensis (insufficient data in TimeTree to place this taxon)
not_in_tree = c("Larus_californicus", "Strix_varia", "Larus_fuscus", "Carduus_crispus", 
                "Sialia_mexicana", "Accipiter_cooperii", "Quercus_petraea", "Primula_elatior",
                "Larus_delawarensis")

resid = resid[which(!resid$sp_name_checked %in% not_in_tree),] 

## change species names that time tree changed
resid$sp_name_checked[which(resid$sp_name_checked == "Grus_canadensis")] <- "Antigone_canadensis"

rownames(resid) = resid$sp_name_checked

resid = select(resid, resid)

## bind tip data to tree
p4 <- phylo4d(tree, tip.data = resid)

## use Pagel's lambda to test the strength of phylogenetic autocorrelation:
physig = phyloSignal(p4, methods = "Lambda")
physig$stat #Lambda = 0.32
physig$pvalue < 0.05 #FALSE



