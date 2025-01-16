## examine model fit 


lm_disp_int <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp + 
                     m_int*DispersalPotentialKmY*ClimVeloTKmY_spp + int,
                   data = data, 
                   start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                   algorithm = "port")

lm_disp_int_cont <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp + int,
                        data = cont, 
                        start = list(int = 0, m = 1, m_cv = 1),
                        algorithm = "port")
summary(lm_disp_int_cont)

AIC(lm_disp_int_cont)





summary(lm_disp_int)

plot(lm_disp_int)

df_disp_int

df_disp_int <- expand.grid(DispersalPotentialKmY = c(0, 1000),
                           ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp),
                                                  by = 0.1))

df_disp_int$pred_lm <- predict(lm_disp_int, se.fit = FALSE, newdata = df_disp_int)

ints_lm_disp_int <- do.call(rbind, lapply(data$ClimVeloTKmY_spp, fun_cv, mod = lm_disp_int))


data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = DispersalPotentialKmY)) +
  geom_smooth(data = df_disp_int, aes(x = ClimVeloTKmY_spp, y = pred_lm, colour = DispersalPotentialKmY,
                                      group = DispersalPotentialKmY),
              inherit.aes = FALSE, alpha = 0.1, method = "lm") +
  geom_point(alpha = 0.7, aes(shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = ClimVeloTKmY_spp, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Mean climate velocity (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Maximum\npotential\ndispersal\nrate (km/y)') +
  geom_ribbon(data = cbind(data, ints_lm_disp_int), aes(x = ClimVeloTKmY_spp, 
                                                        y = Estimate,
                                                        ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data, ints_lm_disp_int), aes(x = ClimVeloTKmY_spp, 
                                                      y = Estimate), 
            colour = "black", inherit.aes = FALSE) 


## plot dispersal 
df_disp_int <- expand.grid(DispersalPotentialKmY = seq(min(data$DispersalPotentialKmY), max(data$DispersalPotentialKmY),
                                                       by = 1),
                           ClimVeloTKmY_spp = c(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp)))

df_disp_int$pred_lm <- predict(lm_disp_int, se.fit = FALSE, newdata = df_disp_int)

ints_lm_disp_int <- do.call(rbind, lapply(data$DispersalPotentialKmY, fun_disp, mod = lm_disp_int))

data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = df_disp_int, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp,
                                      group = ClimVeloTKmY_spp),
              inherit.aes = FALSE) +
  geom_point(alpha = 0.7, aes(shape = group)) +
  # geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
  #            fill = "transparent", pch = 2,
  #            aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  # geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
  #            fill = "transparent", pch = 1,
  #            aes(x = DispersalPotentialKmY, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean climate\nvelocity (km/y)') +
  geom_ribbon(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                                        y = Estimate,
                                                        ymin = Q2.5, ymax = Q97.5), 
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = cbind(data, ints_lm_disp_int), aes(x = DispersalPotentialKmY, 
                                                      y = Estimate), 
            colour = "black", inherit.aes = FALSE) +
  scale_x_log10()

data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, colour = group)) +
  geom_histogram()

data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = DispersalPotentialKmY, colour = group)) +
  geom_point()

data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = DispersalPotentialKmY, colour = group)) +
  geom_point()

data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, colour = ID)) +
  geom_histogram() +
  facet_wrap(~group)


data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ID, scales = "free")

data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY)) +
  geom_point() +
  geom_smooth(method = "lm")

data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ID, scales = "free")

data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY)) +
  geom_point() +
  geom_smooth(method = "lm")+
  scale_x_log10()




data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ID)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) 

data %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, colour = ID)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) 



###############################################################################
####      analyzing influence of outliers / points with high leverage      ####
###############################################################################
## Cook's distance can be used to calculate influence on all parameters at once 
exps <- data

## get predictions from full model
preds_disp_full <- predict(lm_disp, newdata = data.frame(DispersalPotentialKmY = exps$DispersalPotentialKmY))
preds_climvelo_full <- predict(lm_climvelo, newdata = data.frame(ClimVeloTKmY_spp = exps$ClimVeloTKmY_spp))
preds_limrate_full <- predict(lm_limrate, newdata = data.frame(LimitingRate = exps$LimitingRate))
preds_int_full <- predict(lm_disp_int, newdata = data.frame(DispersalPotentialKmY = exps$DispersalPotentialKmY,
                                                            ClimVeloTKmY_spp = exps$ClimVeloTKmY_spp))

## calculate cook's distance
MSE_disp = mean(residuals(lm_disp)^2) ## mean squared error
MSE_climvelo = mean(residuals(lm_climvelo)^2) ## mean squared error
MSE_limrate = mean(residuals(lm_limrate)^2) ## mean squared error
MSE_int = mean(residuals(lm_disp_int)^2) 
p = 2 ## number of estimated parameters in model

## loop through points, removing one by one, fitting model, and getting prediction for point that is removed

## note: must assign new data frame to old 'lat' object, or else weird error occurs 
## so must refresh data each time 
cooks <- c() 
x = 1
while(x <= nrow(exps)) {
  data = exps
  
  ## remove point x
  data <- data[-x,]
  
  ##  fit models to data without the point:
  lm_disp_sub <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                     data = data, 
                     start = list(int = 0, m = 1),
                     algorithm = "port")
  
  lm_climvelo_sub <- nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                         data = data, 
                         start = list(int = 0, m = 1),
                         algorithm = "port")
  
  lm_limrate_sub <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                        data = data, 
                        start = list(int = 0, m = 1),
                        algorithm = "port")
  
  lm_disp_int_sub <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
                       m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
                     data = data, 
                     start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                     algorithm = "port")
  
  ## get model predictions
  preds_disp <- predict(lm_disp_sub, newdata = data.frame(DispersalPotentialKmY = exps$DispersalPotentialKmY))
  preds_climvelo <- predict(lm_climvelo_sub, newdata = data.frame(ClimVeloTKmY_spp = exps$ClimVeloTKmY_spp))
  preds_limrate <- predict(lm_limrate_sub, newdata = data.frame(LimitingRate = exps$LimitingRate))
  preds_int <- predict(lm_disp_int_sub, newdata = data.frame(DispersalPotentialKmY = exps$DispersalPotentialKmY,
                                                               ClimVeloTKmY_spp = exps$ClimVeloTKmY_spp))
  
  
  ## calculate cook' distance 
  cooks_disp <- sum((preds_disp - preds_disp_full)^2)/p*MSE_disp
  cooks_climvelo <- sum((preds_climvelo - preds_climvelo_full)^2)/p*MSE_climvelo
  cooks_limrate <- sum((preds_limrate - preds_limrate_full)^2)/p*MSE_limrate
  cooks_int <- sum((preds_int - preds_int_full)^2)/4*MSE_int
  
  ## save 
  cooks <- rbind(cooks, data.frame(cooks_disp = cooks_disp,
                                   cooks_climvelo = cooks_climvelo,
                                   cooks_limrate = cooks_limrate, 
                                   cooks_int = cooks_int))
  
  print(paste0("On obs. number: ", x))
  
  x = x + 1
}

data <- exps

hist(cooks$cooks_disp)
hist(cooks$cooks_climvelo)
hist(cooks$cooks_limrate)
hist(cooks$cooks_int)

## add cook's distance + coefficients to dataset 
exps <- cbind(exps, cooks)

## plot
exps_influential <- exps %>%
  mutate(is_bigger = cooks_limrate >= 1) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = is_bigger, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  # stat_function(colour = "grey", fun = function(x){x},
  #               linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_y_continuous(limits = c(-6, 25), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = "Cook's distance > 1", 
       shape = "") 

## run models without high influence points 
data_disp_sub <- exps[exps$cooks_disp < 1, ]
data_climvelo_sub <- exps[exps$cooks_climvelo < 1, ]
data_limrate_sub <- exps[exps$cooks_limrate < 1, ]
data_int_sub <- exps[exps$cooks_int < 1, ]

mod_lowlev_disp <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                       data = data_disp_sub, 
                       start = list(int = 0, m = 1),
                       algorithm = "port") 

mod_lowlev_climvelo <- nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                           data = data_climvelo_sub, 
                           start = list(int = 0, m = 1),
                           algorithm = "port") 

mod_lowlev_limrate <- nls(formula = ShiftKmY ~ m*LimitingRate + int,
                          data = data_limrate_sub, 
                          start = list(int = 0, m = 1),
                          algorithm = "port") 
mod_lowlev_int <- nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp +
                         m_int*ClimVeloTKmY_spp*DispersalPotentialKmY + int,
                       data = data_int_sub, 
                       start = list(int = 0, m = 1, m_cv = 1, m_int = 1),
                       algorithm = "port")

summary(mod_lowlev_int)

## make predictions from model
df_limrate_lowlev <- data.frame(LimitingRate = seq(min(data_limrate_sub$LimitingRate), 
                                                   max(data_limrate_sub$LimitingRate),
                                                   by = 0.001))

df_limrate_lowlev$pred_lm <- predict(mod_lowlev_limrate, se.fit = FALSE, newdata = df_limrate_lowlev)
df_limrate_lowlev$facet = "Without influential points"

ints_lm_limrate_lowlev <- predict2_nls(mod_lowlev_limrate, interval = "conf")

df_limrate$facet = "With influential points"
ints_lm_limrate_lowlev$facet = "Without influential points"
ints_lm_limrate$facet = "With influential points"

## plot predictions from model with and without the high leverage points 
facet <- data_disp_sub %>%
  select(-c(cooks_limrate,cooks_disp,cooks_climvelo)) %>%
  mutate(facet = "Without influential points") %>%
  rbind(., mutate(data, facet = "With influential points")) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(data, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(data, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-6, 25), expand = c(0,0.5)) +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_ribbon(data = cbind(data, ints_lm_limrate), aes(x = LimitingRate, y = fitted(lm_limrate),
                                                       ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_ribbon(data = cbind(data_limrate_sub, ints_lm_limrate_lowlev), aes(x = LimitingRate, y = fitted(mod_lowlev_limrate),
                                                                          ymin = Q2.5, ymax = Q97.5),
              fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = df_limrate, aes(x = LimitingRate, y = pred_lm), colour = "black", inherit.aes = FALSE) +
  geom_line(data = df_limrate_lowlev, aes(x = LimitingRate, y = pred_lm), 
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~facet)

ggsave(facet, path = "figures/model_results", filename = "model-predictions_influence_expansions.png", 
       width = 8, height = 4)

