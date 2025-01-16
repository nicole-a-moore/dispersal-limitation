## asymptoic model
library(tidyverse)
library(qpcR)

data <- read.csv("data-processed/model-data_expansions.csv")

## FOR EACH STUDY:
## - compare linear model of range shift ~ dispersal to saturating model of range shift ~ dispersal
## - expect linear model to be better fit (not sampling enough parameter space)
## - expect direction of linear effect to vary across studies 
## - compare linear model of range shift ~ body size to saturating model of range shift ~ body size
## - expect linear model to be better fit 
## - expect direction to vary across studies

## fit all the models for each study
studies <- data %>%
  group_by(ID) %>%
  group_split()

studies = studies[-which(lapply(studies, FUN = nrow) < 5)]
## get rid of studies with less than 5 observations

## fit the linear models to each study
linear_mods <- purrr::map(.x = studies,
          .f = ~nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                    data = .x, 
                    start = list(m = 1, int = 0),
                    algorithm = "port"))

## fit the asymptotic models to each study
asym_mods <- purrr::map(.x = studies,
                          .f = ~nls(formula = ShiftKmY ~ L*ClimVeloTKmY_spp*(1-exp(1)^(-1*(2/ClimVeloTKmY_spp)*(DispersalPotentialKmY))),
                                    data = .x, 
                                    start = list(L = 3),
                                    algorithm = "port"))

linear_mods[[1]]
asym_mods[[1]]

## get parameters from each model
coefs_linear <- purrr::map(.x = linear_mods, 
                    .f = ~coefficients(.x)) %>%
  bind_rows()

coefs_asym <- purrr::map(.x = asym_mods, 
                    .f = ~coefficients(.x)) %>%
  bind_rows()

coefs_asym ## predicted to be close to 1
hist(coefs_asym$L) ## asymptote (cv)

coefs_linear$m
hist(coefs_linear$m) ## slope
## mixed
hist(coefs_linear$int) ## slope

## get info about model fit 
aic_linear <- purrr::map(.x = linear_mods,
                         .f = ~calculate_aic(residuals(.x),
                                             3, 
                                             length(predict(.x)))) %>%
  unlist()


aic_asym <- purrr::map(.x = asym_mods,
                         .f = ~calculate_aic(residuals(.x),
                                             2, 
                                             length(predict(.x)))) %>%
  unlist()

which(aic_linear - aic_asym > 2) ## some linear models are significantly better fit than asym 
which(aic_linear - aic_asym < -2) ## 1 asym models are significantly better fit than linear

get_rsqu <- function(data, type) {
  if(type == "lm") {
    mod <-  nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + int,
                data = data, 
                start = list(m = 1, int = 0),
                algorithm = "port")
  }
  else if(type == "asym") {
    mod <- nls(formula = ShiftKmY ~ L*ClimVeloTKmY_study*(1-exp(1)^(-1*(2/ClimVeloTKmY_study)*(DispersalPotentialKmY))),
               data = data, 
               start = list(L = 3),
               algorithm = "port")
  }

  return(summary(lm(ShiftKmY ~ predict(mod), data = data))$r.squared)
}
    

r2_linear = purrr::map(.x = studies,
                       .f = ~get_rsqu(.x, type = "lm")) %>%
  unlist()
r2_asym = purrr::map(.x = studies,
                       .f = ~get_rsqu(.x, type = "asym")) %>%
  unlist()

r2_linear - r2_asym

## make predictions
new_dfs <- purrr::map(.x = studies,
                      .f = ~expand.grid(DispersalPotentialKmY = seq(min(.x$DispersalPotentialKmY), 
                                                                    max(.x$DispersalPotentialKmY),
                                                                    by = 0.01),
                                        ClimVeloTKmY_spp = unique(.x$ClimVeloTKmY_study))) 

for(i in 1:length(linear_mods)) {
  new_df <- new_dfs[[i]]
  study = studies[[i]]
  new_df$ID = unique(study$ID)
  new_df$pred_lm <- predict(linear_mods[[i]], se.fit = FALSE, newdata = new_df) 
  new_df$pred_asym <- predict(asym_mods[[i]], se.fit = FALSE, newdata = new_df) 
  new_dfs[[i]] <- new_df
}
new_dfs <- bind_rows(new_dfs)

## plot predictions
data %>%
  filter(ID %in% new_dfs$ID) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp,
                                    group = ClimVeloTKmY_spp),
            inherit.aes = FALSE) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_asym,
                                        group = ClimVeloTKmY_spp), colour = "black",
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
  facet_wrap(~ID, scales = "free")

## plot predictions
data %>%
  filter(ID %in% new_dfs$ID) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp,
                                group = ID),
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
       colour = 'Mean climate\nvelocity (km/y)') 

data %>%
  filter(ID %in% new_dfs$ID) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_asym, colour = ClimVeloTKmY_spp,
                                group = ID),
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
       colour = 'Mean climate\nvelocity (km/y)') 

## zoom:
data %>%
  filter(ID %in% new_dfs$ID) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_lm, colour = ClimVeloTKmY_spp,
                                group = ID),
            inherit.aes = FALSE) +
  geom_line(data = new_dfs, aes(x = DispersalPotentialKmY, y = pred_asym, colour = ClimVeloTKmY_spp,
                                group = ID),
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
  scale_y_continuous(limits = c(-6, 20), expand = c(0,0.5)) +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean climate\nvelocity (km/y)') +
  scale_x_continuous(limits = c(0, 30))




## ACROSS ALL STUDIES
## - compare linear model of range shift ~ dispersal + study-level cv to saturating model of range shift ~ dispersal
## - expect saturating model to be a better fit 

## - compare saturating model to limiting rate model + model of cv alone
## - exepct limiting rate model to be a beter fit
## - if saturating model + limiting rate model come out on top, good sign 


## fit linear model
linear <-  nls(formula = ShiftKmY ~ m*DispersalPotentialKmY + m_cv*ClimVeloTKmY_spp + int,
             data = data, 
             start = list(m = 1, m_cv = 1, int = 0),
             algorithm = "port")
summary(linear)
AICc(linear)

## fit asymptotic model 
asym <-  nls(formula = ShiftKmY ~ L*ClimVeloTKmY_spp*(1-exp(1)^(-1*(2/ClimVeloTKmY_spp)*(DispersalPotentialKmY))),
                   data = data, 
                   start = list(L = 3),
                   algorithm = "port")

summary(asym)
AICc(asym)

## fit limiting rate model 
limrate <-  nls(formula = ShiftKmY ~ m*LimitingRate + int,
             data = data, 
             start = list(m = 1, int = 0),
             algorithm = "port")
summary(limrate)
AICc(limrate)

## fit climate velocity model 
climvel <-  nls(formula = ShiftKmY ~ m*ClimVeloTKmY_spp + int,
                data = data, 
                start = list(m = 1, int = 0),
                algorithm = "port")
summary(climvel)
AICc(climvel)

r_squ(linear, data)
r_squ(asym, data)
r_squ(climvel, data)
r_squ(limrate, data)

AICc(linear)
AICc(asym)
AICc(climvel)
AICc(limrate)

summary(linear)
summary(asym)
summary(climvel)
summary(limrate)


## 
vec = c(AICc(limrate), AICc(asym), AICc(linear), AICc(climvel))
vec - min(vec)

## lim rate significantly better than all
## aic of lim rate model is best
## aic of asym model is worst


## plot
df_disp_int <- expand.grid(DispersalPotentialKmY = seq(min(data$DispersalPotentialKmY), max(data$DispersalPotentialKmY),
                                                       by = 1),
                           ClimVeloTKmY_spp = seq(min(data$ClimVeloTKmY_spp), max(data$ClimVeloTKmY_spp), by = 1))

df_disp_int$pred_lm <- predict(linear, se.fit = FALSE, newdata = df_disp_int)

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
       colour = 'Mean climate\nvelocity (km/y)') 


## zoom: 

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
  scale_x_continuous(limits = c(0, 30))





## trash 
x = seq(from = 0, to = 5, by = 0.1)
y = 1.1*2*(1-exp(1)^(-1*(1/2)*(x)))
plot(x = x, y = y)
points(x=c(1:100), y=c(1:100))

x = seq(from = 0, to = 100, by = 0.1)
y = 1.1*50*(1-exp(1)^(-1*(1/50)*(x)))
plot(x = x, y = y)
points(x=c(1:100), y=c(1:100))



preds <- data.frame(x = rep(1:30, 3),
                    cv = rep(c(1,5,10), each = 30),
                    y = c(1*(1-exp(1)^(-1*(2/1)*(1:30))),
                          5*(1-exp(1)^(-1*(2/5)*(1:30))),
                          10*(1-exp(1)^(-1*(2/10)*(1:30)))))

preds %>%
  ggplot(aes(x = x, y = y, colour = cv, group = cv)) + geom_line()


unique(new_dfs$ID)

asym_mods[[10]]
View(studies[[10]])

new_dfs %>%
  filter(ID == "A181_P1") %>%
  ggplot(aes(x = DispersalPotentialKmY, y = pred_lm)) +
  geom_point() +
  geom_point(aes(y = pred_asym), colour = "red")
  
new_dfs %>%
  filter(ID == "A181_P1") %>%
  ggplot(aes(x = DispersalPotentialKmY, y = pred_lm)) +
  geom_point() +
  geom_point(aes(y = pred_asym), colour = "red")

preds <- predict(asym_mods[[10]], newdata = data.frame(DispersalPotentialKmY = seq(0, 
                                                                           300,
                                                                           by = 0.01),
                                               ClimVeloTKmY_spp = unique(studies[[10]]$ClimVeloTKmY_study)))


new_dfs %>%
  filter(ID == "A181_P1") %>%
  ggplot(aes(x = DispersalPotentialKmY, y = pred_lm)) +
  geom_point() +
  geom_point(aes(y = pred_asym), colour = "red") +
  geom_point(data = df)




df = data.frame(pred_lm = preds, DispersalPotentialKmY = seq(0,300,by = 0.01))

## issue: can't fit asymptotic model within studies since species not sampled across both parts of curve in each study
## to make a point, could count # species sampled before and after break?

