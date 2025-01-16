

## for species whose dispersal < rate of climate change, dispersal explains intraspecific differences in the responses of species better than climate velocity 
## i.e., dispersal rate explains why species respond to the same/similar climate velocities differently 
data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  filter(LimitingRate == DispersalPotentialKmY) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ID, group = ID)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  geom_smooth(method = 'lm',  se=FALSE,
              aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ID),
              colour = "black", size = 0.5) +
  scale_x_log10() +
  scale_y_log10()

data %>%
  filter(LimitingRate == DispersalPotentialKmY) %>%
  ggplot(aes(x = ClimVeloTKmY_spp , y = ShiftKmY, colour = ID, group = ID)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  geom_smooth(method = 'lm',  se=FALSE,
              aes(x = ClimVeloTKmY_spp , y = ShiftKmY, colour = ID),
              colour = "black", size = 0.5) 

data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  filter(LimitingRate != DispersalPotentialKmY) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) +
  scale_x_log10()

data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp_max, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp_max),
              colour = "black", size = 0.5) 

data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp + sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Mimimum of potential dispersal rate and\nmaximum rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) 

data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  ggplot(aes(x = ClimVeloTKmY_spp+ sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp + sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = ClimVeloTKmY_spp+ sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) 

data %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp),
              colour = "black", size = 0.5) +
  scale_x_continuous(limits = c(0, 10))

data %>%
  filter(LimitingRate == DispersalPotentialKmY) %>%
  ggplot(aes(x = ClimVeloTKmY_spp, y = ShiftKmY, 
             colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Mean rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = ClimVeloTKmY_spp, y = ShiftKmY, colour = ClimVeloTKmY_spp),
              colour = "black", size = 0.5) +
  scale_x_continuous(limits = c(0, 10))

data %>%
  filter(LimitingRate == DispersalPotentialKmY) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand mean rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp),
              colour = "black", size = 0.5)

data %>%
  filter(LimitingRate != DispersalPotentialKmY) %>%
  ggplot(aes(x = ClimVeloTKmY_spp_max, y = ShiftKmY, colour = ClimVeloTKmY_spp_max, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum rate of climate change (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = ClimVeloTKmY_spp_max, y = ShiftKmY, colour = ClimVeloTKmY_spp_max),
              colour = "black", size = 0.5)
data %>%
  filter(LimitingRate != DispersalPotentialKmY) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp_max, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp_max),
              colour = "black", size = 0.5)


### greater lags?
data %>%
  mutate(shift_lag = ClimVeloTKmY_spp + sd_cv_sppspecific - ShiftKmY) %>%
  ggplot(aes(x = what_is_limiting, y = shift_lag)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") 

data %>%
  mutate(is_lag = ClimVeloTKmY_spp + sd_cv_sppspecific > ShiftKmY) %>%
  ggplot(aes(x = what_is_limiting, colour= is_lag)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") 
## greater proportion of species that lag behind climate velocity 

data %>%
  mutate(shift_lag = ClimVeloTKmY_spp + sd_cv_sppspecific - ShiftKmY) %>%
  group_by(what_is_limiting) %>%
  summarise(mean_lag = mean(shift_lag))



data$sd_cv_sppspecific



data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  filter(LimitingRate !=  ClimVeloTKmY_spp+ sd_cv_sppspecific) %>%
  ggplot(aes(x = ClimVeloTKmY_spp+ sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = ClimVeloTKmY_spp+ sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) +
  scale_x_continuous(limits = c(0, 10))


data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  filter(LimitingRate !=  ClimVeloTKmY_spp+ sd_cv_sppspecific) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) +
  scale_x_continuous(limits = c(0, 10))



data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) 

data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  ggplot(aes(x = ClimVeloTKmY_spp + sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  geom_smooth(method = 'lm', inherit.aes = F, 
              aes(x = ClimVeloTKmY_spp + sd_cv_sppspecific, y = ShiftKmY, colour = ClimVeloTKmY_spp+ sd_cv_sppspecific),
              colour = "black", size = 0.5) 


data$ClimVeloTKmY_spp_max = data$ClimVeloTKmY_spp + data$sd_cv_sppspecific

## compare shift ~ mean vs max clim change 
lm_mean = lm(data = data,
             ShiftKmY ~ ClimVeloTKmY_spp)

lm_max = lm(data = data,
             ShiftKmY ~ ClimVeloTKmY_spp_max)

AIC(lm_mean, lm_max)
AIC(lm_mean) - AIC(lm_max)

## very similar AIC - mean not necessarily a better model than max

## compare shift ~ limiting rate using mean versus max climate change 
data <- data %>%
  mutate(LimitingRate2 = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp_max,
                               DispersalPotentialKmY, ClimVeloTKmY_spp_max)) 

lm_mean_limrate = lm(data = data,
             ShiftKmY ~ LimitingRate)

lm_max_limrate = lm(data = data,
            ShiftKmY ~ LimitingRate2)

AIC(lm_mean_limrate, lm_max_limrate)
AIC(lm_mean_limrate) - AIC(lm_max_limrate)

lm_disp = lm(data = data,
             ShiftKmY ~ DispersalPotentialKmY)
## limiting rate model using max climate velocity is a better fit

## compare lim rate models to cv and dispersal alone
AIC(lm_mean, lm_mean_limrate, lm_disp) ## mean lim rate model is a slightly better fit than cv alone, much better than disp
AIC(lm_mean) - AIC(lm_mean_limrate)

AIC(lm_max, lm_max_limrate, lm_disp) 
AIC(lm_max) - AIC(lm_max_limrate) ## max lim rate model is a much  better fit than cv alone, much better than disp

aics <- AIC(lm_max, lm_mean, lm_max_limrate, lm_mean_limrate, lm_disp) 

arrange(aics,  AIC)


## models of dispersal limited species only
lm_max_disp = lm(data = filter(data, LimitingRate2 == DispersalPotentialKmY),
                     ShiftKmY ~ DispersalPotentialKmY)
lm_max_clim = lm(data = filter(data, LimitingRate2 == DispersalPotentialKmY),
                 ShiftKmY ~ ClimVeloTKmY_spp_max)
lm_max_clim_mean = lm(data = filter(data, LimitingRate2 == DispersalPotentialKmY),
                 ShiftKmY ~ ClimVeloTKmY_spp)

AIC(lm_max_disp, lm_max_clim, lm_max_clim_mean)
AIC(lm_max_disp) - AIC(lm_max_clim)

lm_mean_disp = lm(data = filter(data, LimitingRate == DispersalPotentialKmY),
                 ShiftKmY ~ DispersalPotentialKmY)
lm_mean_clim = lm(data = filter(data, LimitingRate == DispersalPotentialKmY),
                 ShiftKmY ~ ClimVeloTKmY_spp)

## clim velo model is better
AIC(lm_mean_disp, lm_mean_clim, lm_mean_clim_max)
AIC(lm_mean_disp) - AIC(lm_mean_clim)




data %>%
  filter(LimitingRate2 != ClimVeloTKmY_spp_max) %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ShiftKmY, colour = ClimVeloTKmY_spp_max, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Maximum rate of climate change (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  scale_x_log10() +
  scale_y_log10()



data %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp + sd_cv_sppspecific,
                               DispersalPotentialKmY, ClimVeloTKmY_spp+ sd_cv_sppspecific)) %>%
  filter(LimitingRate == DispersalPotentialKmY) %>%
  summarize(mean_disp = mean(LimitingRate), mean_shift = mean(ShiftKmY),
            mean_diff = mean(ShiftKmY - LimitingRate)) 






## traits 
## body size 
data %>%
  ggplot(aes(x = BodySize, y = DispersalPotentialKmY, shape = group, group = group,
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm")

data %>%
  filter(LimitingRate2 != ClimVeloTKmY_spp_max) %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, shape = group, group = group,
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm") 

data %>%
  filter(LimitingRate2 == ClimVeloTKmY_spp_max) %>%
  ggplot(aes(x = BodySize, y = ShiftKmY, shape = group, 
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm") 

## range size 
data %>%
  filter(Area_km2_range > 1e4) %>% ## get rid of outlier
  ggplot(aes(x = Area_km2_range, y = DispersalPotentialKmY, shape = group, group = group,
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm") 

data %>%
  filter(Area_km2_range > 1e4) %>% ## get rid of outlier
  filter(LimitingRate2 == ClimVeloTKmY_spp_max) %>%
  ggplot(aes(x = Area_km2_range, y = ShiftKmY, shape = group, group = group,
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm") 

data %>%
  filter(Area_km2_range > 1e4) %>% ## get rid of outlier
  filter(LimitingRate2 != ClimVeloTKmY_spp_max) %>%
  ggplot(aes(x = Area_km2_range, y = ShiftKmY, shape = group, group = group,
             colour = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10() +
  geom_smooth(method = "lm") 





data %>%
  ggplot(aes(x = DispersalPotentialKmY, y = ClimVeloTKmY_spp_max, colour = ClimVeloTKmY_spp_max, shape = group)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Maximum potential dispersal rate (km/y)",
       y = "Maximum rate of climate change (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black")+
  scale_x_log10() +
  scale_y_log10()
