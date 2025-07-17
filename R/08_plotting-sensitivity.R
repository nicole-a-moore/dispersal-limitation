## making plots to visualize the sensitivty of measurements to which climate velocity scale / which dispersal distance measurement is used 
library(tidyverse)
theme_set(theme_bw())

#############################
##       PREPARE DATA      ##
#############################
## read in data 
dd <- read.csv("data-processed/v3_with-cv.csv")

## filter to leading edge shifts with positive mean climate velocity 
dd <- filter(dd, !Param %in% c("O", "TE") & ClimVeloKmY_RelScale >= 0)

## get rid of range contractions that are father than 1 sd from the mean shift 
dd <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

##########################
##       PLOT DATA      ##
##########################

## violin plot with scale of climate velocity on x axis, measure of climate velocity on y axis 
## colour by what is limiting

## make column for each climate velocity scale that says whether or not max dispersal is slower than climate velocity 
data <- dd %>%
  select(-sdClimVeloKmY_110km, -sdClimVeloKmY_50km, -sdClimVeloKmY_25km, -sdClimVeloKmY_RelScale) %>%
  gather(key = "cv_type", value = "ClimVeloKmY", c(ClimVeloKmY_25km, ClimVeloKmY_50km, 
                                                   ClimVeloKmY_110km, ClimVeloKmY_RelScale)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY <= ClimVeloKmY, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloKmY, NA)) %>%
  mutate(cv_lab = ifelse(cv_type == "ClimVeloKmY_25km", "25km", 
                         ifelse(cv_type == "ClimVeloKmY_50km", "50km", 
                                ifelse(cv_type == "ClimVeloKmY_110km", "110km", 
                                       "Species-relevant scale")))) %>%
  mutate(cv_lab = factor(cv_lab, levels = c("25km", "50km", "110km", "Species-relevant scale"), ordered = TRUE))

## arrange by what is limiting
data <- arrange(data, what_is_limiting)

data %>%
  ggplot(aes(x = cv_lab, y = ClimVeloKmY, fill = ClimVeloKmY, shape = group)) +
  geom_violin(inherit.aes = F, data = data , aes(x = cv_lab, y = ClimVeloKmY)) +
  geom_point(position = position_jitter(),
             colour = ifelse(data$what_is_limiting == "Climate", "transparent", "black")) +
  scale_shape_manual(values = c(21,24,22,23)) +
  scale_fill_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  labs(x = "Spatial resolution of climate data", y = "Mean velocity of climate change (km/y)", 
       shape = "") + 
  guides(size = "legend", fill = "none")

## save
ggsave(path = "figures/visualizing-sensitivity", filename = "sensitivity_cv-scale.png",
       width = 6, height = 3.5)

data %>%
  ggplot(aes(x = cv_lab, y = ShiftKmY, fill = ClimVeloKmY, shape = group)) +
  geom_violin(inherit.aes = F, data = data , aes(x = cv_lab, y = ClimVeloKmY)) +
  geom_point(position = position_jitter(),
             colour = ifelse(data$what_is_limiting == "Climate", "transparent", "black")) +
  scale_shape_manual(values = c(21,24,22,23)) +
  scale_fill_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  labs(x = "Spatial scale of climate velocity", y = "Range shift rate (km/y)", 
       shape = "") + 
  guides(size = "legend", fill = "none")


## plot that shows mean + 3rd quantile and max dispersal
data <- dd %>%
  mutate(facet = ifelse(DispersalPotentialKmY <= q3ClimVeloKmY_RelScale, 
                        "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_median = ifelse(MedianDispersalPotentialKmY <= ClimVeloKmY_RelScale, "Dispersal", "Climate"),
         what_is_limiting_median_upper = ifelse(MedianDispersalPotentialKmY <= q3ClimVeloKmY_RelScale, 
                                                "Dispersal", "Climate"),
         what_is_limiting = ifelse(DispersalPotentialKmY <= ClimVeloKmY_RelScale, "Dispersal", "Climate"),
         what_is_limiting_upper = ifelse(DispersalPotentialKmY <= q3ClimVeloKmY_RelScale, "Dispersal", "Climate")) %>%
  mutate(sensitive = ifelse(what_is_limiting_median == what_is_limiting_median_upper & 
                              what_is_limiting_median_upper ==  what_is_limiting & 
                              what_is_limiting == what_is_limiting_upper,
                            "no", "yes")) %>%
  mutate(cat = ifelse(sensitive == "no" & what_is_limiting == "Dispersal", 
                      "Dispersal always > velocity\nof climate change", 
                ifelse(sensitive == "no" & what_is_limiting == "Climate", 
                       "Dispersal always < velocity\nof climate change", 
                       "Sensitive to type of\nmeasurement used"))) %>%
  gather(key = "disp_measure_type", value = "disp_measure", 
         c(MedianDispersalPotentialKmY, DispersalPotentialKmY)) 

data <- data %>%
  group_by(cat) %>%
  arrange(-ClimVeloKmY_RelScale) %>%
  mutate(rank = 1:n()) 

p <- data %>%
  mutate(disp_measure_type = ifelse(disp_measure_type == "MedianDispersalPotentialKmY", "Median potential dispersal rate",
                                    "Maximum potential dispersal rate")) %>%
  mutate(disp_measure = ifelse(disp_measure >= 15, 15, disp_measure)) %>%
  ggplot(aes(x = rank, y = disp_measure, colour = disp_measure_type)) +
  geom_pointrange(inherit.aes = F, 
                  data = data, aes(x = rank, y = ClimVeloKmY_RelScale, 
                                   ymax = q3ClimVeloKmY_RelScale,
                                   ymin = q1ClimVeloKmY_RelScale),
                  size = 0.1, linewidth = 0.08) +
  labs(x = "Range shift observation", 
       colour = "", 
       y = "Rate (km/y)") +
  facet_wrap(~cat) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        plot.margin = margin(c(0,0,0,0)),
        legend.margin = margin(c(1, 1, 1, 1))) +
  geom_abline(slope = 0, intercept = 0) + 
  geom_point(size = 0.85) 

## hack a legend for point range 
legend <- data %>%
  filter(cat == "Sensitive to type of\nmeasurement used") %>%
  mutate(cat = "q1, mean, and q3 velocity of climate change") %>%
  ggplot(aes(x = rank, y = ClimVeloKmY_RelScale, 
             ymax = q3ClimVeloKmY_RelScale,
             ymin = q1ClimVeloKmY_RelScale, colour = cat)) +
  geom_point() +
  geom_linerange(size = 0.1, linewidth = 0.1) +
  theme(legend.position = "bottom", 
        legend.margin = margin(c(1, 1, 5, 1))) +
  scale_color_manual(values = c("black")) +
  labs(colour = " ") +
  scale_y_continuous(expand = c(0,0))

legend <- ggpubr::get_legend(legend)

plot_grid(p, legend, ncol = 1, nrow = 2, 
          rel_heights = c(15, 1))

## save
ggsave(path = "figures/visualizing-sensitivity", filename = "sensitivity_dispersal-and-cv.png",
       width = 8, height = 3.5)

