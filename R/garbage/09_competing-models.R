## competing models, inspecting model fits and making tables 
library(tidyverse) 
library(gt)
library(AICcmodavg)

## functions to calculate aic, r squared 
r_squ <- function(nls, data) {
  return(summary(lm(ShiftKmY ~ predict(nls), data = data))$r.squared)
}

get_formula <- function(nls) {
  return(summary(nls)$formula)
}

## functions to calculate aic
calculate_aic <- function(residuals, k, n) {
  mse = mean(residuals^2)
  return(2*k + n*log(mse) + n*log(2*pi) + n + 2*k*(k+1)/(n-k-1))
}

## function to calculate bic
calculate_bic <- function(residuals, k, n) {
  ll <- 0.5 * (-1 * n * (log(2 * pi) + 1 - log(n) + log(sum(residuals^2))))
  
  df.ll <- k + 1 
  bic <- -2 * ll + log(n) * df.ll
  return(bic)
}

## r2
calculate_r2 <- function(residuals, response) {
  ss_resid = sum(residuals^2)
  ss_tot = sum((response - mean(response))^2)
  return(1 - ss_resid/ss_tot)
}
## log likelihood
calculate_ll <- function(residuals, k, n) {
  mse = mean(residuals^2)
  return(-n/2*log(2*pi) -n/2*log(mse) - n/2)
}

###################################
##       TESTING EXPANSIONS      ##
###################################
## read in model fits 
exp_models = readRDS("data-processed/modelfits_expansions.rds")
names(exp_models)

## read in data 
exps <- read.csv("data-processed/model-data_expansions.csv")

## version with traits
exp_traits <- filter(exps, !is.na(BodySize), !is.na(Area_km2_range))

## save each model as individual objects 
for (i in 1:length(exp_models)) {
  x <- names(exp_models)[i]
  eval(call("<-", as.name(x), exp_models[[i]]))
}

## model summary
summary(nls_disp)
summary(nls_climvelo)
summary(nls_limrate)

summary(lm_disp)
summary(lm_climvelo)
summary(lm_limrate)
summary(lm_disp_int)

summary(lm_dr)
summary(lm_bs)
summary(lm_rs)

## fit models with fixed slope, fixed intercept
## (calculate the residual error around the 1:1 line)
residuals_disp = exps$ShiftKmY - exps$DispersalPotentialKmY
residuals_climvelo = exps$ShiftKmY - exps$ClimVeloTKmY_spp
residuals_limrate = exps$ShiftKmY - exps$LimitingRate

## calculate AICc, r2, ll for these models
aic_fixed <- lapply(list(residuals_disp, residuals_climvelo, residuals_limrate), 
                    calculate_aic, k = 1, n = nrow(exps))
bic_fixed <- lapply(list(residuals_disp, residuals_climvelo, residuals_limrate), 
                    calculate_bic, k = 1, n = nrow(exps))
r2_fixed <- lapply(list(residuals_disp, residuals_climvelo, residuals_limrate), 
                   calculate_r2, response = exps$ShiftKmY)
ll_fixed <- lapply(list(residuals_disp, residuals_climvelo, residuals_limrate), 
                   calculate_ll,  k = 1, n = nrow(exps))
## note: r2 is negative for first two models, indicating really bad fit 

## make df
fixed = data.frame(Model = c("fixed_disp", "fixed_climvelo", "fixed_limrate"),
           AICc = unlist(aic_fixed),
           BIC = unlist(bic_fixed),
           r_squared = unlist(r2_fixed),
           LL = unlist(ll_fixed),
           n = nrow(exps),
           K = 1,
           Parameter = NA, 
           Estimate = NA,
           Formula = c("DispersalPotentialKmY", "ClimVeloTKmY_spp", "LimitingRate"),
           type = "dispersal")

fixed$`Std. Error` = NA
fixed$`t value` = NA
fixed$`Pr(>|t|)` = NA

## get model equations
formulas <- unlist(lapply(exp_models, get_formula))

## calculate r squared 
exp_rsq <- unlist(list(lapply(exp_models[1:7], FUN = r_squ, data = exps),
                             lapply(exp_models[8:10], FUN = r_squ, data = exp_traits)))

## get n 
n <- c(rep(nrow(exps), 7), rep(nrow(exp_traits), 3))

## get aic
aic_disp <- aictab(cand.set = exp_models[1:7], modnames = names(exp_models)[1:7]) %>%
  data.frame() %>%
  select(Modnames, K, AICc, LL) %>%
  rename("Model" = Modnames) %>%
  mutate(type = "dispersal")

aic_proxy <- aictab(cand.set = exp_models[8:10], modnames = names(exp_models)[8:10]) %>%
  data.frame() %>%
  select(Modnames, K, AICc, LL) %>%
  rename("Model" = Modnames) %>%
  mutate(type = "proxy")

aic_table <- rbind(aic_disp, aic_proxy)

## get BIC 
bic <- data.frame(Model = names(exp_models)[1:10],
           BIC = unlist(lapply(exp_models[1:10], BIC)))

aic_table <- left_join(aic_table, bic) %>%
  select(Model, K, AICc, BIC, LL, type) %>%
  group_by(type) %>%
  arrange(BIC)


## get coefs, join to r squared + se and aic table
coefs <- lapply(exp_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(exp_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = exp_rsq[i]
  cur$n = n[i]
  cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
  coefs[[i]] <- cur
}
coefs <- coefs %>%
  bind_rows() %>%
  left_join(., aic_table) 

## bind to fixed intercept and slope model info
coefs <- rbind(coefs, fixed)

## compete the models, get delta AIC and weights
## dispersal
weights_disp <- qpcR::akaike.weights(unique(coefs$AICc[which(coefs$type == "dispersal")]))
weights_disp <- as.data.frame(weights_disp) %>%
  mutate(AICc = unique(coefs$AICc[which(coefs$type == "dispersal")]))

weights_proxy <- qpcR::akaike.weights(unique(coefs$AICc[which(coefs$type == "proxy")]))
weights_proxy <- as.data.frame(weights_proxy) %>%
  mutate(AICc = unique(coefs$AICc[which(coefs$type == "proxy")]))

coefs <- rbind(weights_disp, weights_proxy) %>%
  left_join(coefs, .)

## sort by model, then type, then rank
coefs <- group_by(coefs, type) %>%
  arrange(deltaAIC, .by_group = TRUE)

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, everything()) 

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, "")))

## make gt table
table_exp_disp <- coefs %>% 
  ungroup() %>%
  filter(type == "dispersal") %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range expansions: dispersal models"
  ) 

table_exp_proxy <- coefs %>% 
  ungroup() %>%
  filter(type == "proxy") %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range expansions: proxy trait models"
  ) 

gtsave(table_exp_disp, path = "figures/model_results", filename = "table_expansions_dispersal.png")
gtsave(table_exp_proxy, path = "figures/model_results", filename = "table_expansions_proxytraits.png")



#####################################
##       TESTING CONTRACTIONS      ##
#####################################
## read in model fits 
cont_models = readRDS("data-processed/modelfits_contractions.rds")
names(cont_models)

## read in data 
contractions <- read.csv("data-processed/model-data_contractions.csv")

## version with traits
contractions_traits <- filter(contractions, !is.na(BodySize), !is.na(Area_km2_range))

## save each model as individual objects 
for (i in 1:length(cont_models)) {
  x <- names(cont_models)[i]
  eval(call("<-", as.name(x), cont_models[[i]]))
}

## model summary
summary(nls_disp_cont)
summary(nls_climvelo_cont)
summary(nls_limrate_cont)

summary(lm_disp_cont)
summary(lm_climvelo_cont)
summary(lm_limrate_cont)
summary(lm_disp_int_cont)

summary(lm_dr_cont)
summary(lm_bs_cont)
summary(lm_rs_cont)

## fit models with fixed slope, fixed intercept
## (calculate the residual error around the 1:1 line)
residuals_disp_cont = contractions$ShiftKmY - contractions$DispersalPotentialKmY
residuals_climvelo_cont = contractions$ShiftKmY - contractions$ClimVeloTKmY_spp
residuals_limrate_cont = contractions$ShiftKmY - contractions$LimitingRate

## calculate AICc, r2, ll for these models
aic_fixed_cont <- lapply(list(residuals_disp_cont, residuals_climvelo_cont, residuals_limrate_cont), 
                    calculate_aic, k = 1, n = nrow(exps))
r2_fixed_cont <- lapply(list(residuals_disp_cont, residuals_climvelo_cont, residuals_limrate_cont), 
                   calculate_r2, response = exps$ShiftKmY)
ll_fixed_cont <- lapply(list(residuals_disp_cont, residuals_climvelo_cont, residuals_limrate_cont), 
                   calculate_ll,  k = 1, n = nrow(exps))
## note: r2 is negative for first two models, indicating really bad fit 

## make df
fixed = data.frame(Model = c("fixed_disp", "fixed_climvelo", "fixed_limrate"),
                   AICc = unlist(aic_fixed_cont), 
                   r_squared = unlist(r2_fixed_cont),
                   LL = unlist(ll_fixed_cont),
                   n = nrow(contractions),
                   K = 1,
                   Parameter = NA, 
                   Estimate = NA,
                   Formula = c("DispersalPotentialKmY", "ClimVeloTKmY_spp", "LimitingRate"),
                   type = "dispersal")

fixed$`Std. Error` = NA
fixed$`t value` = NA
fixed$`Pr(>|t|)` = NA

## get model equations
formulas <- unlist(lapply(cont_models, get_formula))

## calculate r squared 
cont_rsq <- unlist(list(lapply(cont_models[1:7], FUN = r_squ, data = contractions),
                       lapply(cont_models[8:10], FUN = r_squ, data = contractions_traits)))

## get n 
n <- c(rep(nrow(contractions), 7), rep(nrow(contractions_traits), 3))

## get aic
aic_disp <- aictab(cand.set = cont_models[1:7], modnames = names(cont_models)[1:7]) %>%
  data.frame() %>%
  select(Modnames, K, AICc, LL) %>%
  rename("Model" = Modnames) %>%
  mutate(type = "dispersal")

aic_proxy <- aictab(cand.set = cont_models[8:10], modnames = names(cont_models)[8:10]) %>%
  data.frame() %>%
  select(Modnames, K, AICc, LL) %>%
  rename("Model" = Modnames) %>%
  mutate(type = "proxy")

aic_table <- rbind(aic_disp, aic_proxy)

## get coefs, join to r squared + se and aic table
coefs <- lapply(cont_models, FUN = function(x) {summary(x)$coefficients})
for (i in 1:length(coefs)) {
  cur <- as.data.frame(coefs[[i]])
  cur$Model = names(cont_models)[i]
  cur$Parameter = rownames(cur)
  cur$r_squared = exp_rsq[i]
  cur$n = n[i]
  cur$Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[i,2]
  coefs[[i]] <- cur
}
coefs <- coefs %>%
  bind_rows() %>%
  left_join(., aic_table) 

## bind to fixed intercept and slope model info
coefs <- rbind(coefs, fixed)

## compete the models, get delta AIC and weights
## dispersal
weights_disp <- qpcR::akaike.weights(unique(coefs$AICc[which(coefs$type == "dispersal")]))
weights_disp <- as.data.frame(weights_disp) %>%
  mutate(AICc = unique(coefs$AICc[which(coefs$type == "dispersal")]))

weights_proxy <- qpcR::akaike.weights(unique(coefs$AICc[which(coefs$type == "proxy")]))
weights_proxy <- as.data.frame(weights_proxy) %>%
  mutate(AICc = unique(coefs$AICc[which(coefs$type == "proxy")]))

coefs <- rbind(weights_disp, weights_proxy) %>%
  left_join(coefs, .)

## sort by model, then type, then rank
coefs <- group_by(coefs, type) %>%
  arrange(deltaAIC, .by_group = TRUE)

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, `Std. Error`, everything()) 

## replace NA with blank
coefs <- coefs %>% 
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~replace_na(.x, "")))

## make gt table
table_cont_disp <- coefs %>% 
  ungroup() %>%
  filter(type == "dispersal") %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range contractions: dispersal models"
  ) 

table_cont_proxy <- coefs %>% 
  ungroup() %>%
  filter(type == "proxy") %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range contractions: proxy trait models"
  ) 

gtsave(table_cont_disp, path = "figures/model_results", filename = "table_contractions_dispersal.png")
gtsave(table_cont_proxy, path = "figures/model_results", filename = "table_contractions_proxytraits.png")

