## competing models, inspecting model fits and making tables 
library(tidyverse)
library(gt)
library(AICcmodavg)

## function to calculate r squared 
r_squ <- function(nls, data) {
  return(summary(lm(ShiftKmY ~ predict(nls), data = data))$r.squared)
}

get_formula <- function(nls) {
  return(summary(nls)$formula)
}

###################################
##       TESTING EXPANSIONS      ##
###################################
## read in model fits 
expansion_models = readRDS("data-processed/modelfits_expansions.rds")
names(expansion_models)

## read in data 
expansions <- read.csv("data-processed/model-data_expansions.csv")

## version with traits
expansions_traits <- filter(expansions, !is.na(BodySize), !is.na(Area_km2_range))

## save each model as individual objects 
for (i in 1:length(expansion_models)) {
  x <- names(expansion_models)[i]
  eval(call("<-", as.name(x), expansion_models[[i]]))
}

## model summaries
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

## plot residuals 
plot(nls_disp) 
plot(nls_climvelo)
plot(nls_limrate)

plot(lm_disp) 
plot(lm_climvelo)
plot(lm_limrate)
plot(lm_int)

plot(lm_dr)
plot(lm_bs)
plot(lm_rs)

hist(residuals(nls_disp))
hist(residuals(nls_climvelo))
hist(residuals(nls_limrate))

hist(residuals(lm_disp))
hist(residuals(lm_climvelo))
hist(residuals(lm_limrate))
hist(residuals(lm_disp_int))

hist(residuals(lm_dr))
hist(residuals(lm_bs))
hist(residuals(lm_rs))

## get model equations
formulas <- unlist(lapply(expansion_models, get_formula))

## calculate r squared 
expansion_rsq <- unlist(list(lapply(expansion_models[1:7], FUN = r_squ, data = expansions),
                             lapply(expansion_models[8:10], FUN = r_squ, data = expansions_traits)))

## get n 
n <- c(rep(nrow(expansions), 7), rep(nrow(expansions_traits), 3))

## get AIC, delta AIC and rank 
aic_disp <- aictab(cand.set = expansion_models[1:7], modnames = names(expansion_models)[1:7]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(), 
         type = "dispersal") %>% 
  select(-Cum.Wt)

aic_proxy <- aictab(cand.set = expansion_models[8:10], modnames = names(expansion_models)[8:10]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(),
         type = "proxy trait") %>% 
  select(-Cum.Wt)

aic_table <- rbind(aic_disp, aic_proxy)

## get coefs, join to r squared and aic table
coefs <- lapply(expansion_models, coef) %>%
  bind_rows() %>%
  mutate(Model = names(expansion_models), 
         r_squared = expansion_rsq, 
         n = n, 
         Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[,2]) %>%
  gather(key = "Parameter", value = "Estimate", c(int, m, m_cv, m_int)) %>% 
  left_join(., aic_table)

## sort by model, then type, then rank
coefs <- group_by(coefs, type) %>%
  arrange(rank, desc(Parameter), .by_group = TRUE)

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## get rid of slopes that are na
coefs <- coefs[!is.na(coefs$Estimate),]

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, everything()) 

## make gt table
table_expansions <- coefs %>% 
  ungroup() %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range expansions",
  ) %>%
  tab_row_group(
    label = "Proxy trait models",
    rows = 14:25
  ) %>%
  tab_row_group(
    label = "Dispersal limitation models",
    rows = 1:13
  ) 

gtsave(table_expansions, path = "figures/model_results", filename = "table_expansions.png")


#####################################
##       ACCOUNTING FOR ERROR      ##
#####################################
## read in model fits 
error_models = readRDS("data-processed/modelfits_error.rds")
names(error_models)

## read in data 
error <- read.csv("data-processed/model-data_error.csv")

## version with traits
error_traits <- filter(error, !is.na(BodySize), !is.na(Area_km2_range))

## save each model as individual objects 
for (i in 1:length(error_models)) {
  x <- names(error_models)[i]
  eval(call("<-", as.name(x), error_models[[i]]))
}

## model summary
summary(nls_disp_error)
summary(nls_climvelo_error)
summary(nls_limrate_error)

summary(lm_disp_error)
summary(lm_climvelo_error)
summary(lm_limrate_error)
summary(lm_disp_int_error)

summary(lm_dr_error)
summary(lm_bs_error)
summary(lm_rs_error)

## plot residuals 
plot(nls_disp_error) 
plot(nls_climvelo_error)
plot(nls_limrate_error)

plot(lm_disp_error) 
plot(lm_climvelo_error)
plot(lm_limrate_error)
plot(lm_disp_int_error)

plot(lm_dr_error)
plot(lm_bs_error)
plot(lm_rs_error)

hist(residuals(nls_disp_error))
hist(residuals(nls_climvelo_error))
hist(residuals(nls_limrate_error))

hist(residuals(lm_disp_error))
hist(residuals(lm_climvelo_error))
hist(residuals(lm_limrate_error))
hist(residuals(lm_disp_int_error))

hist(residuals(lm_dr_error))
hist(residuals(lm_bs_error))
hist(residuals(lm_rs_error))

## get model equations
formulas <- unlist(lapply(error_models, get_formula))

## calculate r squared 
error_rsq <- unlist(list(lapply(error_models[1:7], FUN = r_squ, data = error),
                             lapply(error_models[8:10], FUN = r_squ, data = error_traits)))

## get n 
n <- c(rep(nrow(error), 7), rep(nrow(error_traits), 3))

## get AIC, delta AIC and rank 
aic_disp <- aictab(cand.set = error_models[1:7], modnames = names(error_models)[1:7]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(), 
         type = "dispersal") %>% 
  select(-Cum.Wt)

aic_proxy <- aictab(cand.set = error_models[8:10], modnames = names(error_models)[8:10]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(),
         type = "proxy trait") %>% 
  select(-Cum.Wt)

aic_table <- rbind(aic_disp, aic_proxy)

## get coefs, join to r squared and aic table
coefs <- lapply(error_models, coef) %>%
  bind_rows() %>%
  mutate(Model = names(error_models), 
         r_squared = error_rsq, 
         n = n, 
         Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[,2]) %>%
  gather(key = "Parameter", value = "Estimate", c(int, m, m_int, m_cv)) %>%
  left_join(., aic_table) 

## sort by model, then type, then rank
coefs <- group_by(coefs, type) %>%
  arrange(rank, desc(Parameter), .by_group = TRUE)

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## get rid of slopes that are na
coefs <- coefs[!is.na(coefs$Estimate),]

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, everything()) 

## make gt table
table_error <- coefs %>% 
  ungroup() %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range expansions + error",
  ) %>%
  tab_row_group(
    label = "Proxy trait models",
    rows = 14:25
  ) %>%
  tab_row_group(
    label = "Dispersal limitation models",
    rows = 1:13
  ) 

gtsave(table_error, path = "figures/model_results", filename = "table_expansions-with-error.png")


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

## plot residuals 
plot(nls_disp_cont) 
plot(nls_climvelo_cont)
plot(nls_limrate_cont)

plot(lm_disp_cont) 
plot(lm_climvelo_cont)
plot(lm_limrate_cont)
plot(lm_disp_int_cont)

plot(lm_dr_cont) 
plot(lm_bs_cont)
plot(lm_rs_cont)

hist(residuals(nls_disp_cont))
hist(residuals(nls_climvelo_cont))
hist(residuals(nls_limrate_cont))

hist(residuals(lm_disp_cont))
hist(residuals(lm_climvelo_cont))
hist(residuals(lm_limrate_cont))
hist(residuals(lm_disp_int_cont))

hist(residuals(lm_dr_cont))
hist(residuals(lm_bs_cont))
hist(residuals(lm_rs_cont))

## get model equations
formulas <- unlist(lapply(cont_models, get_formula))

## calculate r squared 
cont_rsq <- unlist(list(lapply(cont_models[1:7], FUN = r_squ, data = contractions),
                        lapply(cont_models[8:10], FUN = r_squ, data = contractions_traits)))

## get n 
n <- c(rep(nrow(contractions), 7), rep(nrow(contractions_traits), 3))

## get AIC, delta AIC and rank 
aic_disp <- aictab(cand.set = cont_models[1:7], modnames = names(cont_models)[1:7]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(), 
         type = "dispersal") %>% 
  select(-Cum.Wt)

aic_proxy <- aictab(cand.set = cont_models[8:10], modnames = names(cont_models)[8:10]) %>%
  data.frame() %>%
  rename("Model" = Modnames) %>%
  mutate(rank = row_number(),
         type = "proxy trait") %>% 
  select(-Cum.Wt)

aic_table <- rbind(aic_disp, aic_proxy)

## get coefs, join to r squared and aic table
coefs <- lapply(cont_models, coef) %>%
  bind_rows() %>%
  mutate(Model = names(cont_models), 
         r_squared = cont_rsq, 
         n = n, 
         Formula = str_split_fixed(sapply(formulas, toString), "ShiftKmY, ", 2)[,2]) %>%
  gather(key = "Parameter", value = "Estimate", c(int, m, m_cv, m_int)) %>%
  left_join(., aic_table) 

## sort by model, then type, then rank
coefs <- group_by(coefs, type) %>%
  arrange(rank, desc(Parameter), .by_group = TRUE)

## round everything to 2 decimal places 
nums <- unlist(lapply(coefs, is.numeric), use.names = FALSE)  
coefs[,nums] <- round(coefs[,nums], 2)

## get rid of slopes that are na
coefs <- coefs[!is.na(coefs$Estimate),]

## reorder/rename columns 
coefs <- coefs %>%
  rename("R2" = r_squared) %>%
  select(Model, Formula, Parameter, Estimate, everything()) 

## make gt table
table_cont <- coefs %>% 
  ungroup() %>%
  select(-type, -Model) %>%
  gt() %>%
  tab_header(
    title = "Range contractions",
  ) %>%
  tab_row_group(
    label = "Proxy trait models",
    rows = 14:25
  ) %>%
  tab_row_group(
    label = "Dispersal limitation models",
    rows = 1:13
  ) 

gtsave(table_cont, path = "figures/model_results", filename = "table_contractions.png")

