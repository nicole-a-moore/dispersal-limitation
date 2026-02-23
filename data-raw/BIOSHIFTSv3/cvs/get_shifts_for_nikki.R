
sa_cvs <- readr::read_csv("data-processed/final/cv.csv")
sa_cvs

shifts <- get_shifts() %>%add_methods()
shifts %>%
  count(method_id, midpoint_firstperiod, midpoint_secondperiod)
sa_cvs %>%
  filter(article_id == "A007",
         poly_id == "P1") %>%
  left_join(shifts)

sp_csvs <- readr::read_csv("data-processed/final/sp_cv.csv")



sa_cvs %>%
  slice(8) %>%
  left_join(id_key %>% distinct(old_id, article_id, poly_id))

id_key %>%
  slice(8) %>%
  left_join(sa_cvs) %>% glimpse()

sa_cvs %>%
  inner_join(id_key %>% slice(8)) %>%
  select(-method_id) %>%
  distinct() %>% glimpse()


method_key <- get_shifts() %>%
  distinct(id, article_id, poly_id, type, method_id, param, eco) %>%
  add_methods() %>%
  distinct(article_id, poly_id, eco, type, param, method_id, midpoint_firstperiod,
           midpoint_secondperiod, duration)





# make bioshifts key for Nikki --------------------------------------------

orig <- readr::read_csv("data-raw/bioshifts_v3.csv") %>%
  janitor::clean_names()

sp_cvs <- readr::read_csv("data-processed/final/sp_cv.csv")


# filter to lat -----------------------------------------------------------

# filter original V3 to latitudinal shifts
orig <- orig %>% filter(type == "LAT")
orig %>% glimpse()

# filter species' CVs to latitudinal shifts
sp_cvs <- sp_cvs %>% filter(type == "LAT")

# attach species' level CVs to new bioshifts
new_lat <- get_shifts(type = "LAT",
                      group = c("Birds", "Vascular Plants")) %>%
  # add methods so we can match by duration
  add_methods() %>% 
  # add articles so we can match by article name
  add_articles() %>%
  # join species-level CVs
  left_join(sp_cvs)

new_lat %>% glimpse()



# isolate IDs -------------------------------------------------------------
# filter original to columns that we'll use to match -- 
# shift rate and certain methods
orig_filtered <- orig %>% 
  select(id, article, subsp, subsp_or_pop, param, sp_name_publication, sp_name_checked, rate,
         midpoint_firstperiod, midpoint_secondperiod, duration, obs_type,article ) %>%
  # rename ID column so it's different than new "id" column,
  # change rate to "calc_rate" to match new. 
  rename(old_id = id,
         calc_rate = rate)

# join by methods and stuff
joined <- new_lat %>% left_join(orig_filtered) %>%
  relocate(old_id)


# find birds and plants from original bioshifts
old_birds_plants <- orig %>% filter(type == "LAT") %>% 
  filter(class == "Aves" | phylum == "Tracheophyta") %>% distinct(sp_name_checked) %>% pull(sp_name_checked)

# see how many are in new set
(old_birds_plants %in% joined$sp_name_checked) %>% table() # all of them. 

# see if any new ones are not in old set
(joined$sp_name_checked %in% old_birds_plants) %>% table() # all of them. 

# rearrange some columns
joined <- joined %>%
  relocate(c(subsp_or_pop, sublocality_or_pop), .after = subsp) %>% glimpse() %>%
  select(-doi, -id_bioshifts_v1, -id_core) %>% glimpse()

# see if there are missing CVs
joined %>% 
  rowwise() %>%
  mutate(cv_na = sum(is.na(across(cv_temp_mean_res1km:cv_temp_q75_res110km)))) %>%
  ungroup() %>%
  filter(cv_na > 0) %>% glimpse()
# ok yes there's ~500 birds and plant shifts that don't have matching ranges. 
# get rid of those

joined <- joined %>% filter(!is.na(range_source))

# check missing vals again
joined %>%
  rowwise() %>%
  mutate(cv_na = sum(is.na(across(cv_temp_mean_res1km:cv_temp_q75_res110km)))) %>%
  ungroup() %>%
  filter(cv_na > 0) %>% glimpse() %>%
  count(article_id)
# ok there's 13 rows missing some combo of values... I'll figure that out. 


# save
readr::write_csv(joined,
                 "sp_shifts_and_cvs_for_nikki.csv")
