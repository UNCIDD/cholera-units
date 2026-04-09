#-----------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 01_data_prep.R
# Purpose: Initial data import and cleaning
#-----------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize
#-----------------------------------------------------------#

# load packages
library(here)
library(spData)    
library(sf)
library(giscoR)
library(glue)
library(janitor)
library(tidyverse)
library(rlang)

# source project functions
source(here("analysis","00_functions_settings.R"), local = T)

#-----------------------------------------------------------#
# Import & preliminary cleaning----
#-----------------------------------------------------------#

## Sequence data----
sequences <- read_csv(here("data/raw/sequences.csv"))

###WHO cases----
who_cases <- read_csv(here("data/raw/who_annual_case_reports.csv"))

### Phylogeography----

phylo_raw <- read_tsv(here("data/raw/phylogeography/2024-10-04_constant_relaxed_run",
                           "2024-10-03_constant_relaxed_discrete.pruned.log")) 

### Introductions----
#T1-T8 estimates from Weill et al 2017,
#T9-T15 estimates from our phylogeographic analysis
#T16 & T17 lineage estimates from Xiao et al
introductions_raw <- read_csv(here("data/raw","introductions_w.csv")) 

#-----------------------------------------------------------#
# Geography & population sizes----
#-----------------------------------------------------------#

### Subregions----
#Define subregions
subregions <- read_delim(here("data/raw/UN_regions.csv"),
                         delim = ";",
                         guess_max = 1500) |> 
  janitor::clean_names() |> 
  filter(region_name=="Africa") |> 
  select(sub_region_name, intermediate_region_name, country = country_or_area) |> 
  mutate(subregion = str_replace_all(if_else(is.na(intermediate_region_name), 
                                             sub_region_name, intermediate_region_name),
                                     "Middle", "Central")) |> 
  select(subregion, country) |> 
  clean_df() |> 
  mutate(subregion = str_to_lower(subregion))

###Shapefile----
full_africa_sf <- gisco_get_countries(region = "Africa") |> 
  janitor::clean_names() |> 
  select(country = name_engl, geometry) |> 
  clean_df() |> 
  # remove islands not part of analysis
  filter_out(country %in% islands_exclude)
#remove prince edward island & marion island
full_africa_sf$geometry[[which(full_africa_sf$country=="south africa")]][1:2] <- NULL

### Population sizes----
#population sizes, and pull out just 2013 populations
population_sizes <- population |> 
  clean_df() |> 
  filter(country %in% full_africa_sf$country) |> 
  arrange(country, year)
population_sizes_2013 <- population_sizes |> 
  filter(year==2013) |> 
  select(-year)

# Add population size and subregion to shapefile
full_africa_sf <- full_africa_sf |> 
  left_join(population_sizes_2013, by = "country") |> 
  rename(pop = population) |> 
  left_join(subregions, by = "country")

#-----------------------------------------------------------#
# Manipulate introductions df----
#-----------------------------------------------------------#

subregions_w <- subregions |> 
  group_by(subregion) |> 
  mutate(n = row_number()) |> 
  pivot_wider(id_cols = "subregion", values_from = country, 
              names_from = n, names_prefix = "cntry")

# Reshape, exclude countries not in analysis
introductions <- introductions_raw |> 
  pivot_longer(cols = c(starts_with("loc"), starts_with("p")), 
               names_to = c(".value","location_order"), 
               names_pattern = "(.*)([1-5])",
               values_drop_na = TRUE) |> 
  rename(country = loc) |> 
  clean_df() |> 
  left_join(subregions_w, by = c("country" = "subregion")) |> 
  pivot_longer(cols = starts_with("cntry"), names_to = "n", 
               values_to= "country_sr", names_prefix = "cntry") |> 
  select(-n) |> 
  distinct() |> 
  mutate(region = if_else(country %in% subregions$subregion, 1, 0),
         country = if_else(country %in% subregions$subregion, country_sr, country)) |> 
  group_by(te, country) |> 
  mutate(n = n()) |> 
  group_by(te) |> 
  mutate(dup = if_else(n>1 & !is.na(country_sr), 1, 0)) |> 
  ungroup() |> 
  filter(dup==0) |> 
  select(-country_sr, -n, -dup) |> 
  group_by(te, location_order) |> 
  mutate(prob_dist = p/n()) |> 
  ungroup()

#-----------------------------------------------------------#
# Phylogeographic analysis output----
#-----------------------------------------------------------#

# Reshape to long dataset with transition rates
phylo <- phylo_raw |> 
  pivot_longer(cols = starts_with("Location"), 
               names_to = c(".value","countries"),
               names_pattern = c("(Location\\..*\\.)(.*\\..*)")) |> 
  mutate(countries = str_to_title(str_replace(countries, "\\.", "\\, ")),
         countries2 = countries) |> 
  janitor::clean_names() |> 
  separate_wider_delim(countries2, delim = ", ", names = c("country1", "country2")) |> 
  select(country1, country2, countries, starts_with("location"), state)

#-----------------------------------------------------------#
# Combine sequence & case data----
#-----------------------------------------------------------#

cases_lineages <- sequences |> 
  select(country, year, te) |>
  mutate(te = if_else(str_detect(te, "sporadic"), "S", te)) |> 
  full_join(who_cases, by = c("country", "year")) |> 
  arrange(country, year, te)|> 
  mutate(te = factor(te, levels = te_order, ordered = T),
         cases = replace_na(cases, 0)) |> 
  group_by(country, year, te) |> 
  mutate(sample_num = if_else(!is.na(te), row_number(), 0),
         samples = max(sample_num)) |> 
  ungroup() |> 
  filter(sample_num==samples) |> 
  group_by(country, year) |> 
  mutate(tot_samples = sum(samples)) |> 
  ungroup() |> 
  distinct() |> 
  select(-sample_num)


#check
cases_lineages |> 
  group_by(country, te) |> 
  summarize(any_seq = sum(tot_samples)>0) |> 
  ungroup() |> 
  left_join(cases_lineages |> 
              group_by(country) |> 
              summarize(any_cases = sum(cases)>0) |> 
              ungroup(), by = "country") |> 
  distinct(any_seq, any_cases)

#-----------------------------------------------------------#
# Save out----
#-----------------------------------------------------------#

saveRDS(cases_lineages, here("data", "cases_lineages.RDS"))
saveRDS(full_africa_sf, here("data", "full_africa_sf.RDS"))
saveRDS(subregions, here("data", "subregions.RDS"))
saveRDS(introductions, here("data", "introductions.RDS"))
saveRDS(phylo, here("data", "phylogeography.RDS"))




