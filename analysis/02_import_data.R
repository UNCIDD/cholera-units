##import & manipulate case & lineage data

## Initialize ----------------------------------------------------------------
## ---- initialize1
library(tidyverse)
library(spData) 
library(sf)
library(here)

source(here("analysis","00_functions_settings.R"), local = T)


#-----------------------------------------------------------#
# Import----
#-----------------------------------------------------------#

### Shapefiles----
#Shapefile - full continent
full_africa_sf <- readRDS(here("data", "full_africa_sf.RDS")) 

### Subregions----
subregions <- readRDS(here("data", "subregions.RDS")) 

### Lineages & Cases----
#bring in lineage & case data
cases_lineages_full <- readRDS(here("data", "cases_lineages.RDS")) 

#-----------------------------------------------------------#
# Geography----
#-----------------------------------------------------------#

# Shapefile excluding North Africa
africa_sf <- full_africa_sf |> 
  filter_out(country %in% north_africa_exclude)

#### Combine Sudan & South Sudan----
if(sep_sudan==F){
  africa_sf <- africa_sf |> 
    mutate(country = if_else(str_detect(country, "sudan"), "sudan", country),
           subregion = if_else(country=="sudan", "eastern africa", subregion)) |> 
    group_by(country) |> 
    mutate(geometry = st_union(geometry),
           pop = sum(pop)) |> 
    ungroup() |> 
    distinct() 
  subregions <- subregions |> 
    mutate(subregion = if_else(country=="sudan", "eastern africa", subregion))
}

#-----------------------------------------------------------#
# Cases & lineages----
#-----------------------------------------------------------#

if(sep_sudan==F){
  cases_lineages_full <- cases_lineages_full |> 
    mutate(country = if_else(str_detect(country, "sudan"), "sudan", country)) |> 
    group_by(country, year, te) |> 
    mutate(samples = sum(samples),
           tot_samples = sum(tot_samples)) |> 
    distinct() |> 
    group_by(country, year) |> 
    mutate(cases = sum(cases),
           tot_samples = sum(samples)) |> 
    ungroup() |> 
    distinct() |> 
    filter(!is.na(te) | (is.na(te) & tot_samples==0))
}

cases_lineages <- cases_lineages_full |> 
  # subset by year if min_year is updated 
  filter(year>=min_year) |> 
  left_join(subregions, by = "country") |> 
  select(country, subregion, year, cases, te, samples, tot_samples) |> 
  distinct() |> 
  arrange(country, year)

if(test_subset){
  cases_lineages <- cases_lineages |> 
    filter(country %in% c("drc", "kenya", "tanzania", "uganda", "rwanda", "burundi"))
}

#### Country & te lists----
all_countries <- unique(cases_lineages_full$country)
countries_subset <- unique(cases_lineages$country)
te_subset <- te_order[which(te_order %in% as.character(unique((cases_lineages |> filter(!is.na(te)))$te)))] 

#-----------------------------------------------------------#
# Arrivals/Introductions----
#-----------------------------------------------------------#

# Observed dates that lineages are first seen in each country
first_seen <- cases_lineages_full |> 
  select(country, year, te) |> 
  filter_out(is.na(te) | te=="S") |> 
  group_by(te) |> 
  mutate(first_seen = min(year)) |> 
  ungroup() |> 
  filter(year==first_seen) |> 
  arrange(te, country) |> 
  left_join(subregions, by = "country") |> 
  select(te, country, subregion, first_seen) |> 
  distinct() 

# Arrival times of lineages by country
#first, create shell of countries & lineages
arrival_times <- tibble(country = rep(all_countries, 
                                          times = length(te_order[te_order!="S"])),
                            te = rep(te_order[te_order!="S"], 
                                     each = length(all_countries))) |> 
  left_join(cases_lineages_full |> 
              filter_out(is.na(te)) |> 
              select(country, year, te) |> 
              distinct(), by = c("country", "te")) |> 
  mutate(year = replace_na(year, 2099)) |> 
  left_join(first_seen |> select(te, first_seen) |> distinct(), by = "te") |> 
  group_by(country, te) |> 
  mutate(arrival_year = min(year)) |> 
  ungroup() |> 
  select(country, te, arrival_year, first_seen) |> 
  distinct() |> 
  mutate(obs_time_to_arr = if_else(arrival_year!=2099,arrival_year-first_seen,as.numeric(NA))) |> 
  arrange(te, country) 
  

# List of likely seed countries
seed_countries <- unique(c(unique(first_seen$country), "malawi", "south africa", 
                           "zimbabwe", "zambia"))

#introductions

introductions <- readRDS(here("data", "introductions.RDS")) |> 
  mutate_at(c("year_min", "year_max", "year_med"), ~floor(.))

introductions_dist <- readRDS(here("data", "introductions.RDS")) 



