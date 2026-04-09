#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 03_create_partitions.R
# Purpose: Create partitions for downsampling sensitivity analysis
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(here)
library(groupdata2)


seed <- 22

#set.seed(seed)

### Source functions & Read data----
source(here("analysis","00_functions_settings.R"), local = T)
source(here("analysis","02_import_data.R"), local = T)

#-----------------------------------------------------------#
# Create downsample groups----
#-----------------------------------------------------------#

#pull out information on sequenced samples
samples <- cases_lineages |> 
  select(country, year, te, samples) |> 
  distinct() |> 
  #calculate time since first obs
  mutate(time = year-min_year+1) |> 
  #remove sporadic outbreak calls
  filter(!is.na(te) & te!="S") |> 
  #location & strain indicators
  mutate(location = purrr::map_chr(country, ~as.character(which(countries_subset == .x))),
         strain = purrr::map_chr(te, ~as.character(which(te_subset == .x)))) |> 
  arrange(location, time) |> 
  mutate(location = str_c("loc-", location)) |> 
  select(location, time, samples, te, strain) |> 
  arrange(te, location, time, samples) |> select(-strain) 

# Downsample based on specified downsampling fraction (see "00_functions_settings.R")
downsample_groups <- groupdata2::group(samples, randomize = T, n=1/(1-sample_fraction))

#-----------------------------------------------------------#
# Save out ----
#-----------------------------------------------------------#

saveRDS(downsample_groups, here::here("data/downsample_groups.rds"))

