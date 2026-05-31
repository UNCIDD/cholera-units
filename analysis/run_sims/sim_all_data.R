#--------------------------------------------#
#Run outbreak simulations - All data
#--------------------------------------------#

library(tidyverse)
library(here)
source(here("analysis", "00_functions_settings.R"), local = T)

# sampling_options <- c("all_data", "post_1990", "post_2000")

## Global params
stan_chains <- 4
which_obs <- "all_data"
sampling <- which_obs
min_year <- get_min_year(which_obs)

#simulation settings
nsamp <- 2000
years <- 30


#Run simulation

if(!exists("sim_all")){
  sim_all <- T
}

if(sim_all){
  downstream_outbreaks <- T
  secondary_outbreaks <- T
  
  source(here("analysis", "09_sim_outbreaks_poststan.R"), local = T)
}
