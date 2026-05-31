#--------------------------------------------#
#Run outbreak simulations - All data, downstream
#--------------------------------------------#

library(tidyverse)
library(here)
source(here("analysis", "00_functions_settings.R"), local = T)

## Global params
sim_all <- F
source(here("analysis", "run_sims", "sim_all_data.R"), local = T)

#simulation settings
downstream_outbreaks <- T
secondary_outbreaks <- F


#Run simulation

source(here("analysis/", "09_sim_outbreaks_poststan.R"), local = T)


