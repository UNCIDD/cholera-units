# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 05_run_hmm.R
# Purpose: HMM model execution
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(cmdstanr)
library(here)

# test_subset <- T
# min_year <- 1970

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

# sourced from "00_functions_settings.R" if not updated
if (!dir.exists(out_dir)) {dir.create(out_dir)}
if (!dir.exists(xi_dir)) {dir.create(xi_dir)}
if (!dir.exists(str_c(xi_dir,"/draws"))) {dir.create(str_c(xi_dir,"/draws"))}

set.seed(seed)

#-----------------------------------------------------------#
# Import----
#-----------------------------------------------------------#

#import data that is prepped for Stan run
## ---- import 

if(simulate){
  out_dir <- str_c(here(), "/generated_data/stan/simulated_data")
  if (!dir.exists(out_dir)) {dir.create(out_dir)}
  source(here("analysis/hmm_sim_test/04b_hmm_setup_simulate.R"), local = T)
}else{
  source(here("analysis","04_hmm_setup.R"), local = T)
}

#-----------------------------------------------------------#
# Stan data----
#-----------------------------------------------------------#

## ---- stan_data

stan_data <- list(
  #-Indices-#
  M = M, #countries
  V = V, #strains
  `T` = n_time, #years
  K = 2, #states

  N = N, #time*location
  NV = N*V, #time*location*strain
  N_sequences = N_sequences |> 
    select(tot_samples) |> pull(), #total number of observed sequences across all strains (space x time)
  N_gamma = N_gamma, #number of space x time units with introductions
  N_loc_gamma = N_loc_gamma,
  N_yr_gamma = N_yr_gamma,
  
  #-Observations-#
  y = y_obs, #number of observed sequences of each strain
  cases = cases_long |> 
    pull(cases), #observed cases by country & year
  pop = pop_sizes, #population sizes by country
  
  #-Mappings-#
  map_to_location = spatial_objs$map_to_location,
  map_to_time = spatial_objs$map_to_time,
  map_to_n = t(spatial_objs$map_to_n),
  map_from_edge = spatial_objs$map_from_edge,
  N_edges = nrow(spatial_objs$adj_nodes),
  node1_by_dest = spatial_objs$node1_by_dest,
  dist_by_dest = map_dbl(1:nrow(spatial_objs$adj_nodes),
                         ~dist_mat[spatial_objs$adj_nodes$node1[.], spatial_objs$adj_nodes$node2[.]]),
  start_by_dest = spatial_objs$start_by_dest,
  end_by_dest = spatial_objs$end_by_dest,
  map_v_to_gamma = map_v_to_gamma,
  map_gamma_loc = map_gamma_loc,
  map_gamma_yr = map_gamma_yr,
  
  #-Introductions-#
  intro_times0 = intros_mat,
  
  epsilon = epsilon,
  
  #-Priors-#
  mu_delta = -3,
  sd_delta = 0.5,
  p_gamma_loc = t(gamma_loc), #mean parameter for beta distribution
  p_gamma_yr = t(gamma_yr), #mean parameter for beta distribution
  l_gamma = 500, #count parameter for beta distribution
  
  sd_lambda = 5, 
  
  #spatial effect
  p_w = 0.7,
  mu_w = c(-3,0), 
  sd_w = c(1,1),

  #connectivity measures
  mu_kappa = 0.5, #gravity constant
  sd_kappa = 0.1, #gravity constant
  mu_tau = c(0.35, 0.45), #population parameter
  sd_tau = 500, #population parameter
  mu_zeta = 2.25, #distance parameter
  sd_zeta = 0.1, #distance parameter
  
  #case parameters
  mu_e = 0.8, #under-reporting
  sd_e = 100, #under-reporting
  mu_eta = 0.45, #case transformation
  sd_eta = 500 #case transformation
  
)


#-----------------------------------------------------------#
# Run stan model----
#-----------------------------------------------------------#

## ---- stan_model

if(run_model){
 
  #specify model file 
  stan_model <- cmdstan_model(here("analysis","stan","hmm_fb_model_connectivity.stan"), stanc_options = list("O1"))
  
  #fit the model
  model_fit <- stan_model$sample(
    data = stan_data,
    seed = seed,
    init = 0.1, 
    chains = stan_chains,
    parallel_chains = stan_chains,
    iter_warmup = stan_iter_warmup,
    iter_sampling = stan_iter_sample,
    max_treedepth = 14L,
    adapt_delta = .9,
    refresh = 250,
    save_warmup = save_warmup,
    show_messages = T)
  
  # parameter names
  vars <- c("lp__", "eta", "kappa", "zeta", "tau", "log_delta", "delta", "w", 
            "rho","alpha","pi_0","beta","gamma","lambda","lambda_star","prev_lambda",
            "pred_cases","zstar","e_i")
  
  model_objs <- list()
  model_objs$metadata <- model_fit$metadata()
  model_objs$diagnostic_summary <- model_fit$diagnostic_summary()
  model_objs$model_draws <- model_fit$draws(inc_warmup = F, format = "draws_df", variables = vars)
  model_objs$var_summary <- model_fit$summary(variables = c(vars, "xi"))
  model_objs$neffs <- bayesplot::neff_ratio(model_fit, pars = vars)
  model_objs$rhats <- bayesplot::rhat(model_fit, pars = vars)
  
  xi_model_draws <- model_fit$draws(variables = "xi", inc_warmup = F, format = "draws_df")
  
  #save model object
  qs2::qs_save(model_objs, glue::glue("{out_dir}/model_objs_{Sys.Date()}"))
  qs2::qs_save(xi_model_draws, glue::glue("{xi_dir}/draws/xi_model_draws_{Sys.Date()}"))
  
}




