# Simulate outbreak spread from new introductions -------------------------------

## Initialize --------------------------------------------------------------------
## ---- initialize
library(tidyverse)
library(here)
library(furrr)
library(qs2)

# plan(multisession, workers = 4)
plan(list(tweak(multisession, workers = 2), tweak(multisession, workers = I(8))))

# which_obs <- "all_data"
source(here("analysis","00_functions_settings.R"), local = T)


if(!exists("nsamp")){
  nsamp <- 1000
}
if(!exists("yrs")){
  yrs <- 30
}
if(!exists("downstream_outbreaks")){
  downstream_outbreaks <- T
}
if(!exists("secondary_outbreaks")){
  secondary_outbreaks <- T
}
if (!dir.exists(sim_dir)) {dir.create(sim_dir)}

if(!exists("sim_all")){
  sim_all <- T
}

show_progress <- F

### Read data----
source(here("analysis","02_import_data.R"), local = T) 
source(here("analysis","04_hmm_setup.R"), local = T) 


#Bring in connectivity measure (used as weight/metric for clustering) from model output
xi_files <- list.files(paste0(xi_dir,"/draws/"))
xi_run_dates <- str_extract_all(xi_files, "(?<=draws_).*") |> list_c() |> as_date()
xi_most_recent <- xi_files[which(xi_run_dates==max(xi_run_dates))]
connectivity <- qs2::qs_read(here::here(xi_dir,"draws",xi_most_recent))

model_files <- list.files(paste0(out_dir))
model_run_dates <- str_extract_all(model_files, "(?<=objs_).*") |> list_c() |> as_date()
most_recent_run <- str_c("model_objs_", model_run_dates[which(model_run_dates==max(model_run_dates))])
model <- qs2::qs_read(here::here(out_dir,most_recent_run))


### Setup for simulation

## We want to simulate spread based on introductions to different locations
## Start with point introductions by country, 2 strains

#initialize simulation 
M <- length(countries_subset)
n_time <- yrs    # number of years
epsilon <- 1e-7 #small number


#single introduction of a single strain
if(test_subset){
  seed_countries <- seed_countries[seed_countries %in% countries_subset]
}

if(sim_all){
  seed_countries <- countries_subset
}

#Distribution of observed cases by country over time
case_dists <- cases_lineages |> 
  select(country, cases, year) |> 
  #more recent surveillance data
  filter(!is.na(cases),
         year>1989) |> 
  group_by(country) |> 
  mutate(max = max(cases)) |> 
  ungroup() |>
  arrange(country, year) |> 
  filter((cases>0 & max>0) | (cases==0 & max==0)) |> 
  select(-max) |>
  distinct() |> 
  rowwise() |> 
  mutate(country_id = which(countries_subset == country)) |> 
  ungroup()

### Sample parameters----

#Will be using connectivity, persistence, & case transformation from model draws
param_draws <- append(list(pull_params("xi",connectivity)),
                    purrr::map(c("delta", "eta"), pull_params, model$model_draws)) |> 
  set_names(c("xi", "delta", "eta"))
runs <- model$metadata$iter_sampling*stan_chains

## Simulate----

### Downstream----

library(tictoc)
if(downstream_outbreaks){
  tic("downstream")
  sims_downstream <- seed_countries |>
    purrr::set_names() |>
    furrr::future_map(~{
      furrr::future_map(1:nsamp, sim_spread_func,
                        loc_int = which(countries_subset==.x),
                        .options = furrr_options(seed = TRUE))}, 
      .options = furrr_options(seed = TRUE), .progress = show_progress)
  downstream_pred <- furrr::future_map(sims_downstream, ~{
    d <- bind_rows(.x, .id = "draw")
  }, 
  .options = furrr_options(seed = TRUE), .progress = show_progress) |>
    bind_rows(.id = "seed_country")
  toc()
}

### Secondary----

if(secondary_outbreaks){
  tic("secondary")
  sims_secondary <- countries_subset |>
    purrr::set_names() |>
    furrr::future_map(~{
      furrr::future_map(1:nsamp, sim_spread_func,
                        loc_int = which(countries_subset==.x),
                        downstream = F,
                        .options = furrr_options(seed = TRUE))},
      .options = furrr_options(seed = TRUE), .progress = show_progress)
  secondary_pred <- furrr::future_map(sims_secondary, ~{
    d <- bind_rows(.x, .id = "draw") |>
      secondary_func()}, 
    .options = furrr_options(seed = TRUE), .progress = show_progress) |>
    bind_rows(.id = "seed_country")
  toc()
}

## Save data ----

if(downstream_outbreaks){
  saveRDS(sims_downstream, glue::glue("{sim_dir}/sim_downstream.rds"))
  saveRDS(downstream_pred, glue::glue("{sim_dir}/downstream_pred.rds"))
}

if(secondary_outbreaks){
  saveRDS(secondary_pred, glue::glue("{sim_dir}/secondary_pred.rds"))
  saveRDS(sims_secondary, glue::glue("{sim_dir}/sim_secondary.rds"))
}


