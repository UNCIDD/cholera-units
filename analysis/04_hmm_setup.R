#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 04_hmm_setup.R
# Purpose: Data setup to run HMM model in Stan
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(cmdstanr)
library(here)

do_plots <- T  # whether to run intermediate plots

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

# sourced from "00_functions_settings.R" if not updated
if (!dir.exists(out_dir)) {dir.create(out_dir)}

set.seed(seed)

#-----------------------------------------------------------#
# Prep data----
#-----------------------------------------------------------#

### A. Import data ----

#import
## ---- import
source(here("analysis","02_import_data.R"), local = T)

## ---- data_prep
cases_lineages <- cases_lineages |> 
  select(-subregion)


### B. Indices ----
## ---- parameters

# indices for observations
M <- length(countries_subset)
V <- length(te_subset[te_subset!="S"]) # number of strains
n_time <- (max(cases_lineages$year)-min(cases_lineages$year))+1  # number of time slices
N <- M * n_time  # Max number of space*time observations

# small number
epsilon <- 1e-7

#country/year shell
mt_shell <- tibble(country = rep(1:M, each=n_time),
                   time = rep(1:n_time, M)) |> 
  arrange(country, time) |> 
  mutate(location = str_c("loc-", country)) |> 
  select(location, time)

#country/strain shell
m_v_shell <- tibble(m = rep(1:M, each = V),
                    v = rep(1:V, times = M))


### C. Connectivity ----

## ---- spatial
#create spatial mappings
spatial_objs <- create_connections(M, n_time)

#pull centroids for each country in subset & create distance matrix
distances <- st_centroid(africa_sf  |> 
                           filter(country %in% countries_subset)) |> 
  st_distance() |> 
  as.data.frame() |> 
  set_names((africa_sf  |> 
               filter(country %in% countries_subset))$country) |> 
  purrr::map_dfc( ~ {units::set_units(.x, km)}) 
distances <- distances |> 
  mutate(c = names(distances)) |> 
  mutate(m = purrr::map_int(c, ~which(countries_subset == .x))) |> 
  arrange(m) |> 
  dplyr::select(all_of(countries_subset)) |> 
  as.matrix()
row.names(distances) <- 1:length(countries_subset)
colnames(distances) <- 1:length(countries_subset)

#distance matrix
dist_mat <- matrix(distances, nrow = M, ncol = M)
#population sizes
pop_sizes <- purrr::map_vec(countries_subset, ~africa_sf$pop[africa_sf$country==.x])


### D. Cases ----
## ---- cases

#pull out annual case counts
cases <- cases_lineages |> 
  select(country, year, cases) |> 
  distinct() |> 
  mutate(time = year-min_year+1) |> 
  mutate(location = purrr::map_chr(country, ~as.character(which(countries_subset == .x)))) |> 
  arrange(location, time) |> 
  mutate(location = str_c("loc-", location)) |> 
  select(location, cases, time)

#cases with 1 row for each country & year
cases_long <- mt_shell |> 
  left_join(cases) |> 
  mutate(cases = replace_na(cases, 0))

### E. Samples ----
## ---- samples

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

if(downsample_random){
  set.seed(downsample_seed)
  partitioned_samples <- readRDS(here::here("data/downsample_groups.RDS"))
  samples <- samples |> 
    left_join(partitioned_samples, by = c("location", "time", "te", "samples")) |> 
    mutate(samples = if_else(.groups == downsample_run_num, 0, samples)) |> 
    select(-.groups)
}

samples <- samples |>  
  pivot_wider(names_from = te, values_from = samples) |> 
  select(location, time, all_of(te_subset[te_subset!="S"]))

#samples with 1 row for each country & year
samples_long <- mt_shell |> 
  left_join(samples, by = c("location", "time")) |> 
  mutate_at(te_subset[te_subset!="S"], ~replace_na(., 0))

### F. Introductions ----

## ---- intros
intro_m_v <- introductions |> 
  filter(te %in% te_subset[te_subset!="S"],
         country %in% countries_subset) |> 
  mutate(t = year_med-min(min_year)+1,
         m = purrr::map_dbl(country, ~which(countries_subset == .x)),
         v = purrr::map_int(te, ~which(te_subset == .x))) |> 
  select(m, t, v) |> 
  arrange(as.numeric(v), m, t) |> 
  group_by(v, m) |> 
  filter(t == min(t)) |> 
  ungroup() |> 
  right_join(m_v_shell) |> 
  pivot_wider(id_cols = m, names_from = v, values_from = t) |> 
  select(m, str_c(1:length(te_subset[te_subset!="S"]))) |> 
  arrange(m)

intros_mat <- as.matrix(intro_m_v |> 
                     mutate_all(~if_else(is.na(.) | . > n_time | .==1, 0, .)) |> 
                     select(-m))
row.names(intros_mat) <- intro_m_v$m


### G. Prep for Stan ----
## ---- stan_prep

#introductions, using mean timing of introduction
N_gamma <- length(which(colSums(intros_mat)>0))
map_v_to_gamma <- vector("integer", length = V)
g = 0
for(v in 1:V){
  if(v %in% which(colSums(intros_mat)>0)){
    g = g+1
    map_v_to_gamma[v] <- g
  }else{
    map_v_to_gamma[v] <- 0
  }
  
}

# probabilities by year for each lineage
intro_yr_p <- purrr::map(te_subset[-length(te_subset)], ~{
  d <- introductions |> 
    select(te, starts_with("year")) |> 
    filter(te==.x) |> 
    distinct() |> 
    mutate_at(c("year_min","year_med", "year_max"), ~(.-min(cases_lineages$year)+1))
  yrs <- max(cases_lineages$year)-min(cases_lineages$year)+1
  mu <- d$year_med
  sigma <- ifelse(d$year_min!=d$year_med, (d$year_med-d$year_min)/2, 0.1)
  cprobs <- pnorm(seq(1,yrs), mean = mu, sd = sigma)
  probs <- cprobs
  for(i in 2:length(probs)){
    probs[i] <- cprobs[i]-cprobs[i-1]
  }
  p <- probs[(min_year-min(cases_lineages$year)+1):length(probs)]
  return(probs)
})
# probabilities by country for each lineage
intro_loc_p <- purrr::map(te_subset[-length(te_subset)], ~{
  d <- data.frame(country = all_countries) |> 
    left_join(introductions |> 
                filter(te==.x) |> 
                select(country, prob_dist) |> 
                distinct(), by = "country") |> 
    filter(country %in% countries_subset) |> 
    mutate(p = if_else(is.na(prob_dist),0,prob_dist)) |> 
    pull(p)
}) 

N_loc_gamma <- max(purrr::map_int(intro_loc_p,~length(.x[.x>1e-3])))
N_yr_gamma <- max(purrr::map_int(intro_yr_p,~length(.x[.x>1e-3])))

map_gamma_loc <- purrr::map(intro_loc_p, ~{
  names(.x) <- 1:length(.x)
  g <- sort(.x, decreasing=T)[1:N_loc_gamma]
  h <- sort(as.numeric(names(g)))
  a <- vector("integer",length = length(.x))
  i <- 0
  for(m in 1:M){
    if(m %in% h){
      i = i+1
      a[m] <- i
    }else{
      a[m] <- 0
    }
  }
  return(a)
}) |> 
  bind_cols(.name_repair = "unique_quiet") |> 
  as.matrix()
map_gamma_loc[,which(map_v_to_gamma==0)] <- 0
map_gamma_yr <- purrr::map(intro_yr_p, ~{
  names(.x) <- 1:length(.x)
  g <- sort(.x, decreasing=T)[1:N_yr_gamma]
  h <- sort(as.numeric(names(g)))
  a <- vector("integer",length = length(.x))
  i <- 0
  for(t in 1:n_time){
    if(t %in% h){
      i = i+1
      a[t] <- i
    }else{
      a[t] <- 0
    }
  }
  return(a)
}) |> 
  bind_cols(.name_repair = "unique_quiet") |> 
  as.matrix()
map_gamma_yr[,which(map_v_to_gamma==0)] <- 0

intro_loc_p_mat <- intro_loc_p |> 
  bind_cols(.name_repair = "unique_quiet") |> 
  as.matrix()
intro_yr_p_mat <- intro_yr_p |> 
  bind_cols(.name_repair = "unique_quiet") |> 
  as.matrix()


gamma_loc <- matrix(1e-3, nrow = N_loc_gamma, ncol = N_gamma)
gamma_yr <- matrix(1e-3, nrow = N_yr_gamma, ncol = N_gamma)

#priors for location & year introduction probabilities
for(v in 1:V){
  if(map_v_to_gamma[v]>0){
    gamma_loc[,map_v_to_gamma[v]] <- intro_loc_p_mat[which(map_gamma_loc[,v]>0),v]
    gamma_yr[,map_v_to_gamma[v]] <- intro_yr_p_mat[which(map_gamma_yr[,v]>0),v]
  }
}

gamma_loc[which(gamma_loc<1e-3)] <- 1e-3
gamma_yr[which(gamma_yr<1e-3)] <- 1e-3

# Convert cases & samples to matrices
y_obs <- as.matrix(samples_long |> select(-location, -time))


# Prep sequences for Stan
N_samp <- nrow(y_obs)
N_sequences <- cases_lineages |> 
  select(country, year, tot_samples) |> 
  distinct() |> 
  mutate(time = year-min_year+1) |> 
  mutate(location = purrr::map_chr(country, ~as.character(which(countries_subset == .x)))) |> 
  arrange(location, time) |> 
  mutate(location = str_c("loc-", location)) |> 
  select(location, time, tot_samples) |> 
  arrange(location, time, tot_samples)

N_sequences <- mt_shell |> 
  left_join(N_sequences, by = c("location", "time")) |> 
  mutate(tot_samples = replace_na(tot_samples, 0)) 



