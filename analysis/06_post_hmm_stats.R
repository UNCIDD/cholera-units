#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 06_post_hmm_stats.R
# Purpose: Evaluating inferences from HMM
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

## ---- initialize

library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(spData) 
library(sf)
library(posterior)
library(here)
library(cowplot)

options(scipen = 999999999)

# run_model <- F
# test_subset <- T
# min_year <- 1990


# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

if (!dir.exists(out_dir)) {dir.create(out_dir)}
if (!dir.exists(xi_dir)) {dir.create(xi_dir)}

set.seed(seed)


#-----------------------------------------------------------#
# Run model or load prior run----
#-----------------------------------------------------------# 

## ---- run_model 

if(run_model){
  source(here("analysis","05_run_hmm.R"), local = T)
}else{
  print("Loading prior model run")
}

## ---- load_model

if(run_model){
  print("New model run")
  model_seed <- model_objs$metadata$seed
}else{
  #run HMM data prep
  source(here("analysis","04_hmm_setup.R"), local = T)
  #load model output
  model_runs <- file.mtime(paste0(out_dir,"/",list.files(out_dir)))
  model_objs <- qs2::qs_read(here::here(out_dir,list.files(out_dir)[which(model_runs==max(model_runs))]))
  model_seed <- model_objs$metadata$seed
  seed_match <- seed==model_seed
  if(!seed_match){
    seed <- model_seed
    source(here("analysis","04_hmm_setup.R"), local = T)
  }
}

#-----------------------------------------------------------#
# Stan output----
#-----------------------------------------------------------# 

## A. Model Diagnostics ----

## ---- init

# parameter names
pars <- c("lp__", "eta", "kappa", "zeta", "tau[1]", "tau[2]", "delta") 

color_scheme_set("blue")

## ---- diagnostic_summ
model_objs$diagnostic_summary

#Neff ratio
model_objs$neffs[pars[-1]]

#Neffs & rhats for main parameters
mcmc_neff(model_objs$neffs[pars[-1]]) + yaxis_text(hjust = 1)
mcmc_rhat(model_objs$rhats[pars[-1]]) + yaxis_text(hjust = 1)

#Neffs & rhats for omega parameters
w_pars <- model_objs$metadata$model_params[str_which(model_objs$metadata$model_params, "^w\\[")]

w_rhats <- model_objs$rhats[w_pars]
mcmc_rhat(w_rhats) + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank())

if(length(names(which(w_rhats >=1.1)))>0){
  print(c("Rhat >= 1.1:", names(purrr::map(which(w_rhats >=1.1), ~{
    str_c("edge ", .x, ": ", 
          str_c(countries_subset[which(spatial_objs$map_from_edge==.x, arr.ind = T)], 
                collapse = ", "))
  }) |> set_names())))
}else{
  print("All Rhats <1.1")
}

w_neffs <- model_objs$neffs[w_pars]
mcmc_neff(w_neffs) + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank())

if(length(names(which(w_neffs <=0.1)))>0){
  print(c("Neff ratio <= 0.1:", names(purrr::map(which(w_neffs <=0.1), ~{
    str_c("edge ", .x, ": ", 
          str_c(countries_subset[which(spatial_objs$map_from_edge==.x, arr.ind = T)], 
                collapse = ", "))
  }) |> set_names())))
}else{
  print("All Neff Ratios >0.1")
}

#rhats for pi parameters
pi_pars <- model_objs$metadata$model_params[str_which(model_objs$metadata$model_params, "^pi_0\\[")]

pi_rhat <- model_objs$rhats[pi_pars]
mcmc_rhat(pi_rhat) + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank())

#Autocorrelation
mcmc_acf(model_objs$model_draws, pars = pars)

# Trace plots
traceplot(obj = model_objs$model_draws, "lp_")
traceplot(obj = model_objs$model_draws, "kappa")
traceplot(obj = model_objs$model_draws, "zeta")
traceplot(obj = model_objs$model_draws, "^eta") 
traceplot(obj = model_objs$model_draws, "e_i")
traceplot(obj = model_objs$model_draws, "log_delta")
traceplot(obj = model_objs$model_draws, "^delta")
traceplot(obj = model_objs$model_draws, "tau")

# Parameter histograms
mcmc_hist(model_objs$model_draws, pars = "kappa")
mcmc_hist(model_objs$model_draws, pars = "zeta")
mcmc_hist(model_objs$model_draws, pars = "log_delta")
mcmc_hist(model_objs$model_draws, pars = "delta")
mcmc_hist(model_objs$model_draws, pars = "eta")
mcmc_hist(model_objs$model_draws, pars = "tau[1]")
mcmc_hist(model_objs$model_draws, pars = "tau[2]")

# Pair plots
mcmc_pairs(model_objs$model_draws, pars = pars,
           off_diag_args = list(size = 0.75))

## B. Parameter summaries ----
## ---- param_summ

purrr::map(pars[-1], model_fit_summary) |> 
  set_names(pars[-1]) |> 
  bind_rows(.id = "Parameter") |> select(-true) |> 
  print(n=length(pars[-1]))

purrr::map(str_c("e_i[", 1:M, "]"), model_fit_summary) |> 
  set_names(str_c("e_i[", 1:M, "]")) |> 
  bind_rows(.id = "Parameter") |> select(-true) |> 
  print(n=M)

## C. Datasets & Visualizations ----

### Predicted probabilities & prevalence----
## ---- est_dfs
# Manipulate observed data to match stan output for comparison
all_samples <- cases_lineages |> 
  select(country, year, te, samples) |> 
  distinct() |> 
  mutate(time = year-min_year+1) |> 
  filter(!is.na(te)) |> 
  mutate(location = purrr::map_chr(country, ~as.character(which(countries_subset == .x))),
         strain = purrr::map_chr(te, ~as.character(which(te_subset == .x)))) |> 
  arrange(location, time) |> 
  mutate(location = str_c("loc-", location)) |> 
  select(location, time, samples, te, strain) |> 
  arrange(te, location, time, samples) |> select(-strain) |>  
  pivot_wider(names_from = te, values_from = samples) |> 
  select(location, time, all_of(te_subset))

all_samples_long <- mt_shell |> 
  left_join(all_samples, by = c("location", "time")) |> 
  mutate_at(te_subset, ~replace_na(., 0))

te_lookup <- tibble(from = te_subset, to = 1:length(te_subset))

df_obs <- as.data.frame(y_obs) |> 
  mutate(location = str_c("loc-", spatial_objs$map_to_location[row_number()]),
         time = spatial_objs$map_to_time[row_number()]) |> 
  pivot_longer(col = -c(location, time),
               names_to = "te",
               values_to = "samples") |> 
  mutate(strain = recode_values(te, from = te_lookup$from, to = te_lookup$to)) |> 
  select(location, time, strain, samples, te)

# Manipulate stan output for comparison
pred_rho <- poststan_df("rho", state = F)
pred_alpha <- poststan_df("^alpha", state = T)
pred_pi <- poststan_df("pi_0", state = F) |> 
  mutate(location = str_c("loc-", ind),
         time = 1) |> 
  translate_cols() |> 
  select(country, te, year, mean) 
pred_beta <- poststan_df("^beta")
pred_gamma <- poststan_df("gamma", state = F)
pred_lambda <- poststan_df("^lambda\\[", state = F)
pred_lambda_star <- poststan_df("lambda_star", state = F)

#generated quantities
pred_prev_lambda <- poststan_df("prev_lambda", state = F)
pred_cases <- poststan_df("pred_cases", state = F)
pred_z <- poststan_df("zstar", state = F) |> 
  mutate(median = median-1)

cases_obs_pred <- cases_long |> 
  left_join(pred_cases |> 
              select(location, time, strain, median, mean), by = c("location", "time")) |> 
  translate_cols() 

## ---- plot_rho
rho_plot <- postest_plot(pred_rho, title = "Sub-Lineage Presence") + 
  guides(color = guide_legend(nrow = 2), fill =guide_legend(nrow = 2)) +
  theme(text = element_text(size = 8), legend.position = "none")
rho_plot
## ---- plot_alpha
postest_plot(pred_alpha, title = "Alpha")
## ---- plot_beta
postest_plot(pred_beta, title = "Beta")
## ---- plot_pi
pi_plot <- ggplot() +
  geom_point(aes(x = te, y = mean, color = te), data = pred_pi) +
  facet_wrap(~country, ncol = 2, scales = "free_y") +
  plot_colors +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  plottheme
pi_plot

## ---- plot_lstar
postest_plot(pred_lambda_star, title = "Prevalence")  +
  theme(text = element_text(size = 8), legend.position = "none")
## ---- plot_gamma
postest_plot(pred_gamma, title = "Gamma") + 
  theme(text = element_text(size = 8), legend.position = "none")
## ---- plot_lambda
postest_plot(pred_lambda, title = "Lambda", show_cases = F) + 
  theme(text = element_text(size = 8), legend.position = "none")
## ---- plot_prev
prev_plot <- postest_plot(pred_prev_lambda, title = "Lambda Prev")  +
  theme(text = element_text(size = 8), legend.position = "none")
prev_plot

## ---- plot_zstar
zstar_plot <- pred_z |> 
  translate_cols() |> 
  select(country, year, te, median) |> 
  rename(pa = median) |> 
  filter(pa==1) |> 
  mutate(source = "Predicted",
         pa = pa*(as.numeric(te)*20)) |> 
  bind_rows(translate_cols(df_obs) |> 
              mutate(pa = as.double(samples > 0), 
                     source = "Present") %>% 
              dplyr::select(country, year, te, samples, pa, source) %>% 
              filter(pa == 1) %>% 
              distinct() |>  
              arrange(country, year, te, source) |> 
              mutate(pa = pa*(as.numeric(te)*20+5)) |> 
              distinct() |> 
              select(country, year, te, pa, source)) |> 
  mutate(source = factor(source)) |> 
  ggplot() +
  geom_point(aes(x = year, y = pa, color = te, shape = source), size = 1, alpha = 0.8,
             position = position_dodge(width = 0.5)) +
  facet_wrap(~country, ncol = 2) + 
  plot_colors +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_shape_manual(values = c("Predicted" = 8, "Present"=1)) +
  labs(title = "Predicted Presence",
       y = "Presence",
       color = "TE",
       shape = "TE") +
  theme_classic() +
  plottheme +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 6))
zstar_plot

### Cases----
## ---- pred_cases

prev_cases_plot <- ggplot() +
  geom_bar(data = cases_obs_pred |> distinct(country, year, cases), 
           aes(x = year, y = cases), stat = "identity", fill = "white") +
  geom_bar(data = cases_obs_pred |> distinct(country, year, te, mean),
           aes(x = year, y = mean, fill = te),
           show.legend = T, stat = "identity", position = "stack") +
  facet_wrap(~country, ncol = 2, scales = "free_y") +
  guides(fill=guide_legend(nrow=2, title.position = "left")) +
  plot_fillcolors +
  theme_bw() +
  ylab("Cases") +
  theme(legend.position = "bottom",
        legend.title = element_blank())
prev_cases_plot


### Connectivity----

## ---- map_prep
pred_w <- model_objs$var_summary |> 
  filter(str_detect(variable, "^w")) |> 
  mutate(node_n = str_extract(variable, "(?<=\\[)[0-9]+"))
pred_xi <- model_objs$var_summary |> 
  filter(str_detect(variable, "^xi")) |>
  mutate(node_n = str_extract(variable, "(?<=\\[)[0-9]+"))

centroids <- st_centroid(africa_sf  |> 
                           filter(country %in% countries_subset)) 

nodes <- st_centroid(africa_sf  |> 
                       filter(country %in% countries_subset)) |> 
  dplyr::select(country, geometry) |> 
  mutate(long = purrr::map_dbl(geometry, ~.x[1]),
         lat = purrr::map_dbl(geometry, ~.x[2])) 
st_geometry(nodes) <- NULL

edges <- spatial_objs$adj_nodes |> 
  mutate(country1 = purrr::map_chr(node1, ~ {countries_subset[[.x]]}),
         country2 = purrr::map_chr(node2, ~ {countries_subset[[.x]]})) |> 
  left_join(nodes |> rename(xstart = long, ystart = lat), by = c("country1" = "country")) |> 
  left_join(nodes |> rename(xend = long, yend = lat), by = c("country2" = "country")) |> 
  left_join(pred_w |> dplyr::select(mean_w = mean, node_n) |> 
              mutate(node_n = as.integer(node_n)), 
            by = c("n"="node_n")) |> 
  left_join(pred_xi |> dplyr::select(mean_xi = mean, node_n) |> 
              mutate(node_n = as.integer(node_n)), 
            by = c("n"="node_n")) |> 
  filter(node1!=node2)


#Plot w
## ---- map_w
w_map_plot <- ggplot() +
  geom_sf(data = full_africa_sf, fill = "transparent", color = "grey40") + 
  geom_sf(data = africa_sf, aes(fill = pop), alpha = 0.9) + 
  geom_curve(data = (edges |> filter(!is.na(mean_w)) |> distinct() |> arrange(mean_w) |> 
                       filter(mean_w>=(-0.5))),
             aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                 color = (mean_w)),
             angle = 145,
             arrow = arrow(length = unit(0.1,"cm"))) +
  map_options_func(1, "Omega pred", exp = F, color_h = "#CC0000", color_l = "white", color_m = "#FF9999") 
w_map_plot


#Plot Xi
## ---- map_xi

xi_map_plot <- ggplot() + 
  geom_sf(data = full_africa_sf, fill = "transparent", color = "grey40") + 
  geom_sf(data = africa_sf, aes(fill = pop), alpha = 0.9) + 
  geom_curve(aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                 linewidth = log_xi_norm, 
                 color = log_xi_norm), 
             data = edges |> filter(!is.na(mean_xi)) |> 
               distinct() |> arrange(mean_xi) |> 
               filter(mean_xi>=exp(-7)) |> 
               mutate(xi_norm = scale(mean_xi),
                      log_xi_norm = scale(log(mean_xi))),
             angle = 145,
             arrow = arrow(length = unit(0.1,"cm"))) +
  map_options_func(0.5, "Xi", exp = F, color_h = "#CC0000", color_l = "steelblue3", color_m = "white")
xi_map_plot


### Pred & obs lineages----
## ---- pred_obs

lineages_toplot <- cases_lineages |> 
  select(country, year, te) |> 
  mutate(source = "Observed") |> 
  bind_rows(pred_z |> 
              translate_cols() |> 
              filter(median==1) |> 
              select(country, year, te) |> 
              mutate(source = "Inferred")) |> 
  filter(!is.na(te)) |> 
  arrange(country, year, te) |> 
  distinct()

cases_transmission_plot <- cases_lineages_plot(case_data = cases_lineages |>
                                                      filter(!is.na(country)) |>
                                                      mutate(cases = as.numeric(cases)) |>
                                                      distinct(), 
                                                    lineage_data = lineages_toplot)
cases_transmission_plot 


## Save data----
## ---- save_clust

saveRDS(edges, glue::glue("{xi_dir}/edges.rds"))
saveRDS(nodes, glue::glue("{xi_dir}/nodes.rds"))
saveRDS(countries_subset, glue::glue("{xi_dir}/countries.rds"))
saveRDS(lineages_toplot, glue::glue("{xi_dir}/cases_lineages_assign.rds"))
saveRDS(cases_obs_pred, glue::glue("{xi_dir}/cases_prev_pred.rds"))
saveRDS(pred_cases, glue::glue("{xi_dir}/gen_pred_cases.rds"))
saveRDS(pred_xi, glue::glue("{xi_dir}/pred_xi.rds"))
saveRDS(pred_w, glue::glue("{xi_dir}/pred_w.rds"))
saveRDS(pred_rho, glue::glue("{xi_dir}/pred_rho.rds"))
saveRDS(pred_prev_lambda, glue::glue("{xi_dir}/pred_prev_lambda.rds"))
saveRDS(pred_z, glue::glue("{xi_dir}/pred_z_star.rds"))
saveRDS(model_seed, glue::glue("{xi_dir}/model_run_seed.rds"))

if(downsample_random){
  saveRDS(y_obs, glue::glue("{xi_dir}/downsampled_y.rds"))
}

if(save_plots){
  
  if (!dir.exists(fig_dir)) {dir.create(fig_dir)}
  
  ggsave(plot = rho_plot, 
         filename = glue::glue("{fig_dir}/rho_plot.png"), height = 17)
  ggsave(plot = cases_transmission_plot, 
         filename = glue::glue("{fig_dir}/cases_transmission_plot.png"), 
         width = 17, height = 9)
  ggsave(plot = prev_plot, 
         filename = glue::glue("{fig_dir}/prev_plot.png"), height = 17)
  ggsave(plot = zstar_plot, 
         filename = glue::glue("{fig_dir}/zstar_plot.png"), height = 17)
  ggsave(plot = w_map_plot, 
         filename = glue::glue("{fig_dir}/w_map_plot.png"))
  ggsave(plot = xi_map_plot, 
         filename = glue::glue("{fig_dir}/xi_map_plot.png"))
  
  
}