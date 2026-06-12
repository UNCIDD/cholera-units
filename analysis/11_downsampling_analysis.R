#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 11_downsampling_analysis.R
# Purpose: Validation across downsampled partitions, 
#          each with 20% downsampling
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(here)

options(scipen = 9999)

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)
set.seed(seed)

#-----------------------------------------------------------#
# Read data----
#-----------------------------------------------------------#

source(here("analysis","02_import_data.R"), local = T)
downsample_groups <- readRDS(here("data","downsample_groups.RDS"))

downsample_random <- T
sample_fraction <- 0.8

results_list <- purrr::map(c(1:5), ~{
  downsample_run_num <- .x 
  source(here("analysis","00_functions_settings.R"), local = T)
  source(here("analysis","02_import_data.R"), local = T)
  print(downsample_run_num)
  print(xi_dir)
  cases_lineages_hmm <- readRDS(str_c(xi_dir, "/cases_lineages_assign.rds")) |> 
    mutate(country = as.factor(country)) |> 
    clean_df()
  cases_prev_pred <- readRDS(str_c(xi_dir, "/cases_prev_pred.rds"))
  y_obs <- readRDS(str_c(xi_dir, "/downsampled_y.rds"))
  pred_rho <- readRDS(str_c(xi_dir, "/pred_rho.rds")) 
  pred_z <- readRDS(str_c(xi_dir, "/pred_z_star.rds")) 
  results_list <- list(cases_lineages_hmm, cases_prev_pred, y_obs, pred_rho, pred_z) |> 
    set_names(c("cases_lineages_hmm", "cases_prev_pred", "y_obs", "pred_rho", "pred_z"))
  return(results_list)
}) |> 
  set_names(str_c("iter_", c(1:5)))

#-----------------------------------------------------------#
# Spatial objects----
#-----------------------------------------------------------#

#create spatial mappings
countries <- all_countries
if(simulate){
  countries <- sim_countries
}
M <- length(countries) 
n_time <- (max(cases_lineages$year)-min(cases_lineages$year))+1    # number of time slices
N <- M * n_time  # Max number of space*time observations
spatial_objs <- create_connections(M, n_time)

#-----------------------------------------------------------#
# Combine originals for comparison----
#-----------------------------------------------------------#

orig <- downsample_groups |> 
  filter(.groups %in% str_remove_all(names(results_list), "iter_")) |> 
  ungroup()

#-----------------------------------------------------------#
# Identify true/false positive/negative----
#-----------------------------------------------------------#

downsampled <- purrr::map2(results_list, names(results_list), ~{
  samps <- orig |> 
    filter(.groups %in% str_remove_all(.y, "iter_")) |> 
    mutate(present = 1)
  z_samps <- .x$pred_z |> 
    select(location, time, median, strain) |> 
    left_join(samps |> distinct(location, time, present), by = c("location", "time")) |> 
    filter(present==1)
  z <- samps |> 
    mutate(strain = as.character(as.numeric(te))) |> 
    full_join(z_samps |> select(strain, location, time, median),
              by = c("location", "time", "strain")) |> 
    arrange(location, time, strain) |> 
    translate_cols() |> 
    filter(present==1 | median==1) |> 
    select(country, year, te, z = median, present) |> 
    mutate(tp = if_else(!is.na(present) & present==1 & z==1, 1, 0),
           fp_clonal = if_else(is.na(present) & z==1, 1, 0),
           fn = if_else(!is.na(present) & present==1 & z==0, 1, 0)) |> 
    group_by(country, year) |> 
    mutate(fp_nonclonal = if_else(is.na(present) & z==1 & max(tp==1, na.rm = T), 0, fp_clonal)) |> 
    group_by(te) |>
    mutate(te_tp = sum(tp),
           te_fp = sum(fp_clonal),
           te_fn = sum(fn)) |>
    ungroup() |>
    mutate(tp = sum(tp),
           fp_clonal = sum(fp_clonal),
           fp_nonclonal = sum(fp_nonclonal),
           fn = sum(fn)) |>
    select(te, starts_with("te_"), tp, fp_clonal, fp_nonclonal, fn) |>
    distinct() |>
    mutate(recall_tot = tp/(tp+fn),
           recall_te = te_tp/(te_tp+te_fn),
           precision_tot = (tp)/(tp+fp_clonal),
           precision_tot_nonclonal = (tp)/(tp+fp_nonclonal),
           precision_te = (te_tp)/(te_tp+te_fp)) |>
    arrange(te)
  
  return(z)
}) |> 
  bind_rows(.id = "group")


#-----------------------------------------------------------#
# Calculate recall/precision----
#-----------------------------------------------------------#

# Precision & recall overall, assuming both clonal and nonclonal outbreaks
validation_overall <- downsampled |> 
  select(group, recall_tot, precision_tot, precision_tot_nonclonal) |> 
  distinct() |> 
  summarize(mean_recall = mean(recall_tot, na.rm = T),
            sd_recall = sd(recall_tot, na.rm = T),
            mean_precision = mean(precision_tot, na.rm = T),
            sd_precision = sd(precision_tot, na.rm = T),
            mean_precision_nonclonal = mean(precision_tot_nonclonal, na.rm = T),
            sd_precision_nonclonal = sd(precision_tot_nonclonal, na.rm = T))

# Precision & recall by lineage
validation_te <- downsampled |> 
  select(group, te, recall_te, precision_te) |> 
  distinct() |> 
  group_by(te) |> 
  summarize(mean_recall_te = mean(recall_te, na.rm = T),
            sd_recall_te = sd(recall_te, na.rm = T),
            mean_precision_te = mean(precision_te, na.rm = T))

# Overall
validation_overall
# By lineage
validation_te

# True positives
downsampled |> 
  distinct(group, tp) |> 
  summarize(tp_tot = sum(tp))

# False negatives
downsampled |> 
  distinct(group, fn) |> 
  summarize(fn_tot = sum(fn))

# False positives, clonal
downsampled |> 
  distinct(group, fp_clonal) |> 
  summarize(fp_cloncal_tot = sum(fp_clonal))

# False positives, nonclonal
downsampled |> 
  distinct(group, fp_nonclonal) |> 
  summarize(fp_nonclonal_tot = sum(fp_nonclonal))




