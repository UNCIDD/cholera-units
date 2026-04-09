#--------------------------------------------#
#Run results reports - Downsampling, 20%
#--------------------------------------------#

library(here)

source(here("notebooks/generate_reports", "render_report.R"), local = T)
source(here("notebooks/generate_reports", "default_params.R"), local = T)


downsample_run_num <- 5 ##UPDATE FOR DIFFERENT SAMPLES 

## Update default knitting params

downsample <- T
sample_fraction <- 0.8
downsample_seed <- 6843 #seed+4
run_model <- T
stan_chains <- 4
stan_iter_warmup <- 1000
stan_iter_sample <- 2000
save_warmup <- F

title_info <- glue::glue("downsample{sample_fraction*100} {downsample_run_num}")

#Render----

render_report(report = "results")
#render_report(report = "clustering")
render_report(report = "diagnostics")

