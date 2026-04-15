#--------------------------------------------#
#Run results reports - All data, gravity
#--------------------------------------------#

library(here)

source(here("notebooks/generate_reports", "render_report.R"), local = T)
source(here("notebooks/generate_reports", "default_params.R"), local = T)

# sampling_options <- c("all_data", "post_1990", "post_2000")

## Update default knitting params

run_model <- T
# test_subset <- T
# stan_iter_sample <- 2000
# stan_iter_warmup <- 2000
# min_year <- 1990
# save_warmup <- F

source(here("analysis","05_run_hmm.R"), local = T, echo = TRUE, print.eval = TRUE)

#Render----

run_model <- F

render_report(report = "results")
#render_report(report = "clustering")
render_report(report = "diagnostics")

