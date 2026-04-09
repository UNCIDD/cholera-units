#--------------------------------------------#
#Run results reports - All data, gravity
#--------------------------------------------#

library(here)

source(here("notebooks/generate_reports", "render_report.R"), local = T)
source(here("notebooks/generate_reports", "default_params.R"), local = T)

# sampling_options <- c("all_data", "post_1990", "post_2000")

## Update default knitting params

run_model <- T
# save_warmup <- F


#Render----

render_report(report = "results")
#render_report(report = "clustering")
render_report(report = "diagnostics")

