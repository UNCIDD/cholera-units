#-----------------------------------#
#Run results reports
#-----------------------------------#

library(here)
library(tidyverse)

get_min_year <- function(obs = which_obs){
  if(obs=="all_data"){
    yr <- 1970
  }else{
    yr <- as.numeric(str_extract(obs, "(?<=\\_)[:digit:]*"))
  }
  return(yr)
}

render_report <- function(report, ...){
  if(report=="results"){
    md <- "01_hmm_results.Rmd"
  }else if(report=="clustering"){
    md <- "04_clustering_results.Rmd"
  }else if(report=="diagnostics"){
    md <- "02_hmm_diagnostics.Rmd"
  }else if(report=="test"){
      md <- "test.Rmd"
  }
  
  if(!is.null(title_info)){
    title <- str_c("OUCh Model ", str_to_title(str_replace_all(report, "_", " ")), ": ", 
                   str_to_title(str_replace_all(which_obs, "_", " ")),
                   " ", title_info)
    outname <- str_c(glue::glue("HMM_{report}_{which_obs}_", str_replace_all(title_info, " |,", "_"), ".html"))
  }else{
    title <- str_c("OUCh Model ", 
                   str_to_title(str_replace_all(report, "_", " ")), ": ", 
                   str_to_title(str_replace_all(which_obs, "_", " ")))
    outname <- glue::glue("HMM_{report}_{which_obs}.html")
    
  }
  
  params <- list(doc_title = title,
                 downsample_random = downsample_random, 
                 which_obs = which_obs,
                 stan_chains = stan_chains,
                 stan_iter_warmup = stan_iter_warmup,
                 stan_iter_sample = stan_iter_sample,
                 save_plots = save_plots,
                 min_year = get_min_year(which_obs),
                 seed = seed,
                 downsample_seed = downsample_seed,
                 sample_fraction = sample_fraction,
                 sep_sudan = sep_sudan)
  
  if(report=="results"){
    params <- append(params,
                     list(run_model = run_model,
                          save_warmup = save_warmup))
  }
  if(report=="diagnostics"){
    params <- append(params,
                     list(save_warmup = save_warmup))
  }
  
  
  rmarkdown::render(here("notebooks", md),
                    output_file = outname,
                    output_dir = here("notebooks/html", Sys.Date()),
                    params = params, 
                    envir = new.env())
  
}
