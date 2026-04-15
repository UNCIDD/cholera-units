
#----------Functions for data cleaning & analysis-----------------#
#------------Also includes globals & plot themes------------------#


### Data cleaning----

#' function to clean and standardize country names & lineages
#' @param data dataset to modify
#' @param seq whether the dataset contains sequencing information
#' @return returns a tibble

clean_df <- function(data, seq = F) {
  d <- data |> 
    janitor::clean_names() |> 
    #Make country names consistent across files, and remove special characters
    #Fix apostrophes, remove accents, remove strange characters
    mutate(country = str_trim(str_to_lower(country)),
           country = str_replace_all(country, c("’"="'", 
                                                "ô"="o",
                                                "ã"="a",
                                                "é"="e",
                                                "í"="i",
                                                "democratic"="dem\\.",
                                                "republic"="rep\\.",
                                                " \\(.*"="",
                                                "sao tome .*"="sao tome",
                                                ".* south africa"="south africa",
                                                ".* tanzania"="tanzania")),
           country = country |> 
             replace_values(
               "swaziland" ~ "eswatini",
               "congo" ~ "rep. of the congo",
               "dem. rep. of the congo" ~ "drc"
               )) |> 
    filter_out(is.na(country) | country=="?") 
    
  if(seq) {
    if("collection_year" %in% names(d)){
      d <- d |> 
        rename(year = collection_year)
    }
    d <- d |> 
      mutate(across(c(accession, taxa), ~str_replace_all(., "[^0-9A-Za-z///' //-//_]", ""))) |> 
      mutate(te = case_when(
        str_detect(te, "new lineage 1") ~ "T16",
        str_detect(te, "new lineage 2") ~ "T17",
        str_detect(str_to_lower(te), "sporadic") ~ "sporadic",
        str_detect(te, "AFR") ~ str_replace(te,"AFR","T"),
        te=="UND" ~ NA_character_,
        T ~ te)) |> 
      filter_out(is.na(te)) |> 
      select(taxa, accession, country, year, te) |> 
      arrange(country, year, te, accession)
  }
  return(d)
}

#' Function to update country names for consistency
#' @param var variable coding country names

country_rename <- function(var){
    case_when(
      {{ var }} == "drc" ~ "DRC", 
      str_detect({{ var }},"cote") ~ "Cote d'Ivoire",
      T ~ str_to_title({{ var }}))
}

### Stan prep----

# Printing functions for result reports
print_sample_func <- function(){
  if(test_subset){
    print(glue::glue("     Fully observed data in a subsample of {M} countries and {n_time} years."))
  }else{
    if(downsample_random){
      print(glue::glue("     randomly downsample all locations and strains with a sampling fraction of {sample_fraction}."))
    }else{
      print(glue::glue("     Fully observed data."))
    }
  }
}

# print_locations_func <- function(subsample, loc){
#   if(subsample){
#     if(loc=="country"){
#       print(str_c(countries_subset, collapse = ", "))
#     }else if(loc=="region"){
#       print(str_c(unique(subregions$subregion[which(subregions$country %in% countries_subset)]), collapse = ", "))
#     }
#     
#   }else{
#     print("All countries")
#   }
# }

#' Function to create connectivity measures
#' 
#' 
create_connections <- function(n_countries, time_slices) {
  #Create edge matrices
  adj_nodes <- tibble(node1 = rep(1:n_countries, each = n_countries), 
                      node2 = rep(1:n_countries, times = n_countries)) |> 
    filter(node1 != node2) |> 
    arrange(node2) |> 
    mutate(n = row_number()) 
  
  adj_nodes_single <- adj_nodes %>% 
    filter(node1 > node2) 
  
  node1_by_dest <- adj_nodes$node1
  node2_by_dest <- adj_nodes$node2
  
  start_by_dest <- adj_nodes |> group_by(node2) |> slice_min(n) |> pull(n)
  end_by_dest <-  adj_nodes |> group_by(node2) |> slice_max(n) |> pull(n)
  
  map_from_edge <- matrix(0, nrow = n_countries, ncol = n_countries)
  for (i in 1:n_countries) {
    for(j in 1:n_countries) {
      ind <- which(node1_by_dest == i & node2_by_dest == j)
      if (length(ind) > 0) {
        map_from_edge[i, j] = ind
      }
    }
  }
  
  
  #Year to country mapping vector
  map_to_location <- rep(1:n_countries, each = time_slices)
  #Country to year mapping vector
  map_to_time <- rep(1:time_slices, times = n_countries)
  
  map_to_n <- bind_cols(n = 1:(time_slices*n_countries), m = map_to_location, t = map_to_time) |>  
    pivot_wider(id_cols = t, names_from = m, values_from = n) |> 
    dplyr::select(-t) |> 
    as.matrix()
  colnames(map_to_n) <- NULL
  
  connectivity <- list(adj_nodes, node1_by_dest, node2_by_dest, start_by_dest, 
                       end_by_dest, map_from_edge, map_to_location, map_to_time,
                       map_to_n) %>% 
    set_names(c("adj_nodes", "node1_by_dest", 
                "node2_by_dest", "start_by_dest", 
                "end_by_dest", "map_from_edge", 
                "map_to_location", "map_to_time",
                "map_to_n"))
  return(connectivity)
}

# Functions to create node and edge lists

make_nodes <- function(country_set = countries){
  st_centroid(africa_sf  |> 
                filter(country %in% country_set)) |> 
    dplyr::select(country, geometry) |> 
    mutate(long = purrr::map_dbl(geometry, ~.x[1]),
           lat = purrr::map_dbl(geometry, ~.x[2])) |> 
    st_drop_geometry()
}

make_edges <- function(nodes_list, country_set = countries){
  spatial_objs$adj_nodes |> 
    mutate(country1 = purrr::map_chr(node1, ~ {country_set[[.x]]}),
           country2 = purrr::map_chr(node2, ~ {country_set[[.x]]})) |> 
    left_join(nodes_list |> rename(xstart = long, ystart = lat), by = c("country1" = "country")) |> 
    left_join(nodes_list |> rename(xend = long, yend = lat), by = c("country2" = "country")) |> 
    list_countries()
}

make_single_edges <- function(edge_list) {
  edge_list |> 
    select(country1, country2, countries, n) |> 
    group_by(countries) |> 
    mutate(n1 = min(n), n2 = max(n)) |> 
    filter(row_number()==1) |> 
    ungroup() |> 
    select(-n, -countries)
}


### Analysis----

#### Post-stan Functions & settings ----

pd <- position_dodge(width = 1e-6)

#trace plots
traceplot <- function(obj = model_objs$model_draws, pars){
  bayesplot::mcmc_trace(obj, regex_pars = pars)  +
    ggplot2::scale_color_discrete() +
    theme(legend.position = "bottom",
          text = element_text(family = "Times New Roman"))
}

#parameter summaries
model_fit_summary <- function(var, true_val=NULL, name = "mean", ...){
  if(is.null(true_val)){
    true_val <- NA
  }
  df <- model_objs$var_summary |> 
    filter(variable==var) |> 
    mutate(true = true_val) |>  
    dplyr::select(true, mean, q5, q95, rhat, ess_bulk, ess_tail) 
  if(name!="mean"){
    nameq <- sym(str_c("mean_",str_remove_all(var,"\\[.*")))
    df <- df |> 
      rename(!!nameq := mean)
    return(df)
  }
  return(df)
}

translate_cols <- function(df, min = min_year, sim = simulate, locs = countries_subset){
  if(sim){
    countries_subset <- sim_countries
  }
  if("strain" %in% names(df)){
    df <- df |> 
      mutate(country = purrr::map_chr(location, ~ {locs[[as.numeric(str_remove(.x, "loc-"))]]}),
             te = factor(purrr::map_chr(strain, ~{te_subset[[as.numeric(.x)]]}), ordered = T, levels = te_order),
             year = {{min}}+time-1)
  }else{
    df <- df |> 
      mutate(country = purrr::map_chr(location, ~ {locs[[as.numeric(str_remove(.x, "loc-"))]]}),
             year = {{min}}+time-1)
  }
  return(df)
}

#output dataframes
poststan_df <- function(var, state = T){
  if(state){
    d <- model_objs$var_summary |> 
      filter(str_detect(variable, var)) |> 
      mutate(strain = str_extract(variable, "(?<=\\[)[0-9]+(?=,)"),
             state = as.numeric(str_extract(variable, "(?<=,)[0-9]+(?=,)"))-1,
             ind = str_extract(variable, "(?<=,)[0-9]+(?=\\])") |> as.numeric(),
             location = str_c("loc-", spatial_objs$map_to_location[ind]),
             time = spatial_objs$map_to_time[ind]) |> 
      filter(state==1) |> 
      dplyr::select(-state)
  }else{
    d <- model_objs$var_summary |> 
      filter(str_detect(variable, var)) |> 
      mutate(strain = str_extract(variable, "(?<=,)[0-9]+(?=\\])"),
             ind = str_extract(variable, "(?<=\\[)[0-9]+(?=,)") |> as.numeric(),
             location = str_c("loc-", spatial_objs$map_to_location[ind]),
             time = spatial_objs$map_to_time[ind])
  }
  
}

#plot dataframes
postest_plot <- function(est_df, obs = df_obs, title, show_cases = T){

  d_obs <- obs |>  
    translate_cols() |> 
    mutate(pa = as.double(samples > 0)) |> 
    dplyr::select(country, year, te, samples, pa) |>  
    filter(pa == 1) |> 
    distinct() |> 
    arrange(country, year, te) |> 
    group_by(country, year) |> 
    mutate(pa = row_number()/10) |> 
    ungroup()
  
  plot <- ggplot(translate_cols(est_df)) +
    geom_ribbon(aes(x = year, ymin = q5, ymax = q95, fill = te), 
                show.legend = T, alpha = .1, position = position_dodge(0.2)) +
    geom_line(aes(x = year, y = mean, color = te), 
              position = position_dodge(0.2), show.legend = T) +
    facet_wrap(~country, ncol = 2) +
    theme_bw() +
    ggside::geom_xsidepoint(data = d_obs, aes(x = year, y = pa, color = te), 
                            alpha = 0.8, size = 0.9, show.legend = T) +
    ggside::scale_xsidey_continuous(limits = c(0, max(d_obs$pa + 0.1)), minor_breaks = NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(title = title,
         color = "TE",
         fill = "TE") +
    theme_bw() +
    plottheme +
    plot_colors +
    plot_fillcolors +
    guides(fill=guide_legend(nrow=2, title.position = "left"),
           color=guide_legend(nrow=2, title.position = "left")) +
    ggside::theme_ggside_void() +
    theme(legend.position = "bottom",
          ggside.axis.text.y = element_blank(),
          ggside.panel.scale = 0.22)
  return(plot)
}


#### Clustering----

list_countries <- function(data){
  df <- data |> 
    mutate(countries = unlist(purrr::map(purrr::map(str_split(str_c(country1, ", ",
                                                                    country2), ", "),
                                                    sort),
                                         str_c, collapse = ", ")))
  return(df)
}

norm_func <- function(measure, type = c("mean", "robust", "standard")){
  if(type=="mean"){
    x <- ({{measure}}-mean({{measure}}, na.rm = T))/(max({{measure}}, na.rm = T)-min({{measure}}, na.rm = T))
  }else if(type=="robust"){
    x <- ({{measure}}-median({{measure}}, na.rm = T))/(quantile({{measure}}, p = 0.75, na.rm = T)-quantile({{measure}}, p = 0.25, na.rm = T))
  }else{
    x <- scale({{measure}})
  }
}

calc_cluster_metrics <- function(data, measure, metric_name, log_measure = NULL) {
  if(is.null(log_measure)){
    d <- data |> 
      mutate(log_met = log({{measure}}))
  }else{
    d <- data |> 
      rename(log_met = log_measure)
  }
  d <- d |> 
    mutate(inv_met = if_else({{measure}}==0, 1/({{measure}}+1e-7), 1/{{measure}}),
           inv_log_met = 1/log_met,
           #normalize 
           met_norm = norm_func({{measure}}, "standard"),
           met_norm2 = norm_func({{measure}}, "mean"),
           met_norm_robust = norm_func({{measure}}, "robust"),
           log_met_norm = norm_func(log_met, "standard"),
           log_met_norm2 = norm_func(log_met, "mean"),
           log_met_norm_robust = norm_func(log_met, "robust"),
           #inverse metric as dissimilarity measure
           inv_met_norm = norm_func(inv_met, "standard"),
           inv_log_met_norm = norm_func(inv_log_met, "standard"))
  names(d) <- str_replace_all(names(d), "_met", glue("_{metric_name}"))
  names(d) <- str_replace_all(names(d), "met_", glue("{metric_name}_"))
  return(d)
}

xi_metric_matrix_func <- function(data, metric){
  mat <- data |> 
    select(country1, country2, {{metric}}) |> 
    arrange(country1, country2) |> 
    pivot_wider(id_cols = c(country1), names_from = country2, values_from = {{metric}}) |>  
    arrange(country1) |> 
    select(all_of(countries)) |>
    as.matrix()
  mat[which(is.na(mat))] <- 0
  row.names(mat) <- colnames(mat)
  return(mat)
}

diff_metric_matrix_func <- function(data, metric, country_list = countries){
  mat <- data |> 
    select(country1, country2, {{metric}}) |> 
    arrange(country1, country2) |> 
    pivot_wider(id_cols = c(country1), names_from = country2, values_from = {{metric}}) |>  
    arrange(country1) |> 
    select(all_of({{country_list}})) |>
    #select(all_of(countries[countries!="South Sudan"])) |> 
    as.matrix()
  mat[which(is.na(mat))] <- 0
  row.names(mat) <- colnames(mat)
  return(mat)
}

##### Heatmap order----
# heatmap_order <- c("Guinea", "Sierra Leone", "Liberia", "Mauritania", "Guinea-Bissau",
#                    "Senegal", "Cote D'ivoire", "Mali", "Togo",
#                    "Burkina Faso", "Benin", "Ghana", "Nigeria", "Cameroon", "Niger", "Chad", "Equatorial Guinea", "Sao Tome", "Gabon",
#                    "Republic Of The Congo", "Central African Republic", "South Sudan", 
#                    "Sudan", "Angola", "DRC", "Ethiopia", "Somalia",
#                    "Uganda", "Tanzania", "Rwanda", "Kenya", "Zambia", "Madagascar", "Djibouti",
#                    "Namibia", "Burundi", "Comoros", 
#                    "Botswana", "Zimbabwe", "Malawi", "Mozambique", "Eswatini", "Lesotho", "South Africa")

# Used for Figs 3, S1, S3
heatmap_order <-  c("Liberia", "Cote d'Ivoire", "Guinea", "Guinea-Bissau", 
                    "Mali", "Mauritania", "Senegal", "Sierra Leone", "Benin", 
                    "Burkina Faso", "Togo", "Nigeria", "Cameroon", "Ghana",  
                    "Chad", "Equatorial Guinea", "Niger", "Rep. of the Congo", 
                    "Central African Rep.", "Gabon", "Sao Tome", 
                    "Somalia",  "Ethiopia", "Eritrea", "Djibouti", 
                    "Rwanda", "Sudan", "South Sudan", "Uganda", "DRC", 
                    "Burundi", "Kenya", "Tanzania", "Zambia", 
                    "Angola", "Comoros", "Madagascar", "Malawi", 
                    "Mozambique", "Zimbabwe", "South Africa", "Eswatini", 
                    "Lesotho", "Botswana", "Namibia")  

# Combining Sudan & South Sudan (Used for Figs 4-5)
heatmap_order_comb <-  c("Liberia", "Cote d'Ivoire", "Guinea", "Guinea-Bissau", 
                    "Mali", "Mauritania", "Senegal", "Sierra Leone", "Benin", 
                    "Burkina Faso", "Togo", "Nigeria", "Cameroon", "Ghana",  
                    "Chad", "Equatorial Guinea", "Niger", "Rep. of the Congo", 
                    "Central African Rep.", "Gabon", "Sao Tome", 
                    "Somalia",  "Ethiopia", "Eritrea", "Djibouti", 
                    "Rwanda", "Sudan & South Sudan", "Uganda", "DRC", 
                    "Burundi", "Kenya", "Tanzania", "Zambia", 
                    "Angola", "Comoros", "Madagascar", "Malawi", 
                    "Mozambique", "Zimbabwe", "South Africa", "Eswatini", 
                    "Lesotho", "Botswana", "Namibia")  

### Outbreak Simulation functions----

pull_params <- function(var){
  connectivity$post_warmup_draws  |> 
    select(starts_with(var))
}

sim_spread_func <- function(simnum, loc_int, yrs = n_time, n_M = M, 
                            par = param_draws, mu = case_dists, downstream = T,
                            nruns = runs){
  
  iter <- sample(1:nruns, 1)
  
  #sample parameters
  eta <- par$eta[iter,] |> pull()
  delta <- par$delta[iter,] |> pull()
  xi <- t(par$xi)[,iter]
  
  #create spatial mappings
  if(!downstream){
    yrs <- 2
  }
  sp <- create_connections(n_M, yrs)
  
  #initialize matrices
  alpha <- matrix(0, nrow = n_M, ncol = yrs)
  phi <- alpha
  p <- alpha
  cases <- alpha
  
  #initialize introductions
  alpha[loc_int,1] <- 1
  p[loc_int,1] <- 1
  
  #sample cases
  lambda_cases <- mu |> filter(country_id==loc_int) |> slice_sample(n=1) |> pull(cases)
  cases[loc_int,1] <- rpois(1,lambda_cases)
  
  
  if(downstream){
    for (t in 2:(yrs)) {
      for(i in 1:n_M){
        
        #initialize
        phi[i,t] = 1
        
        for (k in sp$start_by_dest[i]:sp$end_by_dest[i]) {
          j = sp$node1_by_dest[k] 
          ji = sp$map_from_edge[j,i]
          
          #update with introduction rate, cases in j, connectivity, & presence in j
          phi[i,t] = phi[i,t] * (1-(1-exp(-(cases[j,t-1]^eta * xi[ji]))))
          
        }#end other country loop
        
        phi[i,t] = phi[i,t] * (1-(1-exp(-(cases[i,t-1]^eta * delta))))
        
        #probability of colonization -> probability of presence
        alpha[i,t] = 1 - phi[i,t] 
        
        #Present/absent
        p[i,t] <- rbinom(1,1,alpha[i,t])
        
        #cases & samples
        lambda_cases <- mu |> filter(country_id==i) |> slice_sample(n=1) |> pull(cases) 
        cases[i,t] <- p[i,t]*rpois(1,lambda_cases)
        
      }#end country loop
    }#end time loop
  }else{
    i <- loc_int
    t <- 2
    #initialize
    phi[i,t] = 1
    #from i to j instead here
    for (k in sp$start_by_dest[i]:sp$end_by_dest[i]) {
      j = sp$node1_by_dest[k] 
      ij = sp$map_from_edge[i,j]
      
      #update with introduction rate, cases in j, connectivity, & presence in j
      phi[j,t] = 1-(1-exp(-(cases[i,t-1]^eta * xi[ij])))
      #probability of colonization -> probability of presence
      alpha[j,t] = 1 - phi[j,t] 
      
      #Present/absent
      p[j,t] <- rbinom(1,1,alpha[j,t])
      
      #cases & samples
      lambda_cases <- mu |> filter(country_id==j) |> slice_sample(n=1) |> pull(cases) 
      cases[j,t] <- p[j,t]*rpois(1,lambda_cases)
      
      
    }#end other country loop
    #probability of colonization -> probability of presence
    alpha[i,t] <- 1 
    
    #Present/absent
    p[i,t] <- 1
    
    #cases & samples
    lambda_cases <- mu |> filter(country_id==i) |> slice_sample(n=1) |> pull(cases) 
    cases[i,t] <- p[i,t]*rpois(1,lambda_cases)
  }
  
  out <- outbreak_df_func(p, "Presence", downstream = downstream)
  return(out)
}

arrival_func <- function(df){
  df |> 
    mutate(arrived = if_else(Presence==1, 1, as.numeric(NA))) |> 
    group_by(country) |> 
    fill(arrived) |> 
    mutate(arrived = if_else(is.na(arrived),0,arrived),
           yr = time-1) |> 
    filter(arrived==1 | max(arrived)==0) |> 
    mutate(arrival_time = if_else(arrived==1,min(yr, na.rm = T),n_time)) |> 
    select(country, arrival_time) |> 
    distinct() |> 
    ungroup()
}

mean_arrival <- function(df, max_time = n_time) {
  df |> 
    group_by(country) |> 
    summarize(arrival_time = mean(arrival_time)) |> 
    mutate(arrival_time = case_when(
      arrival_time==0 ~ 1,
      arrival_time>=max_time ~ max_time,
      T ~ arrival_time)) |> 
    ungroup()
}

secondary_func <- function(df){
  df |> 
    filter(time==2) |> 
    select(-time) |> 
    group_by(country) |> 
    summarize(pct = mean(Presence)) |> 
    ungroup()
}

outbreak_df_func <- function(parameter, label, downstream = T){
  pred <- as.data.frame({{parameter}}) 
  pred <- pred %>% 
    magrittr::set_colnames(str_c("time-", 1:ncol(pred))) %>% 
    mutate(location = str_c("loc-", 1:M)) %>% 
    pivot_longer(contains("time"),
                 names_to = "time",
                 values_to = label) %>% 
    mutate(time = as.numeric(str_remove_all(time, "time-"))) |> 
    translate_func(min = 1) 
  arrival <- arrival_func(pred)
  if(downstream){
    return(arrival)
  }else{
    return(pred)
  }
}

#--------------------------------------------#
### Global settings----
# Set globals for runs/analyses/saving
#--------------------------------------------#

#### initialize lineages---- 
te_order <- c(str_c("T", c(1:13, 15:17)), "S")
pub_te_order <- c(str_c("AFR", c(1:13, 15:17)), "Sporadic Outbreak")

#### excluded countries---- 
north_africa_exclude <- c("algeria", "egypt", "libya", "morocco", "tunisia", "western sahara")
islands_exclude <- c("british indian ocean territory", "french southern and antarctic lands", 
                     "saint helena, ascension and tristan da cunha", "mauritius", "seychelles",
                     "cape verde", "cabo verde")

#parameters

#### Stan model settings---- 
get_min_year <- function(obs = which_obs){
  if(obs=="all_data"){
    yr <- 1970
  }else{
    yr <- as.numeric(str_extract(obs, "(?<=\\_)[:digit:]*"))
  }
  return(yr)
}

if(!exists("run_model")){
  run_model <- F
}
if(!exists("save_warmup")){
  save_warmup <- T
}
if(!exists("stan_iter_warmup")){
  stan_iter_warmup <- 250
}
if(!exists("stan_iter_sample")){
  stan_iter_sample <- 200
}
if(!exists("stan_chains")){
  stan_chains <- 4
}
if(!exists("seed")){
  seed <- 4
}
#whether to combine Sudan and South Sudan
if(!exists("sep_sudan")){
  sep_sudan <- F
}

#sampling schemes
if(!exists("sampling")){
  # sampling_options <- c("all_data", "post_1990", "post_2000")
  sampling <- "all_data"
}
if(!exists("min_year")){
  min_year <- get_min_year(sampling)
}
if(!exists("downsample_random")){
  downsample_random <- F 
}
if(!exists("sample_fraction")){
  sample_fraction <- 0.8 
}
if(!exists("downsample_seed")){
  downsample_seed <- 5 
}
if(downsample_random){
  sampling <- glue::glue("downsample{sample_fraction*100}_{downsample_run_num}")
}
if(!exists("test_subset")){
  #whether to subset
  test_subset <- F
}
#simulation to test HMM 
if(!exists("simulate", mode="character")){
  simulate <- F
}

#whether to save out
if(!exists("save_plots")){
  save_plots <- F
}

#directories
out_dir <- str_c(here(), "/generated_data/stan/", sampling)
if(sep_sudan){
  out_dir <- str_c(out_dir,"_sepSudan")
}

fig_dir <- str_c(here(), "/figures/analysis_output/", sampling)
if(sep_sudan){
  fig_dir <- str_c(fig_dir,"_sepSudan")
}

xi_dir <- str_c(here(), "/generated_data/xi_dfs/", sampling)

if(sep_sudan){
  xi_dir <- str_c(xi_dir,"_sepSudan")
}

sim_dir <- str_c(here(), "/generated_data/sim_dfs/", sampling)
if(sep_sudan){
  sim_dir <- str_c(sim_dir,"_sepSudan")
}

clust_dir <- str_c(here(), "/generated_data/clustering/", sampling)
if(sep_sudan){
  clust_dir <- str_c(clust_dir,"_sepSudan")
}

#--------------------------------------------#
### Plotting themes & functions----
#--------------------------------------------#

#Plotting settings
#### Lineage palette----
te_palette <- c("dodgerblue", "navy", "blue2", 
                "yellowgreen", "green4", 
                "lightcoral", "red2", 
                "gold1", "darkorange1", 
                "purple3", "mediumorchid2", "plum3",   
                "brown4","deeppink3",
                "lightskyblue2", "mediumaquamarine", 
                "grey20")
names(te_palette) <- te_order

pub_te_palette <- te_palette
names(pub_te_palette) <- pub_te_order

plot_colors <- list(scale_color_manual(values = te_palette, drop = F))
plot_fillcolors <- list(scale_fill_manual(values = te_palette, drop = F))

pub_plot_colors <- list(scale_color_manual(values = pub_te_palette, drop = F))
pub_plot_fillcolors <- list(scale_fill_manual(values = pub_te_palette, drop = F))

pub_cluster_palette <- c("1" = "#a1d99b", "2" = "#31a354", "3" = "mediumpurple3", "4" = "mediumpurple1")

#### Theme objects----
##### overall plottheme----
plottheme <-  theme(text = element_text(size = 9, family = "serif"),
                    legend.title = element_text(size = 9),
                    legend.text = element_text(size = 9),
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9))

##### case_lineage_plot_theme----
case_lineage_plot_theme <- list(theme(text = element_text(family = "serif"),
                                      # Axes
                                      axis.title.x = element_text(size = 5.5),
                                      axis.title.y = element_blank(),
                                      axis.text = element_text(size = 5.5),
                                      axis.ticks = element_line(linewidth = 0.01),
                                      axis.ticks.length = unit(0.08,"cm"),
                                      axis.line = element_line(linewidth = 0.08),
                                      # Legend
                                      legend.position = "bottom",
                                      legend.title = element_text(size = 6, face = "bold"), 
                                      legend.title.position = "top",
                                      legend.text = element_text(size = 5.5),
                                      legend.text.position = "right",
                                      legend.key.width = unit(0.1,"cm"),
                                      legend.key.height = unit(0.1,"cm"),
                                      legend.key.spacing = unit(1, "mm"),
                                      # Plot
                                      plot.margin = unit(c(t=0.1, r=0.0, b=0.05, l=0.05), "cm")))

map_options_func <- function(mid, param, exp = F, color_m = "white", 
                             color_l = "dodgerblue3", color_h = "firebrick4"){
  options <- list(
    scale_fill_continuous(low = "grey90", high = "grey50", na.value = "ivory1",
                          labels=scales::comma),
    theme_void(),
    theme(legend.position = "right",
          legend.key.width  = unit(0.25, "cm")),
    scale_linewidth_continuous(range = c(0.001, 1.6)),
    scale_alpha_continuous(range = c(0.3,1)),
    guides(fill = guide_colorbar(
      title = "Population"),
      color = guide_colorbar(
        title = param),
      linewidth = "none", 
      alpha = "none"))
  if(exp){
    options <- append(options, list(scale_color_gradient2(
      low = {{color_l}},
      mid = {{color_m}},
      high = {{color_h}},
      midpoint = {{mid}},
      trans = "exp",
      breaks = scales::pretty_breaks(n=4, min.n = 3),
      expand = 1)))
  }else{
    options <- append(options, list(scale_color_gradient2(
      low = {{color_l}},
      mid = {{color_m}},
      high = {{color_h}},
      midpoint = {{mid}})))
  }
  return(options)
}


#Function to plot prevalence by lineage as donut plot
prev_bar_plot <- function(data, y_var, te_var){
  data |> 
    ggplot() +
    geom_bar(aes(x = year, y = {{y_var}}, fill = {{te_var}}), 
             stat = "identity", position = "stack") +
    scale_fill_manual(values = c(pub_te_palette, "Unsequenced" = "grey50"), 
                      breaks = c(pub_te_order[-length(pub_te_order)], "Unsequenced"),
                      drop = F) +
    ylab("Reported Cholera Cases") +
    labs(color = "Lineage",
         fill = "Lineage") +
    guides(fill=guide_legend(nrow=1, title.position = "left", drop = F),
           color=guide_legend(nrow=1, title.position = "left", drop = F)) +
    theme_classic() +
    theme(text = element_text(family = "serif"),
          legend.position = "bottom",
          legend.text = element_text(size = 5, family = "serif"),
          legend.title = element_blank(), 
          legend.key.width = unit(0.1,"cm"),
          legend.key.height = unit(0.1,"cm"),
          plot.margin = unit(c(t=0.1, r=0.0, b=0.05, l=0.05), "cm"),
          axis.line = element_line(linewidth = 0.08),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.01),
          axis.ticks.length = unit(0.08,"cm"))
}

#Function to plot prevalence by lineage as donut plot
prev_donut_plot <- function(data, te_var = te_assign) {
  data |> 
    mutate(fract = te_count/sum(te_count),
           label = round(fract, 2),
           label = if_else(label>0.02, glue::glue("{as.character(label*100)}%"), "")) |> 
    arrange({{te_var}}) |> 
    mutate(ymax = cumsum(fract),
           ymin = ymax-fract,
           label_pos = (ymax+ymin)/2) |> 
    ggplot() +
    geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = {{te_var}})) +
    geom_label(x=3.6, aes(y=label_pos, label=label), size = 4.5, size.unit = "pt", 
               label.size=NA, fill = "transparent", 
               color = "white") +
    coord_polar(theta = "y") +
    xlim(c(2,4.2)) +
    scale_fill_manual(values = c(pub_te_palette, "Unsequenced" = "grey50"), 
                      breaks = c(te_order[-length(te_order)], "Unsequenced"),
                      drop = F) +
    theme_void() +
    theme(text = element_text(family = "serif"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")
}

#Function to plot cases and detected lineages by country and year
cases_lineages_plot <- function(case_data, lineage_data){
  if(!("source" %in% names(lineage_data))){
    lineage_data <- lineage_data |> 
      mutate(source = "Observed")
  }
  p <- ggplot() +
    geom_point(data = case_data,
               aes(x = year, y = country, size = cases), 
               color = "grey35", alpha = 0.1) +
    geom_point(data = case_data,
               aes(x = year, y = country, size = cases), 
               color = "grey40", shape = 1, alpha = 0.3) +
    geom_point(data = lineage_data |> 
                 filter(te!="S", source == "Observed") |> 
                 arrange(country, year, te),
               aes(x = year, y = country, color = te, shape = source),
               alpha = 0.65, size = 0.9,
               position = position_jitter(width = 0.18, height = 0.22, seed = seed)) +
    geom_point(data = lineage_data |> 
                 filter(te!="S", source == "Inferred") |> 
                 arrange(country, year, te),
               aes(x = year, y = country, color = te, shape = source), 
               alpha = 0.65, size = 0.8, 
               position = position_jitter(width = 0.18, height = 0.22, seed = seed)) +
    geom_point(data = lineage_data |> 
                 filter(te=="S") |> 
                 mutate(te = "Sporadic Outbreak"),
               aes(x = year, y = country, color = te, shape = source), 
               alpha = 0.9, size = 0.5, shape = 16) +
    theme_minimal() +
    xlab("Year") +
    plot_colors +
    scale_shape_manual(values = c("Inferred" = 8, "Observed" = 2)) +
    scale_size_continuous(range = c(0.01,3), breaks = c(100,1000,10000,100000)) +
    labs(size = "Reported Suspected\nCholera Cases", 
         color = "Lineage",
         shape = "Source") +
    case_lineage_plot_theme
  
  if(length(unique(lineage_data$source))<2){
    p <- p + 
      guides(size = guide_legend(order=2, nrow = 1),
             color = guide_legend(nrow = 1, override.aes=list(fill=NA)),
             shape = "none",
             fill = "none") 
  }else{
    p <- p +
      guides(size = guide_legend(order=2, nrow = 1),
             color = guide_legend(nrow = 1, override.aes=list(fill=NA)),
             shape = guide_legend(nrow = 1),
             fill = "none") 
  }
}

secondary_risk_map <- function(country_list, edgelist){
  map(country_order[which(country_order %in% 
                            {{country_list}}$seed_country)], 
      ~ {
        ggplot() + 
          geom_sf(data = full_africa_sf, fill = "grey90") +
          #map risk of outbreaks
          geom_sf(data = secondary_sfs[[which(names(secondary_sfs)==.x)]], 
                  aes(fill = pct), alpha = 1) +
          geom_sf(data = full_africa_sf, fill = "transparent", color = "grey30") +
          #plot connectivity measure
          geom_curve(data = edgelist |> 
                       filter(country1==.x),
                     aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                         linewidth = log_xi_norm,
                         alpha = log_xi_norm, 
                         color = log_xi_norm), 
                     curvature = -0.35,
                     arrow = arrow(length = unit(0.3, "mm"))) +
          theme_void() +
          ggtitle(.x) +
          scale_fill_gradient2(high = "#660000", mid = "#B21807", 
                               low = "white", midpoint = 0.5) +
          scale_color_gradient2(high = "yellow", mid = "#AAF0D1", 
                                low = "grey95", midpoint = 0, 
                                limits = c(min(edgelist$log_xi_norm),
                                           max(edgelist$log_xi_norm))) +
          scale_alpha_continuous(range = c(0.5,0.98)) +
          scale_linewidth_continuous(range = c(0.01, 0.5), 
                                     limits = c(min(edgelist$log_xi_norm),
                                                max(edgelist$log_xi_norm))) +
          guides(linewidth = "none", 
                 fill = guide_colorbar(title = "Risk of outbreak\nin following year", nrow = 1),
                 alpha = "none",
                 color = guide_colorbar(title = "Standardized Connectivity\nMeasure", nrow = 1)) +
          theme(text = element_text(family = "serif"),
                legend.position = "none",
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 7.5),
                legend.key.height = unit(1.5,"mm"),
                legend.spacing = unit(5, "mm"),
                plot.title.position = 'plot', 
                plot.title = element_text(size = 8, hjust = 0.5))
      }) 
}
