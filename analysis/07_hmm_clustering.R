#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 07_hmm_clustering.R
# Purpose: Clustering from HMM output to define transmission units
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

## ---- initialize
library(tidyverse)
library(spData) 
library(sf)
library(here)
library(plotly)
library(cowplot)
library(ggraph)
library(igraph)
library(cluster)
library(factoextra)
library(NbClust)
library(glue)

options(scipen = 999999999)


# Source functions
# which_obs <- "all_data" # update if running sensitivity analyses or testing
source(here("analysis","00_functions_settings.R"), local = T)

if(simulate){
  out_dir <- str_c(here(), "/generated_data/stan/simulated_data")
  if (!dir.exists(out_dir)) {dir.create(out_dir)}
}

# run_model <- F
# test_subset <- T
# min_year <- 1990
# set.seed(seed)
# filter_cutoff <- 0.5

#-----------------------------------------------------------#
# Import----
#-----------------------------------------------------------#
##----import

source(here("analysis","02_import_data.R"), local = T)

#check analysis
sampling

#Bring in connectivity measure (used as weight/metric for clustering) from model output
model_runs <- file.mtime(paste0(xi_dir,"/draws/",list.files(paste0(xi_dir,"/draws/"))))
xi_model_draws <- qs2::qs_read(here::here(xi_dir,"draws",list.files(paste0(xi_dir,"/draws/")))[which(model_runs==max(model_runs))])
xi_model_draws <- xi_model_draws |> 
  select(-c(".chain",".iteration",".draw"))

# number of draws for multiple clustering
nsamp <- nrow(xi_model_draws) 

#-----------------------------------------------------------#
# Spatial objects----
#-----------------------------------------------------------#

# Country list
countries <- all_countries
if(test_subset){
  countries <- countries_subset
}
if(simulate){
  countries <- sim_countries
}

# Indices
M <- length(countries) 
n_time <- 1

# Spatial objects
spatial_objs <- create_connections(M, n_time)

# Distance matrix
distances <- st_centroid(africa_sf  |> 
                           filter(country %in% countries)) |> 
  st_distance() |> 
  as.data.frame() |> 
  set_names((africa_sf  |> 
               filter(country %in% countries))$country) |> 
  purrr::map_dfc( ~ {units::set_units(.x, km)}) 
distances <- distances |> 
  mutate(c = names(distances)) |> 
  mutate(m = purrr::map_int(c, ~which(countries == .x))) |> 
  arrange(m) |> 
  dplyr::select(all_of(countries)) |> 
  as.matrix()

### Edges & Nodes----
##----create_edges

nodes <- make_nodes(countries)

edges <- make_edges(nodes)

if(simulate){
  edges <- edges |> 
    left_join(p_xi |> dplyr::select(node_n, true_xi) |> 
                mutate(node_n = as.integer(node_n)), 
              by = c("n"="node_n")) |> 
    filter(node1!=node2) 
}

#Create symmetric edge list
single_edges <- make_single_edges(edges)

edge_list <- as.matrix(single_edges |> select(country1, country2))

#create graph object for louvain clustering
graph <- igraph::graph_from_edgelist(as.matrix(edge_list), directed = F)

#-----------------------------------------------------------#
# Louvain clustering----
#-----------------------------------------------------------#

## Connectivity distribution---- 
# sample from draws
xi_draws_samp <- xi_model_draws %>% 
  select(-starts_with(".")) %>% 
  mutate(draw = row_number()) %>% 
  pivot_longer(cols = -draw, names_to = "n", names_transform = list(n=as.integer),
               names_pattern = "xi\\[(\\d*+)\\]",
               values_to = "xi") %>% 
  arrange(draw, n) |> 
  left_join(edges |> 
              select(n, country1, country2) |> 
              list_countries(), by = "n") 

xi_draws_samp_w <- xi_draws_samp |> 
  pivot_wider(id_cols = n, names_from = draw, 
              names_prefix = "xi_draw_", values_from = xi) 

#combine with edges
xi_df_list <- purrr::map(1:(length(xi_draws_samp_w)-1), ~{
  varname <- glue::glue("xi_draw_{.x}")
  varname_q <- ensym(varname)
  d <- single_edges |> 
    left_join(xi_draws_samp_w |> 
                select(n, all_of(varname)) |> 
                rename(xi1=!!varname_q), 
              by = c("n1"="n")) |>
    left_join(xi_draws_samp_w |> 
                select(n, all_of(varname)) |> 
                rename(xi2=!!varname_q), 
              by = c("n2"="n")) |> 
    mutate(max_xi = pmax(xi1, xi2)) |> 
    calc_cluster_metrics(measure = max_xi, metric_name = "xi") |> 
    select(-xi1, -xi2) |> 
    arrange(country2, country1)
  return(d)
})

louvain_weights_list <- purrr::map(xi_df_list, ~{
  wt <- .x |> 
    select(max_xi) |> 
    pull()
  return(wt)
})

## Run Louvain clustering----
louvain_list <- purrr::map(1:nsamp, ~{cluster_louvain(graph, weights = louvain_weights_list[[.x]])})

##----louv_hist
## Count optimal cluters----
# the number of optimal clusters from louvain method
cluster_count <- purrr::map_int(louvain_list, ~{
  max(unique(.x$membership))
})
cluster_count <- data.frame(cluster_count) |> 
  mutate(cluster_count = as.factor(cluster_count))
clusters_hist <- ggplot(cluster_count) +
  geom_bar(aes(x = cluster_count), stat = "count") +
  theme_classic()
clusters_hist

## Divide into clusters----
louvain_hmm_subgrps <- purrr::imap(louvain_list, ~{
  k_clust = max(.x$membership)
  members <- data.frame(draw = .y, 
                        cluster = .x$membership,
                        country = .x$names)
  df <- tibble(draw = .y, 
               country1 = rep(members$country, times = length(members$country)),
               country2 = rep(members$country, each = length(members$country))) |> 
    left_join(members |> rename(cluster1 = cluster), by = c("country1"="country")) |> 
    left_join(members |> rename(cluster2 = cluster), by = c("country2"="country")) |> 
    mutate(same = if_else(cluster1==cluster2, 1, 0))
  return(df)}) |> 
  bind_rows() |> 
  group_by(country1, country2) |> 
  summarize(same = sum(same)) |> 
  ungroup() |> 
  mutate(pct = same/nsamp) |> 
  filter(country1!=country2)

louvain_hmm_subgrps2 <- louvain_hmm_subgrps |> 
  select(country1, country2, pct) |> 
  mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
  mutate(country1 = factor(country1, levels = heatmap_order, 
                           ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), 
                           ordered = T)) |> 
  arrange(country2, country1)

## Louvain heatmap----
louvain_heatmap <- ggplot(data = louvain_hmm_subgrps2, 
                          aes(y = country2, x = country1)) +
  geom_tile(aes(fill = pct), color = "white",
            lwd = 0.25,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#FFFFE0", mid = "#79CDCD", high = "#000066", 
                       midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")
louvain_heatmap

## Louvain as distance----
louvain_hmm_subgrps <- louvain_hmm_subgrps |> 
  mutate(inv_pct = 1/pct)
louvain_diff_mat <- diff_metric_matrix_func(louvain_hmm_subgrps, metric = pct, country_list = countries)
louvain_diff_dist <- dist(louvain_diff_mat)

## Consensus clustering----
unique_weights <- length(unique(louvain_hmm_subgrps$pct))
subgrps <- louvain_hmm_subgrps |> 
  list_countries() |> 
  select(countries, pct) |> 
  distinct()

# Start consensus algorithm
iter <- 0
consensus_threshold <- 0.15 
while(unique_weights>2){
  iter <- iter+1
  new_D <- subgrps |> 
    mutate(new_weight = if_else(pct>consensus_threshold, pct, 0)) |> 
    select(countries, new_weight)
  new_weights_list <- purrr::map(xi_df_list, ~{
    wt <- .x |> 
      list_countries() |> 
      left_join(new_D, by = "countries") |> 
      select(country1, country2, new_weight) |> 
      distinct() |> 
      select(new_weight) |> 
      pull()
    return(wt)
  })
  new_louvain_list <- purrr::map(1:nsamp, ~{cluster_louvain(graph, weights = new_weights_list[[.x]])})
  subgrps <- purrr::imap(new_louvain_list, ~{
    k_clust = max(.x$membership)
    members <- data.frame(draw = .y, 
                          cluster = .x$membership,
                          country = .x$names)
    df <- tibble(draw = .y, 
                 country1 = rep(members$country, times = length(members$country)),
                 country2 = rep(members$country, each = length(members$country))) |> 
      left_join(members |> rename(cluster1 = cluster), by = c("country1"="country")) |> 
      left_join(members |> rename(cluster2 = cluster), by = c("country2"="country")) |> 
      mutate(same = if_else(cluster1==cluster2, 1, 0))
    return(df)}) |> 
    bind_rows() |> 
    group_by(country1, country2) |> 
    summarize(same = sum(same)) |> 
    ungroup() |> 
    mutate(pct = same/nsamp) |> 
    filter(country1!=country2) |> 
    list_countries() |> 
    select(countries, pct)
  unique_weights <- length(unique(subgrps$pct))
  print(glue("iteration: {iter}, unique weights: {unique_weights}"))
}

new_louvain_subgrps2 <- louvain_hmm_subgrps |> 
  select(country1, country2) |> 
  list_countries() |> 
  left_join(subgrps |> 
              select(countries, pct) |> distinct(), by = "countries") |> 
  mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
  mutate(country1 = factor(country1, levels = heatmap_order, 
                           ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), 
                           ordered = T)) |> 
  arrange(country2, country1)

# Visualize consensus clustering
new_louvain_heatmap <- ggplot(data = new_louvain_subgrps2, 
                              aes(y = country2, x = country1)) +
  geom_tile(aes(fill = pct), color = "white",
            lwd = 0.25,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#FFFFE0", mid = "#79CDCD", high = "#000066", 
                       midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")
new_louvain_heatmap
# 3 clusters: c("Nigeria", "Djibouti", "South Africa")
# 4 clusters: c("Liberia", "Nigeria", "Djibouti", "South Africa")

consensus_clusters <- purrr::imap(c("Nigeria", "Djibouti", "South Africa"), 
                                  ~{
                                    d <- new_louvain_subgrps2 |> 
                                      filter(country1==.x | 
                                               (country2==.x  & pct==1)) |> 
                                      select(country1) |> 
                                      distinct() |> 
                                      mutate(cluster = .y)
                                    return(d)
                                  }) |> 
  bind_rows()



## Validation: Hierarchical clustering on louvain output----
# Run hierarchical clustering on consensus clustering output
sil_plot <- fviz_nbclust(louvain_diff_mat, FUN = hcut, method = "silhouette")
sil_plot
gapstat_plot <- fviz_nbclust(louvain_diff_mat, FUN = hcut, method = "gap_stat")
gapstat_plot
#3-5 clusters

louvain_k_clust <- 3

# get agglomerative coefficient for each linkage method
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x, d) {
  agnes(d, method = x)$ac
}
purrr::map_dbl(m, ac, d = louvain_diff_dist)
## ALL DATA
# average    single  complete      ward 
# 0.6379071 0.4396745 0.6877312 0.9046040

## Post-1990
# average    single  complete      ward 
# 0.6063893 0.3850411 0.6426837 0.8891913

louvain_diana_clust <- diana(louvain_diff_dist)
# get divisive coefficient for DIANA
louvain_diana_clust$dc
#ALL DATA divisive coefficient: 0.6879017
#post-1990 divisive coefficient: 0.6368392

# Run DIANA clustering
louvain_hclust <- diana(louvain_diff_dist)

dendogram <- fviz_dend(louvain_hclust, cex = 0.8, type = "phylogenic",
                       k = 3, repel = T, phylo_layout = "layout_as_tree",
                       palette = "lancet") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
dendogram

louvain_hclusters <- tibble(country = names(cutree(as.hclust(louvain_hclust), k = 2)),
                            cluster_2k = factor(cutree(as.hclust(louvain_hclust), k = 2)),
                            cluster_3k = factor(cutree(as.hclust(louvain_hclust), k = 3)),
                            cluster_4k = factor(cutree(as.hclust(louvain_hclust), k = 4)),
                            cluster_5k = factor(cutree(as.hclust(louvain_hclust), k = 5)))


# Shapefiles for visualizing 
louvain_clust_sf <- africa_sf |> 
  left_join(louvain_hclusters,
            by = "country")
consensus_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(consensus_clusters,
            by = c("country"="country1")) |> 
  mutate(cluster = as.factor(cluster))

cluster_palette <- c("1" = "#a1d99b", "5" = "#e66101", "4" = "darkorange", 
                     "2" = "darkorange2", "3" = "#31a354")

#map mean clusters
pal <- c("1" = "#a1d99b", "3" = "#31a354", "2" = "mediumpurple1", "4" = "mediumpurple3")

louvain_clust_map_plot <- ggplot() + 
  geom_sf(data = full_africa_sf, fill = "transparent", color = "black") +
  geom_sf(data = louvain_clust_sf, 
          aes(fill = cluster_3k), color = "grey5", alpha = 0.85) +
  geom_sf(data = louvain_clust_sf |> group_by(cluster_3k) |> summarize(),
          fill = "transparent", color = "grey45", linewidth = 0.75) + 
  theme_void() +
  scale_fill_manual(values = pal, na.value = "grey70") +
  theme(legend.position = "none") 
louvain_clust_map_plot



pal2 <- c("3" = "#a1d99b", "2" = "#31a354", "4" = "mediumpurple3", "1" = "mediumpurple1")

consensus_clust_map_plot <- ggplot() + 
  geom_sf(data = full_africa_sf, fill = "transparent", color = "black") +
  geom_sf(data = consensus_clust_sf, 
          aes(fill = cluster), color = "grey5", alpha = 0.85) +
  geom_sf(data = consensus_clust_sf |> group_by(cluster) |> summarize(),
          fill = "transparent", color = "grey45", linewidth = 0.75) + 
  theme_void() +
  scale_fill_manual(values = pal2, na.value = "grey70") +
  theme(legend.position = "none") 
consensus_clust_map_plot


#-----------------------------------------------------------#
# Mean Connectivity (for heat map)----
#-----------------------------------------------------------#

# Mean connectivity measure
xi_summ <- xi_model_draws |>
  summarise_all(mean) |>
  pivot_longer(cols = names(xi_model_draws), names_to = "variable", values_to = "mean") |>
  mutate(node_n = as.integer(str_extract(variable, "(?<=\\[)[0-9]+")))

# Calculate metrics from average xi from model runs
mean_xi_df <- single_edges |>
  left_join(xi_summ |> select(node_n, xi1 = mean), by = c("n1"="node_n")) |>
  left_join(xi_summ |> select(node_n, xi2 = mean), by = c("n2"="node_n")) |>
  mutate(max_xi = pmax(xi1, xi2)) |>
  calc_cluster_metrics(measure = max_xi, metric_name = "xi") |>
  select(-xi1, -xi2)

# Object for all edges
mean_xi_df_alledges <- mean_xi_df |>
  list_countries() |>
  select(countries, contains("xi")) |>
  right_join(edges |> select(country1, country2, countries, n), by = "countries") |>
  select(country1, country2, countries, n, everything()) |>
  arrange(n)

## Distance measure for plotting
xi_dist_plotting <- mean_xi_df_alledges |>
  mutate_at(c("country1", "country2"), ~country_rename(.))  |>
  mutate(country1 = factor(country1, levels = heatmap_order, ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), ordered = T)) |>
  select(country1, country2, log_xi_norm) |>
  arrange(country2, country1)

# Visualize
ggplot(data = xi_dist_plotting,
       aes(y = country2, x = country1)) +
  geom_tile(aes(fill = log_xi_norm, alpha = log_xi_norm), color = "white",
            lwd = 0.3,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#5CB3FF", mid = "grey90", high = "darkorange2",
                       midpoint = 0, na.value = "grey96", limits = c(-2.5, 4)) +
  scale_alpha_continuous(range = c(0.7,0.98)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

# Edges to map
hmm_mapping_edges <- edges |>
  left_join(mean_xi_df_alledges |> select(n, log_xi_norm, xi_norm), by = "n") |>
  group_by(countries) |>
  filter(row_number()==1) |>
  ungroup() |>
  arrange(log_xi_norm) |>
  mutate(country1 = factor(country_rename(country1), levels = heatmap_order, ordered = T),
         country2 = factor(country_rename(country2), levels = heatmap_order, ordered = T))

#-----------------------------------------------------------#
# Validation: HCluster Connectivity Distribution ----
#-----------------------------------------------------------#

xi_df_alledges_list <- purrr::map(xi_df_list, ~{
  d <- .x |> 
    list_countries() |> 
    select(countries, contains("xi")) |> 
    right_join(edges |> select(country1, country2, n) |> 
                 list_countries(), by = "countries") |> 
    select(country1, country2, n, everything()) |> 
    arrange(n) |> 
    mutate_at(c("country2", "country1"), ~case_when(
      .=="lesotho" ~ "AALesotho",
      .=="ghana" ~ "ABGhana",
      .=="kenya" ~ "ACKenya",
      T ~ .
    )) |> 
    arrange(country2, country1)
  return(d)
})



# convert to matrices of dissimilarity measure for DIANA
countries_orig <- countries
countries[countries=="lesotho"] <- "AALesotho"
countries[countries=="ghana"] <- "ABGhana"
countries[countries=="kenya"] <- "ACKenya"
countries <- sort(countries)
inv_xi_mat_list <- purrr::map(xi_df_alledges_list, ~{
  mat <- xi_metric_matrix_func(.x, inv_xi)
  return(mat)
})
inv_xi_dists <- purrr::map(inv_xi_mat_list, as.dist)

## Clustering: Distribution----

### HClustering, 3 clusters----
# run divisive clustering 
clusters_samps <- purrr::map(inv_xi_dists, hclust, method = "complete")

k_clust <- 3
##----clustgroup
# divide into clusters
hclust_sub_grp_samps <- purrr::map(clusters_samps, ~cutree(as.hclust(.x), k = k_clust))

#counts of subgroups for visualization
hclust_subgrps_df <- purrr::imap(hclust_sub_grp_samps, ~{
  #identify those that cluster with Nigeria vs not
  # clusters <- c(if_else(.x[["Nigeria"]]!=1,1,2),
  #                 .x[["Nigeria"]], 
  #                 if_else(.x[["Nigeria"]]!=3,3,2))
  # clusters <- c(1,
  #               if_else(.x[["Nigeria"]]!=2,2,.x[["Nigeria"]]), 
  #               3)
  clusters <- c(1:3)
  
  # assign clusters based on scheme above
  df <- tibble(draw = .y,
               country = names(.x),
               cluster = as.numeric(.x)) |> 
    left_join(purrr::map(.x, ~which(clusters==.x)) |> bind_cols() |> 
                pivot_longer(cols = everything(),
                             names_to = "country",
                             values_to = "cluster_assign"),
              by = "country")
  return(df)
}, k = k_clust) %>% 
  bind_rows() %>% 
  mutate_at(c("country", "cluster"), ~factor(.)) %>%
  # determine frequency of clustering in each grou
  count(cluster_assign, country, .drop = F) %>% 
  rename(count = n) %>% 
  mutate(cluster_assign = as.integer(cluster_assign),
         country = as.character(country)) %>% 
  arrange(country, cluster_assign) |> 
  mutate(country = if_else(str_detect(country,"AA|AB|AC"), 
                           str_remove_all(country, "AA|AB|AC"), country))

### HClustering, 2 clusters----
k_clust2 <- 2
##----clustgroup
# divide into clusters
hclust_sub_grp2_samps <- purrr::map(clusters_samps, ~cutree(as.hclust(.x), k = k_clust2))

#counts of subgroups for visualization
hclust_subgrps2_df <- purrr::imap(hclust_sub_grp2_samps, ~{
  #identify those that cluster with Nigeria vs not
  #clusters <- c(if_else(.x[["Nigeria"]]==2,1,2), .x[["Nigeria"]])
  clusters <- c(1,2)
  # assign clusters based on scheme above
  df <- tibble(draw = .y,
               country = names(.x),
               cluster = as.numeric(.x)) |> 
    left_join(purrr::map(.x, ~which(clusters==.x)) |> bind_cols() |> 
                pivot_longer(cols = everything(),
                             names_to = "country",
                             values_to = "cluster_assign"),
              by = "country")
  return(df)
}) %>% 
  bind_rows() %>% 
  mutate_at(c("country", "cluster_assign"), ~factor(.)) %>%
  # determine frequency of clustering in each grou
  count(cluster_assign, country, .drop = F) %>% 
  rename(count = n) %>% 
  mutate(cluster_assign = as.integer(cluster_assign),
         country = as.character(country)) %>% 
  arrange(country, cluster_assign)  |> 
  mutate(country = if_else(str_detect(country,"AA|AB|AC"), 
                           str_remove_all(country, "AA|AB|AC"), country))

hclust_subgrps_list <- list(hclust_subgrps_df, hclust_subgrps2_df) |> 
  set_names(c("hclust_3k", "hclust_2k"))

### Map Hclusters----

#### Cluster iterations map----
##----clustmap_samp
#map clusters across multiple iterations
iter_clust_sf <- africa_sf |> 
  left_join(hclust_subgrps_df,
            by = "country") |> 
  mutate(cluster_assign = factor(cluster_assign, levels = 1:k_clust))
countries <- countries_orig

clust_map_samps_plotlist <- purrr::map(1:3, ~{
  ggplot() +
    geom_sf(data = full_africa_sf, fill = "transparent", color = "grey40") +
    geom_sf(data = iter_clust_sf |> filter(cluster_assign==.x), 
            aes(fill = cluster_assign, alpha = count)) +
    # geom_sf(data = mean_clust_sf |> group_by(cluster) |> summarize(),
    #         fill = "transparent", color = "grey45", linewidth = 0.75) + 
    # geom_sf(data = mean_clust2_sf |> group_by(cluster) |> summarize(),
    #         fill = "transparent", color = "grey20", linewidth = 0.9) + 
    theme_void() +
    scale_fill_manual(values = pal) + 
    scale_alpha_continuous(range = c(0,1)) +
    theme(legend.position = "none") 
})
clust_map_samps_plotlist
clust_map_samps_plot <- cowplot::plot_grid(clust_map_samps_plotlist[[2]], 
                                           clust_map_samps_plotlist[[3]],
                                           clust_map_samps_plotlist[[1]], nrow = 3)
clust_map_samps_plot


#### Country-pair heat map----
#determine how often each country clusters with other countries 
##----clust_heat
hclust_subgrps_bycountry <- purrr::imap(hclust_sub_grp_samps, ~{
  clust <- tibble(country = names(.x),
                  cluster = .x)
  df <- tibble(country1 = rep(names(.x), times = length(names(.x))),
               country2 = rep(names(.x), each = length(names(.x)))) |> 
    left_join(clust |> rename(cluster1 = cluster), by = c("country1"="country")) |> 
    left_join(clust |> rename(cluster2 = cluster), by = c("country2"="country")) |> 
    mutate(same = if_else(cluster1==cluster2, 1, 0))
  return(df)
}) %>% 
  bind_rows() %>% 
  mutate_at(c("country1", "country2"), ~factor(.)) %>%
  group_by(country1, country2) |> 
  summarize(same = sum(same)) |> 
  ungroup()
hclust_subgrps_bycountry <- hclust_subgrps_bycountry |> 
  mutate(pct = same/nsamp) |> 
  mutate_at(c("country1", "country2"), ~if_else(str_detect(.,"AA|AB|AC|zz"), 
                                                str_remove_all(., "AA|AB|AC|zz"), .)) |> 
  arrange(country1, country2)

hclust_subgrps_bycountry2 <- purrr::imap(hclust_sub_grp2_samps, ~{
  clust <- tibble(country = names(.x),
                  cluster = .x)
  df <- tibble(country1 = rep(names(.x), times = length(names(.x))),
               country2 = rep(names(.x), each = length(names(.x)))) |> 
    left_join(clust |> rename(cluster1 = cluster), by = c("country1"="country")) |> 
    left_join(clust |> rename(cluster2 = cluster), by = c("country2"="country")) |> 
    mutate(same = if_else(cluster1==cluster2, 1, 0))
  return(df)
}) %>% 
  bind_rows() %>% 
  mutate_at(c("country1", "country2"), ~factor(.)) %>%
  group_by(country1, country2) |> 
  summarize(same = sum(same)) |> 
  ungroup()
hclust_subgrps_bycountry2 <- hclust_subgrps_bycountry2 |> 
  mutate(pct_2k = same/nsamp) |> 
  mutate_at(c("country1", "country2"), ~if_else(str_detect(.,"AA|AB|AC|zz"), 
                                                str_remove_all(., "AA|AB|AC|zz"), .)) |> 
  arrange(country1, country2)

hclust_subgrps_bycountry <- hclust_subgrps_bycountry |> 
  left_join(hclust_subgrps_bycountry2 |> select(-same), by = c("country1", "country2"))  |> 
  mutate_at(c("country1", "country2"), ~country_rename(.))

#visualize heatmap of how often individual countries cluster together
hclust_countries_heatmap <- ggplot(data = hclust_subgrps_bycountry |> 
                                    mutate_at(c("country1", "country2"), ~if_else(str_detect(.,"AA|AB|AC"), 
                                                             str_remove_all(., "AA|AB|AC"), .)) |> 
                                    mutate_at(c("country1", "country2"), ~if_else(.=="Democratic Republic Of The Congo", "DRC", .))  |> 
                                    filter(country1!=country2) |> 
                                    mutate(country1 = factor(country1, 
                                                             levels = heatmap_order, 
                                                             ordered = T),
                                           country2 = factor(country2, 
                                                             levels = rev(heatmap_order), 
                                                             ordered = T)), 
                                  aes(x = country1, y = country2)) +
  geom_tile(aes(fill = pct)) +
  theme_classic() +
  scale_fill_gradient2(low = "#ffffd9", mid = "#66cccc", high = "#000066", midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
hclust_countries_heatmap


#### Country clustering maps----
##----clust_country
countries_plotting <- country_rename(countries) 
#map how often countries cluster together, faceted by country
map_list <- purrr::map(countries_plotting, ~{
  ggplot() + 
    geom_sf(data = full_africa_sf |> mutate(country = country_rename(country)), 
            fill = "transparent", color = "grey40") + 
    geom_sf(data = africa_sf |> 
              mutate(country = country_rename(country)) %>% 
              left_join(hclust_subgrps_bycountry |> filter(country1==.x),
                        by = c("country"="country2")), aes(fill = pct), alpha = 0.9) +
    geom_curve(aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                   linewidth = log_xi_norm, 
                   alpha = log_xi_norm,
                   color = log_xi_norm), 
               data = hmm_mapping_edges |> 
                 filter(country1==.x | country2==.x),
               angle = 145,
               arrow = arrow(length = unit(0.05,"cm"))) +
    map_options_func(0.57, "Xi", exp = F, color_l = "lightgoldenrod2", color_h = "darkred", color_m = "white") +
    theme_void() +
    ggtitle(.x) +
    scale_linewidth_continuous(range = c(0.02, 0.7)) + 
                               #limits = c(log(min(max_edges$mean_xi)), log(max(max_edges$mean_xi)))) +
    scale_fill_gradient2(low = "white", mid = "lightskyblue2", high = "dodgerblue4", midpoint = 0.5) +
    theme(legend.position = "none",
          plot.title.position = 'plot', 
          plot.title = element_text(size = 8.5, hjust = 0.5)) }) |> 
  set_names(countries_plotting)
map_list_grid <- cowplot::plot_grid(plotlist = map_list, ncol = 5)
map_list_grid



##Save out----
saveRDS(xi_dist_plotting, str_c(clust_dir, "/xi_dist_plotting.rds"))
saveRDS(louvain_hmm_subgrps2, str_c(clust_dir, "/louvain_hmm_subgrps.rds"))
saveRDS(hmm_mapping_edges, str_c(clust_dir, "/hmm_mapping_edges.rds"))
saveRDS(louvain_hclusters, str_c(clust_dir, "/hmm_louvain_hclusters.rds"))
saveRDS(hclust_subgrps_list, str_c(clust_dir, "/hmm_hclust_subgrps_list.rds"))
saveRDS(cluster_count, str_c(clust_dir, "/hmm_cluster_count.rds"))
saveRDS(consensus_clusters, str_c(clust_dir, "/hmm_consensus_clusters.rds"))

if(save_plots){
  
  ggsave(plot = dist_plot, 
         filename = glue::glue("{fig_dir}/dist_plot.png"))
  ggsave(plot = clust_map_samps_plot, 
         filename = glue::glue("{fig_dir}/cluster_map_{sampling}_{k_clust}clust_{nsamp}iter.png"))
  ggsave(plot = vizclust,
         filename = glue::glue("{fig_dir}/vizclust_{k_clust}clust.png"))
  ggsave(map_list_grid,
         filename = glue::glue("{fig_dir}/cluster_map_bycountry_{sampling}_{k_clust}clust_{nsamp}iter.png"),
         height = 12, width = 17, bg = "white")
  ggsave(plot = clust_map_plot,
         filename = glue::glue("{fig_dir}/cluster_map_{sampling}_{k_clust}clust.png"))
  
}
