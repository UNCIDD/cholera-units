#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 08_phylogeography_clustering.R
# Purpose: Clustering analysis for phylogeographic analysis output
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

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
sampling <- "all_data"

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

set.seed(seed)

#-----------------------------------------------------------#
# Import----
#-----------------------------------------------------------#

##----import

sep_sudan <- T
source(here("analysis","02_import_data.R"), local = T)

phylo <- readRDS(here("data", "phylogeography.RDS"))

#-----------------------------------------------------------#
# Prep phylo data----
#-----------------------------------------------------------#

min_supported_rate <- min((phylo |> filter(location_indicators==1))$location_rates)

phylo2 <- phylo |> 
  filter_out(str_detect(countries, "Other")) |> 
  select(-countries) |> 
  standardize_countries(country_name = country1) |> 
  standardize_countries(country_name = country2) |> 
  list_countries() |> 
  select(country1, country2, countries, everything())
  
#-----------------------------------------------------------#
# Spatial objects----
#-----------------------------------------------------------#

# Country list
phylo_countries <- sort(unique(c(unique(phylo2$country1), unique(phylo2$country2))))

# Indices
M <- length(phylo_countries) 
n_time <- 1

# Spatial objects
spatial_objs <- create_connections(M, n_time)

#pull centroids for each country in subset & create distance matrix
distances <- st_centroid(africa_sf  |> 
                           filter(country %in% phylo_countries)) |> 
  st_distance() |> 
  as.data.frame() |> 
  set_names((africa_sf  |> 
               filter(country %in% phylo_countries))$country) |> 
  purrr::map_dfc( ~ {units::set_units(.x, km)}) 
distances <- distances |> 
  mutate(c = names(distances)) |> 
  mutate(m = purrr::map_int(c, ~which(phylo_countries == .x))) |> 
  arrange(m) |> 
  dplyr::select(all_of(phylo_countries)) |> 
  as.matrix()

# Create nodes and edges
phylo_nodes <- make_nodes(phylo_countries)
phylo_edges <- make_edges(phylo_nodes, country_set = phylo_countries)
  
#-----------------------------------------------------------#
# Phylo mapping objects----
#-----------------------------------------------------------#

#Median transition rate measure for map edges
phylo_mapping <- phylo2 |>
  rename(rate = location_rates) |>
  #for clustering, multiply rate by indicator for edge weights
  mutate(clust_rate = rate*location_indicators,
         log_clust_rate = log(if_else(location_indicators==1,rate*location_indicators,
                                      min_supported_rate/2))) |>
  group_by(countries) |>
  #calculate the percentage of runs where edge is supported
  mutate(ind_pct = sum(location_indicators)/n(),
         clust_rate = median(clust_rate),
         log_clust_rate = median(log_clust_rate)) |>
  ungroup() |>
  #for mapping, weight edges by median among runs with indicator=1
  mutate(rate = if_else(location_indicators==1, rate, as.numeric(NA))) |>
  summarize(rate = median(rate, na.rm = T), .by=c(country1, country2, countries, ind_pct,
                                                  clust_rate, log_clust_rate)) |>
  #for mapping, scale by % supported
  calc_cluster_metrics(measure = rate*ind_pct, metric_name = "rate") |>
  calc_cluster_metrics(measure = clust_rate, metric_name = "clust_rate", log_measure = "log_clust_rate") |>
  left_join(phylo_edges |> select(starts_with(c("n", "country"))),
            by = c("country1", "country2")) |>
  select(countries, starts_with("n"), everything()) |>
  arrange(n) |>
  distinct() |>
  mutate(across(contains("rate"), ~if_else(is.na(.), ind_pct, .)))

#object for all edges
phylo_mapping_alledges <- phylo_mapping |>
  select(country1, country2, countries, ind_pct, contains("rate")) |>
  bind_rows(phylo_mapping |>
              select(country1 = country2, country2 = country1, countries,
                     ind_pct, contains("rate"))) |>
  left_join(phylo_edges |>
              select(country1, country2, n),
            by = c("country1", "country2")) |>
  arrange(n)

phylo_dist_plotting <- phylo_mapping_alledges |>
  mutate_at(c("country1", "country2"), ~country_rename(.)) |>
  mutate(country1 = factor(country1, levels = heatmap_order, ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), ordered = T)) |>
  select(countries, contains("country"), contains("rate")) |>
  arrange(country2, country1)

# Edges to map
phylo_mapping_edges <- phylo_edges |>
  left_join(phylo_mapping_alledges |>
              select(n, ind_pct, clust_rate, log_clust_rate, rate_norm, rate_norm2,
                     log_rate_norm, log_rate_norm2,
                     log_rate, clust_rate_norm, clust_rate_norm2,
                     log_clust_rate_norm, log_clust_rate_norm2),
            by = "n") |>
  group_by(countries) |>
  filter(row_number()==1) |>
  ungroup() |>
  arrange(rate_norm)


#-----------------------------------------------------------#
# Louvain clustering----
#-----------------------------------------------------------#

# sample from all draws
phylo_clust_list <- purrr::map(unique(phylo2$state), ~{
  phylo2 |> 
    filter(state == .x) |> 
    rename(rate = location_rates) |> 
    #for clustering, multiply rate by indicator for edge weights
    mutate(clust_rate = rate*location_indicators,
           log_clust_rate = log(if_else(location_indicators==1,rate*location_indicators,
                                        min_supported_rate/2))) |> 
    group_by(countries) |>
    #calculate the percentage of runs where edge is supported
    mutate(ind_pct = sum(location_indicators)/n()) |>
    ungroup() |> 
    #for mapping, scale by % supported
    calc_cluster_metrics(measure = rate*ind_pct, metric_name = "rate") |> 
    calc_cluster_metrics(measure = clust_rate, metric_name = "clust_rate", log_measure = "log_clust_rate") |>
    left_join(phylo_edges |> select(starts_with(c("n", "country"))),
              by = c("country1", "country2")) |>
    select(countries, starts_with("n"), everything(), -contains("norm")) |>
    arrange(n) |>
    distinct()
})

# edge weights
louvain_phylo_weights_list <- purrr::map(phylo_clust_list, ~{
  wt <- .x |> 
    select(clust_rate) |> 
    pull()
  return(wt)
})


#All samples - Create symmetric edge list
single_edges_samps <- purrr::map(phylo_clust_list, ~{
  phylo_edges |>
    filter(n %in% .x$n) |> 
    make_single_edges()
})

edge_list_samps <- purrr::map(single_edges_samps, ~{
  as.matrix(.x |> select(country1, country2))
})

#create graph object for louvain clustering
graph_samps <- purrr::map(edge_list_samps, graph_from_edgelist, directed = F)

# Run Louvain clustering
louvain_list <- purrr::map(1:length(graph_samps), ~{
  cluster_louvain(graph_samps[[.x]], weights = louvain_phylo_weights_list[[.x]])})

#count the number of optimal clusters from louvain method
cluster_count <- purrr::map_int(louvain_list, ~{
  max(unique(.x$membership))
}) 
cluster_count <- data.frame(cluster_count) |> 
  mutate(cluster_count = as.factor(cluster_count))
clusters_hist <- ggplot(cluster_count) +
  geom_bar(aes(x = cluster_count), stat = "count") +
  theme_classic()
clusters_hist

# Louvain - divide into clusters
louvain_clusters <- purrr::imap(louvain_list, ~{
  members <- data.frame(draw = .y, 
                        cluster = .x$membership,
                        country = .x$names)
})

louvain_phylo_subgrps <- purrr::imap(louvain_list, ~{
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
  summarize(pct = sum(same)/n()) |> 
  ungroup() |> 
  filter(country1!=country2)

#### Visualize louvain clustering----

louvain_heatmap <- ggplot(data = louvain_phylo_subgrps |> 
                            mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
                            mutate(country1_s = factor(country1, 
                                                       levels = heatmap_order, 
                                                       ordered = T),
                                   country2_s = factor(country2, 
                                                       levels = rev(heatmap_order), 
                                                       ordered = T)), 
                            aes(y = country2_s, x = country1_s)) +
  geom_tile(aes(fill = pct), color = "white",
            lwd = 0.25,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#F4FAC4", mid = "#79CDCD", high = "#16256A", 
                       midpoint = 0.5, limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
louvain_heatmap

#### Consensus clustering----

unique_weights <- length(unique(louvain_phylo_subgrps$pct))

subgrps <- louvain_phylo_subgrps |> 
  list_countries() |> 
  select(countries, pct) |> 
  distinct()

# run consensus algorithm
iter <- 0
nsamp <- length(unique(phylo2$state))
consensus_threshold <- 0.15
while(unique_weights>2){
  iter <- iter+1
  new_D <- subgrps |> 
    mutate(new_weight = if_else(pct>consensus_threshold, pct, 0)) |> 
    select(countries, new_weight)
  new_weights_list <- purrr::map(phylo_clust_list, ~{
    wt <- .x |> 
      left_join(new_D, by = "countries") |> 
      distinct(country1, country2, new_weight) |> 
      select(new_weight) |> 
      pull()
    return(wt)
  })
  new_louvain_list <- purrr::map(1:nsamp, ~{cluster_louvain(graph_samps[[.x]], weights = new_weights_list[[.x]])})
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

# Assign groupings 
new_louvain_subgrps2 <- louvain_phylo_subgrps |> 
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

# Visualize
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

# Assign consensu clusters
consensus_clusters <- purrr::imap(c("Guinea", "Benin",
                                    "Sudan", "Malawi"), 
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

#### Hierarchical clustering----

# Using Louvain as distance metric, validate consensus clustering

summary(louvain_phylo_subgrps$pct)
diff_min <- quantile(louvain_phylo_subgrps$pct, 0.25)
louvain_phylo_subgrps <- louvain_phylo_subgrps |> 
  mutate(inv_pct = if_else(pct<diff_min, 1/diff_min, 1/pct))
louvain_diff_mat <- diff_metric_matrix_func(louvain_phylo_subgrps, metric = pct, country_list = phylo_countries)
louvain_diff_dist <- dist(louvain_diff_mat)

# Visualize "distances"
louvain_dist_heatmap <- ggplot(data = louvain_phylo_subgrps |> 
                            mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
                            mutate(country1_s = factor(country1, 
                                                       levels = heatmap_order, 
                                                       ordered = T),
                                   country2_s = factor(country2, 
                                                       levels = rev(heatmap_order), 
                                                       ordered = T)), 
                          aes(y = country2_s, x = country1_s)) +
  geom_tile(aes(fill = log(inv_pct)), color = "white",
            lwd = 0.275,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(high = "#F4FAC4", mid = "#79CDCD", low = "#16256A",
                       midpoint = 2.5,
                       limits = c(0,5.57)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
louvain_dist_heatmap

# optimal k
sil_plot <- fviz_nbclust(louvain_diff_mat, FUN = hcut, method = "silhouette")
sil_plot
gapstat_plot <- fviz_nbclust(louvain_diff_mat, FUN = hcut, method = "gap_stat")
gapstat_plot
#4 clusters


louvain_k_clust<-4

# Run DIANA
louvain_diana_clust <- diana(louvain_diff_dist)
louvain_diana_clust$dc
#divisive coefficient: 0.66

# Visualize clusters
dendogram <- fviz_dend(louvain_diana_clust, cex = 0.8, type = "phylogenic",
                       k = 4, repel = T, phylo_layout = "layout_as_tree",
                       palette = "lancet") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
dendogram

phylo_louvain_hclusters <- tibble(country = names(cutree(as.hclust(louvain_diana_clust), k = 2)),
                                  cluster_2k = factor(cutree(as.hclust(louvain_diana_clust), k = 2)),
                                  cluster_3k = factor(cutree(as.hclust(louvain_diana_clust), k = 3)),
                                  cluster_4k = factor(cutree(as.hclust(louvain_diana_clust), k = 4)),
                                  cluster_5k = factor(cutree(as.hclust(louvain_diana_clust), k = 5)))

# Shapefiles
louvain_clust_sf <- africa_sf |> 
  left_join(phylo_louvain_hclusters,
            by = "country")
consensus_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(consensus_clusters,
            by = c("country"="country1")) |> 
  mutate(cluster = as.factor(cluster))

pal <- c("4" = "#a1d99b", "1" = "#31a354", "3" = "mediumpurple1", "2" = "mediumpurple3")

# Map DIANA clusters
louvain_clust_map_plot <- ggplot() + 
  geom_sf(data = full_africa_sf, fill = "transparent", color = "black") +
  geom_sf(data = louvain_clust_sf, 
          aes(fill = cluster_4k), color = "grey5", alpha = 0.85) +
  geom_sf(data = louvain_clust_sf |> group_by(cluster_4k) |> summarize(),
          fill = "transparent", color = "grey45", linewidth = 0.75) + 
  scale_fill_manual(values = pal, na.value = "grey70") +
  theme_void() +
  theme(legend.position = "none") 
louvain_clust_map_plot

pal <- c("4" = "#a1d99b", "3" = "#31a354", "1" = "mediumpurple1", "2" = "mediumpurple3")

# Map consensus clusters
consensus_clust_map_plot <- ggplot() + 
  geom_sf(data = full_africa_sf, fill = "transparent", color = "black") +
  geom_sf(data = consensus_clust_sf, 
          aes(fill = cluster), color = "grey5", alpha = 0.85) +
  geom_sf(data = consensus_clust_sf |> group_by(cluster) |> summarize(),
          fill = "transparent", color = "grey45", linewidth = 0.75) + 
  scale_fill_manual(values = pal, na.value = "grey70") +
  theme_void() +
  theme(legend.position = "none") 
consensus_clust_map_plot


#-----------------------------------------------------------#
# Save out----
#-----------------------------------------------------------#

sep_sudan <- F
source(here("analysis","00_functions_settings.R"), local = T)

saveRDS(phylo_dist_plotting, str_c(clust_dir, "/phylo_dist_plotting.rds"))
saveRDS(louvain_phylo_subgrps, str_c(clust_dir, "/louvain_phylo_subgrps.rds"))
saveRDS(phylo_mapping_edges, str_c(clust_dir, "/phylo_mapping_edges.rds"))
saveRDS(phylo_louvain_hclusters, str_c(clust_dir, "/phylo_louvain_hclusters.rds"))
saveRDS(cluster_count, str_c(clust_dir, "/phylo_cluster_count.rds"))
saveRDS(consensus_clusters, str_c(clust_dir, "/phylo_consensus_clusters.rds"))


