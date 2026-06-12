#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 13_Figure3.R
# Purpose: This script creates Figure 3 showing clustering results
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
library(ggpattern)

options(scipen = 9999)

# Source functions
sampling <- "all_data"

source(here("analysis","00_functions_settings.R"), local = T)

set.seed(seed)

filter_cutoff <- 0.5

#-----------------------------------------------------------#
# Read data----
#-----------------------------------------------------------#

source(here("analysis","02_import_data.R"), local = T) 

phylo_dist_plotting <- readRDS(str_c(clust_dir, "/phylo_dist_plotting.rds"))
phylo_louvain_hclusters <- readRDS(str_c(clust_dir, "/phylo_louvain_hclusters.rds"))
phylo_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_phylo_subgrps.rds"))
phylo_mapping_edges <- readRDS(str_c(clust_dir, "/phylo_mapping_edges.rds"))
phylo_consensus_clusters <- readRDS(str_c(clust_dir, "/phylo_consensus_clusters.rds"))

hmm_dist_plotting <- readRDS(str_c(clust_dir, "/xi_dist_plotting.rds"))
hmm_louvain_hclusters <- readRDS(str_c(clust_dir, "/hmm_louvain_hclusters.rds"))
hmm_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_hmm_subgrps.rds"))
hmm_mapping_edges <- readRDS(str_c(clust_dir, "/hmm_mapping_edges.rds"))
hmm_consensus_clusters <- readRDS(str_c(clust_dir, "/hmm_consensus_clusters.rds"))

#-----------------------------------------------------------#
# Spatial objects----
#-----------------------------------------------------------#

countries <- all_countries
M <- length(countries) 
n_time <- 1
spatial_objs <- create_connections(M, n_time)

### Edges & Nodes
nodes <- make_nodes()
edges <- make_edges(nodes)

#-----------------------------------------------------------#
# Plotting prep----
#-----------------------------------------------------------#

## Filter edges----
hmm_filter_stat <- quantile(hmm_mapping_edges$log_xi_norm, filter_cutoff)
hmm_filter_stat_mean <- mean(hmm_mapping_edges$log_xi_norm)
phylo_filter_stat <- quantile(hmm_mapping_edges$log_xi_norm, filter_cutoff)

hmm_filtered_edges <- hmm_mapping_edges |> 
  filter(log_xi_norm>hmm_filter_stat_mean)
phylo_filtered_edges <- phylo_mapping_edges |> 
  filter(ind_pct>filter_cutoff)

## Distances----
hmm_dist_plotting <- hmm_dist_plotting |> 
  select(country1, country2, log_xi_norm) |> 
  mutate_at(c("country1", "country2"),  ~as.character(.))
phylo_dist_plotting <- phylo_dist_plotting |> 
  select(country1, country2, log_rate_norm) |> 
  mutate_at(c("country1", "country2"),  ~as.character(.))

## Consistent clustering----
clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Mauritania", "Central African Rep."),
  cluster_4k_assign = factor(c(1:4))) |> 
  left_join(hmm_louvain_hclusters |> select(country, hmm_cluster = cluster_4k),
            by = "country") |> 
  left_join(phylo_louvain_hclusters |> select(country, phylo_cluster = cluster_4k),
            by = "country")

hmm_louvain_hclusters <- hmm_louvain_hclusters |> 
  left_join(clust_reassign |> select(hmm_cluster, cluster_4k_assign),
            by = c("cluster_4k"="hmm_cluster"))
phylo_louvain_hclusters <- phylo_louvain_hclusters |> 
  left_join(clust_reassign |> select(phylo_cluster, cluster_4k_assign),
            by = c("cluster_4k"="phylo_cluster"))

consensus_clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Cameroon", "Mauritania"),
  cluster_phylo4k_assign = factor(c(1:4)),
  cluster_hmm3k_assign = factor(c(1:3, 3)))  |> 
  left_join(hmm_consensus_clusters |> select(country1, hmm_cluster = cluster),
            by = c("country"="country1")) |> 
  left_join(phylo_consensus_clusters |> select(country1, phylo_cluster = cluster),
            by = c("country"="country1"))
  
#Shapefiles for mapping
hmm_clust_sf <- africa_sf |> 
  left_join(hmm_louvain_hclusters,
            by = "country")
phylo_clust_sf <- africa_sf |> 
  left_join(phylo_louvain_hclusters,
            by = "country")


hmm_consensus_clusters <- hmm_consensus_clusters |> 
  left_join(consensus_clust_reassign |> select(contains("hmm")) |> distinct(),
            by = c("cluster" = "hmm_cluster")) |> 
  rename(cluster_assign = cluster_hmm3k_assign) |> 
  mutate(cluster_assign = factor(cluster_assign))
hmm_consensus_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(hmm_consensus_clusters,
            by = c("country"="country1")) 

consensus_clust_reassign |> select(contains("phylo")) |> distinct()
phylo_consensus_clusters <- phylo_consensus_clusters |> 
  left_join(consensus_clust_reassign |> select(contains("phylo")) |> distinct(),
            by = c("cluster" = "phylo_cluster")) |> 
  rename(cluster_assign = cluster_phylo4k_assign) |> 
  mutate(cluster_assign = factor(cluster_assign))
phylo_consensus_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(phylo_consensus_clusters,
            by = c("country"="country1"))

## Background layers for maps----
clust_map_background <- list(
  geom_sf_pattern(data = full_africa_sf |> filter(country %in% north_africa_exclude), 
                  fill = "grey90", color = "grey40", alpha = 0.2,
                  pattern_alpha = 0.2, pattern_spacing = 0.01, pattern_fill = "grey70",
                  pattern_colour = "grey30", pattern_size = 0.02),
  geom_sf(data = full_africa_sf |> filter(country %in% c(countries, "South Sudan")), 
          fill = "white", color = "grey40", alpha = 0.8)
)


## Mapping options----
min_color_val <- floor(min(min(hmm_mapping_edges$log_xi_norm, na.rm = T), 
                           min(phylo_mapping_edges$log_rate_norm, na.rm = T)))
max_color_val <- ceiling(max(max(hmm_mapping_edges$log_xi_norm, na.rm = T), 
                             max(phylo_mapping_edges$log_rate_norm, na.rm = T)))
max_abs_val <- max(abs(min_color_val), abs(max_color_val))

fig3_heatmap_options <- append(list(scale_fill_gradient2(low = "#5CB3FF", 
                                                         mid = "grey90", 
                                                         high = "darkorange2",
                                                         midpoint = 0, 
                                                         na.value = "grey96", 
                                                         limits = c(min_color_val,
                                                                    max_color_val))),
                               fig3_heatmap_options)
fig3_map_cluster_options <- append(list(scale_color_gradient2(low = "#5CB3FF", 
                                                              mid = "grey90", 
                                                              high = "darkorange2",
                                                              midpoint = 0, 
                                                              na.value = "grey60", 
                                                              limits = c(min_color_val,
                                                                         max_color_val))),
                                   fig3_map_cluster_options)

#-----------------------------------------------------------#
# Mapped edges----
#-----------------------------------------------------------#

## Phylo Mapped edges----
figure_3a_phylo_map <- ggplot() + 
  clust_map_background +
  geom_sf(data = phylo_consensus_clust_sf, aes(fill = cluster_assign)) +
  geom_sf(data = phylo_consensus_clust_sf |> filter(!is.na(cluster_assign)) |> 
            group_by(cluster_assign) |> summarize(),
          fill = "transparent", color = "grey10", linewidth = 0.4) + 
  geom_curve(data = phylo_filtered_edges,
             aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                 linewidth = log_rate_norm, 
                 alpha = log_rate_norm,
                 color = log_rate_norm), 
             curvature = -0.35, lineend = "round") +
  geom_point(data = st_centroid((full_africa_sf |> 
                                   filter(country %in% 
                                            (unique(c(phylo_filtered_edges$country1, phylo_filtered_edges$country2)))))$geometry),
             aes(x = unlist(map(geometry,1)), y = unlist(map(geometry,2))),
             size = 0.5, color = "grey20") +
  guides(linewidth = "none", 
         fill = "none",
         alpha = "none",
         color = guide_colorbar(title = "Standardized\nTransition Rates")) +
  fig3_map_cluster_options
figure_3a_phylo_map

## HMM Mapped edges----
figure_3c_hmm_map <- ggplot() + 
  clust_map_background + 
  geom_sf(data = hmm_consensus_clust_sf, aes(fill = cluster_assign)) +
  geom_sf(data = hmm_consensus_clust_sf |> filter(!is.na(cluster_assign)) |> 
            group_by(cluster_assign) |> summarize(),
          fill = "transparent", color = "grey30", linewidth = 0.4) + 
  geom_curve(aes(x = xstart, y = ystart, xend = xend, yend = yend,     
                 linewidth = log_xi_norm,
                 alpha = log_xi_norm, 
                 color = log_xi_norm), 
             data = hmm_filtered_edges,
             curvature = -0.35,
             lineend = "round") +
  geom_point(data = st_centroid(africa_sf$geometry),
             aes(x = unlist(map(geometry,1)), y = unlist(map(geometry,2))),
             size = 0.5, color = "grey20") +
  guides(linewidth = "none", 
         fill = "none",
         alpha = "none",
         color = guide_colorbar(title = "Standardized Connectivity\nMeasure")) +
  fig3_map_cluster_options
figure_3c_hmm_map

#-----------------------------------------------------------#
# Distance heatmaps----
#-----------------------------------------------------------#

countries_df <- tibble(country1 = rep(heatmap_order, times = length(heatmap_order)),
                       country2 = rep(heatmap_order, each = length(heatmap_order))) |> 
  filter(country1!=country2)

### Phylo distance heatmap----

phylo_dist_plotting_ordered <- countries_df |> 
  left_join(phylo_dist_plotting |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order_phylo, ordered = T),
         country2 = factor(country2, 
                           levels = rev(heatmap_order_phylo), ordered = T)) |> 
  arrange(country1, desc(country2))

# Object for outlining clusters
phylo_clust_box <- phylo_consensus_clusters |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order_phylo, ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(x = as.numeric(first(country1)), xend = as.numeric(last(country1))) |> 
  ungroup() |> 
  mutate(country1 = factor(country1, 
                           levels = rev(heatmap_order_phylo), ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(y = as.numeric(first(country1)), yend = as.numeric(last(country1))) |> 
  ungroup() |> 
  select(cluster_assign, x, xend, y, yend) |> 
  distinct() |> 
  arrange(cluster_assign)

figure_3b_phylo_heatmap <- ggplot() +
  geom_tile(data = phylo_dist_plotting_ordered, 
            aes(y = country2, x = country1, fill = log_rate_norm, 
                alpha = log_rate_norm), 
            color = "transparent", lwd = 0.01, linetype = 1) +
  geom_rect(data = phylo_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                        ymin = y-0.5, ymax = yend+0.5, 
                                        color = cluster_assign), fill = NA,
            lwd = 0.35) +
  fig3_heatmap_options +
  scale_color_manual(values = pub_cluster_palette)
figure_3b_phylo_heatmap


### HMM distance heatmap----
hmm_dist_plotting_ordered <- countries_df |> 
  left_join(hmm_dist_plotting |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor((country1), 
                           levels = heatmap_order, ordered = T),
         country2 = factor((country2), 
                           levels = rev(heatmap_order), ordered = T))

# Object for outlining clusters
hmm_clust_box <- hmm_consensus_clusters |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order, ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(x = as.numeric(first(country1)), xend = as.numeric(last(country1))) |> 
  ungroup() |> 
  mutate(country1 = factor(country1, 
                           levels = rev(heatmap_order), ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(y = as.numeric(first(country1)), yend = as.numeric(last(country1))) |> 
  ungroup() |> 
  select(cluster_assign, x, xend, y, yend) |> 
  distinct() |> 
  arrange(cluster_assign)

figure_3d_hmm_heatmap <- ggplot() +
  geom_tile(data = hmm_dist_plotting_ordered, 
            aes(y = country2, x = country1, fill = log_xi_norm, 
                alpha = log_xi_norm), 
            color = "transparent", lwd = 0.01, linetype = 1) +
  geom_rect(data = hmm_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                      ymin = y-0.5, ymax = yend+0.5, 
                   color = cluster_assign), fill = NA, 
            lwd = 0.35) +
  fig3_heatmap_options +
  scale_color_manual(values = pub_cluster_palette)

figure_3d_hmm_heatmap

#-----------------------------------------------------------#
# Save figure panels----
#-----------------------------------------------------------#

## Note: Full figure assembled in Adobe Illustrator
ggsave("Figure_3a_phylo_map.pdf", plot = figure_3a_phylo_map, 
       width = 2.7, height = 2.75, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Figure_3b_phylo_heatmap.pdf", plot = figure_3b_phylo_heatmap + theme(legend.position = "none"), 
       width = 2.4, height = 2.35, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Figure_3c_hmm_map.pdf", plot = figure_3c_hmm_map, 
       width = 2.7, height = 2.75, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Figure_3d_hmm_heatmap.pdf", plot = figure_3d_hmm_heatmap + theme(legend.position = "none"), 
       width = 2.4, height = 2.35, dpi = 400,
       path = here::here("figures/manuscript_figures"))

#-----------------------------------------------------------#
# Statistics for manuscript text----
#-----------------------------------------------------------#

library(coda)

### Phylo: Within vs between cluster analysis----

sep_sudan <- T
sampling <- "all_data"

source(here("analysis","00_functions_settings.R"), local = T)
phylo <- readRDS(here("data", "phylogeography.RDS"))
countries_subset <- c(countries, "south sudan")

distances <- st_centroid(full_africa_sf  |> 
                           filter(country %in% countries_subset)) |> 
  st_distance() |> 
  as.data.frame() |> 
  set_names((full_africa_sf  |> 
               filter(country %in% countries_subset))$country) |> 
  purrr::map_dfc( ~ {units::set_units(.x, km)}) 
distances <- distances |> 
  mutate(c = names(distances)) |> 
  mutate(m = purrr::map_int(c, ~which(countries_subset == .x))) |> 
  arrange(m) |> 
  dplyr::select(all_of(countries_subset)) |> 
  as.matrix()

dist <- data.frame(country2 = colnames(distances)) |>
  bind_cols(distances) |> 
  pivot_longer(cols = -country2, names_to = "country1", values_to = "distance") |> 
  list_countries()  |>  
  distinct()

phylo_dist <- phylo_mapping_edges |> 
  select(country1, country2, countries, log_rate, clust_rate, log_clust_rate, ind_pct) |> 
  filter(country1 %in% phylo_louvain_hclusters$country,
         country2 %in% phylo_louvain_hclusters$country) |> 
  left_join(dist |> select(countries, distance), by = c("countries")) |> 
  mutate(supported = factor(if_else(ind_pct>=0.46, 1, 0)),
         distance_s = distance/1000) |> 
  mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
  mutate(country1 = factor(country1, levels = heatmap_order, ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), ordered = T)) 

#### Distance vs support----
dist_support_model <- glm(supported ~ distance_s, family = "binomial"(link=logit), data = phylo_dist)
broom::tidy(dist_support_model, exponentiate = T, conf.int = T)

#### Distance vs transition rate----
dist_rate_model <- lm(clust_rate ~ distance_s, 
                      data = phylo_dist |> 
                        mutate(distance_s = if_else(distance_s>2.5, 2.5, distance_s)))
broom::tidy(dist_rate_model, exponentiate = F, conf.int = T)


#### Between vs within prep----
phylo_clust_compare <- phylo_dist |> 
  mutate_at(c("country1", "country2"), ~as.character(.)) |> 
  right_join(phylo_consensus_clusters |> select(country1, "country1_clust"=cluster_assign) |> 
              mutate(country1 = as.character(country1),
                     country1_clust = factor(country1_clust,
                                             labels = c("phyS", "phyCE", "phyW2", "phyW1"))), 
            by = "country1") |> 
  right_join(phylo_consensus_clusters |> select(country1, "country2_clust"=cluster_assign) |> 
              mutate(country1 = as.character(country1),
                     country2_clust = factor(country2_clust,
                                             labels = c("phyS", "phyCE", "phyW2", "phyW1"))), 
            by = c("country2"="country1")) |> 
  filter_out(is.na(country1)) |> 
  mutate(same_cluster = if_else(country1_clust==country2_clust, 1, 0),
         cluster_pair_summ = case_when(
           same_cluster==1 & str_detect(country1_clust,"phyW") ~ "Western",
           same_cluster==1 & str_detect(country1_clust,"CE|S") ~ "Eastern",
           (country1_clust=="phyS" | country2_clust=="phyS") &
             (country1_clust=="phyCE" | country2_clust=="phyCE") ~ "Eastern",
           str_detect(country1_clust,"phyW") & str_detect(country2_clust,"phyW") ~ "Western",
           (str_detect(country1_clust,"CE|S") | str_detect(country2_clust,"CE|S")) &
             (str_detect(country1_clust,"phyW") | str_detect(country2_clust,"phyW")) ~
             "Eastern_Western"),
         cluster_pair = case_when(
           same_cluster==1  ~ country1_clust,
            (country1_clust=="phyS" | country2_clust=="phyS") &
             (country1_clust=="phyCE" | country2_clust=="phyCE") ~ "CE_S",
           (country1_clust=="phyS" | country2_clust=="phyS") &
             (country1_clust=="phyW1" | country2_clust=="phyW1") ~ "W1_S",
           (country1_clust=="phyS" | country2_clust=="phyS") &
             (country1_clust=="phyW2" | country2_clust=="phyW2") ~ "W2_S",
           (country1_clust=="phyCE" | country2_clust=="phyCE") &
             (country1_clust=="phyW1" | country2_clust=="phyW1") ~ "W1_CE",
           (country1_clust=="phyCE" | country2_clust=="phyCE") &
             (country1_clust=="phyW2" | country2_clust=="phyW2") ~ "W2_CE",
           str_detect(country1_clust,"phyW") & str_detect(country2_clust,"phyW") ~ "W1_W2"),
         rate = exp(log_rate)) |> 
  distinct() |> 
  group_by(same_cluster) |> 
  mutate(mean_rate_same = mean(rate),
         med_rate_same = median(rate)) |> 
  ungroup()

full_phylo_clust <- phylo |> 
  mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
  group_by(countries) |> 
  mutate(ind_pct = sum(location_indicators)/n()) |> 
  ungroup() |> 
  mutate(rate_ind = location_rates*location_indicators,
         rate = location_rates*ind_pct) |> 
  left_join(phylo_clust_compare |> select(country1, country2, same_cluster,
                                          country1_clust, country2_clust,
                                          cluster_pair_summ), 
            by = c("country1"="country2", "country2"="country1")) |> 
  filter(!str_detect(countries, "Other"))

#### Same vs different clusters----

compare_samevdiff <- full_phylo_clust |> 
  group_by(state, same_cluster) |> 
  summarize(mean_rate = mean(rate)) |> 
  ungroup() |> 
  group_by(state) |> 
  mutate(rate_compare = if_else(same_cluster==1, mean_rate/lag(mean_rate),
                                as.numeric(NA))) |> 
  fill(rate_compare, .direction = "updown") |> 
  ungroup() |> 
  select(state, rate_compare) |> 
  distinct() 

compare_samevdiff_mcmc <- as.mcmc(compare_samevdiff$rate_compare)
compare_samevdiff_hpd <- HPDinterval(compare_samevdiff_mcmc, prob = 0.95)

compare_samevdiff_hpd <- c(mean(compare_samevdiff$rate_compare), compare_samevdiff_hpd)
names(compare_samevdiff_hpd) <- c("mean", "lower", "upper")
compare_samevdiff_hpd

  
#### Between cluster analysis----
  # Comparing western vs CE/S clusters using weighted rates from heatmap

compare_between <- full_phylo_clust |> 
  group_by(state, same_cluster, cluster_pair_summ) |> 
  summarize(mean_rate_pair_summ = mean(rate))  |> 
  ungroup() |> 
  select(state, same_cluster, cluster_pair_summ, mean_rate_pair_summ) |> 
  distinct() |> 
  filter(same_cluster==0) |> 
  select(-same_cluster) |> 
  distinct() |> 
  pivot_wider(id_cols = c(state),
              values_from = mean_rate_pair_summ,
              names_from = cluster_pair_summ) |> 
  fill(c(Eastern_Western, Western, Eastern), .direction = "downup") |> 
  mutate(e_vs_ew = Eastern/Eastern_Western,
         w_vs_ew = Western/Eastern_Western)

# Between phyW1 & phyW2 vs between Eastern (phyCE/phyS) & Western (phyW1/phyW2)
compare_between_w_vs_ew_mcmc <- as.mcmc(compare_between$w_vs_ew)
compare_between_w_vs_ew_hpd <- HPDinterval(compare_between_w_vs_ew_mcmc, prob = 0.95)
compare_between_w_vs_ew_hpd <- round(c(mean(compare_between$w_vs_ew), compare_between_w_vs_ew_hpd),1)
names(compare_between_w_vs_ew_hpd) <- c("mean", "lower", "upper")
compare_between_w_vs_ew_hpd

# Between phyCE & phyS vs between Eastern (phyCE/phyS) & Western (phyW1/phyW2)
compare_between_e_vs_ew_mcmc <- as.mcmc(compare_between$e_vs_ew)
compare_between_e_vs_ew_hpd <- HPDinterval(compare_between_e_vs_ew_mcmc, prob = 0.95)
compare_between_e_vs_ew_hpd <- round(c(mean(compare_between$e_vs_ew), compare_between_e_vs_ew_hpd),1)
names(compare_between_e_vs_ew_hpd) <- c("mean", "lower", "upper")
compare_between_e_vs_ew_hpd


## HMM: Within vs between cluster analysis----

sep_sudan <- F
sampling <- "all_data"

source(here("analysis","00_functions_settings.R"), local = T)

distances <- st_centroid(full_africa_sf  |> 
                           filter(country %in% countries_subset)) |> 
  st_distance() |> 
  as.data.frame() |> 
  set_names((full_africa_sf  |> 
               filter(country %in% countries_subset))$country) |> 
  purrr::map_dfc( ~ {units::set_units(.x, km)}) 
distances <- distances |> 
  mutate(c = names(distances)) |> 
  mutate(m = purrr::map_int(c, ~which(countries_subset == .x))) |> 
  arrange(m) |> 
  dplyr::select(all_of(countries_subset)) |> 
  as.matrix()

dist <- data.frame(country2 = colnames(distances)) |>
  bind_cols(distances) |> 
  pivot_longer(cols = -country2, names_to = "country1", values_to = "distance") |> 
  list_countries()  |>  
  distinct()

hmm_louvain_hclusters_rename <- hmm_louvain_hclusters |> 
  mutate(country = country_rename(country))

hmm_dist <- hmm_mapping_edges |> 
  select(country1, country2, countries, log_xi_norm) |> 
  filter(country1 %in% hmm_louvain_hclusters_rename$country,
         country2 %in% hmm_louvain_hclusters_rename$country) |> 
  left_join(dist |> 
              select(countries, distance), 
            by = c("countries")) |> 
  mutate_at(c("country1", "country2"), ~country_rename(.)) |> 
  mutate(country1 = factor(country1, levels = heatmap_order, ordered = T),
         country2 = factor(country2, levels = rev(heatmap_order), ordered = T)) 


hmm_clust_compare <- hmm_dist |> 
  mutate_at(c("country1", "country2"), ~as.character(.)) |> 
  right_join(hmm_consensus_clusters |> select(country1, "country1_clust"=cluster_assign) |> 
               mutate(country1 = as.character(country1),
                      country1_clust = factor(country1_clust,
                                              labels = c("hmmS", "hmmCE", "hmmW"))), 
             by = "country1") |> 
  right_join(hmm_consensus_clusters |> select(country1, "country2_clust"=cluster_assign) |> 
               mutate(country1 = as.character(country1),
                      country2_clust = factor(country2_clust,
                                              labels = c("hmmS", "hmmCE", "hmmW"))), 
             by = c("country2"="country1")) |> 
  filter(!is.na(country1)) |> 
  mutate(same_cluster = if_else(country1_clust==country2_clust, 1, 0),
         cluster_pair_summ = case_when(
           str_detect(country1_clust,"hmmW") ~ "Western",
           (country1_clust=="hmmS" | country2_clust=="hmmS") &
             (country1_clust=="hmmCE" | country2_clust=="hmmCE") ~ "Eastern",
           (str_detect(country1_clust,"CE|S") | str_detect(country2_clust,"CE|S")) &
             (str_detect(country1_clust,"hmmW") | str_detect(country2_clust,"hmmW")) ~
             "Eastern_Western"),
         cluster_pair = case_when(
           same_cluster==1  ~ country1_clust,
           (country1_clust=="hmmS" | country2_clust=="hmmS") &
             (country1_clust=="hmmCE" | country2_clust=="hmmCE") ~ "CE_S",
           (country1_clust=="hmmS" | country2_clust=="hmmS") &
             (country1_clust=="hmmW" | country2_clust=="hmmW") ~ "W_S",
           (country1_clust=="hmmCE" | country2_clust=="hmmCE") &
             (country1_clust=="hmmW" | country2_clust=="hmmW") ~ "W_CE")) |> 
  group_by(same_cluster) |> 
  mutate(mean_rate_same = mean(log_xi_norm),
         med_rate_same = median(log_xi_norm)) |> 
  ungroup() |> 
  group_by(cluster_pair, same_cluster) |> 
  mutate(mean_rate_pair =  mean(log_xi_norm),
         med_rate_pair =  median(log_xi_norm)) |> 
  ungroup() |> 
  group_by(cluster_pair_summ) |> 
  mutate(mean_rate_pair_summ = mean(log_xi_norm),
         med_rate_pair_summ = median(log_xi_norm)) |> 
  ungroup() |> 
  group_by(cluster_pair_summ, same_cluster) |> 
  mutate(mean_rate_pair_summ2 = mean(log_xi_norm),
         med_rate_pair_summ2 = median(log_xi_norm)) |> 
  ungroup()


hmm_within_between_model <- glm(same_cluster ~ log_xi_norm, 
                                family = "binomial"(link=logit), 
                                data = hmm_clust_compare)
broom::tidy(hmm_within_between_model, exponentiate = T, conf.int = T)




