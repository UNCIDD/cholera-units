#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 15_Figure5.R
# Purpose: This script creates Figure 5 mapping arrival times
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(sf)
library(here)
library(plotly)
library(cowplot)

options(scipen = 9999)

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

set.seed(seed)

#-----------------------------------------------------------#
# Read data----
#-----------------------------------------------------------#

source(here("analysis","02_import_data.R"), local = T) 

downstream_pred <- readRDS(str_c(sim_dir, "/downstream_pred.rds"))
sim_downstream <- readRDS(str_c(sim_dir, "/sim_downstream.rds"))

sim_secondary <- readRDS(str_c(sim_dir, "/sim_secondary.rds")) 

hmm_dist_plotting <- readRDS(str_c(clust_dir, "/xi_dist_plotting.rds"))
hmm_consensus_clusters <- readRDS(str_c(clust_dir, "/hmm_consensus_clusters.rds"))

#-----------------------------------------------------------#
# Spatial objects----
#-----------------------------------------------------------#

countries <- all_countries
M <- length(countries) 
n_time <- 1
spatial_objs <- create_connections(M, n_time)
nodes <- make_nodes()
edges <- make_edges(nodes)

#-----------------------------------------------------------#
# Prep data for plotting----
#-----------------------------------------------------------#

## Countries with highest risks----
# For consistency with Figure 4
secondary_risk <- sim_secondary |> 
  purrr::map_depth(2, ~{
    .x |> 
      filter(year==2) |> 
      select(country, Presence)}) |> 
  purrr::map(~{
    .x |> 
      bind_rows(.id = "draw") |> 
      group_by(country) |> 
      mutate(pct = sum(Presence)/n()) |> 
      ungroup() |> 
      select(country, pct) |> 
      distinct()
  }) |> 
  bind_rows(.id = "seed_country")

secondary_high_risk <- secondary_risk |> 
  group_by(seed_country) |> 
  mutate(risk = if_else(seed_country==country, 0, pct),
         max_risk = max(risk)) |> 
  ungroup() |> 
  distinct(seed_country, max_risk) |> 
  filter(max_risk>0.21) |> 
  distinct(seed_country) |> 
  pull()

## Arrival time distributions----
# downstream_pred |> 
#   filter(seed_country == "ghana") |> 
#   ggplot() +
#   geom_histogram(aes(x = arrival_time)) +
#   facet_wrap(~country) +
#   theme_classic()
arrival_time_summ <- downstream_pred |> 
  mutate(arrival_time = if_else(arrival_time>10,10,arrival_time)) |> 
  group_by(seed_country, country) |> 
  summarize(median_arrival = median(arrival_time),
            mean_arrival = mean(arrival_time)) |> 
  ungroup()


## Shapefiles----
all_downstream_sfs <- countries |> 
  country_rename(sudan_rename = T) |> 
  set_names() |> 
  purrr::map(~{
    africa_sf |> 
      mutate(country = country_rename(country, sudan_rename = T)) |> 
      left_join(arrival_time_summ |> 
                  mutate_at(c("country", "seed_country"), ~country_rename(., sudan_rename = T)) |> 
                  filter(seed_country==.x),
                by = "country")
  })
    
hmm_consensus_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country, sudan_rename = T)) |> 
  left_join(hmm_consensus_clusters |> 
              mutate(country1 = country_rename(country1, sudan_rename = T)),
            by = c("country"="country1")) 

## Consensus clusters----
# Ensure correct ordering
consensus_clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Cameroon"),
  cluster_hmm3k_assign = factor(1:3))  |> 
  left_join(hmm_consensus_clusters |> 
              mutate(country1 = country_rename(country1, sudan_rename = T)) |> 
              select(country1, hmm_cluster = cluster),
            by = c("country"="country1")) 
hmm_consensus_clusters2 <- hmm_consensus_clusters |> 
  left_join(consensus_clust_reassign |> select(contains("hmm")) |> distinct(),
            by = c("cluster" = "hmm_cluster")) |> 
  rename(cluster_assign = cluster_hmm3k_assign) |> 
  mutate(cluster_assign = factor(cluster_assign))

# Object for drawing boxes around clusters in heat map
hmm_clust_box <- hmm_consensus_clusters2 |> 
  mutate(country1 = factor(country_rename(country1, sudan_rename = T), 
                           levels = heatmap_order_comb, ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(x = as.numeric(first(country1)), xend = as.numeric(last(country1))) |> 
  ungroup() |> 
  mutate(country1 = factor(country1, 
                           levels = rev(heatmap_order_comb), ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(y = as.numeric(first(country1)), yend = as.numeric(last(country1))) |> 
  ungroup() |> 
  select(cluster_assign, x, xend, y, yend) |> 
  distinct() |> 
  arrange(cluster_assign)

## Standardize names----
full_africa_sf <- full_africa_sf |> 
  mutate(country = country_rename(country, sudan_rename = T))
secondary_high_risk <- country_rename(secondary_high_risk, sudan_rename = T)

#-----------------------------------------------------------#
# Figure 5 maps ----
#-----------------------------------------------------------#

figure5_maplist <- map(heatmap_order_comb[which(heatmap_order_comb %in% 
                                          secondary_high_risk)], ~
                      {
                        ggplot() + 
                          geom_sf(data = full_africa_sf, fill = "grey80") +
                          geom_sf(data = all_downstream_sfs[[which(names(all_downstream_sfs)==.x)]], 
                                  aes(fill = median_arrival), alpha = 1) +
                          geom_sf(data = full_africa_sf, fill = "transparent", color="grey70") +
                          geom_sf(data = hmm_consensus_clust_sf |> filter(!is.na(cluster)) |> 
                                    group_by(cluster) |> summarize(),
                                  fill = "transparent", color = "black", linewidth = 0.45) + 
                          theme_void() +
                          ggtitle(.x) +
                          scale_fill_gradient2(low = "#15317E", mid = "#A0CFEC", 
                                               high = "white", midpoint = 5,
                                               limits = c(0,10)) +
                          guides(fill = guide_colorbar(title = "Time to Arrival (years)")) +
                          theme(legend.position = "bottom",
                                legend.text = element_text(size = 7),
                                legend.title = element_text(size = 7.5),
                                legend.key.height = unit(0.25, "cm"),
                                plot.title.position = 'plot', 
                                plot.title = element_text(size = 8, hjust = 0.45)) }) |> 
  set_names(heatmap_order_comb[which(heatmap_order_comb %in% 
                                  secondary_high_risk)])

#Get legend
figure5_legend <- get_plot_component(figure5_maplist[[1]], "guide-box-bottom", 
                                     return_all = T)

figure5_maplist <- map(figure5_maplist, ~{.x + theme(legend.position = "none")})


#-----------------------------------------------------------#
# Figure 5 heatmap ----
#-----------------------------------------------------------#


heatmap_options <- list(scale_fill_gradient2(low = "#15317E", mid = "#A0CFEC", 
                                             high = "white",
                                             midpoint = 5,
                                             limits = c(0,10), na.value = "grey80"),
                        theme_classic(),
                        theme(text = element_text(size = 6,
                                                  family = "serif"),
                              axis.title.y = element_text(size = 7),
                              axis.title.x = element_blank(),
                              axis.text.x = element_text(angle = 90, 
                                                         hjust = 1,
                                                         margin = margin(t = 0.01)),
                              axis.text.y = element_text(margin = margin(r = 0.01)),
                              axis.ticks = element_line(linewidth = 0.01),
                              axis.ticks.length = unit(0.05,"cm"),
                              axis.line = element_line(linewidth = 0.05),
                              legend.position = "none",
                              plot.margin = unit(c(t = 0.1, r = 0.1, 
                                                   b = 0.1, l = 0.01), "cm")),
                        coord_fixed())

figure5_spread_heatmap <- tibble(seed_country = rep(heatmap_order_comb, 
                                                    times = length(heatmap_order_comb)),
                                 country = rep(heatmap_order_comb, 
                                               each = length(heatmap_order_comb))) |> 
  filter(seed_country!=country) |> 
  left_join(arrival_time_summ |>  
              mutate_at(c("seed_country", "country"), 
                        ~country_rename(., sudan_rename = T))) |> 
  mutate(seed_country = factor(seed_country, 
                               levels = rev(heatmap_order_comb), ordered = T),
         country = factor(country, 
                          levels = (heatmap_order_comb), ordered = T)) |> 
  ggplot() +
  geom_tile(aes(y = seed_country, x = country, fill = median_arrival), 
            color = "white", lwd = 0.05, linetype = 1) +
  #Outline clusters
  geom_rect(data = hmm_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                      ymin = y-0.5, ymax = yend+0.5, 
                                      color = cluster_assign), fill = NA, 
            lwd = 0.35) +
  theme_classic() +
  heatmap_options +
  ylab("Country with an introduction") +
  scale_color_manual(values = pub_cluster_palette)
figure5_spread_heatmap 


#-----------------------------------------------------------#
# Combine plots ----
#-----------------------------------------------------------#

figure5_grid1 <- cowplot::plot_grid(plotlist = figure5_maplist[1:10], ncol = 5)
figure5_grid2 <- plot_grid(plot_grid(plotlist = figure5_maplist[11:14], ncol = 2),
                           figure5_spread_heatmap + theme(legend.position = "none"), 
                           nrow = 1, rel_widths = c(0.7,1))
figure_5 <- plot_grid(plot_grid(figure5_grid1, figure5_grid2,
                                nrow = 2, rel_heights = c(1,1)),
                      figure5_legend, rel_heights = c(1,0.08), ncol = 1)
figure_5


#-----------------------------------------------------------#
# Supplemental Figure S5 maps ----
#-----------------------------------------------------------#

supp_fig5_maplist <- map(heatmap_order_comb, ~
                           {
                             ggplot() + 
                               geom_sf(data = full_africa_sf, fill = "grey80") +
                               geom_sf(data = all_downstream_sfs[[which(names(all_downstream_sfs)==.x)]], 
                                       aes(fill = median_arrival), alpha = 1) +
                               geom_sf(data = full_africa_sf, fill = "transparent", color="grey70") +
                               geom_sf(data = hmm_consensus_clust_sf |> filter(!is.na(cluster)) |> 
                                         group_by(cluster) |> summarize(),
                                       fill = "transparent", color = "black", linewidth = 0.45) + 
                               theme_void() +
                               ggtitle(.x) +
                               scale_fill_gradient2(low = "#15317E", mid = "#A0CFEC", 
                                                    high = "white", midpoint = 5,
                                                    limits = c(0,10)) +
                               guides(fill = guide_colorbar(title = "Time to Arrival (years)")) +
                               theme(legend.position = "none",
                                     legend.text = element_text(size = 7),
                                     legend.title = element_text(size = 7.5),
                                     legend.key.height = unit(0.25, "cm"),
                                     plot.title.position = 'plot', 
                                     plot.title = element_text(size = 8, hjust = 0.45)) }) |> 
  set_names(heatmap_order_comb)



supp_fig5_grid <- cowplot::plot_grid(plotlist = supp_fig5_maplist, ncol = 6)
supp_fig5 <- plot_grid(supp_fig5_grid,
                      figure5_legend, rel_heights = c(1,0.08), ncol = 1)
supp_fig5


#-----------------------------------------------------------#
# Save plots----
#-----------------------------------------------------------#

ggsave("Figure_5.pdf", plot = figure_5, height = 8, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Figure_5.png", plot = figure_5, height = 8, width = 7, dpi = 400,
       bg = "white", path = here::here("figures/manuscript_figures"))

ggsave("Supplemental_FigS5.pdf", plot = supp_fig5, height = 9, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_Fig5.png", plot = supp_fig5, height = 9, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))

#-----------------------------------------------------------#
# Statistics for manuscript text ----
#-----------------------------------------------------------#

arrival_cluster <- arrival_time_summ |> 
  left_join(hmm_consensus_clusters, by = c("country"="country1")) |> 
  left_join(hmm_consensus_clusters |> rename(seed_cluster = cluster), by = c("seed_country"="country1")) |> 
  mutate_at(c("country", "seed_country"), ~country_rename(., sudan_rename = T))

cluster_borders <- tibble(cluster1 = rep(1:3, each = 2),
                          cluster2 = c(2,3,1,3,1,2),
                          border_countries = c(str_c("Chad", "Central African Rep.", 
                                                     "Rep. of the Congo", sep = ", "),
                                               "Rep. of the Congo",
                                               str_c("Sudan", "DRC", sep = ", "),
                                               str_c("DRC", "Zambia", "Tanzania", sep = ", "),
                                               "Angola",
                                               str_c("Angola", "Namibia", "Botswana", 
                                                     "Zimbabwe", "Mozambique", "Malawi", sep = ", ")))

arrival_cluster2 <- arrival_cluster |> 
  group_by(seed_country) |> 
  mutate(same = if_else(cluster==seed_cluster, 1, 0),
         southeast = if_else(same!=1 & ((cluster==2 & seed_cluster==3)
                                        | (cluster==3 & seed_cluster==2)), 1, 0)) |> 
  ungroup() |> 
  group_by(seed_country, same, southeast) |> 
  mutate(cross_border_time = min(mean_arrival)) |> 
  ungroup() 


hmm_consensus_clusters2 <- hmm_consensus_clusters2 |> 
  mutate(country1 = country_rename(country1))

downstream_clust <- downstream_pred  |> 
  filter(seed_country!=country) |> 
  mutate_at(c("seed_country", "country"), ~country_rename(.)) |> 
  right_join(hmm_consensus_clusters2 |> select(country1, "country1_clust"=cluster_assign) |> 
               mutate(country1 = as.character(country1),
                      country1_clust = factor(country1_clust,
                                              levels = c(1:3),
                                              labels = c("hmmS", "hmmCE", "hmmW"))), 
             by = c("seed_country"="country1"))  |> 
  right_join(hmm_consensus_clusters2 |> select(country1, "country2_clust"=cluster_assign) |> 
               mutate(country1 = as.character(country1),
                      country2_clust = factor(country2_clust,
                                              levels = c(1:3),
                                              labels = c("hmmS", "hmmCE", "hmmW"))), 
             by = c("country"="country1")) |> 
  mutate(same_cluster = (if_else(country1_clust==country2_clust,1,0)),
         cluster = factor(case_when(
           same_cluster==1  ~ country1_clust,
           (country1_clust=="hmmS" | country2_clust=="hmmS") &
             (country1_clust=="hmmCE" | country2_clust=="hmmCE") ~ "CE-S",
           (country1_clust=="hmmS" | country2_clust=="hmmS") &
             (country1_clust=="hmmW" | country2_clust=="hmmW") ~ "W-S",
           (country1_clust=="hmmCE" | country2_clust=="hmmCE") &
             (country1_clust=="hmmW" | country2_clust=="hmmW") ~ "W-CE")))

between_vs_in_arrival <- downstream_clust |> 
  group_by(same_cluster) |> 
  mutate(mean_time = mean(arrival_time)) |> 
  group_by(cluster) |> 
  mutate(mean_time_cluster = mean(arrival_time)) |> 
  ungroup() |> 
  select(same_cluster, cluster, mean_time, mean_time_cluster) |> 
  distinct()
between_vs_in_arrival |> select(same_cluster, mean_time) |> distinct()

# average time to arrival within vs between units
broom::tidy(lm(arrival_time ~ same_cluster, data = downstream_clust), conf.int = T)

cluster_downstream <- downstream_clust |> 
  select(same_cluster, country1_clust, country2_clust, arrival_time) |> 
  group_by(country1_clust, country2_clust) |> 
  mutate(mean_time = mean(arrival_time)) |> 
  ungroup() |> 
  select(-arrival_time) |> 
  distinct()
cluster_downstream

#Timing of spread throughout same transmission unit - Ghana
downstream_clust |> 
  filter(seed_country=="Ghana") |> 
  group_by(country) |> 
  mutate(avg_time = mean(arrival_time)) |> 
  ungroup() |> 
  filter(same_cluster==1) |> 
  distinct(seed_country, country, avg_time) |> 
  arrange(avg_time) |> 
  print(n=21)

summary((downstream_clust |> 
           filter(seed_country=="Ghana") |> 
           group_by(country) |> 
           mutate(avg_time = mean(arrival_time)) |> 
           ungroup() |> 
           filter(same_cluster==1) |> 
           distinct(seed_country, country, avg_time)
           )$avg_time)

#Timing of spread throughout same transmission unit - Cameroon
downstream_clust |> 
  filter(seed_country=="Cameroon") |> 
  group_by(country) |> 
  mutate(avg_time = mean(arrival_time)) |> 
  ungroup() |> 
  filter(same_cluster==1) |> 
  distinct(seed_country, country, avg_time) |> 
  arrange(avg_time) |> 
  print(n=21)

summary((downstream_clust |> 
           filter(seed_country=="Cameroon") |> 
           group_by(country) |> 
           mutate(avg_time = mean(arrival_time)) |> 
           ungroup() |> 
           filter(same_cluster==1) |> 
           distinct(seed_country, country, avg_time)
)$avg_time)
  