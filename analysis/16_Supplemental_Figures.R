#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 16_Figure4.R
# Purpose: This script creates the following supplemental material:
# 
# Fig. S1. Proportion of iterations where country-pairs clustered together
# Fig. S2. Inferred probability of presence of cholera lineages over time by country
# Fig. S3. Inferred cholera transmission units and cross-country transmission in sub-Saharan Africa
# Fig. S6. Lineage prevalence.
# Table S3
#
# Note the files where additional supplemental material is created:
#
# Fig. S4: file 14_Figure4.R
# Fig. S5: file 15_Figure5.R
# Fig. S7: file 10_arrival_times.R
# Table S2: file 14_Figure4.R
#
#-----------------------------------------------------------------#

#-----------------------------------------------------------#
# Initialize----
#-----------------------------------------------------------#

library(tidyverse)
library(sf)
library(here)
library(plotly)
library(cowplot)
library(ggpattern)

options(scipen = 9999)

# Source functions
sampling <- "all_data"

source(here("analysis","00_functions_settings.R"), local = T)
source(here("analysis","02_import_data.R"), local = T)

set.seed(seed)

filter_cutoff <- 0.5

countries <- countries_subset

#-----------------------------------------------------------#
# Fig. S1 - Louvain clustering----
#-----------------------------------------------------------#

### Read data----

phylo_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_phylo_subgrps.rds"))
hmm_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_hmm_subgrps.rds"))

### Plot settings----

louvain_heatmap_options <- list(theme_classic(),
                        theme(axis.text.x = element_text(angle = 90, hjust = 0.98),
                              axis.text.x.bottom = element_text(vjust = 0.5),
                              text = element_text(size = 7,
                                                  family = "serif"),
                              axis.title = element_blank()),
                        coord_fixed())

countries_df <- tibble(country1 = rep(heatmap_order, times = length(heatmap_order)),
                       country2 = rep(heatmap_order, each = length(heatmap_order))) |> 
  filter(country1!=country2)

### Country ordering----
phylo_louvain_subgrps_ordered <- countries_df |> 
  left_join(phylo_louvain_subgrps |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order_phylo, ordered = T),
         country2 = factor(country2, 
                           levels = rev(heatmap_order_phylo), ordered = T)) |> 
  arrange(country1, desc(country2))
hmm_louvain_subgrps_ordered <- countries_df |> 
  left_join(hmm_louvain_subgrps |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order, ordered = T),
         country2 = factor(country2, 
                           levels = rev(heatmap_order), ordered = T)) |> 
  arrange(country1, desc(country2))

### Fig S1a - Phylo Louvain clustering heatmap----
figS1a_louvain_phylo_heatmap <- ggplot(data = phylo_louvain_subgrps_ordered, 
                                aes(y = country2, x = country1)) +
  geom_tile(aes(fill = pct), color = "white",
            lwd = 0.25,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#F4FAC4", mid = "#79CDCD", high = "#16256A",
                       midpoint = 0.5, limits = c(0,1), na.value = "grey90") +
  louvain_heatmap_options +
  guides(fill = guide_colorbar(title = 
                                 "Proportion of\niterations where countries\nclustered together")) +
  theme(legend.position = "bottom")
figS1a_louvain_phylo_heatmap

### Fig S1b - HMM Louvain clustering heatmap----
figS1b_louvain_hmm_heatmap <- ggplot(data = hmm_louvain_subgrps_ordered, 
                              aes(y = country2, x = country1)) +
  geom_tile(aes(fill = pct), color = "white",
            lwd = 0.25,
            linetype = 1) +
  theme_classic() +
  scale_fill_gradient2(low = "#F4FAC4", mid = "#79CDCD", high = "#16256A", 
                       midpoint = 0.5, limits = c(0,1), na.value = "grey90") +
  louvain_heatmap_options
figS1b_louvain_hmm_heatmap

### Combine plots----
louvain_heatmaps_legend <- get_plot_component(figS1a_louvain_phylo_heatmap, "guide-box-bottom", 
                                      return_all = T)

figS1_louvain_heatmaps <- cowplot::plot_grid(figS1a_louvain_phylo_heatmap + 
                                         theme(legend.position = "none"), 
                                         figS1b_louvain_hmm_heatmap + 
                                         theme(legend.position = "none"),
                                       ncol = 2, labels = c("A.", "B."),
                                       label_fontfamily = "serif", 
                                       label_size = 10)
figS1_louvain_heatmaps <- cowplot::plot_grid(figS1_louvain_heatmaps, louvain_heatmaps_legend,
                                       rel_heights = c(1.1, 0.12), nrow = 2)
figS1_louvain_heatmaps

#-----------------------------------------------------------#
# Fig. S2 - Probability of presence----
#-----------------------------------------------------------#

## Read data----
pred_rho <- readRDS(str_c(xi_dir, "/pred_rho.rds")) |> 
  translate_cols() |> 
  select(country, year, te, rho = mean, sd) |> 
  left_join(cases_lineages_full |> select(country, year, te, samples),
            by = c("country", "year", "te")) |> 
  mutate(country = country_rename(country),
         te = factor(te, levels = te_order, labels = pub_te_order,
                     ordered = T),
         present = if_else(!is.na(samples), as.numeric(te)/16, samples))

## Fig S2 - Plot probabilities----
figS2_rho_plot <- ggplot(pred_rho, aes(x = year)) +
  geom_line(aes(y = rho, color = te), position = position_dodge(0.2)) +
  ggside::geom_xsidepoint(aes(y = present, color = te), alpha = 0.8, shape = 16,
                          size = 0.4) +
  facet_wrap(~country, ncol = 5) +
  theme_bw() +
  scale_color_manual(values = pub_te_palette, drop = T) +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = c(0,1.05), breaks = seq(0,1,0.5),
                     minor_breaks = seq(0,1,0.25)) +
  ggside::scale_xsidey_continuous(limits = c(0.1, max(pred_rho$present + 0.1, na.rm = T)), minor_breaks = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  ggside::theme_ggside_void() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(size = 6.5,
                                  margin = margin(0.03,0,0.03,0, "cm")),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 7),
        legend.key.width = unit(0.25,"cm"),
        legend.key.height = unit(0.4,"cm"),
        legend.key.spacing = unit(1, "mm"),
        axis.title = element_text(size = 7.5),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 6),
        legend.margin = margin(0.04, 0.0, 0.02, 0.0, "cm"),
        ggside.axis.text.y = element_blank(),
        ggside.panel.scale = 0.2)  +
  ylab("Probability of Presence")
figS2_rho_plot

#-----------------------------------------------------------#
# Map & heatmap options----
#-----------------------------------------------------------#

phylo_mapping_edges <- readRDS(str_c(clust_dir, "/phylo_mapping_edges.rds"))
hmm_mapping_edges <- readRDS(str_c(clust_dir, "/hmm_mapping_edges.rds"))

min_color_val <- floor(min(min(hmm_mapping_edges$log_xi_norm, na.rm = T), 
                           min(phylo_mapping_edges$log_rate_norm, na.rm = T)))
max_color_val <- ceiling(max(max(hmm_mapping_edges$log_xi_norm, na.rm = T), 
                             max(phylo_mapping_edges$log_rate_norm, na.rm = T)))
max_abs_val <- max(abs(min_color_val), abs(max_color_val))

## Background layers for maps----
clust_map_background <- list(
  geom_sf_pattern(data = full_africa_sf |> filter(country %in% north_africa_exclude), 
                  fill = "grey90", color = "grey40", alpha = 0.2,
                  pattern_alpha = 0.2, pattern_spacing = 0.01, pattern_fill = "grey70",
                  pattern_colour = "grey30", pattern_size = 0.02),
  geom_sf(data = full_africa_sf |> filter(country %in% c(countries, "South Sudan")), 
          fill = "white", color = "grey40", alpha = 0.8)
)

countries_df <- tibble(country1 = rep(heatmap_order, times = length(heatmap_order)),
                       country2 = rep(heatmap_order, each = length(heatmap_order))) |> 
  filter(country1!=country2)

#-----------------------------------------------------------#
# Fig. S3 A&B - Post 1990 Transmission Units----
#-----------------------------------------------------------#

## Read data----
sampling <- "post_1990"

source(here("analysis","00_functions_settings.R"), local = T)

hmm_dist_plotting <- readRDS(str_c(clust_dir, "/xi_dist_plotting.rds"))
hmm_louvain_hclusters <- readRDS(str_c(clust_dir, "/hmm_louvain_hclusters.rds"))
hmm_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_hmm_subgrps.rds"))
hmm_mapping_edges <- readRDS(str_c(clust_dir, "/hmm_mapping_edges.rds"))
hmm_consensus_clusters <- readRDS(str_c(clust_dir, "/hmm_consensus_clusters.rds"))

## Update mapping options----
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

## Spatial objects----
countries <- all_countries
M <- length(countries) 
n_time <- 1
spatial_objs <- create_connections(M, n_time)
nodes <- make_nodes()
edges <- make_nodes(nodes)

## Prep for plotting----
hmm_filter_stat <- quantile(hmm_mapping_edges$log_xi_norm, filter_cutoff)
hmm_filter_stat_mean <- mean(hmm_mapping_edges$log_xi_norm)

hmm_filtered_edges <- hmm_mapping_edges |> 
  filter(log_xi_norm>hmm_filter_stat_mean)

hmm_dist_plotting <- hmm_dist_plotting |> 
  select(country1, country2, log_xi_norm) |> 
  mutate_at(c("country1", "country2"),  ~as.character(.))

## Consistent clustering----
clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Mauritania", "Central African Rep."),
  cluster_4k_assign = factor(c(1:4))) |> 
  left_join(hmm_louvain_hclusters |> select(country, hmm_cluster = cluster_4k),
            by = "country") 

hmm_louvain_hclusters <- hmm_louvain_hclusters |> 
  left_join(clust_reassign |> select(hmm_cluster, cluster_4k_assign),
            by = c("cluster_4k"="hmm_cluster"))
consensus_clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Cameroon", "Mauritania"),
  cluster_hmm3k_assign = factor(c(1:3, 3)))  |> 
  left_join(hmm_consensus_clusters |> select(country1, hmm_cluster = cluster),
            by = c("country"="country1")) 

## Shapefiles for mapping----
hmm_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(hmm_louvain_hclusters,
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

## Fig S3A HMM Mapped edges----
hmm_map_post_1990 <- ggplot() + 
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
hmm_map_post_1990

## Fig S3B HMM distance heatmap----

hmm_dist_plotting_ordered <- countries_df |> 
  left_join(hmm_dist_plotting |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor((country1), 
                           levels = heatmap_order, ordered = T),
         country2 = factor((country2), 
                           levels = rev(heatmap_order), ordered = T))
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

hmm_heatmap_post_1990 <- ggplot() +
  geom_tile(data = hmm_dist_plotting_ordered, 
            aes(y = country2, x = country1, fill = log_xi_norm, 
                alpha = log_xi_norm), 
            color = "transparent", lwd = 0.01, linetype = 1) +
  geom_rect(data = hmm_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                      ymin = y-0.5, ymax = yend+0.5, 
                                      color = cluster_assign), fill = NA, 
            lwd = 0.5) +
  fig3_heatmap_options +
  scale_color_manual(values = pub_cluster_palette)
hmm_heatmap_post_1990

#-----------------------------------------------------------#
# Fig. S6 A - Post 1990 Assigned prevalence----
#-----------------------------------------------------------#

## Bar/donut plot: Assigned prevalence based on sequences

## Read data----
cases_prev_pred <- readRDS(str_c(xi_dir, "/cases_prev_pred.rds")) |> 
  clean_df() |> 
  mutate(te = factor(te, levels = te_order, labels = pub_te_order))

## Bar graph----
prev1990_bar <- cases_prev_pred |> 
  select(country, year, mean, te) |> 
  distinct() |> 
  group_by(year, te) |> 
  summarize(mean = sum(mean), .groups = "keep") |> 
  prev_bar_plot(y_var = mean, te_var = te) +
  xlab("Year") +
  scale_x_continuous(limits = c(1990-0.5,2023+0.5), expand = c(0,0)) +
  theme(plot.margin = unit(c(t=0.1, r=0.05, b=0.05, l=0.1), "cm"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7))
prev1990_bar

## Case lineages assignments - donut plot----
prev1990_donut <- cases_prev_pred |> 
  group_by(te) |> 
  summarize(te_count = sum(mean)) |> 
  ungroup() |> 
  prev_donut_plot(te_var = te)
prev1990_donut

## Inset donut onto bar graph for panels A & B----
prev1990_plot <- ggdraw()+
  draw_plot(prev1990_bar + 
              theme(legend.position = "none"))+
  draw_plot(prev1990_donut, height=0.45,x=0.3,y=0.6) 
prev1990_plot

#-----------------------------------------------------------#
# Fig. S3 C&D - Post 2000 Transmission Units----
#-----------------------------------------------------------#

## Read data----
sampling <- "post_2000"

source(here("analysis","00_functions_settings.R"), local = T)

hmm_dist_plotting <- readRDS(str_c(clust_dir, "/xi_dist_plotting.rds"))
hmm_louvain_hclusters <- readRDS(str_c(clust_dir, "/hmm_louvain_hclusters.rds"))
hmm_louvain_subgrps <- readRDS(str_c(clust_dir, "/louvain_hmm_subgrps.rds"))
hmm_mapping_edges <- readRDS(str_c(clust_dir, "/hmm_mapping_edges.rds"))
hmm_consensus_clusters <- readRDS(str_c(clust_dir, "/hmm_consensus_clusters.rds"))

## Update mapping options----
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

## Prep for plotting----
hmm_filter_stat <- quantile(hmm_mapping_edges$log_xi_norm, filter_cutoff)
hmm_filter_stat_mean <- mean(hmm_mapping_edges$log_xi_norm)

hmm_filtered_edges <- hmm_mapping_edges |> 
  filter(log_xi_norm>hmm_filter_stat_mean)

hmm_dist_plotting <- hmm_dist_plotting |> 
  select(country1, country2, log_xi_norm) |> 
  mutate_at(c("country1", "country2"),  ~as.character(.))

## Consistent clustering
clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Mauritania", "Central African Rep."),
  cluster_4k_assign = factor(c(1:4))) |> 
  left_join(hmm_louvain_hclusters |> select(country, hmm_cluster = cluster_4k),
            by = "country") 

hmm_louvain_hclusters <- hmm_louvain_hclusters |> 
  left_join(clust_reassign |> select(hmm_cluster, cluster_4k_assign),
            by = c("cluster_4k"="hmm_cluster"))
consensus_clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Cameroon", "Mauritania"),
  cluster_hmm3k_assign = factor(c(1:3, 3)))  |> 
  left_join(hmm_consensus_clusters |> select(country1, hmm_cluster = cluster),
            by = c("country"="country1")) 

## Shapefiles for mapping----
hmm_clust_sf <- africa_sf |> 
  mutate(country = country_rename(country)) |> 
  left_join(hmm_louvain_hclusters,
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

## Fig S3C HMM Mapped edges----
hmm_map_post_2000 <- ggplot() + 
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
hmm_map_post_2000

## Fig S3D HMM distance heatmap----

heatmap_order_2000 <-  c("Gambia", "Guinea-Bissau", "Senegal", "Sierra Leone","Guinea", 
                         "Liberia", "Cote d'Ivoire", "Mauritania", "Mali", "Benin", 
                         "Burkina Faso", "Togo", "Ghana", "Nigeria", "Cameroon", 
                         "Chad", "Niger", "Equatorial Guinea", "Gabon", "Sao Tome", 
                         "Rep. of the Congo", "Central African Rep.", "DRC", 
                         "Sudan", "South Sudan", "Somalia",  "Eritrea", "Ethiopia",  
                         "Djibouti", "Rwanda", "Uganda", "Kenya", "Tanzania", "Zambia", 
                         "Burundi", "Angola", "Namibia", "Comoros", "Madagascar",  
                         "Malawi", "Botswana", "Mozambique", "Zimbabwe", "South Africa", 
                         "Eswatini", "Lesotho") 

hmm_dist_plotting_ordered <- countries_df |> 
  left_join(hmm_dist_plotting |>  
              mutate_at(c("country1", "country2"), 
                        ~country_rename(.))) |> 
  mutate(country1 = factor((country1), 
                           levels = heatmap_order_2000, ordered = T),
         country2 = factor((country2), 
                           levels = rev(heatmap_order_2000), ordered = T))
hmm_clust_box <- hmm_consensus_clusters |> 
  mutate(country1 = factor(country1, 
                           levels = heatmap_order_2000, ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(x = as.numeric(first(country1)), xend = as.numeric(last(country1))) |> 
  ungroup() |> 
  mutate(country1 = factor(country1, 
                           levels = rev(heatmap_order_2000), ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(y = as.numeric(first(country1)), yend = as.numeric(last(country1))) |> 
  ungroup() |> 
  select(cluster_assign, x, xend, y, yend) |> 
  distinct() |> 
  arrange(cluster_assign)

hmm_heatmap_post_2000 <- ggplot() +
  geom_tile(data = hmm_dist_plotting_ordered, 
            aes(y = country2, x = country1, fill = log_xi_norm, 
                alpha = log_xi_norm), 
            color = "transparent", lwd = 0.01, linetype = 1) +
  geom_rect(data = hmm_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                      ymin = y-0.5, ymax = yend+0.5, 
                                      color = cluster_assign), fill = NA, 
            lwd = 0.5) +
  fig3_heatmap_options +
  scale_color_manual(values = pub_cluster_palette)

hmm_heatmap_post_2000

#-----------------------------------------------------------#
# Fig. S6 B - Post 2000 Assigned prevalence----
#-----------------------------------------------------------#

## Read data----
cases_prev_pred <- readRDS(str_c(xi_dir, "/cases_prev_pred.rds")) |> 
  clean_df() |> 
  mutate(te = factor(te, levels = te_order, labels = pub_te_order))

## Bar graph----
prev2000_bar <-  cases_prev_pred |> 
  select(country, year, mean, te) |> 
  distinct() |> 
  group_by(year, te) |> 
  summarize(mean = sum(mean), .groups = "keep") |> 
  prev_bar_plot(y_var = mean, te_var = te) +
  xlab("Year") +
  scale_x_continuous(limits = c(1990-0.5,2023+0.5), expand = c(0,0)) +
  theme(plot.margin = unit(c(t=0.1, r=0.05, b=0.05, l=0.1), "cm"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7))
prev2000_bar

## Case lineages assignments - donut plot----
prev2000_donut <- cases_prev_pred |> 
  group_by(te) |> 
  summarize(te_count = sum(mean)) |> 
  ungroup() |> 
  prev_donut_plot(te_var = te)
prev2000_donut

## Inset donut onto bar graph for panels A & B----
prev2000_plot <- ggdraw()+
  draw_plot(prev2000_bar + theme(legend.position = "none"))+
  draw_plot(prev2000_donut, height=0.45,x=0.3,y=0.6) 
prev2000_plot

#-----------------------------------------------------------#
# Fig. S6 - Combine plots----
#-----------------------------------------------------------#

prev_plots_legend <- get_plot_component(prev1990_bar, "guide-box-bottom", 
                                              return_all = T)

prev_plots <- cowplot::plot_grid(prev1990_plot + theme(legend.position = "none"), 
                                 prev2000_plot + theme(legend.position = "none"),
                                       ncol = 2, labels = c("A.", "B."),
                                       label_fontfamily = "serif", 
                                       label_size = 9)
prev_plots <- cowplot::plot_grid(prev_plots, prev_plots_legend,
                                       rel_heights = c(1, 0.1), nrow = 2)
prev_plots

#-----------------------------------------------------------#
# Save figures----
#-----------------------------------------------------------#

ggsave("Supplemental_FigS1_louvain_heatmaps.pdf", plot = figS1_louvain_heatmaps, 
       width = 7.5, height = 4.5, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS1_louvain_heatmaps.png", plot = figS1_louvain_heatmaps, 
       width = 7.5, height = 4.5, dpi = 400,
       path = here::here("figures/manuscript_figures"))

ggsave("Supplemental_FigS2_rho_plot.png", plot = figS2_rho_plot, 
       width = 7.5, height = 8, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS2_rho_plot.pdf", plot = figS2_rho_plot, 
       width = 7.5, height = 9, dpi = 400,
       path = here::here("figures/manuscript_figures"))

ggsave("Supplemental_FigS3A_heatmap_post1990.pdf", plot = hmm_heatmap_post_1990 + theme(legend.position = "none"), 
       width = 2.4, height = 2.35, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS3B_map_post1990.pdf", plot = hmm_map_post_1990, 
       width = 2.7, height = 2.75, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS3C_heatmap_post2000.pdf", plot = hmm_heatmap_post_2000 + theme(legend.position = "none"), 
       width = 2.4, height = 2.35, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS3D_map_post2000.pdf", plot = hmm_map_post_2000, 
       width = 2.7, height = 2.75, dpi = 400,
       path = here::here("figures/manuscript_figures"))

ggsave("Supplemental_FigS6_prev_plots.pdf", plot = prev_plots, 
       width = 7.5, height = 3, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS6_prev_plots.png", plot = prev_plots, 
       width = 7.5, height = 3, dpi = 400,
       path = here::here("figures/manuscript_figures"))

