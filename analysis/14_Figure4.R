#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 14_Figure4.R
# Purpose: This script creates Figure 4 showing outbreak risks
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

options(scipen = 9999)

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

set.seed(seed)

num_sims <- 2000

#-----------------------------------------------------------#
# Read data----
#-----------------------------------------------------------#

source(here("analysis","02_import_data.R"), local = T) 

secondary_pred <- readRDS(str_c(sim_dir, "/secondary_pred.rds"))
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

## Outbreak risks----
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
  bind_rows(.id = "seed_country") |> 
  mutate_at(c("seed_country", "country"), ~country_rename(., sudan_rename = T))

## Countries with highest risks----
secondary_high_risk <- secondary_risk |> 
  group_by(seed_country) |> 
  mutate(risk = if_else(seed_country==country, 0, pct),
         max_risk = max(risk)) |> 
  ungroup() |> 
  distinct(seed_country, max_risk) |> 
  filter(max_risk>0.21) |> 
  mutate(seed_country = country_rename(seed_country, sudan_rename = T))

## Shapefiles----
secondary_sfs <- country_rename(countries, sudan_rename = T) |> 
  set_names() |> 
  purrr::map(~{
    africa_sf |> 
      mutate(country = country_rename(country, sudan_rename = T)) |> 
      left_join(secondary_risk |> 
                  filter(seed_country==.x),
                by = "country")}) 

## Country display order----
country_order <-  heatmap_order_comb

#Country renames
# change_from <- "Sudan"
# change_to <- "Sudan & South Sudan"
# for(i in 1:length(change_from)){
#   names(secondary_sfs)[names(secondary_sfs)==change_from[i]] <- change_to[i]
#   secondary_high_risk$seed_country[secondary_high_risk$seed_country==change_from[i]] <- change_to[i]
#   secondary_risk$seed_country[secondary_risk$seed_country==change_from[i]] <- change_to[i]
# }

## Edges----
hmm_edges <- edges |> 
  mutate_at(c("country1", "country2"), ~country_rename(.))  |> 
  left_join(hmm_dist_plotting, by = c("country1", "country2")) |>
  arrange(log_xi_norm) |> 
  mutate_at(c("country1", "country2"), ~country_rename(., sudan_rename = T)) |> 
  mutate(country1 = factor(country1, levels = country_order, ordered = T),
         country2 = factor(country2, levels = country_order, ordered = T))

## Consensus clusters----
# Ensure correct ordering
consensus_clust_reassign <- tibble(
  country = c("South Africa", "Kenya", "Cameroon"),
  cluster_hmm3k_assign = factor(1:3))  |> 
  left_join(hmm_consensus_clusters |> select(country1, hmm_cluster = cluster),
            by = c("country"="country1")) 
hmm_consensus_clusters <- hmm_consensus_clusters |> 
  left_join(consensus_clust_reassign |> select(contains("hmm")) |> distinct(),
            by = c("cluster" = "hmm_cluster")) |> 
  rename(cluster_assign = cluster_hmm3k_assign) |> 
  mutate(cluster_assign = factor(cluster_assign),
         country1 = country_rename(country1, sudan_rename = T))

# Object for drawing boxes around clusters in heat map
hmm_clust_box <- hmm_consensus_clusters |> 
  mutate(country1 = factor(country1, 
                           levels = country_order, ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(x = as.numeric(first(country1)), xend = as.numeric(last(country1))) |> 
  ungroup() |> 
  mutate(country1 = factor(country1, 
                           levels = rev(country_order), ordered = T)) |> 
  arrange(country1) |> 
  group_by(cluster_assign) |> 
  mutate(y = as.numeric(first(country1)), yend = as.numeric(last(country1))) |> 
  ungroup() |> 
  select(cluster_assign, x, xend, y, yend) |> 
  distinct() |> 
  arrange(cluster_assign)

#-----------------------------------------------------------#
# Figure 4 maps ----
#-----------------------------------------------------------#

#map over countries
figure4_maplist <- secondary_risk_map(secondary_high_risk)

#Get legend
figure4_legend <- get_plot_component(figure4_maplist[[1]], "guide-box-bottom", 
                                      return_all = T)
#Remove legends from individual maps
figure4_maplist <- purrr::map(figure4_maplist, ~{.x + theme(legend.position = "none")})

#-----------------------------------------------------------#
# Figure 4 heatmap ----
#-----------------------------------------------------------#

## Plot options----
heatmap_options <- list(scale_fill_gradient2(high = "#660000", mid = "#B21807", 
                                             low = "white", midpoint = 0.5, limits = c(0,1)),
                        theme_classic(),
                        theme(text = element_text(size = 6,
                                                  family = "serif"),
                              axis.title.y = element_text(size = 7),
                              axis.text.x.bottom = element_text(vjust = 0.5),
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



## Assign clusters---- 
# from hmm_consensus_clusters
secondary_risk_clust <- secondary_risk  |> 
  filter(seed_country!=country) |> 
  mutate_at(c("seed_country", "country"), ~country_rename(., sudan_rename = T)) |> 
  mutate(outbreaks = pct*num_sims,
         nonoutbreaks = (1-pct)*num_sims) |> 
  right_join(hmm_consensus_clusters |> select(country1, "country1_clust"=cluster_assign) |> 
               mutate(country1 = as.character(country1),
                      country1_clust = factor(country1_clust,
                                              levels = c(1:3),
                                              labels = c("hmmS", "hmmCE", "hmmW"))), 
             by = c("seed_country"="country1"))  |> 
  right_join(hmm_consensus_clusters |> select(country1, "country2_clust"=cluster_assign) |> 
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

## Outbreaks by seed----
# Side bar plot
outbreaks_by_seed <- purrr::imap(names(sim_secondary), ~{
  c <- .x
  sim_secondary[[.y]] |> 
    purrr::map(~{.x |> 
        filter(time==2, country!=c) |> 
        select(country, Presence)}) |> 
    bind_rows(.id = "run") |> 
    mutate(seed_country=c) |> 
    select(seed_country, everything()) |> 
    mutate_at(c("seed_country", "country"), ~country_rename(., sudan_rename = T)) |> 
    left_join(secondary_risk_clust |> select(seed_country, country, country1_clust, 
                                             country2_clust, same_cluster),
              by = c("seed_country", "country")) |> 
    group_by(country) |> 
    mutate(country_risk = sum(Presence)/num_sims) |> 
    ungroup() |> 
    select(seed_country, country, same_cluster, country2_clust, country_risk) |> 
    distinct() |> 
    #Any outbreaks 
    mutate(outbreaks = sum(country_risk)) |> 
    #Same cluster
    group_by(same_cluster) |> 
    mutate(outbreaks_same = sum(country_risk)) |> 
    #By cluster
    group_by(country2_clust) |> 
    mutate(outbreaks_byclust = sum(country_risk)) |> 
    ungroup() |> 
    select(seed_country, same_cluster, country2_clust, starts_with("outbreaks")) |> 
    distinct()
}) |> 
  bind_rows() 

## Create heatmap----
figure4_risk_heatmap <- tibble(seed_country = rep(country_order, 
                                                   times = length(country_order)),
                                country = rep(country_order, 
                                              each = length(country_order))) |> 
  filter(seed_country!=country) |> 
  left_join(secondary_pred |>  
              mutate_at(c("seed_country", "country"), 
                        ~country_rename(., sudan_rename = T))) |> 
  mutate(seed_country = factor(seed_country, 
                               levels = rev(country_order), ordered = T),
         country = factor(country, 
                          levels = (country_order), ordered = T)) |> 
  ggplot() +
  geom_tile(aes(y = seed_country, x = country, fill = pct), 
            color = "white", lwd = 0.05, linetype = 1) +
  #Outline clusters
  geom_rect(data = hmm_clust_box, aes(xmin = x-0.5, xmax = xend+0.5, 
                                      ymin = y-0.5, ymax = yend+0.5, 
                                      color = cluster_assign), fill = NA, 
            lwd = 0.35) +
  #theme_classic() +
  heatmap_options +
  ggside::geom_ysidebar(data = outbreaks_by_seed |> 
                          mutate(seed_country = if_else(seed_country=="Sudan", 
                                                        "Sudan & South Sudan", 
                                                        seed_country)) |> 
                          mutate(seed_country = factor(seed_country, 
                                                       levels = rev(country_order), 
                                                       ordered = T),
                                 same_cluster = factor(same_cluster,
                                                       levels = c(0,1),
                                                       labels = c("Different Transmission Unit",
                                                                  "Same Transmission Unit"))) |> 
                          arrange(seed_country) |> 
                          select(seed_country, outbreaks_same, same_cluster) |> 
                          distinct(),
                        aes(x = outbreaks_same, y = seed_country, 
                            alpha = same_cluster), stat = "identity") + 
  ylab("Country with an outbreak") +
  ggside::scale_ysidex_continuous(position = "top", guide = guide_axis("Outbreaks"),
                                  breaks = c(1:3), n.breaks = 3) +
  scale_color_manual(values = pub_cluster_palette) +
  scale_alpha_manual(values = c("Same Transmission Unit"=0.4, 
                                "Different Transmission Unit"=1)) +
  theme(ggside.axis.text.x = element_text(angle = 0))
figure4_risk_heatmap


#-----------------------------------------------------------#
# Combine plots ----
#-----------------------------------------------------------#

figure4_grid1 <- cowplot::plot_grid(plotlist = figure4_maplist[1:10], ncol = 5)
figure4_grid2 <- plot_grid(plot_grid(plotlist = figure4_maplist[11:14], ncol = 2),
                           figure4_risk_heatmap + theme(legend.position = "none"), 
                           nrow = 1, rel_widths = c(0.7,1))
figure_4 <- plot_grid(plot_grid(figure4_grid1, figure4_grid2,
                          nrow = 2, rel_heights = c(1,1)),
                          figure4_legend, rel_heights = c(1,0.08), ncol = 1)
figure_4

#-----------------------------------------------------------#
# Supplemental Figure S4 maps ----
#-----------------------------------------------------------#

supp_fig_maplist <- secondary_risk_map(secondary_risk)
supp_fig_grid <- cowplot::plot_grid(plotlist = supp_fig_maplist, ncol = 6)
supp_fig <- plot_grid(supp_fig_grid,
                      figure4_legend, rel_heights = c(1,0.08), ncol = 1)
supp_fig

#-----------------------------------------------------------#
# Table S2. Outbreaks Table ----
#-----------------------------------------------------------#

outbreaks_by_seed_export <- outbreaks_by_seed |> 
  pivot_wider(id_cols = c(seed_country, outbreaks, 
                          country2_clust, outbreaks_byclust), 
              names_from = same_cluster, values_from = outbreaks_same) |> 
  set_names(c("seed_country", "outbreaks", "country2_clust", 
              "outbreaks_byclust", "outbreaks_diff_clust", "outbreaks_same_clust")) |> 
  pivot_wider(id_cols = c(seed_country, outbreaks, 
                          outbreaks_diff_clust, outbreaks_same_clust), 
              names_from = country2_clust, 
              values_from = outbreaks_byclust, names_prefix = "outbreaks_") |> 
  group_by(seed_country) |> 
  fill(starts_with("outbreak"), .direction = "downup") |> 
  ungroup() |> 
  distinct() |> 
  mutate_if(is.numeric, ~round(.,2)) |> 
  set_names(c("Country with an Outbreak", "Expected Outbreaks within 1 Year",
              "Outbreaks in Different Cluster", "Outbreaks in Same Cluster",
              "hmmW", "hmmS", "hmmCE"))

#-----------------------------------------------------------#
# Save plots & table----
#-----------------------------------------------------------#

ggsave("Figure_4.pdf", plot = figure_4, height = 8, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Figure_4.png", plot = figure_4, height = 8, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))

ggsave("Supplemental_FigS4.pdf", plot = supp_fig, height = 9, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("Supplemental_FigS4.png", plot = supp_fig, height = 9, width = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"))

write_csv(outbreaks_by_seed_export, here::here("figures/manuscript_figures","TableS2_outbreak_risk.csv"))

#-----------------------------------------------------------#
# Statistics for manuscript text----
#-----------------------------------------------------------#

summarize_risk <- function(data){
  data |> 
    summarize(tot_outbreaks = sum(outbreaks),
              tot_nonoutbreaks = sum(nonoutbreaks),
              tot_potential = tot_outbreaks+tot_nonoutbreaks,
              risk = tot_outbreaks/tot_potential) 
}

### Between vs within ----
between_vs_in_risk <- secondary_risk_clust |> 
  group_by(same_cluster) |> 
  mutate(mean_risk = mean(pct)) |> 
  group_by(cluster) |> 
  mutate(mean_risk_cluster = mean(pct)) |> 
  ungroup() |> 
  distinct(same_cluster, cluster, mean_risk, mean_risk_cluster) 
between_vs_in_risk |> distinct(same_cluster, mean_risk) 

within_vs_between <- secondary_risk_clust |> 
  select(same_cluster, outbreaks, nonoutbreaks) |> 
  group_by(same_cluster) |> 
  summarize_risk()
wb_rr <- within_vs_between$risk[2]/within_vs_between$risk[1]
wb_se <- sqrt(((1/within_vs_between$tot_outbreaks[2])+(1/within_vs_between$tot_outbreaks[1])) -
                      ((1/within_vs_between$tot_potential[2]) + (1/within_vs_between$tot_potential[1])))
wb_ci_low <- exp(log(wb_rr)-1.96*wb_se)
wb_ci_hi <- exp(log(wb_rr)+1.96*wb_se)
print(glue::glue("Within vs between RR = {round(wb_rr,1)} 95% CI: {round(wb_ci_low,1)}, {round(wb_ci_hi,1)}"))

### Between by cluster ----
clust_specific <- secondary_risk_clust |> 
  select(country1_clust, country2_clust, outbreaks, nonoutbreaks) |> 
  group_by(country1_clust, country2_clust) |> 
  summarize_risk()
clust_specific_rr <- clust_specific$risk[4]/clust_specific$risk[6]
clust_specific_se <- sqrt(((1/clust_specific$tot_outbreaks[4])+(1/clust_specific$tot_outbreaks[6])) -
                ((1/clust_specific$tot_potential[4]) + (1/clust_specific$tot_potential[6])))
clust_specific_ci_low <- exp(log(clust_specific_rr)-1.96*clust_specific_se)
clust_specific_ci_hi <- exp(log(clust_specific_rr)+1.96*clust_specific_se)
print(glue::glue("hmmCE to hmmS vs to hmmW RR = {round(clust_specific_rr,1)} 95% CI: {round(clust_specific_ci_low,1)}, {round(clust_specific_ci_hi,1)}"))


cluster_risk <- secondary_risk_clust |> 
  select(same_cluster, country1_clust, country2_clust, pct) |> 
  group_by(country1_clust, country2_clust) |> 
  mutate(mean_risk = mean(pct)) |> 
  ungroup() |> 
  select(-pct) |> 
  distinct()
cluster_risk

