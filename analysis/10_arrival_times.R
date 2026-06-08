#-----------------------------------------------------------------#
# Project: Defining epidemiologically relevant transmission
#          units in Africa
# File: 10_arrival_times.R
# Purpose: Analyzing onward spread from outbreak simulations
#         this script also creates Figure S7 for the supplement
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

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

#-----------------------------------------------------------#
# Read data----
#-----------------------------------------------------------#

source(here("analysis","02_import_data.R"), local = T) 

cases_lineages_hmm <- readRDS(str_c(xi_dir, "/cases_lineages_assign.rds")) |> 
  mutate(country=country_rename(country))
downstream_pred <- readRDS(str_c(sim_dir, "/downstream_pred.rds"))
sim_downstream <- readRDS(str_c(sim_dir, "/sim_downstream.rds"))

#-----------------------------------------------------------#
# Standardize names----
#-----------------------------------------------------------#

#Country renames for consistency 
countries <- country_rename(countries_subset)
introductions_dist <- introductions_dist |> 
  mutate(country=country_rename(country))
arrival_times <- arrival_times |> 
  mutate(country=country_rename(country))

#-----------------------------------------------------------#
# Downstream arrival stats----
#-----------------------------------------------------------#

downstream_pred <- downstream_pred |> 
  group_by(seed_country, country) |> 
  summarize(median_arrival_time = median(arrival_time, na.rm = T),
            mean_arrival = mean(arrival_time, na.rm = T),
            sd_arrival = sd(arrival_time, na.rm = T),
            q25_arrival = quantile(arrival_time, 0.25, na.rm = T),
            q75_arrival = quantile(arrival_time, 0.75, na.rm = T)) |> 
  ungroup() |> 
  mutate_at(c("country", "seed_country"), 
            ~(country_rename(.))) |> 
  rename(arrival_time=median_arrival_time)

#-----------------------------------------------------------#
# Seed countries---
#-----------------------------------------------------------#

# Likely seed countries based on phylogeographic analyses/observed data
likely_seeds <- introductions_dist |> 
  filter_out(is.na(prob_dist)) |> 
  mutate(country = if_else(te=="T2" & country == "Liberia", "Uganda", country)) |> 
  group_by(te) |>
  mutate(keep = if_else(prob_dist==max(prob_dist) | prob_dist>0.1, 1, 0)) |>
  ungroup() |>
  mutate(keep = case_when(
    te=="T1" & country %in% c("Benin", "Burkina Faso", "Cameroon", "Cote d'Ivoire", 
                              "Ghana", "Guinea", "Liberia", "Mali", "Niger", 
                              "Nigeria", "Sierra Leone") ~ 1,
    te=="T2" & country %in% c("Uganda", "Nigeria") ~ 1,
    te=="T3" & country == "Ethiopia" ~ 1,
    te=="T4" & country %in% c("Ethiopia", "Somalia", "Djibouti") ~ 1,
    te=="T5" & country %in% c("Somalia", "DRC", "Rwanda") ~ 1,
    te=="T6" & country %in% c("Djibouti", "Somalia") ~ 1,
    te=="T7" & country == "Sierra Leone" ~ 1,
    te=="T8" & country %in% c("Mozambique", "South Africa") ~ 1,
    te %in% str_c("T", 9:17) ~ keep,
    T ~ 0
  )) |> 
  filter(keep == 1) |> 
  select(-keep) |> 
  group_by(te) |> 
  mutate(order = row_number()) |> 
  ungroup() |> 
  select(te, seed = country, order)


#-----------------------------------------------------------#
# Arrival times---
#-----------------------------------------------------------#

# Inferred arrival times from hmm 
arrival_times_hmm <- data.frame(country = rep(countries, times = length(te_order[te_order!="S"])),
                                te = rep(te_order[te_order!="S"], each = length(all_countries))) |> 
  left_join(cases_lineages_hmm |> 
              filter(!is.na(te),
                     source == "Inferred") |> 
              select(country, year, te) |> 
              distinct(), by = c("country", "te")) |> 
  mutate(year = replace_na(year, 2099)) |> 
  group_by(country, te) |> 
  mutate(arrival_year_hmm = min(year)) |> 
  ungroup() |> 
  select(country, te, arrival_year_hmm) |> 
  distinct() |> 
  arrange(te, country) |> 
  filter(arrival_year_hmm<2099)



#Time to arrival, all likely seeds
timings_df <- te_order[te_order!="S"] |> 
  set_names() |> 
  purrr::map(~{
    downstream_pred |> 
      filter(seed_country %in% (likely_seeds |> filter(te == .x))$seed,
             seed_country!=country) |> 
      left_join(arrival_times |> filter(te == .x) |> 
                  select(country, obs_time_to_arr, arrival_year, first_seen), 
                by = "country") |> 
      left_join(arrival_times_hmm |> filter(te == .x) |> 
                  select(country, arrival_year_hmm), 
                by = "country")}) |> 
  bind_rows(.id = "te") |> 
  left_join(introductions_dist |> select(te, year_med, year_min, year_max) |> distinct(), 
            by = "te") |> 
  rename(sim_arrival_time = arrival_time,
         arrival_time_obs = obs_time_to_arr) |> 
  mutate(arrival_year = if_else(arrival_year==2099, as.numeric(NA), arrival_year),
         obs_arrival_time_med = arrival_year-year_med,
         obs_arrival_time_min = arrival_year-year_min,
         obs_arrival_time_max = arrival_year-year_max,
         hmm_arrival_time_med = arrival_year_hmm-year_med,
         hmm_arrival_time_min = arrival_year_hmm-year_min,
         hmm_arrival_time_max = arrival_year_hmm-year_max) |> 
  select(-year_med, -year_min, -year_max, -arrival_year, -first_seen, -arrival_year_hmm,
         -arrival_time_obs) |> 
  pivot_longer(cols = contains("arrival_time_"),
               names_to = c("source", "measure"),
               names_pattern = "(...)\\_(.*)",
               values_to = "arrival_time"
  ) |> 
  left_join(likely_seeds, by = c("te", "seed_country"="seed")) |> 
  mutate(te = factor(te, ordered = T, levels = te_order, labels = pub_te_order)) 

#-----------------------------------------------------------#
# Prep arrival times for plotting---
#-----------------------------------------------------------#

timings_toplot <- timings_df |> 
  rename(Simulated = sim_arrival_time,
         Observed = arrival_time) |> 
  filter(measure %in% c("arrival_time_min", "arrival_time_max", "arrival_time_med")) |>
  mutate(measure = factor(measure, levels = c("arrival_time_min", "arrival_time_med", "arrival_time_max"),
                          labels = c("Earliest", "Median", "Latest")),
         source = factor(source, levels = c("obs", "hmm"), labels = c("Detected", "Inferred"))) |> 
  pivot_wider(id_cols = c("te", "seed_country", "country", "source", "Simulated",
                          "q25_arrival", "q75_arrival", "mean_arrival","sd_arrival"),
              names_from = measure, values_from = Observed) |> 
  # No observed/inferred T3 outside of Ethiopia -> remove
  filter(!is.na(Median),
         te!="AFR3") |> 
  mutate(te_seed = str_c(te, ", Seed: ", seed_country),
         te = factor(te, ordered = T, levels = pub_te_order, labels = pub_te_order)) |> 
  arrange(te, seed) |> 
  mutate(te_seed = factor(te_seed, levels = unique(te_seed), ordered = T)) |> 
  arrange(te, country)

## Assess seed performance by lineage/source----

timings_df |> 
  group_by(source) |> 
  filter_out(is.na(arrival_time)) |> 
  filter(measure=="arrival_time_med") |> 
  group_by(source, te, seed_country) |> 
  summarize(corr = cor(sim_arrival_time, arrival_time)) |> 
  filter(corr == max(corr) | is.na(corr)) |> 
  ungroup() |> 
  print(n = 36)
# hmm: T1 - Ghana, T5 - Rwanda, T6 - Djibouti, T7 - Sierra Leone, T8 - Mozambique, T9 - Guinea, T10 - Kenya, T11 - Mozambique, T12 - Nigeria, T13 - Kenya
# obs: T1 - Ghana, T5 - Somalia, T6 - Djibouti, T7 - Sierra Leone, T8 - , T9 - Guinea, T10 - Kenya, T11 - Mozambique, T12 - Cameroon, T13 - Kenya

timings_df |> 
  filter_out(is.na(arrival_time)) |> 
  filter(measure=="arrival_time_med") |> 
  group_by(te, seed_country) |> 
  summarize(corr = cor(sim_arrival_time, arrival_time)) |> 
  filter(corr == max(corr) | is.na(corr)) |> 
  ungroup() |> 
  print(n = 36)
# comb: T1 - Ghana, T5 - Rwanda, T6 - Djibouti, T7 - Sierra Leone, T8 - South Africa, T9 - Guinea, T10 - Kenya, T11 - Mozambique, T12 - Nigeria, T13 - Kenya

## Assign optimal seed----
optimal_seed <- likely_seeds |> 
  mutate(te = factor(te, ordered = T, levels = te_order, labels = pub_te_order)) |> 
  mutate(best = case_when(
    te %in% c(str_c("AFR",c(7,9:11,13:17))) ~ 1,
    te == "AFR12" & seed == "Nigeria" ~ 1,
    te == "AFR1" & seed == "Ghana" ~ 1,
    te == "AFR2" & seed == "Uganda" ~ 1,
    te == "AFR4" & seed == "Somalia" ~ 1,
    te == "AFR5" & seed == "Rwanda" ~ 1,
    te == "AFR6" & seed == "Djibouti" ~ 1,
    te == "AFR8" & seed == "Mozambique" ~ 1,
    T ~ 0))

final_timings_toplot <- timings_toplot |> 
  left_join(optimal_seed |> select(-order), by = c("te", "seed_country"="seed")) |> 
  filter(best==1) |> 
  arrange(te, seed, country) |> 
  select(-best)

#-----------------------------------------------------------#
# Plot simulated vs observed arrival times---
#-----------------------------------------------------------#

# by lineage (not in manuscript)
final_timings_facet_plot <- final_timings_toplot |> 
  ggplot(aes(y = Simulated, x = Median, color = te, shape = source, alpha = source)) +
  geom_abline(color = "grey70", alpha = 0.4) +
  geom_pointrange(aes(xmin = Earliest, xmax = Latest, color = te), size = 0.25,
                  position = position_jitter(width = 0.5, height = 0.5)) +
  theme_bw() +
  facet_wrap(~te_seed, ncol = 5) +
  scale_color_manual(values = pub_te_palette) +
  scale_shape_manual(values = c("Detected"=16, "Inferred" = 1)) +
  scale_alpha_manual(values = c("Detected"=0.9, "Inferred" = 0.4)) +
  scale_x_continuous(limits = c(-1,30), breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = c(-1,30), breaks = scales::pretty_breaks()) +
  labs(color = "Lineage",
       shape = "Source") +
  ylab("Simulated Arrival") +
  xlab("Observed Arrival") +
  guides(alpha = "none") +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.35, "cm"),
        legend.key.spacing.x = unit(0.1, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 8, family = "serif")) +
  coord_fixed()

## Figure S7 for supplement
fig_S7_timings <- final_timings_toplot |> 
  ggplot(aes(y = Simulated, x = Median, color = te)) +
  geom_abline(color = "grey70", alpha = 0.4) +
  geom_point(shape = 1, position = position_jitter(width = 0.25, height = 0.5,
                                                   seed = 2),
             size = 2.5) +
  geom_linerange(aes(xmin = Earliest, xmax = Latest, color = te),
                 alpha = 0.5, position = position_jitter(width = 0.25, height = 0.5,
                                                         seed = 2)) +
  geom_linerange(aes(ymin = q25_arrival, ymax = q75_arrival, color = te),
                 alpha = 0.3, linetype = "dotted", position = position_jitter(width = 0.25, height = 0.5,
                                                                              seed = 2)) +
  theme_bw() +
  facet_wrap(~source, ncol = 2) +
  scale_color_manual(values = pub_te_palette) +
  scale_alpha_manual(values = c("Detected"=0.9, "Inferred" = 0.4)) +
  scale_x_continuous(limits = c(-1,31), breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = c(-1,31), breaks = scales::pretty_breaks()) +
  labs(color = "Lineage") +
  ylab("Simulated Arrival") +
  xlab("Observed Arrival") +
  guides(color = guide_legend(nrow = 2, title.position = "left"),
         shape = "none",
         alpha = "none") +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_fixed()
fig_S7_timings


#-----------------------------------------------------------#
# Correlations---
#-----------------------------------------------------------#

timing_obs_corr <- final_timings_toplot |> 
  group_by(source) |> 
  mutate(corr = cor(Simulated, Median, method = "pearson"),
         corr_earliest = cor(Simulated, Earliest, method = "pearson"),
         corr_latest = cor(Simulated, Latest, method = "pearson")) |> 
  ungroup() |> 
  select(source, starts_with(c("corr"))) |> 
  distinct()
timing_obs_corr

te_timing_obs_corr <- final_timings_toplot |> 
  group_by(source, te) |> 
  mutate(corr_te = cor(Simulated, Median, method = "pearson"),
         corr_earliest_te = cor(Simulated, Earliest, method = "pearson"),
         corr_latest_te = cor(Simulated, Latest, method = "pearson")) |> 
  ungroup() |> 
  select(te, te_seed, source, starts_with(c("corr"))) |> 
  distinct() |> 
  arrange(source, te) #|> filter_out(is.na(corr_te))
print(te_timing_obs_corr, n = nrow(te_timing_obs_corr))


#-----------------------------------------------------------#
# Save figures---
#-----------------------------------------------------------#

ggsave("figure_S7_arrivals.pdf", plot = fig_S7_timings, width = 7, height = 6, dpi = 400,
       path = here::here("figures/manuscript_figures"))
ggsave("figure_S7_arrivals.png", plot = fig_S7_timings, width = 7, height = 6, dpi = 400,
       path = here::here("figures/manuscript_figures"))


