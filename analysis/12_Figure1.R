#------------- Figure 1 -------------####
#### Observed lineages & cases & attributable prevalence

### Initialize----
library(tidyverse)
library(here)
library(plotly)
library(cowplot)

options(scipen = 9999)

# Source functions
source(here("analysis","00_functions_settings.R"), local = T)

set.seed(seed)

### 1A: Observed data----
#### Read data----
source(here("analysis","02_import_data.R"), local = T) 

#### Prep for plotting----
country_order <- sort(unique(country_rename(cases_lineages$country)), decreasing = T)

cases_lineages <- cases_lineages |> 
  select(-subregion) |> 
  filter(!is.na(country)) |>
  mutate(country = factor(country_rename(country), levels = country_order, ordered = T))

lineages_toplot <- cases_lineages |> 
  select(country, year, te) |> 
  filter(!is.na(te)) |> 
  mutate(te = fct_drop(te)) |>
  arrange(country, year, te) |> 
  distinct() |> 
  mutate(te = factor(te, levels = te_order, labels = pub_te_order))

cases_toplot <- cases_lineages |>
  filter(cases>0) |> 
  select(country, year, cases) |> 
  distinct()

cases_assign <- cases_lineages |> 
  select(country, year, te, cases) |> 
  distinct() |> 
  mutate(assign = if_else(!is.na(te) & !str_detect(te, "S"), 1, 0)) |> 
  group_by(country, year) |> 
  mutate(n_lineages = sum(assign)) |> 
  ungroup() |> 
  mutate(te_count = if_else(n_lineages>1, (cases/n_lineages)*assign, cases),
         te_assign = if_else(assign == 1, te, "Unsequenced"),
         te_assign = factor(te_assign, 
                            levels = c("Unsequenced", te_order[-length(te_order)]),
                            labels = c("Unsequenced", pub_te_order[-length(pub_te_order)]))) |> 
  group_by(year, te_assign) |> 
  mutate(te_count = sum(te_count)) |> 
  select(year, te_assign, te_count) |> 
  ungroup() |> 
  distinct() |> 
  arrange(year, te_assign)

#### Figure 1A top: Assigned prevalence based on sequences----

## Figure 1A top larger panel - bar graph of assigned prevalence
figure1a_top_bar <- cases_assign |> 
  prev_bar_plot(y_var = te_count, te_var = te_assign) +
  scale_x_continuous(limits = c(1970-0.5,2023+0.5), expand = c(0,0))
figure1a_top_bar

##Figure 1A top inset: case lineages assignments - donut plot
figure1a_top_donut <- cases_assign |> 
  group_by(te_assign) |> 
  summarize(te_count = sum(te_count)) |> 
  ungroup() |> 
  prev_donut_plot() 
figure1a_top_donut

#### Figure 1A bottom: observed cases by country with observed lineages----

cases_country_assign <- cases_lineages |> 
  select(country, year, te, cases) |> 
  distinct() |> 
  mutate(assign = if_else(!is.na(te), 1, 0)) |> 
  group_by(country, year) |> 
  mutate(n_lineages = sum(assign)) |> 
  ungroup() |> 
  mutate(te_count = if_else(n_lineages>1, (cases/n_lineages)*assign, cases),
         te_assign = if_else(assign == 1, te, "Unsequenced"),
         te_assign = factor(te_assign, 
                            levels = c(te_order, "Unsequenced"),
                            labels = c(pub_te_order, "Unsequenced"))) |> 
  distinct() |> 
  arrange(country, year, te_assign) |> 
  group_by(country, year) |> 
  mutate(te_num = row_number()) |> 
  ungroup() |> 
  filter(cases>0 | n_lineages>0) |> 
  mutate(cases = if_else(cases==0, n_lineages, cases))


figure_1a_bottom <- ggplot(data = cases_country_assign, 
                           aes(x = year, y = country, color = te_assign)) +
  geom_point(color = "transparent", show.legend = F) +
  geom_point(data = cases_country_assign |> filter(te_num==1 & te_assign=="Unsequenced"),
             aes(size = cases), 
             alpha = 1, show.legend = T, shape = 1) +
  geom_point(data = cases_country_assign |> filter(te_num==1 & te_assign!="Unsequenced"),
             aes(size = cases), 
             alpha = 1, show.legend = T) +
  geom_point(data = cases_country_assign |> filter(te_num==2),
             aes(size = cases/4), 
             alpha = 1, show.legend = F, shape = 1, stroke = 0.5) +
  geom_point(data = cases_country_assign |> filter(te_num==3),
             aes(size = cases/8), 
             alpha = 1, show.legend = F,  shape = 1, stroke = 0.25) +
  geom_point(data = cases_country_assign |> filter(te_num==4),
             aes(size = cases/16), 
             alpha = 1, show.legend = F, shape = 16, stroke = 0) +
  theme_minimal() +
  xlab("Year") +
  scale_color_manual(values = c(pub_te_palette, "Unsequenced"="grey70"), drop = F) +
  scale_size_continuous(range = c(0.01,2.5), breaks = c(100,1000,10000,100000),
                        limits = c(0,max(cases_country_assign$cases)+1)) +
  scale_x_continuous(limits = c(1970-0.5,2023+0.75), expand = c(0,0)) +
  labs(size = "Reported Suspected\nCholera Cases", color = "Lineage") +
  case_lineage_plot_theme +
  guides(size = guide_legend(order=2, nrow = 1),
         fill = "none",
         color = guide_legend(nrow = 1))
figure_1a_bottom


### 1B: Inferred data----
#### Read data----

countries <- all_countries
M <- length(countries) 
n_time <- length(unique(cases_lineages$year))
spatial_objs <- create_connections(M, n_time)

cases_lineages_hmm <- readRDS(str_c(xi_dir, "/cases_lineages_assign.rds"))
cases_prev_pred <- readRDS(str_c(xi_dir, "/cases_prev_pred.rds")) |> 
  clean_df() |> 
  mutate(te = factor(te, levels = te_order, labels = pub_te_order))
gen_cases_prev <- readRDS(str_c(xi_dir, "/gen_pred_cases.rds")) |> 
  translate_cols() |> 
  clean_df() |> 
  mutate(te = factor(te, levels = te_order, labels = pub_te_order))

#### Prep for plotting----

cases_lineages_hmm2 <- cases_lineages_hmm |> 
  select(country, year, te, source) |> 
  distinct() |> 
  mutate(country = factor(country_rename(country), levels = country_order, ordered = T),
         source = factor(source),
         present = if_else(!is.na(source),1,0)) |> 
  filter(!is.na(te)) |> 
  mutate(source_num = if_else(source=="Observed",2,1)) |> 
  group_by(country, year) |> 
  mutate(sampled = if_else(max(source_num)==2,1,0)) |> 
  ungroup() |> 
  filter(source_num==1) |> 
  select(-source, -source_num, -present) |> 
  left_join(cases_country_assign |> 
              select(country, year, cases) |> 
              distinct(), by = c("country", "year")) |> 
  mutate(persistence = if_else(is.na(cases) | cases==0, 1, 0))

cases_filled <- cases_country_assign |> 
  bind_rows(cases_lineages_hmm2 |> 
              mutate(inferred=1) |> 
              filter(sampled==0) |> 
              select(-sampled)) |> 
  arrange(country, year, inferred, te_assign, te) |> 
  group_by(country, year) |> 
  fill(inferred) |> 
  ungroup() |> 
  filter(!is.na(te_assign)) |> 
  filter(te_assign!="Unsequenced" | (te_assign=="Unsequenced" & is.na(inferred))) |> 
  bind_rows(cases_lineages_hmm2 |> 
              mutate(te_assign = factor(as.character(te), 
                                 levels = c("Unsequenced", te_order[-length(te_order)]),
                                 labels = c("Unsequenced", pub_te_order[-length(pub_te_order)])),
                     inferred = 1) |> 
              select(-sampled)) |> 
  arrange(country, year, te_assign, inferred) |> 
  group_by(country, year, te_assign) |> 
  fill(inferred, te_num, persistence, .direction = "downup") |> 
  select(country, year, cases, te_assign, te_num, persistence, inferred) |> 
  distinct() |> 
  arrange(country, year, te_num, te_assign) |> 
  group_by(country, year) |> 
  mutate(te_num = if_else(is.na(te_num) & n()>1, row_number(), te_num),
         te_num = if_else(is.na(te_num) & n()==1, 1, te_num),
         persistence = if_else(is.na(persistence) | persistence==0, "Cases observed", "Inferred persistence"),
         cases = if_else(persistence==1 & is.na(cases),0,cases)) |> 
  ungroup() 


#### Figure 1B top: Assigned prevalence based on sequences----

## Figure 1B top larger panel: bar graph of assigned prevalence
figure1b_top_bar <- cases_prev_pred |> 
  select(country, year, mean, te) |> 
  distinct() |> 
  group_by(year, te) |> 
  summarize(mean = sum(mean), .groups = "keep") |> 
  prev_bar_plot(y_var = mean, te_var = te) +
  scale_x_continuous(limits = c(1970-0.5,2023+0.5), expand = c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(t=0.1, r=0.05, b=0.05, l=0.4), "cm"))
figure1b_top_bar

##Figure 1B top inset: case lineages assignments - donut plot
figure1b_top_donut <- cases_prev_pred |> 
  group_by(te) |> 
  summarize(te_count = sum(mean)) |> 
  ungroup() |> 
  prev_donut_plot(te_var = te)
figure1b_top_donut

#### Figure 1B bottom: observed cases by country with observed & inferred lineages----

figure_1b_bottom <- ggplot(data = cases_filled, 
                           aes(x = year, y = country, color = te_assign)) +
  geom_point(color = "transparent", show.legend = F) +
  geom_point(data = cases_filled |> filter(te_num==1 & te_assign=="Unsequenced" & 
                                             persistence=="Cases observed"),
             aes(size = cases,
                 shape = persistence), stroke = 0.5,  
             alpha = 1, show.legend = T) +
  geom_point(data = cases_filled |> filter(te_num==1 & te_assign!="Unsequenced" & 
                                             persistence=="Cases observed"),
             aes(size = cases),
             alpha = 1, show.legend = F) +
  geom_point(data = cases_filled |> filter(te_num==2 & persistence=="Cases observed"),
             aes(size = cases/4), 
             alpha = 1, show.legend = F, shape = 1, stroke = 0.5) +
  geom_point(data = cases_filled |> filter(te_num==3 & persistence=="Cases observed"),
             aes(size = cases/8), 
             alpha = 1, show.legend = F,  shape = 1, stroke = 0.25) +
  geom_point(data = cases_filled |> filter(te_num==4 & persistence=="Cases observed"),
             aes(size = cases/16), 
             alpha = 1, show.legend = F, shape = 16, stroke = 0) +
  geom_point(data = cases_filled |> filter(te_num==5 & persistence=="Cases observed"),
             aes(size = cases/32), 
             alpha = 1, show.legend = F, shape = 16, stroke = 0) +
  geom_point(data = cases_filled |> filter(persistence=="Inferred persistence",
                                           te_num==1),
             aes(shape = persistence), 
             size = 1, 
             alpha = 0.9, show.legend = T) +
  geom_point(data = cases_filled |> filter(persistence=="Inferred persistence",
                                           te_num==2),
             aes(shape = persistence), 
             size = 0.5, 
             alpha = 0.9, show.legend = F) +
  geom_point(data = cases_filled |> filter(persistence=="Inferred persistence",
                                           te_num==3),
             aes(shape = persistence), 
             size = 0.15, 
             alpha = 0.9, show.legend = F) +
  theme_minimal() +
  xlab("Year") +
  scale_color_manual(values = c(pub_te_palette, "Unsequenced"="grey70"), drop = F) +
  scale_x_continuous(limits = c(1970-0.5,2023+0.75), expand = c(0,0)) +
  scale_shape_manual(values = c("Inferred persistence" = 3, "Cases observed"=1)) +
  scale_size_continuous(range = c(0.01,2.5), breaks = c(100,1000,10000,100000)) +
  labs(size = "Reported Suspected Cholera Cases", color = "Lineage",
       shape = "Case observation") +
  case_lineage_plot_theme +
  guides(size = guide_legend(order=2, nrow = 1),
         fill = "none",
         shape = guide_legend(nrow = 1),
         color = guide_legend(nrow = 2, byrow = T)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(t=0.1, r=0.05, b=0.05, l=0.4), "cm"))
figure_1b_bottom

### Arrange plot components----

#inset donut onto bar graph for panels A & B
figure_1a_top <- ggdraw()+
  draw_plot(figure1a_top_bar + theme(legend.position = "none"))+
  draw_plot(figure1a_top_donut, height=0.92,x=-0.21,y=0.18) +
  theme(plot.margin = margin(l = 1.03, unit = "cm"))
figure_1a_top

figure_1b_top <- ggdraw()+
  draw_plot(figure1b_top_bar + theme(legend.position = "none"))+
  draw_plot(figure1b_top_donut, height=0.92,x=-0.26,y=0.18) 
figure_1b_top

#Combine prevalence & case/te plots
figure_1a <- plot_grid(figure_1a_top, figure_1a_bottom + theme(legend.position = "none"),
                       rel_heights = c(0.2,1), ncol = 1)
figure_1a

figure_1b <- plot_grid(figure_1b_top, figure_1b_bottom + theme(legend.position = "none"),
                       rel_heights = c(0.2,1), ncol = 1)
figure_1b

figure_1_legend <- get_plot_component(figure_1b_bottom, "guide-box-bottom", 
                                      return_all = T)

figure_1 <- plot_grid(figure_1a, figure_1b, ncol = 2, rel_widths = c(1.074,0.926),
                      labels = c("A.", "B. "), label_size = 8,
                      label_fontfamily = "serif")
figure_1 <- plot_grid(figure_1, figure_1_legend, nrow = 2,
                      rel_heights = c(1,0.06))

figure_1

## Save figures----
ggsave("Figure_1a.pdf", plot = figure_1a, width = 3.75, height = 6, dpi = 300,
       path = here::here("figures/manuscript_figures"), bg = "white")
ggsave("Figure_1b.pdf", plot = figure_1b, width = 3.75, height = 6, dpi = 300,
       path = here::here("figures/manuscript_figures"), bg = "white")
ggsave("Figure_1.pdf", plot = figure_1, width = 7.5, height = 7, dpi = 300,
       path = here::here("figures/manuscript_figures"), bg = "white")
ggsave("Figure_1.png", plot = figure_1, width = 7.5, height = 7, dpi = 400,
       path = here::here("figures/manuscript_figures"), bg = "white")

ggsave("Figure_1a_top.pdf", plot = figure_1a_top, width = 8, height = 3, dpi = 300,
       path = here::here("figures/manuscript_figures"), bg = "white", scale = 0.5)
ggsave("Figure_1a_bottom.pdf", plot = figure_1a_bottom, width = 8, height = 5, dpi = 300,
       path = here::here("figures/manuscript_figures"), bg = "white", scale = 1)


## Manuscript text only----

#### 3rd wave cases - attributable cases since 3rd wave start ----
## larger panel - bar graph of assigned prevalence
thirdwave_bar <- ggplot(data = cases_assign) +
  geom_bar(aes(x = year, y = te_count, fill = te_assign), 
           stat = "identity", position = "stack") +
  scale_fill_manual(values = c(pub_te_palette, "Unsequenced" = "grey50"), 
                    breaks = c(te_order[-length(te_order)], "Unsequenced"),
                    drop = F) +
  ylab("Reported Cholera Cases") +
  xlab("Year") +
  labs(fill = "Transmission Event Lineage") +
  guides(fill=guide_legend(nrow=1, title.position = "left"),
         color=guide_legend(nrow=1, title.position = "left")) +
  plottheme +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9.5, face = "bold"), 
        legend.key.width = unit(0.4,"cm"))
thirdwave_bar

##case lineages assignments - donut plot
cases_assign |> 
  filter(year>= min(cases_assign$year[cases_assign$te_assign %in% c(str_c("AFR", 9:17))], na.rm = T)) |> 
  group_by(te_assign) |> 
  summarize(te_count = sum(te_count)) |> 
  ungroup() |> 
  mutate(fract = te_count/sum(te_count),
         label = round(fract, 2),
         label = if_else(label>0.02, as.character(label), "")) |> 
  arrange(te_assign) |> 
  mutate(ymax = cumsum(fract),
         ymin = ymax-fract,
         label_pos = (ymax+ymin)/2) |> 
  ggplot() +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = te_assign)) +
  geom_label(x=3.6, aes(y=label_pos, label=label), size = 5, label.size=NA, fill = "transparent", parse = T) +
  coord_polar(theta = "y") +
  xlim(c(2,4)) +
  scale_fill_manual(values = c(pub_te_palette, "Unsequenced" = "grey50"), 
                    breaks = c(te_order[-length(te_order)], "Unsequenced"),
                    drop = F) +
  labs(fill = "Transmission\n Event Lineage") +
  plottheme +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"), 
        legend.key.width = unit(1,"cm"))


#### Sequence availability before 1989----

cases_lineages_full |> 
  filter(year<1989, 
         (cases>0 & !is.na(cases))) |> 
  mutate(sequence_data = if_else(!is.na(te), 1, 0)) |> 
  distinct(country, year, sequence_data) |> 
  summarize(n_country_years = n(),
            n_sequence_yrs = sum(sequence_data),
            pct_sequence_yrs = n_sequence_yrs/n_country_years)
## 13.7% of country-years with reported cases

#### Sequence availability third wave----

cases_lineages_full |> 
  filter(year>=1989, 
         (cases>0 & !is.na(cases))) |> 
  mutate(third_wave = if_else(!is.na(te) & te %in% c(str_c("T", 9:17)), 1, 0)) |> 
  distinct(country, year, third_wave) |> 
  summarize(n_country_years = n(),
            n_sequence_yrs = sum(third_wave),
            pct_sequence_yrs = n_sequence_yrs/n_country_years)
## 22.3% of country-years with reported cases have third wave sequence data
cases_lineages_full |> 
  filter(year>=1989, 
         (cases>0 & !is.na(cases))) |> 
  mutate(sequence_data = if_else(!is.na(te), 1, 0)) |> 
  distinct(country, year, sequence_data) |> 
  summarize(n_country_years = n(),
            n_sequence_yrs = sum(sequence_data),
            pct_sequence_yrs = n_sequence_yrs/n_country_years)
## 30.8% of country-years with reported cases have any sequence data after 1988
cases_lineages_full |> 
  filter(year>=1989, !is.na(cases) & cases>0) |> 
  mutate(third_wave = if_else(!is.na(te) & te %in% c(str_c("T", 9:17)), 1, 0)) |> 
  mutate(sequence_data = if_else(!is.na(te) & third_wave==1, 1, 0)) |> 
  group_by(country) |> 
  summarize(sequence_data = max(sequence_data)) |> 
  ungroup() |> 
  filter(sequence_data == 0)
## 9 countries with no sequence data after 1988


#### Total cases ----
cases_lineages_full |> 
  filter(!is.na(cases) & cases>0) |> 
  distinct(country, year, cases) |> 
  summarize(tot_cases = sum(cases))
## 5,210,303 cases

#### All sequence availability ----

cases_lineages_full |> 
  filter(!is.na(cases) & cases>0) |> 
  mutate(sequence_data = if_else(!is.na(te), 1, 0)) |> 
  distinct(country, year, sequence_data) |> 
  summarize(n_country_years = n(),
            n_sequence_yrs = sum(sequence_data),
            pct_sequence_yrs = n_sequence_yrs/n_country_years)
## 26.2% of country-years with reported cases


#### Inferred proportions ----

inferred_by_lineage <- gen_cases_prev |> 
  group_by(te) |> 
  mutate(te_cases = sum(mean),
         te_cases_sd = sqrt(sum(sd^2))) |> 
  ungroup() |> 
  mutate(total_cases = sum(mean)) |> 
  select(te, contains("cases")) |> 
  distinct() |> 
  mutate(prev = te_cases/total_cases,
         prev_lower = (te_cases-1.96*te_cases_sd)/total_cases,
         prev_upper = (te_cases+1.96*te_cases_sd)/total_cases)
inferred_by_lineage
inferred_by_lineage |> 
  filter(te %in% c("AFR9", "AFR10")) |> 
  mutate(total_t9t10 = sum(te_cases),
         prev_t9t10 = total_t9t10/total_cases)

## T10: 32.7% [30.5, 34.8]
## T9: 13.5% [12.1, 14.9]

