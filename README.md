# cholera-units
Cholera Transmission Units in Africa

Code accompanying:

*Defining Epidemiologically Relevant Units of Cholera Transmission in sub-Saharan Africa*
**Bethany L. DiPrete, Javier Perez-Saez, Shirlee Wohl, Nathaniel Matteson, Seungwon Kim, Andrew S. Azman, Justin Lessler**
medRxiv 2025.06.06.25329161; doi: https://doi.org/10.1101/2025.06.06.25329161
---

## Overview

This repository contains the data processing, Bayesian modeling, simulation, and visualization workflows used to investigate how the choice of epidemiological unit influences inference about cholera transmission dynamics across sub-Saharan Africa.

The analysis combines:

- WHO annual cholera case reports
- Lineage calls from *Vibrio cholerae* genomic sequence data
- Output from phylogeographic reconstructions
- Bayesian Hidden Markov Models (HMMs) implemented in Stan
- Simulation-based analyses of introductions and transmission

---

## Repository Structure

```text
cholera-units/
├── analysis/
│   ├── 00_functions_settings.R
│   ├── 01_initial_data_prep.R
│   ├── 02_import_data.R
│   ├── 03_create_partitions.R
│   ├── 04_hmm_setup.R
│   ├── 05_run_hmm.R
│   ├── 06_post_hmm_stats.R
│   ├── 07_hmm_clustering.R
│   ├── 08_phylogeography_clustering.R
│   ├── 09_sim_outbreaks_poststan.R
│   ├── 10_arrival_times.R
│   ├── 11_downsampling_analysis.R
│   ├── 12_Figure1.R
│   ├── 13_Figure3.R
│   ├── 14_Figure4.R
│   ├── 15_Figure5.R
│   ├── 16_Supplemental_Figures.R
│   ├── run_sims/
│   └── stan/
│       └── hmm_fb_model_connectivity.stan
│
├── data/
│   └── raw/
│       ├── who_annual_case_reports.csv
│       ├── sequences.csv
│       ├── introductions_w.csv
│       ├── UN_regions.csv
│       └── phylogeography/
│
├── notebooks/
│   ├── 01_hmm_results.Rmd
│   ├── 02_hmm_diagnostics.Rmd
│   └── generate_reports/
│
└── README.md
```

---

## Analysis Workflow

### 1. Data Preparation

The workflow begins by harmonizing surveillance and genomic data.

| Script | Purpose |
|----------|----------|
| `01_initial_data_prep.R` | Initial cleaning and preparation |
| `02_import_data.R` | Import surveillance, sequence, and geographic data |
| `03_create_partitions.R` | Construct partitions used in downsampling analyses |

### 2. HMM Setup and Fitting

A Bayesian Hidden Markov Model is used to estimate cholera transmission dynamics across countries, years, and lineages.

| Script | Purpose |
|----------|----------|
| `04_hmm_setup.R` | Construct Stan inputs and connectivity matrices |
| `05_run_hmm.R` | Fit the model using CmdStanR |
| `analysis/stan/hmm_fb_model_connectivity.stan` | Stan model implementation |

### 3. Post-Processing

Posterior estimates are summarized and compared across epidemiological units.

| Script | Purpose |
|----------|----------|
| `06_post_hmm_stats.R` | Posterior summaries |
| `07_hmm_clustering.R` | Clustering based on HMM results |
| `08_phylogeography_clustering.R` | Clustering based on phylogeographic results |
| `09_sim_outbreaks_poststan.R` | Outbreak simulations |
| `10_arrival_times.R` | Introduction and arrival-time analyses |
| `11_downsampling_analysis.R` | Sensitivity analyses for sampling intensity |

### 4. Figure Generation

| Script | Output |
|----------|----------|
| `12_Figure1.R` | Figure 1 |
| `13_Figure3.R` | Figure 3 |
| `14_Figure4.R` | Figure 4 |
| `15_Figure5.R` | Figure 5 |
| `16_Supplemental_Figures.R` | Supplementary figures |

---

## Data Sources

The repository includes raw inputs used throughout the analysis:

| File | Description |
|----------|----------|
| `who_annual_case_reports.csv` | Annual cholera case reports |
| `sequences.csv` | Curated genomic sequence metadata |
| `introductions_w.csv` | Introduction-weight data used in modeling |
| `UN_regions.csv` | Regional classifications |
| `phylogeography/*` | BEAST phylogeographic outputs and lineage summaries |

Please refer to the manuscript and supplementary materials for complete details on data provenance and inclusion criteria.

---

## Software Requirements

### Core Software

- R 4.5.2
- CmdStan
- CmdStanR 0.9.0

### Key Packages

| Package | Version |
|----------|----------|
| tidyverse | 2.0.0 |
| ggplot2 | 4.0.1 |
| cmdstanr | 0.9.0 |
| posterior | 1.6.1 |
| coda | 0.19-4.1 |
| sf | 1.0-24 |
| cowplot | 1.2.0 |
| ggpattern | 1.2.1 |
| plotly | 4.11.0 |

A complete software manifest, including all attached and namespace-loaded packages, is provided in [`sessionInfo.md`](sessionInfo.md).

---

## Reproducing the Analysis

Clone the repository:

```bash
git clone https://github.com/UNCIDD/cholera-units.git
cd cholera-units
```

Install dependencies and CmdStan, then execute the analysis scripts in numerical order.

For model diagnostics and result summaries, see:

- `notebooks/01_hmm_results.Rmd`
- `notebooks/02_hmm_diagnostics.Rmd`

---

## Computational Environment

Analyses for the manuscript were generated using R 4.5.2 on macOS Sonoma 14.6.1. See [`sessionInfo.md`](sessionInfo.md) for the complete software environment.

---

## Citation

```text
Defining Epidemiologically Relevant Units of Cholera Transmission in sub-Saharan Africa
Bethany L. DiPrete, Javier Perez-Saez, Shirlee Wohl, Nathaniel Matteson, Seungwon Kim, Andrew S. Azman, Justin Lessler
medRxiv 2025.06.06.25329161; doi: https://doi.org/10.1101/2025.06.06.25329161
```

## License

See the repository LICENSE file for licensing information.
