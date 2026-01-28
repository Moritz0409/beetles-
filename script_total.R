# ============================================================
# TWO-STAGE FOREST BIODIVERSITY INDEX
#
# Stage 1: Beetle-based biodiversity value per tree species (BV_i)
# Stage 2: Forest composition biodiversity index
#
# Extended version:
# - W1_linear to W10_linear
# - Sigma = 0.4 (CES), 0.7 (CES), 1.0 (Cobb–Douglas)
#
# Output:
# - CES_results_only.csv
# - Forest_Biodiversity_Index_All_Scenarios_AllWeights_AllSigmas.csv
# ============================================================


#############################################################
# GLOBAL SETUP
#############################################################

setwd("/Users/Jonas/Desktop/Master/3rd Semester/Integration Module/BEETLES/Input_Output")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(tibble)


#############################################################
# STAGE 1 — BEETLE-BASED TREE SPECIES BIODIVERSITY (BV_i)
#############################################################

# -----------------------------------------------------------
# 1. Read beetle input data
# -----------------------------------------------------------

beetle_data <- read.csv(
  "Index Total Mean CES.csv",
  row.names = 1,
  check.names = FALSE
)

# -----------------------------------------------------------
# 2. Select disjoint beetle categories
# -----------------------------------------------------------

ces_data <- beetle_data[c(
  "unique red list beetle species",
  "Non-unique red list beetle species",
  "unique non-red list beetle species",
  "Common beetles"
), ]

rownames(ces_data) <- c(
  "unique_red_list",
  "non_unique_red_list",
  "unique_non_red_list",
  "common"
)

# -----------------------------------------------------------
# 3. Aggregation functions
# -----------------------------------------------------------

ces_general <- function(x, a, sigma) {
  rho <- (sigma - 1) / sigma
  (sum(a * (x ^ rho)))^(1 / rho)
}

ces_cobb_douglas <- function(x, a) {
  exp(sum(a * log(x)))
}

ces_function <- function(x, a, sigma) {
  if (abs(sigma - 1) < 1e-8) {
    ces_cobb_douglas(x, a)
  } else {
    ces_general(x, a, sigma)
  }
}

linear_index <- function(x, a) {
  sum(a * x)
}

# -----------------------------------------------------------
# 4. Normative weight scenarios
# -----------------------------------------------------------

weights_list <- list(
  W1  = c(0.25, 0.25, 0.25, 0.25),
  W2  = c(0.50, 0.30, 0.15, 0.05),
  W3  = c(0.35, 0.25, 0.25, 0.15),
  W4  = c(0.40, 0.10, 0.40, 0.10),
  W5  = c(0.45, 0.05, 0.40, 0.10),
  W6  = c(0.40, 0.30, 0.25, 0.05),
  W7  = c(0.30, 0.25, 0.25, 0.20),
  W8  = c(0.60, 0.35, 0.03, 0.02),
  W9  = c(0.50, 0.05, 0.40, 0.05),
  W10 = c(0.15, 0.15, 0.20, 0.50)
)

# -----------------------------------------------------------
# 5. Compute BV_i (ALL weights, LINEAR ONLY)
# -----------------------------------------------------------

results <- list()

for (w in names(weights_list)) {
  results[[paste0(w, "_linear")]] <- sapply(
    colnames(ces_data),
    function(tree) linear_index(ces_data[, tree], weights_list[[w]])
  )
}

ces_results_df <- as.data.frame(results)
rownames(ces_results_df) <- colnames(ces_data)

write.csv(
  round(ces_results_df, 3),
  "CES_results_only.csv"
)


#############################################################
# STAGE 2a — BASAL AREA CALCULATION FROM PPA SIMULATIONS
#############################################################

# -----------------------------------------------------------
# 1. Species lookup table
# -----------------------------------------------------------

species_lookup <- tibble(
  sp = 1:8,
  species = c(
    "Acer pseudoplatanus",
    "Acer campestre",
    "Fraxinus Excelsior",
    "Carpinus betulus",
    "Acer platanoides",
    "Quercus robur",
    "Ulmus spp.",
    "Tilia spp."
  )
)

# -----------------------------------------------------------
# 2. Mapping of PPA scenarios
# -----------------------------------------------------------

scenarios <- tibble(
  file = c(
    "Results_combined_census_dry_init_dry.csv",
    "Results_combined_census_middle_init_dry.csv",
    "Results_combined_census_moist_init_dry.csv",
    "Results_combined_census_middle_init_intermediate.csv",
    "Results_combined_census_moist_init_intermediate.csv",
    "Results_combined_census_wet_init_intermediate.csv",
    "Results_combined_census_moist_init_moist.csv",
    "Results_combined_census_wet_init_moist.csv",
    "Results_combined_census_v_wet_init_moist.csv"
  ),
  initial_state = c(
    "Dry","Dry","Dry",
    "Intermediate","Intermediate","Intermediate",
    "Moist","Moist","Moist"
  ),
  groundwater_change_m = c(
    0.0,0.5,1.0,
    0.0,0.5,1.0,
    0.0,0.5,1.0
  )
)

# -----------------------------------------------------------
# 3. Compute basal area
# -----------------------------------------------------------

compute_basal_area <- function(file, init, delta_gw) {
  
  read_csv(file, show_col_types = FALSE) %>%
    filter(cl == 1, time %in% c(0,25,50,75,100)) %>%
    mutate(
      dbh_m = dbh / 100,
      basal_area_m2_ha = n * pi * (dbh_m / 2)^2
    ) %>%
    group_by(time, sp) %>%
    summarise(
      basal_area_m2_ha = sum(basal_area_m2_ha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(species_lookup, by = "sp") %>%
    mutate(
      initial_state = init,
      groundwater_change_m = delta_gw
    )
}

all_results <- pmap_dfr(
  list(
    scenarios$file,
    scenarios$initial_state,
    scenarios$groundwater_change_m
  ),
  compute_basal_area
)

#############################################################
# STAGE 2b — FOREST BIODIVERSITY INDEX (ALL WEIGHTS × SIGMAS)
#############################################################

species_mapping <- c(
  "Acer pseudoplatanus" = "Acer",
  "Fraxinus Excelsior"  = "Fraxinus",
  "Quercus robur"       = "Quercus",
  "Tilia spp."          = "Tilia",
  "Ulmus spp."          = "Ulmus"
)

valid_forest_species <- names(species_mapping)
valid_beetle_species <- unname(species_mapping)

# -----------------------------------------------------------
# Forest aggregation function (CES + Cobb–Douglas)
# -----------------------------------------------------------

forest_aggregate <- function(x, sigma) {
  
  a <- rep(1 / length(x), length(x))
  
  if (abs(sigma - 1) < 1e-8) {
    exp(sum(a * log(x)))   # Cobb–Douglas
  } else {
    rho <- (sigma - 1) / sigma
    (sum(a * (x ^ rho)))^(1 / rho)
  }
}

sigma_scenarios <- tibble(
  aggregation = c("CES_sigma_0.4","CES_sigma_0.7","Cobb_Douglas_sigma_1"),
  sigma = c(0.4, 0.7, 1.0)
)

forest_files <- tibble(
  groundwater_scenario = c("no_increase","0.5m_increase","1m_increase"),
  gw_value = c(0.0,0.5,1.0)
)

forest_results <- list()

for (gw in forest_files$gw_value) {
  
  forest_mean <- all_results %>%
    filter(
      groundwater_change_m == gw,
      initial_state %in% c("Dry","Intermediate","Moist"),
      species %in% valid_forest_species
    ) %>%
    group_by(species, time) %>%
    summarise(
      basal_area_m2_ha = mean(basal_area_m2_ha),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = basal_area_m2_ha,
      names_prefix = "year_"
    ) %>%
    column_to_rownames("species")
  
  rownames(forest_mean) <- species_mapping[rownames(forest_mean)]
  forest_mean <- forest_mean[valid_beetle_species, ]
  
  years <- as.numeric(gsub("year_", "", colnames(forest_mean)))
  
  for (w in paste0("W",1:10,"_linear")) {
    
    BV_i <- ces_results_df[valid_beetle_species, w, drop = TRUE]
    
    for (i in 1:nrow(sigma_scenarios)) {
      
      agg <- sigma_scenarios$aggregation[i]
      sigma_val <- sigma_scenarios$sigma[i]
      
      vals <- sapply(forest_mean, function(BA_vec) {
        
        x_i <- BV_i * as.numeric(BA_vec)
        
        if (agg == "linear") {
          mean(x_i)
        } else {
          forest_aggregate(x_i, sigma_val)
        }
      })
      
      forest_results[[length(forest_results) + 1]] <- tibble(
        groundwater_scenario = forest_files$groundwater_scenario[
          forest_files$gw_value == gw
        ],
        year = years,
        weight_scenario = w,
        aggregation = agg,
        forest_biodiversity_index = vals
      )
    }
  }
}

forest_biodiversity_all <- bind_rows(forest_results)

write_csv(
  forest_biodiversity_all,
  "Forest_Biodiversity_Index_All_Scenarios_AllWeights_AllSigmas.csv"
)

forest_biodiversity_all

#############################################################
# STAGE 3 — HUMAN-READABLE COMPARISON TABLES
#############################################################

for (agg_name in unique(forest_biodiversity_all$aggregation)) {
  
  out <- forest_biodiversity_all %>%
    filter(aggregation == agg_name) %>%
    mutate(
      groundwater_scenario = recode(
        groundwater_scenario,
        "no_increase"   = "no_increase",
        "0.5m_increase" = "gw_0.5m",
        "1m_increase"   = "gw_1m"
      )
    ) %>%
    pivot_wider(
      names_from  = groundwater_scenario,
      values_from = forest_biodiversity_index
    ) %>%
    arrange(weight_scenario, year)
  
  write_csv(
    out,
    paste0("Forest_Biodiversity_", agg_name, ".csv")
  )
}

#############################################################
# STAGE 4 — Sensitivity Analysis
#############################################################

# -----------------------------------------------------------
# 1. Beetle Index
# -----------------------------------------------------------

# -----------------------------------------------------------
# Read beetle BV_i results
# -----------------------------------------------------------

beetle_df <- read_csv("CES_results_only.csv", show_col_types = FALSE)

weight_cols <- grep("^W[0-9]+_linear$", colnames(beetle_df), value = TRUE)

# -----------------------------------------------------------
# Compute ranks per weighting
# -----------------------------------------------------------

beetle_ranks <- beetle_df %>%
  select(all_of(weight_cols)) %>%
  mutate(across(everything(), rank, ties.method = "average"))

# -----------------------------------------------------------
# Spearman rank correlation between weightings
# -----------------------------------------------------------

beetle_spearman <- cor(
  beetle_ranks,
  method = "spearman"
)

# ensure rownames are preserved as first column
beetle_spearman_df <- beetle_spearman %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "weight_scenario")

write_csv(
  beetle_spearman_df,
  "Sensitivity_Beetle_Rank_Spearman_Matrix.csv"
)


# -----------------------------------------------------------
# 2. Forest Index Ranking
# -----------------------------------------------------------

rank_forest_csv <- function(input_file, output_file) {
  
  df <- read_csv(input_file, show_col_types = FALSE)
  
  df_ranked <- df %>%
    rowwise() %>%
    mutate(
      ranks = list(
        rank(
          -c(no_increase, gw_0.5m, gw_1m),  # <-- DESCENDING
          ties.method = "average"
        )
      ),
      no_increase = ranks[[1]],
      gw_0.5m     = ranks[[2]],
      gw_1m       = ranks[[3]]
    ) %>%
    select(-ranks) %>%
    ungroup()
  
  write_csv(df_ranked, output_file)
}

rank_forest_csv(
  "Forest_Biodiversity_CES_sigma_0.4.csv",
  "Forest_Biodiversity_CES_sigma_0.4_RANKS.csv"
)

rank_forest_csv(
  "Forest_Biodiversity_CES_sigma_0.7.csv",
  "Forest_Biodiversity_CES_sigma_0.7_RANKS.csv"
)

rank_forest_csv(
  "Forest_Biodiversity_Cobb_Douglas_sigma_1.csv",
  "Forest_Biodiversity_Cobb_Douglas_sigma_1_RANKS.csv"
)
