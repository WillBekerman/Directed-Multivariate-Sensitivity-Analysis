
# ------------------------------------------------------------------
#     Acknowledgment: This data processing code is implemented with the help of ChatGPT.
#     It was prompted to gather and combine dietary + nutrient data from NHANES 1999–2016,
#     and to produce new datasets based on the Healthy Eating Index 2015.
#     The code was then iteratively refined and adapted to our specific needs.
# ------------------------------------------------------------------

library(haven)
library(dplyr)
library(purrr)
library(tidyr)

data_dir <- "./nhanes_diet_data"

# --- helper: detect first match from candidate list ---
pick_name <- function(cols, candidates) {
  out <- candidates[candidates %in% cols]
  if (length(out) > 0) out[1] else NA_character_
}

# --- define cycles and file paths ---
cycles <- list(
  "1999-2000" = list(dr_file = "DRXIFF.XPT",   map_file = "pyr_iff_mped1.sas7bdat", tot_file = "DRXTOT.XPT"),
  "2001-2002" = list(dr_file = "DRXIFF_B.XPT", map_file = "pyr_iff_mped1.sas7bdat", tot_file = "DRXTOT_B.XPT"),
  "2003-2004" = list(dr_file = "DR1IFF_C.XPT", map_file = "pyr_iff_d1_mped2.sas7bdat", tot_file = "DR1TOT_C.XPT"),
  "2005-2006" = list(dr_file = "DR1IFF_D.XPT", map_file = "fped_dr1iff_0506.sas7bdat", tot_file = "DR1TOT_D.XPT"),
  "2007-2008" = list(dr_file = "DR1IFF_E.XPT", map_file = "fped_dr1iff_0708.sas7bdat", tot_file = "DR1TOT_E.XPT"),
  "2009-2010" = list(dr_file = "DR1IFF_F.XPT", map_file = "fped_dr1iff_0910.sas7bdat", tot_file = "DR1TOT_F.XPT"),
  "2011-2012" = list(dr_file = "DR1IFF_G.XPT", map_file = "fped_dr1iff_1112.sas7bdat", tot_file = "DR1TOT_G.XPT"),
  "2013-2014" = list(dr_file = "DR1IFF_H.XPT", map_file = "fped_dr1iff_1314.sas7bdat", tot_file = "DR1TOT_H.XPT"),
  "2015-2016" = list(dr_file = "DR1IFF_I.XPT", map_file = "fped_dr1iff_1516.sas7bdat", tot_file = "DR1TOT_I.XPT")
)

results <- list()

# --- main processing loop ---
for (cy in names(cycles)) {
  message("\n=== Processing ", cy, " ===")
  
  dr_path  <- file.path(data_dir, cycles[[cy]]$dr_file)
  map_path <- file.path(data_dir, cycles[[cy]]$map_file)
  tot_path <- file.path(data_dir, cycles[[cy]]$tot_file)
  
  if (!file.exists(dr_path) || !file.exists(map_path) || !file.exists(tot_path)) {
    warning("Missing DR, MAP, or TOT for ", cy)
    next
  }
  
  dr  <- tryCatch(read_xpt(dr_path),  error = \(e) NULL)
  map <- tryCatch(read_sas(map_path), error = \(e) NULL)
  tot <- tryCatch(read_xpt(tot_path), error = \(e) NULL)
  if (is.null(dr) || is.null(map) || is.null(tot)) {
    warning("Could not read DR/MAP/TOT for ", cy)
    next
  }
  
  # --- normalize column cases early ---
  names(dr)  <- toupper(names(dr))
  names(map) <- toupper(names(map))
  names(tot) <- toupper(names(tot))
  
  # --- detect key + kcal variable ---
  dr_foodkey  <- pick_name(names(dr),  c("DR1IFDCD","DRDIFDCD","DRXIFDCD","FOODCODE"))
  map_foodkey <- pick_name(names(map), c("DR1IFDCD","DRDIFDCD","DRXIFDCD","FOODCODE"))
  kcal_var    <- pick_name(names(dr),  c("DR1IKCAL","DRXIKCAL","KCAL","DR1TKCAL"))
  
  message("Keys detected for ", cy, ": ", dr_foodkey, " (DR) ↔ ", map_foodkey, " (MAP)")
  
  if (any(is.na(c(dr_foodkey, map_foodkey, kcal_var)))) {
    warning("Key variables missing for ", cy, " — skipping cycle.")
    next
  }
  
  # --- harmonize join keys ---
  dr <- dr %>%
    rename(FOODKEY = all_of(dr_foodkey)) %>%
    mutate(FOODKEY = as.character(FOODKEY),
           SEQN = as.numeric(SEQN))
  
  map <- map %>%
    rename(FOODKEY = all_of(map_foodkey)) %>%
    mutate(FOODKEY = as.character(FOODKEY))
  
  # --- collapse map to unique food codes ---
  map_slim <- map %>%
    group_by(FOODKEY) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
  
  # --- join DR with mapping ---
  diet <- dr %>%
    left_join(map_slim, by = "FOODKEY", relationship = "many-to-one")
  
  # --- restore SEQN if lost during join ---
  if (!"SEQN" %in% names(diet)) {
    seqn_col <- grep("^SEQN", names(diet), value = TRUE)
    if (length(seqn_col) == 1) {
      names(diet)[names(diet) == seqn_col] <- "SEQN"
    } else {
      diet$SEQN <- dr$SEQN
    }
  }
  
  # --- Replace NAs in numeric mapping columns ---
  num_map_cols <- setdiff(names(diet)[sapply(diet, is.numeric)], names(dr))
  if (length(num_map_cols) > 0)
    diet[num_map_cols] <- lapply(diet[num_map_cols], \(x) replace_na(x, 0))
  
  # --- Aggregate to person-level ---
  agg <- diet %>%
    group_by(SEQN) %>%
    summarise(across(all_of(c(kcal_var, num_map_cols)), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename(energy = !!kcal_var)
  
  # --- detect nutrient total vars robustly ---
  sodium_var <- pick_name(names(tot), c("DR1TSODI","DRDTSODI"))
  sfat_var   <- pick_name(names(tot), c("DR1TSFAT","DRXTSFAT"))
  mufa_var   <- pick_name(names(tot), c("DR1TMFAT","DRXTMFAT"))
  pufa_var   <- pick_name(names(tot), c("DR1TPFAT","DRXTPFAT"))
  kcal_tot   <- pick_name(names(tot), c("DR1TKCAL","DRXTKCAL","KCAL","DR1IKCAL"))
  
  message("Found vars → sodium: ", sodium_var, 
          ", SFA: ", sfat_var, 
          ", MUFA: ", mufa_var, 
          ", PUFA: ", pufa_var, 
          ", kcal: ", kcal_tot)
  
  # --- pull total nutrients safely ---
  safe_pull <- function(df, var) {
    if (!is.na(var) && var %in% names(df)) df[[var]] else rep(NA_real_, nrow(df))
  }
  
  tot_keep <- tot %>%
    mutate(
      sodium_mg = safe_pull(., sodium_var),
      SFA_g     = safe_pull(., sfat_var),
      MUFA_g    = safe_pull(., mufa_var),
      PUFA_g    = safe_pull(., pufa_var),
      kcal_nutr = safe_pull(., kcal_tot)
    ) %>%
    transmute(SEQN = as.numeric(SEQN), sodium_mg, SFA_g, MUFA_g, PUFA_g, kcal_nutr)
  
  agg <- agg %>%
    left_join(tot_keep, by = "SEQN") %>%
    mutate(cycle = cy)
  
  results[[cy]] <- agg
}

# --- combine all cycles ---
nhanes_diet_combined <- bind_rows(results)

saveRDS(nhanes_diet_combined,
        file.path(data_dir, "nhanes_diet_combined_1999_2016_with_nutrients.rds"))

message("\n✅ Finished processing all available cycles (1999–2016) with total nutrients.")



















library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# -------------------------------------------------------------------
# Load combined NHANES diet dataset
# -------------------------------------------------------------------
data_dir <- "./nhanes_diet_data"
hei <- readRDS(file.path(data_dir, "nhanes_diet_combined_1999_2016_with_nutrients.rds"))

message("✅ Data loaded: ", nrow(hei), " rows")

# helper: pick available variable among candidate names
pick_var <- function(df, candidates) {
    out <- candidates[candidates %in% names(df)]
    if (length(out) == 0) return(rep(NA_real_, nrow(df)))
    df[[out[1]]]
}
hei <- hei %>%
    mutate(
        kcal = energy,
        sodium_g_per_1000 = case_when(
            is.na(sodium_mg) | is.na(kcal) | kcal <= 0 ~ NA_real_,
            TRUE ~ pmin((sodium_mg / 1000) / (kcal / 1000), 6)
        ),
        satfat_pct  = ifelse(kcal > 0, (SFA_g * 9 / kcal) * 100, NA_real_),
        fatty_ratio = ifelse(SFA_g > 0, (MUFA_g + PUFA_g) / SFA_g, NA_real_)
    )

# --- Normalize variable names for MPED vs FPED cycles ---
hei <- hei %>%
    mutate(
        F_TOTAL  = coalesce(F_TOTAL, DR1I_F_TOTAL),
        F_OTHER  = coalesce(F_OTHER, DR1I_F_OTHER),
        V_TOTAL  = coalesce(V_TOTAL, DR1I_V_TOTAL),
        V_DRKGR  = coalesce(V_DRKGR, DR1I_V_DRKGR),
        G_WHL    = coalesce(G_WHL, DR1I_G_WHOLE),
        G_NWHL   = coalesce(G_NWHL, DR1I_G_REFINED),
        D_TOTAL  = coalesce(D_TOTAL, DR1I_D_TOTAL),
        M_MPF    = coalesce(M_MPF, DR1I_PF_MPS_TOTAL),
        M_FISH_HI = coalesce(M_FISH_HI, DR1I_PF_SEAFD_HI),
        ADD_SUG  = coalesce(ADD_SUG, DR1I_ADD_SUGARS)
    )
hei <- hei %>%
    mutate(
        F_TOTAL   = coalesce(F_TOTAL, DR1I_F_TOTAL),
        F_OTHER   = coalesce(F_OTHER, DR1I_F_OTHER),
        V_TOTAL   = coalesce(V_TOTAL, DR1I_V_TOTAL),
        V_DRKGR   = coalesce(V_DRKGR, DR1I_V_DRKGR),
        G_WHL     = coalesce(G_WHL, DR1I_G_WHOLE),
        G_NWHL    = coalesce(G_NWHL, DR1I_G_REFINED),
        D_TOTAL   = coalesce(D_TOTAL, DR1I_D_TOTAL),
        M_MPF     = coalesce(M_MPF, DR1I_PF_TOTAL),
        M_FISH_HI = coalesce(M_FISH_HI, DR1I_PF_SEAFD_HI),
        ADD_SUG   = coalesce(ADD_SUG, DR1I_ADD_SUGARS)
    )

cap01 <- function(x) pmin(pmax(x, 0), 1)

hei <- hei %>%
    mutate(
        # detect MPED vs FPED style names
        F_TOTAL      = pick_var(., c("F_TOTAL", "DR1I_F_TOTAL")),
        F_OTHER      = pick_var(., c("F_OTHER", "DR1I_F_OTHER")),
        V_TOTAL      = pick_var(., c("V_TOTAL", "DR1I_V_TOTAL")),
        V_DRKGR      = pick_var(., c("V_DRKGR", "DR1I_V_DRKGR")),
        G_WHL        = pick_var(., c("G_WHL", "DR1I_G_WHOLE")),
        G_NWHL       = pick_var(., c("G_NWHL", "DR1I_G_REFINED")),
        D_TOTAL      = pick_var(., c("D_TOTAL", "DR1I_D_TOTAL")),
        M_MPF        = pick_var(., c("M_MPF", "DR1I_PF_MPS_TOTAL")),
        M_FISH_HI    = pick_var(., c("M_FISH_HI", "DR1I_PF_SEAFD_HI")),
        ADD_SUG      = pick_var(., c("ADD_SUG", "DR1I_ADD_SUGARS"))
    ) %>%
    mutate(
        # convert to ounce equivalents (for grain vars)
        G_NWHL_ozeq = G_NWHL / 28.35,
        G_WHL_ozeq  = G_WHL  / 28.35
    ) %>%
    mutate(
        # --- Adequacy components (higher = better) ---
        HEI_TotalFruit   = cap01(F_TOTAL / 0.8) * 5,
        HEI_WholeFruit   = cap01(F_OTHER / 0.4) * 5,
        HEI_TotalVeg     = cap01(V_TOTAL / 1.1) * 5,
        HEI_GreensBeans  = cap01(V_DRKGR / 0.2) * 5,
        HEI_WholeGrain   = cap01(G_WHL_ozeq / 1.5) * 10,
        HEI_Dairy        = cap01(D_TOTAL / 1.3) * 10,
        HEI_TotalProt    = cap01(M_MPF / 2.5) * 5,
        HEI_SeaPlantProt = cap01(M_FISH_HI / 0.8) * 5,
        HEI_FattyAcid    = cap01((fatty_ratio - 1.2) / (2.5 - 1.2)) * 10,
        
        # --- Moderation components (lower = better) ---
        HEI_RefinedGrain = (1 - cap01(G_NWHL_ozeq / 1.8)) * 10,
        HEI_Sodium       = (1 - cap01(sodium_g_per_1000 / 2.0)) * 10,
        HEI_AddSug       = (1 - cap01(ADD_SUG / (26 * kcal / 1000))) * 10,
        HEI_SatFat       = (1 - cap01(satfat_pct / 16)) * 10
    ) %>%
    mutate(
        HEI_Total = 
            HEI_TotalFruit + HEI_WholeFruit + HEI_TotalVeg + HEI_GreensBeans +
                HEI_WholeGrain + HEI_Dairy + HEI_TotalProt + HEI_SeaPlantProt +
                HEI_FattyAcid + HEI_RefinedGrain + HEI_Sodium + HEI_AddSug + HEI_SatFat
    )

# # ------------------------------------------------------------------
# # Summary statistics for HEI components
# # ------------------------------------------------------------------
# hei_summary <- hei %>%
#   summarise(across(starts_with("HEI_"), ~mean(.x, na.rm = TRUE))) %>%
#   pivot_longer(everything(), names_to = "component", values_to = "mean_score") %>%
#   arrange(desc(mean_score))
#
# print(hei_summary, n = 20)
#
# # ------------------------------------------------------------------
# # Diagnostics: distribution and component-level stats (fixed)
# # ------------------------------------------------------------------
# hei_diagnostics <- hei %>%
#   summarise(across(
#     starts_with("HEI_"),
#     list(
#       n    = ~sum(!is.na(.x)),
#       mean = ~mean(.x, na.rm = TRUE),
#       sd   = ~sd(.x, na.rm = TRUE),
#       p0   = ~mean(.x == 0, na.rm = TRUE),
#       p10  = ~quantile(.x, 0.1, na.rm = TRUE),
#       min  = ~min(.x, na.rm = TRUE),
#       max  = ~max(.x, na.rm = TRUE)
#     )
#   )) %>%
#   pivot_longer(
#     everything(),
#     names_to = c("component", ".value"),
#     names_pattern = "(HEI_[^_]+)_(.*)"
#   )
#
# print(hei_diagnostics, n = 20)
#
# # ------------------------------------------------------------------
# # Visualization
# # ------------------------------------------------------------------
# ggplot(hei_summary, aes(x = reorder(component, mean_score), y = mean_score)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Mean HEI-2015 Component Scores (Unweighted)",
#        x = NULL, y = "Mean score (0–10 or 0–5)") +
#   theme_minimal(base_size = 12)
#
# # Histogram of selected components
# comp_vars <- c("HEI_Total", "HEI_TotalFruit", "HEI_WholeFruit", "HEI_Dairy",
#                "HEI_WholeGrain", "HEI_AddSug", "HEI_SatFat")
#
# hei %>%
#   select(all_of(comp_vars)) %>%
#   pivot_longer(everything(), names_to = "component", values_to = "score") %>%
#   ggplot(aes(x = score)) +
#   geom_histogram(binwidth = 1, fill = "seagreen3", color = "white") +
#   facet_wrap(~component, scales = "free_y") +
#   labs(title = "HEI-2015 Component Distributions (Unweighted)",
#        x = "Score", y = "Count") +
#   theme_minimal(base_size = 11)
#
# # ------------------------------------------------------------------
# # Nutrient sanity check
# # ------------------------------------------------------------------
# nutr_sanity <- hei %>%
#   summarise(across(
#     c(sodium_g_per_1000, satfat_pct, fatty_ratio),
#     list(mean = mean, sd = sd, min = min, max = max),
#     na.rm = TRUE
#   )) %>%
#   pivot_longer(everything(),
#                names_to = c("variable", ".value"),
#                names_pattern = "(.*)_(mean|sd|min|max)")
#
# print(nutr_sanity)
# cat("\nExpected ranges:\n",
#     "• sodium_g_per_1000 ≈ 1.3–2.8\n",
#     "• satfat_pct ≈ 8–18\n",
#     "• fatty_ratio ≈ 1.3–2.5\n")

# # ------------------------------------------------------------------
# # Total HEI distribution plot
# # ------------------------------------------------------------------
# ggplot(hei, aes(x = HEI_Total)) +
#   geom_histogram(binwidth = 2, fill = "steelblue", color = "white") +
#   geom_vline(xintercept = 47, linetype = "dashed", color = "red") +
#   annotate("text", x = 47, y = 600, label = "USDA mean ≈ 47",
#            color = "red", hjust = -0.1) +
#   labs(title = "Distribution of HEI-2015 Total Scores (Unweighted)",
#        x = "HEI-2015 Total Score", y = "Count") +
#   theme_minimal(base_size = 12)

message("\n✅ HEI-2015 scoring completed (unweighted, validated).")

# hei %>%
#   group_by(cycle) %>%
#   summarise(prop_missing_HEI_Total = mean(is.na(HEI_Total)))


# Define cycle order and map to numeric wave
cycle_order <- c("1999-2000","2001-2002","2003-2004","2005-2006",
                 "2007-2008","2009-2010","2011-2012","2013-2014","2015-2016")

hei <- hei %>%
  mutate(
    cycle = factor(cycle, levels = cycle_order),
    WAVE = as.integer(cycle)
  )

diet_file_final <- hei %>% dplyr::select(SEQN, WAVE, HEI_Total)
saveRDS(diet_file_final,
        file.path(data_dir, "diet_file_final.rds"))

diet_file_final_big <- hei %>% dplyr::select(SEQN, WAVE, HEI_TotalFruit, HEI_WholeFruit, HEI_TotalVeg, HEI_GreensBeans,
                                             HEI_WholeGrain, HEI_Dairy, HEI_TotalProt, HEI_SeaPlantProt,
                                             HEI_FattyAcid, HEI_RefinedGrain, HEI_Sodium, HEI_AddSug, HEI_SatFat, HEI_Total)
saveRDS(diet_file_final_big,
        file.path(data_dir, "diet_file_final_big.rds"))


