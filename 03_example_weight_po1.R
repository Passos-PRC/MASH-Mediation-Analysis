source("R/01_packages.R")
source("R/02_functions.R")

WEIGHT <- readxl::read_excel(
  Sys.getenv("MEDIATOR_EXTRACTION_PATH",
             unset = "data/raw/MEDIATORS_DATA_EXTRACTION.xlsx"),
  sheet = "Weight"
)

source("R/02_functions.R")

study_meta <- tibble::tribble(
  ~study_id,        ~phase3, ~f4_included, ~at_risk_only, ~followup_weeks,
  "IMPACT",         FALSE, FALSE, TRUE,  48,
  "NCT02970942",    FALSE, FALSE, FALSE, 72,
  "ION224-CS2",     FALSE, FALSE, FALSE, 64,
  "TANDEM",         FALSE, TRUE,  FALSE, 48,
  "ATLAS",          FALSE, TRUE,  FALSE, 48,
  "HARMONY",        FALSE, FALSE, TRUE,  96,
  "DEAN",           TRUE,  TRUE,  FALSE, 48,
  "ESSENCE",        TRUE,  FALSE, TRUE,  72,
  "NCT02962297",    FALSE, TRUE,  FALSE, 24,
  "NCT04906421",    FALSE, FALSE, TRUE,  52,
  "MAESTRO-NASH",   TRUE,  FALSE, FALSE, 52,
  "ALPINE 4",       FALSE, TRUE,  FALSE, 48,
  "ENLIVEN",        FALSE, FALSE, TRUE,  24,
  "NCT03987451",    FALSE, TRUE,  FALSE, 48,
  "DESTINY 1",      FALSE, FALSE, FALSE, 38,
  "NATIVE",         FALSE, FALSE, FALSE, 24,
  "NCT02279524",    FALSE, FALSE, FALSE, 65,
  "NCT01068444",    FALSE, TRUE,  FALSE, 36,
  "NCT02912260",    FALSE, FALSE, FALSE, 38,
  "NCT02443116",    FALSE, FALSE, TRUE,  30,
  "NCT02787304",    FALSE, FALSE, FALSE, 52,
  "STELLAR 3",      TRUE,  TRUE,  FALSE, 48,
  "STELLAR 4",      TRUE,  TRUE,  FALSE, 48,
  "CENTAUR",        FALSE, FALSE, FALSE, 52,
  "FALCON 1",       FALSE, FALSE, TRUE,  48,
  "LEAN",           FALSE, TRUE,  FALSE, 60,
  "NCT03976401",    FALSE, TRUE,  FALSE, 20,
  "MT-3995",        FALSE, FALSE, TRUE,  80,
  "SYNERGY",        FALSE, FALSE, TRUE,  52,
  "ALPINE 2/3",     FALSE, FALSE, TRUE,  24,
  "MIRNA",          FALSE, FALSE, TRUE,  48,
  "ENCORE-NF",      FALSE, FALSE, FALSE, 72,
  "EMMINENCE",      FALSE, FALSE, FALSE, 52
) %>%
  mutate(ge_48w = followup_weeks >= 48)

study_meta_mod <- tibble::tribble(
  ~study_id,        ~followup_weeks, ~t2dm_raw,  ~bmi_raw, ~fibrosis_raw,
  "IMPACT",         48,  "43%",     "38.7",  "F2-F3",
  "NCT02970942",    72,  "62.2",    "35.63", "F1-F3",
  "ION224-CS2",     64,  "51%",     "37.82", "F1-F3",
  "TANDEM",         48,  "82,3%",   "34.59", "F1-F4",
  "ATLAS",          48,  "72%",     "34.12", "F3-F4",
  "HARMONY",        96,  "70%",     "38.06", "F2-F3",
  "DEAN",           48,  "45%",     "29.15", "F0-F4",
  "ESSENCE",        72,  "55,90%",  "34.53", "F2-F3",
  "NCT02962297",    24,  "0%",      "26",    "F0-F4",
  "NCT04906421",    52,  "61%",     "35",    "F2-F3",
  "MAESTRO-NASH",   52,  "67%",     "35.66", "F1-F3",
  "ALPINE 4",       48,  "75.6",    "34.86", "F3-F4",
  "ENLIVEN",        24,  "66",      "36.6",  "F2-F3",
  "NCT03987451",    48,  "75",      "34.9",  "F4",
  "DESTINY 1",      38,  "41,03%",  "36.1",  "F1-F3",
  "NATIVE",         24,  "41.70%",  "23.9",  "F0-F3",
  "NCT02279524",    65,  "68,82%",  "32.7",  "F0-F3",
  "NCT01068444",    36,  "23,33%",  "28.9",  "F0-F4",
  "NCT02912260",    38,  "39,20%",  "35.1",  "F1-F3",
  "NCT02443116",    30,  "61,54%",  "36.1",  "F2-F3",
  "NCT02787304",    52,  "43,37%",  "34.5",  "F0-F3",
  "STELLAR 3",      48,  "70.2%",   "32.75", "F3",
  "STELLAR 4",      48,  "76.9%",   "33.41", "F4",
  "CENTAUR",        52,  "50,50%",  "33.9",  "F1-F3",
  "FALCON 1",       48,  "73,70%",  "35.6",  "F3",
  "LEAN",           60,  "32,70%",  "36",    "F2-F4",
  "NCT03976401",    20,  "50,00%",  "37",    "F4",
  "MT-3995",        80,  "42,60%",  "29.8",  "F2-F3",
  "SYNERGY",        52,  "58%",     "36.1",  "F2-F3",
  "ALPINE 2/3",     24,  "49%",     "38.1",  "F2-F3",
  "MIRNA",          48,  "67,50%",  "32.4",  "F2-F3",
  "ENCORE-NF",      72,  "50.6",    "34",    "F1-F3",
  "EMMINENCE",      52,  "52.3",    "35.16", "F1-F3"
) %>%
  mutate(
    study_id    = stringr::str_to_upper(stringr::str_trim(study_id)),
    t2dm_pct    = as.numeric(stringr::str_remove(stringr::str_replace_all(t2dm_raw, ",", "."), "%")),
    bmi         = as.numeric(stringr::str_replace_all(bmi_raw, ",", ".")),
    fibrosis_f4 = as.integer(stringr::str_detect(fibrosis_raw, "F4"))
  )

output_dir <- Sys.getenv("MEDIATOR_OUTPUT_DIR", unset = "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

prep           <- prepare_mediation_data_1(WEIGHT, study_meta, study_meta_mod)
contrasts      <- prep$contrasts
contrasts_long <- prep$contrasts_long %>%
  mutate(group_type = as.character(group_type), group = as.character(group))
all_groups     <- prep$all_groups %>%
  mutate(group_type = as.character(group_type), group = as.character(group))
desired_order  <- prep$desired_order

b_fit <- fit_bpath_bayes(
  df   = contrasts,
  mods = ~ wt_loss_eff + factor(drug_class),
  file = file.path(output_dir, "bfit_weight_po1_cache")
)

print(b_fit$summary)

b_global <- tibble::tibble(
  b               = b_fit$summary$estimate,
  b_se            = b_fit$summary$se,
  b_lower         = b_fit$summary$lower,
  b_upper         = b_fit$summary$upper,
  BF10            = b_fit$summary$BF10,
  b_RD_per5       = b_fit$summary$estimate * 5,
  b_RD_per5_lower = b_fit$summary$lower    * 5,
  b_RD_per5_upper = b_fit$summary$upper    * 5
)

med <- compute_point_estimates_bayes(
  contrasts      = contrasts,
  contrasts_long = contrasts_long,
  all_groups     = all_groups,
  b_fit          = b_fit
)

summary_table <- med %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3)),
    a_ci           = paste0(round(a, 3),           " [", round(a_lower,        3), ", ", round(a_upper,        3), "]"),
    indirect_rd_ci = paste0(round(indirect_rd, 3), " [", round(indirect_lower, 3), ", ", round(indirect_upper, 3), "]"),
    direct_rd_ci   = paste0(round(direct_rd, 3),   " [", round(direct_lower,   3), ", ", round(direct_upper,   3), "]"),
    total_rd_ci    = paste0(round(total_rd, 3),     " [", round(total_lower,    3), ", ", round(total_upper,    3), "]"),
    prop_indirect_pct = round(prop_indirect_pct, 1)
  ) %>%
  dplyr::select(
    group, k, n_studies, total_n,
    a_ci, indirect_rd_ci, prop_indirect_pct,
    total_rd_ci, direct_rd_ci,
    BF10_a, BF10_b, BF10_total, BF10_indirect
  ) %>%
  arrange(desc(abs(prop_indirect_pct)))

print(summary_table, n = Inf)

plot_bundle <- make_forest_plot(med = med, desired_order = desired_order, b_global = b_global)
p_bubble    <- make_bubble_plot(contrasts = contrasts, b_global = b_global)

ggsave(file.path(output_dir, "FOREST_WEIGHT_PO1.pdf"), plot = plot_bundle$p_final, width = 14, height = 6)
ggsave(file.path(output_dir, "BUBBLE_WEIGHT_PO1.pdf"), plot = p_bubble,            width = 9,  height = 6)

write.csv(summary_table, file.path(output_dir, "mediation_weight_po1.csv"), row.names = FALSE)
