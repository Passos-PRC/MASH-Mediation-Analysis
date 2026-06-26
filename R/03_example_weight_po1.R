source("R/01_packages.R")
source("R/02_functions.R")

WEIGHT <- readxl::read_excel("data/MEDIATORS_SYNTHETIC.xlsx", sheet = "Weight")

study_meta <- tibble::tribble(
  ~study_id,      ~phase3, ~f4_included, ~at_risk_only, ~followup_weeks,
  "IMPACT",       FALSE,   FALSE,        TRUE,          48,
  "NCT02970942",  FALSE,   FALSE,        FALSE,         72,
  "ION224-CS2",   FALSE,   FALSE,        FALSE,         64,
  "TANDEM",       FALSE,   TRUE,         FALSE,         48,
  "ATLAS",        FALSE,   TRUE,         FALSE,         48,
  "HARMONY",      FALSE,   FALSE,        TRUE,          96,
  "DEAN",         TRUE,    TRUE,         FALSE,         48,
  "ESSENCE",      TRUE,    FALSE,        TRUE,          72,
  "MAESTRO-NASH", TRUE,    FALSE,        FALSE,         52,
  "CENTAUR",      FALSE,   FALSE,        FALSE,         52
) %>%
  dplyr::mutate(ge_48w = followup_weeks >= 48)

output_dir <- Sys.getenv("MEDIATOR_OUTPUT_DIR", unset = "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

prep           <- prepare_mediation_data_1(WEIGHT, study_meta)
contrasts      <- prep$contrasts
contrasts_long <- prep$contrasts_long %>%
  dplyr::mutate(group_type = as.character(group_type), group = as.character(group))
all_groups     <- prep$all_groups %>%
  dplyr::mutate(group_type = as.character(group_type), group = as.character(group))
desired_order  <- prep$desired_order

b_fit <- fit_bpath_bayes(
  df   = contrasts,
  mods = ~ wt_loss_eff + factor(drug_class),
  file = file.path(output_dir, "bfit_weight_po1_cache")
)

med <- compute_point_estimates_bayes(
  contrasts      = contrasts,
  contrasts_long = contrasts_long,
  all_groups     = all_groups,
  b_fit          = b_fit
)

fmt_ci <- function(m, lo, hi, digits = 3) {
  sprintf("%s [%s, %s]", round(m, digits), round(lo, digits), round(hi, digits))
}

summary_table <- med %>%
  dplyr::transmute(
    group,
    k, n_studies, total_n,
    a_ard_95cri        = fmt_ci(a,        a_lower,        a_upper),
    b_ard_95cri        = fmt_ci(b,        b_lower,        b_upper),
    indirect_ard_95cri = fmt_ci(indirect, indirect_lower, indirect_upper),
    direct_ard_95cri   = fmt_ci(direct,   direct_lower,   direct_upper),
    total_ard_95cri    = fmt_ci(total,    total_lower,    total_upper),
    BF10_a, BF10_b, BF10_indirect, BF10_total
  )

print(summary_table, n = Inf)

b_global <- dplyr::tibble(
  b       = b_fit$summary$estimate,
  b_lower = b_fit$summary$lower,
  b_upper = b_fit$summary$upper
)

p_forest <- make_forest_plot(med = med, desired_order = desired_order, b_global = b_global)
p_bubble <- make_bubble_plot(contrasts = contrasts, b_global = b_global)

ggplot2::ggsave(file.path(output_dir, "FOREST_WEIGHT_PO1.pdf"), plot = p_forest, width = 14, height = 6)
ggplot2::ggsave(file.path(output_dir, "BUBBLE_WEIGHT_PO1.pdf"), plot = p_bubble, width = 9,  height = 6)

write.csv(summary_table, file.path(output_dir, "mediation_weight_po1.csv"), row.names = FALSE)
