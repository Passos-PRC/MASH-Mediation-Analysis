PRIOR_B_SD   <- 0.5
PRIOR_INT_SD <- 0.5
PRIOR_TAU_SD <- 0.25

.default_priors <- function() {
  c(
    set_prior(sprintf("normal(0, %.2f)", PRIOR_B_SD),   class = "b"),
    set_prior(sprintf("normal(0, %.2f)", PRIOR_INT_SD), class = "Intercept"),
    set_prior(sprintf("normal(0, %.2f)", PRIOR_TAU_SD), class = "sd")
  )
}

clean_id <- function(x) {
  x <- str_trim(str_to_upper(as.character(x)))
  case_when(
    str_detect(x, "STELLAR 3") ~ "STELLAR 3",
    str_detect(x, "STELLAR 4") ~ "STELLAR 4",
    TRUE ~ x
  )
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

smry_draws <- function(x) {
  list(
    m = mean(x, na.rm = TRUE),
    l = safe_quantile(x, 0.025),
    u = safe_quantile(x, 0.975)
  )
}

add_within_between <- function(df) {
  df %>%
    group_by(study_id) %>%
    mutate(
      wt_mean   = mean(wt_loss_eff, na.rm = TRUE),
      wt_within = wt_loss_eff - wt_mean
    ) %>%
    ungroup()
}

.classify_arms <- function(arm_upper) {
  drug_class <- case_when(
    str_detect(arm_upper, "PLACEBO") ~ "Placebo",
    str_detect(arm_upper, "SEMAGLUTIDE|TIRZEPATIDE|SURVODUTIDE|PEMVIDUTIDE|LIRAGLUTIDE") ~ "Incretin",
    str_detect(arm_upper, "EFRUXIFERMIN|ALDAFERMIN|PEGOZAFERMIN|PEGBELFERMIN") ~ "FGF19/21 analogue",
    str_detect(arm_upper, "RESMETIRON|RESMETIROM") ~ "THR-beta agonist",
    str_detect(arm_upper, "LANAFIBRANOR|MSDC-0602K") ~ "PPAR agonist",
    TRUE ~ "Other"
  )
  mechanistic_group <- case_when(
    str_detect(arm_upper, paste(
      "SEMAGLUTIDE|TIRZEPATIDE|SURVODUTIDE|PEMVIDUTIDE|LIRAGLUTIDE",
      "EFRUXIFERMIN|ALDAFERMIN|PEGOZAFERMIN|PEGBELFERMIN",
      "RESMETIRON|RESMETIROM|LANAFIBRANOR|MSDC-0602K|DAPAGLIFLOZIN|PIOGLITAZONE|PXL065",
      sep = "|"
    )) ~ "Metabolic/Steatotic Modulator",
    str_detect(arm_upper, paste(
      "ION224|CENICRIVIROC|SELONSERTIB|FIRSOCOSTAT",
      "DENIFANSTAT|EMRICASAN|ARAMCHOL|ERVOGASTAT|CLESACOSTAT|APARARENONE",
      sep = "|"
    )) ~ "Direct Lipogenesis/Fibrogenesis Inhibitor",
    str_detect(arm_upper, "CILOFEXOR|TROPIFEXOR|VOLIXIBAT") ~ "Bile-Acid Modulator",
    str_detect(arm_upper, "VITAMIN E") ~ "Antioxidant",
    str_detect(arm_upper, "PLACEBO") ~ "Placebo",
    TRUE ~ "Other"
  )
  list(drug_class = drug_class, mechanistic_group = mechanistic_group)
}

.fix_stellar <- function(df) {
  df %>% mutate(study_id = case_when(
    str_detect(study_id, "STELLAR 3") ~ "STELLAR 3",
    str_detect(study_id, "STELLAR 4") ~ "STELLAR 4",
    TRUE ~ study_id
  ))
}

.desired_order <- c(
  "All",
  "drug_class | Incretin",
  "drug_class | FGF19/21 analogue",
  "drug_class | THR-beta agonist",
  "drug_class | PPAR agonist",
  "mechanistic_group | Metabolic/Steatotic Modulator",
  "mechanistic_group | Direct Lipogenesis/Fibrogenesis Inhibitor",
  "group_phase3 | Phase III",
  "group_f4 | F4 included",
  "group_risk | At-risk MASH only",
  "group_time | >=48 weeks"
)

prepare_mediation_data <- function(dataset, study_meta,
                                   outcome_col = "event1",
                                   n_col       = "n(PO1)") {
  stopifnot(
    is.data.frame(dataset),
    is.data.frame(study_meta),
    outcome_col %in% names(dataset),
    n_col       %in% names(dataset)
  )

  dat0 <- dataset %>%
    mutate(
      study_id  = clean_id(as.factor(STUDY)),
      arm_id    = as.factor(ARM),
      arm_upper = str_to_upper(as.character(ARM)),
      n         = as.numeric(.data[[n_col]]),
      wt_change = as.numeric(CFB_mean),
      wt_sd     = as.numeric(SD_CFB),
      events    = as.numeric(.data[[outcome_col]])
    ) %>%
    filter(
      !is.na(study_id), !is.na(arm_id), !is.na(n),
      !is.na(wt_change), !is.na(wt_sd), !is.na(events),
      is.finite(n), n > 0
    )

  if (nrow(dat0) == 0)
    stop("No valid rows after filtering. Check column names and data completeness.")

  classes <- .classify_arms(dat0$arm_upper)
  dat0    <- dat0 %>%
    mutate(drug_class = classes$drug_class, mechanistic_group = classes$mechanistic_group)

  study_meta <- study_meta %>% mutate(study_id = clean_id(study_id))
  dat0       <- .fix_stellar(dat0)
  study_meta <- .fix_stellar(study_meta)

  dat0 <- dat0 %>%
    left_join(study_meta, by = "study_id") %>%
    mutate(
      group_phase3 = if_else(phase3,       "Phase III",         "Non-Phase III"),
      group_f4     = if_else(f4_included,  "F4 included",       "No F4"),
      group_risk   = if_else(at_risk_only, "At-risk MASH only", "Mixed population"),
      group_time   = if_else(ge_48w,       ">=48 weeks",        "<48 weeks")
    )

  placebo <- dat0 %>%
    filter(drug_class == "Placebo") %>%
    group_by(study_id) %>%
    summarise(
      n_pl         = sum(n, na.rm = TRUE),
      wt_pl        = weighted.mean(wt_change, w = n),
      wt_sd_pl     = sqrt(sum(wt_sd^2 * n, na.rm = TRUE) / sum(n)),
      events_pl    = sum(events, na.rm = TRUE),
      nonevents_pl = sum(n - events, na.rm = TRUE),
      .groups = "drop"
    )

  if (n_distinct(placebo$study_id) == 0)
    stop("No placebo arms found. Check ARM values and drug_class classification.")

  active <- dat0 %>%
    filter(drug_class != "Placebo") %>%
    transmute(
      study_id, arm_id, drug_class, mechanistic_group,
      n_act = n, wt_act = wt_change, wt_sd_act = wt_sd,
      events_act = events, nonevents_act = n - events
    )

  contrasts <- active %>%
    inner_join(placebo, by = "study_id") %>%
    mutate(
      contrast_id = paste(study_id, arm_id, sep = "__"),
      wt_loss_eff = wt_pl - wt_act,
      var_wt_loss = (wt_sd_act^2 / n_act) + (wt_sd_pl^2 / n_pl),
      total_n     = n_act + n_pl
    ) %>%
    filter(
      is.finite(wt_loss_eff), is.finite(var_wt_loss),
      var_wt_loss > 0, n_act > 0, n_pl > 0
    ) %>%
    add_within_between()

  study_level_vars <- dat0 %>%
    dplyr::select(study_id, group_phase3, group_f4, group_risk, group_time) %>%
    distinct()

  contrasts <- contrasts %>% left_join(study_level_vars, by = "study_id")

  p_act <- contrasts$events_act / contrasts$n_act
  p_pl  <- contrasts$events_pl  / contrasts$n_pl

  contrasts <- contrasts %>%
    mutate(
      rd_resp     = p_act - p_pl,
      var_rd_resp = pmax(
        (pmax(p_act * (1 - p_act), 0) / n_act) + (pmax(p_pl * (1 - p_pl), 0) / n_pl),
        1e-8
      )
    ) %>%
    filter(is.finite(rd_resp), is.finite(var_rd_resp))

  contrasts_long <- contrasts %>%
    pivot_longer(
      cols      = c(drug_class, mechanistic_group, group_phase3, group_f4, group_risk, group_time),
      names_to  = "group_type",
      values_to = "group"
    ) %>%
    filter(!is.na(group), group != "")

  all_groups <- bind_rows(
    contrasts_long %>% distinct(group_type, group),
    tibble(group_type = "ALL", group = "All")
  )

  list(
    contrasts      = contrasts,
    contrasts_long = contrasts_long,
    all_groups     = all_groups,
    desired_order  = .desired_order
  )
}

prepare_mediation_data_1 <- function(dataset, study_meta) {
  prepare_mediation_data(dataset, study_meta, outcome_col = "event1", n_col = "n(PO1)")
}

prepare_mediation_data_2 <- function(dataset, study_meta) {
  prepare_mediation_data(dataset, study_meta, outcome_col = "event2", n_col = "n(PO2)")
}

analytic_re_posterior <- function(yi, vi, tau2 = 0,
                                   prior_mean = 0, prior_sd = PRIOR_B_SD,
                                   n_draws = 4000) {
  stopifnot(length(yi) == length(vi))
  tau2 <- if (is.na(tau2) || !is.finite(tau2) || tau2 < 0) 0 else tau2
  ok   <- is.finite(yi) & is.finite(vi) & vi > 0
  yi   <- yi[ok]; vi <- vi[ok]
  if (length(yi) == 0) {
    warning("analytic_re_posterior: no finite observations.")
    return(rep(NA_real_, n_draws))
  }
  mv         <- vi + tau2
  w          <- 1 / mv
  prior_prec <- 1 / prior_sd^2
  post_prec  <- sum(w) + prior_prec
  post_var   <- 1 / post_prec
  post_mean  <- post_var * (sum(w * yi) + prior_prec * prior_mean)
  rnorm(n_draws, post_mean, sqrt(post_var))
}

savage_dickey_bf10 <- function(post_draws, prior_sd = PRIOR_B_SD) {
  post_draws <- post_draws[is.finite(post_draws)]
  if (length(post_draws) < 10) return(NA_real_)
  prior_dens <- dnorm(0, 0, prior_sd)
  post_dens  <- tryCatch({
    d <- density(post_draws, n = 1024)
    approx(d$x, d$y, xout = 0)$y
  }, error = function(e) NA_real_)
  if (anyNA(c(prior_dens, post_dens)) || post_dens <= 0) return(NA_real_)
  prior_dens / post_dens
}

indirect_bf10 <- function(a_draws, b_draws,
                           prior_a_sd = PRIOR_B_SD, prior_b_sd = PRIOR_B_SD,
                           n_prior = 80000) {
  n_use          <- min(length(a_draws), length(b_draws))
  if (n_use < 10) return(NA_real_)
  indirect_post  <- a_draws[seq_len(n_use)] * b_draws[seq_len(n_use)]
  indirect_prior <- rnorm(n_prior, 0, prior_a_sd) * rnorm(n_prior, 0, prior_b_sd)
  kde_at_0 <- function(x) tryCatch({
    d <- density(x[is.finite(x)], n = 1024)
    approx(d$x, d$y, xout = 0)$y
  }, error = function(e) NA_real_)
  prior_dens <- kde_at_0(indirect_prior)
  post_dens  <- kde_at_0(indirect_post)
  if (anyNA(c(prior_dens, post_dens)) || prior_dens <= 0) return(NA_real_)
  prior_dens / post_dens
}

.empty_bpath <- function(k, n_studies, model_label, fit = NULL) {
  list(
    summary = tibble(
      estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
      BF10 = NA_real_, model = model_label, k = k, n_studies = n_studies,
      tau2 = NA_real_, I2 = NA_real_
    ),
    posterior = rep(NA_real_, 1),
    fit = fit
  )
}

fit_bpath_bayes <- function(df,
                             mods   = ~ wt_loss_eff + drug_class,
                             term   = "wt_within",
                             iter   = 2000, warmup = 1000,
                             chains = 4,   cores  = 4,
                             seed   = 42,  file   = NULL) {
  required_cols <- c("rd_resp", "var_rd_resp", "study_id", "wt_within", "wt_mean")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0)
    stop("fit_bpath_bayes: missing columns: ", paste(missing_cols, collapse = ", "))

  df        <- df %>% filter(is.finite(rd_resp), is.finite(var_rd_resp), var_rd_resp > 0)
  k_studies <- n_distinct(df$study_id)
  if (k_studies < 2) return(.empty_bpath(nrow(df), k_studies, "INSUFFICIENT_STUDIES"))

  rhs_terms <- attr(terms(mods), "term.labels")
  rhs_terms <- ifelse(rhs_terms == "wt_loss_eff", "wt_within + wt_mean", rhs_terms)
  formula   <- as.formula(paste0(
    "rd_resp | se(sqrt(var_rd_resp)) ~ ",
    paste(rhs_terms, collapse = " + "),
    " + (1 | study_id)"
  ))

  fit <- tryCatch(
    brm(
      formula   = formula,
      data      = df,
      prior     = .default_priors(),
      iter      = iter, warmup = warmup, chains = chains, cores = cores,
      seed      = seed, file = file,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.95, max_treedepth = 12),
      silent    = 2, refresh = 0
    ),
    error = function(e) { warning("brm() failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(.empty_bpath(nrow(df), k_studies, "BRMS_FAILED"))

  post_all  <- tryCatch(as_draws_df(fit), error = function(e) NULL)
  if (is.null(post_all)) return(.empty_bpath(nrow(df), k_studies, "DRAWS_EXTRACTION_FAILED", fit))

  brms_term <- paste0("b_", term)
  if (!brms_term %in% names(post_all)) {
    hits      <- grep(term, names(post_all), value = TRUE, fixed = TRUE)
    brms_term <- if (length(hits)) hits[1] else NA_character_
  }
  if (is.na(brms_term) || !brms_term %in% names(post_all)) {
    warning(sprintf("term '%s' not found in posterior.", term))
    return(.empty_bpath(nrow(df), k_studies, "TERM_NOT_FOUND", fit))
  }

  post_draws <- as.numeric(post_all[[brms_term]])
  tau2       <- tryCatch(
    mean(VarCorr(fit, summary = FALSE)$study_id$sd[, "Intercept"]^2, na.rm = TRUE),
    error = function(e) NA_real_
  )
  mean_vi <- mean(df$var_rd_resp, na.rm = TRUE)
  I2      <- if (is.finite(tau2) && is.finite(mean_vi) && (tau2 + mean_vi) > 0)
               tau2 / (tau2 + mean_vi) else NA_real_

  list(
    summary = tibble(
      estimate  = mean(post_draws, na.rm = TRUE),
      se        = sd(post_draws,   na.rm = TRUE),
      lower     = safe_quantile(post_draws, 0.025),
      upper     = safe_quantile(post_draws, 0.975),
      BF10      = savage_dickey_bf10(post_draws),
      model     = "BRMS", k = nrow(df), n_studies = k_studies, tau2 = tau2, I2 = I2
    ),
    posterior = post_draws,
    fit       = fit
  )
}

fit_meta_effect <- function(df, yi_name, vi_name, tau2 = 0, n_draws = 4000) {
  k_studies <- n_distinct(df$study_id)
  k_rows    <- nrow(df)
  empty_row <- function(msg) list(
    summary  = tibble(estimate = NA_real_, se = NA_real_, lower = NA_real_,
                      upper = NA_real_, BF10 = NA_real_, model = msg,
                      k = k_rows, n_studies = k_studies, tau2 = NA_real_, I2 = NA_real_),
    posterior = rep(NA_real_, n_draws)
  )
  if (k_studies < 2) return(empty_row("INSUFFICIENT_STUDIES"))

  yi         <- as.numeric(df[[yi_name]])
  vi         <- as.numeric(df[[vi_name]])
  post_draws <- tryCatch(
    analytic_re_posterior(yi, vi, tau2 = tau2, n_draws = n_draws),
    error = function(e) rep(NA_real_, n_draws)
  )
  if (all(is.na(post_draws))) return(empty_row("ANALYTIC_FAILED"))

  mean_vi <- mean(vi[is.finite(vi)], na.rm = TRUE)
  I2      <- if (is.finite(tau2) && is.finite(mean_vi) && (tau2 + mean_vi) > 0)
               tau2 / (tau2 + mean_vi) else NA_real_

  list(
    summary = tibble(
      estimate  = mean(post_draws, na.rm = TRUE),
      se        = sd(post_draws,   na.rm = TRUE),
      lower     = safe_quantile(post_draws, 0.025),
      upper     = safe_quantile(post_draws, 0.975),
      BF10      = savage_dickey_bf10(post_draws),
      model     = "Analytic-Bayes", k = k_rows, n_studies = k_studies, tau2 = tau2, I2 = I2
    ),
    posterior = post_draws
  )
}

compute_point_estimates_bayes <- function(contrasts, contrasts_long, all_groups,
                                          b_fit, tau2 = NULL, n_draws = 4000) {
  b_draws   <- as.numeric(unlist(b_fit$posterior))
  b_summary <- b_fit$summary
  if (all(is.na(b_draws))) b_draws <- rep(NA_real_, n_draws)
  if (is.null(tau2)) tau2 <- b_summary$tau2[1]
  if (is.na(tau2) || !is.finite(tau2)) tau2 <- 0

  results <- vector("list", nrow(all_groups))

  for (i in seq_len(nrow(all_groups))) {
    gtype  <- as.character(all_groups$group_type[i])
    grp    <- as.character(all_groups$group[i])
    df_sub <- if (gtype == "ALL") contrasts else
      contrasts_long %>% filter(group_type == gtype, group == grp)
    if (nrow(df_sub) == 0) next

    a_fit     <- tryCatch(fit_meta_effect(df_sub, "wt_loss_eff", "var_wt_loss", tau2, n_draws), error = function(e) NULL)
    total_fit <- tryCatch(fit_meta_effect(df_sub, "rd_resp",     "var_rd_resp",  tau2, n_draws), error = function(e) NULL)
    if (is.null(a_fit) || is.null(total_fit)) next

    a_draws     <- as.numeric(a_fit$posterior)
    total_draws <- as.numeric(total_fit$posterior)
    n_use       <- min(length(a_draws), length(b_draws), length(total_draws))
    a_s         <- a_draws[seq_len(n_use)]
    b_s         <- b_draws[seq_len(n_use)]
    tot_s       <- total_draws[seq_len(n_use)]

    indirect_draws <- a_s * b_s
    direct_draws   <- tot_s - indirect_draws

    a_smry   <- smry_draws(a_draws)
    tot_smry <- smry_draws(total_draws)
    ind_smry <- smry_draws(indirect_draws)
    dir_smry <- smry_draws(direct_draws)

    BF10_indirect <- tryCatch(indirect_bf10(a_draws, b_draws), error = function(e) NA_real_)
    group_label   <- if (gtype == "ALL") "All" else paste(gtype, grp, sep = " | ")
    b_est <- b_summary$estimate[1]
    b_lo  <- b_summary$lower[1]
    b_hi  <- b_summary$upper[1]

    results[[i]] <- tibble(
      group     = group_label,
      k         = nrow(df_sub),
      n_studies = n_distinct(df_sub$study_id),
      total_n   = sum(df_sub$total_n, na.rm = TRUE),
      a         = a_smry$m,  a_lower = a_smry$l,  a_upper = a_smry$u,  BF10_a = a_fit$summary$BF10[1],
      b         = b_est,     b_lower = b_lo,       b_upper = b_hi,      BF10_b = b_summary$BF10[1],
      total     = tot_smry$m, total_lower = tot_smry$l, total_upper = tot_smry$u, BF10_total = total_fit$summary$BF10[1],
      indirect  = ind_smry$m, indirect_lower = ind_smry$l, indirect_upper = ind_smry$u, BF10_indirect = BF10_indirect,
      direct    = dir_smry$m, direct_lower   = dir_smry$l, direct_upper   = dir_smry$u
    )
  }
  bind_rows(results)
}

make_forest_plot <- function(med, desired_order, b_global) {
  med <- med %>% filter(group %in% desired_order)
  if (nrow(med) == 0) { warning("make_forest_plot: no rows match desired_order."); return(NULL) }

  b_label <- b_global$b_lower[1]

  forest_data <- med %>%
    dplyr::select(group, total_n,
                  a, a_lower, a_upper,
                  indirect, indirect_lower, indirect_upper,
                  direct,   direct_lower,   direct_upper) %>%
    pivot_longer(cols = c(a, indirect, direct), names_to = "effect", values_to = "estimate") %>%
    mutate(
      lower = case_when(effect == "a" ~ a_lower, effect == "indirect" ~ indirect_lower, TRUE ~ direct_lower),
      upper = case_when(effect == "a" ~ a_upper, effect == "indirect" ~ indirect_upper, TRUE ~ direct_upper),
      effect     = factor(effect, levels = c("a", "indirect", "direct")),
      group      = factor(group, levels = rev(desired_order)),
      point_size = scales::rescale(log(total_n + 1), to = c(2, 8))
    )

  ggplot(forest_data,
      aes(x = estimate, y = group, xmin = lower, xmax = upper, size = point_size)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#56382D", linewidth = 0.8) +
    geom_errorbarh(height = 0.15, linewidth = 1.2, color = "#56382D") +
    geom_point(shape = 21, fill = "#E76D57", color = "black") +
    facet_wrap(~ effect, scales = "free_x", nrow = 1,
      labeller = as_labeller(c(
        a        = "a: Treatment \u2192 Weight loss",
        indirect = "Indirect effect (a\u00d7b)",
        direct   = "Direct effect"
      ))) +
    scale_size_continuous(name = "Total N", range = c(2, 8), labels = scales::comma) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"), panel.spacing = unit(1.5, "lines")) +
    labs(
      title    = "Mediation analysis: weight loss as mediator of histologic response",
      subtitle = paste0("b-path (posterior mean): ARD ", round(b_global$b[1], 4),
                        " [", round(b_lo, 4), ", ", round(b_global$b_upper[1], 4),
                        "] per 1 pp weight loss | 95% CrI"),
      x = "Absolute risk difference (95% CrI)", y = NULL
    )
}

make_bubble_plot <- function(contrasts, b_global) {
  bubble_data <- contrasts %>%
    mutate(weight_scaled = scales::rescale(sqrt(1 / var_rd_resp), to = c(2, 10)))

  ggplot(bubble_data,
      aes(x = wt_loss_eff, y = rd_resp, size = weight_scaled, color = drug_class)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
    geom_smooth(aes(color = drug_class), method = "lm", se = FALSE,
                linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
    scale_size_identity() +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(
      title    = "Trial-level association: weight loss vs histologic response",
      subtitle = paste0("b-path (posterior mean): ARD ", round(b_global$b[1], 4),
                        " per 1 pp weight loss | bubble size \u221d inverse-variance weight"),
      x     = "Weight loss vs placebo (percentage points)",
      y     = "Absolute risk difference of histologic response",
      color = "Drug class"
    )
}
