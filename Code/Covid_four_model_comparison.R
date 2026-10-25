# Four-model SARS-CoV-2 comparison
# --------------------------------
# Combines:
#   1. Lognormal, uncensored
#   2. WH-informed, uncensored
#   3. Lognormal, censored   [new]
#   4. WH-informed, censored
#
# This script does not refit anything. Run Covid_cens_lognormal.R first.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(loo)
  library(tibble)
})

OLD_DIR <- path.expand("~/Within-host-time-delay-framework/Results/Covid/Noncens")
CENS_DIR <- path.expand("~/Within-host-time-delay-framework/Results/Covid/Cens")
NEW_LOG_DIR <- file.path(CENS_DIR, "lognormal")
OUT_DIR <- file.path(CENS_DIR, "four_model_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Known fit paths from your existing scripts.
fits <- tibble::tribble(
  ~model,                    ~strain,     ~path,
  "Lognormal",              "SGTF",      file.path(OLD_DIR, "fit_park_lognormal_fast_sgtf.rds"),
  "Lognormal",              "non-SGTF",  file.path(OLD_DIR, "fit_park_lognormal_fast_nsgtf.rds"),
  "WH-informed",            "SGTF",      file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds"),
  "WH-informed",            "non-SGTF",  file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds"),
  "Lognormal censored",     "SGTF",      file.path(NEW_LOG_DIR, "fit_censored_lognormal_sgtf.rds"),
  "Lognormal censored",     "non-SGTF",  file.path(NEW_LOG_DIR, "fit_censored_lognormal_nsgtf.rds"),
  "WH-informed censored",   "SGTF",      file.path(CENS_DIR, "fit_date_coarsened_joint_obspair_sgtf.rds"),
  "WH-informed censored",   "non-SGTF",  file.path(CENS_DIR, "fit_date_coarsened_joint_obspair_nsgtf.rds")
)

missing <- fits %>% filter(!file.exists(path))
if (nrow(missing) > 0) {
  print(missing)
  stop("One or more fit files are missing. See table above.")
}

extract_gt <- function(path, model, strain) {
  fit <- readRDS(path)
  vars <- posterior::variables(fit$draws(inc_warmup = FALSE))
  need <- c("mean_gen", "sd_gen")
  miss <- setdiff(need, vars)
  if (length(miss)) stop(model, " / ", strain, " is missing: ", paste(miss, collapse = ", "))

  dr <- posterior::as_draws_df(
    fit$draws(variables = need, inc_warmup = FALSE)
  ) %>% as_tibble()

  tibble(
    model = model,
    strain = strain,
    mean_gen_median = median(dr$mean_gen),
    mean_gen_q025 = quantile(dr$mean_gen, 0.025),
    mean_gen_q975 = quantile(dr$mean_gen, 0.975),
    sd_gen_median = median(dr$sd_gen),
    sd_gen_q025 = quantile(dr$sd_gen, 0.025),
    sd_gen_q975 = quantile(dr$sd_gen, 0.975)
  )
}

gt_summary <- fits %>%
  rowwise() %>%
  do(extract_gt(.$path, .$model, .$strain)) %>%
  ungroup() %>%
  mutate(
    strain_label = recode(strain, SGTF = "Omicron", `non-SGTF` = "Delta"),
    model = factor(
      model,
      levels = c("Lognormal", "WH-informed", "Lognormal censored", "WH-informed censored")
    )
  )

write_csv(gt_summary, file.path(OUT_DIR, "generation_time_four_model_summary.csv"))
print(gt_summary)

# Pointwise information criteria from fit objects.
# IMPORTANT: only compare models that use the SAME observation likelihood:
#   uncensored pair: Lognormal vs WH-informed
#   censored pair:   Lognormal censored vs WH-informed censored
extract_ic <- function(path, model, strain) {
  fit <- readRDS(path)
  ll <- posterior::as_draws_matrix(fit$draws(variables = "log_lik", inc_warmup = FALSE))
  loo_fit <- loo::loo(ll)
  waic_fit <- loo::waic(ll)

  tibble(
    model = model,
    strain = strain,
    n = ncol(ll),
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_fit$estimates["elpd_loo", "SE"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    waic = waic_fit$estimates["waic", "Estimate"],
    p_waic = waic_fit$estimates["p_waic", "Estimate"]
  )
}

ic <- fits %>%
  rowwise() %>%
  do(extract_ic(.$path, .$model, .$strain)) %>%
  ungroup() %>%
  mutate(
    strain_label = recode(strain, SGTF = "Omicron", `non-SGTF` = "Delta"),
    likelihood = if_else(grepl("censored", model, fixed = TRUE), "Censored", "Uncensored")
  )

write_csv(ic, file.path(OUT_DIR, "information_criteria_four_models.csv"))
print(ic)

# Matched differences only. Positive delta_LOOIC means WH-informed has lower
# (better) LOOIC than lognormal under the same observation model.
matched <- ic %>%
  mutate(
    family = if_else(grepl("WH-informed", model, fixed = TRUE), "WH-informed", "Lognormal")
  ) %>%
  select(strain, strain_label, likelihood, family, looic, waic, p_loo, p_waic) %>%
  pivot_wider(
    names_from = family,
    values_from = c(looic, waic, p_loo, p_waic),
    names_sep = "__"
  ) %>%
  mutate(
    delta_looic_log_minus_wh = .data[["looic__Lognormal"]] - .data[["looic__WH-informed"]],
    delta_waic_log_minus_wh = .data[["waic__Lognormal"]] - .data[["waic__WH-informed"]]
  )

write_csv(matched, file.path(OUT_DIR, "matched_model_comparisons.csv"))
print(matched)

message("\nInterpret only within-likelihood comparisons:")
message("  Uncensored: Lognormal vs WH-informed")
message("  Censored:   Lognormal censored vs WH-informed censored")
message("Do NOT compare raw WAIC/LOOIC across censored vs uncensored likelihoods.")
