# COVID Table 3 reviewer-requested predictive comparisons
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(dplyr)
  library(tibble)
  library(readr)
})

BASE_DIR <- "~/Within-host-time-delay-paper/Results/Covid"
NONCENS_DIR <- file.path(BASE_DIR, "Noncens")
CENS_DIR <- file.path(BASE_DIR, "Cens")
CENS_LN_DIR <- file.path(CENS_DIR, "lognormal")
OUT_DIR <- file.path(BASE_DIR, "Reviewer_predictive_comparisons")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

load_loo <- function(path) {
  if (!file.exists(path)) stop("Missing LOO file: ", path)
  x <- readRDS(path)
  if (!inherits(x, "loo")) stop("Not a loo object: ", path)
  x
}

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  posterior::as_draws_matrix(fit$draws(variables = variable, inc_warmup = FALSE))
}

pointwise_elpd <- function(loo_obj) {
  if (is.null(loo_obj$pointwise) || !"elpd_loo" %in% colnames(loo_obj$pointwise))
    stop("LOO object lacks pointwise elpd_loo values.")
  as.numeric(loo_obj$pointwise[, "elpd_loo"])
}

paired_loo_comparison <- function(loo_a, loo_b, strain, observation_model,
                                  model_a = "WH-informed", model_b = "Lognormal") {
  a <- pointwise_elpd(loo_a)
  b <- pointwise_elpd(loo_b)
  if (length(a) != length(b)) {
    stop("Pointwise LOO lengths differ for ", strain, " / ", observation_model,
         ": ", length(a), " vs ", length(b))
  }
  d <- a - b
  n <- length(d)
  delta_elpd <- sum(d)
  se_delta_elpd <- sqrt(n * stats::var(d))
  delta_looic <- 2 * delta_elpd
  se_delta_looic <- 2 * se_delta_elpd
  
  tibble(
    strain = strain,
    observation_model = observation_model,
    model_A = model_a,
    model_B = model_b,
    n_pointwise = n,
    elpd_A = unname(loo_a$estimates["elpd_loo", "Estimate"]),
    elpd_B = unname(loo_b$estimates["elpd_loo", "Estimate"]),
    delta_elpd_A_minus_B = delta_elpd,
    se_delta_elpd = se_delta_elpd,
    looic_A = unname(loo_a$estimates["looic", "Estimate"]),
    looic_B = unname(loo_b$estimates["looic", "Estimate"]),
    delta_looic_B_minus_A = delta_looic,
    se_delta_looic = se_delta_looic
  )
}

pareto_summary <- function(loo_obj, strain, observation_model, model) {
  k <- loo::pareto_k_values(loo_obj)
  tibble(
    strain = strain,
    observation_model = observation_model,
    model = model,
    n_obs = length(k),
    max_pareto_k = max(k, na.rm = TRUE),
    n_k_gt_0_5 = sum(k > 0.5, na.rm = TRUE),
    n_k_gt_0_7 = sum(k > 0.7, na.rm = TRUE),
    n_k_gt_1_0 = sum(k > 1.0, na.rm = TRUE)
  )
}

# Uncensored WH-informed
loo_wh_uncens_omicron <- load_loo(file.path(
  NONCENS_DIR,
  "mechanistic_proportional_incubation_updatedV_obspair_sgtf_loo.rds"
))
loo_wh_uncens_delta <- load_loo(file.path(
  NONCENS_DIR,
  "mechanistic_proportional_incubation_updatedV_obspair_nsgtf_loo.rds"
))

# Uncensored lognormal: reconstruct LOO from saved fits
fit_ln_uncens_omicron_path <- file.path(NONCENS_DIR, "fit_park_lognormal_fast_sgtf.rds")
fit_ln_uncens_delta_path <- file.path(NONCENS_DIR, "fit_park_lognormal_fast_nsgtf.rds")
if (!file.exists(fit_ln_uncens_omicron_path)) stop("Missing: ", fit_ln_uncens_omicron_path)
if (!file.exists(fit_ln_uncens_delta_path)) stop("Missing: ", fit_ln_uncens_delta_path)

fit_ln_uncens_omicron <- readRDS(fit_ln_uncens_omicron_path)
fit_ln_uncens_delta <- readRDS(fit_ln_uncens_delta_path)

vars_omicron <- posterior::variables(fit_ln_uncens_omicron$draws(inc_warmup = FALSE))
vars_delta <- posterior::variables(fit_ln_uncens_delta$draws(inc_warmup = FALSE))
if (!any(grepl("^log_lik(\\[|$)", vars_omicron)))
  stop("Omicron lognormal fit has no log_lik generated quantity.")
if (!any(grepl("^log_lik(\\[|$)", vars_delta)))
  stop("Delta lognormal fit has no log_lik generated quantity.")

loglik_ln_uncens_omicron <- extract_log_lik_matrix(fit_ln_uncens_omicron, "log_lik")
loglik_ln_uncens_delta <- extract_log_lik_matrix(fit_ln_uncens_delta, "log_lik")
loo_ln_uncens_omicron <- loo::loo(loglik_ln_uncens_omicron)
loo_ln_uncens_delta <- loo::loo(loglik_ln_uncens_delta)

saveRDS(loo_ln_uncens_omicron, file.path(NONCENS_DIR, "park_lognormal_fast_sgtf_loo.rds"))
saveRDS(loo_ln_uncens_delta, file.path(NONCENS_DIR, "park_lognormal_fast_nsgtf_loo.rds"))

# Censored WH-informed
loo_wh_cens_omicron <- load_loo(file.path(CENS_DIR, "sgtf_loo.rds"))
loo_wh_cens_delta <- load_loo(file.path(CENS_DIR, "nsgtf_loo.rds"))

# Censored lognormal
loo_ln_cens_omicron <- load_loo(file.path(CENS_LN_DIR, "sgtf_loo.rds"))
loo_ln_cens_delta <- load_loo(file.path(CENS_LN_DIR, "nsgtf_loo.rds"))

# Check paired observation counts
length_check <- tibble(
  strain = c("Omicron", "Delta", "Omicron", "Delta"),
  observation_model = c("Uncensored", "Uncensored", "Censored", "Censored"),
  n_WH = c(
    length(pointwise_elpd(loo_wh_uncens_omicron)),
    length(pointwise_elpd(loo_wh_uncens_delta)),
    length(pointwise_elpd(loo_wh_cens_omicron)),
    length(pointwise_elpd(loo_wh_cens_delta))
  ),
  n_lognormal = c(
    length(pointwise_elpd(loo_ln_uncens_omicron)),
    length(pointwise_elpd(loo_ln_uncens_delta)),
    length(pointwise_elpd(loo_ln_cens_omicron)),
    length(pointwise_elpd(loo_ln_cens_delta))
  )
)
print(length_check)
if (any(length_check$n_WH != length_check$n_lognormal))
  stop("At least one comparison does not have aligned pointwise observation counts.")

table3_pairwise <- bind_rows(
  paired_loo_comparison(loo_wh_uncens_omicron, loo_ln_uncens_omicron, "Omicron", "Uncensored"),
  paired_loo_comparison(loo_wh_uncens_delta, loo_ln_uncens_delta, "Delta", "Uncensored"),
  paired_loo_comparison(loo_wh_cens_omicron, loo_ln_cens_omicron, "Omicron", "Censored"),
  paired_loo_comparison(loo_wh_cens_delta, loo_ln_cens_delta, "Delta", "Censored")
)

table3_pareto <- bind_rows(
  pareto_summary(loo_wh_uncens_omicron, "Omicron", "Uncensored", "WH-informed"),
  pareto_summary(loo_ln_uncens_omicron, "Omicron", "Uncensored", "Lognormal"),
  pareto_summary(loo_wh_uncens_delta, "Delta", "Uncensored", "WH-informed"),
  pareto_summary(loo_ln_uncens_delta, "Delta", "Uncensored", "Lognormal"),
  pareto_summary(loo_wh_cens_omicron, "Omicron", "Censored", "WH-informed"),
  pareto_summary(loo_ln_cens_omicron, "Omicron", "Censored", "Lognormal"),
  pareto_summary(loo_wh_cens_delta, "Delta", "Censored", "WH-informed"),
  pareto_summary(loo_ln_cens_delta, "Delta", "Censored", "Lognormal")
)

table3_manuscript <- table3_pairwise %>%
  transmute(
    Strain = strain,
    Observation = observation_model,
    `WH-informed LOOIC` = round(looic_A, 1),
    `Lognormal LOOIC` = round(looic_B, 1),
    `Delta LOOIC (Lognormal - WH)` = round(delta_looic_B_minus_A, 1),
    `SE of Delta LOOIC` = round(se_delta_looic, 1)
  )

table3_pareto_manuscript <- table3_pareto %>%
  transmute(
    Strain = strain,
    Observation = observation_model,
    Model = model,
    `Max Pareto k` = round(max_pareto_k, 3),
    `n(k > 0.7)` = n_k_gt_0_7,
    `n(k > 1)` = n_k_gt_1_0
  )

cat("\n\n=== TABLE 3: paired predictive comparison ===\n")
cat("Positive Delta LOOIC means lower LOOIC for WH-informed.\n\n")
print(table3_manuscript, n = Inf)

cat("\n\n=== TABLE 3: Pareto-k diagnostics ===\n\n")
print(table3_pareto_manuscript, n = Inf)

cat("\n\n=== loo_compare cross-checks ===\n")
cat("\nOmicron / uncensored\n")
print(loo::loo_compare(list(WH_informed = loo_wh_uncens_omicron, Lognormal = loo_ln_uncens_omicron)))
cat("\nDelta / uncensored\n")
print(loo::loo_compare(list(WH_informed = loo_wh_uncens_delta, Lognormal = loo_ln_uncens_delta)))
cat("\nOmicron / censored\n")
print(loo::loo_compare(list(WH_informed = loo_wh_cens_omicron, Lognormal = loo_ln_cens_omicron)))
cat("\nDelta / censored\n")
print(loo::loo_compare(list(WH_informed = loo_wh_cens_delta, Lognormal = loo_ln_cens_delta)))

write_csv(table3_pairwise, file.path(OUT_DIR, "table3_pairwise_loo_full_precision.csv"))
write_csv(table3_manuscript, file.path(OUT_DIR, "table3_pairwise_loo_manuscript.csv"))
write_csv(table3_pareto, file.path(OUT_DIR, "table3_pareto_k_full_precision.csv"))
write_csv(table3_pareto_manuscript, file.path(OUT_DIR, "table3_pareto_k_manuscript.csv"))
write_csv(length_check, file.path(OUT_DIR, "table3_pointwise_length_check.csv"))

cat("\nSaved outputs to:\n", normalizePath(OUT_DIR), "\n")

library(loo)
library(dplyr)
library(tibble)

# ============================================================
# Helpers
# ============================================================

pointwise_elpd <- function(loo_obj) {
  x <- loo_obj$pointwise[, "elpd_loo"]
  
  if (is.null(x)) {
    stop("Could not extract pointwise elpd_loo values.")
  }
  
  as.numeric(x)
}

paired_loo_comparison <- function(
    loo_a,
    loo_b,
    strain,
    observation_model,
    model_a = "WH-informed",
    model_b = "Lognormal"
) {
  
  a <- pointwise_elpd(loo_a)
  b <- pointwise_elpd(loo_b)
  
  if (length(a) != length(b)) {
    stop(
      "Pointwise LOO lengths differ for ",
      strain, " / ", observation_model,
      ": ", length(a), " vs ", length(b)
    )
  }
  
  # Paired pointwise differences in ELPD
  d <- a - b
  n <- length(d)
  
  # ELPD_A - ELPD_B
  delta_elpd <- sum(d)
  
  # SE based on the paired pointwise differences
  se_delta_elpd <- sqrt(n * stats::var(d))
  
  # LOOIC = -2 * ELPD
  #
  # Therefore:
  # LOOIC_B - LOOIC_A
  # = 2 * (ELPD_A - ELPD_B)
  delta_looic <- 2 * delta_elpd
  se_delta_looic <- 2 * se_delta_elpd
  
  tibble(
    strain = strain,
    observation_model = observation_model,
    model_A = model_a,
    model_B = model_b,
    n_pointwise = n,
    
    elpd_A = unname(loo_a$estimates["elpd_loo", "Estimate"]),
    elpd_B = unname(loo_b$estimates["elpd_loo", "Estimate"]),
    
    delta_elpd_A_minus_B = delta_elpd,
    se_delta_elpd = se_delta_elpd,
    
    looic_A = unname(loo_a$estimates["looic", "Estimate"]),
    looic_B = unname(loo_b$estimates["looic", "Estimate"]),
    
    delta_looic_B_minus_A = delta_looic,
    se_delta_looic = se_delta_looic
  )
}

pareto_summary <- function(
    loo_obj,
    strain,
    observation_model,
    model
) {
  
  k <- loo::pareto_k_values(loo_obj)
  
  tibble(
    strain = strain,
    observation_model = observation_model,
    model = model,
    n_obs = length(k),
    
    median_pareto_k = median(k, na.rm = TRUE),
    max_pareto_k = max(k, na.rm = TRUE),
    
    n_k_gt_0_5 = sum(k > 0.5, na.rm = TRUE),
    n_k_gt_0_7 = sum(k > 0.7, na.rm = TRUE),
    n_k_gt_1_0 = sum(k > 1.0, na.rm = TRUE)
  )
}

# ============================================================
# TABLE 3: COVID / SARS-CoV-2 ONLY
# Assumes these objects already exist:
#
# loo_wh_uncens_omicron
# loo_ln_uncens_omicron
# loo_wh_uncens_delta
# loo_ln_uncens_delta
# loo_wh_cens_omicron
# loo_ln_cens_omicron
# loo_wh_cens_delta
# loo_ln_cens_delta
# ============================================================


# ------------------------------------------------------------
# Check paired observation counts
# ------------------------------------------------------------

length_check <- tibble(
  strain = c(
    "Omicron",
    "Delta",
    "Omicron",
    "Delta"
  ),
  observation_model = c(
    "Uncensored",
    "Uncensored",
    "Censored",
    "Censored"
  ),
  n_WH = c(
    length(pointwise_elpd(loo_wh_uncens_omicron)),
    length(pointwise_elpd(loo_wh_uncens_delta)),
    length(pointwise_elpd(loo_wh_cens_omicron)),
    length(pointwise_elpd(loo_wh_cens_delta))
  ),
  n_lognormal = c(
    length(pointwise_elpd(loo_ln_uncens_omicron)),
    length(pointwise_elpd(loo_ln_uncens_delta)),
    length(pointwise_elpd(loo_ln_cens_omicron)),
    length(pointwise_elpd(loo_ln_cens_delta))
  )
)

print(length_check)

if (any(length_check$n_WH != length_check$n_lognormal)) {
  stop("At least one comparison does not have matching pointwise observation counts.")
}


# ------------------------------------------------------------
# Paired predictive comparisons
# ------------------------------------------------------------

table3_pairwise <- bind_rows(
  
  paired_loo_comparison(
    loo_wh_uncens_omicron,
    loo_ln_uncens_omicron,
    strain = "Omicron",
    observation_model = "Uncensored"
  ),
  
  paired_loo_comparison(
    loo_wh_uncens_delta,
    loo_ln_uncens_delta,
    strain = "Delta",
    observation_model = "Uncensored"
  ),
  
  paired_loo_comparison(
    loo_wh_cens_omicron,
    loo_ln_cens_omicron,
    strain = "Omicron",
    observation_model = "Censored"
  ),
  
  paired_loo_comparison(
    loo_wh_cens_delta,
    loo_ln_cens_delta,
    strain = "Delta",
    observation_model = "Censored"
  )
)


# ------------------------------------------------------------
# Pareto-k diagnostics
# ------------------------------------------------------------

table3_pareto <- bind_rows(
  
  pareto_summary(
    loo_wh_uncens_omicron,
    "Omicron",
    "Uncensored",
    "WH-informed"
  ),
  
  pareto_summary(
    loo_ln_uncens_omicron,
    "Omicron",
    "Uncensored",
    "Lognormal"
  ),
  
  pareto_summary(
    loo_wh_uncens_delta,
    "Delta",
    "Uncensored",
    "WH-informed"
  ),
  
  pareto_summary(
    loo_ln_uncens_delta,
    "Delta",
    "Uncensored",
    "Lognormal"
  ),
  
  pareto_summary(
    loo_wh_cens_omicron,
    "Omicron",
    "Censored",
    "WH-informed"
  ),
  
  pareto_summary(
    loo_ln_cens_omicron,
    "Omicron",
    "Censored",
    "Lognormal"
  ),
  
  pareto_summary(
    loo_wh_cens_delta,
    "Delta",
    "Censored",
    "WH-informed"
  ),
  
  pareto_summary(
    loo_ln_cens_delta,
    "Delta",
    "Censored",
    "Lognormal"
  )
)


# ------------------------------------------------------------
# Manuscript-ready predictive comparison table
# ------------------------------------------------------------

table3_manuscript <- table3_pairwise %>%
  transmute(
    Strain = strain,
    Observation = observation_model,
    
    `WH-informed LOOIC` =
      round(looic_A, 1),
    
    `Lognormal LOOIC` =
      round(looic_B, 1),
    
    `Delta ELPD (WH - Lognormal)` =
      round(delta_elpd_A_minus_B, 1),
    
    `SE of Delta ELPD` =
      round(se_delta_elpd, 1),
    
    `Delta LOOIC (Lognormal - WH)` =
      round(delta_looic_B_minus_A, 1),
    
    `SE of Delta LOOIC` =
      round(se_delta_looic, 1)
  )


# ------------------------------------------------------------
# Manuscript-ready Pareto-k diagnostic table
# ------------------------------------------------------------

table3_pareto_manuscript <- table3_pareto %>%
  transmute(
    Strain = strain,
    Observation = observation_model,
    Model = model,
    
    `Median Pareto k` =
      round(median_pareto_k, 2),
    
    `Max Pareto k` =
      round(max_pareto_k, 2),
    
    `n(k > 0.5)` =
      n_k_gt_0_5,
    
    `n(k > 0.7)` =
      n_k_gt_0_7,
    
    `n(k > 1)` =
      n_k_gt_1_0
  )


# ------------------------------------------------------------
# Print reviewer-ready outputs
# ------------------------------------------------------------

cat("\n\n====================================================\n")
cat("TABLE 3: PAIRED PREDICTIVE COMPARISON\n")
cat("====================================================\n\n")

cat("Positive Delta ELPD means better predictive performance for WH-informed.\n")
cat("Positive Delta LOOIC means lower LOOIC for WH-informed.\n\n")

print(table3_manuscript, n = Inf)


cat("\n\n====================================================\n")
cat("TABLE 3: PARETO-k DIAGNOSTICS\n")
cat("====================================================\n\n")

print(table3_pareto_manuscript, n = Inf)


# ------------------------------------------------------------
# loo_compare cross-checks
# ------------------------------------------------------------

cat("\n\n====================================================\n")
cat("loo_compare CROSS-CHECKS\n")
cat("====================================================\n")

cat("\nOmicron / uncensored\n")
print(
  loo::loo_compare(
    list(
      WH_informed = loo_wh_uncens_omicron,
      Lognormal = loo_ln_uncens_omicron
    )
  )
)

cat("\nDelta / uncensored\n")
print(
  loo::loo_compare(
    list(
      WH_informed = loo_wh_uncens_delta,
      Lognormal = loo_ln_uncens_delta
    )
  )
)

cat("\nOmicron / censored\n")
print(
  loo::loo_compare(
    list(
      WH_informed = loo_wh_cens_omicron,
      Lognormal = loo_ln_cens_omicron
    )
  )
)

cat("\nDelta / censored\n")
print(
  loo::loo_compare(
    list(
      WH_informed = loo_wh_cens_delta,
      Lognormal = loo_ln_cens_delta
    )
  )
)


# ------------------------------------------------------------
# Final sanity checks
# ------------------------------------------------------------

stopifnot(
  all(table3_pareto$n_obs > 0),
  
  all(is.finite(table3_pairwise$delta_elpd_A_minus_B)),
  all(is.finite(table3_pairwise$se_delta_elpd)),
  
  all(is.finite(table3_pairwise$delta_looic_B_minus_A)),
  all(is.finite(table3_pairwise$se_delta_looic))
)

cat("\n\n====================================================\n")
cat("SUCCESS\n")
cat("Table 3 paired LOO comparisons and Pareto-k diagnostics completed.\n")
cat("====================================================\n")