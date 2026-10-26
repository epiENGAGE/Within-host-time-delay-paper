library(loo)
library(posterior)
library(dplyr)
library(tibble)

# ============================================================
# TABLE 5: MPOX PAIRED LOO COMPARISONS + PARETO-k DIAGNOSTICS
# ============================================================

BASE_DIR <- "C:/Users/nj7786/Documents/MPOX"

# ------------------------------------------------------------
# Exact fitted-model files
# ------------------------------------------------------------

weak_fit_file <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair",
  "fit_mpox_score9_mechanistic_proportional_incubation_updatedV_obspair.rds"
)

yang_fit_file <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair_yangWHpriors_symrel",
  "fit_mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair.rds"
)

# Find the exact lognormal fitted object recursively
lognormal_fit_hits <- list.files(
  BASE_DIR,
  pattern = "^fit_mpox_lognormal_score9_fastconv\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(lognormal_fit_hits) == 0) {
  stop(
    "Could not find fit_mpox_lognormal_score9_fastconv.rds anywhere under:\n",
    BASE_DIR
  )
}

if (length(lognormal_fit_hits) > 1) {
  cat("\nMultiple lognormal fit objects found:\n")
  print(lognormal_fit_hits)
  stop("More than one exact lognormal fit object exists.")
}

lognormal_fit_file <- lognormal_fit_hits[1]

required_files <- c(
  weak_fit_file,
  yang_fit_file,
  lognormal_fit_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required fitted model file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("\nUsing fitted objects:\n")
cat("Weak-prior WHI:\n", weak_fit_file, "\n\n")
cat("Yang-informed WHI:\n", yang_fit_file, "\n\n")
cat("Lognormal:\n", lognormal_fit_file, "\n\n")


# ============================================================
# Load fitted models
# ============================================================

fit_weak <- readRDS(weak_fit_file)
fit_yang <- readRDS(yang_fit_file)
fit_lognormal <- readRDS(lognormal_fit_file)


# ============================================================
# Extract pointwise log-likelihood matrices
# ============================================================

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  
  vars <- posterior::variables(
    fit$draws(inc_warmup = FALSE)
  )
  
  if (!any(grepl(paste0("^", variable, "(\\[|$)"), vars))) {
    stop(
      "Fit does not contain generated quantity '",
      variable,
      "'."
    )
  }
  
  posterior::as_draws_matrix(
    fit$draws(
      variables = variable,
      inc_warmup = FALSE
    )
  )
}

loglik_weak <- extract_log_lik_matrix(
  fit_weak,
  "log_lik"
)

loglik_yang <- extract_log_lik_matrix(
  fit_yang,
  "log_lik"
)

loglik_lognormal <- extract_log_lik_matrix(
  fit_lognormal,
  "log_lik"
)


# ============================================================
# Check pointwise observation counts
# ============================================================

length_check <- tibble(
  model = c(
    "WH-informed, weak priors",
    "WH-informed, Yang-informed",
    "Lognormal"
  ),
  n_pointwise = c(
    ncol(loglik_weak),
    ncol(loglik_yang),
    ncol(loglik_lognormal)
  )
)

cat("\n====================================================\n")
cat("POINTWISE OBSERVATION COUNTS\n")
cat("====================================================\n\n")

print(length_check)

if (length(unique(length_check$n_pointwise)) != 1) {
  stop(
    "The three models do not have the same number of pointwise observations. ",
    "Paired LOO comparisons require identical observation-level likelihood contributions."
  )
}

if (unique(length_check$n_pointwise) != 34) {
  warning(
    "Expected 34 high-confidence mpox transmission pairs, but found ",
    unique(length_check$n_pointwise),
    " pointwise observations."
  )
}


# ============================================================
# Compute PSIS-LOO
# ============================================================

loo_weak <- loo::loo(loglik_weak)
loo_yang <- loo::loo(loglik_yang)
loo_lognormal <- loo::loo(loglik_lognormal)


# ============================================================
# Helper: pointwise ELPD
# ============================================================

pointwise_elpd <- function(loo_obj) {
  
  x <- loo_obj$pointwise[, "elpd_loo"]
  
  if (is.null(x)) {
    stop("Could not extract pointwise elpd_loo values.")
  }
  
  as.numeric(x)
}


# ============================================================
# Helper: paired LOO comparison
# ============================================================

paired_loo_comparison <- function(
    loo_a,
    loo_b,
    model_a,
    model_b
) {
  
  a <- pointwise_elpd(loo_a)
  b <- pointwise_elpd(loo_b)
  
  if (length(a) != length(b)) {
    stop(
      "Pointwise LOO lengths differ for ",
      model_a,
      " vs ",
      model_b,
      ": ",
      length(a),
      " vs ",
      length(b)
    )
  }
  
  # Paired observation-level ELPD differences
  d <- a - b
  
  n <- length(d)
  
  # ELPD_A - ELPD_B
  delta_elpd <- sum(d)
  
  # Reviewer-requested SE of the paired difference
  se_delta_elpd <- sqrt(
    n * stats::var(d)
  )
  
  # LOOIC = -2 * ELPD
  #
  # Therefore:
  #
  # LOOIC_B - LOOIC_A
  # = 2 * (ELPD_A - ELPD_B)
  delta_looic_b_minus_a <- 2 * delta_elpd
  se_delta_looic <- 2 * se_delta_elpd
  
  tibble(
    model_A = model_a,
    model_B = model_b,
    n_pointwise = n,
    
    elpd_A =
      unname(
        loo_a$estimates[
          "elpd_loo",
          "Estimate"
        ]
      ),
    
    elpd_B =
      unname(
        loo_b$estimates[
          "elpd_loo",
          "Estimate"
        ]
      ),
    
    looic_A =
      unname(
        loo_a$estimates[
          "looic",
          "Estimate"
        ]
      ),
    
    looic_B =
      unname(
        loo_b$estimates[
          "looic",
          "Estimate"
        ]
      ),
    
    delta_elpd_A_minus_B =
      delta_elpd,
    
    se_delta_elpd =
      se_delta_elpd,
    
    delta_looic_B_minus_A =
      delta_looic_b_minus_a,
    
    se_delta_looic =
      se_delta_looic
  )
}


# ============================================================
# All three Table 5 pairwise comparisons
# ============================================================

table5_pairwise <- bind_rows(
  
  # Weak-prior WHI versus lognormal
  paired_loo_comparison(
    loo_a = loo_weak,
    loo_b = loo_lognormal,
    model_a = "WH-informed, weak priors",
    model_b = "Lognormal"
  ),
  
  # Yang-informed WHI versus lognormal
  paired_loo_comparison(
    loo_a = loo_yang,
    loo_b = loo_lognormal,
    model_a = "WH-informed, Yang-informed",
    model_b = "Lognormal"
  ),
  
  # Weak-prior WHI versus Yang-informed WHI
  paired_loo_comparison(
    loo_a = loo_weak,
    loo_b = loo_yang,
    model_a = "WH-informed, weak priors",
    model_b = "WH-informed, Yang-informed"
  )
)


# ============================================================
# Pareto-k diagnostics
# ============================================================

pareto_summary <- function(
    loo_obj,
    model
) {
  
  k <- loo::pareto_k_values(
    loo_obj
  )
  
  tibble(
    model = model,
    n_obs = length(k),
    
    median_pareto_k =
      median(
        k,
        na.rm = TRUE
      ),
    
    max_pareto_k =
      max(
        k,
        na.rm = TRUE
      ),
    
    n_k_gt_0_5 =
      sum(
        k > 0.5,
        na.rm = TRUE
      ),
    
    n_k_gt_0_7 =
      sum(
        k > 0.7,
        na.rm = TRUE
      ),
    
    n_k_gt_1_0 =
      sum(
        k > 1.0,
        na.rm = TRUE
      )
  )
}


table5_pareto <- bind_rows(
  
  pareto_summary(
    loo_lognormal,
    "Lognormal"
  ),
  
  pareto_summary(
    loo_weak,
    "WH-informed, weak priors"
  ),
  
  pareto_summary(
    loo_yang,
    "WH-informed, Yang-informed"
  )
)


# ============================================================
# Reviewer-ready pairwise comparison table
# ============================================================

table5_manuscript <- table5_pairwise %>%
  transmute(
    
    Comparison =
      paste(
        model_A,
        "vs",
        model_B
      ),
    
    `LOOIC model A` =
      round(
        looic_A,
        2
      ),
    
    `Delta LOOIC (B - A)` =
      round(
        delta_looic_B_minus_A,
        2
      ),
    
    `SE of Delta LOOIC` =
      round(
        se_delta_looic,
        2
      ),
    
    `Delta ELPD (A - B)` =
      round(
        delta_elpd_A_minus_B,
        2
      ),
    
    `SE of Delta ELPD` =
      round(
        se_delta_elpd,
        2
      )
  )


# ============================================================
# Reviewer-ready Pareto table
# ============================================================

table5_pareto_manuscript <- table5_pareto %>%
  transmute(
    
    Model =
      model,
    
    `Median Pareto k` =
      round(
        median_pareto_k,
        2
      ),
    
    `Max Pareto k` =
      round(
        max_pareto_k,
        2
      ),
    
    `n(k > 0.5)` =
      n_k_gt_0_5,
    
    `n(k > 0.7)` =
      n_k_gt_0_7,
    
    `n(k > 1)` =
      n_k_gt_1_0
  )


# ============================================================
# Print full outputs
# ============================================================

cat("\n\n====================================================\n")
cat("TABLE 5: PAIRED PREDICTIVE COMPARISONS\n")
cat("====================================================\n\n")

cat(
  "For each row, positive Delta ELPD means model A predicts better.\n"
)

cat(
  "Positive Delta LOOIC (B - A) means model A has lower LOOIC and predicts better.\n\n"
)

print(
  table5_manuscript,
  n = Inf,
  width = Inf
)


cat("\n\n====================================================\n")
cat("TABLE 5: PARETO-k DIAGNOSTICS\n")
cat("====================================================\n\n")

print(
  table5_pareto_manuscript,
  n = Inf,
  width = Inf
)


# ============================================================
# loo_compare cross-checks
# ============================================================

cat("\n\n====================================================\n")
cat("loo_compare CROSS-CHECKS\n")
cat("====================================================\n")


cat("\nWeak-prior WHI vs lognormal\n")

print(
  loo::loo_compare(
    list(
      WH_weak =
        loo_weak,
      
      Lognormal =
        loo_lognormal
    )
  )
)


cat("\nYang-informed WHI vs lognormal\n")

print(
  loo::loo_compare(
    list(
      WH_Yang =
        loo_yang,
      
      Lognormal =
        loo_lognormal
    )
  )
)


cat("\nWeak-prior WHI vs Yang-informed WHI\n")

print(
  loo::loo_compare(
    list(
      WH_weak =
        loo_weak,
      
      WH_Yang =
        loo_yang
    )
  )
)


# ============================================================
# Sanity checks
# ============================================================

stopifnot(
  
  all(
    is.finite(
      table5_pairwise$
        delta_elpd_A_minus_B
    )
  ),
  
  all(
    is.finite(
      table5_pairwise$
        se_delta_elpd
    )
  ),
  
  all(
    is.finite(
      table5_pairwise$
        delta_looic_B_minus_A
    )
  ),
  
  all(
    is.finite(
      table5_pairwise$
        se_delta_looic
    )
  ),
  
  all(
    table5_pareto$n_obs ==
      unique(
        length_check$n_pointwise
      )
  )
)


# ============================================================
# Save reviewer-ready outputs
# ============================================================

OUT_DIR_REVIEWER <- file.path(
  BASE_DIR,
  "table5_reviewer_LOO"
)

dir.create(
  OUT_DIR_REVIEWER,
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  table5_manuscript,
  file.path(
    OUT_DIR_REVIEWER,
    "table5_pairwise_LOO_comparisons.csv"
  ),
  row.names = FALSE
)

write.csv(
  table5_pareto_manuscript,
  file.path(
    OUT_DIR_REVIEWER,
    "table5_pareto_k_diagnostics.csv"
  ),
  row.names = FALSE
)


cat("\n\n====================================================\n")
cat("SUCCESS\n")
cat("Table 5 paired LOO comparisons and Pareto-k diagnostics completed.\n")
cat("Outputs written to:\n")
cat(OUT_DIR_REVIEWER, "\n")

waic_lognormal <- loo::waic(loglik_lognormal)
waic_weak      <- loo::waic(loglik_weak)
waic_yang      <- loo::waic(loglik_yang)

rbind(
  Lognormal = c(
    WAIC = waic_lognormal$estimates["waic", "Estimate"],
    pWAIC = waic_lognormal$estimates["p_waic", "Estimate"]
  ),
  `WHI, weak priors` = c(
    WAIC = waic_weak$estimates["waic", "Estimate"],
    pWAIC = waic_weak$estimates["p_waic", "Estimate"]
  ),
  `WHI, Yang-informed robust prior` = c(
    WAIC = waic_yang$estimates["waic", "Estimate"],
    pWAIC = waic_yang$estimates["p_waic", "Estimate"]
  )
)
cat("====================================================\n")