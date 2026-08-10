# Simulation 3
required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "tidyr", "ggplot2", "readr",
  "stringr", "purrr", "loo", "tibble", "future.apply", "future"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script. For cmdstanr, use install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))."
  )
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(loo)
  library(tibble)
  library(future.apply)
  library(future)
})

# User controls
OUT_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STAN_FILE <- "~/Within-host-time-delay-framework/Code/sim3.stan"

# Main design. Start small for a smoke test, then scale up.
N_REP <- 50
N_PAIRS <- 200
N_TEST <- 5000

CHAINS <- 4
PARALLEL_CHAINS <- 4
ITER_WARMUP <- 300
ITER_SAMPLING <- 300
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15
REFRESH <- 100
SEED <- 52001
SAVE_FITS <- FALSE
FORCE_RECOMPILE <- FALSE

PARALLELIZE_FITS <- FALSE
N_WORKERS <- 1

WRITE_JOB_LEVEL_RDS <- TRUE

GRID_BY_SI <- 0.5
GRID_BY_LU <- 0.5
GRID_BY_IP <- 0.10
N_Z_QUAD <- 7
Z_MAX_SD <- 3.5

SI_MIN <- -10
SI_MAX <- 20
LU_MAX <- 20
GT_MAX <- 2 * LU_MAX
IP_MIN <- 0.05
IP_MAX <- 30

TRUE_R_GROWTH <- 0.0

TRUE_PARAMETERS <- list(
  log_lambda1 = log(0.16),
  log_lambda_I = log(0.16),
  log_r = log(0.27),
  log_d = log(0.62),
  m = 3.5,
  log_kappa = log(1.05),
  log_sigma_z = log(0.35)
)

PARAM_VARS <- names(TRUE_PARAMETERS)
FIT_SUMMARY_VARS <- c(PARAM_VARS, "mean_gen", "sd_gen")
P_SI_VARS <- sprintf("p_si[%s]", seq_len(length(seq(SI_MIN, SI_MAX, by = GRID_BY_SI))))
LOG_LIK_VARS <- sprintf("log_lik_cell[%s]", seq_len(length(seq(SI_MIN, SI_MAX, by = GRID_BY_SI))))

TRAJECTORY_PARAMS <- c("log_r", "log_d", "m", "log_kappa", "log_sigma_z")
ALL_PARAMS <- PARAM_VARS
PRIOR_TARGET_SET <- TRAJECTORY_PARAMS

# Helpers
normalize_density <- function(raw, dx) {
  z <- sum(raw) * dx
  if (is.finite(z) && z > 0) raw / z else rep(1 / (length(raw) * dx), length(raw))
}

normalize_mass <- function(raw) {
  z <- sum(raw)
  if (is.finite(z) && z > 0) raw / z else rep(1 / length(raw), length(raw))
}

cumulative_trapezoid <- function(y, dx) {
  out <- numeric(length(y))
  if (length(y) > 1) {
    for (i in 2:length(y)) out[i] <- out[i - 1] + dx * 0.5 * (y[i - 1] + y[i])
  }
  out
}

viral_shape_vec <- function(t, log_r, log_d, m, log_kappa) {
  rr <- exp(log_r)
  dd <- exp(log_d)
  kappa <- exp(log_kappa)

  a <- log(rr / dd) + kappa * (t - m)
  log1pexp_a <- ifelse(a > 40, a, log1p(exp(a)))
  log_v <- rr * t - ((rr + dd) / kappa) * log1pexp_a
  exp(pmin(log_v, 700))
}

find_grid_index <- function(x, grid, by, tol = 1e-8) {
  idx <- round((x - grid[1]) / by) + 1L
  if (any(idx < 1L | idx > length(grid))) {
    stop("Value outside grid: ", paste(head(x[idx < 1L | idx > length(grid)], 10), collapse = ", "))
  }
  if (any(abs(grid[idx] - x) > tol)) {
    stop("Value not exactly on grid: ", paste(head(x[abs(grid[idx] - x) > tol], 10), collapse = ", "))
  }
  as.integer(idx)
}

make_z_quad <- function(n_z = 7, z_max_sd = 3.5) {
  x <- seq(-z_max_sd, z_max_sd, length.out = n_z)
  dx <- if (length(x) > 1) diff(x)[1] else 1
  w <- dnorm(x) * dx
  list(x = x, w = w / sum(w))
}

make_grids_and_lookups <- function() {
  si_grid <- seq(SI_MIN, SI_MAX, by = GRID_BY_SI)
  l_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
  u_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
  gt_grid <- seq(0, GT_MAX, by = GRID_BY_LU)
  ip_grid <- seq(IP_MIN, IP_MAX, by = GRID_BY_IP)
  t_grid <- seq(0, GT_MAX, by = GRID_BY_LU)

  delta_min <- min(si_grid) - max(gt_grid)
  delta_max <- max(si_grid) - min(gt_grid)
  delta_grid <- seq(delta_min, delta_max, by = GRID_BY_SI)

  stopifnot(abs(GRID_BY_LU - GRID_BY_SI) < 1e-12)
  stopifnot(abs((GRID_BY_SI / GRID_BY_IP) - round(GRID_BY_SI / GRID_BY_IP)) < 1e-12)

  l_t_idx <- find_grid_index(l_grid, t_grid, GRID_BY_LU)

  lu_t_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
  gt_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
  for (jl in seq_along(l_grid)) {
    for (ju in seq_along(u_grid)) {
      lu_t_idx[jl, ju] <- find_grid_index(l_grid[jl] + u_grid[ju], t_grid, GRID_BY_LU)
      gt_idx[jl, ju] <- find_grid_index(l_grid[jl] + u_grid[ju], gt_grid, GRID_BY_LU)
    }
  }

  si_g_delta_idx <- matrix(NA_integer_, nrow = length(si_grid), ncol = length(gt_grid))
  for (ks in seq_along(si_grid)) {
    for (kg in seq_along(gt_grid)) {
      si_g_delta_idx[ks, kg] <- find_grid_index(si_grid[ks] - gt_grid[kg], delta_grid, GRID_BY_SI)
    }
  }

  ip2_idx <- matrix(0L, nrow = length(delta_grid), ncol = length(ip_grid))
  for (dd in seq_along(delta_grid)) {
    target <- delta_grid[dd] + ip_grid
    idx <- round((target - ip_grid[1]) / GRID_BY_IP) + 1L
    ok <- idx >= 1L & idx <= length(ip_grid)
    if (any(ok)) {
      ok_idx <- which(ok)
      ok[ok_idx] <- abs(ip_grid[idx[ok_idx]] - target[ok_idx]) <= 1e-8
    }
    ip2_idx[dd, ok] <- as.integer(idx[ok])
  }

  sparse_idx <- which(ip2_idx > 0L, arr.ind = TRUE)
  ip_sparse_delta <- as.integer(sparse_idx[, 1])
  ip_sparse_i1 <- as.integer(sparse_idx[, 2])
  ip_sparse_i2 <- as.integer(ip2_idx[sparse_idx])

  zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

  list(
    si_grid = si_grid,
    l_grid = l_grid,
    u_grid = u_grid,
    gt_grid = gt_grid,
    ip_grid = ip_grid,
    t_grid = t_grid,
    delta_grid = delta_grid,
    l_t_idx = l_t_idx,
    lu_t_idx = lu_t_idx,
    gt_idx = gt_idx,
    si_g_delta_idx = si_g_delta_idx,
    ip2_idx = ip2_idx,
    N_ip_sparse = length(ip_sparse_delta),
    ip_sparse_delta = ip_sparse_delta,
    ip_sparse_i1 = ip_sparse_i1,
    ip_sparse_i2 = ip_sparse_i2,
    z_std_grid = zq$x,
    z_std_w = zq$w
  )
}

G <- make_grids_and_lookups()
message(
  "Grid sizes: SI=", length(G$si_grid),
  ", l=", length(G$l_grid),
  ", u=", length(G$u_grid),
  ", GT=", length(G$gt_grid),
  ", IP=", length(G$ip_grid),
  ", IP sparse=", G$N_ip_sparse,
  ", z=", length(G$z_std_grid)
)

mechanistic_distribution <- function(par, r_growth, G) {
  lambda1 <- exp(par$log_lambda1)
  lambda_I <- exp(par$log_lambda_I)
  sigma_z <- exp(par$log_sigma_z)

  base_l <- viral_shape_vec(G$l_grid, par$log_r, par$log_d, par$m, par$log_kappa)
  base_t <- viral_shape_vec(G$t_grid, par$log_r, par$log_d, par$m, par$log_kappa)
  base_ip <- viral_shape_vec(G$ip_grid, par$log_r, par$log_d, par$m, par$log_kappa)

  N_l <- length(G$l_grid); N_u <- length(G$u_grid); N_gt <- length(G$gt_grid)
  N_ip <- length(G$ip_grid); N_z <- length(G$z_std_grid); N_delta <- length(G$delta_grid)
  N_si <- length(G$si_grid)

  fG_mass_z <- matrix(0, nrow = N_gt, ncol = N_z)
  fIP_z <- matrix(0, nrow = N_ip, ncol = N_z)

  for (jz in seq_len(N_z)) {
    z <- sigma_z * G$z_std_grid[jz]
    ez <- exp(z)

    hL <- ez * lambda1 * base_l
    HL <- cumulative_trapezoid(hL, GRID_BY_LU)
    normL <- normalize_density(hL * exp(-HL), GRID_BY_LU)

    for (jl in seq_len(N_l)) {
      rawU <- base_t[G$lu_t_idx[jl, ]]
      normU <- normalize_density(rawU, GRID_BY_LU)
      for (ju in seq_len(N_u)) {
        kg <- G$gt_idx[jl, ju]
        fG_mass_z[kg, jz] <- fG_mass_z[kg, jz] + normL[jl] * normU[ju] * GRID_BY_LU * GRID_BY_LU
      }
    }
    fG_mass_z[, jz] <- normalize_mass(fG_mass_z[, jz])

    hI <- ez * lambda_I * base_ip
    HI <- cumulative_trapezoid(hI, GRID_BY_IP)
    fIP_z[, jz] <- normalize_density(hI * exp(-HI), GRID_BY_IP)
  }

  fIP2 <- as.vector(fIP_z %*% G$z_std_w)
  fIP2 <- normalize_density(fIP2, GRID_BY_IP)

  C_delta_z <- matrix(0, nrow = N_delta, ncol = N_z)
  tilt <- exp(-r_growth * G$ip_grid)
  for (jz in seq_len(N_z)) {
    for (qq in seq_len(G$N_ip_sparse)) {
      dd <- G$ip_sparse_delta[qq]
      ii <- G$ip_sparse_i1[qq]
      jj <- G$ip_sparse_i2[qq]
      C_delta_z[dd, jz] <- C_delta_z[dd, jz] + tilt[ii] * fIP_z[ii, jz] * fIP2[jj] * GRID_BY_IP
    }
  }

  unnorm_si <- numeric(N_si)
  for (ks in seq_len(N_si)) {
    acc_s <- 0
    for (kg in seq_len(N_gt)) {
      dd <- G$si_g_delta_idx[ks, kg]
      acc_s <- acc_s + sum(G$z_std_w * fG_mass_z[kg, ] * C_delta_z[dd, ])
    }
    unnorm_si[ks] <- max(acc_s, 1e-300)
  }

  p_si <- normalize_mass(unnorm_si)
  fG_mix_mass <- as.vector(fG_mass_z %*% G$z_std_w)
  fG_mix_mass <- normalize_mass(fG_mix_mass)
  mean_gen <- sum(G$gt_grid * fG_mix_mass)
  sd_gen <- sqrt(sum((G$gt_grid - mean_gen)^2 * fG_mix_mass))

  list(p_si = p_si, fG_mass = fG_mix_mass, mean_gen = mean_gen, sd_gen = sd_gen)
}

TRUE_DIST <- mechanistic_distribution(TRUE_PARAMETERS, TRUE_R_GROWTH, G)
TRUE_SUMMARY <- tibble(
  true_mean_gen = TRUE_DIST$mean_gen,
  true_sd_gen = TRUE_DIST$sd_gen
)
write_csv(TRUE_SUMMARY, file.path(OUT_DIR, "true_generation_time_targets.csv"))
print(TRUE_SUMMARY)

# Prior scenarios
manuscript_prior <- list(
  prior_log_lambda1_mean = log(0.15),
  prior_log_lambda1_sd = 1.0,
  prior_log_lambda_I_mean = log(0.15),
  prior_log_lambda_I_sd = 1.0,
  prior_log_r_mean = log(0.25),
  prior_log_r_sd = 0.75,
  prior_log_d_mean = log(0.60),
  prior_log_d_sd = 0.75,
  prior_m_mean = 3.5,
  prior_m_sd = 1.5,
  prior_log_kappa_mean = 0.0,
  prior_log_kappa_sd = 1.0,
  prior_log_sigma_z_mean = log(0.35),
  prior_log_sigma_z_sd = 0.75
)

strength_sds <- list(
  weak = list(log = 0.75, m = 1.50),
  moderate = list(log = 0.35, m = 0.60),
  strong = list(log = 0.15, m = 0.25)
)

accuracy_offsets <- list(
  accurate = list(log_r = 0, log_d = 0, m = 0, log_kappa = 0, log_sigma_z = 0),
  mildly_wrong = list(log_r = log(1.25), log_d = log(0.85), m = 0.45, log_kappa = log(1.25), log_sigma_z = log(1.20)),
  badly_wrong = list(log_r = log(1.60), log_d = log(0.65), m = 0.90, log_kappa = log(1.80), log_sigma_z = log(1.60))
)

make_prior <- function(strength = "weak", accuracy = "accurate", target_params = PRIOR_TARGET_SET) {
  prior <- manuscript_prior
  sds <- strength_sds[[strength]]
  offsets <- accuracy_offsets[[accuracy]]

  for (v in target_params) {
    mean_name <- paste0("prior_", v, "_mean")
    sd_name <- paste0("prior_", v, "_sd")
    prior[[mean_name]] <- TRUE_PARAMETERS[[v]] + offsets[[v]]
    prior[[sd_name]] <- if (v == "m") sds$m else sds$log
  }

  prior
}

prior_design <- tidyr::expand_grid(
  accuracy = c("accurate", "mildly_wrong", "badly_wrong"),
  strength = c("weak", "moderate", "strong")
) %>%
  mutate(
    prior_scenario = paste(accuracy, strength, sep = "_"),
    target = "trajectory"
  ) %>%
  bind_rows(tibble(
    accuracy = "manuscript_default",
    strength = "default",
    prior_scenario = "manuscript_default",
    target = "trajectory"
  )) %>%
  mutate(prior_order = row_number())

make_prior_for_row <- function(row) {
  if (row$prior_scenario == "manuscript_default") {
    manuscript_prior
  } else {
    make_prior(strength = row$strength, accuracy = row$accuracy, target_params = PRIOR_TARGET_SET)
  }
}

prior_meta <- prior_design %>%
  mutate(prior_list = purrr::pmap(., function(...) make_prior_for_row(tibble(...)))) %>%
  select(-prior_list) %>%
  arrange(prior_order)
write_csv(prior_meta, file.path(OUT_DIR, "prior_scenario_design.csv"))
print(prior_meta)

prior_truth_distance_table <- function(prior_list, scenario_name) {
  tibble(variable = PARAM_VARS) %>%
    mutate(
      prior_mean = map_dbl(variable, ~ prior_list[[paste0("prior_", .x, "_mean")]]),
      prior_sd = map_dbl(variable, ~ prior_list[[paste0("prior_", .x, "_sd")]]),
      true_value = map_dbl(variable, ~ TRUE_PARAMETERS[[.x]]),
      prior_error = prior_mean - true_value,
      prior_abs_error = abs(prior_error),
      prior_error_in_sd = prior_error / prior_sd,
      prior_scenario = scenario_name
    )
}

prior_truth_distances <- purrr::map2_dfr(
  purrr::pmap(prior_design, function(...) make_prior_for_row(tibble(...))),
  prior_design$prior_scenario,
  prior_truth_distance_table
)
write_csv(prior_truth_distances, file.path(OUT_DIR, "prior_truth_distances.csv"))

# Stan data, initial values, summaries
make_stan_data <- function(counts_full, prior_list, r_growth, G) {
  c(
    list(
      N_norm = length(G$si_grid),
      si_grid = as.numeric(G$si_grid),
      dx_si = GRID_BY_SI,
      counts_full = as.integer(counts_full),
      n_total = as.integer(sum(counts_full)),
      N_l = length(G$l_grid),
      l_grid = as.numeric(G$l_grid),
      dx_l = GRID_BY_LU,
      N_u = length(G$u_grid),
      u_grid = as.numeric(G$u_grid),
      dx_u = GRID_BY_LU,
      N_gt = length(G$gt_grid),
      gt_grid = as.numeric(G$gt_grid),
      N_ip = length(G$ip_grid),
      ip_grid = as.numeric(G$ip_grid),
      dx_ip = GRID_BY_IP,
      N_t = length(G$t_grid),
      t_grid = as.numeric(G$t_grid),
      dx_t = GRID_BY_LU,
      N_delta = length(G$delta_grid),
      delta_grid = as.numeric(G$delta_grid),
      l_t_idx = as.integer(G$l_t_idx),
      lu_t_idx = G$lu_t_idx,
      gt_idx = G$gt_idx,
      si_g_delta_idx = G$si_g_delta_idx,
      N_ip_sparse = as.integer(G$N_ip_sparse),
      ip_sparse_delta = as.integer(G$ip_sparse_delta),
      ip_sparse_i1 = as.integer(G$ip_sparse_i1),
      ip_sparse_i2 = as.integer(G$ip_sparse_i2),
      park_weight_ip = as.numeric(exp(-r_growth * G$ip_grid)),
      N_z = length(G$z_std_grid),
      z_std_grid = as.numeric(G$z_std_grid),
      z_std_w = as.numeric(G$z_std_w)
    ),
    prior_list
  )
}

make_inits <- function(prior_list, jitter_sd = 0.15) {
  force(prior_list)
  function() {
    list(
      log_lambda1 = prior_list$prior_log_lambda1_mean + rnorm(1, 0, jitter_sd),
      log_lambda_I = prior_list$prior_log_lambda_I_mean + rnorm(1, 0, jitter_sd),
      log_r = prior_list$prior_log_r_mean + rnorm(1, 0, jitter_sd),
      log_d = prior_list$prior_log_d_mean + rnorm(1, 0, jitter_sd),
      m = min(max(prior_list$prior_m_mean + rnorm(1, 0, 0.20), 0.3), 11.5),
      log_kappa = prior_list$prior_log_kappa_mean + rnorm(1, 0, jitter_sd),
      log_sigma_z = prior_list$prior_log_sigma_z_mean + rnorm(1, 0, jitter_sd)
    )
  }
}

extract_matrix <- function(fit, variables) {
  draws <- fit$draws(variables = variables, format = "draws_matrix")
  as.matrix(draws)
}

get_q <- function(s, variable, col) {
  s %>% filter(.data$variable == variable) %>% pull({{ col }})
}

summarise_fit <- function(fit, prior_scenario, accuracy, strength, rep_id, truth, prior_list, counts_test) {
  s <- fit$summary(variables = FIT_SUMMARY_VARS)
  mean_row <- s %>% filter(variable == "mean_gen")
  sd_row <- s %>% filter(variable == "sd_gen")

  p_draws <- extract_matrix(fit, P_SI_VARS)
  p_draws <- pmax(p_draws, 1e-300)
  p_post_mean <- colMeans(p_draws)
  p_true <- truth$p_si

  train_log_lik <- rowSums(extract_matrix(fit, LOG_LIK_VARS))
  expected_train_lpd <- mean(train_log_lik)

  test_log_score_draws <- as.vector(p_draws %*% counts_test)
  # The line above is intentionally replaced below; matrix multiplication with probabilities
  # is not a log score. Keep this comment as a guard against accidental reuse.
  test_log_score_draws <- as.vector(log(p_draws) %*% counts_test)
  expected_test_lpd <- mean(test_log_score_draws)

  kl_true_to_postmean <- sum(p_true * log(pmax(p_true, 1e-300) / pmax(p_post_mean, 1e-300)))
  tv_true_to_postmean <- 0.5 * sum(abs(p_true - p_post_mean))
  rmse_p_si <- sqrt(mean((p_true - p_post_mean)^2))

  log_lik_matrix <- extract_matrix(fit, LOG_LIK_VARS)
  waic_obj <- tryCatch(loo::waic(log_lik_matrix), error = function(e) NULL)
  loo_obj <- tryCatch(loo::loo(log_lik_matrix), error = function(e) NULL)

  diag_df <- fit$sampler_diagnostics(format = "df")

  tibble(
    prior_scenario = prior_scenario,
    accuracy = accuracy,
    strength = strength,
    rep_id = rep_id,
    n_pairs = N_PAIRS,
    true_mean_gen = truth$mean_gen,
    true_sd_gen = truth$sd_gen,
    mean_gen_median = mean_row$median,
    mean_gen_q05 = mean_row$q5,
    mean_gen_q95 = mean_row$q95,
    mean_gen_error = mean_row$median - truth$mean_gen,
    mean_gen_abs_error = abs(mean_row$median - truth$mean_gen),
    mean_gen_ci_width = mean_row$q95 - mean_row$q5,
    mean_gen_covered = mean_row$q5 <= truth$mean_gen & mean_row$q95 >= truth$mean_gen,
    sd_gen_median = sd_row$median,
    sd_gen_q05 = sd_row$q5,
    sd_gen_q95 = sd_row$q95,
    sd_gen_error = sd_row$median - truth$sd_gen,
    sd_gen_abs_error = abs(sd_row$median - truth$sd_gen),
    sd_gen_ci_width = sd_row$q95 - sd_row$q5,
    sd_gen_covered = sd_row$q5 <= truth$sd_gen & sd_row$q95 >= truth$sd_gen,
    expected_train_lpd = expected_train_lpd,
    expected_test_lpd = expected_test_lpd,
    test_lpd_per_pair = expected_test_lpd / sum(counts_test),
    kl_true_to_postmean = kl_true_to_postmean,
    tv_true_to_postmean = tv_true_to_postmean,
    rmse_p_si = rmse_p_si,
    waic = if (is.null(waic_obj)) NA_real_ else as.numeric(waic_obj$estimates["waic", "Estimate"]),
    looic = if (is.null(loo_obj)) NA_real_ else as.numeric(loo_obj$estimates["looic", "Estimate"]),
    max_pareto_k = if (is.null(loo_obj)) NA_real_ else suppressWarnings(max(loo_obj$diagnostics$pareto_k, na.rm = TRUE)),
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(diag_df$divergent__, na.rm = TRUE),
    max_treedepth_hit = sum(diag_df$treedepth__ >= MAX_TREEDEPTH, na.rm = TRUE)
  )
}

summarise_parameter_recovery <- function(fit, prior_scenario, accuracy, strength, rep_id, true_par, prior_list) {
  truth_tbl <- tibble(variable = names(true_par), true_value = as.numeric(unlist(true_par))) %>%
    mutate(
      prior_mean = map_dbl(variable, ~ prior_list[[paste0("prior_", .x, "_mean")]]),
      prior_sd = map_dbl(variable, ~ prior_list[[paste0("prior_", .x, "_sd")]])
    )

  fit$summary(variables = PARAM_VARS) %>%
    select(variable, median, q5, q95, rhat, ess_bulk, ess_tail) %>%
    left_join(truth_tbl, by = "variable") %>%
    mutate(
      prior_scenario = prior_scenario,
      accuracy = accuracy,
      strength = strength,
      rep_id = rep_id,
      posterior_error = median - true_value,
      posterior_abs_error = abs(posterior_error),
      prior_error = prior_mean - true_value,
      prior_abs_error = abs(prior_error),
      improvement_over_prior = prior_abs_error - posterior_abs_error,
      posterior_moved_toward_truth = posterior_abs_error < prior_abs_error,
      interval_width = q95 - q5,
      covered = q5 <= true_value & q95 >= true_value,
      posterior_shift_from_prior_sd = (median - prior_mean) / prior_sd,
      prior_truth_distance_sd = (prior_mean - true_value) / prior_sd
    ) %>%
    select(
      prior_scenario, accuracy, strength, rep_id, variable,
      true_value, prior_mean, prior_sd, median, q5, q95,
      prior_error, prior_abs_error, posterior_error, posterior_abs_error,
      improvement_over_prior, posterior_moved_toward_truth,
      interval_width, covered, posterior_shift_from_prior_sd,
      prior_truth_distance_sd, rhat, ess_bulk, ess_tail
    )
}

naturalise_param_results <- function(param_results) {
  positive_log_params <- c("log_lambda1", "log_lambda_I", "log_r", "log_d", "log_kappa", "log_sigma_z")

  param_results %>%
    mutate(
      natural_variable = case_when(
        variable == "log_lambda1" ~ "lambda1",
        variable == "log_lambda_I" ~ "lambda_I",
        variable == "log_r" ~ "r",
        variable == "log_d" ~ "d",
        variable == "log_kappa" ~ "kappa",
        variable == "log_sigma_z" ~ "sigma_z",
        TRUE ~ variable
      ),
      natural_true = ifelse(variable %in% positive_log_params, exp(true_value), true_value),
      natural_prior_mean = ifelse(variable %in% positive_log_params, exp(prior_mean), prior_mean),
      natural_median = ifelse(variable %in% positive_log_params, exp(median), median),
      natural_q5 = ifelse(variable %in% positive_log_params, exp(q5), q5),
      natural_q95 = ifelse(variable %in% positive_log_params, exp(q95), q95),
      natural_prior_abs_error = abs(natural_prior_mean - natural_true),
      natural_posterior_abs_error = abs(natural_median - natural_true),
      natural_improvement_over_prior = natural_prior_abs_error - natural_posterior_abs_error,
      natural_interval_width = natural_q95 - natural_q5,
      natural_covered = natural_q5 <= natural_true & natural_q95 >= natural_true
    )
}

# Compile Stan model
Sys.setenv(MAKEFLAGS = paste0("-j", 1))
message("Compiling FAST Stan model: ", STAN_FILE)
mod <- cmdstan_model(STAN_FILE, force_recompile = FORCE_RECOMPILE)
CMDSTAN_OUTPUT_DIR <- file.path(OUT_DIR, "cmdstan_outputs")
dir.create(CMDSTAN_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
JOB_RDS_DIR <- file.path(OUT_DIR, "job_rds")
dir.create(JOB_RDS_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)

simulated_data <- vector("list", N_REP)
for (rep_id in seq_len(N_REP)) {
  simulated_data[[rep_id]] <- list(
    counts_train = as.vector(rmultinom(1, size = N_PAIRS, prob = TRUE_DIST$p_si)),
    counts_test = as.vector(rmultinom(1, size = N_TEST, prob = TRUE_DIST$p_si))
  )
}

job_table <- tidyr::crossing(
  rep_id = seq_len(N_REP),
  pd_i = seq_len(nrow(prior_design))
) %>%
  arrange(rep_id, pd_i) %>%
  mutate(job_id = row_number())

write_csv(job_table, file.path(OUT_DIR, "simulation_job_table.csv"))
message("Total Stan fits to run: ", nrow(job_table), " = ", N_REP, " replicates × ", nrow(prior_design), " prior scenarios")

run_one_job <- function(job_row) {
  # job_row arrives as a one-row tibble/data.frame from lapply/future_lapply.
  rep_id <- as.integer(job_row$rep_id)
  pd_i <- as.integer(job_row$pd_i)
  job_id <- as.integer(job_row$job_id)

  pd <- prior_design[pd_i, ]
  prior_scenario <- pd$prior_scenario
  accuracy <- pd$accuracy
  strength <- pd$strength
  prior_list <- make_prior_for_row(pd)

  counts_train <- simulated_data[[rep_id]]$counts_train
  counts_test <- simulated_data[[rep_id]]$counts_test
  stan_data <- make_stan_data(counts_train, prior_list, TRUE_R_GROWTH, G)

  job_label <- sprintf("job%04d_rep%03d_%s", job_id, rep_id, prior_scenario)
  job_output_dir <- file.path(CMDSTAN_OUTPUT_DIR, job_label)
  dir.create(job_output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Running ", job_label)

  fit_result <- tryCatch({
    fit <- mod$sample(
      data = stan_data,
      chains = CHAINS,
      parallel_chains = PARALLEL_CHAINS,
      iter_warmup = ITER_WARMUP,
      iter_sampling = ITER_SAMPLING,
      seed = SEED + 100000 * rep_id + pd_i,
      init = make_inits(prior_list, jitter_sd = 0.12),
      adapt_delta = ADAPT_DELTA,
      max_treedepth = MAX_TREEDEPTH,
      refresh = REFRESH,
      save_warmup = FALSE,
      output_dir = job_output_dir
    )
    list(ok = TRUE, fit = fit, error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, fit = NULL, error = conditionMessage(e))
  })

  if (!fit_result$ok) {
    out <- list(
      ok = FALSE,
      fit_row = NULL,
      param_row = NULL,
      failure_row = tibble(
        prior_scenario = prior_scenario,
        accuracy = accuracy,
        strength = strength,
        rep_id = rep_id,
        job_id = job_id,
        error = fit_result$error
      )
    )
    if (WRITE_JOB_LEVEL_RDS) saveRDS(out, file.path(JOB_RDS_DIR, paste0(job_label, "_FAILED.rds")))
    return(out)
  }

  fit <- fit_result$fit

  fit_row <- summarise_fit(
    fit = fit,
    prior_scenario = prior_scenario,
    accuracy = accuracy,
    strength = strength,
    rep_id = rep_id,
    truth = TRUE_DIST,
    prior_list = prior_list,
    counts_test = counts_test
  ) %>% mutate(job_id = job_id, .before = rep_id)

  param_row <- summarise_parameter_recovery(
    fit = fit,
    prior_scenario = prior_scenario,
    accuracy = accuracy,
    strength = strength,
    rep_id = rep_id,
    true_par = TRUE_PARAMETERS,
    prior_list = prior_list
  ) %>% mutate(job_id = job_id, .before = rep_id)

  if (SAVE_FITS) {
    fit$save_object(file.path(OUT_DIR, sprintf("fit_%s_rep%03d.rds", prior_scenario, rep_id)))
  }

  out <- list(ok = TRUE, fit_row = fit_row, param_row = param_row, failure_row = NULL)
  if (WRITE_JOB_LEVEL_RDS) saveRDS(out, file.path(JOB_RDS_DIR, paste0(job_label, ".rds")))
  out
}

job_rows <- split(job_table, seq_len(nrow(job_table)))

if (PARALLELIZE_FITS && N_WORKERS > 1L) {
  message("Running independent fits in parallel with N_WORKERS = ", N_WORKERS)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (future::supportsMulticore()) {
    future::plan(future::multicore, workers = N_WORKERS)
  } else {
    future::plan(future::multisession, workers = N_WORKERS)
  }

  job_results <- future.apply::future_lapply(
    job_rows,
    run_one_job,
    future.seed = TRUE,
    future.scheduling = 1
  )
} else {
  message("Running fits serially")
  job_results <- lapply(job_rows, run_one_job)
}

# Reload all completed job-level RDS files after restarting R
JOB_RDS_DIR <- file.path(OUT_DIR, "job_rds")

job_rds_files <- list.files(
  JOB_RDS_DIR,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(job_rds_files) == 0) {
  stop("No job-level RDS files found in: ", JOB_RDS_DIR)
}

job_results_all <- lapply(job_rds_files, readRDS)

fit_results <- bind_rows(lapply(job_results_all, `[[`, "fit_row"))
param_results <- bind_rows(lapply(job_results_all, `[[`, "param_row"))
failures <- bind_rows(lapply(job_results_all, `[[`, "failure_row"))

fit_results <- fit_results %>%
  arrange(job_id) %>%
  distinct(job_id, .keep_all = TRUE)

param_results <- param_results %>%
  arrange(job_id) %>%
  distinct(job_id, variable, .keep_all = TRUE)

if (nrow(failures) > 0) {
  failures <- failures %>%
    arrange(job_id) %>%
    distinct(job_id, .keep_all = TRUE)
}

# Create the natural-scale object that your later code expects
param_results_natural <- naturalise_param_results(param_results)

# Save the combined/reloaded results
write_csv(fit_results, file.path(OUT_DIR, "fit_level_results.csv"))
write_csv(param_results, file.path(OUT_DIR, "parameter_recovery_results.csv"))
write_csv(
  param_results_natural,
  file.path(OUT_DIR, "parameter_recovery_results_natural_scale.csv")
)

if (nrow(failures) > 0) {
  write_csv(failures, file.path(OUT_DIR, "simulation_failures.csv"))
}

# Paper-ready summaries
fit_summary_by_prior <- fit_results %>%
  group_by(prior_scenario, accuracy, strength) %>%
  summarise(
    n_successful_fits = n(),
    failure_rate = 1 - n() / N_REP,
    median_bias_mean_gen = median(mean_gen_error, na.rm = TRUE),
    median_abs_error_mean_gen = median(mean_gen_abs_error, na.rm = TRUE),
    rmse_mean_gen = sqrt(mean(mean_gen_error^2, na.rm = TRUE)),
    coverage_mean_gen = mean(mean_gen_covered, na.rm = TRUE),
    median_ci_width_mean_gen = median(mean_gen_ci_width, na.rm = TRUE),
    median_bias_sd_gen = median(sd_gen_error, na.rm = TRUE),
    median_abs_error_sd_gen = median(sd_gen_abs_error, na.rm = TRUE),
    rmse_sd_gen = sqrt(mean(sd_gen_error^2, na.rm = TRUE)),
    coverage_sd_gen = mean(sd_gen_covered, na.rm = TRUE),
    median_ci_width_sd_gen = median(sd_gen_ci_width, na.rm = TRUE),
    median_test_lpd_per_pair = median(test_lpd_per_pair, na.rm = TRUE),
    median_kl_true_to_postmean = median(kl_true_to_postmean, na.rm = TRUE),
    median_tv_true_to_postmean = median(tv_true_to_postmean, na.rm = TRUE),
    median_rmse_p_si = median(rmse_p_si, na.rm = TRUE),
    median_waic = median(waic, na.rm = TRUE),
    median_looic = median(looic, na.rm = TRUE),
    convergence_ok_rate = mean(max_rhat < 1.05 & n_divergent / (CHAINS * ITER_SAMPLING) < 0.05, na.rm = TRUE),
    median_max_rhat = median(max_rhat, na.rm = TRUE),
    median_min_ess_bulk = median(min_ess_bulk, na.rm = TRUE),
    median_min_ess_tail = median(min_ess_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(accuracy, strength)

write_csv(fit_summary_by_prior, file.path(OUT_DIR, "paper_fit_summary_by_prior.csv"))
print(fit_summary_by_prior)

param_summary_by_prior <- param_results %>%
  group_by(prior_scenario, accuracy, strength, variable) %>%
  summarise(
    n_successful_fits = n(),
    true_value = first(true_value),
    median_prior_abs_error = median(prior_abs_error, na.rm = TRUE),
    median_posterior_abs_error = median(posterior_abs_error, na.rm = TRUE),
    median_improvement_over_prior = median(improvement_over_prior, na.rm = TRUE),
    moved_toward_truth_rate = mean(posterior_moved_toward_truth, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    median_interval_width = median(interval_width, na.rm = TRUE),
    median_shift_from_prior_sd = median(posterior_shift_from_prior_sd, na.rm = TRUE),
    median_abs_shift_from_prior_sd = median(abs(posterior_shift_from_prior_sd), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(prior_scenario, variable)

write_csv(param_summary_by_prior, file.path(OUT_DIR, "paper_parameter_summary_by_prior.csv"))

param_summary_by_prior_natural <- param_results_natural %>%
  group_by(prior_scenario, accuracy, strength, natural_variable) %>%
  summarise(
    n_successful_fits = n(),
    true_value = first(natural_true),
    median_prior_abs_error = median(natural_prior_abs_error, na.rm = TRUE),
    median_posterior_abs_error = median(natural_posterior_abs_error, na.rm = TRUE),
    median_improvement_over_prior = median(natural_improvement_over_prior, na.rm = TRUE),
    coverage = mean(natural_covered, na.rm = TRUE),
    median_interval_width = median(natural_interval_width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(prior_scenario, natural_variable)

write_csv(param_summary_by_prior_natural, file.path(OUT_DIR, "paper_parameter_summary_by_prior_natural_scale.csv"))

# Contrasts against manuscript/default and weak accurate prior.
baseline_tbl <- fit_results %>%
  filter(prior_scenario == "manuscript_default") %>%
  select(rep_id, baseline_test_lpd = test_lpd_per_pair, baseline_kl = kl_true_to_postmean, baseline_mean_abs_error = mean_gen_abs_error)

contrasts <- fit_results %>%
  left_join(baseline_tbl, by = "rep_id") %>%
  mutate(
    delta_test_lpd_vs_default = test_lpd_per_pair - baseline_test_lpd,
    delta_kl_vs_default = kl_true_to_postmean - baseline_kl,
    delta_mean_abs_error_vs_default = mean_gen_abs_error - baseline_mean_abs_error
  ) %>%
  group_by(prior_scenario, accuracy, strength) %>%
  summarise(
    median_delta_test_lpd_vs_default = median(delta_test_lpd_vs_default, na.rm = TRUE),
    q05_delta_test_lpd_vs_default = quantile(delta_test_lpd_vs_default, 0.05, na.rm = TRUE),
    q95_delta_test_lpd_vs_default = quantile(delta_test_lpd_vs_default, 0.95, na.rm = TRUE),
    median_delta_kl_vs_default = median(delta_kl_vs_default, na.rm = TRUE),
    median_delta_mean_abs_error_vs_default = median(delta_mean_abs_error_vs_default, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(contrasts, file.path(OUT_DIR, "paper_prior_contrasts_vs_default.csv"))

# Figures
scenario_levels <- prior_design$prior_scenario
plot_fit <- fit_results %>%
  mutate(prior_scenario = factor(prior_scenario, levels = scenario_levels))

p_lpd <- ggplot(plot_fit, aes(x = prior_scenario, y = test_lpd_per_pair)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.65) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), alpha = 0.20, size = 0.8) +
  coord_flip() +
  labs(
    x = "Prior scenario",
    y = "Held-out log predictive density per pair",
    title = "Predictive fit under accurate and misspecified mechanistic priors",
    subtitle = "Higher values indicate better held-out prediction. Test data are generated from the known truth."
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_1_test_log_score_by_prior.png"), p_lpd, width = 8.8, height = 6.2, dpi = 300)

p_kl <- ggplot(plot_fit, aes(x = prior_scenario, y = kl_true_to_postmean)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.65) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 2), alpha = 0.20, size = 0.8) +
  coord_flip() +
  labs(
    x = "Prior scenario",
    y = "KL divergence: true SI distribution to posterior mean SI distribution",
    title = "Distance between true and fitted serial-interval distributions",
    subtitle = "Lower values indicate a fitted observed-scale distribution closer to the truth."
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_2_kl_by_prior.png"), p_kl, width = 8.8, height = 6.2, dpi = 300)

p_gt <- ggplot(plot_fit, aes(x = prior_scenario, y = mean_gen_median, ymin = mean_gen_q05, ymax = mean_gen_q95)) +
  geom_hline(yintercept = TRUE_DIST$mean_gen, linetype = "dashed", linewidth = 0.7) +
  geom_pointrange(
    position = position_jitter(width = 0.12, height = 0, seed = 3),
    alpha = 0.35,
    linewidth = 0.30
  ) +
  coord_flip() +
  labs(
    x = "Prior scenario",
    y = "Estimated mean generation time, days",
    title = "Generation-time mean recovery under varying prior accuracy and strength",
    subtitle = "Dashed line is the known data-generating value; intervals are 90% posterior intervals."
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_3_mean_gt_recovery_by_prior.png"), p_gt, width = 8.8, height = 6.2, dpi = 300)

plot_param <- param_results %>%
  filter(variable %in% PRIOR_TARGET_SET) %>%
  mutate(prior_scenario = factor(prior_scenario, levels = scenario_levels))

p_shift <- ggplot(plot_param, aes(x = prior_scenario, y = posterior_shift_from_prior_sd)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.20, width = 0.65) +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    x = "Prior scenario",
    y = "Posterior median shift from prior mean, in prior SDs",
    title = "Prior-data conflict and posterior movement",
    subtitle = "Large shifts indicate that the timing data pull the parameter away from the prior mean."
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_DIR, "figure_4_posterior_shift_from_prior.png"), p_shift, width = 11, height = 7.5, dpi = 300)

plot_param_nat <- param_results_natural %>%
  filter(variable %in% PRIOR_TARGET_SET) %>%
  mutate(prior_scenario = factor(prior_scenario, levels = scenario_levels))

p_param_nat <- ggplot(plot_param_nat, aes(x = prior_scenario, y = natural_median, ymin = natural_q5, ymax = natural_q95)) +
  geom_pointrange(position = position_jitter(width = 0.11, height = 0, seed = 4), alpha = 0.25, linewidth = 0.25) +
  geom_hline(aes(yintercept = natural_true), linetype = "dashed", data = plot_param_nat %>% distinct(natural_variable, natural_true)) +
  coord_flip() +
  facet_wrap(~ natural_variable, scales = "free_x") +
  labs(
    x = "Prior scenario",
    y = "Natural-scale parameter estimate",
    title = "Mechanistic-parameter recovery under accurate and misspecified priors",
    subtitle = "Dashed lines are known data-generating values."
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_DIR, "figure_5_natural_parameter_recovery.png"), p_param_nat, width = 11, height = 7.5, dpi = 300)

p_move <- param_results %>%
  filter(variable %in% PRIOR_TARGET_SET) %>%
  group_by(prior_scenario, accuracy, strength, variable) %>%
  summarise(rate = mean(posterior_moved_toward_truth, na.rm = TRUE), .groups = "drop") %>%
  mutate(prior_scenario = factor(prior_scenario, levels = scenario_levels)) %>%
  ggplot(aes(x = prior_scenario, y = rate)) +
  geom_col(width = 0.65) +
  coord_flip() +
  facet_wrap(~ variable) +
  labs(
    x = "Prior scenario",
    y = "Proportion of fits where posterior is closer to truth than prior mean",
    title = "Does the posterior correct dodgy mechanistic priors?"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_DIR, "figure_6_posterior_moves_toward_truth.png"), p_move, width = 11, height = 7.5, dpi = 300)
