# Simulation 1
setwd("~/Measles-county-ranking/GT/New V code/NEW V FIXED")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

OUT_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 1"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

this_file_dir <- "~/Measles-county-ranking/GT/New V code"
STAN_FILE <- "~/Within-host-time-delay-framework/Code/mechanistic_obspair.stan"

# Main simulation setup.
N_REP <- 100
N_SMALL <- c(30, 50, 100)

# MCMC settings for simulation refits.
CHAINS <- 4
PARALLEL_CHAINS <- 4
ITER_WARMUP <- 250
ITER_SAMPLING <- 250
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15
REFRESH <- 100
SEED <- 41001
SAVE_FITS <- TRUE
FORCE_RECOMPILE <- FALSE

# Numerical grid settings.
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

# Single data-generating scenario.
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

# Priors
PRIOR_LIST <- list(
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

# Grid and lookup-table helpers
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

# Mechanistic distribution calculator in R, used to simulate truth
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
      # Observed-pair timing: U is the selected observed transmission event,
      # not the first successful transmission. Conditional on LP = l,
      # f_U(u | l, obs) is proportional to V(l + u), normalized over the u-grid.
      # The transmission-intensity scale cancels and is not identifiable from
      # serial-interval pair timing alone.
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

# True model-implied distribution for the single scenario
TRUE_DIST <- mechanistic_distribution(TRUE_PARAMETERS, TRUE_R_GROWTH, G)
TRUE_SUMMARY <- tibble(
  scenario = "single_scenario",
  true_mean_gen = TRUE_DIST$mean_gen,
  true_sd_gen = TRUE_DIST$sd_gen
)
print(TRUE_SUMMARY)
write_csv(TRUE_SUMMARY, file.path(OUT_DIR, "true_generation_time_targets.csv"))

# Stan data builder and fitting functions
make_stan_data <- function(counts_full, r_growth, G) {
  keep <- which(counts_full > 0)
  c(
    list(
      N_si = length(keep),
      obs_index = as.integer(keep),
      counts = as.integer(counts_full[keep]),
      N_norm = length(G$si_grid),
      si_grid = as.numeric(G$si_grid),
      dx_si = GRID_BY_SI,
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
    PRIOR_LIST
  )
}

make_inits <- function(true_par = NULL, jitter_sd = 0.10) {
  force(true_par)
  function() {
    base <- if (is.null(true_par)) {
      list(
        log_lambda1 = log(0.15), log_lambda_I = log(0.15),
        log_r = log(0.25), log_d = log(0.60), m = 3.5,
        log_kappa = 0, log_sigma_z = log(0.35)
      )
    } else true_par
    list(
      log_lambda1 = base$log_lambda1 + rnorm(1, 0, jitter_sd),
      log_lambda_I = base$log_lambda_I + rnorm(1, 0, jitter_sd),
      log_r = base$log_r + rnorm(1, 0, jitter_sd),
      log_d = base$log_d + rnorm(1, 0, jitter_sd),
      m = min(max(base$m + rnorm(1, 0, 0.15), 0.3), 11.5),
      log_kappa = base$log_kappa + rnorm(1, 0, jitter_sd),
      log_sigma_z = base$log_sigma_z + rnorm(1, 0, jitter_sd)
    )
  }
}

summarise_fit <- function(fit, n_pairs, rep_id, truth) {
  s <- fit$summary(variables = FIT_SUMMARY_VARS)
  mean_row <- s %>% filter(variable == "mean_gen")
  sd_row <- s %>% filter(variable == "sd_gen")
  tibble(
    n_pairs = n_pairs,
    rep_id = rep_id,
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
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(fit$sampler_diagnostics(format = "df")$divergent__, na.rm = TRUE)
  )
}

summarise_parameter_recovery <- function(fit, n_pairs, rep_id, true_par) {
  truth_tbl <- tibble(variable = names(true_par), true_value = as.numeric(unlist(true_par)))
  
  fit$summary(variables = PARAM_VARS) %>%
    select(variable, median, q5, q95, rhat, ess_bulk, ess_tail) %>%
    left_join(truth_tbl, by = "variable") %>%
    mutate(
      n_pairs = n_pairs,
      rep_id = rep_id,
      error = median - true_value,
      abs_error = abs(error),
      interval_width = q95 - q5,
      covered = q5 <= true_value & q95 >= true_value
    ) %>%
    select(n_pairs, rep_id, variable, true_value, median, q5, q95,
           error, abs_error, interval_width, covered, rhat, ess_bulk, ess_tail)
}

# Compile model

# Use all available cores for C++ compilation where supported.
Sys.setenv(MAKEFLAGS = paste0("-j", max(1L, parallel::detectCores(logical = TRUE) - 1L)))

mod <- cmdstan_model(STAN_FILE, force_recompile = FORCE_RECOMPILE)

# Simulation loop
set.seed(SEED)
fit_rows <- list()
param_rows <- list()
failed_rows <- list()
row_id <- 1L
param_row_id <- 1L
fail_id <- 1L

for (n_pairs in N_SMALL) {
  for (rep_id in seq_len(N_REP)) {
    message("\nSimulation n=", n_pairs, ", replicate=", rep_id, "/", N_REP)
    
    counts_full <- as.vector(rmultinom(1, size = n_pairs, prob = TRUE_DIST$p_si))
    stan_data <- make_stan_data(counts_full, TRUE_R_GROWTH, G)
    
    fit_result <- tryCatch({
      fit <- mod$sample(
        data = stan_data,
        chains = CHAINS,
        parallel_chains = PARALLEL_CHAINS,
        iter_warmup = ITER_WARMUP,
        iter_sampling = ITER_SAMPLING,
        seed = SEED + 100000 * rep_id + n_pairs,
        init = make_inits(TRUE_PARAMETERS, jitter_sd = 0.15),
        adapt_delta = ADAPT_DELTA,
        max_treedepth = MAX_TREEDEPTH,
        refresh = REFRESH,
        save_warmup = FALSE
      )
      
      list(ok = TRUE, fit = fit, error = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, fit = NULL, error = conditionMessage(e))
    })
    
    if (!fit_result$ok) {
      failed_rows[[fail_id]] <- tibble(n_pairs = n_pairs, rep_id = rep_id, error = fit_result$error)
      fail_id <- fail_id + 1L
      if (length(failed_rows) > 0) write_csv(bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
      next
    }
    
    fit <- fit_result$fit
    fit_rows[[row_id]] <- summarise_fit(fit, n_pairs, rep_id, TRUE_DIST)
    row_id <- row_id + 1L
    
    param_rows[[param_row_id]] <- summarise_parameter_recovery(fit, n_pairs, rep_id, TRUE_PARAMETERS)
    param_row_id <- param_row_id + 1L
    
    if (SAVE_FITS) {
      fit$save_object(file.path(OUT_DIR, sprintf("fit_n%s_rep%s.rds", n_pairs, rep_id)))
    }
    
    # Incremental writes so a long run can be stopped without losing progress.
    if (length(fit_rows) > 0) write_csv(bind_rows(fit_rows), file.path(OUT_DIR, "simulation_fit_level_results.csv"))
    if (length(param_rows) > 0) write_csv(bind_rows(param_rows), file.path(OUT_DIR, "simulation_parameter_recovery_results.csv"))
    if (length(failed_rows) > 0) write_csv(bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
  }
}

fit_results <- bind_rows(fit_rows)
param_results <- bind_rows(param_rows)
failures <- bind_rows(failed_rows)

write_csv(fit_results, file.path(OUT_DIR, "simulation_fit_level_results.csv"))
write_csv(param_results, file.path(OUT_DIR, "simulation_parameter_recovery_results.csv"))
if (nrow(failures) > 0) write_csv(failures, file.path(OUT_DIR, "simulation_failures.csv"))

# Summaries
recovery_summary <- fit_results %>%
  group_by(n_pairs) %>%
  summarise(
    n_successful_fits = n(),
    failure_rate = 1 - n() / N_REP,
    median_bias_mean_gen = median(mean_gen_error, na.rm = TRUE),
    median_abs_error_mean_gen = median(mean_gen_abs_error, na.rm = TRUE),
    rmse_mean_gen = sqrt(mean(mean_gen_error^2, na.rm = TRUE)),
    median_ci_width_mean_gen = median(mean_gen_ci_width, na.rm = TRUE),
    coverage_mean_gen = mean(mean_gen_covered, na.rm = TRUE),
    median_bias_sd_gen = median(sd_gen_error, na.rm = TRUE),
    median_abs_error_sd_gen = median(sd_gen_abs_error, na.rm = TRUE),
    rmse_sd_gen = sqrt(mean(sd_gen_error^2, na.rm = TRUE)),
    median_ci_width_sd_gen = median(sd_gen_ci_width, na.rm = TRUE),
    coverage_sd_gen = mean(sd_gen_covered, na.rm = TRUE),
    convergence_ok_rate = mean(max_rhat < 1.05 & n_divergent / (CHAINS * ITER_SAMPLING) < 0.05, na.rm = TRUE),
    median_min_ess_bulk = median(min_ess_bulk, na.rm = TRUE),
    median_min_ess_tail = median(min_ess_tail, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(recovery_summary, file.path(OUT_DIR, "paper_recovery_summary_by_sample_size.csv"))
print(recovery_summary)

param_recovery_summary <- param_results %>%
  group_by(n_pairs, variable) %>%
  summarise(
    n_successful_fits = n(),
    true_value = first(true_value),
    median_bias = median(error, na.rm = TRUE),
    median_abs_error = median(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    median_interval_width = median(interval_width, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    median_rhat = median(rhat, na.rm = TRUE),
    median_ess_bulk = median(ess_bulk, na.rm = TRUE),
    median_ess_tail = median(ess_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_pairs, variable)

write_csv(param_recovery_summary, file.path(OUT_DIR, "paper_parameter_recovery_summary_by_sample_size.csv"))
print(param_recovery_summary)

positive_log_params <- c("log_lambda1", "log_lambda_I", "log_r", "log_d", "log_kappa", "log_sigma_z")

param_results_natural <- param_results %>%
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
    natural_median = ifelse(variable %in% positive_log_params, exp(median), median),
    natural_q5 = ifelse(variable %in% positive_log_params, exp(q5), q5),
    natural_q95 = ifelse(variable %in% positive_log_params, exp(q95), q95),
    natural_error = natural_median - natural_true,
    natural_abs_error = abs(natural_error),
    natural_interval_width = natural_q95 - natural_q5,
    natural_covered = natural_q5 <= natural_true & natural_q95 >= natural_true
  )

write_csv(param_results_natural, file.path(OUT_DIR, "simulation_parameter_recovery_results_natural_scale.csv"))

param_recovery_summary_natural <- param_results_natural %>%
  group_by(n_pairs, natural_variable) %>%
  summarise(
    n_successful_fits = n(),
    true_value = first(natural_true),
    median_bias = median(natural_error, na.rm = TRUE),
    median_abs_error = median(natural_abs_error, na.rm = TRUE),
    rmse = sqrt(mean(natural_error^2, na.rm = TRUE)),
    median_interval_width = median(natural_interval_width, na.rm = TRUE),
    coverage = mean(natural_covered, na.rm = TRUE),
    median_rhat = median(rhat, na.rm = TRUE),
    median_ess_bulk = median(ess_bulk, na.rm = TRUE),
    median_ess_tail = median(ess_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_pairs, natural_variable)

write_csv(param_recovery_summary_natural, file.path(OUT_DIR, "paper_parameter_recovery_summary_by_sample_size_natural_scale.csv"))
print(param_recovery_summary_natural)

# Figures
plot_results <- fit_results %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = N_SMALL),
    mean_gen_covered_f = ifelse(mean_gen_covered, "90% interval covers truth", "misses truth"),
    convergence_ok = max_rhat < 1.05 & n_divergent / (CHAINS * ITER_SAMPLING) < 0.05
  )

p_recovery_mean <- ggplot(
  plot_results,
  aes(x = n_pairs_f, y = mean_gen_median, ymin = mean_gen_q05, ymax = mean_gen_q95)
) +
  geom_hline(yintercept = TRUE_DIST$mean_gen, linetype = "dashed", linewidth = 0.7) +
  geom_pointrange(
    aes(shape = mean_gen_covered_f),
    position = position_jitter(width = 0.18, height = 0, seed = 11),
    alpha = 0.55,
    linewidth = 0.35
  ) +
  labs(
    x = "Number of simulated transmission pairs",
    y = "Estimated mean generation time, days",
    title = "Recovery of the true generation-time mean from small samples",
    subtitle = "Each vertical interval is one refit. Dashed line is the known data-generating value; intervals are 90% posterior intervals.",
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "figure_1_mean_generation_recovery_intervals.png"), p_recovery_mean, width = 8.5, height = 5.2, dpi = 300)

error_long <- plot_results %>%
  select(rep_id, n_pairs, mean_gen_abs_error, sd_gen_abs_error) %>%
  pivot_longer(cols = c(mean_gen_abs_error, sd_gen_abs_error), names_to = "quantity", values_to = "absolute_error_days") %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = N_SMALL),
    quantity = recode(quantity, mean_gen_abs_error = "Mean generation time", sd_gen_abs_error = "Generation-time SD")
  )

p_error <- ggplot(error_long, aes(x = n_pairs_f, y = absolute_error_days)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.55) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 12), alpha = 0.25, size = 1) +
  facet_wrap(~ quantity, scales = "free_y") +
  labs(x = "Number of simulated transmission pairs", y = "Absolute error, days", title = "Estimation error by sample size") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_2_absolute_error_by_sample_size.png"), p_error, width = 8.5, height = 4.8, dpi = 300)

width_long <- plot_results %>%
  select(rep_id, n_pairs, mean_gen_ci_width, sd_gen_ci_width) %>%
  pivot_longer(cols = c(mean_gen_ci_width, sd_gen_ci_width), names_to = "quantity", values_to = "interval_width_days") %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = N_SMALL),
    quantity = recode(quantity, mean_gen_ci_width = "Mean generation time", sd_gen_ci_width = "Generation-time SD")
  )

p_width <- ggplot(width_long, aes(x = n_pairs_f, y = interval_width_days)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.55) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 13), alpha = 0.25, size = 1) +
  facet_wrap(~ quantity, scales = "free_y") +
  labs(x = "Number of simulated transmission pairs", y = "90% posterior interval width, days", title = "Posterior uncertainty by sample size") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_3_posterior_uncertainty_by_sample_size.png"), p_width, width = 8.5, height = 4.8, dpi = 300)

calibration_summary <- plot_results %>%
  group_by(n_pairs) %>%
  summarise(
    coverage_mean_gen = mean(mean_gen_covered, na.rm = TRUE),
    coverage_sd_gen = mean(sd_gen_covered, na.rm = TRUE),
    convergence_ok_rate = mean(convergence_ok, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(coverage_mean_gen, coverage_sd_gen, convergence_ok_rate), names_to = "diagnostic", values_to = "rate") %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = N_SMALL),
    diagnostic = recode(diagnostic, coverage_mean_gen = "Coverage: mean GT", coverage_sd_gen = "Coverage: GT SD", convergence_ok_rate = "Convergence OK")
  )

p_calibration <- ggplot(calibration_summary, aes(x = n_pairs_f, y = rate)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ diagnostic) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of simulated transmission pairs", y = "Proportion of simulation replicates", title = "Calibration and convergence across sample sizes") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_4_coverage_and_convergence.png"), p_calibration, width = 8.5, height = 4.8, dpi = 300)

param_plot_results <- param_results %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = N_SMALL),
    variable = factor(variable, levels = PARAM_VARS)
  )

p_param_coverage <- param_plot_results %>%
  group_by(n_pairs, variable) %>%
  summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop") %>%
  mutate(n_pairs_f = factor(n_pairs, levels = N_SMALL)) %>%
  ggplot(aes(x = n_pairs_f, y = coverage)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ variable, ncol = 4) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of simulated transmission pairs", y = "Empirical coverage", title = "Mechanistic-parameter coverage") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_parameter_coverage_by_sample_size.png"), p_param_coverage, width = 10, height = 6, dpi = 300)

p_param_width <- param_plot_results %>%
  ggplot(aes(x = n_pairs_f, y = interval_width)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.55) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), alpha = 0.25, size = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  labs(x = "Number of simulated transmission pairs", y = "90% posterior interval width", title = "Mechanistic-parameter uncertainty") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "figure_parameter_interval_width_by_sample_size.png"), p_param_width, width = 10, height = 6, dpi = 300)
