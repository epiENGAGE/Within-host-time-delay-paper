# Simulation 2
setwd("~/Within-host-time-delay-framework/Code")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

# User controls
OUT_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

this_file_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
STAN_FILE <- file.path(this_file_dir, "mechanistic_obspair2.stan")

# Simulation design.
N_REP <- 5
N_PAIRS <- 200
TARGET_MEANS <- c(3, 10, 20)

SAVE_FITS <- TRUE

# MCMC settings.
CHAINS <- 4
PARALLEL_CHAINS <- 4
ITER_WARMUP <- 500
ITER_SAMPLING <- 500
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15
REFRESH <- 10
SEED <- 52001
FORCE_RECOMPILE <- FALSE

# Numerical grid settings.
GRID_BY_SI <- 1.0
GRID_BY_LU <- 1.0
GRID_BY_IP <- 0.25
N_Z_QUAD <- 5
Z_MAX_SD <- 3.5

SI_MIN <- -30
SI_MAX <- 80
LU_MAX <- 90
GT_MAX <- 2 * LU_MAX
IP_MIN <- 0.05
IP_MAX <- 70

TRUE_R_GROWTH <- 0.0

# Priors.
PRIOR_LIST <- list(
  prior_log_lambda1_mean = log(0.10),
  prior_log_lambda1_sd = 2.0,
  prior_log_lambda_I_mean = log(0.12),
  prior_log_lambda_I_sd = 2.0,
  prior_log_r_mean = log(0.20),
  prior_log_r_sd = 1.25,
  prior_log_d_mean = log(0.30),
  prior_log_d_sd = 1.25,
  prior_m_mean = 8.0,
  prior_m_sd = 6.0,
  prior_log_kappa_mean = log(0.75),
  prior_log_kappa_sd = 1.25,
  prior_log_sigma_z_mean = log(0.35),
  prior_log_sigma_z_sd = 1.0
)

PARAM_VARS <- c(
  "log_lambda1", "log_lambda_I", "log_r", "log_d", "m", "log_kappa", "log_sigma_z"
)
FIT_SUMMARY_VARS <- c(PARAM_VARS, "mean_gen", "sd_gen")

PARAM_BOUNDS <- list(
  log_lambda1 = c(-8, 8),
  log_lambda_I = c(-8, 8),
  log_r = c(-5, 2),
  log_d = c(-5, 2),
  m = c(0.25, 12),
  log_kappa = c(-5, 3),
  log_sigma_z = c(-6, 1.5)
)

clamp <- function(x, lo, hi, eps = 1e-6) {
  min(max(x, lo + eps), hi - eps)
}

clamp_init_list <- function(init) {
  for (nm in names(PARAM_BOUNDS)) {
    init[[nm]] <- clamp(init[[nm]], PARAM_BOUNDS[[nm]][1], PARAM_BOUNDS[[nm]][2])
  }
  init
}

# Helper functions
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
  # Updated manuscript version:
  #   V(t) = exp(r t) * {1 + (r/d) exp[kappa (t - m)]}^(-(r+d)/kappa)
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

make_z_quad <- function(n_z = 5, z_max_sd = 3.5) {
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
message("Grid sizes: SI=", length(G$si_grid),
        ", l=", length(G$l_grid),
        ", u=", length(G$u_grid),
        ", GT=", length(G$gt_grid),
        ", IP=", length(G$ip_grid),
        ", delta=", length(G$delta_grid),
        ", sparse IP terms=", G$N_ip_sparse,
        ", z=", length(G$z_std_grid))

# R mechanistic distribution for simulation truth
mechanistic_distribution <- function(par, r_growth, G) {
  lambda1 <- exp(par$log_lambda1)
  lambda_I <- exp(par$log_lambda_I)
  sigma_z <- exp(par$log_sigma_z)
  
  base_l <- viral_shape_vec(G$l_grid, par$log_r, par$log_d, par$m, par$log_kappa)
  base_t <- viral_shape_vec(G$t_grid, par$log_r, par$log_d, par$m, par$log_kappa)
  base_ip <- viral_shape_vec(G$ip_grid, par$log_r, par$log_d, par$m, par$log_kappa)
  
  N_l <- length(G$l_grid)
  N_u <- length(G$u_grid)
  N_gt <- length(G$gt_grid)
  N_ip <- length(G$ip_grid)
  N_z <- length(G$z_std_grid)
  N_delta <- length(G$delta_grid)
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
      for (jz in seq_len(N_z)) {
        acc_s <- acc_s + G$z_std_w[jz] * fG_mass_z[kg, jz] * C_delta_z[dd, jz]
      }
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

# Regime construction and calibration
seed_regimes <- tibble::tribble(
  ~regime, ~target_mean_gen, ~log_lambda1, ~log_lambda_I, ~log_r, ~log_d, ~m, ~log_kappa, ~log_sigma_z,
  "mean_3",   3,  log(0.16),  log(0.16), log(0.27), log(0.62),  3.5, log(1.05), log(0.35),
  "mean_10", 10,  log(0.06),  log(0.10), log(0.18), log(0.25),  7.5, log(0.70), log(0.35),
  "mean_20", 20,  log(0.025), log(0.08), log(0.12), log(0.13), 11.0, log(0.45), log(0.35)
)

row_to_par <- function(row) {
  as.list(row[, PARAM_VARS, drop = TRUE])
}

shift_hazards <- function(par, shift) {
  par$log_lambda1 <- par$log_lambda1 + shift
  par
}

calibrate_to_target_mean <- function(par, target, G, r_growth = 0, search = c(-8, 8)) {
  mean_for_shift <- function(shift) {
    pp <- shift_hazards(par, shift)
    mechanistic_distribution(pp, r_growth, G)$mean_gen
  }
  
  f <- function(shift) mean_for_shift(shift) - target
  lo <- search[1]
  hi <- search[2]
  flo <- f(lo)
  fhi <- f(hi)
  
  if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) {
    root <- uniroot(f, interval = c(lo, hi), tol = 0.02)$root
    method <- "uniroot"
  } else {
    opt <- optimize(function(s) abs(f(s)), interval = c(lo, hi))
    root <- opt$minimum
    method <- "closest_on_search_interval"
    warning(
      "Could not exactly bracket target mean ", target,
      ". Using closest hazard shift on search interval. Achieved mean may differ from target."
    )
  }
  
  pp <- shift_hazards(par, root)
  dist <- mechanistic_distribution(pp, r_growth, G)
  list(par = pp, dist = dist, shift = root, method = method)
}

message("Calibrating regimes to target mean GTs...")
regime_objects <- vector("list", nrow(seed_regimes))

for (i in seq_len(nrow(seed_regimes))) {
  row <- seed_regimes[i, ]
  par0 <- row_to_par(row)
  cal <- calibrate_to_target_mean(
    par = par0,
    target = row$target_mean_gen,
    G = G,
    r_growth = TRUE_R_GROWTH
  )
  regime_objects[[i]] <- list(
    regime = row$regime,
    target_mean_gen = row$target_mean_gen,
    par = cal$par,
    dist = cal$dist,
    hazard_shift = cal$shift,
    calibration_method = cal$method
  )
  
  message(
    "  ", row$regime,
    ": target=", row$target_mean_gen,
    ", achieved mean=", round(cal$dist$mean_gen, 3),
    ", achieved sd=", round(cal$dist$sd_gen, 3),
    ", hazard shift=", round(cal$shift, 3),
    " [", cal$method, "]"
  )
}

true_targets <- purrr::map_dfr(regime_objects, function(x) {
  tibble(
    regime = x$regime,
    target_mean_gen = x$target_mean_gen,
    true_mean_gen = x$dist$mean_gen,
    true_sd_gen = x$dist$sd_gen,
    hazard_shift = x$hazard_shift,
    calibration_method = x$calibration_method,
    log_lambda1 = x$par$log_lambda1,
    log_lambda_I = x$par$log_lambda_I,
    log_r = x$par$log_r,
    log_d = x$par$log_d,
    m = x$par$m,
    log_kappa = x$par$log_kappa,
    log_sigma_z = x$par$log_sigma_z
  )
})

print(true_targets)
write_csv(true_targets, file.path(OUT_DIR, "true_regime_targets.csv"))

# Fail early if calibration generated a true parameter outside Stan support.
for (nm in names(PARAM_BOUNDS)) {
  vals <- true_targets[[nm]]
  bad <- vals <= PARAM_BOUNDS[[nm]][1] | vals >= PARAM_BOUNDS[[nm]][2]
  if (any(bad)) {
    stop(
      "Calibrated parameter ", nm, " is outside the Stan bounds for regime(s): ",
      paste(true_targets$regime[bad], collapse = ", "),
      ". Widen PARAM_BOUNDS and the Stan parameter declaration, or change the regime calibration."
    )
  }
}

# Stan data and fitting helpers
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
        log_lambda1 = PRIOR_LIST$prior_log_lambda1_mean,
        log_lambda_I = PRIOR_LIST$prior_log_lambda_I_mean,
        log_r = PRIOR_LIST$prior_log_r_mean,
        log_d = PRIOR_LIST$prior_log_d_mean,
        m = PRIOR_LIST$prior_m_mean,
        log_kappa = PRIOR_LIST$prior_log_kappa_mean,
        log_sigma_z = PRIOR_LIST$prior_log_sigma_z_mean
      )
    } else true_par
    
    clamp_init_list(list(
      log_lambda1 = base$log_lambda1 + rnorm(1, 0, jitter_sd),
      log_lambda_I = base$log_lambda_I + rnorm(1, 0, jitter_sd),
      log_r = base$log_r + rnorm(1, 0, jitter_sd),
      log_d = base$log_d + rnorm(1, 0, jitter_sd),
      m = base$m + rnorm(1, 0, 0.25),
      log_kappa = base$log_kappa + rnorm(1, 0, jitter_sd),
      log_sigma_z = base$log_sigma_z + rnorm(1, 0, jitter_sd)
    ))
  }
}

summarise_fit <- function(fit, regime, target_mean_gen, n_pairs, rep_id, truth, true_par) {
  s <- fit$summary(variables = FIT_SUMMARY_VARS)
  mean_row <- s %>% filter(variable == "mean_gen")
  sd_row <- s %>% filter(variable == "sd_gen")
  
  param_summ <- s %>%
    filter(variable %in% PARAM_VARS) %>%
    select(variable, median, q5, q95, rhat, ess_bulk, ess_tail) %>%
    mutate(true_value = as.numeric(unlist(true_par))[match(variable, names(true_par))])
  
  tibble(
    regime = regime,
    target_mean_gen = target_mean_gen,
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
    n_divergent = sum(fit$sampler_diagnostics(format = "df")$divergent__, na.rm = TRUE),
    # Useful coarse check for whether raw parameters are identifiable in each regime.
    median_abs_param_error = median(abs(param_summ$median - param_summ$true_value), na.rm = TRUE)
  )
}

# Compile model
mod <- cmdstan_model(STAN_FILE, force_recompile = FORCE_RECOMPILE)

# Simulation loop
set.seed(SEED)
fit_rows <- list()
failed_rows <- list()
row_id <- 1L
fail_id <- 1L

total_jobs <- length(regime_objects) * N_REP
job_id <- 0L

for (reg_obj in regime_objects) {
  regime <- reg_obj$regime
  target_mean <- reg_obj$target_mean_gen
  true_par <- reg_obj$par
  true_dist <- reg_obj$dist
  
  for (rep_id in seq_len(N_REP)) {
    job_id <- job_id + 1L
    message(
      "\nSimulation ", job_id, "/", total_jobs,
      ": regime=", regime,
      ", target mean=", target_mean,
      ", n=", N_PAIRS,
      ", replicate=", rep_id, "/", N_REP
    )
    
    counts_full <- as.vector(rmultinom(1, size = N_PAIRS, prob = true_dist$p_si))
    stan_data <- make_stan_data(counts_full, TRUE_R_GROWTH, G)
    
    fit_result <- tryCatch({
      fit <- mod$sample(
        data = stan_data,
        chains = CHAINS,
        parallel_chains = PARALLEL_CHAINS,
        iter_warmup = ITER_WARMUP,
        iter_sampling = ITER_SAMPLING,
        seed = SEED + 100000 * rep_id + round(1000 * target_mean),
        init = make_inits(true_par, jitter_sd = 0.15),
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
      failed_rows[[fail_id]] <- tibble(
        regime = regime,
        target_mean_gen = target_mean,
        n_pairs = N_PAIRS,
        rep_id = rep_id,
        error = fit_result$error
      )
      fail_id <- fail_id + 1L
      write_csv(bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
      next
    }
    
    fit <- fit_result$fit
    
    fit_rows[[row_id]] <- summarise_fit(
      fit = fit,
      regime = regime,
      target_mean_gen = target_mean,
      n_pairs = N_PAIRS,
      rep_id = rep_id,
      truth = true_dist,
      true_par = true_par
    )
    row_id <- row_id + 1L
    
    if (SAVE_FITS) {
      fit$save_object(file.path(OUT_DIR, sprintf("fit_%s_n%s_rep%s.rds", regime, N_PAIRS, rep_id)))
    }
    
    # Incremental writes so long runs can be stopped without losing progress.
    if (length(fit_rows) > 0) {
      write_csv(bind_rows(fit_rows), file.path(OUT_DIR, "simulation_fit_level_results.csv"))
    }
    if (length(failed_rows) > 0) {
      write_csv(bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
    }
  }
}

fit_results <- bind_rows(fit_rows)
failures <- bind_rows(failed_rows)

write_csv(fit_results, file.path(OUT_DIR, "simulation_fit_level_results.csv"))
if (nrow(failures) > 0) write_csv(failures, file.path(OUT_DIR, "simulation_failures.csv"))

# Summary
recovery_summary <- fit_results %>%
  group_by(regime, target_mean_gen, true_mean_gen, true_sd_gen, n_pairs) %>%
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
    convergence_ok_rate = mean(max_rhat < 1.05 & n_divergent == 0, na.rm = TRUE),
    median_min_ess_bulk = median(min_ess_bulk, na.rm = TRUE),
    median_min_ess_tail = median(min_ess_tail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(target_mean_gen)

write_csv(recovery_summary, file.path(OUT_DIR, "paper_recovery_summary_by_regime.csv"))
print(recovery_summary)

# Figures
plot_results <- fit_results %>%
  mutate(
    regime_label = paste0("Target mean ", target_mean_gen, " d\nTrue mean ", round(true_mean_gen, 2), " d"),
    regime_label = factor(regime_label, levels = unique(regime_label[order(target_mean_gen)])),
    mean_gen_covered_f = ifelse(mean_gen_covered, "90% interval covers truth", "misses truth"),
    sd_gen_covered_f = ifelse(sd_gen_covered, "90% interval covers truth", "misses truth")
  )

# Main regime-recovery figure
p_recovery_mean <- ggplot(
  plot_results,
  aes(x = regime_label, y = mean_gen_median, ymin = mean_gen_q05, ymax = mean_gen_q95)
) +
  geom_pointrange(
    aes(shape = mean_gen_covered_f),
    position = position_jitter(width = 0.18, height = 0, seed = 21),
    alpha = 0.55,
    linewidth = 0.35
  ) +
  geom_point(
    data = distinct(plot_results, regime_label, true_mean_gen),
    aes(x = regime_label, y = true_mean_gen),
    inherit.aes = FALSE,
    shape = 95,
    size = 12,
    linewidth = 1.2
  ) +
  labs(
    x = "Data-generating disease regime",
    y = "Estimated mean generation time, days",
    title = "Generation-time mean recovery across disease regimes",
    subtitle = paste0("Each interval is one refit to n = ", N_PAIRS, " simulated serial intervals. Horizontal ticks mark true values."),
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUT_DIR, "figure_1_mean_gt_recovery_by_regime.png"),
  p_recovery_mean,
  width = 9,
  height = 5.2,
  dpi = 300
)

# Absolute error figure
error_long <- plot_results %>%
  select(regime_label, target_mean_gen, rep_id, mean_gen_abs_error, sd_gen_abs_error) %>%
  pivot_longer(
    cols = c(mean_gen_abs_error, sd_gen_abs_error),
    names_to = "quantity",
    values_to = "absolute_error_days"
  ) %>%
  mutate(
    quantity = recode(
      quantity,
      mean_gen_abs_error = "Mean generation time",
      sd_gen_abs_error = "Generation-time SD"
    )
  )

p_error <- ggplot(error_long, aes(x = regime_label, y = absolute_error_days)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.55) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 22),
    alpha = 0.25,
    size = 1
  ) +
  facet_wrap(~ quantity, scales = "free_y") +
  labs(
    x = "Data-generating disease regime",
    y = "Absolute error, days",
    title = "GT recovery error across short, medium, and long regimes",
    subtitle = paste0("Each point is one simulated dataset with n = ", N_PAIRS, ". Lower is better.")
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(OUT_DIR, "figure_2_absolute_error_by_regime.png"),
  p_error,
  width = 9,
  height = 4.8,
  dpi = 300
)

# Calibration and convergence figure
calibration_summary <- plot_results %>%
  group_by(regime_label, target_mean_gen) %>%
  summarise(
    coverage_mean_gen = mean(mean_gen_covered, na.rm = TRUE),
    coverage_sd_gen = mean(sd_gen_covered, na.rm = TRUE),
    convergence_ok_rate = mean(max_rhat < 1.05 & n_divergent == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(coverage_mean_gen, coverage_sd_gen, convergence_ok_rate),
    names_to = "diagnostic",
    values_to = "rate"
  ) %>%
  mutate(
    diagnostic = recode(
      diagnostic,
      coverage_mean_gen = "Coverage: mean GT",
      coverage_sd_gen = "Coverage: GT SD",
      convergence_ok_rate = "Convergence OK"
    )
  )

p_calibration <- ggplot(calibration_summary, aes(x = regime_label, y = rate)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ diagnostic) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Data-generating disease regime",
    y = "Proportion of simulation replicates",
    title = "Calibration and convergence across disease regimes",
    subtitle = "Dashed line marks the nominal 90% level for posterior interval coverage."
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(OUT_DIR, "figure_3_coverage_and_convergence_by_regime.png"),
  p_calibration,
  width = 9,
  height = 4.8,
  dpi = 300
)
