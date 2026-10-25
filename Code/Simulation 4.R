# Simulation 4: joint recovery of epidemiological delays
# ------------------------------------------------------
# Scientific question:
#   If serial-interval data are generated from the within-host-informed
#   observed-pair model, can the same fitted model simultaneously recover
#   the latent-period (LP), incubation-period (IP), generation-time (GT),
#   and serial-interval (SI) distributions?
#
# This is adapted directly from Simulation 1. The DGM, priors, grids,
# quadrature and MCMC settings are kept the same; the sample size is n = 150.
# The main change is that all four delay distributions are treated as targets.
# This OPTIMIZED version preserves N_REP=100, n=150, the priors, grids,
# quadrature, chains and MCMC iterations. Speedups come from changing the
# implementation, not weakening the simulation.
#
# Outputs:
#   1) posterior recovery of mean and SD for LP/IP/GT/SI across replicates;
#   2) CRPS, bias, RMSE, interval width and 90% coverage for each target;
#   3) full-distribution recovery using posterior-median parameter fits;
#   4) figures showing simultaneous recovery of all four delays.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# -----------------------------------------------------------------------------
# Paths: edit PROJECT_DIR if your local clone lives somewhere else.
# -----------------------------------------------------------------------------
PROJECT_DIR <- path.expand("~/Within-host-time-delay-framework")
CODE_DIR <- file.path(PROJECT_DIR, "Code")

PROJECT_DIR <- getwd()
CODE_DIR <- getwd()
OUT_DIR <- file.path(PROJECT_DIR, "Results", "Simulations", "Simulation 4")
STAN_FILE <- file.path(CODE_DIR, "sim4_joint_delay_optimized.stan")
GQ_STAN_FILE <- file.path(CODE_DIR, "sim4_joint_delay_gq_optimized.stan")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Main simulation setup
# -----------------------------------------------------------------------------
N_REP <- 100
N_PAIRS <- 150

CHAINS <- 4
PARALLEL_CHAINS <- 4
ITER_WARMUP <- 250
ITER_SAMPLING <- 250
ADAPT_DELTA <- 0.99
MAX_TREEDEPTH <- 15
REFRESH <- 0
SEED <- 44001
SAVE_FITS <- FALSE
FORCE_RECOMPILE <- FALSE
WRITE_EVERY <- 5L  # checkpoint every 5 replicates; does not change inference

# -----------------------------------------------------------------------------
# Numerical grid settings -- identical to Simulation 1
# -----------------------------------------------------------------------------
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

# Same DGM as Simulation 1.
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

# These are produced by sim4_joint_delay.stan.
DELAY_SUMMARY_VARS <- c(
  "mean_lp", "sd_lp",
  "mean_ip", "sd_ip",
  "mean_gen", "sd_gen",
  "mean_si", "sd_si"
)

DELAY_LABELS <- c(
  lp = "Latent period",
  ip = "Incubation period",
  gen = "Generation time",
  si = "Serial interval"
)

# Same priors as Simulation 1.
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

# -----------------------------------------------------------------------------
# Numerical helpers
# -----------------------------------------------------------------------------
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
    for (i in 2:length(y)) {
      out[i] <- out[i - 1] + dx * 0.5 * (y[i - 1] + y[i])
    }
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
  stopifnot(abs((GRID_BY_SI / GRID_BY_IP) -
                  round(GRID_BY_SI / GRID_BY_IP)) < 1e-12)

  # Because l/u/t/GT all use the same spacing, the conditional U vector and
  # its GT destination are contiguous slices. We only need their start cells.
  lu_start_idx <- find_grid_index(l_grid, t_grid, GRID_BY_LU)
  gt_start_idx <- find_grid_index(l_grid, gt_grid, GRID_BY_LU)

  # Compact aligned-slice representation of IP2 = IP1 + delta.
  # This replaces the old ~O(N_ip^2) sparse scatter lookup with one
  # dot-product per delta cell.
  ip_i1_start <- integer(length(delta_grid))
  ip_i2_start <- integer(length(delta_grid))
  ip_overlap_len <- integer(length(delta_grid))

  for (dd in seq_along(delta_grid)) {
    shift <- as.integer(round(delta_grid[dd] / GRID_BY_IP))

    if (shift >= 0L) {
      i1 <- 1L
      i2 <- 1L + shift
      L <- length(ip_grid) - shift
    } else {
      i1 <- 1L - shift
      i2 <- 1L
      L <- length(ip_grid) + shift
    }

    if (L > 0L) {
      ip_i1_start[dd] <- i1
      ip_i2_start[dd] <- i2
      ip_overlap_len[dd] <- L
    } else {
      ip_i1_start[dd] <- 1L
      ip_i2_start[dd] <- 1L
      ip_overlap_len[dd] <- 0L
    }
  }

  # For a fixed SI, delta = SI - GT decreases one delta-grid cell whenever
  # GT increases one cell. Therefore all required delta values are one
  # contiguous block, read in reverse order.
  si_delta_low_idx <- integer(length(si_grid))
  for (ks in seq_along(si_grid)) {
    low_delta <- si_grid[ks] - max(gt_grid)
    si_delta_low_idx[ks] <- find_grid_index(
      low_delta, delta_grid, GRID_BY_SI
    )
  }

  zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

  list(
    si_grid = si_grid,
    l_grid = l_grid,
    u_grid = u_grid,
    gt_grid = gt_grid,
    ip_grid = ip_grid,
    t_grid = t_grid,
    delta_grid = delta_grid,
    lu_start_idx = as.integer(lu_start_idx),
    gt_start_idx = as.integer(gt_start_idx),
    ip_i1_start = as.integer(ip_i1_start),
    ip_i2_start = as.integer(ip_i2_start),
    ip_overlap_len = as.integer(ip_overlap_len),
    si_delta_low_idx = as.integer(si_delta_low_idx),
    z_std_grid = zq$x,
    z_std_w = zq$w
  )
}

G <- make_grids_and_lookups()

# -----------------------------------------------------------------------------
# Mechanistic distribution calculator
#
# This is the Simulation-1 truth calculator extended to return all four delay
# distributions, not only generation time and SI.
# -----------------------------------------------------------------------------
mechanistic_distribution <- function(par, r_growth, G) {
  lambda1 <- exp(par$log_lambda1)
  lambda_I <- exp(par$log_lambda_I)
  sigma_z <- exp(par$log_sigma_z)

  base_l <- viral_shape_vec(
    G$l_grid, par$log_r, par$log_d, par$m, par$log_kappa
  )
  base_t <- viral_shape_vec(
    G$t_grid, par$log_r, par$log_d, par$m, par$log_kappa
  )
  base_ip <- viral_shape_vec(
    G$ip_grid, par$log_r, par$log_d, par$m, par$log_kappa
  )

  N_l <- length(G$l_grid)
  N_u <- length(G$u_grid)
  N_gt <- length(G$gt_grid)
  N_ip <- length(G$ip_grid)
  N_z <- length(G$z_std_grid)
  N_delta <- length(G$delta_grid)
  N_si <- length(G$si_grid)

  fL_z <- matrix(0, nrow = N_l, ncol = N_z)
  fG_mass_z <- matrix(0, nrow = N_gt, ncol = N_z)
  fIP_z <- matrix(0, nrow = N_ip, ncol = N_z)

  for (jz in seq_len(N_z)) {
    ez <- exp(sigma_z * G$z_std_grid[jz])

    # Latent period.
    hL <- ez * lambda1 * base_l
    HL <- cumulative_trapezoid(hL, GRID_BY_LU)
    normL <- normalize_density(hL * exp(-HL), GRID_BY_LU)
    fL_z[, jz] <- normL

    # Generation time. Each U trajectory and each GT update are contiguous
    # slices, so there is no scalar ju scatter loop.
    gcol <- numeric(N_gt)
    for (jl in seq_len(N_l)) {
      t0 <- G$lu_start_idx[jl]
      g0 <- G$gt_start_idx[jl]
      rawU <- base_t[t0:(t0 + N_u - 1L)]
      normU <- normalize_density(rawU, GRID_BY_LU)
      idx <- g0:(g0 + N_u - 1L)
      gcol[idx] <- gcol[idx] +
        normL[jl] * normU * GRID_BY_LU * GRID_BY_LU
    }
    fG_mass_z[, jz] <- normalize_mass(gcol)

    # Incubation period.
    hI <- ez * lambda_I * base_ip
    HI <- cumulative_trapezoid(hI, GRID_BY_IP)
    fIP_z[, jz] <- normalize_density(hI * exp(-HI), GRID_BY_IP)
  }

  # Marginal LP/IP/GT.
  fLP_density <- normalize_density(
    as.vector(fL_z %*% G$z_std_w), GRID_BY_LU
  )
  fLP_mass <- normalize_mass(fLP_density * GRID_BY_LU)

  fIP_density <- normalize_density(
    as.vector(fIP_z %*% G$z_std_w), GRID_BY_IP
  )
  fIP_mass <- normalize_mass(fIP_density * GRID_BY_IP)

  fG_mass <- normalize_mass(as.vector(fG_mass_z %*% G$z_std_w))

  # Difference in incubation periods. One vector dot product per delta/z
  # replaces the old scalar sparse scatter loop.
  C_delta_z <- matrix(0, nrow = N_delta, ncol = N_z)
  tilt <- exp(-r_growth * G$ip_grid)

  for (jz in seq_len(N_z)) {
    weighted_ip1 <- tilt * fIP_z[, jz]

    for (dd in seq_len(N_delta)) {
      L <- G$ip_overlap_len[dd]
      if (L > 0L) {
        i1 <- G$ip_i1_start[dd]
        i2 <- G$ip_i2_start[dd]

        C_delta_z[dd, jz] <- sum(
          weighted_ip1[i1:(i1 + L - 1L)] *
            fIP_density[i2:(i2 + L - 1L)]
        ) * GRID_BY_IP
      }
    }
  }

  # SI = GT + delta. For each SI, use a reversed contiguous C_delta block.
  unnorm_si <- numeric(N_si)
  for (ks in seq_len(N_si)) {
    dlow <- G$si_delta_low_idx[ks]
    cblock <- C_delta_z[
      dlow:(dlow + N_gt - 1L), , drop = FALSE
    ]
    cblock <- cblock[N_gt:1L, , drop = FALSE]

    per_z <- colSums(fG_mass_z * cblock)
    unnorm_si[ks] <- max(sum(G$z_std_w * per_z), 1e-300)
  }
  p_si <- normalize_mass(unnorm_si)

  moment <- function(grid, mass) {
    mass <- normalize_mass(mass)
    mu <- sum(grid * mass)
    s <- sqrt(sum((grid - mu)^2 * mass))
    c(mean = mu, sd = s)
  }

  lp_m <- moment(G$l_grid, fLP_mass)
  ip_m <- moment(G$ip_grid, fIP_mass)
  gt_m <- moment(G$gt_grid, fG_mass)
  si_m <- moment(G$si_grid, p_si)

  list(
    lp_mass = fLP_mass,
    ip_mass = fIP_mass,
    gt_mass = fG_mass,
    si_mass = p_si,
    lp_density = fLP_density,
    ip_density = fIP_density,
    gt_density = fG_mass / GRID_BY_LU,
    si_density = p_si / GRID_BY_SI,
    mean_lp = unname(lp_m["mean"]),
    sd_lp = unname(lp_m["sd"]),
    mean_ip = unname(ip_m["mean"]),
    sd_ip = unname(ip_m["sd"]),
    mean_gen = unname(gt_m["mean"]),
    sd_gen = unname(gt_m["sd"]),
    mean_si = unname(si_m["mean"]),
    sd_si = unname(si_m["sd"])
  )
}

# -----------------------------------------------------------------------------
# True model-implied joint delay targets
# -----------------------------------------------------------------------------
TRUE_DIST <- mechanistic_distribution(TRUE_PARAMETERS, TRUE_R_GROWTH, G)

TRUE_VALUES <- c(
  mean_lp = TRUE_DIST$mean_lp,
  sd_lp = TRUE_DIST$sd_lp,
  mean_ip = TRUE_DIST$mean_ip,
  sd_ip = TRUE_DIST$sd_ip,
  mean_gen = TRUE_DIST$mean_gen,
  sd_gen = TRUE_DIST$sd_gen,
  mean_si = TRUE_DIST$mean_si,
  sd_si = TRUE_DIST$sd_si
)

TRUE_SUMMARY <- tibble(
  variable = names(TRUE_VALUES),
  true_value = as.numeric(TRUE_VALUES)
) %>%
  dplyr::mutate(
    moment = dplyr::if_else(stringr::str_starts(variable, "mean_"), "Mean", "SD"),
    delay_code = stringr::str_remove(variable, "^(mean_|sd_)"),
    delay = dplyr::recode(delay_code, !!!DELAY_LABELS)
  )

print(TRUE_SUMMARY)
readr::write_csv(TRUE_SUMMARY, file.path(OUT_DIR, "true_joint_delay_targets.csv"))

# -----------------------------------------------------------------------------
# Stan data builder
# -----------------------------------------------------------------------------
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

      N_delta = length(G$delta_grid),
      delta_grid = as.numeric(G$delta_grid),

      lu_start_idx = as.integer(G$lu_start_idx),
      gt_start_idx = as.integer(G$gt_start_idx),

      ip_i1_start = as.integer(G$ip_i1_start),
      ip_i2_start = as.integer(G$ip_i2_start),
      ip_overlap_len = as.integer(G$ip_overlap_len),

      si_delta_low_idx = as.integer(G$si_delta_low_idx),

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
        log_lambda1 = log(0.15),
        log_lambda_I = log(0.15),
        log_r = log(0.25),
        log_d = log(0.60),
        m = 3.5,
        log_kappa = 0,
        log_sigma_z = log(0.35)
      )
    } else {
      true_par
    }

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

# -----------------------------------------------------------------------------
# Recovery metrics
# All dplyr verbs below are namespace-qualified so this script is safe even if
# packages such as MASS/plyr are attached in the calling R session.
# -----------------------------------------------------------------------------
crps_empirical <- function(x, truth) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n == 0) return(NA_real_)
  # CRPS = E|X-y| - 1/2 E|X-X'|, evaluated from posterior draws.
  mean(abs(x - truth)) - sum((2 * seq_len(n) - n - 1) * x) / (n^2)
}

summarise_delay_recovery <- function(fit, rep_id) {
  s <- fit$summary(variables = DELAY_SUMMARY_VARS)
  draws <- as_draws_df(fit$draws(variables = DELAY_SUMMARY_VARS))

  out <- lapply(DELAY_SUMMARY_VARS, function(v) {
    row <- s %>% dplyr::filter(variable == v)
    truth <- unname(TRUE_VALUES[v])

    tibble(
      n_pairs = N_PAIRS,
      rep_id = rep_id,
      variable = v,
      true_value = truth,
      median = row$median,
      q05 = row$q5,
      q95 = row$q95,
      error = row$median - truth,
      abs_error = abs(row$median - truth),
      interval_width = row$q95 - row$q5,
      covered = row$q5 <= truth & row$q95 >= truth,
      crps = crps_empirical(draws[[v]], truth)
    )
  })

  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      moment = dplyr::if_else(stringr::str_starts(variable, "mean_"), "Mean", "SD"),
      delay_code = stringr::str_remove(variable, "^(mean_|sd_)"),
      delay = dplyr::recode(delay_code, !!!DELAY_LABELS)
    )
}

summarise_diagnostics <- function(fit, rep_id) {
  s <- fit$summary(variables = PARAM_VARS)
  tibble(
    n_pairs = N_PAIRS,
    rep_id = rep_id,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(fit$sampler_diagnostics(format = "df")$divergent__, na.rm = TRUE)
  )
}

posterior_median_parameters <- function(fit) {
  s <- fit$summary(variables = PARAM_VARS)
  vals <- setNames(s$median, s$variable)
  as.list(vals[PARAM_VARS])
}

safe_kl <- function(p, q, eps = 1e-12) {
  p <- normalize_mass(p)
  q <- normalize_mass(q)
  sum(p * log((p + eps) / (q + eps)))
}

distribution_metrics <- function(fitted_dist, rep_id) {
  specs <- list(
    lp = list(true = TRUE_DIST$lp_mass, fit = fitted_dist$lp_mass),
    ip = list(true = TRUE_DIST$ip_mass, fit = fitted_dist$ip_mass),
    gen = list(true = TRUE_DIST$gt_mass, fit = fitted_dist$gt_mass),
    si = list(true = TRUE_DIST$si_mass, fit = fitted_dist$si_mass)
  )

  dplyr::bind_rows(lapply(names(specs), function(code) {
    p <- normalize_mass(specs[[code]]$true)
    q <- normalize_mass(specs[[code]]$fit)
    tibble(
      n_pairs = N_PAIRS,
      rep_id = rep_id,
      delay_code = code,
      delay = unname(DELAY_LABELS[code]),
      total_variation = 0.5 * sum(abs(p - q)),
      mass_rmse = sqrt(mean((p - q)^2)),
      kl_true_to_fit = safe_kl(p, q)
    )
  }))
}

distribution_long <- function(dist, rep_id) {
  dplyr::bind_rows(
    tibble(rep_id = rep_id, delay_code = "lp", delay = DELAY_LABELS[["lp"]],
           time = G$l_grid, density = dist$lp_density),
    tibble(rep_id = rep_id, delay_code = "ip", delay = DELAY_LABELS[["ip"]],
           time = G$ip_grid, density = dist$ip_density),
    tibble(rep_id = rep_id, delay_code = "gen", delay = DELAY_LABELS[["gen"]],
           time = G$gt_grid, density = dist$gt_density),
    tibble(rep_id = rep_id, delay_code = "si", delay = DELAY_LABELS[["si"]],
           time = G$si_grid, density = dist$si_density)
  )
}

TRUE_DENSITIES <- distribution_long(TRUE_DIST, rep_id = 0L) %>%
  dplyr::select(-rep_id) %>%
  dplyr::rename(true_density = density)

stopifnot(all(c("delay_code", "delay", "time", "true_density") %in% names(TRUE_DENSITIES)))

# -----------------------------------------------------------------------------
# Compile optimized Stan models
#
# Sampling model:
#   - likelihood only (no expensive generated quantities during HMC)
#   - stanc O1 compiler optimizations
#
# Standalone GQ model:
#   - reconstructs LP/IP/GT/SI after sampling
#   - runs without autodiff/gradient construction
#
# Note: this optimized formulation deliberately does NOT force reduce_sum().
# The expensive operations are shared convolutions, not conditionally
# independent likelihood terms; reduce_sum would duplicate those calculations.
# Parallelism is therefore used across chains and standalone-GQ chains, where
# it is useful, while the within-chain speedup comes from eliminating the large
# scalar autodiff loops.
# -----------------------------------------------------------------------------
for (f in c(STAN_FILE, GQ_STAN_FILE)) {
  if (!file.exists(f)) {
    stop(
      "Cannot find ", f, "\n",
      "Put both optimized Stan files in your Code directory or edit the paths above."
    )
  }
}

Sys.setenv(
  MAKEFLAGS = paste0(
    "-j", max(1L, parallel::detectCores(logical = TRUE) - 1L)
  )
)

message("Compiling optimized fitting model with stanc O1...")
mod <- cmdstan_model(
  STAN_FILE,
  stanc_options = list("O1"),
  force_recompile = FORCE_RECOMPILE
)

message("Compiling standalone generated-quantities model with stanc O1...")
gq_mod <- cmdstan_model(
  GQ_STAN_FILE,
  stanc_options = list("O1"),
  force_recompile = FORCE_RECOMPILE
)

# -----------------------------------------------------------------------------
# Simulation loop
# -----------------------------------------------------------------------------
set.seed(SEED)

delay_rows <- list()
diag_rows <- list()
dist_metric_rows <- list()
density_rows <- list()
failed_rows <- list()

delay_id <- 1L
diag_id <- 1L
dist_id <- 1L
density_id <- 1L
fail_id <- 1L

for (rep_id in seq_len(N_REP)) {
  message("\nSimulation 4: n=", N_PAIRS, ", replicate=", rep_id, "/", N_REP)

  # Only SI is supplied to the fit, exactly as in Simulation 1.
  counts_full <- as.vector(rmultinom(1, size = N_PAIRS, prob = TRUE_DIST$si_mass))
  stan_data <- make_stan_data(counts_full, TRUE_R_GROWTH, G)

  fit_result <- tryCatch({
    fit <- mod$sample(
      data = stan_data,
      chains = CHAINS,
      parallel_chains = PARALLEL_CHAINS,
      iter_warmup = ITER_WARMUP,
      iter_sampling = ITER_SAMPLING,
      seed = SEED + 100000 * rep_id + N_PAIRS,
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
    failed_rows[[fail_id]] <- tibble(
      n_pairs = N_PAIRS,
      rep_id = rep_id,
      error = fit_result$error
    )
    fail_id <- fail_id + 1L
    readr::write_csv(dplyr::bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
    next
  }

  fit <- fit_result$fit

  # Expensive delay reconstruction is deliberately outside HMC. Standalone
  # generated quantities uses the exact posterior draws but does not build an
  # autodiff graph or enter the leapfrog trajectory.
  gq_result <- tryCatch({
    gq_fit <- gq_mod$generate_quantities(
      fitted_params = fit,
      data = stan_data,
      seed = SEED + 900000 + rep_id,
      parallel_chains = PARALLEL_CHAINS
    )
    list(ok = TRUE, fit = gq_fit, error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, fit = NULL, error = conditionMessage(e))
  })

  if (!gq_result$ok) {
    failed_rows[[fail_id]] <- tibble(
      n_pairs = N_PAIRS,
      rep_id = rep_id,
      error = paste0("Standalone GQ failed: ", gq_result$error)
    )
    fail_id <- fail_id + 1L
    next
  }

  gq_fit <- gq_result$fit

  # Posterior uncertainty for all 8 scalar delay targets.
  delay_rows[[delay_id]] <- summarise_delay_recovery(gq_fit, rep_id)
  delay_id <- delay_id + 1L

  # Standard MCMC diagnostics.
  diag_rows[[diag_id]] <- summarise_diagnostics(fit, rep_id)
  diag_id <- diag_id + 1L

  # Full fitted delay distributions evaluated at the posterior median parameters.
  med_par <- posterior_median_parameters(fit)
  med_dist <- mechanistic_distribution(med_par, TRUE_R_GROWTH, G)

  dist_metric_rows[[dist_id]] <- distribution_metrics(med_dist, rep_id)
  dist_id <- dist_id + 1L

  density_rows[[density_id]] <- distribution_long(med_dist, rep_id)
  density_id <- density_id + 1L

  if (SAVE_FITS) {
    fit$save_object(file.path(OUT_DIR, sprintf("fit_n%s_rep%s.rds", N_PAIRS, rep_id)))
  }

  # Checkpoint periodically; writing four growing CSVs every replicate can
  # become surprisingly expensive and has no inferential benefit.
  if (rep_id %% WRITE_EVERY == 0L || rep_id == N_REP) {
    readr::write_csv(dplyr::bind_rows(delay_rows), file.path(OUT_DIR, "joint_delay_recovery_by_replicate.csv"))
    readr::write_csv(dplyr::bind_rows(diag_rows), file.path(OUT_DIR, "simulation_diagnostics.csv"))
    readr::write_csv(dplyr::bind_rows(dist_metric_rows), file.path(OUT_DIR, "distribution_recovery_by_replicate.csv"))
    readr::write_csv(dplyr::bind_rows(density_rows), file.path(OUT_DIR, "fitted_delay_densities_by_replicate.csv"))
    if (length(failed_rows) > 0) {
      readr::write_csv(dplyr::bind_rows(failed_rows), file.path(OUT_DIR, "simulation_failures.csv"))
    }
  }
}

# -----------------------------------------------------------------------------
# Aggregate results
# -----------------------------------------------------------------------------
delay_results <- dplyr::bind_rows(delay_rows)
diagnostics <- dplyr::bind_rows(diag_rows)
dist_metrics <- dplyr::bind_rows(dist_metric_rows)
fitted_densities <- dplyr::bind_rows(density_rows)
failures <- dplyr::bind_rows(failed_rows)

readr::write_csv(delay_results, file.path(OUT_DIR, "joint_delay_recovery_by_replicate.csv"))
readr::write_csv(diagnostics, file.path(OUT_DIR, "simulation_diagnostics.csv"))
readr::write_csv(dist_metrics, file.path(OUT_DIR, "distribution_recovery_by_replicate.csv"))
readr::write_csv(fitted_densities, file.path(OUT_DIR, "fitted_delay_densities_by_replicate.csv"))
if (nrow(failures) > 0) {
  readr::write_csv(failures, file.path(OUT_DIR, "simulation_failures.csv"))
}

recovery_summary <- delay_results %>%
  dplyr::group_by(delay_code, delay, moment, variable) %>%
  dplyr::summarise(
    n_successful_fits = dplyr::n(),
    true_value = dplyr::first(true_value),
    median_bias = median(error, na.rm = TRUE),
    median_abs_error = median(abs_error, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mean_crps = mean(crps, na.rm = TRUE),
    median_interval_width = median(interval_width, na.rm = TRUE),
    coverage_90 = mean(covered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(delay_code, levels = c("lp", "ip", "gen", "si")), moment)

readr::write_csv(recovery_summary, file.path(OUT_DIR, "paper_joint_delay_recovery_summary.csv"))
print(recovery_summary)

distribution_summary <- dist_metrics %>%
  dplyr::group_by(delay_code, delay) %>%
  dplyr::summarise(
    n_successful_fits = dplyr::n(),
    median_total_variation = median(total_variation, na.rm = TRUE),
    median_mass_rmse = median(mass_rmse, na.rm = TRUE),
    median_kl_true_to_fit = median(kl_true_to_fit, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(distribution_summary, file.path(OUT_DIR, "paper_full_distribution_recovery_summary.csv"))
print(distribution_summary)

# -----------------------------------------------------------------------------
# Figure 1: main joint-recovery figure
# Posterior median estimates across simulated datasets; dashed line is truth.
# -----------------------------------------------------------------------------
plot_dat <- delay_results %>%
  dplyr::mutate(
    delay = factor(delay, levels = unname(DELAY_LABELS[c("lp", "ip", "gen", "si")])),
    moment = factor(moment, levels = c("Mean", "SD"))
  )

truth_plot <- TRUE_SUMMARY %>%
  dplyr::mutate(
    delay = factor(delay, levels = unname(DELAY_LABELS[c("lp", "ip", "gen", "si")])),
    moment = factor(moment, levels = c("Mean", "SD"))
  )

p_joint <- ggplot(plot_dat, aes(x = 1, y = median)) +
  geom_violin(width = 0.72, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  geom_point(
    data = truth_plot,
    aes(x = 1, y = true_value),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.2,
    size = 3
  ) +
  facet_grid(moment ~ delay, scales = "free_y") +
  scale_x_continuous(breaks = NULL) +
  labs(
    x = NULL,
    y = "Posterior median estimate, days",
    title = "Joint recovery of epidemiological delay summaries",
    subtitle = paste0(
      "n = ", N_PAIRS, " serial intervals per replicate; cross marks the known data-generating value"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  file.path(OUT_DIR, "figure_joint_delay_summary_recovery.png"),
  p_joint,
  width = 11,
  height = 6,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Figure 2: full distribution recovery
# Band = 5th-95th percentile of the fitted density across simulation refits.
# Line = exact DGM density.
# -----------------------------------------------------------------------------
density_summary <- fitted_densities %>%
  dplyr::group_by(delay_code, delay, time) %>%
  dplyr::summarise(
    fitted_q05 = quantile(density, 0.05, na.rm = TRUE),
    fitted_median = median(density, na.rm = TRUE),
    fitted_q95 = quantile(density, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(TRUE_DENSITIES, by = c("delay_code", "delay", "time")) %>%
  dplyr::mutate(delay = factor(delay, levels = unname(DELAY_LABELS[c("lp", "ip", "gen", "si")])))

p_density <- ggplot(density_summary, aes(x = time)) +
  geom_ribbon(aes(ymin = fitted_q05, ymax = fitted_q95), alpha = 0.20) +
  geom_line(aes(y = fitted_median), linewidth = 0.8) +
  geom_line(aes(y = true_density), linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ delay, scales = "free", ncol = 2) +
  labs(
    x = "Days",
    y = "Density",
    title = "Simultaneous recovery of the four delay distributions",
    subtitle = "Solid line: median fitted density across refits; band: 5th-95th percentile; dashed line: true DGM density"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(OUT_DIR, "figure_joint_delay_distribution_recovery.png"),
  p_density,
  width = 9,
  height = 7,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Figure 3: calibration / CRPS summary
# -----------------------------------------------------------------------------
p_crps <- ggplot(
  plot_dat,
  aes(x = delay, y = crps)
) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~ moment, scales = "free_y") +
  labs(
    x = NULL,
    y = "CRPS, days",
    title = "Probabilistic recovery of joint delay summaries"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(OUT_DIR, "figure_joint_delay_crps.png"),
  p_crps,
  width = 9,
  height = 4.8,
  dpi = 300
)

coverage_dat <- recovery_summary %>%
  dplyr::mutate(
    delay = factor(delay, levels = unname(DELAY_LABELS[c("lp", "ip", "gen", "si")])),
    moment = factor(moment, levels = c("Mean", "SD"))
  )

p_coverage <- ggplot(coverage_dat, aes(x = delay, y = coverage_90)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.90, linetype = "dashed") +
  facet_wrap(~ moment) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Empirical coverage",
    title = "Coverage of 90% posterior intervals"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(OUT_DIR, "figure_joint_delay_coverage.png"),
  p_coverage,
  width = 9,
  height = 4.8,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Diagnostics summary
# -----------------------------------------------------------------------------
diag_summary <- diagnostics %>%
  dplyr::summarise(
    n_successful_fits = dplyr::n(),
    failure_rate = 1 - dplyr::n() / N_REP,
    convergence_ok_rate = mean(
      max_rhat < 1.05 & n_divergent / (CHAINS * ITER_SAMPLING) < 0.05,
      na.rm = TRUE
    ),
    median_max_rhat = median(max_rhat, na.rm = TRUE),
    median_min_ess_bulk = median(min_ess_bulk, na.rm = TRUE),
    median_min_ess_tail = median(min_ess_tail, na.rm = TRUE),
    median_divergences = median(n_divergent, na.rm = TRUE)
  )

readr::write_csv(diag_summary, file.path(OUT_DIR, "simulation_diagnostic_summary.csv"))
print(diag_summary)

message("\nSimulation 4 complete. Results written to: ", OUT_DIR)

# Save them directly to the plotting directory
readr::write_csv(
  TRUE_DENSITIES,
  file.path(SIM_DIR, "true_joint_delay_densities.csv")
)
