# ============================================================
# Mpox SI model with flexible Yang-informed WH priors
#
# This script:
#   1) Fits Yang et al. Figure 2 mpox viral-load data using the manuscript V(t).
#   2) Converts posterior draws for theta = [log_r, log_d, m, log_kappa]
#      into infection-time priors for the Miura mpox SI model.
#   3) Refits the mpox observed-pair SI model using one of:
#        prior_mode = 1: independent Normal priors
#        prior_mode = 2: multivariate Normal Yang prior
#        prior_mode = 3: robust mixture prior, Yang MVN + broad MVN
#
# Why this improves the previous script:
#   - It preserves posterior covariance among log_r, log_d, m, and log_kappa.
#   - It allows prior-strength inflation.
#   - It carries uncertainty in the symptom-onset-to-infection-time shift for m.
#   - It allows a robust mixture prior so Yang data can inform the fit without
#     forcing compatibility if the SI data disagree.
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(readxl)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(loo)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

this_file_dir <- getwd()

YANG_STAN_FILE <- file.path(this_file_dir, "yang_mpox_wh_viral_load_model.stan")
if (!file.exists(YANG_STAN_FILE)) YANG_STAN_FILE <- "yang_mpox_wh_viral_load_model.stan"
if (!file.exists(YANG_STAN_FILE)) YANG_STAN_FILE <- file.path("/mnt/data", "yang_mpox_wh_viral_load_model.stan")
if (!file.exists(YANG_STAN_FILE)) stop("Could not find yang_mpox_wh_viral_load_model.stan")

MPOX_STAN_FILE <- file.path(this_file_dir, "mpox_mechanistic_obspair_flexible_WHpriors.stan")
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- "mpox_mechanistic_obspair_flexible_WHpriors.stan"
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- file.path("/mnt/data", "mpox_mechanistic_obspair_flexible_WHpriors.stan")
if (!file.exists(MPOX_STAN_FILE)) stop("Could not find mpox_mechanistic_obspair_flexible_WHpriors.stan")

YANG_SOURCE_XLSX <- Sys.getenv("YANG_SOURCE_XLSX", unset = "")
if (!nzchar(YANG_SOURCE_XLSX)) {
  candidate_paths <- c(
    "Yang_2024_mpox_source_data.xlsx",
    file.path(this_file_dir, "Yang_2024_mpox_source_data.xlsx"),
    file.path("/mnt/data", "Yang_2024_mpox_source_data.xlsx")
  )
  hit <- candidate_paths[file.exists(candidate_paths)]
  if (length(hit) == 0) {
    stop("Could not find Yang_2024_mpox_source_data.xlsx. Put it in the working directory or set YANG_SOURCE_XLSX.")
  }
  YANG_SOURCE_XLSX <- hit[[1]]
}

MPOX_DATA_URL <- "https://raw.githubusercontent.com/fmiura/MpxSI_2022/main/data/Anonym_All_data_score_indexed.csv"

OUT_DIR <- "mpox_score9_flexible_Yang_WH_prior_outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Controls
# ------------------------------------------------------------

# Yang viral-load fit controls
chains_yang <- 4
parallel_chains_yang <- 4
iter_warmup_yang <- 1000
iter_sampling_yang <- 1000
adapt_delta_yang <- 0.995
max_treedepth_yang <- 20
seed_yang <- 4401
refresh_yang <- 25

TREAT_Y_EQ_3_AS_LEFT_CENSORED <- TRUE
Y_LOD_LOG10 <- 3.0

# Choose viral-load source used to construct the prior.
# Options:
#   "shared_all_sites", "skin_lesion", "rectal_swab", "saliva", "oropharyngeal_swab"
YANG_PRIOR_SOURCE <- "shared_all_sites"

# Weak priors used only to fit Yang viral-load curves.
yang_fit_prior_list <- list(
  prior_log_r_mean = log(0.35),
  prior_log_r_sd = 1.0,
  prior_log_d_mean = log(0.20),
  prior_log_d_sd = 1.0,
  prior_m_mean = 6.0,          # Yang time scale: days post symptom onset
  prior_m_sd = 4.0,
  prior_log_kappa_mean = 0.0,
  prior_log_kappa_sd = 1.0
)

# Time-origin conversion: Yang m is post-symptom-onset; SI model m is post-infection.
INCUBATION_MEAN_FOR_SHIFT <- 8.0
INCUBATION_SD_FOR_SHIFT <- 3.0

# Prior strength. Values >1 make Yang prior more conservative.
YANG_PRIOR_SD_MULTIPLIER <- 2.0

# prior_mode sent to Stan:
#   1 = independent Normal priors, same information as old script
#   2 = Yang multivariate Normal prior
#   3 = robust mixture prior: Yang MVN + broad MVN
PRIOR_MODE <- 3

# In robust mixture mode, this is the prior weight on the Yang-informed MVN.
# Smaller means more robust / more willing to fall back to the broad prior.
YANG_MIX_WEIGHT <- 0.80

# Event-hazard and heterogeneity priors remain weak:
# Yang Figure 2 does not identify these SI-model event-hazard scales.
prior_log_lambda1_mean <- log(0.08)
prior_log_lambda1_sd <- 1.5
prior_log_lambda_I_mean <- log(0.08)
prior_log_lambda_I_sd <- 1.5
prior_log_sigma_z_mean <- log(0.50)
prior_log_sigma_z_sd <- 1.0

# Broad WH trajectory prior used for independent mode and mixture fallback.
broad_log_r_mean <- log(0.15)
broad_log_r_sd <- 1.0
broad_log_d_mean <- log(0.20)
broad_log_d_sd <- 1.0
broad_m_mean <- 8.0
broad_m_sd <- 5.0
broad_log_kappa_mean <- 0.0
broad_log_kappa_sd <- 1.25

# Mpox SI model controls
r_growth <- 0.0
GRID_BY_SI <- 0.25
GRID_BY_LU <- 0.25
GRID_BY_IP <- 0.05
SI_MIN <- -15
SI_MAX <- 40
LU_MAX <- 45
GT_MAX <- 2 * LU_MAX
IP_MIN <- 0.05
IP_MAX <- 60
N_Z_QUAD <- 11
Z_MAX_SD <- 3.5

chains_mpox <- 4
parallel_chains_mpox <- 4
iter_warmup_mpox <- 500
iter_sampling_mpox <- 500
adapt_delta_mpox <- 0.95
max_treedepth_mpox <- 12
seed_mpox <- 9001
refresh_mpox <- 5

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

extract_figure2_block <- function(raw, cols, site_label) {
  out <- raw[-c(1, 2), cols, drop = FALSE]
  names(out) <- c("hiv_status", "day_post_onset", "viral_load_log10")

  out %>%
    mutate(
      site = site_label,
      hiv_status = as.character(hiv_status),
      day_post_onset = suppressWarnings(as.numeric(day_post_onset)),
      viral_load_log10 = suppressWarnings(as.numeric(viral_load_log10))
    ) %>%
    filter(!is.na(day_post_onset), !is.na(viral_load_log10), !is.na(hiv_status)) %>%
    mutate(
      hiv = if_else(hiv_status == "HIV", 1L, 0L),
      is_left_censored = if_else(
        TREAT_Y_EQ_3_AS_LEFT_CENSORED & viral_load_log10 <= Y_LOD_LOG10 + 1e-9,
        1L,
        0L
      )
    )
}

make_yang_stan_data <- function(df) {
  df <- df %>% mutate(site = droplevels(factor(site)))
  c(
    list(
      N = nrow(df),
      t = as.numeric(df$day_post_onset),
      y_log10 = as.numeric(df$viral_load_log10),
      is_left_censored = as.integer(df$is_left_censored),
      y_lod_log10 = Y_LOD_LOG10,
      K_site = nlevels(df$site),
      site_id = as.integer(df$site),
      hiv = as.integer(df$hiv)
    ),
    yang_fit_prior_list
  )
}

make_spd <- function(S, jitter = 1e-6) {
  S <- as.matrix(S)
  S <- 0.5 * (S + t(S))
  eig <- eigen(S, symmetric = TRUE)
  min_eig <- min(eig$values)
  if (min_eig <= jitter) {
    S <- S + diag(jitter - min_eig + jitter, nrow(S))
  }
  S
}

construct_theta_prior <- function(yang_draws) {
  # Base theta on Yang time scale first.
  theta_yang <- as.matrix(yang_draws[, c("log_r", "log_d", "m", "log_kappa")])
  colnames(theta_yang) <- c("log_r", "log_d", "m", "log_kappa")

  # Convert m from post-symptom-onset to post-infection by shifting the mean.
  theta_infection <- theta_yang
  theta_infection[, "m"] <- theta_infection[, "m"] + INCUBATION_MEAN_FOR_SHIFT

  theta_mean <- colMeans(theta_infection)
  theta_cov <- cov(theta_infection)

  # Inflate Yang posterior uncertainty and add independent incubation-shift uncertainty to m.
  theta_cov <- (YANG_PRIOR_SD_MULTIPLIER^2) * theta_cov
  theta_cov["m", "m"] <- theta_cov["m", "m"] + INCUBATION_SD_FOR_SHIFT^2
  theta_cov <- make_spd(theta_cov)

  broad_mean <- c(
    log_r = broad_log_r_mean,
    log_d = broad_log_d_mean,
    m = broad_m_mean,
    log_kappa = broad_log_kappa_mean
  )
  broad_cov <- diag(c(
    broad_log_r_sd,
    broad_log_d_sd,
    broad_m_sd,
    broad_log_kappa_sd
  )^2)
  rownames(broad_cov) <- colnames(broad_cov) <- names(broad_mean)
  broad_cov <- make_spd(broad_cov)

  list(
    theta_mean = as.numeric(theta_mean),
    theta_cov = theta_cov,
    broad_mean = as.numeric(broad_mean),
    broad_cov = broad_cov,
    theta_corr_yang_post_onset = cor(theta_yang),
    theta_corr_infection_prior = cov2cor(theta_cov),
    marginal_prior = tibble(
      parameter = c("log_r", "log_d", "m", "log_kappa"),
      prior_mean = as.numeric(theta_mean),
      prior_sd = sqrt(diag(theta_cov)),
      broad_mean = as.numeric(broad_mean),
      broad_sd = sqrt(diag(broad_cov)),
      prior_mode = PRIOR_MODE,
      yang_prior_sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
      incubation_mean_shift = INCUBATION_MEAN_FOR_SHIFT,
      incubation_sd_shift = INCUBATION_SD_FOR_SHIFT,
      yang_mix_weight = YANG_MIX_WEIGHT
    )
  )
}

find_grid_index <- function(x, grid, by, tol = 1e-8) {
  idx <- round((x - grid[1]) / by) + 1L

  if (any(idx < 1L | idx > length(grid))) {
    bad <- x[idx < 1L | idx > length(grid)]
    stop("Value outside grid: ", paste(head(bad, 10), collapse = ", "))
  }

  if (any(abs(grid[idx] - x) > tol)) {
    bad <- x[abs(grid[idx] - x) > tol]
    stop("Value not exactly on grid: ", paste(head(bad, 10), collapse = ", "))
  }

  as.integer(idx)
}

expand_weighted_obs_index <- function(obs_index, counts) {
  rep.int(as.integer(obs_index), as.integer(counts))
}

make_l_t_idx <- function(l_grid, t_grid, by) find_grid_index(l_grid, t_grid, by)

make_lu_t_idx <- function(l_grid, u_grid, t_grid, by) {
  out <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
  for (jl in seq_along(l_grid)) {
    for (ju in seq_along(u_grid)) {
      out[jl, ju] <- find_grid_index(l_grid[jl] + u_grid[ju], t_grid, by)
    }
  }
  out
}

make_gt_idx <- function(l_grid, u_grid, gt_grid, by) {
  out <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
  for (jl in seq_along(l_grid)) {
    for (ju in seq_along(u_grid)) {
      out[jl, ju] <- find_grid_index(l_grid[jl] + u_grid[ju], gt_grid, by)
    }
  }
  out
}

make_si_g_delta_idx <- function(si_grid, gt_grid, delta_grid, by) {
  out <- matrix(NA_integer_, nrow = length(si_grid), ncol = length(gt_grid))
  for (ks in seq_along(si_grid)) {
    for (kg in seq_along(gt_grid)) {
      delta <- si_grid[ks] - gt_grid[kg]
      out[ks, kg] <- find_grid_index(delta, delta_grid, by)
    }
  }
  out
}

make_ip_sparse <- function(delta_grid, ip_grid, by_ip, tol = 1e-8) {
  sparse_delta <- integer(0)
  sparse_i1 <- integer(0)
  sparse_i2 <- integer(0)

  for (dd in seq_along(delta_grid)) {
    target <- delta_grid[dd] + ip_grid
    idx <- round((target - ip_grid[1]) / by_ip) + 1L
    ok <- idx >= 1L & idx <= length(ip_grid)

    if (any(ok)) {
      ok_idx <- which(ok)
      ok[ok_idx] <- abs(ip_grid[idx[ok_idx]] - target[ok_idx]) <= tol
    }

    if (any(ok)) {
      sparse_delta <- c(sparse_delta, rep.int(dd, sum(ok)))
      sparse_i1 <- c(sparse_i1, which(ok))
      sparse_i2 <- c(sparse_i2, idx[ok])
    }
  }

  list(
    N_ip_sparse = length(sparse_delta),
    ip_sparse_delta = as.integer(sparse_delta),
    ip_sparse_i1 = as.integer(sparse_i1),
    ip_sparse_i2 = as.integer(sparse_i2)
  )
}

make_z_quad <- function(n_z = 11, z_max_sd = 3.5) {
  x <- seq(-z_max_sd, z_max_sd, length.out = n_z)
  dx <- if (length(x) > 1L) diff(x)[1L] else 1
  w_raw <- dnorm(x) * dx
  w <- w_raw / sum(w_raw)
  list(x = x, w = w)
}

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

# ------------------------------------------------------------
# 1. Extract and fit Yang Figure 2 viral-load data
# ------------------------------------------------------------

raw_fig2 <- read_xlsx(YANG_SOURCE_XLSX, sheet = "Figure 2", col_names = FALSE)

yang_vl_dat <- bind_rows(
  extract_figure2_block(raw_fig2, cols = 2:4,   site_label = "skin_lesion"),
  extract_figure2_block(raw_fig2, cols = 6:8,   site_label = "rectal_swab"),
  extract_figure2_block(raw_fig2, cols = 10:12, site_label = "saliva"),
  extract_figure2_block(raw_fig2, cols = 14:16, site_label = "oropharyngeal_swab")
) %>%
  mutate(site = factor(site, levels = c("skin_lesion", "rectal_swab", "saliva", "oropharyngeal_swab")))

write_csv(yang_vl_dat, file.path(OUT_DIR, "yang_figure2_viral_load_long_used_for_priors.csv"))

if (YANG_PRIOR_SOURCE != "shared_all_sites") {
  yang_fit_dat <- yang_vl_dat %>% filter(site == YANG_PRIOR_SOURCE)
} else {
  yang_fit_dat <- yang_vl_dat
}

if (nrow(yang_fit_dat) < 10) stop("Yang prior source has too few observations: ", YANG_PRIOR_SOURCE)

message("\nYang viral-load source: ", YANG_PRIOR_SOURCE)
print(yang_fit_dat %>% count(site, hiv_status, is_left_censored))
message("Yang time variable is day_post_onset: days post symptom onset.")

mod_yang <- cmdstan_model(YANG_STAN_FILE, force_recompile = FALSE)
yang_fit <- mod_yang$sample(
  data = make_yang_stan_data(yang_fit_dat),
  chains = chains_yang,
  parallel_chains = parallel_chains_yang,
  iter_warmup = iter_warmup_yang,
  iter_sampling = iter_sampling_yang,
  seed = seed_yang,
  adapt_delta = adapt_delta_yang,
  max_treedepth = max_treedepth_yang,
  refresh = refresh_yang,
  save_warmup = FALSE
)

yang_fit$save_object(file.path(OUT_DIR, paste0("fit_yang_wh_", YANG_PRIOR_SOURCE, ".rds")))

yang_summary <- yang_fit$summary(
  variables = c("r", "d", "m", "kappa", "log_r", "log_d", "log_kappa", "alpha_site", "beta_hiv", "sigma_y")
) %>%
  mutate(fit = YANG_PRIOR_SOURCE) %>%
  select(fit, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

write_csv(yang_summary, file.path(OUT_DIR, paste0("summary_yang_wh_", YANG_PRIOR_SOURCE, ".csv")))
print(yang_summary)

yang_draws <- as_draws_df(yang_fit$draws(variables = c("log_r", "log_d", "m", "log_kappa", "r", "d", "kappa")))
write_csv(as.data.frame(yang_draws), file.path(OUT_DIR, paste0("draws_yang_wh_", YANG_PRIOR_SOURCE, ".csv")))

theta_prior <- construct_theta_prior(yang_draws)
write_csv(theta_prior$marginal_prior, file.path(OUT_DIR, "constructed_theta_prior_marginals.csv"))
write.csv(theta_prior$theta_cov, file.path(OUT_DIR, "constructed_theta_prior_cov.csv"), row.names = TRUE)
write.csv(theta_prior$broad_cov, file.path(OUT_DIR, "constructed_theta_broad_cov.csv"), row.names = TRUE)
write.csv(theta_prior$theta_corr_yang_post_onset, file.path(OUT_DIR, "yang_theta_posterior_correlation_post_onset.csv"), row.names = TRUE)
write.csv(theta_prior$theta_corr_infection_prior, file.path(OUT_DIR, "constructed_theta_prior_correlation_infection_time.csv"), row.names = TRUE)

message("\nConstructed theta prior marginals:")
print(theta_prior$marginal_prior)
message("\nYang posterior correlation on post-onset scale:")
print(round(theta_prior$theta_corr_yang_post_onset, 3))
message("\nConstructed prior correlation on infection-time scale:")
print(round(theta_prior$theta_corr_infection_prior, 3))

# ------------------------------------------------------------
# 2. Build mpox SI data and lookup tables
# ------------------------------------------------------------

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

mpox_raw <- read_csv(MPOX_DATA_URL, show_col_types = FALSE)
mpox_score9 <- mpox_raw %>%
  filter(score == 9, Conf_pair == "high", Conf_sym == "high") %>%
  mutate(serial = as.numeric(SI))

if (nrow(mpox_score9) != 34) {
  warning("Expected 34 score==9 high/high pairs, but found ", nrow(mpox_score9), ". Check upstream data.")
}

N_PAIRS <- nrow(mpox_score9)

serial_w_mpox <- mpox_score9 %>%
  count(serial, name = "n") %>%
  arrange(serial)

write_csv(mpox_score9, file.path(OUT_DIR, "mpox_score9_high_confidence_pairs.csv"))
write_csv(serial_w_mpox, file.path(OUT_DIR, "mpox_score9_weighted_serial_intervals.csv"))

message("\nMpox high-confidence serial data:")
message("  unique SI values = ", nrow(serial_w_mpox), ", total pairs = ", sum(serial_w_mpox$n))
message("  SI range = ", min(serial_w_mpox$serial), " to ", max(serial_w_mpox$serial))

message("\nBuilding lookup tables...")
l_t_idx <- make_l_t_idx(l_grid, t_grid, GRID_BY_LU)
lu_t_idx <- make_lu_t_idx(l_grid, u_grid, t_grid, GRID_BY_LU)
gt_idx <- make_gt_idx(l_grid, u_grid, gt_grid, GRID_BY_LU)
si_g_delta_idx <- make_si_g_delta_idx(si_grid, gt_grid, delta_grid, GRID_BY_SI)
ip_sparse <- make_ip_sparse(delta_grid, ip_grid, GRID_BY_IP)
zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

make_stan_data_mech <- function(serial_df, r_growth) {
  obs_index <- find_grid_index(
    x = as.numeric(serial_df$serial),
    grid = si_grid,
    by = GRID_BY_SI
  )

  obs_index_pt <- expand_weighted_obs_index(obs_index, serial_df$n)

  # For prior_mode == 1, these independent priors are used.
  # For prior_mode == 2/3, theta MVN/mixture priors are used instead for WH trajectory.
  independent_prior <- list(
    prior_log_r_mean = theta_prior$marginal_prior$prior_mean[theta_prior$marginal_prior$parameter == "log_r"],
    prior_log_r_sd = theta_prior$marginal_prior$prior_sd[theta_prior$marginal_prior$parameter == "log_r"],
    prior_log_d_mean = theta_prior$marginal_prior$prior_mean[theta_prior$marginal_prior$parameter == "log_d"],
    prior_log_d_sd = theta_prior$marginal_prior$prior_sd[theta_prior$marginal_prior$parameter == "log_d"],
    prior_m_mean = theta_prior$marginal_prior$prior_mean[theta_prior$marginal_prior$parameter == "m"],
    prior_m_sd = theta_prior$marginal_prior$prior_sd[theta_prior$marginal_prior$parameter == "m"],
    prior_log_kappa_mean = theta_prior$marginal_prior$prior_mean[theta_prior$marginal_prior$parameter == "log_kappa"],
    prior_log_kappa_sd = theta_prior$marginal_prior$prior_sd[theta_prior$marginal_prior$parameter == "log_kappa"]
  )

  c(
    list(
      N_si = nrow(serial_df),
      obs_index = as.integer(obs_index),
      counts = as.integer(serial_df$n),

      N_obs = length(obs_index_pt),
      obs_index_pt = as.integer(obs_index_pt),

      N_norm = length(si_grid),
      si_grid = as.numeric(si_grid),
      dx_si = GRID_BY_SI,

      N_l = length(l_grid),
      l_grid = as.numeric(l_grid),
      dx_l = GRID_BY_LU,

      N_u = length(u_grid),
      u_grid = as.numeric(u_grid),
      dx_u = GRID_BY_LU,

      N_gt = length(gt_grid),
      gt_grid = as.numeric(gt_grid),

      N_ip = length(ip_grid),
      ip_grid = as.numeric(ip_grid),
      dx_ip = GRID_BY_IP,

      N_t = length(t_grid),
      t_grid = as.numeric(t_grid),
      dx_t = GRID_BY_LU,

      N_delta = length(delta_grid),
      delta_grid = as.numeric(delta_grid),

      l_t_idx = as.integer(l_t_idx),
      lu_t_idx = lu_t_idx,
      gt_idx = gt_idx,
      si_g_delta_idx = si_g_delta_idx,

      N_ip_sparse = ip_sparse$N_ip_sparse,
      ip_sparse_delta = ip_sparse$ip_sparse_delta,
      ip_sparse_i1 = ip_sparse$ip_sparse_i1,
      ip_sparse_i2 = ip_sparse$ip_sparse_i2,

      park_weight_ip = as.numeric(exp(-r_growth * ip_grid)),

      N_z = N_Z_QUAD,
      z_std_grid = as.numeric(zq$x),
      z_std_w = as.numeric(zq$w),

      prior_mode = PRIOR_MODE,
      prior_theta_mean = theta_prior$theta_mean,
      prior_theta_cov = theta_prior$theta_cov,
      prior_theta_broad_mean = theta_prior$broad_mean,
      prior_theta_broad_cov = theta_prior$broad_cov,
      prior_theta_mix_weight = YANG_MIX_WEIGHT
    ),
    list(
      prior_log_lambda1_mean = prior_log_lambda1_mean,
      prior_log_lambda1_sd = prior_log_lambda1_sd,
      prior_log_lambda_I_mean = prior_log_lambda_I_mean,
      prior_log_lambda_I_sd = prior_log_lambda_I_sd,
      prior_log_sigma_z_mean = prior_log_sigma_z_mean,
      prior_log_sigma_z_sd = prior_log_sigma_z_sd
    ),
    independent_prior
  )
}

data_mpox <- make_stan_data_mech(serial_w_mpox, r_growth)

# ------------------------------------------------------------
# 3. Fit mpox SI model
# ------------------------------------------------------------

mod_mpox <- cmdstan_model(MPOX_STAN_FILE, force_recompile = TRUE)

make_inits_mech <- function() {
  function() {
    list(
      log_lambda1 = prior_log_lambda1_mean,
      log_lambda_I = prior_log_lambda_I_mean,
      log_r = data_mpox$prior_theta_mean[1],
      log_d = data_mpox$prior_theta_mean[2],
      m = min(max(data_mpox$prior_theta_mean[3], 0.5), 30),
      log_kappa = data_mpox$prior_theta_mean[4],
      log_sigma_z = prior_log_sigma_z_mean
    )
  }
}

message("\nFitting mpox score==9 SI model with flexible WH prior...")
message("  PRIOR_MODE = ", PRIOR_MODE)
message("  YANG_PRIOR_SOURCE = ", YANG_PRIOR_SOURCE)
message("  YANG_PRIOR_SD_MULTIPLIER = ", YANG_PRIOR_SD_MULTIPLIER)
message("  YANG_MIX_WEIGHT = ", YANG_MIX_WEIGHT)

fit_mpox <- mod_mpox$sample(
  data = data_mpox,
  chains = chains_mpox,
  parallel_chains = parallel_chains_mpox,
  iter_warmup = 150,
  iter_sampling = 150,
  seed = seed_mpox,
  init = make_inits_mech(),
  adapt_delta = adapt_delta_mpox,
  max_treedepth = max_treedepth_mpox,
  refresh = refresh_mpox,
  save_warmup = FALSE
)

fit_slug <- paste0(
  "source_", YANG_PRIOR_SOURCE,
  "_mode_", PRIOR_MODE,
  "_sdmult_", gsub("\\.", "p", YANG_PRIOR_SD_MULTIPLIER),
  "_mix_", gsub("\\.", "p", YANG_MIX_WEIGHT)
)

fit_mpox$save_object(file.path(OUT_DIR, paste0("fit_mpox_score9_", fit_slug, ".rds")))

pars_summary <- c(
  "log_lambda1", "log_lambda_I", "log_r", "log_d", "m", "log_kappa",
  "log_sigma_z", "mean_gen", "sd_gen"
)

summ_mpox <- fit_mpox$summary(variables = pars_summary) %>%
  mutate(
    dataset = "mpox_score9",
    fit = fit_slug,
    prior_mode = PRIOR_MODE,
    yang_prior_source = YANG_PRIOR_SOURCE,
    yang_prior_sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
    yang_mix_weight = YANG_MIX_WEIGHT
  ) %>%
  select(dataset, fit, prior_mode, yang_prior_source, yang_prior_sd_multiplier,
         yang_mix_weight, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_mpox)
write_csv(summ_mpox, file.path(OUT_DIR, paste0("summary_mpox_score9_", fit_slug, ".csv")))

gen_vars_mpox <- grep("^generation_prob_grid\\[", posterior::variables(fit_mpox$draws()), value = TRUE)
if (length(gen_vars_mpox) == 0) stop("generation_prob_grid was not found in the fitted object.")
message("Confirmed generation_prob_grid exists: entries = ", length(gen_vars_mpox))

log_lik_mat <- extract_log_lik_matrix(fit_mpox, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit_mpox, "log_lik_grouped")

if (ncol(log_lik_mat) != N_PAIRS) {
  stop("Pointwise log_lik has ", ncol(log_lik_mat), " columns, expected ", N_PAIRS, ". LOO/WAIC not comparable.")
}

saveRDS(log_lik_mat, file.path(OUT_DIR, paste0("log_lik_mpox_score9_", fit_slug, ".rds")))
saveRDS(log_lik_grouped_mat, file.path(OUT_DIR, paste0("log_lik_grouped_mpox_score9_", fit_slug, ".rds")))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)

saveRDS(loo_fit, file.path(OUT_DIR, paste0("loo_mpox_score9_", fit_slug, ".rds")))
saveRDS(waic_fit, file.path(OUT_DIR, paste0("waic_mpox_score9_", fit_slug, ".rds")))

loo_waic_mpox <- tibble(
  dataset = "mpox_score9",
  fit = fit_slug,
  prior_mode = PRIOR_MODE,
  yang_prior_source = YANG_PRIOR_SOURCE,
  yang_prior_sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
  incubation_mean_shift = INCUBATION_MEAN_FOR_SHIFT,
  incubation_sd_shift = INCUBATION_SD_FOR_SHIFT,
  yang_mix_weight = YANG_MIX_WEIGHT,
  n_pointwise_obs = ncol(log_lik_mat),
  elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
  se_elpd_loo = loo_fit$estimates["elpd_loo", "SE"],
  looic = loo_fit$estimates["looic", "Estimate"],
  se_looic = loo_fit$estimates["looic", "SE"],
  p_loo = loo_fit$estimates["p_loo", "Estimate"],
  elpd_waic = waic_fit$estimates["elpd_waic", "Estimate"],
  se_elpd_waic = waic_fit$estimates["elpd_waic", "SE"],
  waic = waic_fit$estimates["waic", "Estimate"],
  se_waic = waic_fit$estimates["waic", "SE"],
  p_waic = waic_fit$estimates["p_waic", "Estimate"]
)

print(loo_waic_mpox)
write_csv(loo_waic_mpox, file.path(OUT_DIR, paste0("metrics_mpox_score9_", fit_slug, ".csv")))

message("\nDiagnostics:")
diag_txt <- capture.output(fit_mpox$cmdstan_diagnose())
writeLines(diag_txt, file.path(OUT_DIR, paste0("cmdstan_diagnose_mpox_score9_", fit_slug, ".txt")))
cat(paste(diag_txt, collapse = "\n"))

draws_mpox <- as_draws_df(fit_mpox$draws(variables = pars_summary)) %>%
  mutate(
    dataset = "mpox_score9",
    fit = fit_slug,
    prior_mode = PRIOR_MODE,
    yang_prior_source = YANG_PRIOR_SOURCE,
    yang_prior_sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
    yang_mix_weight = YANG_MIX_WEIGHT
  )

write_csv(as.data.frame(draws_mpox), file.path(OUT_DIR, paste0("draws_mpox_score9_", fit_slug, ".csv")))

p_mean_gen <- ggplot(draws_mpox, aes(x = mean_gen)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = paste0("Mpox score==9 flexible WH prior: ", fit_slug)
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, paste0("mean_gen_density_mpox_score9_", fit_slug, ".png")),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)

message("\nDone. Outputs saved in: ", OUT_DIR)
