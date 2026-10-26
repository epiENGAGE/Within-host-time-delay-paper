# ============================================================
# Full mpox WH-prior analysis
# 1) Fit Yang et al. mpox viral-load data using your V(t) shape.
# 2) Convert the Yang posterior into priors for r, d, kappa, and m.
# 3) Refit the Miura mpox score==9 transmission-pair SI model with those priors.
#
# IMPORTANT TIME-SCALE FACT:
#   Yang Figure 2 x-axis is day_post_onset = days post symptom onset.
#   The transmission-pair model's V(t) uses t = days since infection.
#   Therefore the Yang posterior m is a peak time AFTER SYMPTOM ONSET.
#   The SI model prior for m must be shifted to infection time by adding an
#   incubation-period offset. Set INCUBATION_MEAN_FOR_SHIFT and
#   INCUBATION_SD_FOR_SHIFT below.
#
# This script saves to a NEW output directory and will not overwrite the weak-prior fit.
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

MPOX_STAN_FILE <- file.path(this_file_dir, "mpox_mechanistic_prop_updatedV_fast_genprob_obspair_mmix.stan")
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- "mpox_mechanistic_prop_updatedV_fast_genprob_obspair_mmix.stan"
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- file.path("mpox_mechanistic_prop_updatedV_fast_genprob_obspair_mmix.stan")
if (!file.exists(MPOX_STAN_FILE)) stop("Could not find mpox_mechanistic_prop_updatedV_fast_genprob_obspair.stan")

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

OUT_DIR_YANG <- "yang_mpox_wh_parameter_outputs_for_priors"
OUT_DIR_MPOX <- "stan_outputs_mpox_score9_obspair_yangWHpriors"
dir.create(OUT_DIR_YANG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_MPOX, showWarnings = FALSE, recursive = TRUE)

message("Using Yang Stan file: ", YANG_STAN_FILE)
message("Using mpox SI Stan file: ", MPOX_STAN_FILE)
message("Using Yang workbook: ", YANG_SOURCE_XLSX)
message("Yang output directory: ", OUT_DIR_YANG)
message("Mpox SI output directory: ", OUT_DIR_MPOX)

# ------------------------------------------------------------
# Controls
# ------------------------------------------------------------

# Yang viral-load fitting controls
chains_yang <- 4
parallel_chains_yang <- 4
iter_warmup_yang <- 1000
iter_sampling_yang <- 1000
adapt_delta_yang <- 0.95
max_treedepth_yang <- 12
seed_yang <- 4401
refresh_yang <- 25

TREAT_Y_EQ_3_AS_LEFT_CENSORED <- TRUE
Y_LOD_LOG10 <- 3.0

# Use all four Figure 2 body sites with one shared r/d/m/kappa shape and site-specific intercepts.
# This produces one prior source for r, d, m_post_onset, and kappa.
# If you want a lesion-only prior, filter yang_vl_dat below before fitting.
YANG_PRIOR_LABEL <- "shared_all_sites"

# Priors used only to fit Yang viral-load curves. These are weak regularising priors.
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

# Time-origin conversion for m.
# Yang m is measured in days post symptom onset; the SI model m is days since infection.
# Do NOT use a hard +8 shift. Instead, build the prior from posterior draws:
#   m_infection_draw = m_yang_post_onset_draw + incubation_draw
# and approximate that empirical convolution by a finite normal mixture in Stan.
#
# Use either:
#   INCUBATION_PRIOR_MODE = "lognormal_quantiles"
#   INCUBATION_PRIOR_MODE = "csv_draws" with a CSV containing a column named incubation_days
INCUBATION_PRIOR_MODE <- "lognormal_quantiles"
INCUBATION_DRAWS_CSV <- ""

# Edit these to the incubation distribution you want. Quantiles are better than a Normal mean/sd
# because incubation times are positive and usually right-skewed.
INCUBATION_MEDIAN_DAYS <- 8.5
INCUBATION_Q05_DAYS <- 4.0
INCUBATION_Q95_DAYS <- 21.0

# Number of mixture components used to approximate the shifted m prior.
# 5 is usually enough for a smooth 1D prior without making sampling annoying.
M_PRIOR_N_MIX <- 5
M_PRIOR_MIN_COMPONENT_SD <- 0.20

# Strength of Yang WH priors in the SI model.
# 1.0 = use Yang posterior uncertainty directly; >1 inflates Yang uncertainty before building priors.
YANG_PRIOR_SD_MULTIPLIER <- 1.0

# Leave lambda1/lambdaI/sigma_z weak: Yang Figure 2 does not identify these SI-model event-hazard scales.
prior_log_lambda1_mean <- log(0.08)
prior_log_lambda1_sd <- 1.5
prior_log_lambda_I_mean <- log(0.08)
prior_log_lambda_I_sd <- 1.5
prior_log_sigma_z_mean <- log(0.50)
prior_log_sigma_z_sd <- 1.0

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
# 1. Extract and fit Yang Figure 2 viral-load data
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

raw_fig2 <- read_xlsx(YANG_SOURCE_XLSX, sheet = "Figure 2", col_names = FALSE)

# Exact mapping used here from the workbook:
#   cols 2:4   = skin lesion block
#   cols 6:8   = rectal swab block
#   cols 10:12 = saliva block
#   cols 14:16 = oropharyngeal swab block
# Within each block: HIV status, day post symptom onset, viral load log10.
yang_vl_dat <- bind_rows(
  extract_figure2_block(raw_fig2, cols = 2:4,   site_label = "skin_lesion"),
  extract_figure2_block(raw_fig2, cols = 6:8,   site_label = "rectal_swab"),
  extract_figure2_block(raw_fig2, cols = 10:12, site_label = "saliva"),
  extract_figure2_block(raw_fig2, cols = 14:16, site_label = "oropharyngeal_swab")
) %>%
  mutate(site = factor(site, levels = c("skin_lesion", "rectal_swab", "saliva", "oropharyngeal_swab")))

write_csv(yang_vl_dat, file.path(OUT_DIR_YANG, "yang_figure2_viral_load_long_used_for_priors.csv"))

message("\nYang Figure 2 viral-load data used for WH prior fit:")
print(yang_vl_dat %>% count(site, hiv_status, is_left_censored))
message("Yang time variable is day_post_onset: days post symptom onset.")

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

mod_yang <- cmdstan_model(YANG_STAN_FILE, force_recompile = TRUE)
yang_fit <- mod_yang$sample(
  data = make_yang_stan_data(yang_vl_dat),
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

yang_fit$save_object(file.path(OUT_DIR_YANG, paste0("fit_yang_wh_", YANG_PRIOR_LABEL, ".rds")))

yang_summary <- yang_fit$summary(variables = c("r", "d", "m", "kappa", "log_r", "log_d", "log_kappa", "alpha_site", "beta_hiv", "sigma_y")) %>%
  mutate(fit = YANG_PRIOR_LABEL) %>%
  select(fit, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)
write_csv(yang_summary, file.path(OUT_DIR_YANG, paste0("summary_yang_wh_", YANG_PRIOR_LABEL, ".csv")))
print(yang_summary)

# ------------------------------------------------------------
# 2. Convert Yang WH posterior to SI-model priors
# ------------------------------------------------------------

yang_draws <- as_draws_df(yang_fit$draws(variables = c("log_r", "log_d", "m", "log_kappa", "r", "d", "kappa")))
write_csv(as.data.frame(yang_draws), file.path(OUT_DIR_YANG, paste0("draws_yang_wh_", YANG_PRIOR_LABEL, ".csv")))

posterior_mean_sd <- function(x) c(mean = mean(x), sd = sd(x))

# Convert a desired lognormal median and 5th/95th percentile into lognormal parameters.
# This avoids pretending the incubation shift is Normal.
lnorm_params_from_median_q05_q95 <- function(median, q05, q95) {
  if (median <= 0 || q05 <= 0 || q95 <= 0 || q05 >= q95) {
    stop("Incubation quantiles must be positive with q05 < q95.")
  }
  meanlog <- log(median)
  sdlog <- (log(q95) - log(q05)) / (qnorm(0.95) - qnorm(0.05))
  list(meanlog = meanlog, sdlog = sdlog)
}

make_incubation_draws <- function(n) {
  if (INCUBATION_PRIOR_MODE == "csv_draws") {
    if (!nzchar(INCUBATION_DRAWS_CSV) || !file.exists(INCUBATION_DRAWS_CSV)) {
      stop("INCUBATION_PRIOR_MODE='csv_draws' but INCUBATION_DRAWS_CSV does not exist.")
    }
    inc_df <- read_csv(INCUBATION_DRAWS_CSV, show_col_types = FALSE)
    if (!("incubation_days" %in% names(inc_df))) {
      stop("Incubation CSV must contain a column named incubation_days.")
    }
    inc_pool <- inc_df$incubation_days
    inc_pool <- inc_pool[is.finite(inc_pool) & inc_pool > 0]
    if (length(inc_pool) < 100) stop("Need at least 100 positive incubation draws.")
    return(sample(inc_pool, size = n, replace = TRUE))
  }

  if (INCUBATION_PRIOR_MODE == "lognormal_quantiles") {
    inc_par <- lnorm_params_from_median_q05_q95(
      median = INCUBATION_MEDIAN_DAYS,
      q05 = INCUBATION_Q05_DAYS,
      q95 = INCUBATION_Q95_DAYS
    )
    return(rlnorm(n, meanlog = inc_par$meanlog, sdlog = inc_par$sdlog))
  }

  stop("Unknown INCUBATION_PRIOR_MODE: ", INCUBATION_PRIOR_MODE)
}

# Equal-frequency normal mixture approximation to an empirical 1D prior.
# This preserves skewness and tail shape much better than collapsing to one mean/sd.
make_equal_frequency_normal_mixture <- function(x, K = 5, min_sd = 0.20) {
  x <- sort(x[is.finite(x)])
  if (length(x) < K * 20) stop("Too few draws to build mixture prior.")
  groups <- cut(seq_along(x), breaks = K, labels = FALSE)
  out <- tibble(
    component = seq_len(K),
    weight = as.numeric(tabulate(groups, nbins = K)) / length(x),
    mean = as.numeric(tapply(x, groups, mean)),
    sd = pmax(as.numeric(tapply(x, groups, sd)), min_sd)
  )
  out
}

# Inflate Yang posterior draws if requested. For m this happens before the incubation shift.
inflate_draws <- function(x, multiplier) {
  if (multiplier <= 1) return(x)
  x + rnorm(length(x), mean = 0, sd = sqrt(multiplier^2 - 1) * sd(x))
}

set.seed(seed_yang + 99)

yang_log_r_prior_draws <- inflate_draws(yang_draws$log_r, YANG_PRIOR_SD_MULTIPLIER)
yang_log_d_prior_draws <- inflate_draws(yang_draws$log_d, YANG_PRIOR_SD_MULTIPLIER)
yang_log_kappa_prior_draws <- inflate_draws(yang_draws$log_kappa, YANG_PRIOR_SD_MULTIPLIER)
yang_m_post_onset_prior_draws <- inflate_draws(yang_draws$m, YANG_PRIOR_SD_MULTIPLIER)

incubation_prior_draws <- make_incubation_draws(length(yang_m_post_onset_prior_draws))
m_infection_prior_draws <- yang_m_post_onset_prior_draws + incubation_prior_draws

# Keep impossible/extreme m prior draws from driving mixture tails outside the Stan parameter range.
# The Stan parameter itself is constrained to [0.25, 35].
m_infection_prior_draws <- m_infection_prior_draws[
  is.finite(m_infection_prior_draws) &
    m_infection_prior_draws > 0.25 &
    m_infection_prior_draws < 35
]

if (length(m_infection_prior_draws) < 500) {
  stop("Too few valid shifted m prior draws after filtering. Check incubation prior.")
}

m_prior_mix <- make_equal_frequency_normal_mixture(
  m_infection_prior_draws,
  K = M_PRIOR_N_MIX,
  min_sd = M_PRIOR_MIN_COMPONENT_SD
)

prior_log_r_mean <- mean(yang_log_r_prior_draws)
prior_log_r_sd <- max(sd(yang_log_r_prior_draws), 0.05)

prior_log_d_mean <- mean(yang_log_d_prior_draws)
prior_log_d_sd <- max(sd(yang_log_d_prior_draws), 0.05)

prior_log_kappa_mean <- mean(yang_log_kappa_prior_draws)
prior_log_kappa_sd <- max(sd(yang_log_kappa_prior_draws), 0.05)

# These are only used for initialization and summaries; Stan uses m_prior_mix for the actual prior.
prior_m_mean <- mean(m_infection_prior_draws)
prior_m_sd <- max(sd(m_infection_prior_draws), 0.25)

write_csv(
  tibble(m_infection_prior_draw = m_infection_prior_draws),
  file.path(OUT_DIR_MPOX, "yangWH_m_infection_prior_draws_empirical_convolution.csv")
)
write_csv(
  m_prior_mix,
  file.path(OUT_DIR_MPOX, "yangWH_m_infection_prior_normal_mixture.csv")
)

constructed_priors <- tibble(
  parameter = c("log_r", "log_d", "log_kappa", "m"),
  transmission_model_time_scale = c("days^-1 on infection-time V(t)", "days^-1 on infection-time V(t)", "transition sharpness on infection-time V(t)", "days since infection"),
  yang_source_time_scale = c("same rate parameter from day-post-onset V(t)", "same rate parameter from day-post-onset V(t)", "same shape parameter from day-post-onset V(t)", "Yang m is days post symptom onset"),
  prior_construction = c(
    "Yang posterior draws, optionally inflated",
    "Yang posterior draws, optionally inflated",
    "Yang posterior draws, optionally inflated",
    "Empirical convolution: Yang m posterior draws + positive incubation-period draws, approximated in Stan by a finite normal mixture"
  ),
  prior_mean = c(prior_log_r_mean, prior_log_d_mean, prior_log_kappa_mean, prior_m_mean),
  prior_sd = c(prior_log_r_sd, prior_log_d_sd, prior_log_kappa_sd, prior_m_sd),
  yang_posterior_mean = c(mean(yang_draws$log_r), mean(yang_draws$log_d), mean(yang_draws$log_kappa), mean(yang_draws$m)),
  yang_posterior_sd = c(sd(yang_draws$log_r), sd(yang_draws$log_d), sd(yang_draws$log_kappa), sd(yang_draws$m)),
  sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
  incubation_prior_mode = c(NA, NA, NA, INCUBATION_PRIOR_MODE),
  incubation_median = c(NA, NA, NA, INCUBATION_MEDIAN_DAYS),
  incubation_q05 = c(NA, NA, NA, INCUBATION_Q05_DAYS),
  incubation_q95 = c(NA, NA, NA, INCUBATION_Q95_DAYS),
  m_prior_mixture_components = c(NA, NA, NA, M_PRIOR_N_MIX)
)

write_csv(constructed_priors, file.path(OUT_DIR_MPOX, "yangWH_constructed_priors_used_in_mpox_SI_fit.csv"))
message("\nConstructed Yang-WH priors for mpox SI model:")
print(constructed_priors)
message("\nNormal mixture approximation used for m prior:")
print(m_prior_mix)

# ------------------------------------------------------------
# 3. Refit mpox score==9 SI model using Yang-informed priors
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
  warning("Expected 34 score==9 high/high pairs, but found ", nrow(mpox_score9), ". Check the upstream data file.")
}

serial_w_mpox <- mpox_score9 %>%
  count(serial, name = "n") %>%
  arrange(serial)

write.csv(mpox_score9, file.path(OUT_DIR_MPOX, "mpox_score9_high_confidence_pairs.csv"), row.names = FALSE)
write.csv(serial_w_mpox, file.path(OUT_DIR_MPOX, "mpox_score9_weighted_serial_intervals.csv"), row.names = FALSE)

message("\nMpox high-confidence serial data:")
message("  unique SI values = ", nrow(serial_w_mpox), ", total pairs = ", sum(serial_w_mpox$n))
message("  SI range = ", min(serial_w_mpox$serial), " to ", max(serial_w_mpox$serial))

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

message("\nBuilding lookup tables for mpox SI fit...")
l_t_idx <- make_l_t_idx(l_grid, t_grid, GRID_BY_LU)
lu_t_idx <- make_lu_t_idx(l_grid, u_grid, t_grid, GRID_BY_LU)
gt_idx <- make_gt_idx(l_grid, u_grid, gt_grid, GRID_BY_LU)
si_g_delta_idx <- make_si_g_delta_idx(si_grid, gt_grid, delta_grid, GRID_BY_SI)
ip_sparse <- make_ip_sparse(delta_grid, ip_grid, GRID_BY_IP)
zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

message("Lookup tables built.")
message("  N SI grid       = ", length(si_grid))
message("  N l grid        = ", length(l_grid))
message("  N u grid        = ", length(u_grid))
message("  N gt grid       = ", length(gt_grid))
message("  N ip grid       = ", length(ip_grid))
message("  N t grid        = ", length(t_grid))
message("  N delta grid    = ", length(delta_grid))
message("  N z quad        = ", N_Z_QUAD)
message("  N ip sparse     = ", ip_sparse$N_ip_sparse)

make_stan_data_mech <- function(serial_df, r_growth) {
  obs_index <- find_grid_index(
    x = as.numeric(serial_df$serial),
    grid = si_grid,
    by = GRID_BY_SI
  )

  obs_index_pt <- expand_weighted_obs_index(obs_index, serial_df$n)

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
      z_std_w = as.numeric(zq$w)
    ),
    list(
      prior_log_lambda1_mean = prior_log_lambda1_mean,
      prior_log_lambda1_sd = prior_log_lambda1_sd,
      prior_log_lambda_I_mean = prior_log_lambda_I_mean,
      prior_log_lambda_I_sd = prior_log_lambda_I_sd,
      prior_log_r_mean = prior_log_r_mean,
      prior_log_r_sd = prior_log_r_sd,
      prior_log_d_mean = prior_log_d_mean,
      prior_log_d_sd = prior_log_d_sd,
      prior_m_mean = prior_m_mean,
      prior_m_sd = prior_m_sd,
      K_m_mix = nrow(m_prior_mix),
      m_mix_log_w = as.numeric(log(m_prior_mix$weight)),
      m_mix_mu = as.numeric(m_prior_mix$mean),
      m_mix_sd = as.numeric(m_prior_mix$sd),
      prior_log_kappa_mean = prior_log_kappa_mean,
      prior_log_kappa_sd = prior_log_kappa_sd,
      prior_log_sigma_z_mean = prior_log_sigma_z_mean,
      prior_log_sigma_z_sd = prior_log_sigma_z_sd
    )
  )
}

data_mpox <- make_stan_data_mech(serial_w_mpox, r_growth)

mod_mpox <- cmdstan_model(MPOX_STAN_FILE, force_recompile = TRUE)

make_inits_mech <- function() {
  function() {
    list(
      log_lambda1 = prior_log_lambda1_mean,
      log_lambda_I = prior_log_lambda_I_mean,
      log_r = prior_log_r_mean,
      log_d = prior_log_d_mean,
      m = min(max(prior_m_mean, 0.5), 30),
      log_kappa = prior_log_kappa_mean,
      log_sigma_z = prior_log_sigma_z_mean
    )
  }
}

message("\nFitting mpox score==9 SI model with Yang-informed WH priors...")
fit_mpox <- mod_mpox$sample(
  data = data_mpox,
  chains = chains_mpox,
  parallel_chains = parallel_chains_mpox,
  iter_warmup = iter_warmup_mpox,
  iter_sampling = iter_sampling_mpox,
  seed = seed_mpox,
  init = make_inits_mech(),
  adapt_delta = adapt_delta_mpox,
  max_treedepth = max_treedepth_mpox,
  refresh = refresh_mpox,
  save_warmup = FALSE
)

pars_summary <- c(
  "log_lambda1", "log_lambda_I", "log_r", "log_d", "m", "log_kappa", "log_sigma_z", "mean_gen", "sd_gen"
)

summ_mpox <- fit_mpox$summary(variables = pars_summary) %>%
  mutate(dataset = "mpox_score9_yangWHpriors") %>%
  select(dataset, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)
print(summ_mpox)

write.csv(
  summ_mpox,
  file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair_summary.csv"),
  row.names = FALSE
)

fit_mpox$save_object(file.path(OUT_DIR_MPOX, "fit_mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair.rds"))

# Confirm generated quantities needed downstream exist.
gen_vars_mpox <- grep("^generation_prob_grid\\[", posterior::variables(fit_mpox$draws()), value = TRUE)
if (length(gen_vars_mpox) == 0) stop("generation_prob_grid was not found in the fitted object. Check the Stan file.")
message("Confirmed generation_prob_grid exists: entries = ", length(gen_vars_mpox))

extract_log_lik_matrix <- function(fit, variable = "log_lik") as_draws_matrix(fit$draws(variables = variable))

log_lik_mat <- extract_log_lik_matrix(fit_mpox, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit_mpox, "log_lik_grouped")

saveRDS(log_lik_mat, file = file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_log_lik_matrix.rds"))
write.csv(log_lik_mat, file = file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_log_lik_matrix.csv"), row.names = FALSE)
saveRDS(log_lik_grouped_mat, file = file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_log_lik_grouped_matrix.rds"))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)

saveRDS(loo_fit, file = file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_loo.rds"))
saveRDS(waic_fit, file = file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_waic.rds"))

loo_waic_mpox <- data.frame(
  dataset = "mpox_score9_yangWHpriors",
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
  p_waic = waic_fit$estimates["p_waic", "Estimate"],
  row.names = NULL
)
print(loo_waic_mpox)
write.csv(loo_waic_mpox, file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_loo_waic.csv"), row.names = FALSE)

message("\nDiagnostics: mpox score==9 Yang-WH-prior fit")
print(fit_mpox$cmdstan_diagnose())

draws_mpox <- as_draws_df(fit_mpox$draws(variables = pars_summary)) %>%
  mutate(dataset = "mpox_score9_yangWHpriors")

write.csv(
  as.data.frame(draws_mpox),
  file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair_draws.csv"),
  row.names = FALSE
)

p_mean_gen <- ggplot(draws_mpox, aes(x = mean_gen)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = "Mpox score==9 observed-pair model with Yang WH-informed priors"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR_MPOX, "mpox_score9_yangWHpriors_mean_gen_density.png"),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)

message("\nDone. Outputs saved in:")
message("  Yang WH prior outputs: ", OUT_DIR_YANG, "/")
message("  Mpox SI Yang-WH-prior outputs: ", OUT_DIR_MPOX, "/")
message("Fit object:")
message("  ", file.path(OUT_DIR_MPOX, "fit_mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair.rds"))
