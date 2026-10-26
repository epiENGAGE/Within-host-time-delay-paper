# ============================================================
# Mpox high-confidence-pair lognormal generation-interval model
# FAST exact-on-grid convolution implementation
#
# This replaces the slow quad_terms implementation with a discrete
# incubation-difference convolution:
#   SI = GT + IP2 - IP1
#   delta = IP2 - IP1
#   f_SI(s) = integral f_GT(s - delta) f_delta(delta) d_delta
#
# This is algebraically equivalent to the brute-force quadrature on the
# same grids when rho = 0. The script stops if rho != 0 so that the model
# is not silently changed.
#
# Fits Miura et al. mpox data restricted to score == 9 pairs:
#   Conf_pair == high and Conf_sym == high
#
# Paired Stan file:
#   mpox_lognormal_generation_fastconv.stan
#
# Outputs:
#   stan_outputs_mpox_lognormal_score9_fastconv/
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
  library(loo)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

this_file_dir <- getwd()
stan_file <- file.path(this_file_dir, "mpox_lognormal_generation_fastconv.stan")
if (!file.exists(stan_file)) stan_file <- "mpox_lognormal_generation_fastconv.stan"
if (!file.exists(stan_file)) stop("Could not find mpox_lognormal_generation_fastconv.stan. Put it in the working directory or edit stan_file.")

DATA_URL <- "https://raw.githubusercontent.com/fmiura/MpxSI_2022/main/data/Anonym_All_data_score_indexed.csv"
OUT_DIR <- "stan_outputs_mpox_lognormal_score9_fastconv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Using Stan file:  ", stan_file)
message("Using mpox data:  ", DATA_URL)
message("Output directory: ", OUT_DIR)

# ------------------------------------------------------------
# Controls
# ------------------------------------------------------------

# Serial-interval normalization grid. Keep this the same as the slow script
# if you want exactly comparable numbers.
GRID_BY_SI <- 0.25
SI_MIN <- -30
SI_MAX <- 60

# Incubation-period quadrature grid. For exact-on-grid convolution with the
# current Stan file, use the same step as the serial grid.
GRID_BY_IP <- 0.25
IP_MIN <- 0.05
IP_MAX <- 45

if (abs(GRID_BY_SI - GRID_BY_IP) > 1e-12) {
  stop("For this fast exact-on-grid implementation, set GRID_BY_SI == GRID_BY_IP.")
}

# Park correction epidemic growth rate. This is not a within-host growth rate.
# Default 0 means no Park correction is applied.
park_growth_rate <- 0.0

# Fixed incubation-period distribution used by the lognormal comparator.
# These are deliberately editable and are NOT estimated from within-host data here.
inc_mean_days <- 8.5
inc_sd_days <- 4.5

# The fast convolution is exactly equivalent to the original brute-force
# quadrature only for rho = 0, because nonzero rho makes the generation-time
# density depend on IP1 before marginalizing over IP2 - IP1.
rho <- 0.0
if (abs(rho) > 1e-12) {
  stop("Fast convolution version is exact only for rho = 0. Use the slow quad_terms version for rho != 0.")
}

# Priors on lognormal generation-interval parameters. These are on the log-scale
# parameters logmean_gen and logsd_gen, not directly on mean/sd in days.
prior_gen_mean_days <- 12
prior_gen_sd_days <- 7
prior_logmean_gen_sd <- 1.0
prior_logsd_gen_sd <- 0.75

chains <- 4
parallel_chains <- 4
iter_warmup <- 500
iter_sampling <- 500
adapt_delta <- 0.95
max_treedepth <- 12
seed <- 4001
refresh <- 5

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

lnorm_params_from_mean_sd <- function(mean, sd) {
  if (mean <= 0 || sd <= 0) stop("mean and sd must be positive")
  sigma2 <- log(1 + (sd^2 / mean^2))
  sigma <- sqrt(sigma2)
  mu <- log(mean) - 0.5 * sigma2
  list(mu = mu, sigma = sigma)
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

make_serial_score9 <- function(df) {
  df %>%
    filter(score == 9, Conf_pair == "high", Conf_sym == "high") %>%
    group_by(SI) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(SI) %>%
    rename(serial = SI)
}

make_incubation_difference_density <- function(ip_grid, dx_ip, logmean_inc, logsd_inc, park_growth_rate) {
  # delta = IP2 - IP1. The density is approximated by
  # C(delta) = integral exp(-r * IP1) f_IP(IP1) f_IP(IP1 + delta) dIP1.
  # The final serial likelihood is normalized over the SI grid, so this does
  # not need to integrate to 1 when a Park correction or truncation is present.

  delta_grid <- seq(min(ip_grid) - max(ip_grid), max(ip_grid) - min(ip_grid), by = dx_ip)
  f_ip <- dlnorm(ip_grid, meanlog = logmean_inc, sdlog = logsd_inc)
  park_w <- exp(-park_growth_rate * ip_grid)

  C_delta_density <- numeric(length(delta_grid))

  for (dd in seq_along(delta_grid)) {
    target_ip2 <- ip_grid + delta_grid[dd]
    idx2 <- round((target_ip2 - ip_grid[1]) / dx_ip) + 1L
    ok <- idx2 >= 1L & idx2 <= length(ip_grid)

    if (any(ok)) {
      # Guard against floating point index slips.
      ok_idx <- which(ok)
      ok[ok_idx] <- abs(ip_grid[idx2[ok_idx]] - target_ip2[ok_idx]) < 1e-8
    }

    if (any(ok)) {
      C_delta_density[dd] <- sum(park_w[ok] * f_ip[ok] * f_ip[idx2[ok]]) * dx_ip
    }
  }

  list(
    delta_grid = as.numeric(delta_grid),
    C_delta_density = as.numeric(C_delta_density)
  )
}

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

# ------------------------------------------------------------
# Build grids and fixed distributions
# ------------------------------------------------------------

si_norm_grid <- seq(SI_MIN, SI_MAX, by = GRID_BY_SI)
ip_grid <- seq(IP_MIN, IP_MAX, by = GRID_BY_IP)

inc_par <- lnorm_params_from_mean_sd(inc_mean_days, inc_sd_days)
prior_gen_par <- lnorm_params_from_mean_sd(prior_gen_mean_days, prior_gen_sd_days)

message("Fixed incubation distribution:")
message("  mean days = ", inc_mean_days)
message("  sd days   = ", inc_sd_days)
message("  meanlog   = ", signif(inc_par$mu, 4))
message("  sdlog     = ", signif(inc_par$sigma, 4))
message("rho = ", rho)
message("park_growth_rate = ", park_growth_rate)

message("Building incubation-difference convolution...")
conv <- make_incubation_difference_density(
  ip_grid = ip_grid,
  dx_ip = GRID_BY_IP,
  logmean_inc = inc_par$mu,
  logsd_inc = inc_par$sigma,
  park_growth_rate = park_growth_rate
)
message("  N_delta = ", length(conv$delta_grid))
message("  delta range = ", min(conv$delta_grid), " to ", max(conv$delta_grid))
message("  nonzero C_delta entries = ", sum(conv$C_delta_density > 0))

# ------------------------------------------------------------
# Load mpox serial-interval data
# ------------------------------------------------------------

mpox_raw <- read.csv(DATA_URL, stringsAsFactors = FALSE)
serial_w <- make_serial_score9(mpox_raw)

message("Mpox score == 9 data:")
message("  unique SI values = ", nrow(serial_w))
message("  total pairs      = ", sum(serial_w$n))
message("  observed SI range = ", min(serial_w$serial), " to ", max(serial_w$serial))

obs_index <- find_grid_index(as.numeric(serial_w$serial), si_norm_grid, GRID_BY_SI)
obs_index_obs <- expand_weighted_obs_index(obs_index, serial_w$n)

# ------------------------------------------------------------
# Stan data
# ------------------------------------------------------------

stan_data <- list(
  N_norm = length(si_norm_grid),
  si_norm_grid = as.numeric(si_norm_grid),
  si_norm_step = GRID_BY_SI,

  N_si = nrow(serial_w),
  obs_index = as.integer(obs_index),
  counts = as.integer(serial_w$n),

  N_obs = length(obs_index_obs),
  obs_index_obs = as.integer(obs_index_obs),

  N_delta = length(conv$delta_grid),
  delta_grid = as.numeric(conv$delta_grid),
  delta_step = GRID_BY_IP,
  C_delta_density = as.numeric(conv$C_delta_density),

  prior_logmean_gen_mean = prior_gen_par$mu,
  prior_logmean_gen_sd = prior_logmean_gen_sd,
  prior_logsd_gen_mean = log(prior_gen_par$sigma),
  prior_logsd_gen_sd = prior_logsd_gen_sd
)

# ------------------------------------------------------------
# Compile and fit
# ------------------------------------------------------------

mod <- cmdstan_model(stan_file, force_recompile = TRUE)

make_inits <- function() {
  function() {
    list(
      logmean_gen = prior_gen_par$mu,
      logsd_gen = prior_gen_par$sigma
    )
  }
}

message("\nFitting mpox lognormal fast-convolution comparator...")
fit <- mod$sample(
  data = stan_data,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  seed = seed,
  init = make_inits(),
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = refresh,
  save_warmup = FALSE
)

# ------------------------------------------------------------
# Outputs: summaries
# ------------------------------------------------------------

pars_summary <- c(
  "logmean_gen",
  "logsd_gen",
  "mean_gen",
  "sd_gen"
)

summ <- fit$summary(variables = pars_summary) %>%
  mutate(model = "Lognormal", dataset = "mpox_score9") %>%
  select(dataset, model, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ)

write.csv(
  summ,
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_summary.csv"),
  row.names = FALSE
)

fit$save_object(file.path(OUT_DIR, "fit_mpox_lognormal_score9_fastconv.rds"))

# ------------------------------------------------------------
# LOO / WAIC
# ------------------------------------------------------------

log_lik_mat <- extract_log_lik_matrix(fit, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit, "log_lik_grouped")

saveRDS(log_lik_mat, file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_log_lik_matrix.rds"))
write.csv(log_lik_mat, file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_log_lik_matrix.csv"), row.names = FALSE)
saveRDS(log_lik_grouped_mat, file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_log_lik_grouped_matrix.rds"))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)

saveRDS(loo_fit, file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_loo.rds"))
saveRDS(waic_fit, file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_waic.rds"))

loo_waic <- data.frame(
  dataset = "mpox_score9",
  model = "Lognormal",
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

print(loo_waic)

write.csv(
  loo_waic,
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_loo_waic.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

message("\nDiagnostics: mpox lognormal fast-convolution comparator")
print(fit$cmdstan_diagnose())

# ------------------------------------------------------------
# Draws for compact plotting / tables
# ------------------------------------------------------------

draws <- as_draws_df(fit$draws(variables = pars_summary)) %>%
  mutate(model = "Lognormal", dataset = "mpox_score9")

write.csv(
  as.data.frame(draws),
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_draws.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Posterior fitted serial density grid
# ------------------------------------------------------------

serial_density_draws <- as_draws_matrix(fit$draws(variables = "serial_density_norm_grid"))
serial_density_mean <- colMeans(serial_density_draws)
serial_density_q05 <- apply(serial_density_draws, 2, quantile, probs = 0.05)
serial_density_q50 <- apply(serial_density_draws, 2, quantile, probs = 0.50)
serial_density_q95 <- apply(serial_density_draws, 2, quantile, probs = 0.95)

serial_fit_grid <- data.frame(
  serial = si_norm_grid,
  mean_density = serial_density_mean,
  q05_density = serial_density_q05,
  median_density = serial_density_q50,
  q95_density = serial_density_q95
)

write.csv(
  serial_fit_grid,
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_serial_density_grid.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Quick plots
# ------------------------------------------------------------

p_mean_gen <- ggplot(draws, aes(x = mean_gen)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = "Mpox score == 9 lognormal generation-interval comparator"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_mean_gen_density.png"),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)

obs_plot_df <- serial_w %>% mutate(freq = n)
serial_fit_plot <- serial_fit_grid %>%
  filter(serial >= min(serial_w$serial) - 10, serial <= max(serial_w$serial) + 20)

p_fit <- ggplot() +
  geom_col(data = obs_plot_df, aes(x = serial, y = freq), alpha = 0.35, width = 0.8) +
  geom_line(
    data = serial_fit_plot,
    aes(x = serial, y = median_density * sum(serial_w$n) * GRID_BY_SI),
    linewidth = 0.8
  ) +
  geom_ribbon(
    data = serial_fit_plot,
    aes(
      x = serial,
      ymin = q05_density * sum(serial_w$n) * GRID_BY_SI,
      ymax = q95_density * sum(serial_w$n) * GRID_BY_SI
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Serial interval, days",
    y = "Frequency",
    title = "Observed and fitted mpox serial intervals: lognormal comparator"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mpox_lognormal_score9_fastconv_serial_fit.png"),
  p_fit,
  width = 7,
  height = 4.5,
  dpi = 300
)

message("\nDone. Outputs saved in ", OUT_DIR, "/")
message("Fit object:")
message("  ", file.path(OUT_DIR, "fit_mpox_lognormal_score9_fastconv.rds"))
