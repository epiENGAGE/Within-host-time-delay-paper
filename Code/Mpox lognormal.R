# Mpox lognormal analysis
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
  library(loo)
})


# Small local helper used in path detection.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Paths
this_file_dir <- "~/Within-host-time-delay-framework/Code"
stan_file <- file.path(this_file_dir, "mpox_lognormal.stan")

DATA_URL <- "https://raw.githubusercontent.com/fmiura/MpxSI_2022/main/data/Anonym_All_data_score_indexed.csv"
OUT_DIR <- "~/Within-host-time-delay-framework/Results/Mpox"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Controls
GRID_BY_SI <- 0.25
SI_MIN <- -30
SI_MAX <- 60

# Incubation-period quadrature grid.
GRID_BY_IP <- 0.25
IP_MIN <- 0.05
IP_MAX <- 45

r_growth <- 0.0

inc_mean_days <- 8.5
inc_sd_days <- 4.5

rho <- 0.0

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

# Helper functions
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

make_quadrature_terms <- function(si_grid, ip_grid, dx_ip, logmean_inc, logsd_inc, r_growth) {
  s_id <- integer(0)
  const_lp <- numeric(0)
  log_inc1 <- numeric(0)
  log_gen <- numeric(0)

  log_dx2 <- 2 * log(dx_ip)

  for (ks in seq_along(si_grid)) {
    s <- si_grid[ks]

    for (i1 in seq_along(ip_grid)) {
      inc1 <- ip_grid[i1]
      log_f_inc1 <- dlnorm(inc1, meanlog = logmean_inc, sdlog = logsd_inc, log = TRUE)

      for (i2 in seq_along(ip_grid)) {
        inc2 <- ip_grid[i2]
        gen <- s + inc1 - inc2

        if (is.finite(gen) && gen > 0) {
          s_id <- c(s_id, ks)
          log_f_inc2 <- dlnorm(inc2, meanlog = logmean_inc, sdlog = logsd_inc, log = TRUE)

          # alpha1 = -inc1, alpha2 = s - inc2, generation = alpha2 - alpha1 = s + inc1 - inc2.
          # Park correction contribution is r * alpha1 = -r * inc1.
          const_lp <- c(const_lp, -r_growth * inc1 + log_f_inc1 + log_f_inc2 + log_dx2)
          log_inc1 <- c(log_inc1, log(inc1))
          log_gen <- c(log_gen, log(gen))
        }
      }
    }
  }

  list(
    N_terms = length(s_id),
    s_id = as.integer(s_id),
    const_lp = as.numeric(const_lp),
    log_inc1 = as.numeric(log_inc1),
    log_gen = as.numeric(log_gen)
  )
}

# Build grids
si_norm_grid <- seq(SI_MIN, SI_MAX, by = GRID_BY_SI)
ip_grid <- seq(IP_MIN, IP_MAX, by = GRID_BY_IP)

inc_par <- lnorm_params_from_mean_sd(inc_mean_days, inc_sd_days)
prior_gen_par <- lnorm_params_from_mean_sd(prior_gen_mean_days, prior_gen_sd_days)

# Load mpox serial-interval data
mpox_raw <- read.csv(DATA_URL, stringsAsFactors = FALSE)
serial_w <- make_serial_score9(mpox_raw)

obs_index <- find_grid_index(as.numeric(serial_w$serial), si_norm_grid, GRID_BY_SI)
obs_index_obs <- expand_weighted_obs_index(obs_index, serial_w$n)

# Precompute quadrature terms
quad_terms <- make_quadrature_terms(
  si_grid = si_norm_grid,
  ip_grid = ip_grid,
  dx_ip = GRID_BY_IP,
  logmean_inc = inc_par$mu,
  logsd_inc = inc_par$sigma,
  r_growth = r_growth
)

stan_data <- c(
  list(
    N_norm = length(si_norm_grid),
    si_norm_grid = as.numeric(si_norm_grid),
    si_norm_step = GRID_BY_SI,

    N_si = nrow(serial_w),
    obs_index = as.integer(obs_index),
    counts = as.integer(serial_w$n),

    N_obs = length(obs_index_obs),
    obs_index_obs = as.integer(obs_index_obs),

    logmean_inc = inc_par$mu,
    logsd_inc = inc_par$sigma,
    rho = rho,

    prior_logmean_gen_mean = prior_gen_par$mu,
    prior_logmean_gen_sd = prior_logmean_gen_sd,
    prior_logsd_gen_mean = prior_gen_par$sigma,
    prior_logsd_gen_sd = prior_logsd_gen_sd
  ),
  quad_terms
)

# Compile and fit
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

inits <- function() {
  list(
    logmean_gen = prior_gen_par$mu,
    logsd_gen = prior_gen_par$sigma
  )
}

fit_mpox <- mod$sample(
  data = stan_data,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  seed = seed,
  init = inits,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = refresh,
  save_warmup = FALSE
)

# Outputs: summary
pars_summary <- c(
  "logmean_gen",
  "logsd_gen",
  "mean_gen",
  "sd_gen"
)

summ <- fit_mpox$summary(variables = pars_summary) %>%
  mutate(dataset = "mpox_score9") %>%
  select(dataset, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ)

write.csv(
  summ,
  file.path(OUT_DIR, "mpox_lognormal_generation_score9_summary.csv"),
  row.names = FALSE
)

fit_mpox$save_object(file.path(OUT_DIR, "fit_mpox_lognormal_generation_score9.rds"))

# LOO / WAIC
extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

log_lik_mat <- extract_log_lik_matrix(fit_mpox, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit_mpox, "log_lik_grouped")

saveRDS(log_lik_mat, file.path(OUT_DIR, "mpox_lognormal_generation_score9_log_lik_matrix.rds"))
write.csv(log_lik_mat, file.path(OUT_DIR, "mpox_lognormal_generation_score9_log_lik_matrix.csv"), row.names = FALSE)
saveRDS(log_lik_grouped_mat, file.path(OUT_DIR, "mpox_lognormal_generation_score9_log_lik_grouped_matrix.rds"))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)

saveRDS(loo_fit, file.path(OUT_DIR, "mpox_lognormal_generation_score9_loo.rds"))
saveRDS(waic_fit, file.path(OUT_DIR, "mpox_lognormal_generation_score9_waic.rds"))

loo_waic <- data.frame(
  dataset = "mpox_score9",
  model = "lognormal_generation",
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
write.csv(loo_waic, file.path(OUT_DIR, "mpox_lognormal_generation_score9_loo_waic.csv"), row.names = FALSE)

# Diagnostics
print(fit_mpox$cmdstan_diagnose())

# Draws and plots
draws <- as_draws_df(fit_mpox$draws(variables = pars_summary)) %>%
  mutate(dataset = "mpox_score9", model = "lognormal_generation")

write.csv(
  as.data.frame(draws),
  file.path(OUT_DIR, "mpox_lognormal_generation_score9_draws.csv"),
  row.names = FALSE
)

p_mean_gen <- ggplot(draws, aes(x = mean_gen)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = "Mpox score == 9 lognormal generation-interval model"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mpox_lognormal_generation_score9_mean_gen_density.png"),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)

serial_density_draws <- as_draws_matrix(fit_mpox$draws(variables = "serial_density_norm_grid"))
serial_density_summary <- data.frame(
  si = si_norm_grid,
  median = apply(serial_density_draws, 2, median),
  q025 = apply(serial_density_draws, 2, quantile, probs = 0.025),
  q975 = apply(serial_density_draws, 2, quantile, probs = 0.975)
)
write.csv(serial_density_summary, file.path(OUT_DIR, "mpox_lognormal_generation_score9_serial_density_summary.csv"), row.names = FALSE)

obs_plot_df <- serial_w %>% mutate(prop = n / sum(n))

p_si <- ggplot() +
  geom_col(data = obs_plot_df, aes(x = serial, y = prop), width = 0.85, alpha = 0.45) +
  geom_line(data = serial_density_summary, aes(x = si, y = median * GRID_BY_SI), linewidth = 0.8) +
  geom_ribbon(data = serial_density_summary, aes(x = si, ymin = q025 * GRID_BY_SI, ymax = q975 * GRID_BY_SI), alpha = 0.2) +
  labs(
    x = "Serial interval, days",
    y = "Probability mass per grid cell",
    title = "Mpox score == 9 observed serial intervals and fitted lognormal model"
  ) +
  coord_cartesian(xlim = c(min(serial_w$serial) - 5, max(serial_w$serial) + 10)) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mpox_lognormal_generation_score9_serial_fit.png"),
  p_si,
  width = 7,
  height = 4.5,
  dpi = 300
)
