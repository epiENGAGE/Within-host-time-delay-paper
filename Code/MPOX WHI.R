# Mox within-host-informed model
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(loo)
})

# Paths
this_file_dir <- setwd("~/Within-host-time-delay-framework/Code")
stan_file <- file.path(this_file_dir, "mpox_whi.stan")
if (!file.exists(stan_file)) stan_file <- "mpox_whi.stan"
if (!file.exists(stan_file)) stop("Could not find mpox_whi.stan")

DATA_URL <- "https://raw.githubusercontent.com/fmiura/MpxSI_2022/main/data/Anonym_All_data_score_indexed.csv"
OUT_DIR <- "~/Within-host-time-delay-framework/Results/Mpox"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

chains <- 4
parallel_chains <- 4
iter_warmup <- 500
iter_sampling <- 500
adapt_delta <- 0.95
max_treedepth <- 12
seed <- 7001
refresh <- 5

# Grids
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

# Load and filter mpox serial interval data
mpox_raw <- read_csv(DATA_URL, show_col_types = FALSE)

mpox_score9 <- mpox_raw %>%
  filter(score == 9, Conf_pair == "high", Conf_sym == "high") %>%
  mutate(serial = as.numeric(SI))

serial_w_mpox <- mpox_score9 %>%
  count(serial, name = "n") %>%
  arrange(serial)

write.csv(mpox_score9, file.path(OUT_DIR, "mpox_score9_high_confidence_pairs.csv"), row.names = FALSE)
write.csv(serial_w_mpox, file.path(OUT_DIR, "mpox_score9_weighted_serial_intervals.csv"), row.names = FALSE)

# Helper functions
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

make_l_t_idx <- function(l_grid, t_grid, by) {
  find_grid_index(l_grid, t_grid, by)
}

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

# Build lookup tables
l_t_idx <- make_l_t_idx(l_grid, t_grid, GRID_BY_LU)
lu_t_idx <- make_lu_t_idx(l_grid, u_grid, t_grid, GRID_BY_LU)
gt_idx <- make_gt_idx(l_grid, u_grid, gt_grid, GRID_BY_LU)
si_g_delta_idx <- make_si_g_delta_idx(si_grid, gt_grid, delta_grid, GRID_BY_SI)
ip_sparse <- make_ip_sparse(delta_grid, ip_grid, GRID_BY_IP)
zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

# Stan data builder
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
      # Weak, non-WH-data-derived priors. These are deliberately broad and only
      # intended to put the model on a plausible acute-infection time scale.
      prior_log_lambda1_mean = log(0.08),
      prior_log_lambda1_sd = 1.5,

      prior_log_lambda_I_mean = log(0.08),
      prior_log_lambda_I_sd = 1.5,

      prior_log_r_mean = log(0.15),
      prior_log_r_sd = 1.0,

      prior_log_d_mean = log(0.20),
      prior_log_d_sd = 1.0,

      prior_m_mean = 8.0,
      prior_m_sd = 5.0,

      prior_log_kappa_mean = 0.0,
      prior_log_kappa_sd = 1.25,

      prior_log_sigma_z_mean = log(0.50),
      prior_log_sigma_z_sd = 1.0
    )
  )
}

data_mpox <- make_stan_data_mech(serial_w_mpox, r_growth)

stopifnot("gt_idx" %in% names(data_mpox))
stopifnot("si_g_delta_idx" %in% names(data_mpox))
stopifnot("N_ip_sparse" %in% names(data_mpox))
stopifnot("park_weight_ip" %in% names(data_mpox))

# Compile
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

# Initial values
make_inits_mech <- function() {
  function() {
    list(
      log_lambda1 = log(0.08),
      log_lambda_I = log(0.08),
      log_r = log(0.15),
      log_d = log(0.20),
      m = 8.0,
      log_kappa = 0,
      log_sigma_z = log(0.50)
    )
  }
}

inits_mech <- make_inits_mech()

# Fit
fit_mpox <- mod$sample(
  data = data_mpox,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  seed = seed,
  init = inits_mech,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = refresh,
  save_warmup = FALSE
)

# Outputs
pars_summary <- c(
  "log_lambda1",
  "log_lambda_I",
  "log_r",
  "log_d",
  "m",
  "log_kappa",
  "log_sigma_z",
  "mean_gen",
  "sd_gen"
)

summ_mpox <- fit_mpox$summary(variables = pars_summary) %>%
  mutate(dataset = "mpox_score9") %>%
  select(dataset, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_mpox)

write.csv(
  summ_mpox,
  file.path(OUT_DIR, "mpox_score9_mechanistic_proportional_incubation_updatedV_obspair_summary.csv"),
  row.names = FALSE
)

fit_mpox$save_object(file.path(OUT_DIR, "fit_mpox_score9_mechanistic_proportional_incubation_updatedV_obspair.rds"))

# Confirm generated quantities needed downstream exist.
gen_vars_mpox <- grep(
  "^generation_prob_grid\\[",
  posterior::variables(fit_mpox$draws()),
  value = TRUE
)

# LOO / WAIC
extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

log_lik_mat <- extract_log_lik_matrix(fit_mpox, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit_mpox, "log_lik_grouped")

saveRDS(log_lik_mat, file = file.path(OUT_DIR, "mpox_score9_log_lik_matrix.rds"))
write.csv(log_lik_mat, file = file.path(OUT_DIR, "mpox_score9_log_lik_matrix.csv"), row.names = FALSE)
saveRDS(log_lik_grouped_mat, file = file.path(OUT_DIR, "mpox_score9_log_lik_grouped_matrix.rds"))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)

saveRDS(loo_fit, file = file.path(OUT_DIR, "mpox_score9_loo.rds"))
saveRDS(waic_fit, file = file.path(OUT_DIR, "mpox_score9_waic.rds"))

loo_waic_mpox <- data.frame(
  dataset = "mpox_score9",
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
write.csv(loo_waic_mpox, file.path(OUT_DIR, "mpox_score9_loo_waic.csv"), row.names = FALSE)

# Diagnostics
print(fit_mpox$cmdstan_diagnose())

# Draws for compact plotting / tables
draws_mpox <- as_draws_df(fit_mpox$draws(variables = pars_summary)) %>%
  mutate(dataset = "mpox_score9")

write.csv(
  as.data.frame(draws_mpox),
  file.path(OUT_DIR, "mpox_score9_mechanistic_proportional_incubation_updatedV_obspair_draws.csv"),
  row.names = FALSE
)

p_mean_gen <- ggplot(draws_mpox, aes(x = mean_gen)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = "Mpox score==9 observed-pair updated-V(t) proportional-incubation model"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mpox_score9_mean_gen_density.png"),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)
