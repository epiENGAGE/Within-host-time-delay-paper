# Noncensored SARS-Cov-2 analysis
setwd("~/Within-host-time-delay-framework/Code")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(loo)
})

# Paths
this_file_dir <- getwd()

stan_file <- file.path(this_file_dir, "covid_noncens.stan")

DATA_DIR <- "~/Within-host-time-delay-framework/Data"
SERIAL_FILE <- file.path(DATA_DIR, "serial-netherlands.xlsx")

OUT_DIR <- "~/Within-host-time-delay-framework/Results/Covid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Controls
r_sgtf <- 0.15
r_nsgtf <- -0.05

GRID_BY_SI <- 0.25
GRID_BY_LU <- 0.25
GRID_BY_IP <- 0.05

SI_MIN <- -10
SI_MAX <- 20

LU_MAX <- 20
GT_MAX <- 2 * LU_MAX

IP_MIN <- 0.05
IP_MAX <- 30

N_Z_QUAD <- 11
Z_MAX_SD <- 3.5

chains <- 4
parallel_chains <- 4
iter_warmup <- 500
iter_sampling <- 500
adapt_delta <- 0.95
max_treedepth <- 12
seed <- 3001
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

# Load serial interval data
make_weighted_serial <- function(df, strain_value, household_value = "within") {
  df %>%
    filter(strain == strain_value, household == household_value) %>%
    group_by(serial) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    arrange(serial)
}

expand_weighted_obs_index <- function(obs_index, counts) {
  rep.int(as.integer(obs_index), as.integer(counts))
}

serialdata <- read_xlsx(SERIAL_FILE)

serial_w_sgtf <- make_weighted_serial(serialdata, "SGTF", "within")
serial_w_nsgtf <- make_weighted_serial(serialdata, "non-SGTF", "within")

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
      prior_log_lambda1_mean = log(0.15),
      prior_log_lambda1_sd = 1.0,
      
      # Proportional-incubation symptom-onset hazard scale:
      # h_IP(t | z) = exp(z) * lambda_I * V(t)
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
  )
}

data_sgtf <- make_stan_data_mech(serial_w_sgtf, r_sgtf)
data_nsgtf <- make_stan_data_mech(serial_w_nsgtf, r_nsgtf)

stopifnot("gt_idx" %in% names(data_sgtf))
stopifnot("si_g_delta_idx" %in% names(data_sgtf))
stopifnot("N_ip_sparse" %in% names(data_sgtf))
stopifnot("park_weight_ip" %in% names(data_sgtf))

# Compile
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

# Initial values
make_inits_mech <- function() {
  function() {
    list(
      log_lambda1 = log(0.15),
      log_lambda_I = log(0.15),
      log_r = log(0.25),
      log_d = log(0.60),
      m = 3.5,
      log_kappa = 0,
      log_sigma_z = log(0.35)
    )
  }
}

inits_mech <- make_inits_mech()

# Fit function
fit_one <- function(data, this_seed, label) {
  message("\nFitting ", label, "...")
  mod$sample(
    data = data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = this_seed,
    init = inits_mech,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = refresh,
    save_warmup = FALSE
  )
}

fit_sgtf <- fit_one(data_sgtf, seed, "SGTF / Omicron")
fit_nsgtf <- fit_one(data_nsgtf, seed + 1, "non-SGTF / Delta")

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

summ_sgtf <- fit_sgtf$summary(variables = pars_summary) %>%
  mutate(strain = "SGTF")

summ_nsgtf <- fit_nsgtf$summary(variables = pars_summary) %>%
  mutate(strain = "non-SGTF")

summ_all <- bind_rows(summ_sgtf, summ_nsgtf) %>%
  select(strain, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_all)

write.csv(
  summ_all,
  file.path(OUT_DIR, "mechanistic_proportional_incubation_updatedV_obspair_summary.csv"),
  row.names = FALSE
)

fit_sgtf$save_object(file.path(OUT_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds"))
fit_nsgtf$save_object(file.path(OUT_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds"))

# Confirm that generated quantities needed downstream exist.
gen_vars_sgtf <- grep(
  "^generation_prob_grid\\[",
  posterior::variables(fit_sgtf$draws()),
  value = TRUE
)

gen_vars_nsgtf <- grep(
  "^generation_prob_grid\\[",
  posterior::variables(fit_nsgtf$draws()),
  value = TRUE
)

# LOO / WAIC
extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

compute_and_save_loo_waic <- function(fit, strain_slug, strain_label) {
  log_lik_mat <- extract_log_lik_matrix(fit, "log_lik")
  log_lik_grouped_mat <- extract_log_lik_matrix(fit, "log_lik_grouped")
  
  saveRDS(
    log_lik_mat,
    file = file.path(OUT_DIR, paste0("mechanistic_proportional_incubation_updatedV_obspair_", strain_slug, "_log_lik_matrix.rds"))
  )
  
  write.csv(
    log_lik_mat,
    file = file.path(OUT_DIR, paste0("mechanistic_proportional_incubation_updatedV_obspair_", strain_slug, "_log_lik_matrix.csv")),
    row.names = FALSE
  )
  
  saveRDS(
    log_lik_grouped_mat,
    file = file.path(OUT_DIR, paste0("mechanistic_proportional_incubation_updatedV_obspair_", strain_slug, "_log_lik_grouped_matrix.rds"))
  )
  
  loo_fit <- loo(log_lik_mat)
  waic_fit <- waic(log_lik_mat)
  
  saveRDS(
    loo_fit,
    file = file.path(OUT_DIR, paste0("mechanistic_proportional_incubation_updatedV_obspair_", strain_slug, "_loo.rds"))
  )
  
  saveRDS(
    waic_fit,
    file = file.path(OUT_DIR, paste0("mechanistic_proportional_incubation_updatedV_obspair_", strain_slug, "_waic.rds"))
  )
  
  data.frame(
    strain = strain_label,
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
}

loo_waic_sgtf <- compute_and_save_loo_waic(fit_sgtf, "sgtf", "SGTF")
loo_waic_nsgtf <- compute_and_save_loo_waic(fit_nsgtf, "nsgtf", "non-SGTF")
loo_waic_all <- bind_rows(loo_waic_sgtf, loo_waic_nsgtf)

print(loo_waic_all)

write.csv(
  loo_waic_all,
  file.path(OUT_DIR, "mechanistic_proportional_incubation_updatedV_obspair_loo_waic.csv"),
  row.names = FALSE
)

# Draws for compact plotting / tables
draws_sgtf <- as_draws_df(fit_sgtf$draws(variables = pars_summary)) %>%
  mutate(strain = "SGTF")

draws_nsgtf <- as_draws_df(fit_nsgtf$draws(variables = pars_summary)) %>%
  mutate(strain = "non-SGTF")

draws_all <- bind_rows(draws_sgtf, draws_nsgtf)

write.csv(
  as.data.frame(draws_all),
  file.path(OUT_DIR, "mechanistic_proportional_incubation_updatedV_obspair_draws.csv"),
  row.names = FALSE
)

# Posterior mean generation interval plot
p_mean_gen <- ggplot(draws_all, aes(x = mean_gen, fill = strain)) +
  geom_density(alpha = 0.35) +
  labs(
    x = "Mean generation interval, days",
    y = "Posterior density",
    title = "Observed-pair updated-V(t) proportional-incubation mechanistic model"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "mechanistic_proportional_incubation_updatedV_obspair_mean_gen_density.png"),
  p_mean_gen,
  width = 7,
  height = 4.5,
  dpi = 300
)
