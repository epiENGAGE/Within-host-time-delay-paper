# Covid analysis accounting for censoring
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readxl)
  library(loo)
})

# Paths
stan_file <- "~/Within-host-time-delay-framework/Code/covid_cens.stan"

serial_file <- "~/Within-host-time-delay-framework/Data/serial-netherlands.xlsx"

out_dir <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

# Controls

# Same growth correction values used in the previous COVID application.
r_sgtf <- 0.15
r_nsgtf <- -0.05

GRID_BY_LU <- 0.25
GRID_BY_IP <- 0.05

LU_MAX <- 20
GT_MAX <- 2 * LU_MAX
IP_MIN <- 0.05
IP_MAX <- 30

SI_INT_MIN <- -10
SI_INT_MAX <- 20
serial_cat <- SI_INT_MIN:SI_INT_MAX

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

# Helpers
find_grid_index <- function(x, grid, by, tol = 1e-8) {
  idx <- round((x - grid[1]) / by) + 1L
  if (any(idx < 1L | idx > length(grid))) {
    bad <- x[idx < 1L | idx > length(grid)]
    stop("Value outside grid: ", paste(head(bad, 10), collapse = ", "))
  }
  if (any(abs(grid[idx] - x) > tol)) {
    bad <- x[abs(grid[idx] - x) > tol]
    stop("Value not on grid: ", paste(head(bad, 10), collapse = ", "))
  }
  as.integer(idx)
}

make_z_quad <- function(n_z = 11, z_max_sd = 3.5) {
  x <- seq(-z_max_sd, z_max_sd, length.out = n_z)
  dx <- if (length(x) > 1L) diff(x)[1L] else 1
  w_raw <- dnorm(x) * dx
  list(x = x, w = w_raw / sum(w_raw))
}

make_weighted_serial <- function(df, strain_value, household_value = "within") {
  df %>%
    filter(.data$strain == strain_value, .data$household == household_value) %>%
    group_by(.data$serial) %>%
    dplyr::summarize(n = sum(.data$n), .groups = "drop") %>%
    arrange(.data$serial)
}

make_counts_on_support <- function(serial_df, serial_cat) {
  out <- integer(length(serial_cat))
  idx <- match(as.integer(serial_df$serial), serial_cat)
  if (any(is.na(idx))) {
    bad <- serial_df$serial[is.na(idx)]
    stop("Observed serial interval outside support: ", paste(bad, collapse = ", "))
  }
  out[idx] <- as.integer(serial_df$n)
  out
}

expand_obs_cat <- function(counts) {
  rep.int(seq_along(counts), counts)
}

# Grids and lookup tables
l_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
u_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
gt_grid <- seq(0, GT_MAX, by = GRID_BY_LU)
ip_grid <- seq(IP_MIN, IP_MAX, by = GRID_BY_IP)
t_grid <- seq(0, GT_MAX, by = GRID_BY_LU)

delta_grid <- seq(min(ip_grid) - max(ip_grid), max(ip_grid) - min(ip_grid), by = GRID_BY_IP)

lu_t_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
gt_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
for (jl in seq_along(l_grid)) {
  for (ju in seq_along(u_grid)) {
    lu <- l_grid[jl] + u_grid[ju]
    lu_t_idx[jl, ju] <- find_grid_index(lu, t_grid, GRID_BY_LU)
    gt_idx[jl, ju] <- find_grid_index(lu, gt_grid, GRID_BY_LU)
  }
}

# All IP1/IP2 pairs; this is the explicit marginalisation over both onset times.
pair_delta <- integer(length(ip_grid) * length(ip_grid))
pair_i1 <- integer(length(pair_delta))
pair_i2 <- integer(length(pair_delta))
qq <- 0L
for (ii in seq_along(ip_grid)) {
  for (jj in seq_along(ip_grid)) {
    qq <- qq + 1L
    delta <- ip_grid[jj] - ip_grid[ii]
    pair_delta[qq] <- find_grid_index(delta, delta_grid, GRID_BY_IP)
    pair_i1[qq] <- ii
    pair_i2[qq] <- jj
  }
}

coarse_cat <- integer(0)
coarse_g <- integer(0)
coarse_delta <- integer(0)
coarse_w <- numeric(0)
for (kc in seq_along(serial_cat)) {
  y <- serial_cat[kc]
  for (kg in seq_along(gt_grid)) {
    x_no_delta <- gt_grid[kg]
    w <- pmax(0, 1 - abs(x_no_delta + delta_grid - y))
    ok <- which(w > 0)
    if (length(ok) > 0) {
      coarse_cat <- c(coarse_cat, rep.int(kc, length(ok)))
      coarse_g <- c(coarse_g, rep.int(kg, length(ok)))
      coarse_delta <- c(coarse_delta, ok)
      coarse_w <- c(coarse_w, w[ok])
    }
  }
}

zq <- make_z_quad(N_Z_QUAD, Z_MAX_SD)

# Data
serialdata <- read_xlsx(serial_file)

serial_w_sgtf <- make_weighted_serial(serialdata, "SGTF", "within")
serial_w_nsgtf <- make_weighted_serial(serialdata, "non-SGTF", "within")

make_stan_data <- function(serial_df, r_growth) {
  counts <- make_counts_on_support(serial_df, serial_cat)
  obs_cat_pt <- expand_obs_cat(counts)
  
  list(
    N_cat = length(serial_cat),
    counts = as.integer(counts),
    serial_cat = as.numeric(serial_cat),
    
    N_obs = length(obs_cat_pt),
    obs_cat_pt = as.integer(obs_cat_pt),
    
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
    
    N_delta = length(delta_grid),
    delta_grid = as.numeric(delta_grid),
    
    lu_t_idx = lu_t_idx,
    gt_idx = gt_idx,
    
    N_ip_pair = length(pair_delta),
    pair_delta = as.integer(pair_delta),
    pair_i1 = as.integer(pair_i1),
    pair_i2 = as.integer(pair_i2),
    
    N_coarse = length(coarse_w),
    coarse_cat = as.integer(coarse_cat),
    coarse_g = as.integer(coarse_g),
    coarse_delta = as.integer(coarse_delta),
    coarse_w = as.numeric(coarse_w),
    
    park_weight_ip = as.numeric(exp(-r_growth * ip_grid)),
    
    N_z = N_Z_QUAD,
    z_std_grid = as.numeric(zq$x),
    z_std_w = as.numeric(zq$w),
    
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
}

data_sgtf <- make_stan_data(serial_w_sgtf, r_sgtf)
data_nsgtf <- make_stan_data(serial_w_nsgtf, r_nsgtf)

# Compile and sample
mod <- cmdstan_model(stan_file, force_recompile = TRUE)

make_inits <- function() {
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

fit_one <- function(data, this_seed, label) {
  message("\nFitting ", label, "...")
  mod$sample(
    data = data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = this_seed,
    init = make_inits(),
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

summ_sgtf <- fit_sgtf$summary(variables = pars_summary) %>% mutate(strain = "SGTF")
summ_nsgtf <- fit_nsgtf$summary(variables = pars_summary) %>% mutate(strain = "non-SGTF")

summ_all <- bind_rows(summ_sgtf, summ_nsgtf) %>%
  select(strain, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_all)
write.csv(summ_all, file.path(out_dir, "date_coarsened_joint_obspair_summary.csv"), row.names = FALSE)

fit_sgtf$save_object(file.path(out_dir, "fit_date_coarsened_joint_obspair_sgtf.rds"))
fit_nsgtf$save_object(file.path(out_dir, "fit_date_coarsened_joint_obspair_nsgtf.rds"))

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

compute_and_save_loo_waic <- function(fit, strain_slug, strain_label) {
  log_lik_mat <- extract_log_lik_matrix(fit, "log_lik")
  log_lik_grouped_mat <- extract_log_lik_matrix(fit, "log_lik_grouped")
  
  saveRDS(log_lik_mat, file.path(out_dir, paste0(strain_slug, "_log_lik_matrix.rds")))
  saveRDS(log_lik_grouped_mat, file.path(out_dir, paste0(strain_slug, "_log_lik_grouped_matrix.rds")))
  
  loo_fit <- loo(log_lik_mat)
  waic_fit <- waic(log_lik_mat)
  
  saveRDS(loo_fit, file.path(out_dir, paste0(strain_slug, "_loo.rds")))
  saveRDS(waic_fit, file.path(out_dir, paste0(strain_slug, "_waic.rds")))
  
  data.frame(
    strain = strain_label,
    n_pointwise_obs = ncol(log_lik_mat),
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_fit$estimates["elpd_loo", "SE"],
    looic = loo_fit$estimates["looic", "Estimate"],
    se_looic = loo_fit$estimates["looic", "SE"],
    elpd_waic = waic_fit$estimates["elpd_waic", "Estimate"],
    se_elpd_waic = waic_fit$estimates["elpd_waic", "SE"],
    waic = waic_fit$estimates["waic", "Estimate"],
    se_waic = waic_fit$estimates["waic", "SE"]
  )
}

ic_all <- bind_rows(
  compute_and_save_loo_waic(fit_sgtf, "sgtf", "SGTF"),
  compute_and_save_loo_waic(fit_nsgtf, "nsgtf", "non-SGTF")
)
print(ic_all)
write.csv(ic_all, file.path(out_dir, "date_coarsened_joint_obspair_loo_waic.csv"), row.names = FALSE)
