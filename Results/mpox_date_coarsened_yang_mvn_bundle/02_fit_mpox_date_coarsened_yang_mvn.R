# ============================================================
# 02_fit_mpox_date_coarsened_yang_mvn.R
#
# Fits the mpox high-confidence observed-pair serial-interval model using:
#   - a COVID-style calendar-date coarsened observation likelihood;
#   - latent LP, U, IP1, IP2 integration;
#   - Yang-informed MVN prior ONLY for theta = [log_r, log_d, m, log_kappa].
#
# It does NOT use the previous Yang/broad robust mixture prior.
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readr)
  library(loo)
})

this_file_dir <- getwd()
MPOX_STAN_FILE <- file.path(this_file_dir, "mpox_mechanistic_date_coarsened_yang_mvn.stan")
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- "mpox_mechanistic_date_coarsened_yang_mvn.stan"
if (!file.exists(MPOX_STAN_FILE)) MPOX_STAN_FILE <- file.path("/mnt/data", "mpox_mechanistic_date_coarsened_yang_mvn.stan")
if (!file.exists(MPOX_STAN_FILE)) stop("Could not find mpox_mechanistic_date_coarsened_yang_mvn.stan")

MPOX_DATA_URL <- "https://raw.githubusercontent.com/fmiura/MpxSI_2022/main/data/Anonym_All_data_score_indexed.csv"
OUT_DIR <- "mpox_score9_date_coarsened_Yang_MVN_outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Read Yang-informed prior built by 01_fit_yang_construct_prior.R
prior_marginal_file <- file.path(OUT_DIR, "constructed_theta_prior_marginals.csv")
prior_cov_file <- file.path(OUT_DIR, "constructed_theta_prior_cov.csv")
if (!file.exists(prior_marginal_file) || !file.exists(prior_cov_file)) {
  stop("Missing constructed Yang prior files. Run 01_fit_yang_construct_prior.R first.")
}

theta_marg <- read_csv(prior_marginal_file, show_col_types = FALSE)
theta_cov <- as.matrix(read.csv(prior_cov_file, row.names = 1, check.names = FALSE))

theta_order <- c("log_r", "log_d", "m", "log_kappa")
theta_marg <- theta_marg %>% slice(match(theta_order, parameter))
prior_theta_mean <- theta_marg$prior_mean
prior_theta_cov <- theta_cov[theta_order, theta_order]

# Event-hazard and heterogeneity priors remain weak: Yang viral-load data do not identify these.
prior_log_lambda1_mean <- log(0.08)
prior_log_lambda1_sd <- 1.5
prior_log_lambda_I_mean <- log(0.08)
prior_log_lambda_I_sd <- 1.5
prior_log_sigma_z_mean <- log(0.50)
prior_log_sigma_z_sd <- 1.0

# Mpox SI model controls
r_growth <- 0.0
GRID_BY_LU <- 0.25
GRID_BY_IP <- 0.05
SI_INT_MIN <- -15
SI_INT_MAX <- 40
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

make_z_quad <- function(n_z = 11, z_max_sd = 3.5) {
  x <- seq(-z_max_sd, z_max_sd, length.out = n_z)
  dx <- if (length(x) > 1L) diff(x)[1L] else 1
  w_raw <- dnorm(x) * dx
  list(x = x, w = w_raw / sum(w_raw))
}

extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  as_draws_matrix(fit$draws(variables = variable))
}

# ------------------------------------------------------------
# Build grids and COVID-style date-coarsening lookup tables
# ------------------------------------------------------------

serial_cat <- SI_INT_MIN:SI_INT_MAX
l_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
u_grid <- seq(0, LU_MAX, by = GRID_BY_LU)
gt_grid <- seq(0, GT_MAX, by = GRID_BY_LU)
ip_grid <- seq(IP_MIN, IP_MAX, by = GRID_BY_IP)
t_grid <- seq(0, GT_MAX, by = GRID_BY_LU)

delta_grid <- seq(min(ip_grid) - max(ip_grid), max(ip_grid) - min(ip_grid), by = GRID_BY_IP)

message("Building date-coarsened lookup tables...")

lu_t_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
gt_idx <- matrix(NA_integer_, nrow = length(l_grid), ncol = length(u_grid))
for (jl in seq_along(l_grid)) {
  for (ju in seq_along(u_grid)) {
    lu <- l_grid[jl] + u_grid[ju]
    lu_t_idx[jl, ju] <- find_grid_index(lu, t_grid, GRID_BY_LU)
    gt_idx[jl, ju] <- find_grid_index(lu, gt_grid, GRID_BY_LU)
  }
}

# Explicit marginalisation over both incubation periods.
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

# Sparse triangular calendar-date coarsening kernel.
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

message("Lookup tables built.")
message("  N integer SI categories = ", length(serial_cat))
message("  N l grid                = ", length(l_grid))
message("  N u grid                = ", length(u_grid))
message("  N gt grid               = ", length(gt_grid))
message("  N ip grid               = ", length(ip_grid))
message("  N delta grid            = ", length(delta_grid))
message("  N IP pairs              = ", length(pair_delta))
message("  N coarsening terms      = ", length(coarse_w))

# ------------------------------------------------------------
# Mpox high-confidence observed-pair data
# ------------------------------------------------------------

mpox_raw <- read_csv(MPOX_DATA_URL, show_col_types = FALSE)
mpox_score9 <- mpox_raw %>%
  filter(score == 9, Conf_pair == "high", Conf_sym == "high") %>%
  mutate(serial = as.integer(SI))

if (nrow(mpox_score9) != 34) {
  warning("Expected 34 score==9 high/high pairs, but found ", nrow(mpox_score9), ". Check upstream data.")
}
N_PAIRS <- nrow(mpox_score9)

serial_w_mpox <- mpox_score9 %>%
  count(serial, name = "n") %>%
  arrange(serial)

write_csv(mpox_score9, file.path(OUT_DIR, "mpox_score9_high_confidence_pairs.csv"))
write_csv(serial_w_mpox, file.path(OUT_DIR, "mpox_score9_weighted_serial_intervals.csv"))
write_csv(data.frame(serial = serial_cat), file.path(OUT_DIR, "mpox_serial_category_support.csv"))
write_csv(data.frame(gt = gt_grid), file.path(OUT_DIR, "mpox_generation_time_grid.csv"))

message("\nMpox high-confidence serial data:")
message("  unique SI values = ", nrow(serial_w_mpox), ", total pairs = ", sum(serial_w_mpox$n))
message("  SI range = ", min(serial_w_mpox$serial), " to ", max(serial_w_mpox$serial))

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

    prior_log_lambda1_mean = prior_log_lambda1_mean,
    prior_log_lambda1_sd = prior_log_lambda1_sd,
    prior_log_lambda_I_mean = prior_log_lambda_I_mean,
    prior_log_lambda_I_sd = prior_log_lambda_I_sd,
    prior_log_sigma_z_mean = prior_log_sigma_z_mean,
    prior_log_sigma_z_sd = prior_log_sigma_z_sd,

    prior_theta_mean = as.numeric(prior_theta_mean),
    prior_theta_cov = prior_theta_cov
  )
}

data_mpox <- make_stan_data(serial_w_mpox, r_growth)

# ------------------------------------------------------------
# Fit date-coarsened mpox SI model
# ------------------------------------------------------------

mod_mpox <- cmdstan_model(MPOX_STAN_FILE, force_recompile = TRUE)

make_inits <- function() {
  function() {
    list(
      log_lambda1 = prior_log_lambda1_mean,
      log_lambda_I = prior_log_lambda_I_mean,
      log_r = data_mpox$prior_theta_mean[1],
      log_d = data_mpox$prior_theta_mean[2],
      m = min(max(data_mpox$prior_theta_mean[3], 0.5), 35),
      log_kappa = data_mpox$prior_theta_mean[4],
      log_sigma_z = prior_log_sigma_z_mean
    )
  }
}

fit_slug <- paste0(
  "date_coarsened_Yang_MVN_",
  theta_marg$yang_prior_source[1],
  "_sdmult_", gsub("\\.", "p", theta_marg$yang_prior_sd_multiplier[1])
)

message("\nFitting mpox score==9 date-coarsened SI model...")
message("  prior = Yang MVN only")
message("  fit_slug = ", fit_slug)

fit_mpox <- mod_mpox$sample(
  data = data_mpox,
  chains = chains_mpox,
  parallel_chains = parallel_chains_mpox,
  iter_warmup = iter_warmup_mpox,
  iter_sampling = iter_sampling_mpox,
  seed = seed_mpox,
  init = make_inits(),
  adapt_delta = adapt_delta_mpox,
  max_treedepth = max_treedepth_mpox,
  refresh = refresh_mpox,
  save_warmup = FALSE
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
    prior = "Yang_MVN_only",
    yang_prior_source = theta_marg$yang_prior_source[1],
    yang_prior_sd_multiplier = theta_marg$yang_prior_sd_multiplier[1]
  ) %>%
  select(dataset, fit, prior, yang_prior_source, yang_prior_sd_multiplier,
         variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_mpox)
write_csv(summ_mpox, file.path(OUT_DIR, paste0("summary_mpox_score9_", fit_slug, ".csv")))

# Explicit figure-ready draws/objects.
gen_prob_mat <- as_draws_matrix(fit_mpox$draws(variables = "generation_prob_grid"))
pi_cat_mat <- as_draws_matrix(fit_mpox$draws(variables = "pi_cat_save"))
log_lik_mat <- extract_log_lik_matrix(fit_mpox, "log_lik")
log_lik_grouped_mat <- extract_log_lik_matrix(fit_mpox, "log_lik_grouped")

if (ncol(log_lik_mat) != N_PAIRS) {
  stop("Pointwise log_lik has ", ncol(log_lik_mat), " columns, expected ", N_PAIRS, ". LOO/WAIC not comparable.")
}

saveRDS(gen_prob_mat, file.path(OUT_DIR, paste0("generation_prob_grid_mpox_score9_", fit_slug, ".rds")))
saveRDS(pi_cat_mat, file.path(OUT_DIR, paste0("pi_cat_date_coarsened_mpox_score9_", fit_slug, ".rds")))
saveRDS(log_lik_mat, file.path(OUT_DIR, paste0("log_lik_mpox_score9_", fit_slug, ".rds")))
saveRDS(log_lik_grouped_mat, file.path(OUT_DIR, paste0("log_lik_grouped_mpox_score9_", fit_slug, ".rds")))

loo_fit <- loo(log_lik_mat)
waic_fit <- waic(log_lik_mat)
saveRDS(loo_fit, file.path(OUT_DIR, paste0("loo_mpox_score9_", fit_slug, ".rds")))
saveRDS(waic_fit, file.path(OUT_DIR, paste0("waic_mpox_score9_", fit_slug, ".rds")))

loo_waic_mpox <- tibble(
  dataset = "mpox_score9",
  fit = fit_slug,
  prior = "Yang_MVN_only",
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

diag_txt <- capture.output(fit_mpox$cmdstan_diagnose())
writeLines(diag_txt, file.path(OUT_DIR, paste0("cmdstan_diagnose_mpox_score9_", fit_slug, ".txt")))
cat(paste(diag_txt, collapse = "\n"))

draws_mpox <- as_draws_df(fit_mpox$draws(variables = pars_summary)) %>%
  mutate(
    dataset = "mpox_score9",
    fit = fit_slug,
    prior = "Yang_MVN_only"
  )
write_csv(as.data.frame(draws_mpox), file.path(OUT_DIR, paste0("draws_mpox_score9_", fit_slug, ".csv")))

message("\nDone. Outputs saved in: ", OUT_DIR)
