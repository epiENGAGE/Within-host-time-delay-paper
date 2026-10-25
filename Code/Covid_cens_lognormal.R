# SARS-CoV-2 censored lognormal comparator
#
# Purpose
# -------
# Fit the fourth model requested by the reviewer:
#   lognormal generation interval + day-level onset-date coarsening/censoring.
#
# The observation model is matched to Covid cens.R: reported symptom-onset dates
# are treated as day-level observations, giving a triangular kernel for the
# difference of two independent nearest-day recording errors.
#
# The epidemiological structure follows the Park-style lognormal comparator:
#   SI = G + IP2 - IP1,
# with a fixed incubation distribution, an epidemic-growth correction on the
# infector incubation period, and correlation rho between log(IP1) and log(G).
#
# IMPORTANT:
# The original Park et al. analysis is maximum-likelihood rather than Bayesian.
# The priors below are deliberately broad and common across variants. If your old
# uncensored Bayesian lognormal fit used different priors, replace only the four
# PRIOR_* controls below so censored and uncensored Bayesian comparators match.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(readxl)
  library(loo)
  library(tibble)
  library(readr)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
CODE_DIR <- path.expand("~/Within-host-time-delay-paper/Code")
DATA_DIR <- path.expand("~/Within-host-time-delay-paper/Data")
OUT_DIR  <- path.expand("~/Within-host-time-delay-paper/Results/Covid/Cens/lognormal")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STAN_FILE <- file.path(CODE_DIR, "covid_cens_lognormal.stan")
SERIAL_FILE <- file.path(DATA_DIR, "serial-netherlands.xlsx")

if (!file.exists(STAN_FILE)) {
  stop("Could not find Stan file: ", STAN_FILE,
       "\nCopy covid_cens_lognormal.stan into ~/Within-host-time-delay-paper/Code/")
}
if (!file.exists(SERIAL_FILE)) stop("Could not find serial data: ", SERIAL_FILE)

# -----------------------------------------------------------------------------
# Epidemiological settings
# -----------------------------------------------------------------------------
# Same growth corrections as your current SARS-CoV-2 analyses.
r_sgtf  <-  0.15
r_nsgtf <- -0.05

# Park et al. use rho = 0.75 as the main assumption for correlation between
# log incubation period of the infector and log generation interval.
RHO <- 0.75

# Backer et al. backward-incubation summaries for this Dutch dataset.
# These are used only to reconstruct the Park-style forward incubation
# distribution under the specified growth rate.
BACKWARD_INCUBATION <- tibble::tribble(
  ~strain,    ~mean_days, ~sd_days,
  "SGTF",           3.2,      2.2,
  "non-SGTF",       4.4,      2.5
)

# -----------------------------------------------------------------------------
# Numerical settings
# -----------------------------------------------------------------------------
# Reported SI categories exactly match the censored WH analysis.
SI_CAT_MIN <- -10L
SI_CAT_MAX <-  20L
serial_cat <- SI_CAT_MIN:SI_CAT_MAX

# Latent SI grid. We use 0.25 d to match the original lognormal-comparator
# discretisation and keep the convolution computationally manageable.
DX_SI <- 0.25
SI_LATENT_MIN <- SI_CAT_MIN - 1
SI_LATENT_MAX <- SI_CAT_MAX + 1
si_latent_grid <- seq(SI_LATENT_MIN, SI_LATENT_MAX, by = DX_SI)

# Incubation quadrature. 0.25 d matches the existing lognormal implementation.
DX_IP <- 0.25
IP_MIN <- 0.05
IP_MAX <- 30
ip_grid <- seq(IP_MIN, IP_MAX, by = DX_IP)

# -----------------------------------------------------------------------------
# Broad Bayesian priors for lognormal generation interval
# -----------------------------------------------------------------------------
# Common prior center: mean 4 d, SD 2.5 d. The prior SDs are broad on the
# lognormal parameter scale and are intentionally the same for both variants.
PRIOR_GEN_MEAN_DAYS <- 4.0
PRIOR_GEN_SD_DAYS   <- 2.5
PRIOR_LOGMEAN_GEN_SD <- 1.0
PRIOR_LOGSD_GEN_SD   <- 0.75

# -----------------------------------------------------------------------------
# MCMC settings
# -----------------------------------------------------------------------------
CHAINS <- 4
PARALLEL_CHAINS <- 4
ITER_WARMUP <- 500
ITER_SAMPLING <- 500
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12
SEED <- 3301
REFRESH <- 10

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
lnorm_params_from_mean_sd <- function(mean, sd) {
  if (!is.finite(mean) || !is.finite(sd) || mean <= 0 || sd <= 0) {
    stop("mean and sd must be positive finite values")
  }
  sigma2 <- log1p((sd / mean)^2)
  list(
    mu = log(mean) - 0.5 * sigma2,
    sigma = sqrt(sigma2)
  )
}

# Reconstruct the forward incubation distribution from the backward lognormal
# by exponential tilting: f_forward(t) proportional to exp(r*t) f_backward(t).
# Park et al. then approximate the corrected forward distribution by a lognormal
# matched on its mean and SD; we do the same here.
forward_incubation_from_backward <- function(mean_back, sd_back, r_growth,
                                             grid_max = 40, dx = 0.0025) {
  p <- lnorm_params_from_mean_sd(mean_back, sd_back)
  x <- seq(dx, grid_max, by = dx)
  logw <- dlnorm(x, meanlog = p$mu, sdlog = p$sigma, log = TRUE) + r_growth * x
  logw <- logw - max(logw)
  w <- exp(logw)
  z <- sum(w) * dx
  if (!is.finite(z) || z <= 0) stop("Failed to normalize forward incubation distribution")
  f <- w / z
  m1 <- sum(x * f) * dx
  m2 <- sum(x^2 * f) * dx
  s <- sqrt(max(m2 - m1^2, 0))
  approx_par <- lnorm_params_from_mean_sd(m1, s)
  list(mean = m1, sd = s, mu = approx_par$mu, sigma = approx_par$sigma)
}

make_weighted_serial <- function(df, strain_value, household_value = "within") {
  df %>%
    filter(.data$strain == strain_value, .data$household == household_value) %>%
    group_by(.data$serial) %>%
    summarize(n = sum(.data$n), .groups = "drop") %>%
    arrange(.data$serial)
}

make_counts_on_support <- function(serial_df, support) {
  out <- integer(length(support))
  idx <- match(as.integer(serial_df$serial), support)
  if (anyNA(idx)) {
    stop("Observed serial interval outside [", min(support), ", ", max(support), "]: ",
         paste(serial_df$serial[is.na(idx)], collapse = ", "))
  }
  out[idx] <- as.integer(serial_df$n)
  out
}

expand_obs_cat <- function(counts) {
  rep.int(seq_along(counts), counts)
}

# Sparse Park-style convolution terms for a single strain.
# const_lp contains the fixed parts of the quadrature:
#   -r * IP1 + log f(IP1) + log f(IP2) + 2 log(dx_ip)
# The conditional generation density is evaluated in Stan because its parameters
# are inferred.
make_quadrature_terms <- function(si_grid, ip_grid, dx_ip,
                                  logmean_inc, logsd_inc, r_growth) {
  s_id <- integer(0)
  const_lp <- numeric(0)
  log_inc1 <- numeric(0)
  log_gen <- numeric(0)
  
  log_f_ip <- dlnorm(ip_grid, meanlog = logmean_inc, sdlog = logsd_inc, log = TRUE)
  log_dx2 <- 2 * log(dx_ip)
  
  # This is intentionally explicit and transparent. At the current grid sizes it
  # creates roughly 1--2 million sparse terms per strain; construction takes some
  # time once, while Stan then reuses them for every posterior evaluation.
  for (ks in seq_along(si_grid)) {
    s <- si_grid[ks]
    
    for (i1 in seq_along(ip_grid)) {
      gen_vec <- s + ip_grid[i1] - ip_grid
      ok <- which(is.finite(gen_vec) & gen_vec > 0)
      if (length(ok) == 0L) next
      
      n_ok <- length(ok)
      s_id <- c(s_id, rep.int(ks, n_ok))
      const_lp <- c(
        const_lp,
        -r_growth * ip_grid[i1] + log_f_ip[i1] + log_f_ip[ok] + log_dx2
      )
      log_inc1 <- c(log_inc1, rep.int(log(ip_grid[i1]), n_ok))
      log_gen <- c(log_gen, log(gen_vec[ok]))
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

# Same nearest-day coarsening kernel used conceptually in Covid cens.R.
# If latent SI=s and each onset date is rounded to the nearest day independently,
# the difference in rounding errors is triangular on [-1,1].
make_coarsening_matrix <- function(serial_cat, si_latent_grid) {
  outer(
    serial_cat,
    si_latent_grid,
    FUN = function(y, s) pmax(0, 1 - abs(s - y))
  )
}

prior_center <- lnorm_params_from_mean_sd(PRIOR_GEN_MEAN_DAYS, PRIOR_GEN_SD_DAYS)
coarse_w <- make_coarsening_matrix(serial_cat, si_latent_grid)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
serialdata <- read_xlsx(SERIAL_FILE)
serial_w_sgtf  <- make_weighted_serial(serialdata, "SGTF", "within")
serial_w_nsgtf <- make_weighted_serial(serialdata, "non-SGTF", "within")

make_stan_data <- function(serial_df, strain, r_growth) {
  inc_back <- BACKWARD_INCUBATION %>% filter(.data$strain == .env$strain)
  if (nrow(inc_back) != 1L) stop("Incubation assumptions missing for ", strain)
  
  inc_fwd <- forward_incubation_from_backward(
    mean_back = inc_back$mean_days,
    sd_back = inc_back$sd_days,
    r_growth = r_growth
  )
  
  message(
    strain, " forward incubation approximation: mean=",
    sprintf("%.3f", inc_fwd$mean), " d, SD=", sprintf("%.3f", inc_fwd$sd), " d"
  )
  
  terms <- make_quadrature_terms(
    si_grid = si_latent_grid,
    ip_grid = ip_grid,
    dx_ip = DX_IP,
    logmean_inc = inc_fwd$mu,
    logsd_inc = inc_fwd$sigma,
    r_growth = r_growth
  )
  
  counts <- make_counts_on_support(serial_df, serial_cat)
  obs_cat_pt <- expand_obs_cat(counts)
  
  c(
    list(
      N_cat = length(serial_cat),
      counts = as.integer(counts),
      serial_cat = as.numeric(serial_cat),
      N_obs = length(obs_cat_pt),
      obs_cat_pt = as.integer(obs_cat_pt),
      N_si_latent = length(si_latent_grid),
      si_latent_grid = as.numeric(si_latent_grid),
      dx_si = DX_SI,
      coarse_w = coarse_w,
      logmean_inc = inc_fwd$mu,
      logsd_inc = inc_fwd$sigma,
      rho = RHO,
      prior_logmean_gen_mean = prior_center$mu,
      prior_logmean_gen_sd = PRIOR_LOGMEAN_GEN_SD,
      prior_logsd_gen_mean = prior_center$sigma,
      prior_logsd_gen_sd = PRIOR_LOGSD_GEN_SD
    ),
    terms
  )
}

message("Building SGTF / Omicron quadrature...")
data_sgtf <- make_stan_data(serial_w_sgtf, "SGTF", r_sgtf)
message("Building non-SGTF / Delta quadrature...")
data_nsgtf <- make_stan_data(serial_w_nsgtf, "non-SGTF", r_nsgtf)

# Fail early if the R/Stan interface is incomplete.
REQUIRED_STAN_DATA <- c(
  "N_cat", "counts", "serial_cat", "N_obs", "obs_cat_pt",
  "N_si_latent", "si_latent_grid", "dx_si",
  "N_terms", "s_id", "const_lp", "log_inc1", "log_gen",
  "coarse_w", "logmean_inc", "logsd_inc", "rho",
  "prior_logmean_gen_mean", "prior_logmean_gen_sd",
  "prior_logsd_gen_mean", "prior_logsd_gen_sd"
)
check_stan_data <- function(x, label) {
  miss <- setdiff(REQUIRED_STAN_DATA, names(x))
  if (length(miss)) stop(label, " is missing Stan data: ", paste(miss, collapse = ", "))
  if (!is.matrix(x$coarse_w)) stop(label, ": coarse_w must be a matrix")
  if (!identical(dim(x$coarse_w), c(x$N_cat, x$N_si_latent))) {
    stop(label, ": coarse_w has dimensions ", paste(dim(x$coarse_w), collapse = " x "),
         "; expected ", x$N_cat, " x ", x$N_si_latent)
  }
  invisible(TRUE)
}
check_stan_data(data_sgtf, "SGTF")
check_stan_data(data_nsgtf, "non-SGTF")
message("Stan-data checks passed for both variants.")

# -----------------------------------------------------------------------------
# Compile and fit
# -----------------------------------------------------------------------------
mod <- cmdstan_model(STAN_FILE, force_recompile = TRUE)

make_inits <- function(mean_days, sd_days) {
  p <- lnorm_params_from_mean_sd(mean_days, sd_days)
  function() list(logmean_gen = p$mu, logsd_gen = p$sigma)
}

fit_one <- function(data, this_seed, label, init_mean, init_sd) {
  message("\nFitting censored lognormal: ", label, "...")
  mod$sample(
    data = data,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    seed = this_seed,
    init = make_inits(init_mean, init_sd),
    adapt_delta = ADAPT_DELTA,
    max_treedepth = MAX_TREEDEPTH,
    refresh = REFRESH,
    save_warmup = FALSE
  )
}

# Initial values close to the published Park within-household estimates; these are
# only initial states, not priors.
fit_sgtf <- fit_one(data_sgtf, SEED,     "SGTF / Omicron",    3.0, 2.0)
fit_nsgtf <- fit_one(data_nsgtf, SEED+1, "non-SGTF / Delta", 3.8, 2.0)

# -----------------------------------------------------------------------------
# Save fits and posterior summaries
# -----------------------------------------------------------------------------
fit_sgtf$save_object(file.path(OUT_DIR, "fit_censored_lognormal_sgtf.rds"))
fit_nsgtf$save_object(file.path(OUT_DIR, "fit_censored_lognormal_nsgtf.rds"))

pars_summary <- c("logmean_gen", "logsd_gen", "mean_gen", "sd_gen")

summ_all <- bind_rows(
  fit_sgtf$summary(variables = pars_summary) %>% mutate(strain = "SGTF", .before = 1),
  fit_nsgtf$summary(variables = pars_summary) %>% mutate(strain = "non-SGTF", .before = 1)
) %>%
  select(strain, variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

print(summ_all)
write_csv(summ_all, file.path(OUT_DIR, "censored_lognormal_summary.csv"))

# -----------------------------------------------------------------------------
# LOO / WAIC
# -----------------------------------------------------------------------------
extract_log_lik_matrix <- function(fit, variable = "log_lik") {
  posterior::as_draws_matrix(fit$draws(variables = variable, inc_warmup = FALSE))
}

compute_ic <- function(fit, strain_slug, strain_label) {
  ll <- extract_log_lik_matrix(fit, "log_lik")
  ll_grouped <- extract_log_lik_matrix(fit, "log_lik_grouped")
  
  saveRDS(ll, file.path(OUT_DIR, paste0(strain_slug, "_log_lik_matrix.rds")))
  saveRDS(ll_grouped, file.path(OUT_DIR, paste0(strain_slug, "_log_lik_grouped_matrix.rds")))
  
  loo_fit <- loo::loo(ll)
  waic_fit <- loo::waic(ll)
  
  saveRDS(loo_fit, file.path(OUT_DIR, paste0(strain_slug, "_loo.rds")))
  saveRDS(waic_fit, file.path(OUT_DIR, paste0(strain_slug, "_waic.rds")))
  
  tibble(
    strain = strain_label,
    model = "Lognormal censored",
    n_pointwise_obs = ncol(ll),
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_fit$estimates["elpd_loo", "SE"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    se_looic = loo_fit$estimates["looic", "SE"],
    elpd_waic = waic_fit$estimates["elpd_waic", "Estimate"],
    se_elpd_waic = waic_fit$estimates["elpd_waic", "SE"],
    p_waic = waic_fit$estimates["p_waic", "Estimate"],
    waic = waic_fit$estimates["waic", "Estimate"],
    se_waic = waic_fit$estimates["waic", "SE"]
  )
}

ic_all <- bind_rows(
  compute_ic(fit_sgtf, "sgtf", "SGTF"),
  compute_ic(fit_nsgtf, "nsgtf", "non-SGTF")
)

print(ic_all)
write_csv(ic_all, file.path(OUT_DIR, "censored_lognormal_loo_waic.csv"))

# -----------------------------------------------------------------------------
# Save posterior fitted reported-SI probabilities for plotting
# -----------------------------------------------------------------------------
extract_prob_summary <- function(fit, strain) {
  d <- posterior::as_draws_matrix(
    fit$draws(variables = "serial_prob_category", inc_warmup = FALSE)
  )
  
  tibble(
    strain = strain,
    serial = serial_cat,
    prob_median = apply(d, 2, median),
    prob_q025 = apply(d, 2, quantile, probs = 0.025),
    prob_q975 = apply(d, 2, quantile, probs = 0.975)
  )
}

prob_summary <- bind_rows(
  extract_prob_summary(fit_sgtf, "SGTF"),
  extract_prob_summary(fit_nsgtf, "non-SGTF")
)
write_csv(prob_summary, file.path(OUT_DIR, "censored_lognormal_serial_probabilities.csv"))

message("\nDONE. Outputs written to: ", OUT_DIR)
message("Next run the four-model comparison script.")