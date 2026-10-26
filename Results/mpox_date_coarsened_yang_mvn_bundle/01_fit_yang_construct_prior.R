# ============================================================
# 01_fit_yang_construct_prior.R
#
# Fits Yang et al. Figure 2 mpox viral-load data using the manuscript
# V(t), then constructs a Yang-informed infection-time MVN prior for
# theta = [log_r, log_d, m, log_kappa].
#
# This script intentionally constructs ONLY the Yang-informed MVN prior.
# It does not construct or use a robust mixture prior.
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(readxl)
  library(dplyr)
  library(readr)
})

this_file_dir <- getwd()
YANG_STAN_FILE <- file.path(this_file_dir, "yang_mpox_wh_viral_load_model.stan")
if (!file.exists(YANG_STAN_FILE)) YANG_STAN_FILE <- "yang_mpox_wh_viral_load_model.stan"
if (!file.exists(YANG_STAN_FILE)) YANG_STAN_FILE <- file.path("/mnt/data", "yang_mpox_wh_viral_load_model.stan")
if (!file.exists(YANG_STAN_FILE)) stop("Could not find yang_mpox_wh_viral_load_model.stan")

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

OUT_DIR <- "mpox_score9_date_coarsened_Yang_MVN_outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

# Options: "shared_all_sites", "skin_lesion", "rectal_swab", "saliva", "oropharyngeal_swab"
YANG_PRIOR_SOURCE <- "shared_all_sites"

yang_fit_prior_list <- list(
  prior_log_r_mean = log(0.35),
  prior_log_r_sd = 1.0,
  prior_log_d_mean = log(0.20),
  prior_log_d_sd = 1.0,
  prior_m_mean = 6.0,          # Yang scale: days post symptom onset
  prior_m_sd = 4.0,
  prior_log_kappa_mean = 0.0,
  prior_log_kappa_sd = 1.0
)

# Convert Yang m from post-symptom-onset to post-infection.
INCUBATION_MEAN_FOR_SHIFT <- 8.0
INCUBATION_SD_FOR_SHIFT <- 3.0

# Values >1 make the Yang prior weaker while preserving posterior covariance.
YANG_PRIOR_SD_MULTIPLIER <- 2.0

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
  theta_yang <- as.matrix(yang_draws[, c("log_r", "log_d", "m", "log_kappa")])
  colnames(theta_yang) <- c("log_r", "log_d", "m", "log_kappa")

  theta_infection <- theta_yang
  theta_infection[, "m"] <- theta_infection[, "m"] + INCUBATION_MEAN_FOR_SHIFT

  theta_mean <- colMeans(theta_infection)
  theta_cov <- cov(theta_infection)
  theta_cov <- (YANG_PRIOR_SD_MULTIPLIER^2) * theta_cov
  theta_cov["m", "m"] <- theta_cov["m", "m"] + INCUBATION_SD_FOR_SHIFT^2
  theta_cov <- make_spd(theta_cov)

  list(
    theta_mean = as.numeric(theta_mean),
    theta_cov = theta_cov,
    theta_corr_yang_post_onset = cor(theta_yang),
    theta_corr_infection_prior = cov2cor(theta_cov),
    marginal_prior = tibble(
      parameter = c("log_r", "log_d", "m", "log_kappa"),
      prior_mean = as.numeric(theta_mean),
      prior_sd = sqrt(diag(theta_cov)),
      prior_type = "Yang_MVN_only",
      yang_prior_source = YANG_PRIOR_SOURCE,
      yang_prior_sd_multiplier = YANG_PRIOR_SD_MULTIPLIER,
      incubation_mean_shift = INCUBATION_MEAN_FOR_SHIFT,
      incubation_sd_shift = INCUBATION_SD_FOR_SHIFT
    )
  )
}

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
write.csv(theta_prior$theta_corr_yang_post_onset, file.path(OUT_DIR, "yang_theta_posterior_correlation_post_onset.csv"), row.names = TRUE)
write.csv(theta_prior$theta_corr_infection_prior, file.path(OUT_DIR, "constructed_theta_prior_correlation_infection_time.csv"), row.names = TRUE)

message("\nConstructed Yang MVN theta prior marginals:")
print(theta_prior$marginal_prior)
message("\nDone. Outputs saved in: ", OUT_DIR)
