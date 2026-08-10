# Covid R0 plot
CENS_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

CENS_SGTF_RDS <- file.path(
  CENS_DIR,
  "fit_date_coarsened_joint_obspair_sgtf.rds"
)

CENS_NSGTF_RDS <- file.path(
  CENS_DIR,
  "fit_date_coarsened_joint_obspair_nsgtf.rds"
)

file.exists(CENS_SGTF_RDS)
file.exists(CENS_NSGTF_RDS)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(vroom)
  library(zoo)
  library(mgcv)
  library(mvtnorm)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
})

ROOT <- "~/Within-host-time-delay-framework/Data"

OLD_OUT_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Noncens"

NEW_OUT_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

FIG_DIR <- Sys.getenv(
  "FIG_DIR",
  unset = file.path(NEW_OUT_DIR, "figure5_three_model_rt_advantage")
)

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

LOGNORMAL_SGTF_RDS <- Sys.getenv(
  "LOGNORMAL_SGTF_RDS",
  unset = file.path(OLD_OUT_DIR, "fit_park_lognormal_fast_sgtf.rds")
)

LOGNORMAL_NSGTF_RDS <- Sys.getenv(
  "LOGNORMAL_NSGTF_RDS",
  unset = file.path(OLD_OUT_DIR, "fit_park_lognormal_fast_nsgtf.rds")
)

MECH_SGTF_RDS <- Sys.getenv(
  "MECH_SGTF_RDS",
  unset = file.path(OLD_OUT_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds")
)

MECH_NSGTF_RDS <- Sys.getenv(
  "MECH_NSGTF_RDS",
  unset = file.path(OLD_OUT_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds")
)

CASE_FILE <- Sys.getenv(
  "CASE_FILE",
  unset = file.path(ROOT, "COVID-19_aantallen_gemeente_per_dag.csv")
)

VARIANT_FILE <- Sys.getenv(
  "VARIANT_FILE",
  unset = file.path(ROOT, "variant-netherlands.xlsx")
)

OUT_PREFIX <- file.path(FIG_DIR, "figure5_three_model_rt_advantage")

NSAMPLE_TARGET <- as.integer(Sys.getenv("NSAMPLE_TARGET", unset = "400"))

MAX_GT_DAYS <- as.numeric(Sys.getenv("MAX_GT_DAYS", unset = "14"))
GT_BIN_BY <- as.numeric(Sys.getenv("GT_BIN_BY", unset = "0.5"))

SHOW_RIBBONS <- as.logical(Sys.getenv("SHOW_RIBBONS", unset = "TRUE"))

set.seed(101)

# Find censored fits
find_censored_fit <- function(strain_slug) {
  env_name <- if (strain_slug == "sgtf") "CENS_SGTF_FIT" else "CENS_NSGTF_FIT"
  env_path <- Sys.getenv(env_name, unset = "")
  
  if (nzchar(env_path) && file.exists(env_path)) {
    return(env_path)
  }
  
  hits <- list.files(
    NEW_OUT_DIR,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  hits <- hits[
    !grepl("_loo\\.rds$|_waic\\.rds$|loo\\.rds$|waic\\.rds$", basename(hits), ignore.case = TRUE)
  ]
  
  hits <- hits[
    grepl("fit|cens|censored", basename(hits), ignore.case = TRUE)
  ]
  
  if (strain_slug == "sgtf") {
    hits <- hits[
      grepl("sgtf", basename(hits), ignore.case = TRUE) &
        !grepl("nsgtf|non[-_ ]?sgtf", basename(hits), ignore.case = TRUE)
    ]
  } else {
    hits <- hits[
      grepl("nsgtf|non[-_ ]?sgtf", basename(hits), ignore.case = TRUE)
    ]
  }
  
  if (length(hits) != 1) {
    message("\nCandidate censored fit files for ", strain_slug, ":")
    print(hits)
    stop(
      "\nCould not uniquely identify censored ", strain_slug, " fit file.\n",
      "Set manually before running:\n",
      "Sys.setenv(CENS_SGTF_FIT='/full/path/to/censored_sgtf_fit.rds')\n",
      "Sys.setenv(CENS_NSGTF_FIT='/full/path/to/censored_nsgtf_fit.rds')\n"
    )
  }
  
  hits[1]
}

CENS_SGTF_RDS <- Sys.getenv(
  "CENS_SGTF_RDS",
  unset = find_censored_fit("sgtf")
)

CENS_NSGTF_RDS <- Sys.getenv(
  "CENS_NSGTF_RDS",
  unset = find_censored_fit("nsgtf")
)

# File checks
needed_files <- c(
  LOGNORMAL_SGTF_RDS,
  LOGNORMAL_NSGTF_RDS,
  MECH_SGTF_RDS,
  MECH_NSGTF_RDS,
  CENS_SGTF_RDS,
  CENS_NSGTF_RDS,
  CASE_FILE,
  VARIANT_FILE
)

missing_files <- needed_files[!file.exists(needed_files)]

# Utility functions
normalize_prob_matrix <- function(x, eps = 1e-12) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  
  rs <- rowSums(x)
  
  bad <- !is.finite(rs) | rs <= 0
  if (any(bad)) {
    x[bad, ] <- eps
    rs <- rowSums(x)
  }
  
  x / rs
}

subset_rows_evenly <- function(x, n) {
  x <- as.matrix(x)
  
  if (nrow(x) <= n) return(x)
  
  idx <- unique(round(seq(1, nrow(x), length.out = n)))
  x[idx, , drop = FALSE]
}

fit_variables <- function(fit_obj) {
  tryCatch(
    posterior::variables(fit_obj$draws(inc_warmup = FALSE)),
    error = function(e) {
      tryCatch(
        fit_obj$metadata()$model_params,
        error = function(e2) character()
      )
    }
  )
}

variable_regex <- function(variable) {
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", variable)
  paste0("^", escaped, "(\\[|$)")
}

has_base_variable <- function(fit_obj, variable) {
  vars <- fit_variables(fit_obj)
  any(stringr::str_detect(vars, variable_regex(variable)))
}

draws_matrix_one <- function(fit_obj, variable) {
  if (!has_base_variable(fit_obj, variable)) {
    return(NULL)
  }
  
  tryCatch(
    posterior::as_draws_matrix(
      fit_obj$draws(variables = variable, inc_warmup = FALSE)
    ),
    error = function(e) NULL
  )
}

draws_df_vars <- function(fit_obj, variables) {
  have <- variables[
    vapply(variables, function(v) has_base_variable(fit_obj, v), logical(1))
  ]
  
  if (length(have) == 0) return(tibble())
  
  posterior::as_draws_df(
    fit_obj$draws(variables = have, inc_warmup = FALSE)
  ) %>%
    as_tibble()
}

infer_gt_grid <- function(n_cols, model_name) {
  candidates <- list(
    seq(0, 20, by = 0.1),
    seq(0, 20, by = 0.2),
    seq(0, 20, by = 0.25),
    seq(0, 20, by = 0.5),
    seq(0, 30, by = 0.1),
    seq(0, 30, by = 0.2),
    seq(0, 30, by = 0.25),
    seq(0, 30, by = 0.5),
    seq(0, 40, by = 0.1),
    seq(0, 40, by = 0.2),
    seq(0, 40, by = 0.25),
    seq(0, 40, by = 0.5),
    seq(-10, 20, by = 0.1),
    seq(-10, 20, by = 0.2),
    seq(-10, 20, by = 0.25),
    seq(-10, 20, by = 0.5)
  )
  
  lens <- vapply(candidates, length, integer(1))
  idx <- which(lens == n_cols)
  
  if (length(idx) == 1) {
    return(candidates[[idx]])
  }
  
  stop(
    "Could not infer generation-time grid for ",
    model_name,
    ". Matrix has ",
    n_cols,
    " columns. Add the correct grid to infer_gt_grid()."
  )
}

density_grid_to_bin_mass <- function(density_mat, grid, breaks) {
  density_mat <- as.matrix(density_mat)
  dx <- median(diff(grid))
  
  out <- matrix(
    0,
    nrow = nrow(density_mat),
    ncol = length(breaks) - 1
  )
  
  for (j in seq_len(length(breaks) - 1)) {
    keep <- grid >= breaks[j] & grid < breaks[j + 1]
    
    if (any(keep)) {
      out[, j] <- rowSums(density_mat[, keep, drop = FALSE]) * dx
    }
  }
  
  normalize_prob_matrix(out)
}

lognormal_params_from_mean_sd <- function(mean, sd) {
  mean <- pmax(mean, 1e-8)
  sd <- pmax(sd, 1e-8)
  
  sdlog <- sqrt(log(1 + (sd / mean)^2))
  meanlog <- log(mean) - 0.5 * sdlog^2
  
  tibble(
    meanlog = meanlog,
    sdlog = sdlog
  )
}

lognormal_bin_mass_from_params <- function(meanlog, sdlog, breaks) {
  out <- t(vapply(seq_along(meanlog), function(i) {
    p <- diff(plnorm(
      breaks,
      meanlog = meanlog[i],
      sdlog = sdlog[i]
    ))
    
    if (!is.finite(sum(p)) || sum(p) <= 0) {
      p <- rep(1 / (length(breaks) - 1), length(breaks) - 1)
    } else {
      p <- p / sum(p)
    }
    
    p
  }, numeric(length(breaks) - 1)))
  
  out
}

extract_generation_mass_cmdstan <- function(fit_obj,
                                            model_name,
                                            strain_name,
                                            nsample,
                                            max_days = 14,
                                            by = 0.5) {
  breaks <- seq(0, max_days, by = by)
  
  density_vars <- c(
    "generation_density_norm_grid",
    "gen_density_norm_grid",
    "gt_density_norm_grid",
    "generation_density_grid",
    "gen_density_grid",
    "gt_density_grid",
    "generation_time_density_norm_grid",
    "generation_interval_density_norm_grid"
  )
  
  prob_vars <- c(
    "generation_prob_grid",
    "gen_prob_grid",
    "gt_prob_grid",
    "p_gen",
    "p_gt",
    "generation_interval_prob"
  )
  
  log_unnorm_vars <- c(
    "log_unnorm_gen",
    "log_unnorm_gt",
    "log_unnorm_generation",
    "log_generation_unnorm",
    "log_gt_unnorm"
  )
  
  log_norm_vars <- c(
    "log_norm_gen",
    "log_norm_gt",
    "log_norm_generation"
  )
  
  for (v in density_vars) {
    mat <- draws_matrix_one(fit_obj, v)
    
    if (!is.null(mat)) {
      mat <- subset_rows_evenly(mat, nsample)
      mat[!is.finite(mat) | mat < 0] <- 0
      
      grid <- infer_gt_grid(ncol(mat), paste(model_name, strain_name, v))
      
      message(
        "Using generation density variable `", v, "` for ",
        model_name, " / ", strain_name,
        " with ", nrow(mat), " draws and ", ncol(mat), " grid points."
      )
      
      return(density_grid_to_bin_mass(mat, grid, breaks))
    }
  }
  
  for (v in prob_vars) {
    mat <- draws_matrix_one(fit_obj, v)
    
    if (!is.null(mat)) {
      mat <- subset_rows_evenly(mat, nsample)
      
      if (ncol(mat) == length(breaks) - 1) {
        message(
          "Using generation probability variable `", v, "` directly for ",
          model_name, " / ", strain_name,
          " with ", nrow(mat), " draws and ", ncol(mat), " bins."
        )
        
        return(normalize_prob_matrix(mat))
      }
      
      grid <- infer_gt_grid(ncol(mat), paste(model_name, strain_name, v))
      
      message(
        "Using generation probability-grid variable `", v, "` for ",
        model_name, " / ", strain_name,
        " with ", nrow(mat), " draws and ", ncol(mat), " grid points."
      )
      
      mat <- normalize_prob_matrix(mat)
      return(density_grid_to_bin_mass(mat, grid, breaks))
    }
  }
  
  for (lv in log_unnorm_vars) {
    logu <- draws_matrix_one(fit_obj, lv)
    
    if (is.null(logu)) next
    
    for (nv in log_norm_vars) {
      logn <- draws_matrix_one(fit_obj, nv)
      
      if (is.null(logn)) next
      
      logu <- subset_rows_evenly(logu, nsample)
      logn <- subset_rows_evenly(logn, nsample)
      
      if (nrow(logu) != nrow(logn)) {
        next
      }
      
      dens <- exp(sweep(logu, 1, as.numeric(logn[, 1]), FUN = "-"))
      dens[!is.finite(dens) | dens < 0] <- 0
      
      grid <- infer_gt_grid(ncol(dens), paste(model_name, strain_name, lv, nv))
      
      message(
        "Using `", lv, " - ", nv, "` for ",
        model_name, " / ", strain_name,
        " with ", nrow(dens), " draws and ", ncol(dens), " grid points."
      )
      
      return(density_grid_to_bin_mass(dens, grid, breaks))
    }
  }
  
  param_df <- draws_df_vars(
    fit_obj,
    c(
      "logmean_gen",
      "logsd_gen",
      "meanlog_gen",
      "sdlog_gen",
      "mean_gen",
      "sd_gen"
    )
  )
  
  if (nrow(param_df) > 0) {
    param_df <- as.data.frame(param_df)
    param_df <- param_df[
      round(seq(1, nrow(param_df), length.out = min(nrow(param_df), nsample))),
      ,
      drop = FALSE
    ]
    
    if (all(c("logmean_gen", "logsd_gen") %in% names(param_df))) {
      message(
        "Using lognormal generation parameters `logmean_gen`, `logsd_gen` for ",
        model_name, " / ", strain_name, "."
      )
      
      return(lognormal_bin_mass_from_params(
        meanlog = param_df$logmean_gen,
        sdlog = param_df$logsd_gen,
        breaks = breaks
      ))
    }
    
    if (all(c("meanlog_gen", "sdlog_gen") %in% names(param_df))) {
      message(
        "Using lognormal generation parameters `meanlog_gen`, `sdlog_gen` for ",
        model_name, " / ", strain_name, "."
      )
      
      return(lognormal_bin_mass_from_params(
        meanlog = param_df$meanlog_gen,
        sdlog = param_df$sdlog_gen,
        breaks = breaks
      ))
    }
    
    if (all(c("mean_gen", "sd_gen") %in% names(param_df))) {
      warning(
        "No explicit generation-density grid found for ",
        model_name, " / ", strain_name,
        ". Falling back to lognormal approximation from mean_gen/sd_gen."
      )
      
      pars <- lognormal_params_from_mean_sd(
        mean = param_df$mean_gen,
        sd = param_df$sd_gen
      )
      
      return(lognormal_bin_mass_from_params(
        meanlog = pars$meanlog,
        sdlog = pars$sdlog,
        breaks = breaks
      ))
    }
  }
  
  vars <- fit_variables(fit_obj)
  
}

# Build case trajectories
build_case_trajectory_objects <- function(nsample) {
  message("Building Park-style case trajectories from local data files...")
  
  datevec <- seq(as.Date("2021-11-28"), as.Date("2022-01-30"), by = 1)
  weekbreak <- as.Date("2021-12-19") + 7 * (-4:6)
  
  cases <- vroom::vroom(CASE_FILE, show_col_types = FALSE)
  variant <- readxl::read_xlsx(VARIANT_FILE)
  
  cases_all <- cases %>%
    dplyr::group_by(Date_of_publication) %>%
    dplyr::summarize(cases = sum(Total_reported), .groups = "drop") %>%
    dplyr::arrange(Date_of_publication) %>%
    dplyr::mutate(rollmean = zoo::rollmean(cases, 7, fill = NA)) %>%
    dplyr::filter(
      Date_of_publication >= "2021-11-22",
      Date_of_publication <= "2022-01-30"
    )
  
  cases_all_weekly <- cases_all %>%
    dplyr::mutate(group = cut(Date_of_publication, weekbreak + 1)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      cases = sum(cases),
      date = max(as.Date(Date_of_publication)),
      .groups = "drop"
    )
  
  variant2 <- data.frame(
    omicron = unlist(variant[9, -c(1:2, 13)]),
    delta   = unlist(variant[5, -c(1:2, 13)]),
    total   = unlist(variant[1, -c(1:2, 13)]),
    week    = c(4, 3, 2, 1, 52, 51, 50, 49, 48, 47),
    year    = c(2022, 2022, 2022, 2022, 2021, 2021, 2021, 2021, 2021, 2021),
    date    = rev(cases_all_weekly$date)
  )
  
  variant3 <- variant2 %>%
    dplyr::group_by(year, week) %>%
    dplyr::mutate(
      omicron_prop = omicron / total,
      omicron_prop_lwr = binom.test(omicron, total)[[4]][1],
      omicron_prop_upr = binom.test(omicron, total)[[4]][2],
      delta_prop = delta / total,
      delta_prop_lwr = binom.test(delta, total)[[4]][1],
      delta_prop_upr = binom.test(delta, total)[[4]][2]
    ) %>%
    merge(cases_all_weekly) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(time = 1:dplyr::n())
  
  variant4 <- variant3 %>%
    dplyr::mutate(
      omicron_cases = omicron_prop * cases,
      delta_cases = delta_prop * cases
    )
  
  gfit1 <- mgcv::gam(
    log(delta_cases) ~ s(time, bs = "cs"),
    data = variant4,
    method = "REML"
  )
  
  gfit2 <- mgcv::gam(
    log(omicron_cases) ~ s(time, bs = "cs"),
    data = variant4,
    method = "REML"
  )
  
  time_grid <- seq(1, 10, by = 1 / 14)
  
  gfit1_p1 <- predict(
    gfit1,
    newdata = data.frame(time = time_grid),
    type = "lpmatrix"
  )
  
  gfit2_p1 <- predict(
    gfit2,
    newdata = data.frame(time = time_grid),
    type = "lpmatrix"
  )
  
  gfit1_sim <- mvtnorm::rmvnorm(
    nsample,
    mean = coef(gfit1),
    sigma = vcov(gfit1)
  )
  
  gfit2_sim <- mvtnorm::rmvnorm(
    nsample,
    mean = coef(gfit2),
    sigma = vcov(gfit2)
  )
  
  list(
    datevec = datevec,
    weekbreak = weekbreak,
    time_grid = time_grid,
    cases = cases,
    variant = variant,
    cases_all = cases_all,
    cases_all_weekly = cases_all_weekly,
    variant2 = variant2,
    variant3 = variant3,
    variant4 = variant4,
    gfit1 = gfit1,
    gfit2 = gfit2,
    gfit1_p1 = gfit1_p1,
    gfit2_p1 = gfit2_p1,
    gfit1_sim = gfit1_sim,
    gfit2_sim = gfit2_sim
  )
}

# Renewal calculation
compute_rt_draws_from_gen <- function(
    gen_delta,
    gen_omicron,
    i1_mat,
    i2_mat,
    time_grid,
    model_name
) {
  n_lag <- ncol(gen_delta)
  
  out <- lapply(seq_len(nrow(gen_delta)), function(x) {
    i1 <- c(i1_mat[, x])
    i2 <- c(i2_mat[, x])
    
    g1 <- gen_delta[x, ]
    g2 <- gen_omicron[x, ]
    
    R_delta <- tail(i1, -n_lag) /
      sapply(seq_len(length(i1) - n_lag), function(y) {
        sum(i1[y:(y + n_lag - 1)] * rev(g1))
      })
    
    R_omicron <- tail(i2, -n_lag) /
      sapply(seq_len(length(i2) - n_lag), function(y) {
        sum(i2[y:(y + n_lag - 1)] * rev(g2))
      })
    
    data.frame(
      model = model_name,
      sim = x,
      time = tail(time_grid, -n_lag),
      Delta = R_delta,
      Omicron = R_omicron,
      Advantage = R_omicron / R_delta
    )
  }) %>%
    dplyr::bind_rows()
  
  out
}

summarize_rt <- function(draw_df) {
  draw_df %>%
    tidyr::pivot_longer(
      cols = c("Delta", "Omicron", "Advantage"),
      names_to = "quantity",
      values_to = "value"
    ) %>%
    dplyr::group_by(model, quantity, time) %>%
    dplyr::summarize(
      median = median(value, na.rm = TRUE),
      lwr = quantile(value, 0.025, na.rm = TRUE),
      upr = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# Load models
logfit_sgtf <- readRDS(LOGNORMAL_SGTF_RDS)

logfit_nsgtf <- readRDS(LOGNORMAL_NSGTF_RDS)

mechfit_sgtf <- readRDS(MECH_SGTF_RDS)

mechfit_nsgtf <- readRDS(MECH_NSGTF_RDS)

censfit_sgtf <- readRDS(CENS_SGTF_RDS)

censfit_nsgtf <- readRDS(CENS_NSGTF_RDS)

n_draws_available <- c(
  nrow(posterior::as_draws_matrix(logfit_sgtf$draws(inc_warmup = FALSE))),
  nrow(posterior::as_draws_matrix(logfit_nsgtf$draws(inc_warmup = FALSE))),
  nrow(posterior::as_draws_matrix(mechfit_sgtf$draws(inc_warmup = FALSE))),
  nrow(posterior::as_draws_matrix(mechfit_nsgtf$draws(inc_warmup = FALSE))),
  nrow(posterior::as_draws_matrix(censfit_sgtf$draws(inc_warmup = FALSE))),
  nrow(posterior::as_draws_matrix(censfit_nsgtf$draws(inc_warmup = FALSE)))
)

ns <- min(NSAMPLE_TARGET, n_draws_available)

# Build case trajectories
case_objs <- build_case_trajectory_objects(nsample = ns)

i1_mat <- exp(case_objs$gfit1_p1 %*% t(case_objs$gfit1_sim)) / 7
i2_mat <- exp(case_objs$gfit2_p1 %*% t(case_objs$gfit2_sim)) / 7

# Build generation-interval masses
gen_delta_log <- extract_generation_mass_cmdstan(
  fit_obj = logfit_nsgtf,
  model_name = "Lognormal",
  strain_name = "Delta/non-SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

gen_omicron_log <- extract_generation_mass_cmdstan(
  fit_obj = logfit_sgtf,
  model_name = "Lognormal",
  strain_name = "Omicron/SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

gen_delta_mech <- extract_generation_mass_cmdstan(
  fit_obj = mechfit_nsgtf,
  model_name = "WH-informed",
  strain_name = "Delta/non-SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

gen_omicron_mech <- extract_generation_mass_cmdstan(
  fit_obj = mechfit_sgtf,
  model_name = "WH-informed",
  strain_name = "Omicron/SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

gen_delta_cens <- extract_generation_mass_cmdstan(
  fit_obj = censfit_nsgtf,
  model_name = "Censored WH-informed",
  strain_name = "Delta/non-SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

gen_omicron_cens <- extract_generation_mass_cmdstan(
  fit_obj = censfit_sgtf,
  model_name = "Censored WH-informed",
  strain_name = "Omicron/SGTF",
  nsample = ns,
  max_days = MAX_GT_DAYS,
  by = GT_BIN_BY
)

# Compute Rt and advantage
draws_log <- compute_rt_draws_from_gen(
  gen_delta = gen_delta_log,
  gen_omicron = gen_omicron_log,
  i1_mat = i1_mat,
  i2_mat = i2_mat,
  time_grid = case_objs$time_grid,
  model_name = "Lognormal"
)

draws_mech <- compute_rt_draws_from_gen(
  gen_delta = gen_delta_mech,
  gen_omicron = gen_omicron_mech,
  i1_mat = i1_mat,
  i2_mat = i2_mat,
  time_grid = case_objs$time_grid,
  model_name = "WH-informed"
)

draws_cens <- compute_rt_draws_from_gen(
  gen_delta = gen_delta_cens,
  gen_omicron = gen_omicron_cens,
  i1_mat = i1_mat,
  i2_mat = i2_mat,
  time_grid = case_objs$time_grid,
  model_name = "Censored WH-informed"
)

all_draws <- dplyr::bind_rows(draws_mech, draws_cens, draws_log)
all_summary <- summarize_rt(all_draws)

readr::write_csv(all_draws, paste0(OUT_PREFIX, "_draws.csv"))
readr::write_csv(all_summary, paste0(OUT_PREFIX, "_summary.csv"))

# Plot
plot_df <- all_summary %>%
  dplyr::mutate(
    quantity = factor(
      quantity,
      levels = c("Delta", "Omicron", "Advantage")
    ),
    model = factor(
      model,
      levels = c(
        "WH-informed",
        "Censored WH-informed",
        "Lognormal"
      )
    )
  )

date_breaks <- 2:10
date_labels <- c(
  "Nov 28", "Dec 5", "Dec 12", "Dec 19", "Dec 26",
  "Jan 2", "Jan 9", "Jan 16", "Jan 23"
)

model_linetypes <- c(
  "WH-informed" = "solid",
  "Censored WH-informed" = "longdash",
  "Lognormal" = "dotted"
)

quantity_cols <- c(
  Delta = "black",
  Omicron = "orange",
  Advantage = "purple"
)

p <- ggplot(
  plot_df,
  aes(
    x = time,
    y = median,
    colour = quantity,
    linetype = model,
    group = interaction(quantity, model)
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey45"
  )

if (SHOW_RIBBONS) {
  p <- p +
    geom_ribbon(
      data = plot_df %>%
        dplyr::filter(model == "WH-informed"),
      aes(
        x = time,
        ymin = lwr,
        ymax = upr,
        fill = quantity,
        group = quantity
      ),
      inherit.aes = FALSE,
      alpha = 0.10,
      colour = NA,
      show.legend = FALSE
    )
}

p <- p +
  geom_line(linewidth = 1.05) +
  scale_x_continuous(
    name = "Date",
    breaks = date_breaks,
    labels = date_labels
  ) +
  scale_y_log10(
    name = "Reproduction number",
    breaks = c(0.5, 1, 2, 4),
    limits = c(0.5, 4),
    expand = c(0, 0)
  ) +
  scale_colour_manual(
    name = "Quantity",
    values = quantity_cols
  ) +
  scale_fill_manual(
    values = quantity_cols
  ) +
  scale_linetype_manual(
    name = "Model",
    values = model_linetypes,
    breaks = c(
      "WH-informed",
      "Censored WH-informed",
      "Lognormal"
    ),
    labels = c(
      "WH-informed",
      "WH-informed, censored",
      "Lognormal"
    )
  ) +
  guides(
    fill = "none",
    linetype = guide_legend(
      order = 1,
      override.aes = list(
        colour = "grey20",
        linewidth = 1.5
      )
    ),
    colour = guide_legend(
      order = 2,
      override.aes = list(
        linetype = "solid",
        linewidth = 1.3
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Reproduction number"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(2.2, "cm"),
    legend.key.height = grid::unit(0.55, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = paste0(OUT_PREFIX, ".pdf"),
  plot = p,
  width = 8.2,
  height = 4.8
)

ggsave(
  filename = paste0(OUT_PREFIX, ".png"),
  plot = p,
  width = 8.2,
  height = 4.8,
  dpi = 300
)
