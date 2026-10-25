# ============================================================
# SIMULATION 3
# FAST WINDOWS CRPS EXTRACTION + FIGURE 4
#
# - DOES NOT RERUN STAN
# - Windows-safe parallelisation using future.apply
# - Reads only mean_gen and sd_gen from CmdStan CSVs
# - Reads each chain separately
# - Ignores empty / broken / stale chain files
# - Explicit dplyr::select() everywhere
#
# Panel A:
#   GT mean CRPS
#   GT SD CRPS
#   KL distance to true SI
#   Held-out LPD per pair
#
# Panel B:
#   Posterior movement toward truth
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "future",
  "future.apply",
  "ggplot2",
  "scales",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


suppressPackageStartupMessages({
  
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  
  library(future)
  library(future.apply)
  
  library(ggplot2)
  library(scales)
  library(patchwork)
})


# ============================================================
# 1. PATHS
# ============================================================

OUT_DIR <-
  "~/Within-host-time-delay-framework/Results/Simulations/Simulation 3"

CMDSTAN_OUTPUT_DIR <-
  file.path(
    OUT_DIR,
    "cmdstan_outputs"
  )

FIG_DIR <- OUT_DIR


# ============================================================
# 2. WINDOWS PARALLEL SETTINGS
# ============================================================
#
# Multisession works on Windows.
#
# Do NOT automatically use every CPU core because this task is
# largely disk I/O. 4-6 workers is usually plenty.
#
# Change manually if wanted.
# ============================================================

AVAILABLE_CORES <-
  future::availableCores()

N_WORKERS <-
  max(
    1L,
    min(
      6L,
      AVAILABLE_CORES - 1L
    )
  )


cat(
  "\nAvailable cores:",
  AVAILABLE_CORES,
  "\n"
)

cat(
  "Workers being used:",
  N_WORKERS,
  "\n\n"
)


# ============================================================
# 3. LOAD EXISTING SIMULATION OUTPUT
# ============================================================

fit_results <-
  readr::read_csv(
    file.path(
      OUT_DIR,
      "fit_level_results.csv"
    ),
    show_col_types = FALSE
  )


param_results <-
  readr::read_csv(
    file.path(
      OUT_DIR,
      "parameter_recovery_results.csv"
    ),
    show_col_types = FALSE
  )


true_targets <-
  readr::read_csv(
    file.path(
      OUT_DIR,
      "true_generation_time_targets.csv"
    ),
    show_col_types = FALSE
  )


# ============================================================
# 4. CONSTANTS
# ============================================================

N_REP <- 50L

EXPECTED_CHAINS <- 4L

EXPECTED_DRAWS_PER_CHAIN <- 300L


PRIOR_TARGET_SET <- c(
  "log_r",
  "log_d",
  "m",
  "log_kappa",
  "log_sigma_z"
)


TRUE_MEAN_GEN <-
  true_targets$true_mean_gen[1]

TRUE_SD_GEN <-
  true_targets$true_sd_gen[1]


cat(
  "True GT mean:",
  TRUE_MEAN_GEN,
  "\n"
)

cat(
  "True GT SD:",
  TRUE_SD_GEN,
  "\n\n"
)


# ============================================================
# 5. CRPS FUNCTION
# ============================================================
#
# Empirical CRPS:
#
# E|X - y| - 0.5 E|X - X'|
#
# Efficient O(n log n) implementation.
# ============================================================

crps_draws <- function(
    draws,
    truth
) {
  
  x <-
    as.numeric(
      draws
    )
  
  x <-
    x[
      is.finite(x)
    ]
  
  n <-
    length(x)
  
  
  if (n == 0L) {
    
    return(
      NA_real_
    )
  }
  
  
  x <-
    sort(x)
  
  
  # E|X - truth|
  term1 <-
    mean(
      abs(
        x - truth
      )
    )
  
  
  # 0.5 E|X-X'|
  weights <-
    2 * seq_len(n) -
    n -
    1
  
  
  term2 <-
    sum(
      weights * x
    ) /
    n^2
  
  
  term1 -
    term2
}


# ============================================================
# 6. FIND CMDSTAN CSV FILES FOR ONE FIT
# ============================================================

get_job_csv_files <- function(
    job_id,
    rep_id,
    prior_scenario
) {
  
  job_label <-
    sprintf(
      "job%04d_rep%03d_%s",
      as.integer(job_id),
      as.integer(rep_id),
      prior_scenario
    )
  
  
  job_dir <-
    file.path(
      CMDSTAN_OUTPUT_DIR,
      job_label
    )
  
  
  if (!dir.exists(job_dir)) {
    
    return(
      character(0)
    )
  }
  
  
  files <-
    list.files(
      job_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    )
  
  
  # Exclude diagnostic files
  files <-
    files[
      !grepl(
        "diagnostic",
        basename(files),
        ignore.case = TRUE
      )
    ]
  
  
  files
}


# ============================================================
# 7. FAST READ OF ONE CMDSTAN CHAIN
# ============================================================
#
# CmdStan puts metadata on lines starting with #.
#
# read_csv(comment = "#") skips those.
#
# Only two columns are read:
#
#   mean_gen
#   sd_gen
#
# ============================================================

read_chain_fast <- function(
    file
) {
  
  tryCatch({
    
    dat <-
      readr::read_csv(
        file,
        comment = "#",
        col_select = c(
          "mean_gen",
          "sd_gen"
        ),
        show_col_types = FALSE,
        progress = FALSE,
        lazy = FALSE
      )
    
    
    if (nrow(dat) == 0L) {
      
      return(
        NULL
      )
    }
    
    
    if (!all(
      c(
        "mean_gen",
        "sd_gen"
      ) %in% names(dat)
    )) {
      
      return(
        NULL
      )
    }
    
    
    dat
    
    
  }, error = function(e) {
    
    NULL
    
  })
}


# ============================================================
# 8. PROCESS ONE FIT
# ============================================================

process_one_fit <- function(
    row
) {
  
  job_id <-
    as.integer(
      row$job_id
    )
  
  rep_id <-
    as.integer(
      row$rep_id
    )
  
  prior_scenario <-
    as.character(
      row$prior_scenario
    )
  
  
  # ----------------------------------------------------------
  # Find candidate CSV files
  # ----------------------------------------------------------
  
  files <-
    get_job_csv_files(
      job_id = job_id,
      rep_id = rep_id,
      prior_scenario = prior_scenario
    )
  
  
  if (length(files) == 0L) {
    
    return(
      tibble::tibble(
        
        job_id =
          job_id,
        
        rep_id =
          rep_id,
        
        prior_scenario =
          prior_scenario,
        
        n_csv_found =
          0L,
        
        n_valid_chains =
          0L,
        
        n_draws =
          0L,
        
        mean_gen_crps =
          NA_real_,
        
        sd_gen_crps =
          NA_real_,
        
        crps_ok =
          FALSE
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Sort newest first
  #
  # Useful if old/stale files remain after a rerun.
  # ----------------------------------------------------------
  
  file_info <-
    file.info(
      files
    )
  
  
  files <-
    files[
      order(
        file_info$mtime,
        decreasing = TRUE
      )
    ]
  
  
  # ----------------------------------------------------------
  # Read each chain
  # ----------------------------------------------------------
  
  chain_data <-
    lapply(
      files,
      read_chain_fast
    )
  
  
  valid <-
    !vapply(
      chain_data,
      is.null,
      logical(1)
    )
  
  
  files <-
    files[
      valid
    ]
  
  chain_data <-
    chain_data[
      valid
    ]
  
  
  if (length(chain_data) == 0L) {
    
    return(
      tibble::tibble(
        
        job_id =
          job_id,
        
        rep_id =
          rep_id,
        
        prior_scenario =
          prior_scenario,
        
        n_csv_found =
          length(files),
        
        n_valid_chains =
          0L,
        
        n_draws =
          0L,
        
        mean_gen_crps =
          NA_real_,
        
        sd_gen_crps =
          NA_real_,
        
        crps_ok =
          FALSE
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Prefer chains with exactly 300 posterior draws
  #
  # This matches the original Simulation 3:
  #
  # ITER_SAMPLING = 300
  # ----------------------------------------------------------
  
  chain_n <-
    vapply(
      chain_data,
      nrow,
      integer(1)
    )
  
  
  correct_length <-
    chain_n ==
    EXPECTED_DRAWS_PER_CHAIN
  
  
  if (
    sum(correct_length) >=
    EXPECTED_CHAINS
  ) {
    
    files <-
      files[
        correct_length
      ]
    
    chain_data <-
      chain_data[
        correct_length
      ]
  }
  
  
  # ----------------------------------------------------------
  # If there are old duplicate valid chains, keep newest four.
  #
  # Since files were already ordered newest -> oldest above,
  # just retain first four.
  # ----------------------------------------------------------
  
  if (
    length(chain_data) >
    EXPECTED_CHAINS
  ) {
    
    chain_data <-
      chain_data[
        seq_len(
          EXPECTED_CHAINS
        )
      ]
    
    files <-
      files[
        seq_len(
          EXPECTED_CHAINS
        )
      ]
  }
  
  
  # ----------------------------------------------------------
  # Combine posterior draws
  # ----------------------------------------------------------
  
  draws <-
    dplyr::bind_rows(
      chain_data
    )
  
  
  n_draws <-
    nrow(
      draws
    )
  
  
  if (n_draws == 0L) {
    
    return(
      tibble::tibble(
        
        job_id =
          job_id,
        
        rep_id =
          rep_id,
        
        prior_scenario =
          prior_scenario,
        
        n_csv_found =
          length(files),
        
        n_valid_chains =
          length(chain_data),
        
        n_draws =
          0L,
        
        mean_gen_crps =
          NA_real_,
        
        sd_gen_crps =
          NA_real_,
        
        crps_ok =
          FALSE
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # CRPS
  # ----------------------------------------------------------
  
  mean_gen_crps <-
    crps_draws(
      draws$mean_gen,
      TRUE_MEAN_GEN
    )
  
  
  sd_gen_crps <-
    crps_draws(
      draws$sd_gen,
      TRUE_SD_GEN
    )
  
  
  tibble::tibble(
    
    job_id =
      job_id,
    
    rep_id =
      rep_id,
    
    prior_scenario =
      prior_scenario,
    
    n_csv_found =
      length(files),
    
    n_valid_chains =
      length(chain_data),
    
    n_draws =
      n_draws,
    
    mean_gen_crps =
      mean_gen_crps,
    
    sd_gen_crps =
      sd_gen_crps,
    
    crps_ok =
      is.finite(mean_gen_crps) &
      is.finite(sd_gen_crps)
  )
}


# ============================================================
# 9. CREATE JOB LIST
# ============================================================

jobs <-
  fit_results %>%
  
  dplyr::select(
    job_id,
    rep_id,
    prior_scenario
  ) %>%
  
  dplyr::distinct() %>%
  
  dplyr::arrange(
    job_id
  )


job_rows <-
  split(
    jobs,
    seq_len(
      nrow(jobs)
    )
  )


cat(
  "Fits to process:",
  length(job_rows),
  "\n\n"
)


# ============================================================
# 10. TEST ONE FIT FIRST
# ============================================================

cat(
  "Testing first fit...\n"
)


test_result <-
  process_one_fit(
    job_rows[[1]]
  )


print(
  test_result
)


if (!test_result$crps_ok) {
  
  stop(
    "First fit failed CRPS extraction. ",
    "Stopping before launching parallel workers."
  )
}


# ============================================================
# 11. WINDOWS PARALLEL PROCESSING
# ============================================================
#
# multisession = separate R processes
#
# Works on Windows.
# ============================================================

old_plan <-
  future::plan()


on.exit(
  future::plan(old_plan),
  add = TRUE
)


future::plan(
  future::multisession,
  workers = N_WORKERS
)


cat(
  "\nStarting Windows parallel extraction...\n\n"
)


start_time <-
  Sys.time()


crps_list <-
  future.apply::future_lapply(
    
    job_rows,
    
    process_one_fit,
    
    future.seed = TRUE,
    
    future.scheduling = 2
  )


end_time <-
  Sys.time()


future::plan(
  future::sequential
)


cat(
  "\nExtraction time:",
  round(
    as.numeric(
      difftime(
        end_time,
        start_time,
        units = "mins"
      )
    ),
    2
  ),
  "minutes\n\n"
)


# ============================================================
# 12. COMBINE RESULTS
# ============================================================

crps_results <-
  dplyr::bind_rows(
    crps_list
  )


# ============================================================
# 13. EXTRACTION CHECKS
# ============================================================

cat(
  "CRPS extraction success:\n"
)


print(
  crps_results %>%
    
    dplyr::count(
      crps_ok
    )
)


cat(
  "\nNumber of posterior draws per fit:\n"
)


print(
  crps_results %>%
    
    dplyr::count(
      n_draws
    ) %>%
    
    dplyr::arrange(
      n_draws
    )
)


cat(
  "\nValid chains per fit:\n"
)


print(
  crps_results %>%
    
    dplyr::count(
      n_valid_chains
    )
)


# ============================================================
# 14. STOP IF EXTRACTION FAILED
# ============================================================

success_rate <-
  mean(
    crps_results$crps_ok
  )


cat(
  "\nSuccess rate:",
  round(
    100 * success_rate,
    1
  ),
  "%\n"
)


if (success_rate < 0.95) {
  
  readr::write_csv(
    crps_results,
    file.path(
      OUT_DIR,
      "crps_extraction_diagnostics_WINDOWS.csv"
    )
  )
  
  
  stop(
    "Less than 95% of fits produced CRPS. ",
    "Check crps_extraction_diagnostics_WINDOWS.csv."
  )
}


# ============================================================
# 15. SAVE CRPS RESULTS
# ============================================================

readr::write_csv(
  
  crps_results,
  
  file.path(
    OUT_DIR,
    "crps_results_WINDOWS_FAST.csv"
  )
)


# ============================================================
# 16. MERGE WITH EXISTING FIT RESULTS
# ============================================================

fit_results_crps <-
  fit_results %>%
  
  dplyr::select(
    -dplyr::any_of(
      c(
        "mean_gen_crps",
        "sd_gen_crps"
      )
    )
  ) %>%
  
  dplyr::left_join(
    
    crps_results %>%
      
      dplyr::select(
        job_id,
        mean_gen_crps,
        sd_gen_crps,
        n_draws,
        n_valid_chains
      ),
    
    by =
      "job_id"
  )


# ============================================================
# 17. CHECK ACTUAL CRPS VALUES
# ============================================================

cat(
  "\nGT mean CRPS:\n"
)

print(
  summary(
    fit_results_crps$mean_gen_crps
  )
)


cat(
  "\nGT SD CRPS:\n"
)

print(
  summary(
    fit_results_crps$sd_gen_crps
  )
)


if (all(
  is.na(
    fit_results_crps$mean_gen_crps
  )
)) {
  
  stop(
    "All mean GT CRPS values are NA."
  )
}


if (all(
  is.na(
    fit_results_crps$sd_gen_crps
  )
)) {
  
  stop(
    "All GT SD CRPS values are NA."
  )
}


# ============================================================
# 18. SAVE FIT-LEVEL RESULTS
# ============================================================

readr::write_csv(
  
  fit_results_crps,
  
  file.path(
    OUT_DIR,
    "fit_level_results_with_crps.csv"
  )
)


# ============================================================
# 19. PRIOR-SCENARIO SUMMARY
# ============================================================

fit_summary_by_prior <-
  fit_results_crps %>%
  
  dplyr::group_by(
    prior_scenario,
    accuracy,
    strength
  ) %>%
  
  dplyr::summarise(
    
    n_successful_fits =
      dplyr::n(),
    
    
    # --------------------------------------------------------
    # CRPS
    # --------------------------------------------------------
    
    median_crps_mean_gen =
      median(
        mean_gen_crps,
        na.rm = TRUE
      ),
    
    
    median_crps_sd_gen =
      median(
        sd_gen_crps,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Existing Simulation 3 metrics
    # --------------------------------------------------------
    
    median_test_lpd_per_pair =
      median(
        test_lpd_per_pair,
        na.rm = TRUE
      ),
    
    
    median_kl_true_to_postmean =
      median(
        kl_true_to_postmean,
        na.rm = TRUE
      ),
    
    
    .groups =
      "drop"
  )


cat(
  "\nPrior-level results:\n\n"
)


print(
  fit_summary_by_prior,
  n = Inf
)


# ============================================================
# 20. SAVE SUMMARY
# ============================================================

readr::write_csv(
  
  fit_summary_by_prior,
  
  file.path(
    OUT_DIR,
    "paper_fit_summary_by_prior_CRPS.csv"
  )
)


# ============================================================
# 21. PANEL A DATA
# ============================================================

plot_heat_data <-
  fit_summary_by_prior %>%
  
  dplyr::filter(
    prior_scenario !=
      "manuscript_default"
  ) %>%
  
  dplyr::mutate(
    
    accuracy = factor(
      
      accuracy,
      
      levels = c(
        "accurate",
        "mildly_wrong",
        "badly_wrong"
      ),
      
      labels = c(
        "Accurate prior",
        "Mildly wrong prior",
        "Badly wrong prior"
      )
    ),
    
    
    strength = factor(
      
      strength,
      
      levels = c(
        "weak",
        "moderate",
        "strong"
      ),
      
      labels = c(
        "Weak",
        "Moderate",
        "Strong"
      )
    )
  ) %>%
  
  dplyr::select(
    
    accuracy,
    
    strength,
    
    median_crps_mean_gen,
    
    median_crps_sd_gen,
    
    median_kl_true_to_postmean,
    
    median_test_lpd_per_pair
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      median_crps_mean_gen,
      median_crps_sd_gen,
      median_kl_true_to_postmean,
      median_test_lpd_per_pair
    ),
    
    names_to =
      "metric_raw",
    
    values_to =
      "value"
  ) %>%
  
  dplyr::group_by(
    accuracy,
    metric_raw
  ) %>%
  
  dplyr::mutate(
    
    weak_value =
      value[
        strength == "Weak"
      ][1],
    
    
    improvement_vs_weak =
      dplyr::case_when(
        
        metric_raw ==
          "median_test_lpd_per_pair" ~
          
          value -
          weak_value,
        
        TRUE ~
          
          weak_value -
          value
      ),
    
    
    percent_improvement_vs_weak =
      dplyr::case_when(
        
        # Higher LPD = better
        metric_raw ==
          "median_test_lpd_per_pair" ~
          
          100 *
          (
            value -
              weak_value
          ) /
          abs(
            weak_value
          ),
        
        
        # Lower CRPS / KL = better
        TRUE ~
          
          100 *
          (
            1 -
              value /
              weak_value
          )
      )
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    metric =
      dplyr::recode(
        
        metric_raw,
        
        median_crps_mean_gen =
          "GT mean CRPS",
        
        median_crps_sd_gen =
          "GT SD CRPS",
        
        median_kl_true_to_postmean =
          "KL distance to true SI",
        
        median_test_lpd_per_pair =
          "Held-out LPD per pair"
      ),
    
    
    metric =
      factor(
        
        metric,
        
        levels = c(
          "GT mean CRPS",
          "GT SD CRPS",
          "KL distance to true SI",
          "Held-out LPD per pair"
        )
      ),
    
    
    label =
      dplyr::case_when(
        
        metric ==
          "Held-out LPD per pair" ~
          
          sprintf(
            "%.5f\n%+.3f%%",
            value,
            percent_improvement_vs_weak
          ),
        
        
        TRUE ~
          
          sprintf(
            "%.4f\n%+.1f%%",
            value,
            percent_improvement_vs_weak
          )
      ),
    
    
    panel =
      "A",
    
    parameter =
      NA_character_,
    
    moved_toward_truth_rate =
      NA_real_
  ) %>%
  
  dplyr::select(
    
    panel,
    
    accuracy,
    
    strength,
    
    metric,
    
    value,
    
    weak_value,
    
    improvement_vs_weak,
    
    percent_improvement_vs_weak,
    
    label,
    
    parameter,
    
    moved_toward_truth_rate
  )


# ============================================================
# 22. PANEL B DATA
# ============================================================

plot_move_data <-
  param_results %>%
  
  dplyr::filter(
    
    variable %in%
      PRIOR_TARGET_SET,
    
    prior_scenario !=
      "manuscript_default",
    
    accuracy %in%
      c(
        "mildly_wrong",
        "badly_wrong"
      )
  ) %>%
  
  dplyr::group_by(
    accuracy,
    strength,
    variable
  ) %>%
  
  dplyr::summarise(
    
    moved_toward_truth_rate =
      mean(
        posterior_moved_toward_truth,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::mutate(
    
    accuracy = factor(
      
      accuracy,
      
      levels = c(
        "mildly_wrong",
        "badly_wrong"
      ),
      
      labels = c(
        "Mildly wrong prior",
        "Badly wrong prior"
      )
    ),
    
    
    strength = factor(
      
      strength,
      
      levels = c(
        "weak",
        "moderate",
        "strong"
      ),
      
      labels = c(
        "Weak",
        "Moderate",
        "Strong"
      )
    ),
    
    
    parameter =
      dplyr::recode(
        
        variable,
        
        log_r =
          "r",
        
        log_d =
          "d",
        
        m =
          "m",
        
        log_kappa =
          "κ",
        
        log_sigma_z =
          "σ_z"
      ),
    
    
    parameter =
      factor(
        
        parameter,
        
        levels = c(
          "r",
          "d",
          "m",
          "κ",
          "σ_z"
        )
      ),
    
    
    panel =
      "B",
    
    metric =
      NA_character_,
    
    value =
      NA_real_,
    
    weak_value =
      NA_real_,
    
    improvement_vs_weak =
      NA_real_,
    
    percent_improvement_vs_weak =
      NA_real_,
    
    label =
      NA_character_
  ) %>%
  
  dplyr::select(
    
    panel,
    
    accuracy,
    
    strength,
    
    metric,
    
    value,
    
    weak_value,
    
    improvement_vs_weak,
    
    percent_improvement_vs_weak,
    
    label,
    
    parameter,
    
    moved_toward_truth_rate
  )


# ============================================================
# 23. COMBINE FIGURE DATA
# ============================================================

sim3_panel_AB_plot_data <-
  dplyr::bind_rows(
    
    plot_heat_data,
    
    plot_move_data
  )


readr::write_csv(
  
  sim3_panel_AB_plot_data,
  
  file.path(
    OUT_DIR,
    "sim3_panel_AB_plot_data_CRPS.csv"
  ),
  
  na = ""
)


# ============================================================
# 24. PREPARE PLOT OBJECTS
# ============================================================

plot_heat <-
  sim3_panel_AB_plot_data %>%
  
  dplyr::filter(
    panel == "A"
  )


plot_move <-
  sim3_panel_AB_plot_data %>%
  
  dplyr::filter(
    panel == "B"
  ) %>%
  
  dplyr::mutate(
    
    parameter_plotmath =
      dplyr::recode(
        
        as.character(
          parameter
        ),
        
        "r" =
          "r",
        
        "d" =
          "d",
        
        "m" =
          "m",
        
        "κ" =
          "kappa",
        
        "σ_z" =
          "sigma[z]"
      )
  )


# ============================================================
# 25. PANEL A
# ============================================================

p_heat <-
  ggplot2::ggplot(
    
    plot_heat,
    
    ggplot2::aes(
      x = strength,
      y = accuracy,
      fill =
        percent_improvement_vs_weak
    )
  ) +
  
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(
      label = label
    ),
    
    size = 4.0,
    lineheight = 0.9
  ) +
  
  ggplot2::facet_wrap(
    ~ metric,
    nrow = 2
  ) +
  
  ggplot2::scale_fill_gradient2(
    
    low =
      "firebrick3",
    
    mid =
      "white",
    
    high =
      "steelblue4",
    
    midpoint =
      0,
    
    name =
      "Improvement\nvs weak prior (%)"
  ) +
  
  ggplot2::labs(
    
    x =
      "Prior strength",
    
    y =
      "Prior accuracy"
  ) +
  
  ggplot2::theme_bw(
    base_size = 15
  ) +
  
  ggplot2::theme(
    
    panel.grid =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 15
      ),
    
    axis.title =
      ggplot2::element_text(
        size = 17
      ),
    
    axis.text =
      ggplot2::element_text(
        size = 14
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1,
        size = 14
      ),
    
    legend.title =
      ggplot2::element_text(
        size = 15
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 13
      ),
    
    plot.margin =
      ggplot2::margin(
        8,
        10,
        8,
        8
      )
  )


# ============================================================
# 26. PANEL B
# ============================================================

p_move_clean <-
  ggplot2::ggplot(
    
    plot_move,
    
    ggplot2::aes(
      x = strength,
      y = moved_toward_truth_rate
    )
  ) +
  
  ggplot2::geom_hline(
    
    yintercept =
      0.5,
    
    linetype =
      "dashed",
    
    linewidth =
      0.5
  ) +
  
  ggplot2::geom_col(
    
    width =
      0.65,
    
    fill =
      "grey35"
  ) +
  
  ggplot2::facet_grid(
    
    accuracy ~ parameter_plotmath,
    
    labeller =
      ggplot2::labeller(
        parameter_plotmath =
          ggplot2::label_parsed
      )
  ) +
  
  ggplot2::scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  
  ggplot2::labs(
    
    x =
      "Prior strength",
    
    y =
      "Posterior closer to truth than prior mean"
  ) +
  
  ggplot2::theme_bw(
    base_size = 15
  ) +
  
  ggplot2::theme(
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 15
      ),
    
    axis.title =
      ggplot2::element_text(
        size = 17
      ),
    
    axis.text =
      ggplot2::element_text(
        size = 14
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1,
        size = 14
      ),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    plot.margin =
      ggplot2::margin(
        8,
        10,
        8,
        8
      )
  )


# ============================================================
# 27. COMBINE
# ============================================================

p_AB <-
  p_heat /
  p_move_clean +
  
  patchwork::plot_layout(
    
    heights = c(
      1.15,
      1.0
    )
  ) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  
  ggplot2::theme(
    
    plot.tag =
      ggplot2::element_text(
        face = "bold",
        size = 24
      ),
    
    plot.tag.position =
      c(
        0.01,
        0.98
      )
  )


print(
  p_AB
)


# ============================================================
# 28. SAVE FIGURE
# ============================================================

ggplot2::ggsave(
  
  filename =
    file.path(
      FIG_DIR,
      "figure_sim3_panel_A_B_prior_sensitivity_CRPS.png"
    ),
  
  plot =
    p_AB,
  
  width =
    14,
  
  height =
    14,
  
  dpi =
    400
)


ggplot2::ggsave(
  
  filename =
    file.path(
      FIG_DIR,
      "figure_sim3_panel_A_B_prior_sensitivity_CRPS.pdf"
    ),
  
  plot =
    p_AB,
  
  width =
    14,
  
  height =
    14,
  
  device =
    cairo_pdf
)


# ============================================================
# 29. DONE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "DONE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nCRPS results:\n",
  file.path(
    OUT_DIR,
    "crps_results_WINDOWS_FAST.csv"
  ),
  "\n"
)

cat(
  "\nFit-level results:\n",
  file.path(
    OUT_DIR,
    "fit_level_results_with_crps.csv"
  ),
  "\n"
)

cat(
  "\nPrior summary:\n",
  file.path(
    OUT_DIR,
    "paper_fit_summary_by_prior_CRPS.csv"
  ),
  "\n"
)

cat(
  "\nFigure:\n",
  file.path(
    FIG_DIR,
    "figure_sim3_panel_A_B_prior_sensitivity_CRPS.png"
  ),
  "\n\n"
)