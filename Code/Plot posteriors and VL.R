# ============================================================
# MPOX SIX-PANEL FIGURE
#
# PANELS:
#   A  B
#   C  D
#
#   [shared legend for A-D]
#
#   E  F
#
# A-D:
#   - NO descriptive panel titles/subtitles
#   - ONLY A, B, C, D
#   - EXACTLY 4 major grid cells horizontally
#   - EXACTLY 4 major grid cells vertically
#     (= 5 major tick positions on each axis)
#
# Rows 1-2 colours:
#   RED   = Weak-prior posterior
#   GREEN = Yang-informed prior
#   BLUE  = Yang-informed posterior
#
# E-F:
#   - neutral grey violins
#   - no boxplots
#   - no points
#   - short dashed median line exactly matches
#     violin width at the median height
#   - model names shown on x-axis
#
# Panel E:
#   y-axis 4 to 16 days
#   ticks every 2 days
#
# Panel F:
#   y-axis 0 to 12 days
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)
})


# ============================================================
# 2. BASE DIRECTORY
# ============================================================

BASE_DIR <- "C:/Users/nj7786/Documents/MPOX"

if (!dir.exists(BASE_DIR)) {
  stop(
    "\nBASE_DIR does not exist:\n",
    BASE_DIR
  )
}


# ============================================================
# 3. EXACT INPUT FILES
# ============================================================

comparison_file <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair_yangWHpriors_symrel",
  "WH_prior_posterior_comparison_logr_logd_logkappa_m_plot_data.csv"
)

weak_draw_file <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair",
  "mpox_score9_mechanistic_proportional_incubation_updatedV_obspair_draws.csv"
)

yang_post_draw_file <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair_yangWHpriors_symrel",
  "mpox_score9_yangWHpriors_mechanistic_proportional_incubation_updatedV_obspair_draws.csv"
)


# ============================================================
# 4. FIND EXACT LOGNORMAL DRAW FILE
# ============================================================

lognormal_hits <- list.files(
  BASE_DIR,
  pattern = "^mpox_lognormal_score9_fastconv_draws\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(lognormal_hits) == 0) {
  
  stop(
    "\nCould not find:\n",
    "mpox_lognormal_score9_fastconv_draws.csv\n\n",
    "anywhere below:\n",
    BASE_DIR
  )
}

if (length(lognormal_hits) > 1) {
  
  cat(
    "\nMultiple exact lognormal draw files found:\n"
  )
  
  print(
    lognormal_hits
  )
  
  stop(
    "\nMore than one exact lognormal draw file exists."
  )
}

lognormal_draw_file <- lognormal_hits[1]


# ============================================================
# 5. OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR <- file.path(
  BASE_DIR,
  "stan_outputs_mpox_score9_obspair_yangWHpriors_symrel"
)

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

out_png <- file.path(
  OUTPUT_DIR,
  "WH_parameter_and_generation_time_6panel_ABCDEF.png"
)

out_pdf <- file.path(
  OUTPUT_DIR,
  "WH_parameter_and_generation_time_6panel_ABCDEF.pdf"
)


# ============================================================
# 6. CHECK INPUT FILES
# ============================================================

required_files <- c(
  comparison_file,
  weak_draw_file,
  yang_post_draw_file,
  lognormal_draw_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "\nThe following required files do not exist:\n\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ============================================================
# 7. PARAMETER NAME CLEANER
# ============================================================

clean_parameter <- function(x) {
  
  x <- as.character(x)
  
  z <- tolower(
    trimws(x)
  )
  
  out <- x
  
  out[
    z %in% c(
      "log_r",
      "log r",
      "logr"
    )
  ] <- "log_r"
  
  out[
    z %in% c(
      "log_d",
      "log d",
      "logd"
    )
  ] <- "log_d"
  
  out[
    z %in% c(
      "log_kappa",
      "log kappa",
      "logkappa"
    )
  ] <- "log_kappa"
  
  out[
    z == "m"
  ] <- "m"
  
  out
}


# ============================================================
# 8. SOURCE NAME CLEANER
# ============================================================

clean_source <- function(x) {
  
  x <- as.character(x)
  
  out <- x
  
  weak_hit <-
    grepl(
      "weak|broad",
      x,
      ignore.case = TRUE
    ) &
    grepl(
      "post",
      x,
      ignore.case = TRUE
    )
  
  prior_hit <-
    grepl(
      "yang",
      x,
      ignore.case = TRUE
    ) &
    grepl(
      "prior",
      x,
      ignore.case = TRUE
    )
  
  yang_post_hit <-
    grepl(
      "yang",
      x,
      ignore.case = TRUE
    ) &
    grepl(
      "post",
      x,
      ignore.case = TRUE
    )
  
  out[weak_hit] <-
    "Weak-prior posterior"
  
  out[prior_hit] <-
    "Yang-informed prior"
  
  out[yang_post_hit] <-
    "Yang-informed posterior"
  
  out
}


# ============================================================
# 9. READ SAVED PARAMETER-COMPARISON DATA
# ============================================================

dat <- utils::read.csv(
  comparison_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lower_names <- tolower(
  names(dat)
)


find_column <- function(candidates) {
  
  for (candidate in candidates) {
    
    hit <- which(
      lower_names ==
        tolower(candidate)
    )
    
    if (length(hit) > 0) {
      
      return(
        names(dat)[hit[1]]
      )
    }
  }
  
  NA_character_
}


parameter_col <- find_column(
  c(
    "parameter",
    "variable",
    "param",
    "par"
  )
)

source_col <- find_column(
  c(
    "source",
    "distribution",
    "group",
    "type",
    "fit"
  )
)

value_col <- find_column(
  c(
    "value",
    "draw",
    "x"
  )
)

if (
  is.na(parameter_col) ||
  is.na(source_col) ||
  is.na(value_col)
) {
  
  stop(
    "\nCould not identify parameter/source/value columns.\n\n",
    "Columns present:\n",
    paste(
      names(dat),
      collapse = "\n"
    )
  )
}


# ============================================================
# 10. BUILD PARAMETER DATA
# ============================================================

parameter_dat <- data.frame(
  
  parameter =
    clean_parameter(
      dat[[parameter_col]]
    ),
  
  source =
    clean_source(
      dat[[source_col]]
    ),
  
  value =
    as.numeric(
      dat[[value_col]]
    ),
  
  stringsAsFactors = FALSE
)

theta_names <- c(
  "log_r",
  "log_d",
  "log_kappa",
  "m"
)

parameter_dat <- parameter_dat[
  parameter_dat$parameter %in% theta_names &
    is.finite(parameter_dat$value),
  ,
  drop = FALSE
]

parameter_dat$parameter <- factor(
  parameter_dat$parameter,
  levels = c(
    "log_r",
    "log_d",
    "log_kappa",
    "m"
  )
)

parameter_dat$source <- factor(
  parameter_dat$source,
  levels = c(
    "Weak-prior posterior",
    "Yang-informed prior",
    "Yang-informed posterior"
  )
)


# ============================================================
# 11. UPPER-PANEL COLOURS
# ============================================================

parameter_cols <- c(
  
  "Weak-prior posterior" =
    "#F8766D",
  
  "Yang-informed prior" =
    "#00BA38",
  
  "Yang-informed posterior" =
    "#619CFF"
)


# ============================================================
# 12. GET CONSISTENT GRID BREAKS FOR PARAMETER PANEL
#
# Every A-D panel receives:
#
#   5 major x-axis ticks = 4 horizontal grid cells
#   5 major y-axis ticks = 4 vertical grid cells
#
# The numerical range is still calculated separately for
# each parameter.
# ============================================================

get_parameter_grid <- function(
    panel_dat,
    n_cells = 4
) {
  
  # ----------------------------------------------------------
  # X RANGE
  # ----------------------------------------------------------
  
  x_raw <- range(
    panel_dat$value,
    finite = TRUE
  )
  
  x_span <- diff(
    x_raw
  )
  
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }
  
  
  # Small padding so the density is not touching the edges
  x_min <- x_raw[1] - 0.03 * x_span
  x_max <- x_raw[2] + 0.03 * x_span
  
  
  x_breaks <- seq(
    from = x_min,
    to = x_max,
    length.out = n_cells + 1
  )
  
  
  # ----------------------------------------------------------
  # DETERMINE MAXIMUM DENSITY ACROSS ALL THREE SOURCES
  # ----------------------------------------------------------
  
  source_levels <- levels(
    parameter_dat$source
  )
  
  
  density_maxima <- numeric(0)
  
  
  for (src in source_levels) {
    
    x <- panel_dat$value[
      as.character(
        panel_dat$source
      ) == src
    ]
    
    x <- x[
      is.finite(x)
    ]
    
    
    if (length(x) >= 2) {
      
      dens <- stats::density(
        x,
        na.rm = TRUE
      )
      
      density_maxima <- c(
        density_maxima,
        max(
          dens$y,
          na.rm = TRUE
        )
      )
    }
  }
  
  
  if (length(density_maxima) == 0) {
    
    stop(
      "Could not calculate density range for parameter panel."
    )
  }
  
  
  y_max_raw <- max(
    density_maxima,
    na.rm = TRUE
  )
  
  
  # Add slight headroom above the tallest density
  y_max <- y_max_raw * 1.05
  
  
  # ----------------------------------------------------------
  # EXACTLY FIVE Y TICKS
  # ----------------------------------------------------------
  
  y_breaks <- seq(
    from = 0,
    to = y_max,
    length.out = n_cells + 1
  )
  
  
  list(
    x_limits = c(
      x_min,
      x_max
    ),
    x_breaks = x_breaks,
    y_limits = c(
      0,
      y_max
    ),
    y_breaks = y_breaks
  )
}


# ============================================================
# 13. PARAMETER PANEL FUNCTION
#
# A-D ALL HAVE THE SAME NUMBER OF GRID CELLS.
# ============================================================

make_parameter_panel <- function(
    parameter_name,
    panel_letter
) {
  
  panel_dat <- parameter_dat[
    parameter_dat$parameter ==
      parameter_name,
    ,
    drop = FALSE
  ]
  
  
  grid_info <- get_parameter_grid(
    panel_dat = panel_dat,
    n_cells = 4
  )
  
  
  ggplot(
    panel_dat,
    aes(
      x = value,
      colour = source,
      fill = source
    )
  ) +
    
    geom_density(
      alpha = 0.20,
      linewidth = 1.15,
      na.rm = TRUE
    ) +
    
    scale_colour_manual(
      values = parameter_cols,
      breaks = c(
        "Weak-prior posterior",
        "Yang-informed prior",
        "Yang-informed posterior"
      ),
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = parameter_cols,
      breaks = c(
        "Weak-prior posterior",
        "Yang-informed prior",
        "Yang-informed posterior"
      ),
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = grid_info$x_breaks,
      labels = scales::label_number(
        accuracy = 0.1,
        trim = TRUE
      ),
      expand = expansion(
        mult = 0
      )
    ) +
    
    scale_y_continuous(
      breaks = grid_info$y_breaks,
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(
        mult = 0
      )
    ) +
    
    coord_cartesian(
      xlim = grid_info$x_limits,
      ylim = grid_info$y_limits
    ) +
    
    labs(
      title = NULL,
      subtitle = NULL,
      tag = panel_letter,
      x = NULL,
      y = "Density",
      colour = NULL,
      fill = NULL
    ) +
    
    theme_minimal(
      base_size = 16
    ) +
    
    theme(
      
      plot.title =
        element_blank(),
      
      plot.subtitle =
        element_blank(),
      
      plot.tag = element_text(
        face = "bold",
        size = 20
      ),
      
      plot.tag.position = c(
        0.01,
        0.99
      ),
      
      panel.grid.major = element_line(
        linewidth = 0.6
      ),
      
      panel.grid.minor =
        element_blank(),
      
      axis.text =
        element_text(
          size = 13
        ),
      
      axis.title.y =
        element_text(
          size = 15
        ),
      
      axis.title.x =
        element_blank(),
      
      legend.title =
        element_blank(),
      
      legend.text =
        element_text(
          size = 14
        ),
      
      legend.key.width =
        grid::unit(
          1.35,
          "cm"
        )
    )
}


# ============================================================
# 14. READ GENERATION-TIME POSTERIOR DRAWS
# ============================================================

weak_draws <- utils::read.csv(
  weak_draw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

yang_post_draws <- utils::read.csv(
  yang_post_draw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lognormal_draws <- utils::read.csv(
  lognormal_draw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 15. CHECK REQUIRED GT COLUMNS
# ============================================================

required_gt_cols <- c(
  "mean_gen",
  "sd_gen"
)

check_gt_columns <- function(
    x,
    object_name
) {
  
  missing_cols <- setdiff(
    required_gt_cols,
    names(x)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "\n",
      object_name,
      " is missing:\n",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
}

check_gt_columns(
  weak_draws,
  "weak_draws"
)

check_gt_columns(
  lognormal_draws,
  "lognormal_draws"
)

check_gt_columns(
  yang_post_draws,
  "yang_post_draws"
)


# ============================================================
# 16. BUILD GENERATION-TIME DATA
# ============================================================

lognormal_gt <- data.frame(
  
  model =
    "Lognormal",
  
  mean_gen =
    as.numeric(
      lognormal_draws$mean_gen
    ),
  
  sd_gen =
    as.numeric(
      lognormal_draws$sd_gen
    )
)

weak_gt <- data.frame(
  
  model =
    "Weak-prior WH",
  
  mean_gen =
    as.numeric(
      weak_draws$mean_gen
    ),
  
  sd_gen =
    as.numeric(
      weak_draws$sd_gen
    )
)

yang_gt <- data.frame(
  
  model =
    "Yang-informed WH",
  
  mean_gen =
    as.numeric(
      yang_post_draws$mean_gen
    ),
  
  sd_gen =
    as.numeric(
      yang_post_draws$sd_gen
    )
)

gt_dat <- rbind(
  lognormal_gt,
  weak_gt,
  yang_gt
)

gt_dat <- gt_dat[
  is.finite(gt_dat$mean_gen) &
    is.finite(gt_dat$sd_gen),
  ,
  drop = FALSE
]

gt_dat$model <- factor(
  gt_dat$model,
  levels = c(
    "Lognormal",
    "Weak-prior WH",
    "Yang-informed WH"
  )
)


# ============================================================
# 17. MEAN AND SD DATA
# ============================================================

gt_mean_dat <- data.frame(
  model = gt_dat$model,
  value = gt_dat$mean_gen
)

gt_sd_dat <- data.frame(
  model = gt_dat$model,
  value = gt_dat$sd_gen
)


# ============================================================
# 18. NUMERICAL SUMMARY
# ============================================================

gt_summary <- gt_dat %>%
  group_by(
    model
  ) %>%
  summarise(
    
    mean_median =
      median(
        mean_gen
      ),
    
    mean_q025 =
      quantile(
        mean_gen,
        0.025
      ),
    
    mean_q975 =
      quantile(
        mean_gen,
        0.975
      ),
    
    sd_median =
      median(
        sd_gen
      ),
    
    sd_q025 =
      quantile(
        sd_gen,
        0.025
      ),
    
    sd_q975 =
      quantile(
        sd_gen,
        0.975
      ),
    
    .groups = "drop"
  )

print(
  gt_summary
)


# ============================================================
# 19. EXACT VIOLIN WIDTH AT MEDIAN
# ============================================================

get_median_violin_segments <- function(
    data,
    violin_width = 0.82
) {
  
  median_dat <- data %>%
    group_by(
      model
    ) %>%
    summarise(
      median_value =
        median(
          value,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
  
  
  violin_only <- ggplot(
    data,
    aes(
      x = model,
      y = value
    )
  ) +
    
    geom_violin(
      width = violin_width,
      trim = FALSE,
      scale = "area"
    )
  
  
  built <- ggplot_build(
    violin_only
  )
  
  
  violin_geom <- built$data[[1]]
  
  
  result <- vector(
    "list",
    length(levels(data$model))
  )
  
  
  model_levels <- levels(
    data$model
  )
  
  
  for (i in seq_along(model_levels)) {
    
    this_model <- model_levels[i]
    
    this_median <- median_dat$median_value[
      as.character(
        median_dat$model
      ) == this_model
    ]
    
    
    vg <- violin_geom[
      violin_geom$group == i,
      ,
      drop = FALSE
    ]
    
    
    if (nrow(vg) == 0) {
      
      stop(
        "Could not recover violin geometry for ",
        this_model
      )
    }
    
    
    vg <- vg[
      order(vg$y),
      ,
      drop = FALSE
    ]
    
    
    width_at_median <- stats::approx(
      x = vg$y,
      y = vg$violinwidth,
      xout = this_median,
      rule = 2
    )$y
    
    
    half_width <-
      (violin_width / 2) *
      width_at_median
    
    
    result[[i]] <- data.frame(
      
      model =
        this_model,
      
      median_value =
        this_median,
      
      x =
        i,
      
      x_start =
        i - half_width,
      
      x_end =
        i + half_width,
      
      stringsAsFactors = FALSE
    )
  }
  
  
  do.call(
    rbind,
    result
  )
}


# ============================================================
# 20. BOTTOM-PANEL FUNCTION
# ============================================================

make_gt_panel <- function(
    data,
    ylim,
    panel_letter,
    y_breaks
) {
  
  VIOLIN_WIDTH <- 0.82
  
  
  median_segments <- get_median_violin_segments(
    data = data,
    violin_width = VIOLIN_WIDTH
  )
  
  
  ggplot(
    data,
    aes(
      x = model,
      y = value
    )
  ) +
    
    geom_violin(
      fill = "grey72",
      colour = "grey30",
      width = VIOLIN_WIDTH,
      linewidth = 0.75,
      trim = FALSE,
      scale = "area"
    ) +
    
    geom_segment(
      data = median_segments,
      aes(
        x = x_start,
        xend = x_end,
        y = median_value,
        yend = median_value
      ),
      inherit.aes = FALSE,
      linetype = "11",
      linewidth = 0.8,
      colour = "black",
      lineend = "butt"
    ) +
    
    scale_y_continuous(
      breaks = y_breaks
    ) +
    
    coord_cartesian(
      ylim = ylim
    ) +
    
    labs(
      title = NULL,
      subtitle = NULL,
      tag = panel_letter,
      x = NULL,
      y = "Days"
    ) +
    
    theme_minimal(
      base_size = 16
    ) +
    
    theme(
      
      plot.title =
        element_blank(),
      
      plot.subtitle =
        element_blank(),
      
      plot.tag = element_text(
        face = "bold",
        size = 20
      ),
      
      plot.tag.position = c(
        0.01,
        0.99
      ),
      
      panel.grid.minor =
        element_blank(),
      
      axis.text.x =
        element_text(
          size = 14,
          angle = 0,
          hjust = 0.5
        ),
      
      axis.text.y =
        element_text(
          size = 13
        ),
      
      axis.title.y =
        element_text(
          size = 15
        ),
      
      axis.title.x =
        element_blank(),
      
      legend.position =
        "none"
    )
}


# ============================================================
# 21. CREATE PANELS A-D
#
# ALL FOUR HAVE:
#   5 x ticks -> 4 grid cells
#   5 y ticks -> 4 grid cells
# ============================================================

p_A <- make_parameter_panel(
  parameter_name = "log_r",
  panel_letter = "A"
)

p_B <- make_parameter_panel(
  parameter_name = "log_d",
  panel_letter = "B"
)

p_C <- make_parameter_panel(
  parameter_name = "log_kappa",
  panel_letter = "C"
)

p_D <- make_parameter_panel(
  parameter_name = "m",
  panel_letter = "D"
)


# ============================================================
# 22. PANEL E
#
# More y-axis labels:
# 4, 6, 8, 10, 12, 14, 16
# ============================================================

p_E <- make_gt_panel(
  data = gt_mean_dat,
  ylim = c(
    4,
    16
  ),
  panel_letter = "E",
  y_breaks = seq(
    4,
    16,
    by = 2
  )
)


# ============================================================
# 23. PANEL F
# ============================================================

p_F <- make_gt_panel(
  data = gt_sd_dat,
  ylim = c(
    0,
    12
  ),
  panel_letter = "F",
  y_breaks = seq(
    0,
    12,
    by = 2
  )
)


# ============================================================
# 24. FINAL LAYOUT
#
# A B
# C D
# G G  <- shared legend
# E F
# ============================================================

layout_design <- "
AB
CD
GG
EF
"


final_plot <- patchwork::wrap_plots(
  
  A = p_A,
  B = p_B,
  C = p_C,
  D = p_D,
  
  G = patchwork::guide_area(),
  
  E = p_E,
  F = p_F,
  
  design = layout_design,
  
  guides = "collect",
  
  heights = c(
    1,
    1,
    0.18,
    1.10
  )
  
) &
  
  theme(
    
    legend.position =
      "bottom",
    
    legend.title =
      element_blank(),
    
    legend.text =
      element_text(
        size = 14
      ),
    
    legend.key.width =
      grid::unit(
        1.35,
        "cm"
      )
  )


# ============================================================
# 25. DISPLAY
# ============================================================

print(
  final_plot
)


# ============================================================
# 26. SAVE PNG
# ============================================================

ggsave(
  filename = out_png,
  plot = final_plot,
  width = 13,
  height = 13,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 27. SAVE PDF
# ============================================================

ggsave(
  filename = out_pdf,
  plot = final_plot,
  width = 13,
  height = 13,
  bg = "white"
)


# ============================================================
# 28. SAVE NUMERICAL SUMMARY
# ============================================================

utils::write.csv(
  gt_summary,
  file.path(
    OUTPUT_DIR,
    "generation_time_three_model_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 29. DONE
# ============================================================

cat(
  "\n====================================\n",
  "DONE\n",
  "====================================\n\n",
  "PNG:\n",
  out_png,
  "\n\n",
  "PDF:\n",
  out_pdf,
  "\n"
)

# ============================================================
# MANUSCRIPT-READY GENERATION-TIME SUMMARIES
# ============================================================

cat("\n============================================\n")
cat("GENERATION-TIME POSTERIOR SUMMARIES\n")
cat("============================================\n\n")

gt_summary %>%
  arrange(model) %>%
  rowwise() %>%
  do({
    x <- .
    
    cat(
      as.character(x$model), "\n",
      "  Mean generation time: ",
      sprintf("%.2f", x$mean_median),
      " days (95% CrI: ",
      sprintf("%.2f", x$mean_q025),
      "--",
      sprintf("%.2f", x$mean_q975),
      ")\n",
      "  SD generation time:   ",
      sprintf("%.2f", x$sd_median),
      " days (95% CrI: ",
      sprintf("%.2f", x$sd_q025),
      "--",
      sprintf("%.2f", x$sd_q975),
      ")\n\n",
      sep = ""
    )
    
    x
  }) %>%
  invisible()