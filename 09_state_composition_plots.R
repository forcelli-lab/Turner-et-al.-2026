## Group 2 state composition for Turner et al.
## Isolation / Proximity / Social Contact are mutually exclusive and exhaust the
## scored window, so they are true proportions that sum to 1 and their three
## paired differences sum to 0.
##
## Produces two figures:
##   <prefix>condition.png  - condition means, stacked bars + ribbons,
##                            with the nested permutation difference below
##   <prefix>byanimal.png   - the same composition per animal

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(egg)
})

graphics.off()
set.seed(42)

## ---------------- CONFIG ----------------
data_dir         <- "01_measures"

## One entry per dataset; each is analyzed and plotted independently, into its
## own subdirectory.
datasets         <- c(ACCs = "ACCs_BMI.csv",
                      ACCg = "ACCg_BMI.csv")

out_dir          <- "09_state_composition"
subject_col      <- "Animal"
hemisphere_col   <- "Hemisphere"
treatment_col    <- "Treatment"
reference_level  <- "Saline"
comparison_level <- "BMI"

## ---- The three Group 2 states: column in the sheet -> label on the figure ----
states <- c(Alone         = "Isolation",
            Proximity     = "Proximity",
            SocialContact = "Social contact")

## Figure labels: two-letter code per animal. Keys are Animal column values.
animal_abbrev <- c(Basquiat = "BA",
                   Cobain   = "CO",
                   Hendrix  = "HE",
                   Joplin   = "JO")

## Short forms for the difference-panel x axis; the legend keeps the full names.
states_short <- c(Alone         = "Isolation",
                  Proximity     = "Proximity",
                  SocialContact = "Contact")


## ---- STYLE: all appearance settings ----
STYLE <- list(

  ## ---- Font ----
  font_family = "Arial",
  font_size   = 14,

  ## ---- Panel dimensions (inches) ----
  ## Set target_figure_width_in to 1.5x the width the estimation scripts print,
  ## and animal_panel_width_in to (3 x that width - 2 x axis margin) / 7.
  target_figure_width_in = 6.15,
  panel_height_top_in    = 3.5,
  panel_height_bot_in    = 2.6,
  animal_panel_width_in  = 1.45,
  animal_panel_height_in = 2.8,
  fig_margin_w_in      = 1.6,   # axis title, tick labels and legend
  fig_margin_h_in      = 2.0,
  fig_dpi              = 600,

  ## ---- State colors (bars and ribbons) ----
  state_colors = c("Isolation"      = "#BFBFBF",
                   "Proximity"      = "#9EC9E2",
                   "Social contact" = "#E8705F"),
  ribbon_alpha = 0.40,          # ribbons are the same hue, lighter
  bar_alpha    = 1.00,

  ## ---- Bar geometry (x axis is arbitrary units; bars sit at 0 and 1) ----
  ## fraction of the Saline-to-BMI distance; the gap between bars is
  ## 1 - 2 * bar_half_width
  bar_half_width = 0.30,
  x_pad          = 0.04,        # blank space outside the outer bar edges
  ribbon_n       = 200,         # points along each ribbon; higher = smoother
  bar_outline    = NA,

  ## ---- Difference panel ----
  violin_max_halfwidth = 0.30,
  violin_bulk_color    = "#FAE2AB",   # center of the null
  violin_tail_color    = "#F38D7D",   # rejection tails, test significant
  violin_tail_color_ns = "#C9C9C9",   # rejection tails, test not significant
  alpha_level          = 0.05,
  violin_offset        = 0.16,
  ci_line_width        = 0.66,
  ci_color             = "black",
  point_fill           = "#F5C457",
  point_color          = "black",
  point_size           = 3,
  point_shape          = 21,
  point_stroke         = 1.1,
  zero_line_color      = "grey50",
  zero_line_type       = "dashed",

  ## ---- Stats annotation ----
  # headroom above the tallest violin / CI, as a fraction of the data range,
  # so the two-line label is never clipped
  label_headroom_frac = 0.22,
  label_pad_frac      = 0.05,
  stats_text_size   = 4,
  stats_text_color  = "black",
  stats_p_sig_color = "red",
  stats_label_fill  = "white",

  ## ---- Axes ----
  axis_line_color = "black",
  n_y_breaks      = 6,

  ## ---- CI search ----
  # proportions live on 0-1, so the bracketing step has to start small
  ci_step_min = 0.01,
  ci_tol      = 0.0005
)

## =====================================================================
## One dataset: the unit structure, the exact enumeration and both figures all
## depend on the dataset, so everything below runs per site.
## =====================================================================

run_site <- function(tag, data_path, out_dir) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## ---------------- LOAD DATA & BUILD STRUCTURE ----------------
  df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  unit_col <- "unit"
  df[[unit_col]] <- paste(df[[subject_col]], df[[hemisphere_col]], sep = "_")

  # proportion of the scored window, not of 1800 -- a truncated session is then
  # on the same footing as a full one
  window_s <- rowSums(df[names(states)])
  for (s in names(states)) df[[paste0("p_", s)]] <- df[[s]] / window_s

  units    <- unique(df[[unit_col]])
  subjects <- unique(df[[subject_col]])
  levelA <- reference_level; levelB <- comparison_level

  subject_units <- lapply(subjects, function(s) unique(df[[unit_col]][df[[subject_col]] == s]))
  names(subject_units) <- subjects

  unit_idx <- lapply(units, function(u) which(df[[unit_col]] == u))
  names(unit_idx) <- units

  combo_list <- lapply(units, function(u) {
    idx <- unit_idx[[u]]
    combn(length(idx), sum(df[[treatment_col]][idx] == levelA))
  })
  names(combo_list) <- units

  ## ---------------- nested statistic and exact permutation ----------------
  # unit difference, averaged within animal, then across animals
  compute_stat <- function(vals, treat) {
    unit_stat <- sapply(units, function(u) {
      idx <- unit_idx[[u]]
      mean(vals[idx][treat[idx] == levelB]) - mean(vals[idx][treat[idx] == levelA])
    })
    names(unit_stat) <- units
    mean(sapply(subjects, function(s) mean(unit_stat[subject_units[[s]]])))
  }

  null_distribution <- function(vals) {
    delta_options <- lapply(units, function(u) {
      idx <- unit_idx[[u]]; combos <- combo_list[[u]]
      sapply(seq_len(ncol(combos)), function(k) {
        A_idx <- idx[combos[, k]]; B_idx <- idx[-combos[, k]]
        mean(vals[B_idx]) - mean(vals[A_idx])
      })
    })
    names(delta_options) <- units
    grid <- do.call(expand.grid, delta_options)
    subject_cols <- sapply(subjects, function(s) {
      units_s <- subject_units[[s]]
      if (length(units_s) == 1) grid[[units_s]] else rowMeans(grid[units_s])
    })
    rowMeans(subject_cols)
  }

  run_exact <- function(col) {
    vals <- df[[col]]; treat <- df[[treatment_col]]
    observed <- compute_stat(vals, treat)
    null_vec <- null_distribution(vals)
    list(observed = observed, null = null_vec,
         p = mean(abs(null_vec) >= abs(observed) - 1e-9))
  }

  # test inversion: subtract a candidate effect from the BMI sessions and re-test
  p_for_shift <- function(col, shift) {
    vals <- df[[col]]; treat <- df[[treatment_col]]
    vals_adj <- vals
    vals_adj[treat == levelB] <- vals[treat == levelB] - shift
    null_vec <- null_distribution(vals_adj)
    mean(abs(null_vec) >= abs(compute_stat(vals_adj, treat)) - 1e-9)
  }

  find_ci_bound <- function(col, observed, direction, alpha = 0.05, max_expand = 40) {
    step <- max(abs(observed) * 0.25, STYLE$ci_step_min)
    inner_t <- 0; outer_t <- step; expand <- 0
    while (p_for_shift(col, observed + direction * outer_t) >= alpha) {
      inner_t <- outer_t; outer_t <- outer_t * 2; expand <- expand + 1
      if (expand > max_expand) stop("Failed to bracket CI bound for ", col)
    }
    while (outer_t - inner_t > STYLE$ci_tol) {
      mid_t <- (inner_t + outer_t) / 2
      if (p_for_shift(col, observed + direction * mid_t) >= alpha) inner_t <- mid_t else outer_t <- mid_t
    }
    observed + direction * ((inner_t + outer_t) / 2)
  }

  ## ---------------- nested condition means ----------------
  # unit mean -> animal mean -> condition mean, matching how the statistic is built
  nested_means <- function(group_by_animal = FALSE) {
    pcols <- paste0("p_", names(states))
    um <- df %>%
      group_by(.data[[subject_col]], .data[[unit_col]], .data[[treatment_col]]) %>%
      summarise(across(all_of(pcols), mean), .groups = "drop")
    am <- um %>%
      group_by(.data[[subject_col]], .data[[treatment_col]]) %>%
      summarise(across(all_of(pcols), mean), .groups = "drop")
    if (group_by_animal) return(am)
    am %>%
      group_by(.data[[treatment_col]]) %>%
      summarise(across(all_of(pcols), mean), .groups = "drop")
  }

  # long format with cumulative band edges, ready for stacked bars
  stack_bands <- function(wide, extra_key = NULL) {
    keys <- c(extra_key, treatment_col)
    out <- do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
      p <- unlist(wide[i, paste0("p_", names(states))])
      p <- p / sum(p)                 # guard against accumulated rounding
      top <- cumsum(p); top[length(top)] <- 1
      bottom <- c(0, top[-length(top)])
      data.frame(wide[i, keys, drop = FALSE], row.names = NULL,
                 state = factor(unname(states), levels = unname(states)),
                 ymin = unname(bottom), ymax = unname(top),
                 stringsAsFactors = FALSE)
    }))
    out$x <- ifelse(out[[treatment_col]] == levelA, 0, 1)
    out
  }

  # smoothstep ribbon linking each state's Saline band to its BMI band
  ribbon_df <- function(bands, extra_key = NULL) {
    hw <- STYLE$bar_half_width
    tt <- seq(0, 1, length.out = STYLE$ribbon_n)
    smooth <- tt * tt * (3 - 2 * tt)
    grps <- if (is.null(extra_key)) list(bands) else split(bands, bands[[extra_key]])
    do.call(rbind, lapply(grps, function(b) {
      do.call(rbind, lapply(levels(b$state), function(st) {
        a <- b[b$state == st & b[[treatment_col]] == levelA, ]
        z <- b[b$state == st & b[[treatment_col]] == levelB, ]
        d <- data.frame(x = hw + tt * (1 - 2 * hw),
                        ymin = a$ymin + smooth * (z$ymin - a$ymin),
                        ymax = a$ymax + smooth * (z$ymax - a$ymax),
                        state = factor(st, levels = levels(b$state)),
                        stringsAsFactors = FALSE)
        if (!is.null(extra_key)) d[[extra_key]] <- a[[extra_key]]
        d
      }))
    }))
  }

  base_theme <- theme_classic(base_size = STYLE$font_size, base_family = STYLE$font_family) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = STYLE$axis_line_color),
          axis.ticks = element_line(color = STYLE$axis_line_color),
          text = element_text(family = STYLE$font_family, size = STYLE$font_size),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.margin = margin(t = 0, b = 0),
          legend.key.size = unit(0.32, "cm"),
          legend.spacing.x = unit(0.15, "cm"))

  ## ---------------- FIGURE 1: condition means ----------------
  bands  <- stack_bands(nested_means())
  ribs   <- ribbon_df(bands)
  hw     <- STYLE$bar_half_width

  p_top <- ggplot() +
    geom_ribbon(data = ribs, aes(x = x, ymin = ymin, ymax = ymax, fill = state),
                alpha = STYLE$ribbon_alpha) +
    geom_rect(data = bands,
              aes(xmin = x - hw, xmax = x + hw, ymin = ymin, ymax = ymax, fill = state),
              alpha = STYLE$bar_alpha, color = STYLE$bar_outline) +
    scale_fill_manual(values = STYLE$state_colors) +
    scale_x_continuous(breaks = c(0, 1), labels = c(levelA, levelB),
                       limits = c(-hw - STYLE$x_pad, 1 + hw + STYLE$x_pad)) +
    scale_y_continuous(name = "Proportion of scored time",
                       breaks = scales::pretty_breaks(n = STYLE$n_y_breaks),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL) + base_theme +
    guides(fill = guide_legend(nrow = 1))

  ## difference panel: one column per state
  res <- lapply(names(states), function(s) run_exact(paste0("p_", s)))
  names(res) <- names(states)
  ci  <- lapply(names(states), function(s) {
    o <- res[[s]]$observed
    c(find_ci_bound(paste0("p_", s), o, -1), find_ci_bound(paste0("p_", s), o, +1))
  })
  names(ci) <- names(states)

  xpos <- setNames(seq_along(states), names(states))
  viol <- do.call(rbind, lapply(names(states), function(s) {
    d <- density(res[[s]]$null, n = 512)
    hwid <- d$y / max(d$y) * STYLE$violin_max_halfwidth
    # rejection region, colored by verdict, per state
    cuts <- quantile(res[[s]]$null,
                     c(STYLE$alpha_level / 2, 1 - STYLE$alpha_level / 2), names = FALSE)
    extreme <- d$x < cuts[1] | d$x > cuts[2]
    tail_col <- if (res[[s]]$p < STYLE$alpha_level)
      STYLE$violin_tail_color else STYLE$violin_tail_color_ns
    data.frame(y = d$x, xmin = xpos[[s]] + STYLE$violin_offset,
               xmax = xpos[[s]] + STYLE$violin_offset + hwid,
               extreme = extreme,
               band_fill = ifelse(extreme, tail_col, STYLE$violin_bulk_color),
               seg = paste(s, cumsum(c(1, diff(as.integer(extreme)) != 0))),
               state = factor(states[[s]], levels = unname(states)),
               stringsAsFactors = FALSE)
  }))
  pts <- data.frame(x = unname(xpos),
                    obs = sapply(names(states), function(s) res[[s]]$observed),
                    lo  = sapply(names(states), function(s) ci[[s]][1]),
                    hi  = sapply(names(states), function(s) ci[[s]][2]),
                    p   = sapply(names(states), function(s) res[[s]]$p),
                    state = factor(unname(states), levels = unname(states)),
                    stringsAsFactors = FALSE)
  pts$label <- sprintf("%+.1f pp\np = %.4f", 100 * pts$obs, pts$p)

  # room for the violins, the CI whiskers and the two-line label above them
  span      <- range(c(viol$y, pts$lo, pts$hi, 0))
  y_floor   <- span[1] - STYLE$label_pad_frac * diff(span)
  y_ceiling <- span[2] + STYLE$label_headroom_frac * diff(span)

  p_bot <- ggplot() +
    geom_hline(yintercept = 0, linetype = STYLE$zero_line_type, color = STYLE$zero_line_color) +
    geom_ribbon(data = viol, aes(y = y, xmin = xmin, xmax = xmax, group = seg, fill = band_fill),
                orientation = "y", color = NA) +
    scale_fill_identity() +
    geom_segment(data = pts, aes(x = x, xend = x, y = lo, yend = hi),
                 color = STYLE$ci_color, linewidth = STYLE$ci_line_width) +
    geom_point(data = pts, aes(x = x, y = obs), fill = STYLE$point_fill,
               color = STYLE$point_color, size = STYLE$point_size,
               shape = STYLE$point_shape, stroke = STYLE$point_stroke) +
    geom_text(data = pts, aes(x = x, y = hi, label = label,
                              color = p < 0.05, fontface = ifelse(p < 0.05, "bold", "plain")),
              vjust = -0.4, size = STYLE$stats_text_size, show.legend = FALSE) +
    scale_color_manual(values = c(`FALSE` = STYLE$stats_text_color, `TRUE` = STYLE$stats_p_sig_color),
                       guide = "none") +
    scale_x_continuous(breaks = unname(xpos), labels = unname(states_short),
                       limits = c(0.5, length(states) + 0.9)) +
    scale_y_continuous(name = paste0("Difference in proportion (", levelB, " - ", levelA, ")"),
                       breaks = scales::pretty_breaks(n = STYLE$n_y_breaks)) +
    coord_cartesian(ylim = c(y_floor, y_ceiling), clip = "off") +
    labs(x = NULL) + base_theme + theme(legend.position = "none")

  # panel heights come from the heights ratio; the saved width is fixed 
  f1 <- ggarrange(p_top, p_bot, ncol = 1,
                  heights = c(STYLE$panel_height_top_in, STYLE$panel_height_bot_in),
                  draw = FALSE)

  f1_name <- file.path(out_dir, "condition.png")
  ggsave(f1_name, f1,
         width  = STYLE$target_figure_width_in,
         height = STYLE$panel_height_top_in + STYLE$panel_height_bot_in + STYLE$fig_margin_h_in,
         dpi = STYLE$fig_dpi, units = "in")
  cat(sprintf("saved %s (%.2f x %.2f in)\n", f1_name,
              STYLE$target_figure_width_in,
              STYLE$panel_height_top_in + STYLE$panel_height_bot_in + STYLE$fig_margin_h_in))

  ## ---------------- FIGURE 2: per animal ----------------
  bands_a <- stack_bands(nested_means(group_by_animal = TRUE), extra_key = subject_col)
  ribs_a  <- ribbon_df(bands_a, extra_key = subject_col)

  p_animal <- ggplot() +
    geom_ribbon(data = ribs_a, aes(x = x, ymin = ymin, ymax = ymax, fill = state),
                alpha = STYLE$ribbon_alpha) +
    geom_rect(data = bands_a,
              aes(xmin = x - hw, xmax = x + hw, ymin = ymin, ymax = ymax, fill = state),
              alpha = STYLE$bar_alpha, color = STYLE$bar_outline) +
    facet_wrap(as.formula(paste("~", subject_col)), nrow = 1,
               labeller = labeller(.default = function(x) animal_abbrev[x])) +
    scale_fill_manual(values = STYLE$state_colors) +
    scale_x_continuous(breaks = c(0, 1), labels = c(levelA, levelB),
                       limits = c(-hw - STYLE$x_pad, 1 + hw + STYLE$x_pad)) +
    scale_y_continuous(name = "Proportion of scored time",
                       breaks = scales::pretty_breaks(n = STYLE$n_y_breaks),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL) + base_theme +
    theme(strip.background = element_blank(),
          strip.text = element_text(family = STYLE$font_family, size = STYLE$font_size))

  f2 <- set_panel_size(p_animal,
                       width  = unit(STYLE$animal_panel_width_in, "in"),
                       height = unit(STYLE$animal_panel_height_in, "in"))
  f2_name <- file.path(out_dir, "by_animal.png")
  ggsave(f2_name, f2,
         width  = grid::convertWidth(sum(f2$widths), "in", valueOnly = TRUE),
         height = grid::convertHeight(sum(f2$heights), "in", valueOnly = TRUE),
         dpi = STYLE$fig_dpi, units = "in")
  cat("saved", f2_name, "\n")

  ## ---------------- numbers behind the figures ----------------
  summary_tbl <- data.frame(
    State      = unname(states),
    Saline     = sapply(names(states), function(s) nested_means()[[paste0("p_", s)]][nested_means()[[treatment_col]] == levelA]),
    BMI        = sapply(names(states), function(s) nested_means()[[paste0("p_", s)]][nested_means()[[treatment_col]] == levelB]),
    Difference = sapply(names(states), function(s) res[[s]]$observed),
    CI_low     = sapply(names(states), function(s) ci[[s]][1]),
    CI_high    = sapply(names(states), function(s) ci[[s]][2]),
    P          = sapply(names(states), function(s) res[[s]]$p),
    row.names = NULL)
  print(summary_tbl, row.names = FALSE)
  write.csv(summary_tbl, file.path(out_dir, "results.csv"), row.names = FALSE)
  stopifnot(abs(sum(summary_tbl$Difference)) < 1e-8)

  invisible(summary_tbl)
}

## ---------------- RUN ALL DATASETS ----------------

for (tag in names(datasets)) {
  cat("\n==== ", tag, " ====\n", sep = "")
  run_site(tag,
           file.path(data_dir, datasets[[tag]]),
           file.path(out_dir, tag))
}
