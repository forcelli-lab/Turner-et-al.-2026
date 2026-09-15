## Plots for Turner et al. supplement -- ACCs vs ACCg comparison, BMI condition
## only, collapsed across hemisphere because it is somewhat unbalanced (unit => Animal). 
## Shows individual infusions, animal mean spaghetti, condition (site) mean, and permutation
## based mean difference estimate.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(egg)
  library(cowplot)
  library(randomizr)
})

graphics.off()

set.seed(42)

## ---------------- CONFIG ----------------
data_dir         <- "01_measures"
data_path        <- file.path(data_dir, "ACCgvss_BMI.csv")

out_dir          <- "08_plots_GvsS"
subject_col      <- "Animal"

filter_col       <- "Treatment"
filter_value     <- "BMI"          # restrict to BMI sessions only

treatment_col    <- "Site"
reference_level  <- "ACCg"
comparison_level <- "ACCs"

max_exact_assignments <- 1e6
comparison_tolerance  <- sqrt(.Machine$double.eps)

## ---- STYLE: all appearance settings ----
STYLE <- list(
  
  ## ---- Font ----
  font_family = "Arial",
  font_size   = 14,
  
  ## ---- Panel dimensions (inches) ----
  panel_width_in  = 2.5,
  panel_height_in = 3.5,
  fig_dpi         = 600,
  
  ### -------- ANIMAL SPECIFICATIONS --------
  ## ---- Per-animal colors; Joplin has no gyrus BMI sessions ----
  animal_colors = c(
    Basquiat = "#1DB100",  # green
    Cobain   = "#0351FF",  # blue
    Hendrix  = "#F42300"   # red
  ),
  
  ## ---- Per-animal symbols ----
  animal_pch = c(
    Basquiat = 5,   # diamond
    Cobain   = 1,   # circle
    Hendrix  = 0    # square
  ),
  ### -------- END OF ANIMAL SPECIFICATIONS --------
  
  ## ---- Raw session points (outermost layer) ----
  raw_point_size    = 2.6,
  raw_point_stroke  = 0.55,
  raw_point_alpha   = 0.45,
  raw_jitter_width  = 0.16,
  
  ## ---- Animal-level means (now the only spaghetti layer) ----
  animal_offset       = 0.60,
  animal_line_width   = 0.75,
  animal_point_size   = 3.25,
  animal_point_stroke = 1,
  animal_point_shape  = 21,      # pch 21: white fill + colored outline
  animal_legend_key_size = 3,
  
  ## ---- Condition (group) means (innermost layer) ----
  condition_offset     = 0.825,
  condition_color      = "#797979",
  condition_line_width = 0.9,
  condition_line_type  = "dotted",
  condition_point_size = 2.5,
  condition_point_stroke = 1,
  condition_point_shape = 8,     # asterisk
  
  ## ---- Column spacing (Reference site / Comparison site / Delta) ----
  ref_comp_gap     = 2.6,   # needs > 2*condition_offset = 1.65 to avoid collision
  comp_delta_gap   = 0.85,  # gap between comparison-site and Delta columns
  xlim_left_pad    = 0.25,
  xlim_right_pad   = 0.5,
  
  ## ---- Delta-panel violin (permutation null, two-tone p-value shading) ----
  violin_max_halfwidth = 0.45,
  violin_bulk_color       = "#FAE2AB",  # center of the null
  violin_tail_color       = "#F38D7D",  # rejection tails, test significant
  violin_tail_color_ns    = "#C9C9C9",  # rejection tails, test not significant
  alpha_level             = 0.05,
  violin_offset_from_delta = 0.15,
  
  ## ---- Mean-difference CI + point (in the Delta panel) ----
  ci_line_width      = 0.66,
  ci_color           = "black",
  point_fill         = "#F5C457", 
  point_color        = "black",
  point_size         = 3,
  point_shape        = 21,
  point_stroke       = 1.1,
  
  ## ---- Reference / guide lines ----
  dashed_line_color = "grey50",
  dashed_line_type  = "dashed",
  
  ## ---- Stats annotation text ----
  stats_text_color        = "black",
  stats_p_sig_color       = "red", # only red when p<0.05
  stats_text_size         = 4,
  stats_text_y_pad_frac   = 0.03,
  stats_label_fill        = "white",
  stats_label_border      = 0,
  stats_label_padding     = 0.3,
  stats_label_right_margin_in = 0.9,  # extra margin so the label doesn't get cut off
  
  ## ---- Axis ----
  axis_line_color  = "black",
  n_y_breaks       = 10,
  n_mdiff_breaks   = 12,
  
  ## ---- white fill behind the animal points so the spaghetti line stops at the shape's edge ----
  mask_color = "white"
)

## Figures label animals by their two-letter code, both letters capitalized,
## not the full name stored in the data. Keys are the values in the Animal
## column.
animal_abbrev <- c(Basquiat = "BA",
                   Cobain   = "CO",
                   Hendrix  = "HE",
                   Joplin   = "JO")

## ---------------- LOAD DATA & BUILD STRUCTURE ----------------
df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)

## Restrict to the BMI sessions
if (!filter_col %in% names(df)) {
  stop("Filter column '", filter_col, "' not found in the data.")
}
df <- df[df[[filter_col]] == filter_value, ]
if (nrow(df) == 0) {
  stop("No rows remain after filtering ", filter_col, " == '", filter_value, "'. Check filter_value/spelling.")
}

## Unit = Animal only (hemisphere pooled/ignored for grouping)
unit_col <- "unit"
df[[unit_col]] <- df[[subject_col]]

units    <- unique(df[[unit_col]])
levelA <- reference_level; levelB <- comparison_level

unit_idx <- lapply(units, function(u) which(df[[unit_col]] == u))
names(unit_idx) <- units

## ---------------- EXACT PERMUTATION TEST ----------------
## Reproduces the same permutation test done in the NestedPermutationTest_GvsS.R script
units_sorted <- sort(unique(as.character(df[[unit_col]])))
unit_counts <- table(
  factor(df[[unit_col]], levels = units_sorted),
  factor(df[[treatment_col]], levels = c(levelA, levelB))
)
if (any(unit_counts == 0L))
  stop("Every animal must contain both sites (", levelA, " and ", levelB, ") among its BMI sessions.")

design        <- declare_ra(blocks = as.character(df[[unit_col]]),
                            block_m = as.numeric(unit_counts[, levelB]))
n_assignments <- obtain_num_permutations(design)
if (!is.finite(n_assignments) || n_assignments > max_exact_assignments)
  stop("Exact space (", format(n_assignments, big.mark = ","), ") exceeds max_exact_assignments.")

assignment_matrix <- obtain_permutation_matrix(design, maximum_permutations = n_assignments)
if (ncol(assignment_matrix) != n_assignments)
  stop("randomizr did not return the full exact enumeration.")

# weight matrix: one row per assignment; permutation_weights %*% values = site stat (1 = ACCs)
permutation_weights <- matrix(0, nrow = n_assignments, ncol = nrow(df))
unit_weight <- 1 / length(units)
for (u in units) {
  idx <- unit_idx[[u]]
  a   <- assignment_matrix[idx, , drop = FALSE]
  permutation_weights[, idx] <- t(ifelse(
    a == 1L,
    unit_weight / unit_counts[u, levelB],
    -unit_weight / unit_counts[u, levelA]
  ))
}

## site statistic: per-animal (ACCs - ACCg) mean difference, averaged equally across animals.
compute_stat <- function(vals, treat) {
  unit_stat <- sapply(units, function(u) {
    idx <- unit_idx[[u]]
    mean(vals[idx][treat[idx] == levelB]) - mean(vals[idx][treat[idx] == levelA])
  })
  names(unit_stat) <- units
  mean(unit_stat)
}

# same test as the stats script (randomizr null), and returns the null vector for the violin
run_exact_with_null <- function(col) {
  vals <- df[[col]]; treat <- df[[treatment_col]]
  observed <- compute_stat(vals, treat)
  null_vec <- as.vector(permutation_weights %*% vals)
  tol <- comparison_tolerance * max(1, abs(observed))
  list(observed = observed, null = null_vec, p = mean(abs(null_vec) >= abs(observed) - tol))
}

## ---------------- permutation CI by test inversion ----------------
# CI is the range of shifts that, subtracted from the ACCs sessions, leave
# the exact test unrejected (p >= 0.05).
p_for_shift <- function(col, shift) {
  vals <- df[[col]]; treat <- df[[treatment_col]]
  vals_adj <- vals
  vals_adj[treat == levelB] <- vals[treat == levelB] - shift
  observed_adj <- compute_stat(vals_adj, treat)
  null_vec <- as.vector(permutation_weights %*% vals_adj)
  tol <- comparison_tolerance * max(1, abs(observed_adj))
  mean(abs(null_vec) >= abs(observed_adj) - tol)
}

find_ci_bound <- function(col, observed, direction, alpha = 0.05,
                          tol = 0.01, max_expand = 40) {
  step <- max(abs(observed) * 0.25, 1)
  inner_t <- 0; outer_t <- step
  expand_count <- 0
  while (p_for_shift(col, observed + direction * outer_t) >= alpha) {
    inner_t <- outer_t; outer_t <- outer_t * 2
    expand_count <- expand_count + 1
    if (expand_count > max_expand) stop("Failed to bracket CI bound for ", col)
  }
  while (outer_t - inner_t > tol) {
    mid_t <- (inner_t + outer_t) / 2
    if (p_for_shift(col, observed + direction * mid_t) >= alpha) inner_t <- mid_t else outer_t <- mid_t
  }
  observed + direction * ((inner_t + outer_t) / 2)
}

get_permutation_ci <- function(col, observed, alpha = 0.05) {
  list(lower = find_ci_bound(col, observed, -1, alpha),
       upper = find_ci_bound(col, observed, +1, alpha))
}

## ---------------- half-violin, tails shaded by verdict ----------------
make_half_violin <- function(null_vals, reference_mean, observed, p_value, xpos = 3,
                             max_halfwidth = 0.55, side = 1, alpha = 0.05) {
  # density is not drawn past the most extreme arrangement
  d <- density(null_vals, n = 512, from = min(null_vals), to = max(null_vals))
  v  <- d$x
  wd <- d$y / max(d$y) * max_halfwidth

  # Rejection region: the outer alpha/2 of the null in each direction.
  cuts <- quantile(null_vals, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  lo <- max(min(v), cuts[1])
  hi <- min(max(v), cuts[2])

  # each band carries its own end points, width interpolated there, so
  # neighboring bands abut exactly
  seg <- function(y0, y1, region) {
    if (!is.finite(y0) || !is.finite(y1) || y1 <= y0) return(NULL)
    keep <- v > y0 & v < y1
    yy <- c(y0, v[keep], y1)
    ww <- c(approx(v, wd, y0)$y, wd[keep], approx(v, wd, y1)$y)
    data.frame(y = yy, xmin = xpos, xmax = xpos + side * ww,
               extreme = region == "extreme",
               seg_id = paste(region, y0, sep = "_"))
  }

  rbind(seg(min(v), lo, "extreme"),
        seg(lo, hi, "central"),
        seg(hi, max(v), "extreme"))
}

## ---------------- BUILD ONE FIGURE ----------------
## show_legend = TRUE only when harvesting the shared legend.
make_plot <- function(col, axis_label, show_legend = FALSE) {
  
  result <- run_exact_with_null(col)
  observed <- result$observed
  
  # Named x-positions for the three columns: reference site, comparison
  # site, Delta.
  x_ref    <- 1
  x_comp   <- x_ref + STYLE$ref_comp_gap
  x_delta  <- x_comp + STYLE$comp_delta_gap
  
  raw_df <- data.frame(Value = df[[col]], Treatment = df[[treatment_col]],
                       Subject = df[[subject_col]])
  raw_df$xpos <- ifelse(raw_df$Treatment == levelA, x_ref, x_comp)
  
  animal_offset    <- STYLE$animal_offset
  condition_offset <- STYLE$condition_offset
  
  # animal mean = average of that animal's BMI sessions at that site
  animal_means <- raw_df %>%
    group_by(Subject, Treatment) %>%
    summarise(MeanValue = mean(Value), .groups = "drop") %>%
    mutate(xpos = ifelse(Treatment == levelA, x_ref + animal_offset, x_comp - animal_offset))
  
  # site mean = average of the animal means, equal weight per animal
  group_means <- animal_means %>%
    group_by(Treatment) %>%
    summarise(MeanValue = mean(MeanValue), .groups = "drop") %>%
    mutate(xpos = ifelse(Treatment == levelA, x_ref + condition_offset, x_comp - condition_offset))
  
  # this is the "zero effect" reference and what the delta axis is rescaled against
  reference_mean <- group_means$MeanValue[group_means$Treatment == levelA]
  null_rescaled  <- reference_mean + result$null
  obs_rescaled   <- reference_mean + observed
  
  ci <- get_permutation_ci(col, observed)
  ci_lower_rescaled <- reference_mean + ci$lower
  ci_upper_rescaled <- reference_mean + ci$upper
  
  # floor at 0 normally, but drops lower if the violin or CI actually needs it
  y_floor <- min(0, min(null_rescaled, na.rm = TRUE), ci_lower_rescaled)
  
  violin_df <- make_half_violin(null_rescaled, reference_mean, observed, p_value = result$p,
                                xpos = x_delta, max_halfwidth = STYLE$violin_max_halfwidth, side = 1)
  
  ci_xpos <- x_delta - STYLE$violin_offset_from_delta
  
  # plot label centered over the whole plot
  label_xpos <- (x_ref + x_delta) / 2
  
  y_ceiling <- max(c(raw_df$Value, null_rescaled, ci_upper_rescaled), na.rm = TRUE)
  
  # figure out how much vertical room (in data units) the 3-line stats box
  # needs, based on its actual physical size in inches, so it doesn't
  # collide with the plot on measures with a small y-range
  n_stats_lines   <- 3
  line_height_in  <- (STYLE$stats_text_size * 1.35) / 25.4
  box_padding_in  <- (STYLE$stats_label_padding * STYLE$stats_text_size * 1.2) / 25.4
  needed_label_height_in <- n_stats_lines * line_height_in + 2 * box_padding_in
  
  data_units_per_in <- (y_ceiling - y_floor) / STYLE$panel_height_in
  headroom_data      <- needed_label_height_in * data_units_per_in * 1.15
  gap_above_data     <- 0.05 * (y_ceiling - y_floor)
  
  label_y_pos          <- y_ceiling + gap_above_data
  y_ceiling_with_label <- label_y_pos + headroom_data
  
  # fixes the geom_jitter scatter so a panel redraws identically
  # tail color marks the verdict: red if the observed effect lands in the
  # rejection region, grey otherwise
  tail_fill <- if (result$p < STYLE$alpha_level)
    STYLE$violin_tail_color else STYLE$violin_tail_color_ns

  set.seed(1)
  p <- ggplot() +
    # raw sessions
    geom_jitter(data = raw_df, aes(x = xpos, y = Value, color = Subject, shape = Subject),
                width = STYLE$raw_jitter_width, height = 0, size = STYLE$raw_point_size,
                alpha = STYLE$raw_point_alpha, stroke = STYLE$raw_point_stroke) +
    # bold lines
    geom_line(data = animal_means,
              aes(x = xpos, y = MeanValue, group = Subject, color = Subject),
              linewidth = STYLE$animal_line_width) +
    geom_point(data = animal_means,
               aes(x = xpos, y = MeanValue, color = Subject),
               fill = STYLE$mask_color, size = STYLE$animal_point_size,
               shape = STYLE$animal_point_shape, stroke = STYLE$animal_point_stroke) +
    # Condition (site) mean spaghetti. Dotted connecting line to be visually distinct from the solid animal spaghetti lines.
    geom_line(data = group_means, aes(x = xpos, y = MeanValue, group = 1),
              color = STYLE$condition_color, linewidth = STYLE$condition_line_width,
              linetype = STYLE$condition_line_type) +
    geom_point(data = group_means, aes(x = xpos, y = MeanValue),
               fill = STYLE$condition_color, color = STYLE$condition_color,
               size = STYLE$condition_point_size, shape = STYLE$condition_point_shape,
               stroke = STYLE$condition_point_stroke) +
    scale_color_manual(values = STYLE$animal_colors,
                       labels = function(x) animal_abbrev[x]) +
    scale_shape_manual(values = STYLE$animal_pch,
                       labels = function(x) animal_abbrev[x]) +
    geom_hline(yintercept = reference_mean, linetype = STYLE$dashed_line_type,
               color = STYLE$dashed_line_color) +
    # dotted guide from the ACCs site mean to the observed effect in the delta panel
    geom_segment(aes(x = group_means$xpos[group_means$Treatment == levelB],
                     xend = ci_xpos, y = obs_rescaled, yend = obs_rescaled),
                 linetype = STYLE$dashed_line_type, color = STYLE$dashed_line_color) +
    geom_ribbon(data = violin_df,
                aes(y = y, xmin = xmin, xmax = xmax, fill = extreme, group = seg_id),
                orientation = "y", color = NA) +
    geom_segment(aes(x = ci_xpos, xend = ci_xpos, y = ci_lower_rescaled, yend = ci_upper_rescaled),
                 color = STYLE$ci_color, linewidth = STYLE$ci_line_width) +
    geom_point(aes(x = ci_xpos, y = obs_rescaled),
               fill = STYLE$point_fill, color = STYLE$point_color,
               size = STYLE$point_size, shape = STYLE$point_shape, stroke = STYLE$point_stroke) +
    # one 3-line stats box (p / obs / CI); whole box goes red+bold if p<0.05
    geom_label(data = data.frame(x = label_xpos, y = label_y_pos,
                                 label = sprintf("p = %.4f\nobs. = %.1f\n95%% CI [%.1f, %.1f]",
                                                 result$p, observed, ci$lower, ci$upper)),
               aes(x = x, y = y, label = label),
               color = if (result$p < 0.05) STYLE$stats_p_sig_color else STYLE$stats_text_color,
               fontface = if (result$p < 0.05) "bold" else "plain",
               hjust = 0.5, vjust = 0, size = STYLE$stats_text_size,
               fill = STYLE$stats_label_fill, linewidth = STYLE$stats_label_border,
               label.padding = unit(STYLE$stats_label_padding, "lines")) +
    scale_x_continuous(breaks = c(x_ref, x_comp, x_delta), labels = c(levelA, levelB, "\u0394"),
                       limits = c(x_ref - STYLE$xlim_left_pad, x_delta + STYLE$xlim_right_pad)) +
    # clip="off" so wide CI text can spill into the right margin instead of getting cut
    coord_cartesian(ylim = c(y_floor, y_ceiling_with_label), clip = "off") +
    scale_fill_manual(values = c(`FALSE` = STYLE$violin_bulk_color, `TRUE` = tail_fill),
                      labels = c("bulk of null", "rejection region (\u03b1 = 0.05)"),
                      name = NULL) +
    scale_y_continuous(
      name = axis_label,
      breaks = scales::pretty_breaks(n = STYLE$n_y_breaks),
      sec.axis = sec_axis(~ . - reference_mean,
                          name = paste0("Mean difference (", levelB, " - ", levelA, ")"),
                          breaks = scales::extended_breaks(n = STYLE$n_mdiff_breaks))
    ) +
    labs(x = NULL) +
    theme_classic(base_size = STYLE$font_size, base_family = STYLE$font_family) +
    theme(legend.position = if (show_legend) "bottom" else "none",
          legend.box = "vertical",
          panel.grid = element_blank(),
          axis.line = element_line(color = STYLE$axis_line_color),
          axis.ticks = element_line(color = STYLE$axis_line_color),
          text = element_text(family = STYLE$font_family, size = STYLE$font_size),
          axis.title = element_text(family = STYLE$font_family, size = STYLE$font_size),
          axis.text = element_text(family = STYLE$font_family, size = STYLE$font_size),
          legend.title = element_text(family = STYLE$font_family, size = STYLE$font_size),
          legend.text = element_text(family = STYLE$font_family, size = STYLE$font_size),
          plot.margin = margin(t = 0.15, r = STYLE$stats_label_right_margin_in,
                               b = 0.05, l = 0.05, unit = "in")) +
    guides(color = guide_legend(title = "Animal", order = 1,
                                override.aes = list(size = STYLE$animal_legend_key_size)),
           shape = guide_legend(title = "Animal", order = 1,
                                override.aes = list(size = STYLE$animal_legend_key_size)),
           fill  = guide_legend(order = 2))
  
  p
}

## ---------------- extract & save legend once, separately ----------------
# same legend on every measure, so grab it once instead of repeating it
extract_legend <- function(p_with_legend) {
  g <- ggplotGrob(p_with_legend)
  is_legend <- sapply(g$grobs, function(x) x$name) == "guide-box"
  g$grobs[[which(is_legend)]]
}

save_legend_png <- function(legend_grob, fname, dpi = 600, pad_in = 0.15) {
  width_in  <- grid::convertWidth(grid::grobWidth(legend_grob), "in", valueOnly = TRUE) + 2 * pad_in
  height_in <- grid::convertHeight(grid::grobHeight(legend_grob), "in", valueOnly = TRUE) + 2 * pad_in
  
  png(fname, width = width_in, height = height_in, units = "in", res = dpi)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(1, "npc") - grid::unit(2 * pad_in, "in"),
                                    height = grid::unit(1, "npc") - grid::unit(2 * pad_in, "in")))
  grid::grid.draw(legend_grob)
  grid::popViewport()
  dev.off()
  cat("saved", fname, sprintf("(%.2f x %.2f in)\n", width_in, height_in))
}

## ---------------- generate plots ----------------
measures_to_plot <- c(
  "Locomotion"          = "Locomotion (sec)",
  "Object Manipulation" = "Manipulation (sec)",
  "Passive Alone"       = "Passive Alone (sec)",
  "SelfDirected"        = "Self-Directed (sec)",
  "SocialContact"       = "Social Contact (sec)",
  "SocBoutCount"        = "Social Bouts (count)",
  "SocBoutDur"          = "Social Bout Duration, median (sec)",
  "Grooming"            = "Grooming (sec)",
  "Passive in Cont."    = "Passive in Contact (sec)",
  "Proximity"           = "Time in Proximity (sec)",
  "Alone"               = "Time Alone (sec)",
  "Passive"             = "Passive Total (sec)"
)

dir.create(out_dir, showWarnings = FALSE)

legend_source_plot <- make_plot(names(measures_to_plot)[1], measures_to_plot[[1]],
                                show_legend = TRUE)
save_legend_png(extract_legend(legend_source_plot), "estimationplot_legend_SiteBMI.png")

panel_gtables <- list()
for (i in seq_along(measures_to_plot)) {
  col        <- names(measures_to_plot)[i]
  axis_label <- measures_to_plot[[i]]

  plt <- make_plot(col, axis_label, show_legend = FALSE)
  # fixed panel size regardless of how wide this measure's axis labels are to make for easier montaging later
  panel_gtables[[col]] <- set_panel_size(plt,
                                       width  = unit(STYLE$panel_width_in, "in"),
                                       height = unit(STYLE$panel_height_in, "in"))
}

aligned_gtables <- align_plots(plotlist = panel_gtables, align = "hv", axis = "tblr")
names(aligned_gtables) <- names(panel_gtables)

for (i in seq_along(measures_to_plot)) {
  col        <- names(measures_to_plot)[i]
  axis_label <- measures_to_plot[[i]]

  g <- aligned_gtables[[col]]
  fname <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9]", "_", axis_label), ".png"))
  
  total_width_in  <- grid::convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
  total_height_in <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  
  ggsave(fname, g, width = total_width_in, height = total_height_in,
         dpi = STYLE$fig_dpi, units = "in")
  
  cat("saved", fname, sprintf("(%.2f x %.2f in)\n", total_width_in, total_height_in))
}