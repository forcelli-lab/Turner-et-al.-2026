## =====================================================================
## Exact nested permutation test for Turner et al.
##
##   Animal x Hemisphere = randomization block ("unit").
##   Becuase the number of sessions for each treatment varies between
##   hemispheres, the number of comparison-level sessions is kept contstant
##   during the shuffles.
##
##   Note - sessions are averaged within a hemisphere, and then hemispheres
##   are averaged within animal. Effects are then averaged across animals.
##
##   Input file: Column for Animal, Session, Hemisphere, Treatment, and 
##   the dependent variables. All should be specified in the config section
##
##   Statistic: mean over animals of
##     (mean over that animal's hemispheres of (BMI mean - Saline mean)).
##
##   Output: P values are exact and two sided.
##   Dependent measures: uncorrected p values reported
##
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(randomizr)
})

## ---- CONFIG ----
## Each dataset is analyzed independently; a combined table is also written.
data_dir <- "01_measures"

datasets <- c(
  ACCs = "ACCs_BMI.csv",
  ACCg = "ACCg_BMI.csv"
)

out_dir  <- "03_permutation"

subject_col      <- "Animal"
hemisphere_col   <- "Hemisphere"
treatment_col    <- "Treatment"
reference_level  <- "Saline"
comparison_level <- "BMI"

measure_cols <- c(
  "Locomotion", "Object Manipulation", "Passive Alone", "SelfDirected",
  "SocialContact", "SocBoutCount", "SocBoutDur",
  "Grooming", "Passive in Cont.", "Proximity", "Alone", "Passive"
)

max_exact_assignments <- 1e6 # this is to stop things from blowing up if the space is big.
#this version of the script is all exaact (no monte carlo) so, if too big, stop.
comparison_tolerance  <- sqrt(.Machine$double.eps) 

combined_results_file <- "combined.csv"

## ---- ONE DATASET ----

analyze_dataset <- function(data_path) {

  ## ---- LOAD & CHECK COLUMNS ----
  ## Read the sheet and confirm the expected columns and treatment labels are there.
  df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)

  missing_cols <- setdiff(c(subject_col, hemisphere_col, treatment_col, measure_cols), names(df))
  if (length(missing_cols) > 0L)
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))

  if (!setequal(unique(as.character(df[[treatment_col]])), c(reference_level, comparison_level)))
    stop("Treatment column must contain exactly ", reference_level, " and ", comparison_level, ".")

  ## ---- NA CHECK ----
  ## Checks for missing values.
  missing_counts <- sapply(measure_cols, function(m) sum(is.na(df[[m]])))
  if (any(missing_counts > 0)) {
    bad_measures <- names(missing_counts)[missing_counts > 0]
    stop("Missing values detected in: ", paste(bad_measures, collapse = ", "), 
         ". Matrix multiplication will fail. Drop these rows before running.")
  }

  ## ---- UNITS ----
  ## Build Animal x Hemisphere units and the row/subject lookups; require both levels per unit.
  df$unit <- interaction(df[[subject_col]], df[[hemisphere_col]], drop = TRUE, lex.order = TRUE)

  units    <- sort(unique(as.character(df$unit))) #sort these to keep things aligned for randmoization
  subjects <- unique(as.character(df[[subject_col]]))

  unit_rows     <- setNames(lapply(units, function(u) which(as.character(df$unit) == u)), units)
  subject_units <- setNames(
    lapply(subjects, function(s) unique(as.character(df$unit[as.character(df[[subject_col]]) == s]))),
    subjects
  )

  unit_counts <- table(
    factor(as.character(df$unit), levels = units),
    factor(as.character(df[[treatment_col]]), levels = c(reference_level, comparison_level))
  )

  if (any(unit_counts == 0L))
    stop("Every Animal x Hemisphere unit must contain both treatment levels.")

  ## ---- GENERATE ALL POSSIBLE SHUFFLES ----
  ## Specify all possible within-unit relabeling via randomizr.
  design        <- declare_ra(blocks = as.character(df$unit),
                              block_m = as.numeric(unit_counts[, comparison_level])) #presorted from above
  n_assignments <- obtain_num_permutations(design)

  if (!is.finite(n_assignments) || n_assignments > max_exact_assignments)
    stop("Exact space (", format(n_assignments, big.mark = ","),
         ") exceeds max_exact_assignments.")

  assignment_matrix <- obtain_permutation_matrix(design, maximum_permutations = n_assignments)
  if (ncol(assignment_matrix) != n_assignments)
    stop("randomizr did not return the full exact enumeration.")

  ## ---- PERMUTATION-WEIGHT MATRIX ----
  ## Turn each assignment into weights so weights %*% values gives the nested statistic (1 = BMI).
  permutation_weights <- matrix(0, nrow = n_assignments, ncol = nrow(df))

  for (s in subjects) {
    unit_weight <- 1 / length(subjects) / length(subject_units[[s]])
    for (u in subject_units[[s]]) {
      idx <- unit_rows[[u]]
      a   <- assignment_matrix[idx, , drop = FALSE]
      permutation_weights[, idx] <- t(ifelse(
        a == 1L,
        unit_weight / unit_counts[u, comparison_level],
        -unit_weight / unit_counts[u, reference_level]
      ))
    }
  }

  ## ---- TEST STATISTIC ----
  ## The real treatment labels, plus the direct statistic used for the observed value.
  observed_assignment <- as.integer(as.character(df[[treatment_col]]) == comparison_level)

  nested_statistic <- function(values, assignment) {
    unit_diff <- vapply(units, function(u) {
      idx <- unit_rows[[u]]
      mean(values[idx][assignment[idx] == 1L]) - mean(values[idx][assignment[idx] == 0L])
    }, numeric(1))
    mean(vapply(subjects, function(s) mean(unit_diff[subject_units[[s]]]), numeric(1)))
  }

  ## ---- RUN ----
  ## Observed statistic and two-sided exact p per measure.
  run_one <- function(measure) {
    v        <- df[[measure]]
    observed <- nested_statistic(v, observed_assignment)
    null     <- as.vector(permutation_weights %*% v)
    tol      <- comparison_tolerance * max(1, abs(observed))
    data.frame(
      Measure = measure, Observed = observed,
      P = mean(abs(null) >= abs(observed) - tol),
      N_arrangements = n_assignments, stringsAsFactors = FALSE
    )
  }

  do.call(rbind, lapply(measure_cols, run_one))
}

## ---- RUN ALL DATASETS ----

dir.create(out_dir, showWarnings = FALSE)

all_results <- lapply(names(datasets), function(tag) {
  cat("\n==== ", tag, " (", datasets[[tag]], ") ====\n", sep = "")

  results <- analyze_dataset(file.path(data_dir, datasets[[tag]]))
  results <- results[order(results$P), ]
  rownames(results) <- NULL

  print(results, row.names = FALSE, digits = 6)

  out_file <- file.path(out_dir, paste0(tag, ".csv"))
  write.csv(results, out_file, row.names = FALSE)
  cat("Saved:", out_file, "\n")

  cbind(Dataset = tag, results)
})

combined <- do.call(rbind, all_results)
write.csv(combined, file.path(out_dir, combined_results_file), row.names = FALSE)
cat("\nSaved combined:", file.path(out_dir, combined_results_file), "\n")
