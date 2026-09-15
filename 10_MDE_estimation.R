## ============================================================================
## Minimum Detectable Effect
## nested permutation test + simulation-based power
##
## The analysis:
##   1. Resamples pooled Saline and BMI observations within each
##      Animal x Hemisphere unit to generate null data sets.
##   2. Modulates BMI values across a range of candidate effect sizes.
##   3. Applies the permutation test to each simulated data set.
##   4. Finds the smallest effect with a lower 95% Wilson power bound >= 0.80.
##
## Input file:
##  an .xlsx file that contains Animal, hemipshere, session, and treatment columns
##  as well as dependent measure columns
##
## Study statistic
##   1. BMI mean - Saline mean within each hemisphere.
##   2. Average (equal weight) hemispheres within animal.
##   3. Average (equal-weight) animals.
##
## Required packages
##   install.packages(c("dplyr", "DeclareDesign", "randomizr", "dqrng"))
##
## For controlling multithreading
##   install.packages("RhpcBLASctl")
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(DeclareDesign)
  library(randomizr)
  library(dqrng)
  library(parallel)
})

## ---- CONFIGURATION ----

data_dir         <- "01_measures"

## Each dataset writes its own results table and progress log.
datasets         <- c(ACCs = "ACCs_BMI.csv",
                      ACCg = "ACCg_BMI.csv")

out_dir          <- "10_MDE"

subject_col    <- "Animal"
hemisphere_col <- "Hemisphere"
session_col    <- "Session"
treatment_col  <- "Treatment"

reference_level  <- "Saline"
comparison_level <- "BMI"

#update if data format changes
measure_cols <- c(
  "Locomotion", "Object Manipulation", "Passive Alone", "SelfDirected",
  "SocialContact", "SocBoutCount", "SocBoutDur",
  "Grooming", "Passive in Cont.", "Proximity", "Alone", "Passive"
)

integer_measure_cols <- "SocBoutCount"

duration_measure_cols <- setdiff(
  measure_cols,
  integer_measure_cols
)

duration_ceiling <- 1800 #sets ceiling for duration-based behaviors based on the total test duration

## Statistical targets
alpha        <- 0.05
target_power <- 0.80

## Lower-confidence-bound criterion
## lower 95% Wilson confidence bound >= .80.
confirm <- list(
  max_sims = 50000L,
  min_sims = 20000L,
  step = 0.25,
  precision = 0.001
)

maximum_increase_pct <- 1000
confidence_level <- 0.95

## Computation
parallel_workers <- 6L
simulation_batch_size <- 250L
permutation_row_chunk_size <- 50000L
progress_every_sims <- 1000L
max_exact_assignments <- 1e6
random_seed <- 20260718L
comparison_tolerance <- sqrt(.Machine$double.eps)

## Output
dir.create(out_dir, showWarnings = FALSE)

## One BLAS thread inside each parallel worker.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(if (parallel_workers > 1L) 1L else 5L)
}

## Live progress logging. Forked workers write to the same file, so the path and
## the print position live in one environment rather than in globals.
log_state <- new.env(parent = emptyenv())

start_log <- function(path) {
  log_state$file          <- path
  log_state$lines_printed <- 0L
  if (file.exists(path)) file.remove(path)
  file.create(path)
  invisible(path)
}

status <- function(..., measure = NULL) {
  prefix <- sprintf("[%s] [pid %d]", format(Sys.time(), "%H:%M:%S"), Sys.getpid())
  if (!is.null(measure)) prefix <- paste0(prefix, " [", measure, "]")
  line <- paste0(prefix, " ", paste0(..., collapse = ""))
  
  con <- file(log_state$file, open = "a")
  on.exit(close(con), add = TRUE)
  writeLines(line, con = con, sep = "\n", useBytes = TRUE)
  flush(con)
  invisible(line)
}

print_new_status <- function(final = FALSE) {
  if (!file.exists(log_state$file)) return(invisible(NULL))
  lines <- readLines(log_state$file, warn = FALSE)
  if (length(lines) > log_state$lines_printed) {
    new_lines <- lines[(log_state$lines_printed + 1L):length(lines)]
    cat(paste0(new_lines, collapse = "\n"), "\n", sep = "")
    flush.console()
    log_state$lines_printed <- length(lines)
  }
  if (final) flush.console()
  invisible(NULL)
}

run_fork_pool <- function(indices, worker_fun, max_workers, poll_seconds = 0.25) {
  results <- vector("list", length(indices))
  names(results) <- as.character(indices)
  pending <- as.list(indices)
  running <- list()
  pid_to_index <- integer(0)
  
  launch_one <- function(index) {
    job <- parallel::mcparallel(
      worker_fun(index),
      mc.set.seed = FALSE,
      silent = TRUE
    )
    running[[as.character(job$pid)]] <<- job
    pid_to_index[as.character(job$pid)] <<- index
  }
  
  while (length(pending) > 0L && length(running) < max_workers) {
    launch_one(pending[[1L]])
    pending <- pending[-1L]
  }
  
  while (length(running) > 0L) {
    print_new_status()
    
    collected <- parallel::mccollect(running, wait = FALSE)
    if (!is.null(collected)) {
      for (pid in names(collected)) {
        index <- unname(pid_to_index[[pid]])
        results[[as.character(index)]] <- collected[[pid]]
        running[[pid]] <- NULL
        pid_to_index <- pid_to_index[names(pid_to_index) != pid]
        
        if (length(pending) > 0L) {
          launch_one(pending[[1L]])
          pending <- pending[-1L]
        }
      }
    }
    
    if (length(running) > 0L) Sys.sleep(poll_seconds)
  }
  
  print_new_status(final = TRUE)
  unname(results)
}

## ============================================================================
## ONE DATASET
## ============================================================================

run_mde <- function(tag, file) {
  
  data_path    <- file.path(data_dir, file)
  output_stub  <- file.path(out_dir, tag)
  results_file <- paste0(output_stub, "_results.csv")
  
  start_log(paste0(output_stub, "_progress.log"))
  
  ## ============================================================================
  ## 1. LOAD AND VALIDATE DATA
  ## ============================================================================
  
  df <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  required_cols <- c(
    subject_col,
    hemisphere_col,
    treatment_col,
    measure_cols,
    if (!is.null(session_col) && nzchar(session_col)) session_col
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  for (col in c(subject_col, hemisphere_col, treatment_col)) {
    if (is.character(df[[col]])) df[[col]] <- trimws(df[[col]])
  }
  
  if (anyNA(df[required_cols])) {
    missing_counts <- vapply(df[required_cols], function(x) sum(is.na(x)), integer(1))
    bad <- names(missing_counts)[missing_counts > 0L]
    stop(
      "Missing values found in: ",
      paste(sprintf("%s (%d)", bad, missing_counts[bad]), collapse = ", ")
    )
  }
  
  if (!setequal(unique(df[[treatment_col]]), c(reference_level, comparison_level))) {
    stop(
      "Treatment labels must be exactly '", reference_level,
      "' and '", comparison_level, "'."
    )
  }
  
  non_numeric <- measure_cols[!vapply(df[measure_cols], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric outcome column(s): ", paste(non_numeric, collapse = ", "))
  }
  
  invalid_outcomes <- vapply(
    df[measure_cols],
    function(x) sum(!is.finite(x) | x < 0),
    integer(1)
  )
  if (any(invalid_outcomes > 0L)) {
    bad <- names(invalid_outcomes)[invalid_outcomes > 0L]
    stop(
      "Non-finite or negative outcome values in: ",
      paste(sprintf("%s (%d)", bad, invalid_outcomes[bad]), collapse = ", ")
    )
  }
  
  for (col in duration_measure_cols) {
    if (any(df[[col]] > duration_ceiling + 1e-10)) {
      stop(col, " contains values above ", duration_ceiling, " seconds.")
    }
  }
  
  for (col in integer_measure_cols) {
    if (any(abs(df[[col]] - round(df[[col]])) > 1e-8)) {
      stop(col, " is configured as a count but contains noninteger values.")
    }
  }
  
  if (!is.null(session_col) && nzchar(session_col)) {
    duplicates <- df %>%
      count(across(all_of(c(subject_col, hemisphere_col, session_col))), name = "n") %>%
      filter(n > 1L)
    
    if (nrow(duplicates)) {
      print(duplicates)
      stop("Duplicate Animal x Hemisphere x Session identifiers found.")
    }
  }
  
  ## ============================================================================
  ## Nested Design
  ## ============================================================================
  
  df$unit <- interaction(
    df[[subject_col]],
    df[[hemisphere_col]],
    drop = TRUE,
    lex.order = TRUE
  )
  
  units <- sort(unique(as.character(df$unit)))
  subjects <- unique(as.character(df[[subject_col]]))
  
  unit_rows <- setNames(
    lapply(units, function(u) which(as.character(df$unit) == u)),
    units
  )
  
  subject_units <- setNames(
    lapply(subjects, function(s) unique(as.character(df$unit[df[[subject_col]] == s]))),
    subjects
  )
  
  unit_treatment_counts <- table(
    factor(as.character(df$unit), levels = units),
    factor(df[[treatment_col]], levels = c(reference_level, comparison_level))
  )
  
  bad_units <- units[apply(unit_treatment_counts, 1L, function(x) any(x == 0L))]
  if (length(bad_units)) {
    stop("Units missing one treatment: ", paste(bad_units, collapse = ", "))
  }
  
  comparison_rows <- df[[treatment_col]] == comparison_level
  
  ## ============================================================================
  ## Generate all possible permutations of labels
  ## ============================================================================
  
  ## randomizr treats BMI as treatment = 1 and Saline as control = 0.
  assignment_declaration <- declare_ra(
    blocks = as.character(df$unit),
    block_m = as.numeric(unit_treatment_counts[, comparison_level])
  )
  
  n_permutations <- obtain_num_permutations(assignment_declaration)
  
  if (!is.finite(n_permutations) || n_permutations > max_exact_assignments) {
    stop(
      "Exact assignment space is ", format(n_permutations, big.mark = ","),
      "; max_exact_assignments is ", format(max_exact_assignments, big.mark = ","), "."
    )
  }
  
  assignment_matrix <- obtain_permutation_matrix(
    assignment_declaration,
    maximum_permutations = n_permutations
  )
  
  if (ncol(assignment_matrix) != n_permutations) {
    stop("randomizr sampled assignments instead of returning the full exact space.")
  }
  
  cat(
    "Exact legal assignments:", format(n_permutations, big.mark = ","), "\n",
    "Theoretical p-value floor:", format(1 / n_permutations, scientific = TRUE), "\n"
  )
  
  ## ============================================================================
  ## Test Statistic
  ## ============================================================================
  
  permutation_weights <- matrix(0, nrow = n_permutations, ncol = nrow(df))
  observed_weights <- numeric(nrow(df))
  
  for (subject in subjects) {
    subject_unit_names <- subject_units[[subject]]
    unit_weight <- 1 / length(subjects) / length(subject_unit_names)
    
    for (u in subject_unit_names) {
      i <- unit_rows[[u]]
      n_ref <- unit_treatment_counts[u, reference_level]
      n_cmp <- unit_treatment_counts[u, comparison_level]
      
      z <- assignment_matrix[i, , drop = FALSE]
      permutation_weights[, i] <- t(ifelse(z == 1L, unit_weight / n_cmp, -unit_weight / n_ref))
      
      observed_weights[i[df[[treatment_col]][i] == reference_level]] <- -unit_weight / n_ref
      observed_weights[i[df[[treatment_col]][i] == comparison_level]] <-  unit_weight / n_cmp
    }
  }
  
  rm(assignment_matrix)
  gc(verbose = FALSE)
  
  ## ============================================================================
  ## SUPPORT FUNCTIONS
  ## ============================================================================
  
  make_seed <- function(measure_index, direction, percent_effect = 0) {
    direction_code <- match(direction, c("increase", "decrease"))
    
    value <- random_seed +
      measure_index * 1000003 +
      direction_code * 10007 +
      round(percent_effect * 100) * 17
    
    as.integer(value %% (.Machine$integer.max - 1L) + 1L)
  }
  
  wilson_interval <- function(successes, trials, confidence = confidence_level) {
    p <- successes / trials
    z <- qnorm(1 - (1 - confidence) / 2)
    denominator <- 1 + z^2 / trials
    center <- (p + z^2 / (2 * trials)) / denominator
    half_width <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) / denominator
    
    c(lower = max(0, center - half_width), upper = min(1, center + half_width))
  }
  
  stochastic_round <- function(x) {
    lower <- floor(x)
    lower + matrix(
      rbinom(length(x), 1L, as.vector(x - lower)),
      nrow = nrow(x),
      ncol = ncol(x)
    )
  }
  
  ## Pooled empirical-null distributions, cached once for every measure.
  null_pools <- setNames(lapply(measure_cols, function(measure) {
    setNames(lapply(units, function(u) df[[measure]][unit_rows[[u]]]), units)
  }), measure_cols)
  
  ## ============================================================================
  ## GENERATE EMPIRICAL NULL DATA SETS
  ## ============================================================================
  
  ## Draw once for each outcome and direction. All candidate effect sizes
  ## within that search use the same simulated null datasets
  draw_base_null <- function(measure, n_sims, seed) {
    dqrng::dqset.seed(seed)
    
    values <- matrix(NA_real_, nrow = nrow(df), ncol = n_sims)
    
    for (u in units) {
      i <- unit_rows[[u]]
      
      ## Pool saline and BMI observations within this hemisphere and resample
      ## from that distribution. This breaks any treatment
      ## effect while preserving the underlying data structure.
      values[i, ] <- matrix(
        dqrng::dqsample(
          null_pools[[measure]][[u]],
          length(i) * n_sims,
          replace = TRUE
        ),
        nrow = length(i),
        ncol = n_sims
      )
    }
    
    values
  }
  
  ## Apply one requested effect to a shared null matrix. Sequential batches
  ## then take column slices from this transformed matrix.
  apply_effect <- function(base_null, measure, percent_effect, direction, rounding_seed) {
    values <- base_null
    before <- values[comparison_rows, , drop = FALSE]
    
    ## Effects are imposed only on rows assigned to BMI. The requested effect is
    ## multiplicative: e.g., 25% increase multiplies BMI values by 1.25, whereas
    ## 25% decrease multiplies them by 0.75.
    multiplier <- if (direction == "increase") {
      1 + percent_effect / 100
    } else {
      1 - percent_effect / 100
    }
    
    after <- before * multiplier
    
    ## Duration outcomes cannot exceed the 30-minute observation ceiling.
    ## Because this can make the achieved effect smaller than the
    ## requested multiplier, achieved effects are tracked separately below.
    if (measure %in% duration_measure_cols && direction == "increase") {
      after <- pmin(after, duration_ceiling)
    }
    
    values[comparison_rows, ] <- after
    
    if (measure %in% integer_measure_cols) {
      set.seed(rounding_seed)
      values <- stochastic_round(values)
      after <- values[comparison_rows, , drop = FALSE]
    }
    
    before_mean <- colMeans(before)
    after_mean <- colMeans(after)
    achieved <- rep(NA_real_, ncol(values))
    positive <- before_mean > 0
    achieved[positive] <- 100 * (after_mean[positive] - before_mean[positive]) / before_mean[positive]
    achieved[before_mean == 0 & after_mean == 0] <- 0
    
    list(
      outcomes = values,
      achieved_effect = achieved
    )
  }
  
  ## DeclareDesign model handler: expose one sequential column block from the
  ## already-generated effect data.
  serve_effect_batch <- function(effect_data, first_col, last_col) {
    columns <- first_col:last_col
    
    data.frame(
      outcomes = I(list(effect_data$outcomes[, columns, drop = FALSE])),
      achieved_effect = I(list(effect_data$achieved_effect[columns])),
      stringsAsFactors = FALSE
    )
  }
  
  ## ============================================================================
  ## PERMUTATION TEST
  ## ============================================================================
  
  classify_exact_batch <- function(simulated) {
    ## Each column of simulated is one complete simulated data set. For each one,
    ## calculate the observed nested statistic and compare it with the statistic
    ## under every possible assignment of labels. 
    observed <- as.vector(crossprod(observed_weights, simulated))
    threshold <- abs(observed) - comparison_tolerance * pmax(1, abs(observed))
    rejection_limit <- floor(alpha * n_permutations)
    
    extreme <- integer(ncol(simulated))
    decided <- logical(ncol(simulated))
    rejected <- logical(ncol(simulated))
    
    for (start in seq(1L, n_permutations, by = permutation_row_chunk_size)) {
      active <- which(!decided)
      if (!length(active)) break
      
      end <- min(start + permutation_row_chunk_size - 1L, n_permutations)
      rows <- start:end
      
      permuted <- permutation_weights[rows, , drop = FALSE] %*%
        simulated[, active, drop = FALSE]
      
      ## Vectorized column-major comparison. Each threshold is repeated for one
      ## complete permutation-statistic column.
      threshold_vector <- rep(threshold[active], each = nrow(permuted))
      extreme[active] <- extreme[active] +
        colSums(abs(permuted) >= threshold_vector)
      
      remaining <- n_permutations - end
      
      no_reject <- active[extreme[active] > rejection_limit]
      decided[no_reject] <- TRUE
      
      active <- which(!decided)
      yes_reject <- active[extreme[active] + remaining <= rejection_limit]
      decided[yes_reject] <- TRUE
      rejected[yes_reject] <- TRUE
    }
    
    unresolved <- which(!decided)
    rejected[unresolved] <- extreme[unresolved] <= rejection_limit
    rejected
  }
  
  exact_batch_test <- function(data) {
    rejected <- classify_exact_batch(data$outcomes[[1]])
    achieved <- data$achieved_effect[[1]]
    
    data.frame(
      p.value = mean(rejected),
      n_rejected = sum(rejected),
      n_sim = length(rejected),
      achieved_sum = sum(achieved, na.rm = TRUE),
      achieved_n = sum(is.finite(achieved)),
      stringsAsFactors = FALSE
    )
  }
  
  ## ============================================================================
  ## SIMULATION DESIGN
  ## ============================================================================
  
  mde_designer <- function(effect_data, first_col, last_col) {
    declare_model(
      handler = serve_effect_batch,
      effect_data = effect_data,
      first_col = first_col,
      last_col = last_col,
      label = "Shared empirical-null model and effect"
    ) +
      declare_test(
        handler = label_test(exact_batch_test),
        label = "Exact nested permutation test"
      )
  }
  
  ## ============================================================================
  ## ESTIMATE POWER AT ONE EFFECT SIZE
  ## ============================================================================
  
  diagnose_effect <- function(
    measure, measure_index, percent_effect, direction, base_null
  ) {
    settings <- confirm
    
    effect_data <- apply_effect(
      base_null = base_null,
      measure = measure,
      percent_effect = percent_effect,
      direction = direction,
      rounding_seed = make_seed(measure_index, direction, percent_effect)
    )
    
    n <- rejected <- 0L
    achieved_sum <- achieved_n <- 0
    reason <- "maximum simulations reached"
    next_progress <- progress_every_sims
    
    status(
      direction, " ", percent_effect, "% started; max ",
      settings$max_sims, " simulations",
      measure = measure
    )
    
    ## Simulations are processed in batches. After the minimum number has been
    ## reached, the loop may stop early when the Wilson interval is above or below the target.
    while (n < settings$max_sims) {
      first_col <- n + 1L
      last_col <- min(n + simulation_batch_size, settings$max_sims)
      
      design <- mde_designer(effect_data, first_col, last_col)
      run <- run_design(design)
      
      n <- n + run$n_sim
      rejected <- rejected + run$n_rejected
      achieved_sum <- achieved_sum + run$achieved_sum
      achieved_n <- achieved_n + run$achieved_n
      
      if (n >= next_progress || n == settings$max_sims) {
        current_power <- rejected / n
        current_ci <- wilson_interval(rejected, n)
        status(
          direction, " ", percent_effect, "%: n=", n,
          ", power=", sprintf("%.3f", current_power),
          " [", sprintf("%.3f", current_ci["lower"]), ", ",
          sprintf("%.3f", current_ci["upper"]), "]",
          measure = measure
        )
        while (next_progress <= n) next_progress <- next_progress + progress_every_sims
      }
      
      if (n >= settings$min_sims) {
        ci <- wilson_interval(rejected, n)
        half_width <- diff(ci) / 2
        
        if (ci["lower"] >= target_power) {
          reason <- "lower bound above target"
          break
        }
        if (ci["upper"] < target_power) {
          reason <- "upper bound below target"
          break
        }
        if (half_width <= settings$precision) {
          reason <- "precision target reached"
          break
        }
      }
    }
    
    power <- rejected / n
    ci <- wilson_interval(rejected, n)
    
    status(
      direction, " ", percent_effect, "% finished: power=",
      sprintf("%.3f", power), " [", sprintf("%.3f", ci["lower"]), ", ",
      sprintf("%.3f", ci["upper"]), "], n=", n, ", ", reason,
      measure = measure
    )
    
    data.frame(
      Direction = direction,
      Requested_Percent_Effect = percent_effect,
      Mean_Achieved_Percent_Effect =
        if (achieved_n) achieved_sum / achieved_n else NA_real_,
      Power_Raw = power,
      Power_CI_Lower = ci["lower"],
      Power_CI_Upper = ci["upper"],
      N_Simulations = n,
      row.names = NULL
    )
  }
  
  
  ## ============================================================================
  ## SEARCH FOR THE MDE
  ## ============================================================================
  
  search_direction <- function(measure, measure_index, direction) {
    direction_start <- proc.time()[["elapsed"]]
    status("starting ", direction, " MDE search (maximum screen -> bisection)",
           measure = measure)
    
    ## Search maximum: 100% for decreases, maximum_increase_pct for increases.
    maximum_effect <- if (direction == "increase") maximum_increase_pct else 100
    
    ## every candidate in this direction is evaluated on
    ## the same empirical-null simulations.
    confirm_base_null <- draw_base_null(
      measure,
      confirm$max_sims,
      make_seed(measure_index, direction, 0)
    )
    
    evaluated <- NULL
    
    evaluate_once <- function(pct) {
      existing <- which(
        !is.null(evaluated) &
          abs(evaluated$Requested_Percent_Effect - pct) < 1e-10
      )
      if (length(existing)) return(evaluated[existing[1], , drop = FALSE])
      
      estimate <- diagnose_effect(
        measure, measure_index, pct, direction, confirm_base_null
      )
      evaluated <<- bind_rows(evaluated, estimate) %>%
        arrange(Requested_Percent_Effect)
      estimate
    }
    
    is_confirmed <- function(x) {
      is.finite(x$Power_CI_Lower) && x$Power_CI_Lower >= target_power
    }
    
    ## Screen the maximum effect first; if it fails, no smaller effect can pass.
    status(
      "screening maximum allowable ", direction, " effect (",
      maximum_effect, "%)",
      measure = measure
    )
    
    maximum_result <- evaluate_once(maximum_effect)
    
    if (!is_confirmed(maximum_result)) {
      if (maximum_result$Power_CI_Upper < target_power) {
        result_status <- paste0(
          "No confirmed ", direction,
          " MDE: even the maximum allowable effect has an upper 95% power ",
          "bound below the target."
        )
      } else {
        result_status <- paste0(
          "No confirmed ", direction,
          " MDE: the maximum allowable effect did not achieve a lower 95% ",
          "power bound at or above the target."
        )
      }
      
      status(result_status, measure = measure)
      return(list(
        mde = NA_real_,
        confirmed = NULL
      ))
    }
    
    ## Bisect to the first confirmed crossing.
    low <- 0
    high <- maximum_effect
    
    status(
      "maximum passed; beginning bisection in [", low, ", ", high, "]",
      measure = measure
    )
    
    while ((high - low) > confirm$step) {
      midpoint <- floor((low + high) / 2 / confirm$step) * confirm$step
      
      ## Protect against a midpoint equal to a bound because of integer rounding.
      if (midpoint <= low) midpoint <- low + confirm$step
      if (midpoint >= high) midpoint <- high - confirm$step
      
      midpoint_result <- evaluate_once(midpoint)
      
      if (is_confirmed(midpoint_result)) {
        high <- midpoint
        decision <- "confirmed; lowering upper bound"
      } else {
        low <- midpoint
        decision <- "not confirmed; raising lower bound"
      }
      
      status(
        "bisection ", direction, ": ", midpoint, "% ", decision,
        "; bracket now [", low, ", ", high, "]",
        measure = measure
      )
    }
    
    provisional <- high
    
    ## Require two consecutive confirmed integer-percent effects.
    ##
    ## Bisection guarantees that provisional is confirmed. Evaluate upward only
    ## ------------------------------------------------------------------------
    candidate <- provisional
    confirmed <- NULL
    
    while (candidate <= maximum_effect) {
      first <- evaluate_once(candidate)
      second_pct <- candidate + confirm$step
      
      if (second_pct > maximum_effect) break
      
      second <- evaluate_once(second_pct)
      
      if (is_confirmed(first) && is_confirmed(second)) {
        confirmed <- first
        break
      }
      
      candidate <- candidate + confirm$step
    }
    
    if (is.null(confirmed)) {
      result_status <- paste0(
        "No strict confirmed crossing: two consecutive effects did not meet ",
        "the lower-95%-confidence-bound criterion within the search range."
      )
      
      status(result_status, measure = measure)
      return(list(
        mde = NA_real_,
        confirmed = NULL
      ))
    }
    
    elapsed <- proc.time()[["elapsed"]] - direction_start
    result_status <-
      "Confirmed by two consecutive lower 95% confidence bounds after bisection."
    
    status(
      direction, " MDE confirmed at ",
      confirmed$Requested_Percent_Effect, "% in ",
      sprintf("%.1f", elapsed), " seconds",
      measure = measure
    )
    
    list(
      mde = confirmed$Requested_Percent_Effect,
      confirmed = confirmed
    )
  }
  
  ## ============================================================================
  ## SUMMARIZE ONE OUTCOME
  ## ============================================================================
  
  nested_reference_mean <- function(measure) {
    unit_mean <- vapply(units, function(u) {
      i <- unit_rows[[u]]
      mean(df[[measure]][i][df[[treatment_col]][i] == reference_level])
    }, numeric(1))
    names(unit_mean) <- units
    
    mean(vapply(subjects, function(s) mean(unit_mean[subject_units[[s]]]), numeric(1)))
  }
  
  value_or_na <- function(row, column) {
    if (is.null(row)) NA else row[[column]][1]
  }
  
  analyze_measure <- function(measure_index) {
    measure <- measure_cols[measure_index]
    measure_start <- proc.time()[["elapsed"]]
    status("outcome started (", measure_index, "/", length(measure_cols), ")", measure = measure)
    
    increase <- search_direction(measure, measure_index, "increase")
    decrease <- search_direction(measure, measure_index, "decrease")
    reference_mean <- nested_reference_mean(measure)
    
    achieved_increase <- value_or_na(increase$confirmed, "Mean_Achieved_Percent_Effect")
    achieved_decrease <- value_or_na(decrease$confirmed, "Mean_Achieved_Percent_Effect")
    
    selected <- data.frame(
      Measure = measure,
      Saline_Mean = reference_mean,
      
      Increase_Achieved_Pct = achieved_increase,
      Increase_Raw_Change = reference_mean * achieved_increase / 100,
      Increase_Power = value_or_na(increase$confirmed, "Power_Raw"),
      Increase_CI_Lower = value_or_na(increase$confirmed, "Power_CI_Lower"),
      Increase_CI_Upper = value_or_na(increase$confirmed, "Power_CI_Upper"),
      Increase_N = value_or_na(increase$confirmed, "N_Simulations"),
      
      Decrease_Achieved_Pct = abs(achieved_decrease),
      Decrease_Raw_Change = reference_mean * abs(achieved_decrease) / 100,
      Decrease_Power = value_or_na(decrease$confirmed, "Power_Raw"),
      Decrease_CI_Lower = value_or_na(decrease$confirmed, "Power_CI_Lower"),
      Decrease_CI_Upper = value_or_na(decrease$confirmed, "Power_CI_Upper"),
      Decrease_N = value_or_na(decrease$confirmed, "N_Simulations"),
      
      stringsAsFactors = FALSE
    )
    
    elapsed_minutes <- (proc.time()[["elapsed"]] - measure_start) / 60
    status(
      "outcome finished in ", sprintf("%.1f", elapsed_minutes), " min; ",
      "increase MDE=", ifelse(is.na(increase$mde), "not found", paste0(increase$mde, "%")),
      ", decrease MDE=", ifelse(is.na(decrease$mde), "not found", paste0(decrease$mde, "%")),
      measure = measure
    )
    
    selected
  }
  
  safe_analyze_measure <- function(i) {
    tryCatch(
      list(ok = TRUE, measure = measure_cols[i], result = analyze_measure(i), error = NA_character_),
      error = function(e) {
        status("ERROR: ", conditionMessage(e), measure = measure_cols[i])
        list(
          ok = FALSE, measure = measure_cols[i], result = NULL,
          error = conditionMessage(e)
        )
      }
    )
  }
  
  
  ## ============================================================================
  ## RUN ANALYSIS AND SAVE RESULTS
  ## ============================================================================
  
  worker_count <- min(parallel_workers, length(measure_cols))
  status(
    "analysis started: ", length(measure_cols), " outcomes, ", worker_count,
    " workers, progress every ", progress_every_sims, " simulations"
  )
  
  print_new_status()
  
  if (.Platform$OS.type == "unix" && worker_count > 1L) {
    wrapped <- run_fork_pool(
      indices = seq_along(measure_cols),
      worker_fun = safe_analyze_measure,
      max_workers = worker_count
    )
  } else {
    wrapped <- vector("list", length(measure_cols))
    for (i in seq_along(measure_cols)) {
      wrapped[[i]] <- safe_analyze_measure(i)
      print_new_status()
    }
  }
  
  failed <- !vapply(wrapped, `[[`, logical(1), "ok")
  if (any(failed)) {
    failures <- data.frame(
      Measure = vapply(wrapped[failed], `[[`, character(1), "measure"),
      Error = vapply(wrapped[failed], `[[`, character(1), "error")
    )
    print(failures)
    stop(sum(failed), " outcome(s) failed.")
  }
  
  results <- lapply(wrapped, `[[`, "result")
  selected_results <- bind_rows(results)
  
  ## Results table.
  results_table <- selected_results %>%
    transmute(
      Measure,
      Saline_Mean = round(Saline_Mean, 1),
      
      Increase_MDE_Pct = round(Increase_Achieved_Pct, 1),
      Increase_Raw_Change = round(Increase_Raw_Change, 1),
      Increase_Power = round(Increase_Power, 3),
      Increase_Power_CI = ifelse(
        is.na(Increase_CI_Lower),
        NA_character_,
        sprintf("%.3f-%.3f", Increase_CI_Lower, Increase_CI_Upper)
      ),
      Increase_Simulations = as.integer(Increase_N),
      
      Decrease_MDE_Pct = round(Decrease_Achieved_Pct, 1),
      Decrease_Raw_Change = round(Decrease_Raw_Change, 1),
      Decrease_Power = round(Decrease_Power, 3),
      Decrease_Power_CI = ifelse(
        is.na(Decrease_CI_Lower),
        NA_character_,
        sprintf("%.3f-%.3f", Decrease_CI_Lower, Decrease_CI_Upper)
      ),
      Decrease_Simulations = as.integer(Decrease_N)
    ) %>%
    arrange(Increase_MDE_Pct)
  
  write.csv(results_table, results_file, row.names = FALSE)
  
  status("analysis complete; writing final output")
  cat("\n=== MDE RESULTS ===\n")
  print(results_table, row.names = FALSE)
  
  cat(
    "\nPermutations evaluated within each simulated data set: ",
    format(n_permutations, big.mark = ","), "\n",
    "Saved results table: ", results_file, "\n",
    "Progress log:        ", log_state$file, "\n",
    sep = ""
  )
  invisible(results_table)
}

## ============================================================================
## RUN ALL DATASETS
## ============================================================================

for (tag in names(datasets)) {
  cat("\n########  ", tag, "  ########\n", sep = "")
  run_mde(tag, datasets[[tag]])
}