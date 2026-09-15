# ============================================================
# Exact nested permutation test: BMI versus SAL
#
# Outcomes:
#   1. Pro-social transition probability
#   2. Social outcome following Approach
#
# Statistic:
#   - Pool counts across sessions within Animal x Hemisphere x Treatment.
#   - Compute BMI - SAL within each retained hemisphere.
#   - Average retained hemispheres equally within animal.
#   - Average retained animals equally.
#
# Zero denominators:
#   For each outcome, a hemisphere is retained only when BOTH treatment
#   denominators are positive under EVERY legal within-hemisphere assignment.
#   The retained hemisphere set is therefore fixed across the observed and
#   permuted statistics.
#
# Multiplicity:
#   Raw two-sided exact p-values and BH FDR across the two outcomes
#   within each region.
#
# Input: 
# #  2 files (ACCs_transitions.xlsx, ACCg_transitions.xlsx) which are output by the behavior transition
#  processing script. Animal, hemisphere and treatment are joined from
#  MANIFEST.csv by that script; no manual annotation step.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(writexl)
  library(randomizr)
})

reference_level <- "SAL"
comparison_level <- "BMI"
assignment_chunk_size <- 50000L
comparison_tolerance <- 1e-12

# written by 02_ACC_behavioral_batch_analysis.R
data_dir <- "02_transitions"
out_dir  <- "05_transition_BMI_vs_SAL"

expected_files <- c(
  "ACCs_transitions.xlsx",
  "ACCg_transitions.xlsx"
)

outcome_definitions <- tribble(
  ~outcome, ~numerator, ~denominator,
  "Pro-social behavioral transition",
  "n_pro_social_behavior_transitions",
  "n_from_non_social_behavior",
  
  "Social outcome following Approach",
  "n_approach_social_success",
  "n_classifiable_approaches"
)


# ---- Input -----------------------------------------------------------------

input_directory <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)

input_files <- file.path(
  input_directory,
  expected_files
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files)) {
  stop(
    "Missing required workbook(s): ",
    paste(
      basename(missing_files),
      collapse = ", "
    ),
    call. = FALSE
  )
}


# ---- Input helpers ---------------------------------------------------------

infer_region <- function(file) {
  sub(
    "_transitions\\.xlsx$",
    "",
    basename(file)
  )
}

normalize_treatment <- function(x) {
  value <- toupper(
    trimws(as.character(x))
  )
  
  case_when(
    value %in% c("BMI", "DRUG") ~
      comparison_level,
    value %in%
      c("SAL", "SALINE", "VEHICLE", "CONTROL") ~
      reference_level,
    TRUE ~ value
  )
}

read_region_data <- function(file) {
  message("Reading: ", basename(file))
  
  raw <- read_excel(
    file,
    sheet = "key_outputs"
  )
  
  required_columns <- unique(c(
    "input_file",
    "Animal",
    "Hemisphere",
    "Treatment",
    outcome_definitions$numerator,
    outcome_definitions$denominator
  ))
  
  missing_columns <- setdiff(
    required_columns,
    names(raw)
  )
  
  if (length(missing_columns)) {
    stop(
      basename(file),
      " key_outputs is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if ("processing_error" %in% names(raw)) {
    bad_rows <- raw |>
      filter(
        !is.na(processing_error),
        trimws(as.character(processing_error)) != ""
      )
    
    if (nrow(bad_rows)) {
      stop(
        basename(file),
        " contains ",
        nrow(bad_rows),
        " row(s) with processing errors.",
        call. = FALSE
      )
    }
  }
  
  if ("processing_status" %in% names(raw)) {
    bad_rows <- raw |>
      filter(
        !is.na(processing_status),
        trimws(as.character(processing_status)) != "",
        toupper(trimws(as.character(processing_status))) != "OK"
      )
    
    if (nrow(bad_rows)) {
      stop(
        basename(file),
        " contains ",
        nrow(bad_rows),
        " row(s) with non-OK processing status.",
        call. = FALSE
      )
    }
  }
  
  data <- raw |>
    transmute(
      Region = infer_region(file),
      input_file =
        as.character(input_file),
      Animal =
        toupper(trimws(as.character(Animal))),
      Hemisphere =
        tools::toTitleCase(
          trimws(as.character(Hemisphere))
        ),
      Treatment =
        normalize_treatment(Treatment),
      across(
        all_of(unique(c(
          outcome_definitions$numerator,
          outcome_definitions$denominator
        ))),
        as.numeric
      )
    ) |>
    mutate(
      unit = interaction(
        Animal,
        Hemisphere,
        drop = TRUE,
        sep = "__"
      ),
      session_row = row_number()
    ) |>
    arrange(
      Animal,
      Hemisphere,
      session_row
    )
  
  if (
    anyNA(data$Animal) ||
    any(data$Animal == "")
  ) {
    stop(
      basename(file),
      " contains missing or blank Animal values.",
      call. = FALSE
    )
  }
  
  if (
    anyNA(data$Hemisphere) ||
    any(data$Hemisphere == "")
  ) {
    stop(
      basename(file),
      " contains missing or blank Hemisphere values.",
      call. = FALSE
    )
  }
  
  invalid_treatments <- setdiff(
    unique(data$Treatment),
    c(
      reference_level,
      comparison_level
    )
  )
  
  if (length(invalid_treatments)) {
    stop(
      basename(file),
      " contains unrecognized Treatment value(s): ",
      paste(
        invalid_treatments,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  measure_columns <- unique(c(
    outcome_definitions$numerator,
    outcome_definitions$denominator
  ))
  
  missing_measure_values <- vapply(
    data[measure_columns],
    function(x) anyNA(x),
    logical(1)
  )
  
  if (any(missing_measure_values)) {
    stop(
      basename(file),
      " contains missing values in: ",
      paste(
        names(missing_measure_values)[
          missing_measure_values
        ],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    any(
      as.matrix(data[measure_columns]) < 0
    )
  ) {
    stop(
      basename(file),
      " contains negative counts.",
      call. = FALSE
    )
  }
  
  incomplete_units <- data |>
    distinct(unit, Treatment) |>
    count(
      unit,
      name = "n_treatments"
    ) |>
    filter(n_treatments != 2L)
  
  if (nrow(incomplete_units)) {
    stop(
      basename(file),
      " contains Animal x Hemisphere units without both BMI and SAL.",
      call. = FALSE
    )
  }
  
  data
}


# ---- Exact assignment mechanism -------------------------------------------

make_local_assignments <- function(
    n_sessions,
    n_comparison
) {
  combinations <- combn(
    n_sessions,
    n_comparison
  )
  
  assignments <- matrix(
    FALSE,
    nrow = ncol(combinations),
    ncol = n_sessions
  )
  
  for (i in seq_len(ncol(combinations))) {
    assignments[
      i,
      combinations[, i]
    ] <- TRUE
  }
  
  assignments
}

enumerate_unit_assignments <- function(data) {
  unit_rows <- split(
    seq_len(nrow(data)),
    data$unit
  )
  
  local_assignments <- map(
    unit_rows,
    function(rows) {
      make_local_assignments(
        n_sessions = length(rows),
        n_comparison = sum(
          data$Treatment[rows] ==
            comparison_level
        )
      )
    }
  )
  
  local_counts <- map_int(
    local_assignments,
    nrow
  )
  
  assignment_index_grid <- expand.grid(
    lapply(
      local_counts,
      seq_len
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(assignment_index_grid) <-
    names(unit_rows)
  
  list(
    unit_rows = unit_rows,
    local_assignments =
      local_assignments,
    assignment_index_grid =
      assignment_index_grid,
    n_exact =
      nrow(assignment_index_grid)
  )
}

declare_block_design <- function(data) {
  blocks <- as.character(data$unit)
  block_levels <- sort(unique(blocks))
  
  block_m <- vapply(
    block_levels,
    function(block_name) {
      sum(
        data$Treatment[
          blocks == block_name
        ] == comparison_level
      )
    },
    integer(1)
  )
  
  declaration <- declare_ra(
    blocks = blocks,
    block_m = block_m
  )
  
  list(
    declaration = declaration,
    n_exact =
      obtain_num_permutations(
        declaration
      )
  )
}

build_assignment_chunk <- function(
    assignment_index_grid,
    row_indices,
    unit_rows,
    local_assignments,
    n_sessions
) {
  assignment_matrix <- matrix(
    FALSE,
    nrow = length(row_indices),
    ncol = n_sessions
  )
  
  for (unit_name in names(unit_rows)) {
    assignment_matrix[
      ,
      unit_rows[[unit_name]]
    ] <- local_assignments[[unit_name]][
      assignment_index_grid[
        row_indices,
        unit_name
      ],
      ,
      drop = FALSE
    ]
  }
  
  assignment_matrix
}


# ---- Fixed estimable units -------------------------------------------------

find_stable_units <- function(
    data,
    unit_rows,
    local_assignments,
    denominator_name
) {
  denominator <- data[[denominator_name]]
  
  stable <- vapply(
    names(unit_rows),
    function(unit_name) {
      rows <- unit_rows[[unit_name]]
      local <- local_assignments[[unit_name]]
      values <- denominator[rows]
      
      comparison_denominator <-
        as.vector(local %*% values)
      
      reference_denominator <-
        sum(values) -
        comparison_denominator
      
      all(
        comparison_denominator > 0 &
          reference_denominator > 0
      )
    },
    logical(1)
  )
  
  names(stable)[stable]
}


# ---- Statistic -------------------------------------------------------------

row_mean_complete <- function(x) {
  if (!ncol(x)) {
    return(
      rep(
        NA_real_,
        nrow(x)
      )
    )
  }
  
  rowMeans(x)
}

compute_outcome_statistics <- function(
    comparison_assignment,
    data,
    unit_rows,
    animal_units,
    numerator_name,
    denominator_name
) {
  numerator <- data[[numerator_name]]
  denominator <- data[[denominator_name]]
  
  n_assignments <-
    nrow(comparison_assignment)
  
  unit_effects <- matrix(
    NA_real_,
    nrow = n_assignments,
    ncol = length(unit_rows)
  )
  
  colnames(unit_effects) <-
    names(unit_rows)
  
  for (unit_index in seq_along(unit_rows)) {
    rows <- unit_rows[[unit_index]]
    
    unit_assignment <-
      comparison_assignment[
        ,
        rows,
        drop = FALSE
      ]
    
    unit_numerator <- numerator[rows]
    unit_denominator <- denominator[rows]
    
    comparison_numerator <-
      as.vector(
        unit_assignment %*%
          unit_numerator
      )
    
    comparison_denominator <-
      as.vector(
        unit_assignment %*%
          unit_denominator
      )
    
    reference_numerator <-
      sum(unit_numerator) -
      comparison_numerator
    
    reference_denominator <-
      sum(unit_denominator) -
      comparison_denominator
    
    if (
      any(comparison_denominator <= 0) ||
      any(reference_denominator <= 0)
    ) {
      stop(
        "Internal error: an unstable unit reached statistic computation.",
        call. = FALSE
      )
    }
    
    unit_effects[, unit_index] <-
      comparison_numerator /
      comparison_denominator -
      reference_numerator /
      reference_denominator
  }
  
  animal_effects <- matrix(
    NA_real_,
    nrow = n_assignments,
    ncol = length(animal_units)
  )
  
  for (
    animal_index in
    seq_along(animal_units)
  ) {
    animal_effects[
      ,
      animal_index
    ] <- row_mean_complete(
      unit_effects[
        ,
        animal_units[[animal_index]],
        drop = FALSE
      ]
    )
  }
  
  row_mean_complete(animal_effects)
}

observed_assignment_vector <- function(data) {
  matrix(
    data$Treatment ==
      comparison_level,
    nrow = 1L
  )
}


# ---- Region analysis -------------------------------------------------------

run_exact_region_analysis <- function(file) {
  data <- read_region_data(file)
  region <- unique(data$Region)
  
  design <- declare_block_design(data)
  exact_design <-
    enumerate_unit_assignments(data)
  
  if (
    exact_design$n_exact !=
    design$n_exact
  ) {
    stop(
      "Exact assignment count does not match randomizr for ",
      region,
      ".",
      call. = FALSE
    )
  }
  
  n_exact <- exact_design$n_exact
  
  message(
    "Exact assignments for ",
    region,
    ": ",
    format(
      n_exact,
      big.mark = ","
    )
  )
  
  outcome_results <- vector(
    "list",
    nrow(outcome_definitions)
  )
  
  observed_unit_details <- vector(
    "list",
    nrow(outcome_definitions)
  )
  
  for (
    outcome_index in
    seq_len(nrow(outcome_definitions))
  ) {
    outcome_name <-
      outcome_definitions$outcome[
        outcome_index
      ]
    
    numerator_name <-
      outcome_definitions$numerator[
        outcome_index
      ]
    
    denominator_name <-
      outcome_definitions$denominator[
        outcome_index
      ]
    
    stable_unit_names <- find_stable_units(
      data = data,
      unit_rows =
        exact_design$unit_rows,
      local_assignments =
        exact_design$local_assignments,
      denominator_name =
        denominator_name
    )
    
    if (!length(stable_unit_names)) {
      stop(
        region,
        ": no permutation-stable units remain for ",
        outcome_name,
        ".",
        call. = FALSE
      )
    }
    
    stable_unit_rows <-
      exact_design$unit_rows[
        stable_unit_names
      ]
    
    unit_to_animal <- data |>
      filter(
        as.character(unit) %in%
          stable_unit_names
      ) |>
      distinct(
        unit,
        Animal
      )
    
    animal_units <- split(
      match(
        as.character(
          unit_to_animal$unit
        ),
        names(stable_unit_rows)
      ),
      unit_to_animal$Animal
    )
    
    observed_statistic <-
      compute_outcome_statistics(
        comparison_assignment =
          observed_assignment_vector(data),
        data = data,
        unit_rows =
          stable_unit_rows,
        animal_units =
          animal_units,
        numerator_name =
          numerator_name,
        denominator_name =
          denominator_name
      )[[1]]
    
    exact_statistics <- numeric(n_exact)
    
    chunk_starts <- seq.int(
      1L,
      n_exact,
      by = assignment_chunk_size
    )
    
    for (
      chunk_number in
      seq_along(chunk_starts)
    ) {
      chunk_start <-
        chunk_starts[[chunk_number]]
      
      chunk_end <- min(
        chunk_start +
          assignment_chunk_size -
          1L,
        n_exact
      )
      
      chunk_rows <-
        chunk_start:chunk_end
      
      assignment_chunk <-
        build_assignment_chunk(
          assignment_index_grid =
            exact_design$assignment_index_grid,
          row_indices =
            chunk_rows,
          unit_rows =
            exact_design$unit_rows,
          local_assignments =
            exact_design$local_assignments,
          n_sessions = nrow(data)
        )
      
      exact_statistics[chunk_rows] <-
        compute_outcome_statistics(
          comparison_assignment =
            assignment_chunk,
          data = data,
          unit_rows =
            stable_unit_rows,
          animal_units =
            animal_units,
          numerator_name =
            numerator_name,
          denominator_name =
            denominator_name
        )
    }
    
    permutation_center <-
      mean(exact_statistics)
    
    centered_observed <-
      observed_statistic -
      permutation_center
    
    raw_exact_p <- mean(
      abs(
        exact_statistics -
          permutation_center
      ) >=
        abs(centered_observed) -
        comparison_tolerance
    )
    
    outcome_results[[outcome_index]] <-
      tibble(
        Region = region,
        outcome = outcome_name,
        numerator = numerator_name,
        denominator =
          denominator_name,
        observed_effect =
          observed_statistic,
        observed_effect_percentage_points =
          100 * observed_statistic,
        ## Mean of the exact null; the test compares |observed - center| to
        ## |null - center|.
        exact_permutation_center =
          permutation_center,
        raw_exact_p =
          raw_exact_p,
        n_exact_assignments =
          n_exact,
        n_randomization_units =
          length(stable_unit_names),
        n_animals =
          length(animal_units),
        direction = case_when(
          observed_statistic > 0 ~
            "BMI higher",
          observed_statistic < 0 ~
            "BMI lower",
          TRUE ~ "No difference"
        )
      )
    
    observed_unit_details[[outcome_index]] <-
      data |>
      filter(
        as.character(unit) %in%
          stable_unit_names
      ) |>
      group_by(
        Region,
        Animal,
        Hemisphere,
        Treatment
      ) |>
      summarise(
        outcome = outcome_name,
        n_sessions = n(),
        numerator_total = sum(
          .data[[numerator_name]]
        ),
        denominator_total = sum(
          .data[[denominator_name]]
        ),
        unit_rate =
          numerator_total /
          denominator_total,
        .groups = "drop"
      ) |>
      select(
        Region,
        outcome,
        Animal,
        Hemisphere,
        Treatment,
        n_sessions,
        numerator_total,
        denominator_total,
        unit_rate
      )
  }
  
  results <- bind_rows(
    outcome_results
  ) |>
    mutate(
      bh_fdr_p = p.adjust(
        raw_exact_p,
        method = "BH"
      ),
      .after = raw_exact_p
    )
  
  list(
    results = results,
    hemisphere_details =
      bind_rows(observed_unit_details)
  )
}


# ---- Run and export --------------------------------------------------------

all_region_results <- map(
  input_files,
  run_exact_region_analysis
)

exact_results <- map_dfr(
  all_region_results,
  "results"
)

hemisphere_details <- map_dfr(
  all_region_results,
  "hemisphere_details"
)

## ---- Values behind the transition diagram ------------------------------------
## Weighted as the test: ratio of sums within unit, then equal weight to each
## hemisphere and each animal.

figure_values <- hemisphere_details |>
  group_by(Region, outcome, Treatment, Animal) |>
  summarise(
    animal_rate = mean(unit_rate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(Region, outcome, Treatment) |>
  summarise(
    n_animals   = n(),
    probability = mean(animal_rate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    arrow = if_else(
      outcome == "Pro-social behavioral transition",
      "Non-social -> Social",
      "Approach -> Social"
    )
  )

## second arrow out of Approach is the complement of the first
figure_values <- bind_rows(
  figure_values,
  figure_values |>
    filter(arrow == "Approach -> Social") |>
    mutate(
      arrow       = "Approach -> Non-social",
      probability = 1 - probability
    )
) |>
  mutate(percent = 100 * probability) |>
  arrange(Region, Treatment, arrow) |>
  select(Region, Treatment, arrow, outcome, n_animals, probability, percent)


dir.create(out_dir, showWarnings = FALSE)

output_file <- file.path(out_dir, "results.xlsx")

write_xlsx(
  list(
    exact_results = exact_results,
    hemisphere_details =
      hemisphere_details,
    figure_values = figure_values
  ),
  output_file
)

cat(
  "\nExact pooled permutation analyzes complete\n",
  "Output: ",
  normalizePath(
    output_file,
    winslash = "/"
  ),
  "\n",
  sep = ""
)