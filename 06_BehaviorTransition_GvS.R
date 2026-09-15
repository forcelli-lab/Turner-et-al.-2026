# ============================================================
# Exact permutation test: ACCs versus ACCg, BMI sessions only
#
# Input: 
#  2 files (ACCs_transitions.xlsx, ACCg_transitions.xlsx) which are output by the behavior transition
#  processing script. Animal, hemisphere and treatment are joined from
#  MANIFEST.csv by that script; no manual annotation step.
#
# Outcomes:
#   1. Pro-social transition probability
#   2. Social outcome following Approach
#
# Statistic:
#   - Retain BMI sessions only.
#   - Ignore/collapse hemispher out of necessity
#   - Pool counts across BMI sessions within Animal x Site.
#   - Compute ACCs - ACCg within each animal.
#   - Average across animals
#
#   Note: every animal must have positive pooled denominators
#   at BOTH sites under EVERY within-animal assignment for
#   both outcomes. The script stops if this condition fails 
#   (prevents weird divide by zeros)
#
# Multiplicity:
#   Raw two-sided exact p-values and BH FDR across the two outcomes.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(writexl)
  library(randomizr)
})

comparison_level <- "ACCs"
bmi_level <- "BMI"

assignment_chunk_size <- 50000L
comparison_tolerance <- 1e-12

# written by 02_ACC_behavioral_batch_analysis.R
data_dir <- "02_transitions"
out_dir  <- "06_transition_GvS"

input_files <- file.path(
  data_dir,
  c(ACCs = "ACCs_transitions.xlsx",
    ACCg = "ACCg_transitions.xlsx")
)
names(input_files) <- c("ACCs", "ACCg")

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files)) {
  stop(
    "Missing required workbook(s): ",
    paste(basename(missing_files), collapse = ", "),
    ". Run 02_ACC_behavioral_batch_analysis.R first.",
    call. = FALSE
  )
}

outcome_definitions <- tribble(
  ~outcome, ~numerator, ~denominator,
  
  "Pro-social behavioral transition",
  "n_pro_social_behavior_transitions",
  "n_from_non_social_behavior",
  
  "Social outcome following Approach",
  "n_approach_social_success",
  "n_classifiable_approaches"
)


# ---- Input files -----------------------------------------------------------

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files)) {
  stop(
    "Missing required workbook(s) in the working directory: ",
    paste(unname(missing_files), collapse = ", "),
    call. = FALSE
  )
}


# ---- Input helpers ---------------------------------------------------------

normalize_treatment <- function(x) {
  value <- toupper(
    trimws(as.character(x))
  )
  
  case_when(
    value %in% c("BMI", "DRUG") ~
      bmi_level,
    value %in%
      c("SAL", "SALINE", "VEHICLE", "CONTROL") ~
      "SAL",
    TRUE ~ value
  )
}

read_site_data <- function(file, site) {
  message("Reading: ", file, " as ", site)
  
  
  raw <- read_excel(
    file,
    sheet = "key_outputs"
  )
  
  required_columns <- unique(c(
    "input_file",
    "Animal",
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
      Site = site,
      input_file =
        as.character(input_file),
      Animal =
        toupper(trimws(as.character(Animal))),
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
    filter(Treatment == bmi_level)
  
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

enumerate_animal_assignments <- function(data) {
  animal_rows <- split(
    seq_len(nrow(data)),
    data$Animal
  )
  
  local_assignments <- map(
    animal_rows,
    function(rows) {
      make_local_assignments(
        n_sessions = length(rows),
        n_comparison = sum(
          data$Site[rows] ==
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
    names(animal_rows)
  
  list(
    animal_rows = animal_rows,
    local_assignments =
      local_assignments,
    assignment_index_grid =
      assignment_index_grid,
    n_exact =
      nrow(assignment_index_grid)
  )
}

declare_block_design <- function(data) {
  blocks <- as.character(data$Animal)
  block_levels <- sort(unique(blocks))
  
  block_m <- vapply(
    block_levels,
    function(animal_name) {
      sum(
        data$Site[
          blocks == animal_name
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
    animal_rows,
    local_assignments,
    n_sessions
) {
  assignment_matrix <- matrix(
    FALSE,
    nrow = length(row_indices),
    ncol = n_sessions
  )
  
  for (animal_name in names(animal_rows)) {
    assignment_matrix[
      ,
      animal_rows[[animal_name]]
    ] <- local_assignments[[animal_name]][
      assignment_index_grid[
        row_indices,
        animal_name
      ],
      ,
      drop = FALSE
    ]
  }
  
  assignment_matrix
}


# ---- Fixed estimable animals ----------------------------------------------

find_stable_animals <- function(
    data,
    animal_rows,
    local_assignments,
    denominator_name
) {
  denominator <- data[[denominator_name]]
  
  stable <- vapply(
    names(animal_rows),
    function(animal_name) {
      rows <- animal_rows[[animal_name]]
      local <- local_assignments[[animal_name]]
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

compute_outcome_statistics <- function(
    comparison_assignment,
    data,
    animal_rows,
    numerator_name,
    denominator_name
) {
  numerator <- data[[numerator_name]]
  denominator <- data[[denominator_name]]
  
  n_assignments <-
    nrow(comparison_assignment)
  
  animal_effects <- matrix(
    NA_real_,
    nrow = n_assignments,
    ncol = length(animal_rows)
  )
  
  colnames(animal_effects) <-
    names(animal_rows)
  
  for (
    animal_index in
    seq_along(animal_rows)
  ) {
    rows <- animal_rows[[animal_index]]
    
    animal_assignment <-
      comparison_assignment[
        ,
        rows,
        drop = FALSE
      ]
    
    animal_numerator <- numerator[rows]
    animal_denominator <- denominator[rows]
    
    comparison_numerator <-
      as.vector(
        animal_assignment %*%
          animal_numerator
      )
    
    comparison_denominator <-
      as.vector(
        animal_assignment %*%
          animal_denominator
      )
    
    reference_numerator <-
      sum(animal_numerator) -
      comparison_numerator
    
    reference_denominator <-
      sum(animal_denominator) -
      comparison_denominator
    
    if (
      any(comparison_denominator <= 0) ||
      any(reference_denominator <= 0)
    ) {
      stop(
        "Internal error: an unstable animal reached statistic computation.",
        call. = FALSE
      )
    }
    
    animal_effects[, animal_index] <-
      comparison_numerator /
      comparison_denominator -
      reference_numerator /
      reference_denominator
  }
  
  rowMeans(animal_effects)
}

observed_assignment_vector <- function(data) {
  matrix(
    data$Site ==
      comparison_level,
    nrow = 1L
  )
}


# ---- Analysis --------------------------------------------------------------

run_exact_site_analysis <- function(files) {
  data <- map2_dfr(
    .x = unname(files),
    .y = names(files),
    .f = read_site_data
  ) |>
    arrange(
      Animal,
      Site,
      input_file
    )
  
  site_counts <- data |>
    distinct(
      Animal,
      Site
    ) |>
    count(
      Animal,
      name = "n_sites"
    )
  
  animals_with_both_sites <-
    site_counts |>
    filter(n_sites == 2L) |>
    pull(Animal)
  
  excluded_animals <- setdiff(
    unique(data$Animal),
    animals_with_both_sites
  )
  
  if (length(excluded_animals)) {
    message(
      "Excluding animals without BMI sessions at both sites: ",
      paste(
        excluded_animals,
        collapse = ", "
      )
    )
  }
  
  data <- data |>
    filter(
      Animal %in%
        animals_with_both_sites
    )
  
  if (!nrow(data)) {
    stop(
      "No animals have BMI sessions at both ACCs and ACCg.",
      call. = FALSE
    )
  }
  
  design <- declare_block_design(data)
  
  exact_design <-
    enumerate_animal_assignments(data)
  
  if (
    exact_design$n_exact !=
    design$n_exact
  ) {
    stop(
      "Exact assignment count does not match randomizr.",
      call. = FALSE
    )
  }
  
  n_exact <- exact_design$n_exact
  
  message(
    "Exact assignments: ",
    format(
      n_exact,
      big.mark = ","
    )
  )
  
  outcome_results <- vector(
    "list",
    nrow(outcome_definitions)
  )
  
  observed_animal_details <- vector(
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
    
    eligible_animals <- names(
      exact_design$animal_rows
    )
    
    stable_animals <- find_stable_animals(
      data = data,
      animal_rows =
        exact_design$animal_rows,
      local_assignments =
        exact_design$local_assignments,
      denominator_name =
        denominator_name
    )
    
    unstable_animals <- setdiff(
      eligible_animals,
      stable_animals
    )
    
    if (length(unstable_animals)) {
      stop(
        "The fixed-animal denominator requirement failed for ",
        outcome_name,
        ". Animal(s) with a zero pooled site denominator under ",
        "at least one legal assignment: ",
        paste(unstable_animals, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    stable_animal_rows <-
      exact_design$animal_rows
    
    observed_statistic <-
      compute_outcome_statistics(
        comparison_assignment =
          observed_assignment_vector(data),
        data = data,
        animal_rows =
          stable_animal_rows,
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
          animal_rows =
            exact_design$animal_rows,
          local_assignments =
            exact_design$local_assignments,
          n_sessions = nrow(data)
        )
      
      exact_statistics[chunk_rows] <-
        compute_outcome_statistics(
          comparison_assignment =
            assignment_chunk,
          data = data,
          animal_rows =
            stable_animal_rows,
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
        comparison = "ACCs minus ACCg",
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
        n_animals =
          length(eligible_animals),
        direction = case_when(
          observed_statistic > 0 ~
            "ACCs higher",
          observed_statistic < 0 ~
            "ACCs lower",
          TRUE ~ "No difference"
        )
      )
    
    observed_animal_details[[outcome_index]] <-
      data |>
      filter(
        Animal %in%
          stable_animals
      ) |>
      group_by(
        Animal,
        Site
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
        outcome,
        Animal,
        Site,
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
    animal_details =
      bind_rows(observed_animal_details)
  )
}


# ---- Run and export --------------------------------------------------------

analysis <- run_exact_site_analysis(
  input_files
)

## ---- Values behind the transition diagram ------------------------------------
## Weighted as the test: ratio of sums within animal, then equal weight to
## each animal.

figure_values <- analysis$animal_details |>
  group_by(outcome, Site) |>
  summarise(
    n_animals   = n(),
    probability = mean(unit_rate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    arrow = if_else(
      outcome == "Pro-social behavioral transition",
      "Non-social -> Social",
      "Approach -> Social"
    )
  )

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
  arrange(Site, arrow) |>
  select(Site, arrow, outcome, n_animals, probability, percent)


dir.create(out_dir, showWarnings = FALSE)

output_file <- file.path(out_dir, "results.xlsx")

write_xlsx(
  list(
    exact_results =
      analysis$results,
    animal_details =
      analysis$animal_details,
    figure_values = figure_values
  ),
  output_file
)

cat(
  "\nExact pooled ACCs-versus-ACCg BMI permutation analyzes complete\n",
  "Output: ",
  normalizePath(
    output_file,
    winslash = "/"
  ),
  "\n",
  sep = ""
)