# ============================================================
# ACC behavioral batch analysis  
# Input: Directory containing CSV files downloaded from behaviorcloud
# Purpose: Cleans the files to calculate event durations and transition counts
# and to output the transitions of interest for next level analysis
#
# Analysis 1: Pro-social transitions
#   FROM any Group 1 state, Isolation, or Proximity
#   TO   Social Contact or any Group 3 state
#
# Analysis 2: First meaningful outcome after Approach
#   SOCIAL:
#     Social Contact, Mounting, Solicit Grooming,
#     Receive Groom, Give Groom
#   NON-SOCIAL / FAILURE:
#     Passive, Locomotion, Manipulation, Self Directed,
#     Isolation, Proximity, Withdraw
#
# Yawn and procedural markers are ignored.
# Stop markers terminate states only; they are never outcomes.
#
# Output sheets:
#   key_outputs
#   behavioral_transition_counts
#   approach_outcomes
#   label_summary
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(writexl)
})

session_duration_seconds <- 30 * 60

# Session identity is taken from the deposit manifest.
data_dir      <- "behavior_sequences"
manifest_file <- "MANIFEST.csv"

# One workbook per site, named so the contents are obvious. Hendrix 3 is a
# bilateral saline session used in both datasets and is written to both.
output_files <- c(ACCs = "ACCs_transitions.xlsx",
                  ACCg = "ACCg_transitions.xlsx")
tol <- 1e-6


# ---- Input -----------------------------------------------------------------

input_directory <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)

manifest_path <- file.path(input_directory, manifest_file)

if (!file.exists(manifest_path)) {
  stop(
    "Manifest not found: ",
    manifest_path,
    call. = FALSE
  )
}

manifest <- read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

input_files <- file.path(
  input_directory,
  manifest$new_filename
)

missing_sessions <- input_files[!file.exists(input_files)]

if (length(missing_sessions)) {
  stop(
    "Listed in the manifest but not on disk: ",
    paste(basename(missing_sessions), collapse = ", "),
    call. = FALSE
  )
}

if (!length(input_files)) {
  stop(
    "Manifest lists no sessions: ",
    manifest_path,
    call. = FALSE
  )
}

output_directory <- "02_transitions"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)



# ---- Event normalization ---------------------------------------------------

normalize_label_key <- function(x) {
  x |>
    as.character() |>
    str_squish() |>
    str_to_lower() |>
    str_replace_all("[^[:alnum:]]+", " ") |>
    str_squish()
}

event_aliases <- c(
  "passive" = "Passive",
  "locomotion" = "Locomotion",
  "manipulation" = "Manipulation",
  "self directed" = "Self Directed",
  "stop group 1" = "Stop Group 1",
  
  "social contact" = "Social Contact",
  "isolation" = "Isolation",
  "proximity" = "Proximity",
  
  "mounting" = "Mounting",
  "solicit grooming" = "Solicit Grooming",
  "receive groom" = "Receive Groom",
  "give groom" = "Give Groom",
  "stop groom" = "Stop Groom",
  
  "approach" = "Approach",
  "withdraw" = "Withdraw",
  "yawn" = "Yawn",
  
  "no response" = "No Response",
  "response" = "Response",
  "response latency" = "Response Latency",
  "stop resp latency" = "Stop Response Latency",
  "end 15 minutes" = "End 15 Minutes"
)

canonicalize_event <- function(x) {
  unname(event_aliases[normalize_label_key(x)])
}


# ---- Behavior dictionary ---------------------------------------------------

behavior_dictionary <- tribble(
  ~event,                  ~event_group,     ~event_type,    ~state_group,
  
  "Passive",               "Group 1",        "State",        "Group 1",
  "Locomotion",            "Group 1",        "State",        "Group 1",
  "Manipulation",          "Group 1",        "State",        "Group 1",
  "Self Directed",         "Group 1",        "State",        "Group 1",
  "Stop Group 1",          "Group 1",        "Stop marker",  "Group 1",
  
  "Social Contact",        "Group 2",        "State",        "Group 2",
  "Isolation",             "Group 2",        "State",        "Group 2",
  "Proximity",             "Group 2",        "State",        "Group 2",
  
  "Mounting",              "Group 3",        "State",        "Group 3",
  "Solicit Grooming",      "Group 3",        "State",        "Group 3",
  "Receive Groom",         "Group 3",        "State",        "Group 3",
  "Give Groom",            "Group 3",        "State",        "Group 3",
  "Stop Groom",            "Group 3",        "Stop marker",  "Group 3",
  
  "Approach",              "Discrete event", "Event",        NA_character_,
  "Withdraw",              "Discrete event", "Event",        NA_character_,
  
  "Yawn",                  "Ignored",        "Ignored",      NA_character_,
  "No Response",           "Procedural",     "Ignored",      NA_character_,
  "Response",              "Procedural",     "Ignored",      NA_character_,
  "Response Latency",      "Procedural",     "Ignored",      NA_character_,
  "Stop Response Latency", "Procedural",     "Ignored",      NA_character_,
  "End 15 Minutes",        "Procedural",     "Ignored",      NA_character_
)

stopifnot(
  !anyDuplicated(behavior_dictionary$event),
  setequal(
    unique(unname(event_aliases)),
    behavior_dictionary$event
  )
)

group1_states <- c(
  "Passive",
  "Locomotion",
  "Manipulation",
  "Self Directed"
)

non_social_states <- c(
  group1_states,
  "Isolation",
  "Proximity"
)

social_states <- c(
  "Social Contact",
  "Mounting",
  "Solicit Grooming",
  "Receive Groom",
  "Give Groom"
)

transition_states <- c(
  non_social_states,
  social_states
)

approach_outcomes_allowed <- c(
  non_social_states,
  social_states,
  "Withdraw"
)



# ---- Input -----------------------------------------------------------------

read_events <- function(file) {
  raw <- read_csv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  required_columns <- c("Event", "Start Time", "Duration")
  
  missing_columns <- setdiff(
    required_columns,
    names(raw)
  )
  
  if (length(missing_columns)) {
    stop(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  label_rows <- raw |>
    transmute(
      raw_event = str_squish(as.character(Event)),
      label_key = normalize_label_key(raw_event),
      event = canonicalize_event(raw_event)
    ) |>
    filter(
      !is.na(raw_event),
      raw_event != ""
    ) |>
    count(
      raw_event,
      label_key,
      event,
      name = "n_rows"
    ) |>
    mutate(
      input_file = basename(file),
      recognized = !is.na(event),
      .before = 1
    )
  
  unknown_labels <- label_rows |>
    filter(!recognized)
  
  if (nrow(unknown_labels)) {
    stop(
      "Unrecognized event label(s): ",
      paste(
        sort(unique(unknown_labels$raw_event)),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  events <- raw |>
    transmute(
      source_row = row_number(),
      raw_event = str_squish(as.character(Event)),
      event = canonicalize_event(raw_event),
      event_time = ymd_hms(
        `Start Time`,
        quiet = TRUE,
        tz = "UTC"
      ),
      duration = suppressWarnings(as.numeric(Duration))
    ) |>
    filter(
      !is.na(raw_event),
      raw_event != "",
      !is.na(event_time)
    ) |>
    left_join(
      behavior_dictionary,
      by = "event"
    ) |>
    arrange(event_time, source_row)
  
  if (
    anyNA(events$event_group) ||
    anyNA(events$event_type)
  ) {
    stop(
      "Internal behavior-dictionary join failed.",
      call. = FALSE
    )
  }
  
  list(
    n_input_rows = nrow(raw),
    events = events,
    label_rows = label_rows
  )
}


# ---- State bouts -----------------------------------------------------------

construct_state_bouts <- function(
    events,
    session_start,
    session_end
) {
  # zero-duration rows carry no time and are not transitions
  state_controls <- events |>
    filter(
      event_time >= session_start,
      event_time < session_end,
      event_type %in% c("State", "Stop marker"),
      is.na(duration) | duration > 0
    ) |>
    arrange(
      state_group,
      event_time,
      source_row
    ) |>
    group_by(
      state_group,
      event_time
    ) |>
    # Last-coded event wins within a group at tied onset times.
    slice_tail(n = 1) |>
    ungroup() |>
    arrange(
      state_group,
      event_time,
      source_row
    ) |>
    group_by(state_group) |>
    mutate(
      end_time = lead(
        event_time,
        default = session_end
      )
    ) |>
    ungroup()
  
  state_bouts <- state_controls |>
    filter(event_type == "State") |>
    transmute(
      source_row,
      event,
      event_group,
      state_group,
      start_time = event_time,
      end_time
    ) |>
    filter(end_time > start_time + tol) |>
    arrange(
      state_group,
      start_time,
      source_row
    ) |>
    group_by(state_group) |>
    mutate(
      new_bout =
        row_number() == 1L |
        event != lag(event) |
        abs(
          as.numeric(
            difftime(
              start_time,
              lag(end_time),
              units = "secs"
            )
          )
        ) > tol,
      bout_id =
        cumsum(replace_na(new_bout, TRUE))
    ) |>
    group_by(
      state_group,
      bout_id
    ) |>
    summarise(
      source_row = first(source_row),
      event = first(event),
      event_group = first(event_group),
      start_time = first(start_time),
      end_time = last(end_time),
      .groups = "drop"
    ) |>
    arrange(start_time, source_row)
  
  if (
    any(
      state_bouts$event %in%
      c("Stop Group 1", "Stop Groom")
    )
  ) {
    stop(
      "Internal error: stop marker emitted as a state bout.",
      call. = FALSE
    )
  }
  
  state_bouts
}


# ---- Pro-social transitions ------------------------------------------------

summarize_prosocial_transitions <- function(state_bouts) {
  sequence <- state_bouts |>
    filter(event %in% transition_states) |>
    transmute(
      source_row,
      event,
      event_group,
      start_time
    ) |>
    arrange(start_time, source_row) |>
    # Merge only immediately repeated identical labels.
    mutate(
      starts_new_element =
        row_number() == 1L |
        event != lag(event),
      sequence_number =
        cumsum(
          replace_na(starts_new_element, TRUE)
        )
    ) |>
    group_by(sequence_number) |>
    slice(1) |>
    ungroup() |>
    arrange(start_time, source_row)
  
  transitions <- sequence |>
    transmute(
      from_sequence_number = sequence_number,
      to_sequence_number = lead(sequence_number),
      from_event = event,
      to_event = lead(event)
    ) |>
    filter(
      !is.na(from_event),
      !is.na(to_event)
    ) |>
    mutate(
      from_class = if_else(
        from_event %in% non_social_states,
        "Non-social",
        "Social"
      ),
      to_class = if_else(
        to_event %in% social_states,
        "Social",
        "Non-social"
      ),
      pro_social =
        from_class == "Non-social" &
        to_class == "Social"
    )
  
  transition_counts <- transitions |>
    count(
      from_event,
      to_event,
      from_class,
      to_class,
      pro_social,
      name = "n_transitions"
    ) |>
    group_by(from_event) |>
    mutate(
      n_transitions_from_event =
        sum(n_transitions),
      proportion_from_event =
        n_transitions /
        n_transitions_from_event
    ) |>
    ungroup()
  
  n_from_non_social <- sum(
    transitions$from_class == "Non-social"
  )
  
  n_pro_social <- sum(
    transitions$pro_social
  )
  
  list(
    session_summary = tibble(
      n_from_non_social_behavior =
        n_from_non_social,
      n_pro_social_behavior_transitions =
        n_pro_social
    ),
    transition_counts = transition_counts
  )
}


# ---- First outcome after Approach -----------------------------------------

summarize_approach_outcomes <- function(
    events,
    state_bouts,
    session_start,
    session_end
) {
  state_onsets <- state_bouts |>
    filter(event %in% approach_outcomes_allowed) |>
    transmute(
      source_row,
      event,
      event_time = start_time
    )
  
  discrete_events <- events |>
    filter(
      event_time >= session_start,
      event_time < session_end,
      event %in% c("Approach", "Withdraw")
    ) |>
    transmute(
      source_row,
      event,
      event_time
    )
  
  event_stream <- bind_rows(
    state_onsets,
    discrete_events
  ) |>
    arrange(event_time, source_row)
  
  approaches <- event_stream |>
    filter(event == "Approach") |>
    arrange(event_time, source_row) |>
    mutate(
      approach_number = row_number(),
      next_approach_time = lead(
        event_time,
        default = session_end
      )
    )
  
  classify_one <- function(
    approach_number,
    approach_time,
    approach_source_row,
    next_approach_time
  ) {
    candidates <- event_stream |>
      filter(
        event != "Approach",
        event %in% approach_outcomes_allowed,
        
        # Later than the current Approach; allow tied timestamps only
        # when the candidate was coded later in source-row order.
        (
          event_time > approach_time + tol |
            (
              abs(
                as.numeric(
                  difftime(
                    event_time,
                    approach_time,
                    units = "secs"
                  )
                )
              ) <= tol &
                source_row > approach_source_row
            )
        ),
        
        # Stop at the next Approach. Events exactly at the next
        # Approach timestamp are conservatively assigned to neither.
        event_time < next_approach_time - tol,
        
        event_time <= session_end + tol
      ) |>
      arrange(event_time, source_row)
    
    if (!nrow(candidates)) {
      return(
        tibble(
          approach_number = approach_number,
          approach_time = approach_time,
          approach_source_row =
            approach_source_row,
          next_outcome_event =
            NA_character_,
          next_outcome_time =
            as.POSIXct(NA, tz = "UTC"),
          next_outcome_class =
            "Unresolved",
          social_success = FALSE,
          classifiable = FALSE
        )
      )
    }
    
    next_event <- candidates |>
      slice(1)
    
    next_event_name <- next_event$event[[1]]
    
    next_class <- if (
      next_event_name %in% social_states
    ) {
      "Social"
    } else {
      "Non-social or Withdraw"
    }
    
    tibble(
      approach_number = approach_number,
      approach_time = approach_time,
      approach_source_row =
        approach_source_row,
      next_outcome_event =
        next_event_name,
      next_outcome_time =
        next_event$event_time[[1]],
      next_outcome_class =
        next_class,
      social_success =
        next_class == "Social",
      classifiable = TRUE
    )
  }
  
  approach_outcomes <- if (nrow(approaches)) {
    pmap_dfr(
      list(
        approaches$approach_number,
        approaches$event_time,
        approaches$source_row,
        approaches$next_approach_time
      ),
      classify_one
    )
  } else {
    tibble(
      approach_number = integer(),
      approach_time =
        as.POSIXct(character(), tz = "UTC"),
      approach_source_row = integer(),
      next_outcome_event = character(),
      next_outcome_time =
        as.POSIXct(character(), tz = "UTC"),
      next_outcome_class = character(),
      social_success = logical(),
      classifiable = logical()
    )
  }
  
  n_approaches <- nrow(approach_outcomes)
  
  n_classifiable <- sum(
    approach_outcomes$classifiable
  )
  
  n_social <- sum(
    approach_outcomes$social_success,
    na.rm = TRUE
  )
  
  n_non_social <- sum(
    approach_outcomes$next_outcome_class ==
      "Non-social or Withdraw",
    na.rm = TRUE
  )
  
  n_unresolved <- sum(
    approach_outcomes$next_outcome_class ==
      "Unresolved",
    na.rm = TRUE
  )
  
  list(
    session_summary = tibble(
      n_approaches_total = n_approaches,
      n_classifiable_approaches =
        n_classifiable,
      n_approach_social_success =
        n_social,
      n_approach_non_social_or_withdraw =
        n_non_social,
      n_approach_unresolved =
        n_unresolved
    ),
    approach_outcomes = approach_outcomes
  )
}


# ---- Analyze one file ------------------------------------------------------

analyze_file <- function(file) {
  imported <- read_events(file)
  events <- imported$events
  
  behavioral_events <- events |>
    filter(
      event_type %in%
        c("State", "Stop marker", "Event")
    )
  
  if (!nrow(behavioral_events)) {
    stop(
      "No recognized behavioral rows found.",
      call. = FALSE
    )
  }
  
  session_start <- min(
    behavioral_events$event_time
  )
  
  session_end <-
    session_start +
    seconds(session_duration_seconds)
  
  state_bouts <- construct_state_bouts(
    events,
    session_start,
    session_end
  )
  
  transition_results <-
    summarize_prosocial_transitions(
      state_bouts
    )
  
  approach_results <-
    summarize_approach_outcomes(
      events,
      state_bouts,
      session_start,
      session_end
    )
  
  session_summary <- bind_cols(
    tibble(
      input_file = basename(file),
      processing_status = "OK",
      processing_error = NA_character_,
      session_start = session_start,
      session_end = session_end,
      n_input_rows = imported$n_input_rows
    ),
    transition_results$session_summary,
    approach_results$session_summary
  )
  
  transition_counts <-
    transition_results$transition_counts |>
    mutate(
      input_file = basename(file),
      .before = 1
    )
  
  approach_outcomes <-
    approach_results$approach_outcomes |>
    mutate(
      input_file = basename(file),
      .before = 1
    )
  
  label_summary <- imported$label_rows |>
    left_join(
      behavior_dictionary,
      by = "event"
    ) |>
    mutate(
      analysis_role = case_when(
        event_type == "State" ~
          "Behavioral state",
        event_type == "Stop marker" ~
          "State terminator only",
        event_type == "Event" ~
          "Discrete behavior event",
        TRUE ~ "Ignored"
      )
    )
  
  list(
    session_summary = session_summary,
    transition_counts = transition_counts,
    approach_outcomes = approach_outcomes,
    label_summary = label_summary
  )
}


# ---- Batch run -------------------------------------------------------------

empty_transition_counts <- tibble(
  input_file = character(),
  from_event = character(),
  to_event = character(),
  from_class = character(),
  to_class = character(),
  pro_social = logical(),
  n_transitions = integer(),
  n_transitions_from_event = integer(),
  proportion_from_event = numeric()
)

empty_approach_outcomes <- tibble(
  input_file = character(),
  approach_number = integer(),
  approach_time =
    as.POSIXct(character(), tz = "UTC"),
  approach_source_row = integer(),
  next_outcome_event = character(),
  next_outcome_time =
    as.POSIXct(character(), tz = "UTC"),
  next_outcome_class = character(),
  social_success = logical(),
  classifiable = logical()
)

empty_label_summary <- tibble(
  input_file = character(),
  raw_event = character(),
  label_key = character(),
  event = character(),
  n_rows = integer(),
  recognized = logical(),
  event_group = character(),
  event_type = character(),
  state_group = character(),
  analysis_role = character()
)

results <- map(
  input_files,
  \(file) {
    message(
      "Processing: ",
      basename(file)
    )
    
    tryCatch(
      analyze_file(file),
      error = function(e) {
        list(
          session_summary = tibble(
            input_file = basename(file),
            processing_status = "ERROR",
            processing_error =
              conditionMessage(e),
            session_start =
              as.POSIXct(NA, tz = "UTC"),
            session_end =
              as.POSIXct(NA, tz = "UTC"),
            n_input_rows = NA_integer_
          ),
          transition_counts =
            empty_transition_counts,
          approach_outcomes =
            empty_approach_outcomes,
          label_summary =
            empty_label_summary
        )
      }
    )
  }
)

key_outputs <- map_dfr(
  results,
  "session_summary"
) |>
  arrange(input_file)

behavioral_transition_counts <- map_dfr(
  results,
  "transition_counts"
) |>
  arrange(
    input_file,
    from_event,
    to_event
  )

approach_outcomes <- map_dfr(
  results,
  "approach_outcomes"
) |>
  arrange(
    input_file,
    approach_number
  )

label_summary <- map_dfr(
  results,
  "label_summary"
) |>
  group_by(
    raw_event,
    event,
    event_group,
    event_type,
    state_group,
    analysis_role
  ) |>
  summarise(
    n_rows = sum(n_rows),
    n_files = n_distinct(input_file),
    .groups = "drop"
  ) |>
  arrange(event, raw_event)


# ---- Identity from the manifest -------------------------------------------
# Used In selects the workbook; a session listed in both is written to both.

manifest <- manifest |>
  transmute(
    input_file = new_filename,
    Animal     = animal,
    Session    = session,
    `Video Code` = video_code,
    Hemisphere = hemisphere,
    Treatment  = treatment,
    used_in    = used_in
  )

unlisted <- setdiff(key_outputs$input_file, manifest$input_file)

if (length(unlisted)) {
  stop(
    "Not in the manifest: ",
    paste(unlisted, collapse = ", "),
    call. = FALSE
  )
}

key_outputs <- manifest |>
  inner_join(key_outputs, by = "input_file") |>
  arrange(Animal, Session)


# ---- Export ---------------------------------------------------------------

for (site in names(output_files)) {
  site_files <- key_outputs$input_file[grepl(site, key_outputs$used_in, fixed = TRUE)]

  write_xlsx(
    list(
      key_outputs = key_outputs |>
        filter(input_file %in% site_files) |>
        select(-used_in),
      behavioral_transition_counts =
        behavioral_transition_counts |>
        filter(input_file %in% site_files),
      approach_outcomes = approach_outcomes |>
        filter(input_file %in% site_files),
      label_summary = label_summary
    ),
    file.path(output_directory, output_files[[site]])
  )

  cat(
    sprintf(
      "%-4s %2d sessions -> %s\n",
      site,
      length(site_files),
      output_files[[site]]
    )
  )
}

cat(
  "\nBatch analysis complete\n",
  "Files: ",
  length(input_files),
  " | Successful: ",
  sum(key_outputs$processing_status == "OK"),
  " | Errors: ",
  sum(key_outputs$processing_status == "ERROR"),
  "\nOutput directory: ",
  normalizePath(
    output_directory,
    winslash = "/"
  ),
  "\n",
  sep = ""
)