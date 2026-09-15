## Behavioral measures for Turner et al., computed from the individual
## behavior-sequence exports in behavior_sequences/.
##
## Parameters set below:
##   - first 1800 s of each session
##   - unscored lead-in (usually only a few seconds) assigned to the first scored behavior
##   - physical grooming scored outside Social Contact is reassigned to it
##   - Social Contact bouts: any scored interruption starts a new bout
##   - zero-duration events contribute neither time nor a bout
##   - Grooming = Give Groom + Receive Groom only

## ---------------- CONFIG ----------------
data_dir       <- "behavior_sequences"
manifest_file  <- "MANIFEST.csv"          # supplies animal / session / video / hemisphere / treatment
out_dir        <- "01_measures"
out_measures   <- "measures.csv"
out_qc         <- "qc.csv"

# Per-site files used for estimation-plot, permutation and MDE scripts.
out_by_site    <- c(ACCs = "ACCs_BMI.csv",
                    ACCg = "ACCg_BMI.csv")

# Region comparison (supplement): BMI sessions only, restricted to animals that
# received BMI at both sites, with a Site column in place of Treatment.
out_region     <- "ACCgvss_BMI.csv"

dir.create(out_dir, showWarnings = FALSE)

## Column names in the source files ------------------------------------------
col_event      <- "Event"
col_group      <- "Group"
col_start      <- "Latency"               # seconds from the file's time origin
col_duration   <- "Duration"              # seconds

## Group labels --------------------------------------------------------------
grp1           <- "Behavior Group 1"
grp2           <- "Behavior Group 2"
grp3           <- "Behavior Group 3"

## Event vocabulary ----------------------------------------------------------
# Coders varied spacing, hyphenation and case; whitespace and hyphens are
# collapsed and case ignored, so only genuine wording differences need an entry.
event_aliases  <- c(
  "self directed"      = "Self Directed",
  "receive groom"      = "Receive Groom",
  "receive grooming"   = "Receive Groom",
  "give groom"         = "Give Groom",
  "give grooming"      = "Give Groom",
  "solicit grooming"   = "Solicit Grooming",
  "stop group 1"       = "Stop Group 1",
  "stop groom"         = "Stop Groom",
  "end 15 minutes"     = "END 15 Minutes"
)

# Non-scored codes: present in the files. Anything
# outside this list and the measures below is flagged in the QC table.
ev_ignore      <- c("Approach", "Withdraw", "Yawn", "END 15 Minutes",
                    "No Response", "Response", "Response Latency",
                    "STOP RESP. LATENCY")

## Measures ------------------------------------------------------------------
# Duration measures: output column -> (group, events summed)
dur_measures <- list(
  locomotion       = list(grp1, "Locomotion"),
  manipulation     = list(grp1, "Manipulation"),
  self_directed    = list(grp1, "Self Directed"),
  passive          = list(grp1, "Passive"),
  stop_group1      = list(grp1, "Stop Group 1"),
  isolation        = list(grp2, "Isolation"),
  proximity        = list(grp2, "Proximity"),
  social_contact   = list(grp2, "Social Contact"),
  give_groom       = list(grp3, "Give Groom"),
  receive_groom    = list(grp3, "Receive Groom"),
  mounting         = list(grp3, "Mounting"),
  solicit_grooming = list(grp3, "Solicit Grooming")
)

# Point-event counts: output column -> event name
cnt_measures <- list(approach_n = "Approach", withdraw_n = "Withdraw")

# Derived quantities.
ev_passive     <- "Passive"
ev_alone       <- c("Isolation", "Proximity")          # Group 2 states = not in contact
ev_contact     <- "Social Contact"
ev_groom_excl  <- c("Mounting", "Solicit Grooming",    # Group 3 states disqualifying "alone"
                    "Receive Groom", "Give Groom")
ev_groom_sum   <- c("Give Groom", "Receive Groom")     # what the analysis sheets call "Grooming"

## Contact implied by physical grooming ---------------------------------------
# Give Groom, Receive Groom and Mounting can't occur without physical contact.
# Where one is scored while Group 2 says Isolation or Proximity the streams
# contradict each other; that time is reassigned to Social Contact and taken
# off the other state, so the session still sums to 1800s. Solicit Grooming
# is excluded -- soliciting from proximity is not a contradiction.
ev_contact_implied <- c("Give Groom", "Receive Groom", "Mounting")

## Lead-in handling ----------------------------------------------------------
# Some files begin scoring after the latency clock starts. The first state in
# these streams is extended back to time 0. Group 3 is excluded: its late start
# means "no grooming yet".
lead_in_groups <- c(grp1, grp2)

## Analysis window -----------------------------------------------------------
window_start   <- 0
window_end     <- 1800                    # seconds; NA = last event in file

## Segment handling ----------------------------------------------------------
# Group streams are scored independently, so transitions are offset by up to a few
# tenths of a second. A nonzero threshold breaks the exact decomposition
# passive_alone + passive_contact = passive, so it stays 0.
min_fragment   <- 0

## Bout definition -----------------------------------------------------------
# Two contact episodes count as one bout only if the break between them is no
# longer than this. 0 = any scored interruption starts a new bout. Affects
# SocBoutCount and the bout durations only -- never a duration total.
bout_gap_tol   <- 0

## Reporting -----------------------------------------------------------------
digits         <- 2

## Output column names --------------------------------------------------------
# Internal names are used throughout the code; these are the headers written to
# the CSV.
out_names <- c(
  animal           = "Animal",
  session          = "Session",
  video_code       = "Video Code",
  hemisphere       = "Hemisphere",
  treatment        = "Treatment",
  locomotion       = "Locomotion",
  manipulation     = "Object Manipulation",
  passive_alone    = "Passive Alone",
  self_directed    = "SelfDirected",
  social_contact   = "SocialContact",
  soc_bout_n       = "SocBoutCount",
  soc_bout_med     = "SocBoutDur",
  grooming         = "Grooming",
  passive_contact  = "Passive in Cont.",
  proximity        = "Proximity",
  isolation        = "Alone",
  passive          = "Passive",
  stop_group1      = "Stop Group 1",
  give_groom       = "Give Groom",
  receive_groom    = "Receive Groom",
  mounting         = "Mounting",
  solicit_grooming = "Solicit Grooming",
  alone_total      = "Not In Contact",
  approach_n       = "Approach",
  withdraw_n       = "Withdraw",
  window_end       = "Window (s)",
  used_in          = "Used In",
  file             = "Source File"
)

# ============================================================================
# INTERVAL HELPERS
# ============================================================================
# Intervals are two-column matrices of [start, end), sorted by start.

iv_new <- function(start, end) {
  m <- cbind(start = start, end = end)
  m[order(m[, 1]), , drop = FALSE]
}

iv_merge <- function(a, tol = 0) {
  if (nrow(a) == 0) return(a)
  a <- a[order(a[, 1]), , drop = FALSE]
  keep <- a[1, , drop = FALSE]
  for (i in seq_len(nrow(a))[-1]) {
    last <- nrow(keep)
    if (a[i, 1] <= keep[last, 2] + tol) {
      keep[last, 2] <- max(keep[last, 2], a[i, 2])
    } else {
      keep <- rbind(keep, a[i, , drop = FALSE])
    }
  }
  keep
}

# Intervals shorter than this are treated as empty. Stream boundaries are
# scored to 0.1 s but stored as start + duration, so an intersection can come
# out a fraction of a nanosecond long purely from floating-point error.
iv_eps <- 1e-6

iv_intersect <- function(a, b) {
  if (nrow(a) == 0 || nrow(b) == 0) return(iv_new(numeric(0), numeric(0)))
  out <- vector("list", nrow(a))
  for (i in seq_len(nrow(a))) {
    lo <- pmax(a[i, 1], b[, 1])
    hi <- pmin(a[i, 2], b[, 2])
    ok <- hi > lo + iv_eps
    out[[i]] <- cbind(start = lo[ok], end = hi[ok])
  }
  iv_new(unlist(lapply(out, function(x) x[, 1])),
         unlist(lapply(out, function(x) x[, 2])))
}

iv_complement <- function(a, lo, hi) {
  a <- iv_merge(a)
  a <- a[a[, 2] > lo & a[, 1] < hi, , drop = FALSE]
  m <- matrix(c(lo, as.vector(t(a)), hi), ncol = 2, byrow = TRUE)
  m[m[, 2] > m[, 1] + iv_eps, , drop = FALSE]
}

iv_clip  <- function(a, lo, hi) iv_intersect(a, iv_new(lo, hi))
iv_total <- function(a) sum(a[, 2] - a[, 1])

# ============================================================================
# LOADING
# ============================================================================

canonical <- function(x) {
  key <- tolower(gsub("[[:space:]-]+", " ", trimws(x)))
  hit <- event_aliases[key]
  ifelse(is.na(hit), trimws(gsub("[[:space:]]+", " ", x)), hit)
}

read_session <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(
    event = canonical(d[[col_event]]),
    group = trimws(ifelse(is.na(d[[col_group]]), "", d[[col_group]])),
    start = as.numeric(d[[col_start]]),
    end   = as.numeric(d[[col_start]]) + as.numeric(d[[col_duration]]),
    stringsAsFactors = FALSE
  )
}

pick <- function(d, group, events) {
  s <- d[d$group == group & d$event %in% events, ]
  iv_new(s$start, s$end)
}

# ============================================================================
# IMPLIED CONTACT
# ============================================================================
# Split any Isolation/Proximity interval overlapping physical grooming and
# relabel the overlapping piece as Social Contact.

apply_implied_contact <- function(d) {
  gi <- iv_merge(pick(d, grp3, ev_contact_implied))
  if (nrow(gi) == 0) return(d)
  rows <- which(d$group == grp2 & d$event %in% ev_alone)
  add <- NULL; drop <- integer(0)
  for (i in rows) {
    span <- iv_new(d$start[i], d$end[i])
    ov <- iv_intersect(span, gi)
    if (nrow(ov) == 0) next
    drop <- c(drop, i)
    keep <- iv_intersect(span, iv_complement(gi, d$start[i], d$end[i]))
    if (nrow(keep) > 0)
      add <- rbind(add, data.frame(event = d$event[i], group = grp2,
                                   start = keep[, 1], end = keep[, 2],
                                   stringsAsFactors = FALSE))
    add <- rbind(add, data.frame(event = ev_contact, group = grp2,
                                 start = ov[, 1], end = ov[, 2],
                                 stringsAsFactors = FALSE))
  }
  if (length(drop) > 0) d <- d[-drop, ]
  if (!is.null(add)) d <- rbind(d, add)
  d[order(d$start), ]
}

# ============================================================================
# QC
# ============================================================================
# Each stream should tile its own extent with no gaps and no overlaps.
# Violations are reported in the qc file

# QC covers the analysis window only.
stream_qc <- function(d, group, file, t0, t1) {
  s <- d[d$group == group & !(d$event %in% ev_ignore), ]
  a <- if (nrow(s)) iv_clip(iv_new(s$start, s$end), t0, t1) else iv_new(numeric(0), numeric(0))
  if (nrow(a) == 0) {
    return(data.frame(file = file, group = group, n = 0, covered = 0,
                      first = NA, last = NA, overlap = 0, gap = 0))
  }
  reach <- cummax(a[, 2])                 # running max, so nesting is not read as a gap
  gaps <- a[-1, 1] - reach[-nrow(a)]
  data.frame(file = file, group = group, n = nrow(a), covered = iv_total(a),
             first = a[1, 1], last = max(a[, 2]),
             overlap = -sum(pmin(gaps, 0)), gap = sum(pmax(gaps, 0)))
}

# ============================================================================
# PER-SESSION EXTRACTION
# ============================================================================

analyze <- function(path) {
  d <- read_session(path)

  for (g in lead_in_groups) {
    i <- which(d$group == g & !(d$event %in% ev_ignore))
    j <- i[which.min(d$start[i])]
    d$start[j] <- min(0, d$start[j])
  }

  d <- apply_implied_contact(d)

  t0 <- if (is.na(window_start)) min(d$start) else window_start
  t1 <- if (is.na(window_end))   max(d$end)   else window_end

  vals <- lapply(dur_measures, function(m)
    iv_total(iv_clip(pick(d, m[[1]], m[[2]]), t0, t1)))

  cnts <- lapply(cnt_measures, function(e)
    sum(d$event == e & d$start >= t0 & d$start < t1))

  passive <- iv_clip(pick(d, grp1, ev_passive), t0, t1)
  alone   <- iv_merge(iv_clip(pick(d, grp2, ev_alone), t0, t1))
  contact <- iv_merge(iv_clip(pick(d, grp2, ev_contact), t0, t1))
  groom   <- iv_merge(iv_clip(pick(d, grp3, ev_groom_excl), t0, t1))

  seg <- iv_intersect(iv_intersect(passive, alone), iv_complement(groom, t0, t1))
  seg <- seg[(seg[, 2] - seg[, 1]) >= min_fragment, , drop = FALSE]
  cseg <- iv_intersect(passive, contact)
  cseg <- cseg[(cseg[, 2] - cseg[, 1]) >= min_fragment, , drop = FALSE]

  sc <- iv_clip(pick(d, grp2, ev_contact), t0, t1)
  sc_bouts <- iv_merge(sc, tol = bout_gap_tol)
  sc_len   <- if (nrow(sc_bouts)) sc_bouts[, 2] - sc_bouts[, 1] else numeric(0)

  known <- c(unlist(lapply(dur_measures, `[[`, 2)), unlist(cnt_measures),
             ev_ignore, unname(event_aliases))
  unknown <- sort(unique(d$event[!(d$event %in% known)]))

  list(
    session = data.frame(
      file        = basename(path),
      window_end  = t1,
      as.data.frame(vals),
      alone_total     = iv_total(alone),
      grooming        = iv_total(iv_clip(pick(d, grp3, ev_groom_sum), t0, t1)),
      passive_alone   = iv_total(seg),
      passive_contact = iv_total(cseg),
      soc_bout_n      = nrow(sc_bouts),
      soc_bout_dur    = if (nrow(sc_bouts)) iv_total(sc) / nrow(sc_bouts) else 0,
      soc_bout_med    = if (length(sc_len)) median(sc_len) else 0,
      as.data.frame(cnts),
      stringsAsFactors = FALSE
    ),
    qc = cbind(
      do.call(rbind, lapply(c(grp1, grp2, grp3), function(g)
        stream_qc(d, g, sub("\\.csv$", "", basename(path)), t0, t1))),
      unknown_events = paste(unknown, collapse = "; ")
    )
  )
}


# ============================================================================
# RUN
# ============================================================================

manifest <- read.csv(file.path(data_dir, manifest_file), stringsAsFactors = FALSE)
res <- lapply(file.path(data_dir, manifest$new_filename), analyze)

measures <- do.call(rbind, lapply(res, `[[`, "session"))
qc       <- do.call(rbind, lapply(res, `[[`, "qc"))

# identity columns come from the manifest
measures <- cbind(
  manifest[, c("animal", "session", "video_code", "hemisphere", "treatment",
               "used_in")],
  measures)

measures$grooming <- measures$give_groom + measures$receive_groom

# identities that must hold exactly; checked before rounding
stopifnot(
  abs(measures$isolation + measures$proximity + measures$social_contact -
        measures$window_end) < 0.05 |
    measures$isolation + measures$proximity + measures$social_contact <
      measures$window_end,
  abs(measures$locomotion + measures$manipulation + measures$self_directed +
        measures$passive + measures$stop_group1 - measures$window_end) < 0.05,
  measures$passive - measures$passive_alone - measures$passive_contact > -0.01,
  abs(measures$soc_bout_n * measures$soc_bout_dur - measures$social_contact) < 0.01,
  !nzchar(qc$unknown_events),
  abs(qc$gap) < 0.01, abs(qc$overlap) < 0.01
)

num <- vapply(measures, is.numeric, logical(1))
measures[num] <- lapply(measures[num], round, digits)
qcn <- vapply(qc, is.numeric, logical(1))
qc[qcn] <- lapply(qc[qcn], round, digits)

measures <- measures[, names(out_names)]
names(measures) <- unname(out_names)

write.csv(measures, file.path(out_dir, out_measures), row.names = FALSE)
write.csv(qc, file.path(out_dir, out_qc), row.names = FALSE)

# one file per site, for the downstream analysis scripts
for (site in names(out_by_site)) {
  rows <- grepl(site, measures$`Used In`, fixed = TRUE)
  write.csv(measures[rows, ], file.path(out_dir, out_by_site[[site]]), row.names = FALSE)
  cat(sprintf("%-4s %2d sessions -> %s\n", site, sum(rows), out_by_site[[site]]))
}

# region comparison
bmi <- measures[measures$Treatment == "BMI", ]
bmi$Site <- bmi$`Used In`
both <- names(which(tapply(bmi$Site, bmi$Animal, function(x) length(unique(x))) == 2))
bmi <- bmi[bmi$Animal %in% both, ]
bmi <- bmi[order(bmi$Animal, bmi$Session), ]
write.csv(bmi, file.path(out_dir, out_region), row.names = FALSE)
cat(sprintf("%-4s %2d sessions (%s) -> %s\n", "GvS", nrow(bmi),
            paste(sprintf("%s %d", names(table(bmi$Site)), table(bmi$Site)), collapse = ", "),
            out_region))
