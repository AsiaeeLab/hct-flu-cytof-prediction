# ==============================================================================
# Multi-Outcome Analysis with Covariates
# ==============================================================================
# Functions for handling multiple virus strains and adjusting for covariates

ensure_p_id_yr <- function(outcome_obj) {
  if (!is.null(outcome_obj$p_id_yr)) {
    return(outcome_obj)
  }

  if (!is.null(outcome_obj$id_yr)) {
    outcome_obj$p_id_yr <- sub("_(\\d)$", "_yr\\1", outcome_obj$id_yr)
    return(outcome_obj)
  }

  if (!is.null(outcome_obj$pid) && !is.null(outcome_obj$yr)) {
    outcome_obj$p_id_yr <- paste0(outcome_obj$pid, "_yr", outcome_obj$yr)
    return(outcome_obj)
  }

  stop("Cannot construct p_id_yr: need one of {p_id_yr}, {id_yr}, or {pid, yr} columns.")
}

# HAI titers come from the REDCap exports (experiments/E17-07-redcap-titers.R).
# They replace the titer columns of whichever counts or metadata file the caller
# read; episodes without a REDCap record get missing titers.
use_redcap_titers <- function(outcome_obj) {
  hai_file <- Sys.getenv("FLU_HAI_FILE", file.path(paths$raw, "hct_hai_redcap_2026-10-29.csv"))
  rc <- read.csv(hai_file, stringsAsFactors = FALSE)
  outcome_obj <- ensure_p_id_yr(outcome_obj)
  idx <- match(sub("_yr", "_", outcome_obj$p_id_yr), rc$id_yr)
  cols <- grep("^geo_mean_", names(outcome_obj), value = TRUE)
  stopifnot(all(cols %in% names(rc)))
  for (col in cols) outcome_obj[[col]] <- rc[[col]][idx]
  outcome_obj
}

#' Prepare outcomes for all virus strains
#'
#' Extracts log fold-change outcomes for all viruses (H1N1, H3N2, B/Victoria, B/Yamagata)
#'
#' @param outcome_obj Data frame containing outcome measurements
#' @param p_id Vector of patient IDs from the cell data
#' @param exclude_viruses Vector of virus names to exclude (e.g., "bvic")
#' @return List with outcomes for each virus
prepare_all_outcomes <- function(outcome_obj, p_id, exclude_viruses = c()) {

  outcome_obj <- use_redcap_titers(outcome_obj)

  viruses <- c("h1n1", "h3n2", "bvic", "byam")
  virus_names <- c("H1N1", "H3N2", "B/Victoria", "B/Yamagata")
  names(virus_names) <- viruses

  # Remove excluded viruses
  if (length(exclude_viruses) > 0) {
    viruses <- setdiff(viruses, exclude_viruses)
  }

  all_outcomes <- list()

  for (virus in viruses) {
    baseline_col <- paste0("geo_mean_", virus, "_v1")
    target_col <- paste0("geo_mean_", virus, "_v3")

    # Get outcomes
    outcome_v1 <- outcome_obj[[baseline_col]]
    outcome_v3 <- outcome_obj[[target_col]]
    outcome_v1[!is.na(outcome_v1) & outcome_v1 < 5] <- NA  # below 5 marks a missing measurement
    outcome_v3[!is.na(outcome_v3) & outcome_v3 < 5] <- NA

    # Calculate response as log(v3/v1) = log(v3) - log(v1)
    outcome_y <- log(outcome_v3) - log(outcome_v1)
    outcome_x <- log(outcome_v1)  # Baseline covariate

    # Check for completeness
    has_both <- !is.na(outcome_y) & !is.na(outcome_x) &
                is.finite(outcome_x) & is.finite(outcome_y)
    has_cell <- outcome_obj$p_id_yr %in% p_id
    complete_index <- has_both & has_cell

    # Extract complete cases
    outcome_y <- outcome_y[complete_index]
    outcome_x <- outcome_x[complete_index]
    patient_ids <- outcome_obj$p_id_yr[complete_index]

    patient_outcome_y <- setNames(outcome_y, patient_ids)
    patient_outcome_x <- setNames(outcome_x, patient_ids)

    all_outcomes[[virus]] <- list(
      y = patient_outcome_y,
      x = patient_outcome_x,
      name = virus_names[virus],
      n_patients = length(unique(patient_ids))
    )
  }

  return(all_outcomes)
}


#' Merge covariates from metadata
#'
#' Adds dose (HD/SD), age, and sex information to patient outcomes
#'
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata data frame
#' @return Merged data frame with covariates
load_episode_timing <- function() {
  raw_dir <- if (exists("paths") && !is.null(paths$raw)) paths$raw else
    path.expand("~/datasets/flu/raw")
  path <- Sys.getenv("FLU_EPISODE_METADATA", unset =
    file.path(raw_dir, "hsct_metadata_15JUL2025.csv"))
  if (!file.exists(path)) stop("Episode-specific transplant timing is required: ", path)
  ensure_p_id_yr(read.csv(path, stringsAsFactors = FALSE))
}

merge_covariates <- function(outcome_obj, meta_obj, episode_meta = NULL) {

  outcome_obj <- use_redcap_titers(outcome_obj)

  # Convert metadata study_id to outcome pid format.
  # outcome_obj has pid like "ASCT_<site>-<number>" (adult), "PSCT_..." (pediatric) or "ASOT_...".
  convert_study_id_to_pid <- function(study_id) {
    if (is.na(study_id) || !nzchar(study_id)) {
      return(NA_character_)
    }
    if (grepl("^ASCTF", study_id)) return(sub("^ASCTF", "ASCT_", study_id))
    if (grepl("^ASCT_", study_id)) return(study_id)
    if (grepl("^ASOTF", study_id)) return(sub("^ASOTF", "ASOT_", study_id))
    if (grepl("^ASOT_", study_id)) return(study_id)
    if (grepl("^PSCTF", study_id)) return(sub("^PSCTF", "PSCT_", study_id))
    if (grepl("^PSCT_", study_id)) return(study_id)
    paste0("PSCT_", study_id)
  }

  meta_obj$pid <- vapply(meta_obj$study_id, convert_study_id_to_pid, character(1))

  # Select relevant covariates and merge by pid (metadata is patient-level, not year-specific).
  covariates <- meta_obj[, c("pid", "intervention_recode",
                             "age_consent_years", "transplant_vaccine_months")]
  colnames(covariates) <- c("pid", "dose", "age", "time_post_transplant")

  merged <- merge(outcome_obj, covariates, by = "pid", all.x = TRUE)
  merged <- ensure_p_id_yr(merged)

  # Timing belongs to the vaccination episode, not the recipient. Keep the
  # enrollment value only for provenance; never approximate a repeat as +12 mo.
  if (is.null(episode_meta)) episode_meta <- load_episode_timing()
  episode_meta <- ensure_p_id_yr(episode_meta)
  stopifnot(!anyDuplicated(episode_meta$p_id_yr), "mo_post" %in% names(episode_meta))
  merged$time_post_transplant_enrollment <- merged$time_post_transplant
  merged$time_post_transplant <- episode_meta$mo_post[
    match(merged$p_id_yr, episode_meta$p_id_yr)]

  # Also add sex from outcome_obj if available
  if ("sex" %in% colnames(outcome_obj)) {
    merged$sex <- outcome_obj$sex[match(merged$id_yr, outcome_obj$id_yr)]
  }

  return(merged)
}


#' Extract covariate matrix for patients
#'
#' Creates a design matrix with covariates (dose, age, sex) for prediction
#'
#' @param merged_data Merged outcome + metadata data frame
#' @param patient_ids Vector of patient IDs to extract
#' @param include_dose Include dose (HD/SD) as covariate
#' @param include_age Include age as covariate
#' @param include_sex Include sex as covariate
#' @param include_time Include time post-transplant (for future use)
#' @return Matrix with covariates (patients × covariates)
extract_covariate_matrix <- function(merged_data, patient_ids,
                                     include_dose = TRUE,
                                     include_age = TRUE,
                                     include_sex = TRUE,
                                     include_time = FALSE) {

  # Filter to requested patients
  data_subset <- merged_data[match(patient_ids, merged_data$p_id_yr), ]

  # Initialize covariate matrix
  covariate_list <- list()

  # Dose (HD=1, SD=0)
  if (include_dose && "dose" %in% colnames(data_subset)) {
    dose_hd <- as.numeric(data_subset$dose == "HD")
    dose_hd[is.na(dose_hd)] <- 0  # Replace NA with 0
    covariate_list$dose_HD <- dose_hd
  }

  # Age (centered and scaled)
  if (include_age && "age" %in% colnames(data_subset)) {
    age <- data_subset$age
    # Handle missing values and zero variance
    if (any(is.na(age))) {
      age[is.na(age)] <- mean(age, na.rm = TRUE)
    }
    age_sd <- sd(age, na.rm = TRUE)
    if (age_sd > 0) {
      age_scaled <- scale(age)[,1]  # Center and scale
    } else {
      age_scaled <- age - mean(age, na.rm = TRUE)  # Just center if no variance
    }
    # Replace any remaining NAs with 0
    age_scaled[is.na(age_scaled)] <- 0
    covariate_list$age <- age_scaled
  }

  # Sex (Male=1, Female=0)
  if (include_sex && "sex" %in% colnames(data_subset)) {
    sex_male <- as.numeric(data_subset$sex == "Male")
    sex_male[is.na(sex_male)] <- 0  # Replace NA with 0
    covariate_list$sex_Male <- sex_male
  }

  # Time post-transplant (for future)
  if (include_time && "time_post_transplant" %in% colnames(data_subset)) {
    time <- data_subset$time_post_transplant
    # Fail on absent episode timing instead of silently imputing enrollment.
    if (any(!is.finite(time))) stop("Missing episode-specific transplant timing for: ",
      paste(patient_ids[!is.finite(time)], collapse = ", "))
    time_sd <- sd(time, na.rm = TRUE)
    if (time_sd > 0) {
      time_scaled <- scale(time)[,1]
    } else {
      time_scaled <- time - mean(time, na.rm = TRUE)
    }
    time_scaled[is.na(time_scaled)] <- 0
    covariate_list$time_post_transplant <- time_scaled
  }

  # Combine into matrix
  if (length(covariate_list) > 0) {
    covariate_matrix <- do.call(cbind, covariate_list)
    rownames(covariate_matrix) <- patient_ids
  } else {
    covariate_matrix <- matrix(nrow = length(patient_ids), ncol = 0)
    rownames(covariate_matrix) <- patient_ids
  }

  return(covariate_matrix)
}


#' Prepare outcome with covariates
#'
#' Extended version of prepare_outcome that includes clinical covariates
#'
#' @param outcome_obj Data frame containing outcome measurements
#' @param meta_obj Metadata data frame with covariates
#' @param outcome_col_1 Column name for baseline measurement (v1)
#' @param outcome_col_2 Column name for target measurement (v3)
#' @param p_id Vector of patient IDs from the cell data
#' @param include_dose Include dose covariate
#' @param include_age Include age covariate
#' @param include_sex Include sex covariate
#' @return List with y, x (baseline), and covariates matrix
prepare_outcome_with_covariates <- function(outcome_obj, meta_obj,
                                            outcome_col_1, outcome_col_2, p_id,
                                            include_dose = TRUE,
                                            include_age = TRUE,
                                            include_sex = TRUE,
                                            include_time = FALSE) {

  # Merge covariates first
  merged_data <- merge_covariates(outcome_obj, meta_obj)

  # Prepare basic outcomes
  outcome_v1 <- merged_data[[outcome_col_1]]
  outcome_v3 <- merged_data[[outcome_col_2]]
  # Titers below detection are recorded as 5; anything lower marks a missing
  # measurement (the older counts file coded those as 1).
  outcome_v1[!is.na(outcome_v1) & outcome_v1 < 5] <- NA
  outcome_v3[!is.na(outcome_v3) & outcome_v3 < 5] <- NA

  # Calculate response
  outcome_y <- log(outcome_v3) - log(outcome_v1)
  outcome_x <- log(outcome_v1)

  # Check completeness
  has_both <- !is.na(outcome_y) & !is.na(outcome_x) &
              is.finite(outcome_x) & is.finite(outcome_y)
  has_cell <- merged_data$p_id_yr %in% p_id
  complete_index <- has_both & has_cell

  # Extract complete cases
  outcome_y <- outcome_y[complete_index]
  outcome_x <- outcome_x[complete_index]
  patient_ids <- merged_data$p_id_yr[complete_index]

  patient_outcome_y <- setNames(outcome_y, patient_ids)
  patient_outcome_x <- setNames(outcome_x, patient_ids)

  # Extract covariates for these patients
  covariate_matrix <- extract_covariate_matrix(
    merged_data, patient_ids,
    include_dose = include_dose,
    include_age = include_age,
    include_sex = include_sex,
    include_time = include_time
  )

  return(list(
    y = patient_outcome_y,
    x = patient_outcome_x,
    covariates = covariate_matrix
  ))
}


#' Prepare outcome with covariates - ALTERNATIVE: predict log(v3)
#'
#' Alternative formulation where response is log(v3) and log(v1) is used as feature.
#'
#' @param outcome_obj Data frame containing outcome measurements
#' @param meta_obj Metadata data frame with covariates
#' @param outcome_col_1 Column name for baseline measurement (v1)
#' @param outcome_col_2 Column name for target measurement (v3)
#' @param p_id Vector of patient IDs from the cell data
#' @param include_dose Include dose covariate
#' @param include_age Include age covariate
#' @param include_sex Include sex covariate
#' @param include_time Include time covariate
#' @return List with y (log(v3)), x (log(v1)), and covariates matrix
prepare_outcome_with_covariates_absolute <- function(outcome_obj, meta_obj,
                                                     outcome_col_1, outcome_col_2, p_id,
                                                     include_dose = TRUE,
                                                     include_age = TRUE,
                                                     include_sex = TRUE,
                                                     include_time = FALSE) {

  # Merge covariates first
  merged_data <- merge_covariates(outcome_obj, meta_obj)

  # Prepare basic outcomes
  outcome_v1 <- merged_data[[outcome_col_1]]
  outcome_v3 <- merged_data[[outcome_col_2]]
  # Titers below detection are recorded as 5; anything lower marks a missing
  # measurement (the older counts file coded those as 1).
  outcome_v1[!is.na(outcome_v1) & outcome_v1 < 5] <- NA
  outcome_v3[!is.na(outcome_v3) & outcome_v3 < 5] <- NA

  # Response is log(v3) directly
  outcome_y <- log(outcome_v3)
  # Baseline is log(v1) - used as feature
  outcome_x <- log(outcome_v1)

  # Check completeness
  has_both <- !is.na(outcome_y) & !is.na(outcome_x) &
              is.finite(outcome_x) & is.finite(outcome_y)
  has_cell <- merged_data$p_id_yr %in% p_id
  complete_index <- has_both & has_cell

  # Extract complete cases
  outcome_y <- outcome_y[complete_index]
  outcome_x <- outcome_x[complete_index]
  patient_ids <- merged_data$p_id_yr[complete_index]

  patient_outcome_y <- setNames(outcome_y, patient_ids)
  patient_outcome_x <- setNames(outcome_x, patient_ids)

  # Extract covariates for these patients
  covariate_matrix <- extract_covariate_matrix(
    merged_data, patient_ids,
    include_dose = include_dose,
    include_age = include_age,
    include_sex = include_sex,
    include_time = include_time
  )

  return(list(
    y = patient_outcome_y,
    x = patient_outcome_x,
    covariates = covariate_matrix
  ))
}


#' Build patient features with covariates
#'
#' Extended version that adds clinical covariates to cluster features
#'
#' @param p_id Vector of patient IDs for each cell
#' @param clusters Vector of cluster assignments for each cell
#' @param baseline_outcome Named vector of baseline outcomes per patient
#' @param covariates Matrix of covariates (patients × covariates)
#' @param cluster_names Optional vector of cluster names
#' @return Matrix with patients as rows, clusters + baseline + covariates as columns
build_patient_features_with_covariates <- function(p_id, clusters, baseline_outcome,
                                                   covariates = NULL,
                                                   cluster_names = NULL) {

  # Build cluster features
  feature_df <- as.data.frame.matrix(table(p_id, clusters))
  feature_df <- log(feature_df + 1)

  if (!is.null(cluster_names)) {
    colnames(feature_df) <- cluster_names
  } else {
    colnames(feature_df) <- paste0("cluster_", 1:ncol(feature_df))
  }

  # Convert to matrix
  feature_matrix <- as.matrix(feature_df)
  patient_ids <- rownames(feature_matrix)

  # Add baseline
  feature_matrix <- cbind(feature_matrix, baseline = baseline_outcome[patient_ids])

  # Add covariates if provided
  if (!is.null(covariates) && nrow(covariates) > 0) {
    # Ensure covariate rows match feature matrix rows
    matched_covariates <- covariates[patient_ids, , drop = FALSE]
    feature_matrix <- cbind(feature_matrix, matched_covariates)
  }

  return(feature_matrix)
}


#' Summary of outcomes across viruses
#'
#' @param all_outcomes List from prepare_all_outcomes
#' @return Data frame summarizing each virus
summarize_outcomes <- function(all_outcomes) {
  summary_df <- data.frame(
    virus = names(all_outcomes),
    virus_name = sapply(all_outcomes, function(x) x$name),
    n_patients = sapply(all_outcomes, function(x) x$n_patients),
    mean_log_fold_change = sapply(all_outcomes, function(x) mean(x$y, na.rm = TRUE)),
    sd_log_fold_change = sapply(all_outcomes, function(x) sd(x$y, na.rm = TRUE))
  )

  return(summary_df)
}
