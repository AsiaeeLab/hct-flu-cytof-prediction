#!/usr/bin/env Rscript
# E14-02: Shared interpretation atlas built from FIRST vaccination episodes only.
#
# Adapted from methods/40-build-final-atlas-production.R. Differences:
#   - k-means micro-clusters, Leiden graph and the model fit use the 262
#     first-season episodes only; second-season episodes are never used to
#     build the atlas or fit a model.
#   - All 343 episodes are mapped onto the finished atlas (counts, CLR,
#     proportions), so second-season compositions can be compared with first.
#   - The number of populations is not fixed at 33. TARGET_K is set from the
#     first-episode nested CV (E14-03 rule) and passed in for stage 2.
#   - STAGE=1 builds the geometry, Leiden sweep, cell assignment and cell sample
#     (independent of k); STAGE=2 selects the resolution for TARGET_K and
#     finishes (meta counts, CLR, model fit, transfer, marker profiles).
suppressPackageStartupMessages({
  library(here)
  library(digest)
  library(data.table)
  library(RANN)
  library(glmnet)
})

get_num_env <- function(var, default) {
  val <- Sys.getenv(var, unset = "")
  if (val == "") return(default)
  as.numeric(val)
}

N_CORES      <- as.integer(get_num_env("N_CORES", 5))
TARGET_K     <- as.integer(get_num_env("TARGET_K", NA))
STAGE        <- as.integer(get_num_env("STAGE", 1))
N_MICRO      <- as.integer(get_num_env("N_MICRO", 5000))
KMEANS_FIT   <- as.integer(get_num_env("KMEANS_FIT_CELLS", 1000000))
SAMPLE_CELLS <- as.integer(get_num_env("SAMPLE_CELLS", 200000))
FORCE        <- as.integer(get_num_env("FORCE", 0))
SEED_BASE    <- 2027L
ALPHA        <- 0.9
LAMBDA_RULE  <- "lambda.min"
VIRUSES      <- c("h1n1", "h3n2", "bvic")
N_NEIGHBORS  <- 30L
CHUNK_SIZE   <- 300000L
RES_GRID     <- seq(0.3, 6.0, by = 0.05)  # Fine grid to hit target k

cat("\n============================================================\n")
cat("BUILD FINAL GLOBAL ATLAS (PRODUCTION)\n")
cat("============================================================\n")
cat("Target k:       ", TARGET_K, "\n", sep = "")
cat("N micro:        ", N_MICRO, "\n", sep = "")
cat("K-means fit:    ", KMEANS_FIT, "\n", sep = "")
cat("Sample cells:   ", SAMPLE_CELLS, "\n", sep = "")
cat("Cores:          ", N_CORES, "\n", sep = "")
cat("Alpha:          ", ALPHA, "\n", sep = "")
cat("Lambda rule:    ", LAMBDA_RULE, "\n", sep = "")
cat("Viruses:        ", paste(VIRUSES, collapse = ", "), "\n", sep = "")
cat("Force:          ", FORCE, "\n\n", sep = "")

tic_total <- Sys.time()

# ==============================================================================
# Source helpers
# ==============================================================================
source(here("methods", "00-paths.R"))
source(here("methods", "10-two-stage-clustering.R"))
source(here("methods", "35-per-fold-atlas-cv.R"))

# ==============================================================================
# Output path
# ==============================================================================
out_dir <- Sys.getenv("FLU_SHARED_ATLAS_DIR", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14/shared_atlas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stage1_file <- file.path(out_dir, "stage1_geometry.RData")

param_list <- list(
  target_k = TARGET_K, n_micro = N_MICRO, kmeans_fit = KMEANS_FIT,
  sample_cells = SAMPLE_CELLS, alpha = ALPHA, lambda_rule = LAMBDA_RULE,
  viruses = VIRUSES, n_neighbors = N_NEIGHBORS, seed_base = SEED_BASE
)
PARAM_HASH <- substr(digest::digest(param_list, algo = "md5"), 1, 8)
out_file <- file.path(out_dir, paste0("final_atlas_", PARAM_HASH, ".RData"))

if (STAGE == 2 && file.exists(out_file) && FORCE == 0) {
  cat("Found existing output (skip): ", out_file, "\n", sep = "")
  quit(status = 0)
}

# ==============================================================================
# Load data
# ==============================================================================
cat("Loading cell data...\n")
load(file.path(paths$raw, "hsct_transplant_object_harmony_final_small.RData"))
full_data <- obj$batch_correction$data
full_p_id <- obj$source
marker_names <- colnames(full_data)
rm(obj); gc()

cat(sprintf("Cells: %d x %d markers\n", nrow(full_data), ncol(full_data)))
cat(sprintf("Unique p_id_yr: %d\n", length(unique(full_p_id))))

cat("Loading outcomes and metadata...\n")
outcome_obj <- read.csv(file.path(paths$raw, "transplant_cytof_counts_hai_12MAY2024.csv"),
                        stringsAsFactors = FALSE)
meta_obj <- read.csv(file.path(paths$raw, "hct_limited-dataset_2025-07-09.csv"),
                     stringsAsFactors = FALSE)

# ==============================================================================
# Prepare outcomes (same as production pipeline)
# ==============================================================================
cat("Preparing outcomes...\n")
if (!exists("prepare_outcome_with_covariates")) {
  source(here("methods", "08-multi-outcome-covariates.R"))
}

outcomes <- list()
for (v in VIRUSES) {
  baseline_col <- paste0("geo_mean_", v, "_v1")
  target_col <- paste0("geo_mean_", v, "_v3")
  outcomes[[v]] <- prepare_outcome_with_covariates(
    outcome_obj, meta_obj,
    outcome_col_1 = baseline_col,
    outcome_col_2 = target_col,
    p_id = full_p_id,
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = TRUE
  )
}

common_p_ids <- Reduce(intersect, lapply(outcomes, function(o) names(o$y)))
cat(sprintf("Common observations: %d\n", length(common_p_ids)))
all_p_ids <- common_p_ids
episode_meta <- read.csv(file.path(paths$raw, "hsct_metadata_15JUL2025.csv"), stringsAsFactors = FALSE)
episode_meta$p_id_yr <- sub("_([0-9]+)$", "_yr\\1", episode_meta$id_yr)
visit_type <- setNames(episode_meta$visit_type, episode_meta$p_id_yr)[all_p_ids]
stopifnot(!anyNA(visit_type))
common_p_ids <- all_p_ids[visit_type == "first"]
cat(sprintf("First-season episodes used to build the atlas: %d (of %d)\n", length(common_p_ids), length(all_p_ids)))
cohort_all <- infer_cohort_from_p_id(all_p_ids)

cohort_vec <- infer_cohort_from_p_id(common_p_ids)
covariates <- outcomes[[VIRUSES[1]]]$covariates[common_p_ids, , drop = FALSE]

# Outcome matrix (n x 3)
Y_mat <- sapply(VIRUSES, function(v) outcomes[[v]]$y[common_p_ids])
colnames(Y_mat) <- VIRUSES
rownames(Y_mat) <- common_p_ids

if (STAGE == 1) {
# ==============================================================================
# Step 1: Build atlas on ALL patients
# ==============================================================================
cat("\n--- Step 1: Build atlas (k-means on all patients) ---\n")
tic_atlas <- Sys.time()

atlas <- build_atlas_for_fold(
  full_data = full_data,
  full_p_id = full_p_id,
  train_patient_ids = common_p_ids,
  n_micro = N_MICRO,
  n_neighbors = N_NEIGHBORS,
  kmeans_fit_n_cells = KMEANS_FIT,
  kmeans_nstart = 1,
  kmeans_iter_max = 50,
  seed = SEED_BASE,
  verbose = TRUE
)

centroids <- atlas$centroids
scale_center <- atlas$scale_center
scale_scale <- atlas$scale_scale
centroid_graph <- atlas$centroid_graph

atlas_s <- as.numeric(difftime(Sys.time(), tic_atlas, units = "secs"))
cat(sprintf("Atlas built in %.1f seconds\n", atlas_s))

# ==============================================================================
# Step 2: Sweep Leiden resolutions to find target k
# ==============================================================================
cat("\n--- Step 2: Leiden resolution sweep (target k=", TARGET_K, ") ---\n", sep = "")
tic_leiden <- Sys.time()

leiden_results <- lapply(RES_GRID, function(res) {
  two_stage_leiden_micro_to_meta(centroid_graph, resolution = res, seed = 123)
})
k_found_vec <- sapply(leiden_results, function(m) m$k_found)
names(k_found_vec) <- RES_GRID

leiden_s <- as.numeric(difftime(Sys.time(), tic_leiden, units = "secs"))
cat(sprintf("Leiden sweep: %.1f seconds\n", leiden_s))

# ==============================================================================
# Step 3: Assign ALL cells to micro-clusters, compute patient counts
# ==============================================================================
cat("\n--- Step 3: Assign all cells & compute patient counts ---\n")
tic_assign <- Sys.time()

micro_counts_all <- compute_patient_micro_counts(
  full_data = full_data,
  full_p_id = full_p_id,
  patient_ids = all_p_ids,
  centroids = centroids,
  scale_center = scale_center,
  scale_scale = scale_scale,
  chunk_size = CHUNK_SIZE,
  verbose = TRUE
)

assign_s <- as.numeric(difftime(Sys.time(), tic_assign, units = "secs"))
cat(sprintf("Cell assignment: %.1f seconds\n", assign_s))

# ==============================================================================
# Step 6: Sample cells for downstream visualization
# ==============================================================================
cat("\n--- Step 6: Sample cells for visualization ---\n")
tic_sample <- Sys.time()

sample_dt <- NULL
X_raw_sample <- NULL
if (SAMPLE_CELLS > 0) {
  set.seed(20260222)
  idx_all <- which(full_p_id %in% common_p_ids)
  n_sample <- min(SAMPLE_CELLS, length(idx_all))
  sel <- sample(idx_all, size = n_sample)

  X_raw_sample <- full_data[sel, , drop = FALSE]
  X_scaled <- scale(X_raw_sample, center = scale_center, scale = scale_scale)
  micro_assign <- RANN::nn2(data = centroids, query = X_scaled, k = 1)$nn.idx[, 1]

  sample_dt <- data.table(
    patient_id = full_p_id[sel],
    cohort = cohort_vec[full_p_id[sel]],
    micro = micro_assign,
    meta = NA_integer_
  )
  # Attach raw marker values as columns
  for (j in seq_along(marker_names)) {
    set(sample_dt, j = marker_names[j], value = X_raw_sample[, j])
  }

  cat(sprintf("Sampled %d cells from %d patients\n",
              n_sample, length(unique(sample_dt$patient_id))))
}

sample_s <- as.numeric(difftime(Sys.time(), tic_sample, units = "secs"))
cat(sprintf("Cell sampling: %.1f seconds\n", sample_s))

save(centroids, scale_center, scale_scale, centroid_graph, leiden_results, k_found_vec, micro_counts_all,
     sample_dt, marker_names, all_p_ids, common_p_ids, atlas_s, leiden_s, assign_s, sample_s, file = stage1_file)
cat("Stage 1 saved: ", stage1_file, "\n"); quit(status = 0)
}
load(stage1_file)
# The geometry (stage 1) is built from cells only; outcome-based steps use the
# episodes with complete titers (methods/08 treats titers coded as 1 as missing).
common_p_ids <- intersect(common_p_ids, rownames(Y_mat))
cat(sprintf("Episodes with complete titers for the outcome-based steps: %d\n", length(common_p_ids)))
stopifnot(is.finite(TARGET_K))
# Find resolution closest to target k
dist_to_target <- abs(k_found_vec - TARGET_K)
best_idx <- which.min(dist_to_target)
selected_res <- RES_GRID[best_idx]
selected_k <- k_found_vec[best_idx]
micro_to_meta <- leiden_results[[best_idx]]$micro_to_meta

cat(sprintf("Resolution sweep: %d settings tested\n", length(RES_GRID)))
cat(sprintf("k range found: %d to %d\n", min(k_found_vec), max(k_found_vec)))
cat(sprintf("Selected: resolution=%.2f -> k=%d (target was %d)\n",
            selected_res, selected_k, TARGET_K))

# Show the k landscape around the target
near_target <- which(abs(k_found_vec - TARGET_K) <= 5)
if (length(near_target) > 0) {
  cat("Resolutions near target:\n")
  for (i in near_target) {
    flag <- if (i == best_idx) " <-- SELECTED" else ""
    cat(sprintf("  res=%.2f -> k=%d%s\n", RES_GRID[i], k_found_vec[i], flag))
  }
}

# Optional relabeling: population IDs carried from the all-episode atlas
# (best marker-profile match keeps the old ID; additional splits get 34+).
RELABEL_FILE <- file.path(out_dir, "relabel_map.csv")
if (file.exists(RELABEL_FILE)) {
  rl <- read.csv(RELABEL_FILE)
  stopifnot(setequal(rl$from, unique(micro_to_meta)))
  micro_to_meta <- rl$to[match(micro_to_meta, rl$from)]
  cat("Applied population relabeling from", RELABEL_FILE, "\n")
}
# Column j of the population count matrix must be population j. The shared
# helper orders populations by first appearance (which equals their labels
# for unrelabeled Leiden output); after relabeling, order them by label.
two_stage_collapse_micro_to_meta_counts <- function(micro_counts, micro_to_meta) {
  mc <- t(rowsum(t(micro_counts), group = micro_to_meta, reorder = TRUE))
  storage.mode(mc) <- "integer"
  colnames(mc) <- paste0("cluster_", sort(unique(micro_to_meta)))
  stopifnot(identical(sort(unique(micro_to_meta)), seq_along(unique(micro_to_meta))))
  mc
}
sample_dt[, meta := micro_to_meta[micro]]
# Collapse to meta-cluster counts
meta_counts_all <- two_stage_collapse_micro_to_meta_counts(micro_counts_all, micro_to_meta)
cat(sprintf("Meta-cluster counts: %d patients x %d clusters\n",
            nrow(meta_counts_all), ncol(meta_counts_all)))

# ==============================================================================
# Step 4: CLR transform
# ==============================================================================
cat("\n--- Step 4: CLR feature transform ---\n")
clr_features <- transform_counts(meta_counts_all, method = "clr")
cat(sprintf("CLR features: %d patients x %d features\n",
            nrow(clr_features), ncol(clr_features)))

# Also compute proportions for reference
meta_props_all <- sweep(meta_counts_all, 1, pmax(1, rowSums(meta_counts_all)), "/")

# ==============================================================================
# Step 5: Fit final EN model on ALL patients
# ==============================================================================
cat("\n--- Step 5: Fit final multitask EN model ---\n")
tic_model <- Sys.time()

# Build feature matrices (matching production pipeline)
base_train_mat <- sapply(VIRUSES, function(v) outcomes[[v]]$x[common_p_ids])
colnames(base_train_mat) <- paste0("baseline_", VIRUSES)

X_base <- as.matrix(cbind(base_train_mat, covariates))
X_clust <- as.matrix(cbind(clr_features[common_p_ids, , drop = FALSE], X_base))

# OLS baseline (deterministic)
cat("  Fitting OLS baseline...\n")
baseline_fit <- fit_ols_baseline_multi(X_base, Y_mat, X_base)
baseline_pred <- baseline_fit$pred_train  # in-sample predictions

# Multitask EN with cv.glmnet for lambda selection
cat("  Fitting multitask EN (alpha=", ALPHA, ")...\n", sep = "")
if (N_CORES > 1) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(N_CORES)
  doParallel::registerDoParallel(cl)
}

# Grouped foldid for lambda CV
foldid <- compute_grouped_foldid(common_p_ids, n_folds = 5L, seed = SEED_BASE + 999)

set.seed(SEED_BASE + 1)
cv_fit <- cv.glmnet(
  x = X_clust,
  y = Y_mat,
  family = "mgaussian",
  alpha = ALPHA,
  standardize = TRUE,
  foldid = foldid,
  parallel = N_CORES > 1
)

s_val <- if (LAMBDA_RULE == "lambda.min") cv_fit$lambda.min else cv_fit$lambda.1se

# Extract coefficients
en_coefs <- coef(cv_fit, s = s_val)
en_pred <- as_pred_matrix(
  predict(cv_fit, newx = X_clust, s = s_val),
  VIRUSES
)

if (N_CORES > 1) {
  try(parallel::stopCluster(cl), silent = TRUE)
}

model_s <- as.numeric(difftime(Sys.time(), tic_model, units = "secs"))
cat(sprintf("Model fitting: %.1f seconds\n", model_s))

# Summarize coefficients
cat("\n  Non-zero coefficients per virus:\n")
for (v in VIRUSES) {
  beta <- en_coefs[[v]]
  nz <- sum(beta[-1] != 0)  # exclude intercept
  cat(sprintf("    %s: %d non-zero (of %d)\n", v, nz, length(beta) - 1))
}

# In-sample performance (for reference only — unbiased estimate comes from nested CV)
cat("\n  In-sample RMSE (for reference only, NOT unbiased):\n")
for (v in VIRUSES) {
  base_rmse <- sqrt(mean((baseline_pred[, v] - Y_mat[, v])^2))
  en_rmse <- sqrt(mean((en_pred[, v] - Y_mat[, v])^2))
  cat(sprintf("    %s: baseline=%.3f, EN=%.3f (improvement=%.1f%%)\n",
              v, base_rmse, en_rmse, 100 * (1 - en_rmse / base_rmse)))
}

# ==============================================================================
# Step 7: Cross-cohort transfer diagnostic
# ==============================================================================
cat("\n--- Step 7: Cross-cohort transfer diagnostic ---\n")
tic_transfer <- Sys.time()

transfer_dt <- NULL
tryCatch({
  peds_ids <- common_p_ids[cohort_vec[common_p_ids] == "peds"]
  adult_ids <- common_p_ids[cohort_vec[common_p_ids] == "adult"]

  cat(sprintf("  Peds: %d, Adult: %d\n", length(peds_ids), length(adult_ids)))

  if (length(peds_ids) >= 30 && length(adult_ids) >= 30) {
    transfer_results <- list()

    for (direction in c("peds_to_adult", "adult_to_peds")) {
      if (direction == "peds_to_adult") {
        train_ids <- peds_ids; test_ids <- adult_ids
      } else {
        train_ids <- adult_ids; test_ids <- peds_ids
      }

      for (v in VIRUSES) {
        # Baseline OLS
        X_base_tr <- as.matrix(cbind(
          baseline = outcomes[[v]]$x[train_ids],
          covariates[train_ids, , drop = FALSE]
        ))
        X_base_te <- as.matrix(cbind(
          baseline = outcomes[[v]]$x[test_ids],
          covariates[test_ids, , drop = FALSE]
        ))
        y_tr <- outcomes[[v]]$y[train_ids]
        y_te <- outcomes[[v]]$y[test_ids]

        base_fit_xfer <- fit_ols_baseline(X_base_tr, y_tr, X_base_te)
        base_mse <- mean((base_fit_xfer$pred_test - y_te)^2)

        # Cluster model (single-task EN for simplicity in transfer)
        X_clust_tr <- as.matrix(cbind(
          clr_features[train_ids, , drop = FALSE],
          X_base_tr
        ))
        X_clust_te <- as.matrix(cbind(
          clr_features[test_ids, , drop = FALSE],
          X_base_te
        ))

        set.seed(SEED_BASE + 500)
        cv_xfer <- cv.glmnet(X_clust_tr, y_tr, alpha = ALPHA, nfolds = 5)
        s_xfer <- cv_xfer$lambda.min
        pred_xfer <- predict(cv_xfer, newx = X_clust_te, s = s_xfer)
        clust_mse <- mean((as.numeric(pred_xfer) - y_te)^2)

        improve_pct <- 100 * (clust_mse / base_mse - 1)

        transfer_results[[length(transfer_results) + 1]] <- data.table(
          direction = direction,
          virus = v,
          baseline_mse = base_mse,
          clusters_mse = clust_mse,
          improve_pct = improve_pct,
          n_train = length(train_ids),
          n_test = length(test_ids)
        )
      }
    }

    transfer_dt <- rbindlist(transfer_results)
    cat("  Transfer results:\n")
    print(transfer_dt[, .(direction, virus,
                          improve_pct = round(improve_pct, 1),
                          n_train, n_test)])
  } else {
    cat("  Skipping transfer: insufficient cohort sizes\n")
  }
}, error = function(e) {
  cat("  Transfer diagnostic failed: ", conditionMessage(e), "\n")
})

transfer_s <- as.numeric(difftime(Sys.time(), tic_transfer, units = "secs"))
cat(sprintf("Transfer diagnostic: %.1f seconds\n", transfer_s))

# ==============================================================================
# Step 8: Compute per-cluster marker profiles (centroid means in original space)
# ==============================================================================
cat("\n--- Step 8: Cluster marker profiles ---\n")

# The centroids are in scaled space; transform back to original
centroid_original <- sweep(centroids, 2, scale_scale, "*")
centroid_original <- sweep(centroid_original, 2, scale_center, "+")
colnames(centroid_original) <- marker_names

# Per meta-cluster mean (weighted by micro-cluster size)
micro_sizes <- colSums(micro_counts_all[common_p_ids, , drop = FALSE])  # first-season cells per micro-cluster
meta_cluster_ids <- sort(unique(micro_to_meta))
n_meta <- length(meta_cluster_ids)

marker_profiles <- matrix(0, nrow = n_meta, ncol = ncol(centroids))
colnames(marker_profiles) <- marker_names
rownames(marker_profiles) <- paste0("cluster_", meta_cluster_ids)

for (mc in seq_along(meta_cluster_ids)) {
  micro_idx <- which(micro_to_meta == meta_cluster_ids[mc])
  w <- micro_sizes[micro_idx]
  if (sum(w) > 0) {
    marker_profiles[mc, ] <- colSums(centroid_original[micro_idx, , drop = FALSE] * w) / sum(w)
  }
}

cat(sprintf("Marker profiles: %d clusters x %d markers\n",
            nrow(marker_profiles), ncol(marker_profiles)))

# ==============================================================================
# Save everything
# ==============================================================================
total_s <- as.numeric(difftime(Sys.time(), tic_total, units = "secs"))

timings <- list(
  atlas_s = atlas_s,
  leiden_s = leiden_s,
  assign_s = assign_s,
  model_s = model_s,
  sample_s = sample_s,
  transfer_s = transfer_s,
  total_s = total_s
)

cat(sprintf("\n=== TOTAL TIME: %.1f seconds (%.1f minutes) ===\n",
            total_s, total_s / 60))

# Patient-level matrices: rows = first-season episodes (as downstream scripts
# expect); the mapping of all 343 episodes is kept as *_every.
meta_counts_every <- meta_counts_all; meta_props_every <- meta_props_all; clr_every <- clr_features
for (m in list(meta_counts_every, meta_props_every, clr_every)) if (!is.null(rownames(m))) stopifnot(identical(rownames(m), all_p_ids))
rownames(meta_counts_every) <- rownames(meta_props_every) <- rownames(clr_every) <- all_p_ids
meta_counts_all <- meta_counts_every[common_p_ids, , drop = FALSE]
meta_props_all <- meta_props_every[common_p_ids, , drop = FALSE]
clr_features <- clr_every[common_p_ids, , drop = FALSE]
cat("Saving to:", out_file, "\n")

save(
  # Parameters
  param_list, timings,
  VIRUSES, ALPHA, LAMBDA_RULE, TARGET_K,
  selected_res, selected_k,

  # Atlas geometry
  centroids, scale_center, scale_scale, centroid_graph,
  micro_to_meta, marker_names,

  # Resolution sweep
  RES_GRID, k_found_vec,

  # Patient-level data
  common_p_ids, cohort_vec, covariates, Y_mat, outcomes,
  all_p_ids, cohort_all, visit_type, meta_counts_every, meta_props_every, clr_every,
  micro_counts_all, meta_counts_all, meta_props_all, clr_features,

  # Model outputs
  cv_fit, en_coefs, en_pred, baseline_pred,
  X_base, X_clust,

  # Visualization data
  sample_dt,

  # Marker profiles
  marker_profiles,

  # Transfer diagnostic
  transfer_dt,

  file = out_file
)

cat("Done.\n")
cat(sprintf("Output: %s (%.1f MB)\n",
            out_file,
            file.info(out_file)$size / 1e6))
