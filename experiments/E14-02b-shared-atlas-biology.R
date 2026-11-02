#!/usr/bin/env Rscript
# ==============================================================================
# E14-02b: cluster biology for the first-episode shared atlas (adapted from
# methods/41-cluster-biology-production.R). Patient-level matrices keep the
# first-season episodes as rows (as downstream scripts expect); the mapping of
# all 343 episodes onto the atlas is saved as *_every with all_p_ids.
# ==============================================================================
# Cluster Biology Prep (Production Configuration)
# ==============================================================================
# Loads the final global atlas from script 40 and computes:
#   - PCA embedding on sampled cells
#   - UMAP embedding on sampled cells
#   - Per-cluster marker means (scaled + original space)
#   - Per-patient cluster abundance data
#
# Usage:
#   Rscript methods/41-cluster-biology-production.R [atlas_file]
#
# Env overrides:
#   UMAP_N_NEIGHBORS (30)
#   UMAP_MIN_DIST (0.3)
#   UMAP_N_THREADS (5)
#   PER_CLUSTER_MAX (8000)
#   FORCE (0)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)

get_num_env <- function(var, default) {
  val <- Sys.getenv(var, unset = "")
  if (val == "") return(default)
  as.numeric(val)
}

UMAP_N_NEIGHBORS <- as.integer(get_num_env("UMAP_N_NEIGHBORS", 30))
UMAP_MIN_DIST <- get_num_env("UMAP_MIN_DIST", 0.3)
UMAP_N_THREADS <- as.integer(get_num_env("UMAP_N_THREADS", 5))
PER_CLUSTER_MAX <- as.integer(get_num_env("PER_CLUSTER_MAX", 8000))
FORCE <- as.integer(get_num_env("FORCE", 0))

source(here("methods", "00-paths.R"))

# Find atlas file
atlas_dir <- file.path(paths$results, "final_atlas_production")
if (length(args) >= 1 && file.exists(args[1])) {
  atlas_file <- args[1]
} else {
  atlas_dir <- Sys.getenv("FLU_SHARED_ATLAS_DIR", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14/shared_atlas")
  atlas_files <- list.files(atlas_dir, pattern = "^final_atlas_.*\\.RData$", full.names = TRUE)
  if (length(atlas_files) == 0) stop("No atlas file found in ", atlas_dir)
  atlas_file <- atlas_files[which.max(file.info(atlas_files)$mtime)]
}

cat("\n============================================================\n")
cat("CLUSTER BIOLOGY PREP (PRODUCTION)\n")
cat("============================================================\n")
cat("Atlas:           ", atlas_file, "\n", sep = "")
cat("UMAP neighbors:  ", UMAP_N_NEIGHBORS, "\n", sep = "")
cat("UMAP min_dist:   ", UMAP_MIN_DIST, "\n", sep = "")
cat("UMAP threads:    ", UMAP_N_THREADS, "\n", sep = "")
cat("Per-cluster max: ", PER_CLUSTER_MAX, "\n\n", sep = "")

tic_total <- Sys.time()

# Output path
out_dir <- Sys.getenv("FLU_SHARED_ATLAS_DIR", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14/shared_atlas")
param_list <- list(
  atlas = basename(atlas_file),
  umap_nn = UMAP_N_NEIGHBORS,
  umap_md = UMAP_MIN_DIST,
  per_cluster_max = PER_CLUSTER_MAX
)
PARAM_HASH <- substr(digest::digest(param_list, algo = "md5"), 1, 8)
out_file <- file.path(out_dir, paste0("biology_production_", PARAM_HASH, ".RData"))

if (file.exists(out_file) && FORCE == 0) {
  cat("Found existing output (skip): ", out_file, "\n")
  quit(status = 0)
}

# Load atlas
cat("Loading atlas...\n")
load(atlas_file)

if (is.null(sample_dt) || nrow(sample_dt) == 0) {
  stop("No sampled cells found in atlas file.")
}

cat(sprintf("Sampled cells: %d\n", nrow(sample_dt)))

# Extract marker columns from sample_dt
marker_cols <- marker_names
X_raw <- as.matrix(sample_dt[, ..marker_cols])
X_scaled <- scale(X_raw, center = scale_center, scale = scale_scale)

# Balance sample by meta cluster
cat("Balancing sample by meta cluster (max ", PER_CLUSTER_MAX, " per cluster)...\n", sep = "")
set.seed(20260222)
sample_dt[, keep := seq_len(.N) <= PER_CLUSTER_MAX, by = meta]
keep_idx <- which(sample_dt$keep)
dt_keep <- sample_dt[keep_idx]
X_scaled_keep <- X_scaled[keep_idx, , drop = FALSE]
X_raw_keep <- X_raw[keep_idx, , drop = FALSE]

cat(sprintf("After balancing: %d cells (from %d)\n", nrow(dt_keep), nrow(sample_dt)))

# PCA
cat("Computing PCA...\n")
pca <- prcomp(X_scaled_keep, center = FALSE, scale. = FALSE)
pca_dt <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# UMAP
umap_dt <- NULL
have_uwot <- requireNamespace("uwot", quietly = TRUE)
if (have_uwot) {
  cat("Computing UMAP...\n")
  set.seed(20260222)
  um <- uwot::umap(
    X_scaled_keep,
    n_neighbors = UMAP_N_NEIGHBORS,
    min_dist = UMAP_MIN_DIST,
    n_threads = UMAP_N_THREADS,
    verbose = TRUE
  )
  umap_dt <- data.table(UMAP1 = um[, 1], UMAP2 = um[, 2])
  cat("UMAP complete.\n")
} else {
  cat("Package 'uwot' not installed; skipping UMAP.\n")
}

# Per-cluster marker means (scaled space, for heatmap)
cat("Computing per-cluster marker means...\n")
means_scaled <- as.data.table(
  rowsum(X_scaled_keep, group = dt_keep$meta) /
    as.numeric(table(dt_keep$meta))
)
means_scaled[, meta := as.integer(rownames(means_scaled))]
setcolorder(means_scaled, c("meta", marker_cols))
means_long <- melt(means_scaled, id.vars = "meta",
                   variable.name = "marker", value.name = "mean_scaled")

# Per-cluster marker means (original space)
means_original <- as.data.table(
  rowsum(X_raw_keep, group = dt_keep$meta) /
    as.numeric(table(dt_keep$meta))
)
means_original[, meta := as.integer(rownames(means_original))]
setcolorder(means_original, c("meta", marker_cols))

timings_bio <- list(
  total_s = as.numeric(difftime(Sys.time(), tic_total, units = "secs"))
)

cat(sprintf("\nTotal time: %.1f seconds\n", timings_bio$total_s))

save(
  meta_counts_every, meta_props_every, clr_every, all_p_ids, cohort_all, visit_type,
  param_list,
  timings_bio,
  atlas_file,
  dt_keep,
  marker_names,
  X_scaled_keep,
  pca_dt,
  umap_dt,
  means_long,
  means_original,
  meta_counts_all,
  meta_props_all,
  clr_features,
  marker_profiles,
  common_p_ids,
  cohort_vec,
  selected_k,
  selected_res,
  file = out_file
)

cat("Saved: ", out_file, "\n", sep = "")
cat("Done.\n")
