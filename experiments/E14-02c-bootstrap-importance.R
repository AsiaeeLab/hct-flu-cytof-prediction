#!/usr/bin/env Rscript
# E14-02c: bootstrap importance on the first-episode shared atlas (copy of methods/42 with a new output directory).
# ==============================================================================
# Bootstrap EN Coefficient Stability (Production Configuration)
# ==============================================================================
# Loads the final global atlas from script 40 and bootstraps:
#   - EN coefficient stability (which clusters are consistently selected)
#   - Selection frequency and sign consistency per feature
#
# No RF (dropped per analysis showing EN dominates).
#
# Usage:
#   Rscript methods/42-bootstrap-importance-production.R [atlas_file]
#
# Env overrides:
#   N_BOOT (200)
#   SEED_BOOT (202701)
#   N_CORES (5)
#   FORCE (0)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(glmnet)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)

get_num_env <- function(var, default) {
  val <- Sys.getenv(var, unset = "")
  if (val == "") return(default)
  as.numeric(val)
}

N_BOOT <- as.integer(get_num_env("N_BOOT", 200))
SEED_BOOT <- as.integer(get_num_env("SEED_BOOT", 202701))
N_CORES <- as.integer(get_num_env("N_CORES", 5))
FORCE <- as.integer(get_num_env("FORCE", 0))

source(here("methods", "00-paths.R"))
source(here("methods", "35-per-fold-atlas-cv.R"))  # for compute_grouped_foldid, fit helpers

# Find atlas file
atlas_dir <- file.path(paths$results, "final_atlas_production")
if (length(args) >= 1 && file.exists(args[1])) {
  atlas_file <- args[1]
} else if (nzchar(Sys.getenv("FLU_PROD_ATLAS"))) {
  atlas_file <- Sys.getenv("FLU_PROD_ATLAS")
  if (!file.exists(atlas_file)) stop("FLU_PROD_ATLAS does not exist: ", atlas_file)
} else if (nzchar(Sys.getenv("FLU_RESULTS_ROOT"))) {
  # A results root without an explicit atlas would otherwise fall back to an
  # older production atlas built on a different episode set.
  stop("FLU_RESULTS_ROOT is set but FLU_PROD_ATLAS is not")
} else {
  atlas_files <- list.files(atlas_dir, pattern = "^final_atlas_.*\\.RData$", full.names = TRUE)
  if (length(atlas_files) == 0) stop("No atlas file found in ", atlas_dir)
  atlas_file <- atlas_files[which.max(file.info(atlas_files)$mtime)]
}

cat("\n============================================================\n")
cat("BOOTSTRAP EN IMPORTANCE (PRODUCTION)\n")
cat("============================================================\n")
cat("Atlas:      ", atlas_file, "\n", sep = "")
cat("Bootstraps: ", N_BOOT, "\n", sep = "")
cat("Seed:       ", SEED_BOOT, "\n", sep = "")
cat("Cores:      ", N_CORES, "\n", sep = "")
cat("Force:      ", FORCE, "\n\n", sep = "")

tic_total <- Sys.time()

# Output path
out_dir <- Sys.getenv("FLU_SHARED_ATLAS_DIR", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14/shared_atlas")
# OUTCOME=composite fits a single-outcome elastic net to the composite outcome
# (mean of the three antigen log fold-rises); coefficients are stored with
# virus = "comp". The default is the three-antigen multitask model.
OUTCOME <- Sys.getenv("OUTCOME", "multitask")
param_out <- list(
  atlas = basename(atlas_file),
  n_boot = N_BOOT,
  seed_boot = SEED_BOOT,
  resampling_unit = "base_recipient"
)
if (OUTCOME == "composite") param_out$outcome <- "composite"
PARAM_HASH <- substr(digest::digest(param_out, algo = "md5"), 1, 8)
out_file <- file.path(out_dir, paste0("bootstrap_importance_", PARAM_HASH, ".RData"))

if (file.exists(out_file) && FORCE == 0) {
  cat("Found existing output (skip): ", out_file, "\n")
  quit(status = 0)
}

# Load atlas
cat("Loading atlas...\n")
load(atlas_file)

VIRUSES_BOOT <- VIRUSES
if (OUTCOME == "composite") {
  Y_mat <- matrix(rowMeans(Y_mat), ncol = 1, dimnames = list(rownames(Y_mat), "comp"))
  VIRUSES_BOOT <- "comp"
}
FAMILY <- if (OUTCOME == "composite") "gaussian" else "mgaussian"
YFIT <- function(Y) if (OUTCOME == "composite") Y[, 1] else Y
patient_ids <- common_p_ids
n_pat <- length(patient_ids)
base_ids <- base_pid_from_p_id_yr(patient_ids)
recipient_rows <- split(seq_len(n_pat), base_ids)

# First: select lambda per virus via cv.glmnet on all patients (multitask)
cat("Selecting EN lambda (multitask cv.glmnet)...\n")
if (N_CORES > 1) {
  suppressPackageStartupMessages(library(doParallel))
  cl <- parallel::makeCluster(N_CORES)
  doParallel::registerDoParallel(cl)
}

foldid <- compute_grouped_foldid(patient_ids, n_folds = 5L, seed = SEED_BOOT)

set.seed(SEED_BOOT)
cvfit_ref <- cv.glmnet(
  x = X_clust,
  y = YFIT(Y_mat),
  family = FAMILY,
  alpha = ALPHA,
  standardize = TRUE,
  foldid = foldid,
  parallel = N_CORES > 1
)
lambda_ref <- cvfit_ref$lambda.min
cat(sprintf("Reference lambda.min = %.6f\n", lambda_ref))

# Bootstrap loop
cat("\nRunning ", N_BOOT, " bootstraps...\n", sep = "")
coef_rows <- list()
perf_rows <- list()

set.seed(SEED_BOOT)
for (b in seq_len(N_BOOT)) {
  if (b %% max(1, floor(N_BOOT / 10)) == 0) {
    cat("  bootstrap ", b, " / ", N_BOOT, "\n", sep = "")
  }

  draw <- sample(names(recipient_rows), length(recipient_rows), replace = TRUE)
  idx_boot <- unlist(recipient_rows[draw], use.names = FALSE)
  inbag <- rep(FALSE, n_pat)
  inbag[idx_boot] <- TRUE
  oob_idx <- which(!inbag)
  if (length(oob_idx) < 10) next

  # Fit multitask EN with fixed lambda
  X_boot <- X_clust[idx_boot, , drop = FALSE]
  Y_boot <- Y_mat[idx_boot, , drop = FALSE]
  X_oob <- X_clust[oob_idx, , drop = FALSE]
  Y_oob <- Y_mat[oob_idx, , drop = FALSE]

  fit <- glmnet(
    x = X_boot, y = YFIT(Y_boot),
    family = FAMILY,
    alpha = ALPHA,
    standardize = TRUE,
    lambda = lambda_ref
  )

  # Extract coefficients per virus
  for (v in VIRUSES_BOOT) {
    beta <- if (OUTCOME == "composite") as.matrix(coef(fit))[, 1] else as.matrix(coef(fit)[[v]])[, 1]
    coef_rows[[length(coef_rows) + 1]] <- data.table(
      bootstrap = b,
      virus = v,
      feature = names(beta),
      coef = as.numeric(beta)
    )
  }

  # OOB performance
  pred_oob <- predict(fit, newx = X_oob, s = lambda_ref)
  if (length(dim(pred_oob)) == 3) {
    pred_oob <- pred_oob[, , 1]
  }
  if (OUTCOME == "composite") pred_oob <- matrix(pred_oob, ncol = 1, dimnames = list(NULL, "comp"))

  for (v in VIRUSES_BOOT) {
    mse_oob <- mean((pred_oob[, v] - Y_oob[, v])^2)
    perf_rows[[length(perf_rows) + 1]] <- data.table(
      bootstrap = b,
      virus = v,
      en_oob_mse = mse_oob,
      en_lambda = lambda_ref
    )
  }
}

if (N_CORES > 1) {
  try(parallel::stopCluster(cl), silent = TRUE)
}

coef_dt <- rbindlist(coef_rows, fill = TRUE)
perf_dt <- rbindlist(perf_rows, fill = TRUE)

# Summary statistics
cat("\n--- Coefficient stability summary ---\n")
stability <- coef_dt[feature != "(Intercept)", .(
  mean_coef = mean(coef),
  sd_coef = sd(coef),
  selection_freq = mean(coef != 0),
  sign_consistency = abs(mean(sign(coef[coef != 0])))
), by = .(virus, feature)][order(virus, -selection_freq)]

cat("Top 5 features per virus by selection frequency:\n")
for (v in VIRUSES_BOOT) {
  cat(sprintf("  %s:\n", v))
  top5 <- stability[virus == v][order(-selection_freq)][1:min(5, .N)]
  for (i in seq_len(nrow(top5))) {
    cat(sprintf("    %s: freq=%.2f, mean=%.4f, sign=%.2f\n",
                top5$feature[i], top5$selection_freq[i],
                top5$mean_coef[i], top5$sign_consistency[i]))
  }
}

timings_boot <- list(
  total_s = as.numeric(difftime(Sys.time(), tic_total, units = "secs"))
)

cat(sprintf("\nTotal time: %.1f seconds (%.1f minutes)\n",
            timings_boot$total_s, timings_boot$total_s / 60))

save(
  param_out,
  timings_boot,
  atlas_file,
  lambda_ref,
  VIRUSES_BOOT,
  ALPHA,
  coef_dt,
  perf_dt,
  stability,
  file = out_file
)

cat("Saved: ", out_file, "\n", sep = "")
cat("Done.\n")
