#!/usr/bin/env Rscript
# ============================================================================
# E3: Recipient-clustered paired bootstrap for prediction-error comparisons
# ----------------------------------------------------------------------------
# Addresses Andrew's review concern: the original bootstrap in
# experiments/E2-posthoc.R::compute_bootstrap_delta resamples observation
# rows. Some recipients contribute multiple vaccination episodes (yr1/yr2),
# whose residuals are correlated. Here we resample the BASE RECIPIENT
# (the part of patient_id before "_yr<n>"), keeping each recipient's bundle
# of 1-2 episodes x 3 antigens together. We also report a TWO-SIDED p-value.
# ============================================================================

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf(fmt, ...)))
  flush.console()
}

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
OUTPUT_ROOT <- FLU_ROOT
N_MICRO     <- Sys.getenv("N_MICRO", "5000")
FOLD_DIR    <- file.path(OUTPUT_ROOT, paste0("per_fold_atlas_", N_MICRO))
OUT_DIR     <- file.path(OUTPUT_ROOT, paste0("post_hoc_", N_MICRO), "recipient_bootstrap")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_RDATA <- file.path(OUT_DIR, "E3_recipient_boot.RData")
OUT_CSV   <- file.path(OUT_DIR, "E3_recipient_boot.csv")

VIRUSES <- c("h1n1", "h3n2", "bvic")
N_FOLDS <- 5L
B       <- 2000L
SEED    <- 7000L
N_CORES <- as.integer(Sys.getenv("N_CORES", "12"))

# ----------------------------------------------------------------------------
# Load per-fold predictions and build long table
# ----------------------------------------------------------------------------
log_msg("Loading %d fold artifacts...", N_FOLDS)

long_rows <- vector("list", N_FOLDS)
for (i in seq_len(N_FOLDS)) {
  fp <- file.path(FOLD_DIR, sprintf("fold_%d", i),
                  sprintf("fold_%d_all.RData", i))
  if (!file.exists(fp)) stop("Missing fold artifact: ", fp)

  e <- new.env(parent = emptyenv())
  load(fp, envir = e)

  pids <- e$test_p_ids
  Yt   <- e$Y_test
  p1   <- e$pred_M1_en_test
  p2   <- e$pred_M2_en_test
  p3   <- e$pred_M3_en_test
  p4   <- e$pred_M4_en_test

  if (!all(colnames(Yt) == VIRUSES)) {
    stop("Unexpected column order in Y_test for fold ", i)
  }
  if (nrow(Yt) != length(pids)) {
    stop("Row count mismatch for fold ", i)
  }

  # Stack across antigens
  fold_dt <- rbindlist(lapply(seq_along(VIRUSES), function(v) {
    data.table(
      fold       = i,
      patient_id = pids,
      virus      = VIRUSES[v],
      Y_true     = Yt[, v],
      pred_M1    = p1[, v],
      pred_M2    = p2[, v],
      pred_M3    = p3[, v],
      pred_M4    = p4[, v]
    )
  }))
  long_rows[[i]] <- fold_dt
}

dt <- rbindlist(long_rows)
dt[, base_recipient := sub("_yr[0-9]+$", "", patient_id)]
dt[, virus := factor(virus, levels = VIRUSES)]

n_obs_total       <- nrow(dt)                              # 343 * 3 = 1029
n_episodes        <- length(unique(dt$patient_id))         # 343
n_base_recipients <- length(unique(dt$base_recipient))     # ~263

log_msg("Long table: %d rows | %d episodes | %d base recipients",
        n_obs_total, n_episodes, n_base_recipients)

stopifnot(n_obs_total == n_episodes * 3L)

# Squared errors for each model
for (k in 1:4) {
  set(dt, j = paste0("err_M", k),
      value = (dt[[paste0("pred_M", k)]] - dt$Y_true)^2)
}

# ----------------------------------------------------------------------------
# Pre-compute structures for fast recipient resampling
# ----------------------------------------------------------------------------
# Map: base_recipient -> integer row indices into `dt`
setkey(dt, base_recipient)
recipient_ids <- unique(dt$base_recipient)
n_p <- length(recipient_ids)

recipient_to_rows <- split(seq_len(nrow(dt)), dt$base_recipient)
# Ensure order matches recipient_ids
recipient_to_rows <- recipient_to_rows[recipient_ids]

# Per-virus row index lookup for fast subsetting per scope
virus_row_idx <- list(
  overall = seq_len(nrow(dt)),
  h1n1    = which(dt$virus == "h1n1"),
  h3n2    = which(dt$virus == "h3n2"),
  bvic    = which(dt$virus == "bvic")
)

# ----------------------------------------------------------------------------
# Single (comparison, scope) job
# ----------------------------------------------------------------------------
COMPARISONS <- data.frame(
  comparison = c("M1_to_M2", "M1_to_M3", "M1_to_M4",
                 "M2_to_M3", "M2_to_M4", "M3_to_M4"),
  model_a    = c("M1", "M1", "M1", "M2", "M2", "M3"),
  model_b    = c("M2", "M3", "M4", "M3", "M4", "M4"),
  stringsAsFactors = FALSE
)

SCOPES <- c("overall", VIRUSES)

jobs <- expand.grid(
  cmp_idx   = seq_len(nrow(COMPARISONS)),
  scope_idx = seq_along(SCOPES),
  stringsAsFactors = FALSE
)

run_one <- function(job_row) {
  cmp_i   <- job_row$cmp_idx
  scope_i <- job_row$scope_idx
  cmp     <- COMPARISONS[cmp_i, ]
  scope   <- SCOPES[scope_i]

  err_a_all <- dt[[paste0("err_", cmp$model_a)]]
  err_b_all <- dt[[paste0("err_", cmp$model_b)]]

  # Scope filter: precomputed virus row indices, intersected per recipient
  scope_rows <- virus_row_idx[[scope]]
  in_scope   <- logical(nrow(dt))
  in_scope[scope_rows] <- TRUE

  # For each recipient, store the scope-filtered row indices
  rec_scope_rows <- lapply(recipient_to_rows, function(ix) ix[in_scope[ix]])

  # Observed (full data) MSEs
  obs_a <- mean(err_a_all[scope_rows])
  obs_b <- mean(err_b_all[scope_rows])
  observed_delta_pct <- (obs_a - obs_b) / obs_a * 100

  # Bootstrap: distinct seed per (comparison, scope) for reproducibility
  set.seed(SEED + cmp_i * 1000L + scope_i)
  boot_delta <- numeric(B)
  for (b in seq_len(B)) {
    samp_rec <- sample.int(n_p, size = n_p, replace = TRUE)
    rows_b   <- unlist(rec_scope_rows[samp_rec], use.names = FALSE)
    if (length(rows_b) == 0L) {
      boot_delta[b] <- NA_real_
      next
    }
    mse_a_b <- mean(err_a_all[rows_b])
    mse_b_b <- mean(err_b_all[rows_b])
    boot_delta[b] <- (mse_a_b - mse_b_b) / mse_a_b * 100
  }

  ci <- as.numeric(stats::quantile(boot_delta,
                                   probs = c(0.025, 0.975),
                                   names = FALSE, na.rm = TRUE))
  p_two_sided <- 2 * min(mean(boot_delta <= 0, na.rm = TRUE),
                         mean(boot_delta >= 0, na.rm = TRUE))
  p_two_sided <- min(p_two_sided, 1)

  list(
    comparison  = cmp$comparison,
    model_a     = cmp$model_a,
    model_b     = cmp$model_b,
    scope       = scope,
    delta_pct   = observed_delta_pct,
    ci_low_pct  = ci[1],
    ci_high_pct = ci[2],
    p_two_sided = p_two_sided,
    boot_draws  = boot_delta
  )
}

# ----------------------------------------------------------------------------
# Run jobs in parallel
# ----------------------------------------------------------------------------
log_msg("Running %d jobs (B=%d) on %d cores...",
        nrow(jobs), B, N_CORES)

job_list <- split(jobs, seq_len(nrow(jobs)))

t0 <- Sys.time()
results <- mclapply(
  job_list,
  function(jr) run_one(jr),
  mc.cores    = N_CORES,
  mc.set.seed = FALSE  # we set seed deterministically inside run_one
)
t1 <- Sys.time()
log_msg("Bootstrap finished in %.1f sec", as.numeric(difftime(t1, t0, units = "secs")))

# Assemble table
results_df <- rbindlist(lapply(results, function(r) {
  data.table(
    comparison    = r$comparison,
    model_a       = r$model_a,
    model_b       = r$model_b,
    scope         = r$scope,
    delta_pct     = r$delta_pct,
    ci_low_pct    = r$ci_low_pct,
    ci_high_pct   = r$ci_high_pct,
    p_two_sided   = r$p_two_sided,
    n_patients    = n_base_recipients,
    n_observations = n_obs_total
  )
}))

# Order: comparison then scope (overall, h1n1, h3n2, bvic)
scope_order <- c("overall", "h1n1", "h3n2", "bvic")
cmp_order   <- COMPARISONS$comparison
results_df[, comparison := factor(comparison, levels = cmp_order)]
results_df[, scope      := factor(scope, levels = scope_order)]
setorder(results_df, comparison, scope)
results_df[, comparison := as.character(comparison)]
results_df[, scope      := as.character(scope)]

# Bootstrap draws keyed by "comparison|scope"
boot_draws <- setNames(
  lapply(results, function(r) r$boot_draws),
  vapply(results, function(r) paste(r$comparison, r$scope, sep = "|"),
         character(1))
)

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
recipient_boot_meta <- list(
  B                = B,
  seed             = SEED,
  n_cores          = N_CORES,
  n_base_recipients = n_base_recipients,
  n_episodes       = n_episodes,
  n_observations   = n_obs_total,
  comparisons      = COMPARISONS,
  scopes           = SCOPES,
  long_table_cols  = c("fold", "patient_id", "base_recipient", "virus",
                       "Y_true", paste0("pred_M", 1:4), paste0("err_M", 1:4))
)

save(
  results_df,
  boot_draws,
  recipient_boot_meta,
  dt,
  file = OUT_RDATA
)
fwrite(results_df, OUT_CSV)

log_msg("Saved RData: %s", OUT_RDATA)
log_msg("Saved CSV:   %s", OUT_CSV)

# ----------------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------------
cat("\n================ E3 Recipient-Clustered Paired Bootstrap ================\n")
cat(sprintf("n_base_recipients = %d | n_episodes = %d | n_observations = %d | B = %d\n\n",
            n_base_recipients, n_episodes, n_obs_total, B))

print_df <- as.data.frame(results_df)
print_df$delta_pct   <- sprintf("%+6.2f", print_df$delta_pct)
print_df$ci_low_pct  <- sprintf("%+6.2f", print_df$ci_low_pct)
print_df$ci_high_pct <- sprintf("%+6.2f", print_df$ci_high_pct)
print_df$p_two_sided <- sprintf("%.4f", print_df$p_two_sided)
print(print_df[, c("comparison", "scope", "delta_pct",
                   "ci_low_pct", "ci_high_pct", "p_two_sided")],
      row.names = FALSE)

log_msg("Done.")
