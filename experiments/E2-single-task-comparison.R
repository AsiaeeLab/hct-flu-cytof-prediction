#!/usr/bin/env Rscript
# E2-single-task-comparison.R
# Compare multitask (mgaussian) vs single-task (gaussian) elastic net
# Uses saved per-fold X/Y matrices — no raw data needed.

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})

source(here::here("methods", "35-per-fold-atlas-cv.R"))

cat("=== Single-task vs Multitask Comparison ===\n")

PER_FOLD_ROOT <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "per_fold_atlas_5000")
OUT_DIR <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "post_hoc_5000/figures"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SEED_BASE <- 7100
ALPHA <- 0.9
VIRUSES <- c("h1n1", "h3n2", "bvic")
n_folds <- 5

# --- Load aggregated data (for multitask results) ---
load(file.path(PER_FOLD_ROOT, "aggregated", "E2_aggregated.RData"))
oof_multi_M3 <- oof_pred_M3_en
oof_multi_M4 <- oof_pred_M4_en
oof_Y_all <- oof_Y

# --- Run single-task models per fold ---
cat("\nFitting single-task (per-virus) elastic net...\n")
oof_single_M3 <- matrix(NA, nrow = nrow(oof_Y_all), ncol = 3, dimnames = dimnames(oof_Y_all))
oof_single_M4 <- matrix(NA, nrow = nrow(oof_Y_all), ncol = 3, dimnames = dimnames(oof_Y_all))

for (f in 1:n_folds) {
  cat(sprintf("  Fold %d...\n", f))
  env <- new.env()
  load(file.path(PER_FOLD_ROOT, paste0("fold_", f), paste0("fold_", f, "_all.RData")), envir = env)

  # Get fold test indices in the aggregated oof matrix
  test_ids <- env$test_p_ids
  test_idx <- match(test_ids, common_p_ids)

  for (v in VIRUSES) {
    # --- M3: titers + clusters (single-task) ---
    X_tr <- env$X_M3_train
    X_te <- env$X_M3_test
    y_tr <- env$Y_train[, v]
    y_te <- env$Y_test[, v]

    set.seed(SEED_BASE + f * 100 + which(VIRUSES == v))
    cvfit <- cv.glmnet(x = X_tr, y = y_tr, family = "gaussian",
                       alpha = ALPHA, standardize = TRUE, foldid = compute_grouped_foldid(env$train_p_ids, 5L, SEED_BASE + f * 100))
    pred <- as.numeric(predict(cvfit, newx = X_te, s = "lambda.min"))
    oof_single_M3[test_idx, v] <- pred

    # --- M4: full (single-task) ---
    X_tr4 <- env$X_M4_train
    X_te4 <- env$X_M4_test

    set.seed(SEED_BASE + f * 100 + which(VIRUSES == v) + 10)
    cvfit4 <- cv.glmnet(x = X_tr4, y = y_tr, family = "gaussian",
                        alpha = ALPHA, standardize = TRUE, foldid = compute_grouped_foldid(env$train_p_ids, 5L, SEED_BASE + f * 100))
    pred4 <- as.numeric(predict(cvfit4, newx = X_te4, s = "lambda.min"))
    oof_single_M4[test_idx, v] <- pred4
  }
}

# --- Compute MSE ---
compute_mse <- function(Y_true, Y_pred, label) {
  per_virus <- sapply(VIRUSES, function(v) mean((Y_true[, v] - Y_pred[, v])^2))
  overall <- mean((Y_true - Y_pred)^2)
  data.frame(model = label, virus = c(VIRUSES, "overall"),
             mse = c(per_virus, overall), stringsAsFactors = FALSE)
}

# Load M1 multitask for reference
mse_df <- rbind(
  compute_mse(oof_Y_all, oof_pred_M1_en, "M1_multi"),
  compute_mse(oof_Y_all, oof_multi_M3, "M3_multi"),
  compute_mse(oof_Y_all, oof_multi_M4, "M4_multi"),
  compute_mse(oof_Y_all, oof_single_M3, "M3_single"),
  compute_mse(oof_Y_all, oof_single_M4, "M4_single")
)

cat("\n=== MSE Summary ===\n")
print(reshape(mse_df, direction = "wide", idvar = "virus", timevar = "model"))

# --- Compute improvements ---
m1_overall <- mse_df$mse[mse_df$model == "M1_multi" & mse_df$virus == "overall"]
m3_multi_o <- mse_df$mse[mse_df$model == "M3_multi" & mse_df$virus == "overall"]
m3_single_o <- mse_df$mse[mse_df$model == "M3_single" & mse_df$virus == "overall"]
m4_multi_o <- mse_df$mse[mse_df$model == "M4_multi" & mse_df$virus == "overall"]
m4_single_o <- mse_df$mse[mse_df$model == "M4_single" & mse_df$virus == "overall"]

cat(sprintf("\nM3 improvement over M1:\n"))
cat(sprintf("  Multitask:  %.2f%%\n", 100 * (1 - m3_multi_o / m1_overall)))
cat(sprintf("  Single-task: %.2f%%\n", 100 * (1 - m3_single_o / m1_overall)))
cat(sprintf("\nM4 improvement over M1:\n"))
cat(sprintf("  Multitask:  %.2f%%\n", 100 * (1 - m4_multi_o / m1_overall)))
cat(sprintf("  Single-task: %.2f%%\n", 100 * (1 - m4_single_o / m1_overall)))

# --- Figure ---
plot_df <- mse_df[mse_df$virus == "overall" & mse_df$model != "M1_multi", ]
plot_df$task <- ifelse(grepl("multi", plot_df$model), "Multitask", "Single-task")
plot_df$model_base <- gsub("_.*", "", plot_df$model)
plot_df$improvement <- 100 * (1 - plot_df$mse / m1_overall)

p1 <- ggplot(plot_df, aes(x = model_base, y = improvement, fill = task)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", improvement)),
            position = position_dodge(width = 0.6), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Multitask" = "#2166ac", "Single-task" = "#92c5de"),
                    name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "MSE improvement (%) over M1") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Per-virus comparison for M3
pv_df <- mse_df[mse_df$model %in% c("M3_multi", "M3_single") & mse_df$virus != "overall", ]
pv_df$task <- ifelse(grepl("multi", pv_df$model), "Multitask", "Single-task")
pv_m1 <- mse_df[mse_df$model == "M1_multi" & mse_df$virus != "overall", ]
pv_df <- merge(pv_df, pv_m1[, c("virus", "mse")], by = "virus", suffixes = c("", "_m1"))
pv_df$improvement <- 100 * (1 - pv_df$mse / pv_df$mse_m1)
pv_df$virus_label <- factor(pv_df$virus, levels = c("h1n1", "h3n2", "bvic"),
                             labels = c("H1N1", "H3N2", "B/Vic"))

p2 <- ggplot(pv_df, aes(x = virus_label, y = improvement, fill = task)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", improvement)),
            position = position_dodge(width = 0.6), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Multitask" = "#2166ac", "Single-task" = "#92c5de"),
                    name = "") +
  labs(x = NULL, y = "MSE improvement (%) over M1") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

fig <- (p1 | p2) + plot_annotation(tag_levels = "A", tag_suffix = ".")
fig <- fig & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(OUT_DIR, "E2_multitask_vs_single.pdf"), fig, width = 12, height = 5)
ggsave(file.path(OUT_DIR, "E2_multitask_vs_single.png"), fig, width = 12, height = 5, dpi = 300)
file.copy(file.path(OUT_DIR, "E2_multitask_vs_single.png"),
          "/home/amir/repo/flu-cellular-prediction/paper/images/E2_multitask_vs_single.png",
          overwrite = TRUE)

# Save CSV
write.csv(mse_df, file.path(OUT_DIR, "E2_multitask_vs_single.csv"), row.names = FALSE)

# Save RData
save(oof_single_M3, oof_single_M4, mse_df,
     file = file.path(PER_FOLD_ROOT, "aggregated", "E2_single_task.RData"))

cat("\nFigure saved to:", file.path(OUT_DIR, "E2_multitask_vs_single.png"), "\n")
cat("Done.\n")
