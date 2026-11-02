# =============================================================================
# Regenerate Figure 4: M4 standardized feature importance
# Uses PRODUCTION ATLAS (33 clusters with consistent IDs matching annotations)
# Phase 1 fixes from 2026-04-28 meeting:
#   1.1 Cluster ordering: numeric (1..33), not lexicographic
#   1.2 Marker shorthand on cluster labels (from cluster_annotations_corrected.csv)
#   1.3 Sign of EN coefficients: signed importance, red/blue for neg/pos
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BIOLOGY_FILE  <- Sys.getenv("FLU_PROD_ATLAS", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/final_atlas_timing_corrected.RData")
BOOT_FILE     <- Sys.getenv("FLU_BOOT_FILE", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/bootstrap_importance_3f0f3ae8.RData")
ANN_FILE      <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
OUT_PNG       <- "/home/amir/repo/flu-cellular-prediction/paper/images/E2_all_feature_importance.png"
OUT_PDF       <- "/home/amir/repo/flu-cellular-prediction/paper/images/E2_all_feature_importance.pdf"
VIRUSES       <- c("h1n1", "h3n2", "bvic")

# Load production-atlas data ---------------------------------------------------
e_bio <- new.env(parent = emptyenv())
load(BIOLOGY_FILE, envir = e_bio)

e_boot <- new.env(parent = emptyenv())
load(BOOT_FILE, envir = e_boot)

clr_features <- e_bio$clr_features  # 343 patients x 33 clusters
coef_dt      <- as.data.table(e_boot$coef_dt)  # bootstrap, virus, feature, coef

# Standardize every coefficient using its actual production feature SD.
feature_sd <- apply(e_bio$X_clust, 2, sd)
sd_dt <- data.table(feature = names(feature_sd), sd = as.numeric(feature_sd))

# Annotations -------------------------------------------------------------------
ann <- fread(ANN_FILE)
ann[, cluster_id := as.integer(cluster_id)]

short_label <- function(annotation) {
  out <- gsub("\\s*\\(.*?\\)", "", annotation)
  out <- gsub("\\s+", " ", out)
  trimws(out)
}
ann[, short := short_label(annotation)]

# Average signed coef across bootstraps ----------------------------------------
agg <- coef_dt[feature != "(Intercept)",
               .(mean_coef = mean(coef), sd_coef = sd(coef)),
               by = .(virus, feature)]

# Multiply by feature SD to get standardized signed importance -----------------
agg <- merge(agg, sd_dt, by = "feature", all.x = TRUE)
# For covariates without a known SD (not in clr_features), use sd = 1
stopifnot(!anyNA(agg$sd))
agg[, signed_imp := mean_coef * sd]

# Build display label and ordering ---------------------------------------------
agg[, is_cluster := grepl("^cluster_", feature)]
agg[, cluster_id := ifelse(is_cluster,
                            as.integer(sub("cluster_", "", feature)),
                            NA_integer_)]

covariate_pretty <- c(
  "baseline_h1n1"        = "Baseline H1N1 titer",
  "baseline_h3n2"        = "Baseline H3N2 titer",
  "baseline_bvic"        = "Baseline B/Vic titer",
  "dose_HD"              = "Dose: HD",
  "age"                  = "Age",
  "sex_Male"             = "Sex: Male",
  "time_post_transplant" = "Time post-transplant"
)

agg[, display_label := ifelse(
  is_cluster,
  paste0("C", cluster_id, ": ",
         ann[match(agg$cluster_id, cluster_id), short]),
  ifelse(feature %in% names(covariate_pretty),
         covariate_pretty[feature], feature)
)]

# Order: covariates first (in fixed natural order), then clusters numerically
covariate_order <- c("baseline_h1n1", "baseline_h3n2", "baseline_bvic",
                      "dose_HD", "age", "sex_Male", "time_post_transplant")

agg[, sort_key := ifelse(is_cluster,
                          1000L + cluster_id,
                          match(feature, covariate_order))]
# Drop any unmapped feature (shouldn't happen but safety)
agg <- agg[!is.na(sort_key)]

# Factor levels: coord_flip reverses, so reverse order in factor
display_levels <- unique(agg[order(-sort_key), display_label])
agg[, display_label := factor(display_label, levels = display_levels)]

# Plot --------------------------------------------------------------------------
agg[, virus := factor(virus, levels = c("h1n1", "h3n2", "bvic"))]
p <- ggplot(agg, aes(x = display_label, y = signed_imp, fill = signed_imp > 0)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ virus, ncol = 3, scales = "free_x",
             labeller = labeller(virus = c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria"))) +
  scale_fill_manual(values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
                    labels = c("TRUE" = "Positive predictor", "FALSE" = "Negative predictor"),
                    name = NULL) +
  labs(
    x = NULL,
    y = "Average signed importance across 200 bootstraps"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(OUT_PNG, p, width = 7.2, height = 6.6, dpi = 600)
ggsave(OUT_PDF, p, width = 7.2, height = 6.6)
cat("Saved:", OUT_PNG, "\n")
cat("Saved:", OUT_PDF, "\n")
