# Marker-profile heatmap using CATALYST-style expression scaling.
#
# Replaces the earlier z-scored / arithmetic-mean heatmap. Each marker channel is
# independently rescaled to the 0-1 range after clipping at its 1st and 99th
# percentiles across all cells, and each cluster is then summarized by the MEDIAN
# of the scaled values. This is the convention used by CATALYST's plotExprHeatmap
# (via .scale_exprs) and by FCSimple::fcs_cluster_heatmap.
#
# Display updated 2026-09-05 to match Josh Simmons's FCSimple palette:
# https://github.com/jsim91/FCSimple/blob/main/R/fcs_cluster_heatmap.R
# The reversed 11-color RdYlBu palette runs blue -> pale yellow -> red.
# Columns follow tmp_josh/heatmap_original-1.png for visual comparison;
# cluster IDs remain those of this paper's 33-cluster atlas. Color is a
# relative display scale, not a positivity gate or evidence of exact matching.

FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

BIO_FILE <- FLU_BIO_FILE
ANN_FILE <- FLU_ANN_FILE
FIG_DIR  <- "/home/amir/repo/flu-cellular-prediction/paper/images"

SCALING_QUANTILE <- 0.01   # clip at the 1st and 99th percentile of each channel

e <- new.env(parent = emptyenv())
load(BIO_FILE, envir = e)
dt_keep      <- as.data.table(e$dt_keep)
marker_names <- e$marker_names

annotation_dt <- fread(ANN_FILE)
setorder(annotation_dt, cluster_id)

# ---------------------------------------------------------------------------
# Per-channel scaling: clip at the quantile bounds, then rescale to 0-1
# ---------------------------------------------------------------------------
scale_channel <- function(x, q = SCALING_QUANTILE) {
  bounds <- quantile(x, probs = c(q, 1 - q), na.rm = TRUE)
  lo <- bounds[1]; hi <- bounds[2]
  if (!is.finite(hi - lo) || (hi - lo) == 0) return(rep(0, length(x)))
  (pmin(pmax(x, lo), hi) - lo) / (hi - lo)
}

scaled <- dt_keep[, lapply(.SD, scale_channel), .SDcols = marker_names]
scaled[, meta := dt_keep$meta]

# Cluster summary is the MEDIAN of the scaled values, not the mean
med <- scaled[, lapply(.SD, median, na.rm = TRUE), by = meta, .SDcols = marker_names]
setorder(med, meta)

mat <- as.matrix(med[, ..marker_names])
rownames(mat) <- paste0("C", med$meta)

# ---------------------------------------------------------------------------
# Cluster rows; fix columns to the collaborator's original marker order.
# ---------------------------------------------------------------------------
row_hc <- hclust(dist(mat), method = "ward.D2")
reference_marker_order <- c(
  "CD4", "CD45RO", "CD3", "CD27", "CD8", "CD38", "CD127", "CD28",
  "HLA-DR", "CD19", "CD20", "CXCR5", "CD21", "CCR2", "CD86", "CD138",
  "CCR7", "CD137", "OX40", "CCR6", "BTLA", "CXCR3", "PD-1", "ICOS",
  "CCR4", "Ki-67", "TIGIT", "CD69", "CD57", "CD25"
)
stopifnot(setequal(reference_marker_order, colnames(mat)))

lab <- merge(med[, .(meta)], annotation_dt[, .(cluster_id, annotation)],
             by.x = "meta", by.y = "cluster_id", all.x = TRUE)
setorder(lab, meta)
lab[is.na(annotation), annotation := "Unannotated"]
row_labels <- paste0("C", lab$meta, " ", lab$annotation)

long <- melt(as.data.table(mat, keep.rownames = "row"),
             id.vars = "row", variable.name = "marker", value.name = "scaled")
long[, idx := as.integer(sub("^C", "", row))]
long <- merge(long, data.table(idx = lab$meta, label = row_labels), by = "idx")

long[, label  := factor(label,  levels = rev(row_labels[row_hc$order]))]
long[, marker := factor(marker, levels = reference_marker_order)]

p <- ggplot(long, aes(x = marker, y = label, fill = scaled)) +
  geom_tile(color = "grey30", linewidth = 0.12) +
  scale_fill_gradientn(
    colours = grDevices::colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    limits = c(0, 1), oob = scales::squish,
    breaks = seq(0, 1, 0.2), name = "Median scaled\nexpression"
  ) +
  labs(x = "Marker", y = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
    axis.text.y = element_text(size = 5.8, color = "black"),
    panel.grid  = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(FIG_DIR, "fig3b_marker_heatmap.png"), p,
       width = 7.2, height = 6.6, dpi = 600, bg = "white")
fwrite(med, file.path(FIG_DIR, "figS4_scaled_marker_medians.csv"))

cat("Wrote", file.path(FIG_DIR, "fig3b_marker_heatmap.png"), "\n")
cat("Cells used:", nrow(dt_keep), "| clusters:", nrow(med), "\n\n")

# Markers most affected by the change in scaling, reported so the annotation
# review can focus on the clusters whose profiles moved.
cat("Scaled expression of selected markers on the follicular-helper candidates:\n")
sel <- c("CXCR5", "PD-1", "TIGIT", "CD4", "CD45RO", "Ki-67", "CD19", "CD38", "CD3")
sel <- intersect(sel, marker_names)
print(med[meta %in% c(6, 10, 16, 25, 28, 29, 32), c("meta", sel), with = FALSE])
