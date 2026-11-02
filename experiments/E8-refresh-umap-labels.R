#!/usr/bin/env Rscript
# The faceted UMAP used legacy cell-type names. Use atlas IDs consistently;
# the current cell-type interpretation is in S4 and Supplementary Table S1.
# Cached coordinates and assignments are retained; no UMAP or atlas is refit.
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
e <- new.env(parent = emptyenv())
load(FLU_BIO_FILE, envir = e)
d <- cbind(as.data.table(e$dt_keep)[, .(meta = as.integer(meta))],
           as.data.table(e$umap_dt))
K <- max(d$meta)
if (!nzchar(Sys.getenv("FLU_BIO_FILE"))) stopifnot(nrow(d) == 173353, K == 33)
stopifnot(identical(sort(unique(d$meta)), seq_len(K)))
d[, cluster_label := factor(paste0("C", meta), levels = paste0("C", seq_len(K)))]
background <- d[, .(UMAP1, UMAP2)]
p <- ggplot(d, aes(UMAP1, UMAP2)) +
  geom_point(data = background, color = "gray85", size = 0.05, alpha = 0.2) +
  geom_point(color = "#d62728", size = 0.15, alpha = 0.7) +
  facet_wrap(~ cluster_label, ncol = if (K > 36) 7 else 6) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold", size = 13),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave("paper/images/fig3_umap_faceted_annotated.png", p,
       width = 16, height = if (K > 36) 15 else 14, dpi = 200, bg = "white")
