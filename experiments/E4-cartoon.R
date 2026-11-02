# Conceptual cartoon: early reconstitution -> mature response capacity
# Using geom_rect for explicit per-stage bar widths (fixes the empty-bar bug).

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

FIG_DIR <- "/home/amir/repo/flu-cellular-prediction/paper/images"

cat_colors <- c(
  "Mature B / memory B"   = "#08519c",
  "Activated follicular B"= "#3182bd",
  "Plasmablast"           = "#6baed6",
  "CD4+ memory T"         = "#fd8d3c",
  "CD8+ memory T"         = "#fb6a4a",
  "Activated APC"         = "#fdae6b",
  "Immature / DN"         = "#737373",
  "Lineage-negative"      = "#bdbdbd"
)
cat_levels <- c("Immature / DN", "Lineage-negative",
                "Activated APC", "Activated follicular B",
                "Mature B / memory B", "Plasmablast",
                "CD4+ memory T", "CD8+ memory T")

# Three stages with (center position, half-width, total composition vector)
stages <- list(
  list(label = "Early reconstitution\n(non-responder phenotype)",
       x = 1, half_w = 0.20,
       comp = c("Immature / DN" = 28, "Lineage-negative" = 22,
                "Activated APC" = 14, "Activated follicular B" = 6,
                "Mature B / memory B" = 8, "Plasmablast" = 2,
                "CD4+ memory T" = 12, "CD8+ memory T" = 8)),
  list(label = "Intermediate",
       x = 2, half_w = 0.35,
       comp = c("Immature / DN" = 18, "Lineage-negative" = 14,
                "Activated APC" = 14, "Activated follicular B" = 12,
                "Mature B / memory B" = 14, "Plasmablast" = 4,
                "CD4+ memory T" = 14, "CD8+ memory T" = 10)),
  list(label = "Mature response capacity\n(responder phenotype)",
       x = 3, half_w = 0.50,
       comp = c("Immature / DN" = 8, "Lineage-negative" = 6,
                "Activated APC" = 12, "Activated follicular B" = 16,
                "Mature B / memory B" = 22, "Plasmablast" = 6,
                "CD4+ memory T" = 18, "CD8+ memory T" = 12))
)

# Build a long-format rectangle data frame
rect_dt <- rbindlist(lapply(stages, function(s) {
  ord <- match(cat_levels, names(s$comp))
  pct <- s$comp[ord]
  ymin <- cumsum(c(0, pct[-length(pct)]))
  ymax <- cumsum(pct)
  data.table(
    stage = s$label, x = s$x,
    xmin = s$x - s$half_w, xmax = s$x + s$half_w,
    ymin = ymin, ymax = ymax,
    category = cat_levels
  )
}))
rect_dt[, category := factor(category, levels = cat_levels)]

# Labels for stages
stage_labels <- data.table(
  x = c(1, 2, 3),
  label = sapply(stages, function(s) s$label)
)

# Arrows between stages
arrows_dt <- data.table(
  x1 = c(1.25, 2.4), x2 = c(1.65, 2.5),
  y = c(108, 108),
  label = c("year 1 → year 2 in many\nNR → R repeaters",
            "expansion of mature B-cell\nand T-memory populations")
)

p <- ggplot() +
  geom_rect(data = rect_dt,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax, fill = category),
            color = "white", linewidth = 0.3) +
  geom_segment(data = arrows_dt,
               aes(x = x1, xend = x2, y = y, yend = y),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "gray30", linewidth = 0.5) +
  geom_text(data = arrows_dt,
            aes(x = (x1 + x2)/2, y = y + 6, label = label),
            color = "gray25", size = 2.7, lineheight = 0.95) +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = stage_labels$label,
                     limits = c(0.4, 3.6)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)),
                     limits = c(0, 124)) +
  scale_fill_manual(values = cat_colors, name = NULL,
                    breaks = cat_levels) +
  labs(x = NULL, y = "Composition (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 9, lineheight = 1.0, vjust = 1),
        plot.margin = margin(15, 6, 6, 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(file.path(FIG_DIR, "E4_cartoon_reconstitution.png"), p,
       width = 10, height = 5.5, dpi = 300)
ggsave(file.path(FIG_DIR, "E4_cartoon_reconstitution.pdf"), p,
       width = 10, height = 5.5, device = cairo_pdf)
cat("Cartoon saved with explicit rectangles.\n")
