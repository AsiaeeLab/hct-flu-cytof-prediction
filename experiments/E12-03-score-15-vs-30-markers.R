#!/usr/bin/env Rscript
# Cellular score computed from 15 versus 30 markers.
#
# Uses the reduced-marker cell reassignments of E2-marker-reduction-perfold.R:
# for each fold, a fixed sample of up to 2,000 cells per episode was reassigned
# to that fold's atlas using either all 30 markers or the fold's 15 top-ranked
# markers (ranked by between-population variance of the atlas centroids, an
# outcome-free criterion). Both versions use the same cells, so differences
# between them come from the marker reduction alone. The full-data score of the
# primary analysis (all cells, 30 markers) is shown as the reference.
#
# For each version, M3 (cellular features plus the three baseline titers) and
# its composite analogue are refitted in each training fold with the primary
# grouped inner folds, and held-out episodes are scored as in E6-01.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({library(data.table); library(glmnet); library(here);
  library(ggplot2); library(patchwork)})
source(here("methods", "35-per-fold-atlas-cv.R"))

ROOT <- FLU_ROOT
FOLD_ROOT <- file.path(ROOT, "per_fold_atlas_5000")
MR <- file.path(ROOT, "post_hoc_5000/marker_reduction")
OUT <- file.path(ROOT, "post_hoc_5000/score_15_markers"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FIG_DIR <- "/home/amir/repo/flu-cellular-prediction/paper/images"
V <- c("h1n1","h3n2","bvic"); B <- 2000L; RESP <- log(4) - 1e-12

comp <- readRDS(file.path(ROOT, "post_hoc_5000/composite/composite_episode_table.rds"))

score_version <- function(f, e, counts) {
  tr <- e$train_p_ids; te <- e$test_p_ids
  clr <- transform_counts(counts, "clr")
  cl <- colnames(e$clr_test)
  stopifnot(setequal(colnames(clr), cl))
  clr <- clr[, cl, drop=FALSE]
  Xtr <- cbind(clr[tr,,drop=FALSE], e$X_M1_train); Xte_c <- clr[te,,drop=FALSE]
  fid <- compute_grouped_foldid(tr, 5L, 7100 + f*1000 + 99)
  fit <- cv.glmnet(Xtr, e$Y_train, family="mgaussian", alpha=0.9, foldid=fid, standardize=TRUE)
  cf <- coef(fit, s="lambda.min")
  sa <- sapply(seq_along(V), function(i) as.vector(Xte_c %*% as.matrix(cf[[i]])[cl,1]))
  cfit <- cv.glmnet(Xtr, rowMeans(e$Y_train), family="gaussian", alpha=0.9, foldid=fid, standardize=TRUE)
  sc <- as.vector(Xte_c %*% as.matrix(coef(cfit, s="lambda.min"))[cl,1])
  data.table(episode_id=te, score_h1n1=sa[,1], score_h3n2=sa[,2], score_bvic=sa[,3], score_comp=sc)
}

rows <- list()
for (f in 1:5) {
  e <- new.env(); load(file.path(FOLD_ROOT, paste0("fold_",f), sprintf("fold_%d_all.RData", f)), e)
  for (b in c(30L, 15L)) {
    cnt <- readRDS(file.path(MR, sprintf("counts_fold%d_budget%d.rds", f, b)))
    rows[[length(rows)+1]] <- score_version(f, e, cnt)[, `:=`(fold=f, version=paste0("sample_", b))]
  }
}
sc <- rbindlist(rows)
full <- comp[, .(episode_id, fold, score_h1n1, score_h3n2, score_bvic, score_comp)][, version := "full_30"]
allv <- rbind(full, sc, use.names=TRUE)
outcome_cols <- comp[, .(episode_id, base_id, quartile, y_h1n1, y_h3n2, y_bvic, y_comp)]
allv <- merge(allv, outcome_cols, by="episode_id")
outs <- c(V, "comp")

# Agreement between the 15- and 30-marker scores on the same cells
w <- dcast(allv, episode_id + fold + quartile ~ version, value.var=paste0("score_", outs))
agree <- rbindlist(lapply(outs, function(o) {
  s30 <- w[[paste0("score_",o,"_sample_30")]]; s15 <- w[[paste0("score_",o,"_sample_15")]]
  sfull <- w[[paste0("score_",o,"_full_30")]]
  cen <- function(x) x - ave(x, w$fold)
  g <- function(x) ave(x, w$quartile, FUN=function(z) z > median(z))
  data.table(outcome=o,
    r_15_vs_30=cor(cen(s15), cen(s30)), r_30sample_vs_full=cor(cen(s30), cen(sfull)),
    r_15_vs_full=cor(cen(s15), cen(sfull)),
    same_group_15_vs_30=mean(g(s15)==g(s30)), same_group_15_vs_full=mean(g(s15)==g(sfull)),
    same_group_30sample_vs_full=mean(g(s30)==g(sfull)))
}))
cat("Score agreement (fold-centred Pearson r; share of episodes in the same score group):\n"); print(agree)
fwrite(agree, file.path(OUT, "score_agreement.csv"))

# Pooled Q1-Q3 fold-rise ratios for each version (recipient bootstrap, E6 seed)
base_levels <- sort(unique(comp$base_id))
ord <- comp$episode_id                           # E6 row order
rows_by_base <- split(seq_along(ord), factor(comp$base_id, levels=base_levels))
set.seed(20260728L)
boot_index <- lapply(seq_len(B), function(b)
  unlist(rows_by_base[sample.int(length(base_levels), length(base_levels), replace=TRUE)], use.names=FALSE))
ratio_rows <- list()
for (ver in c("full_30","sample_30","sample_15")) {
  d <- allv[version==ver][match(ord, episode_id)]
  for (o in outs) {
    s <- d[[paste0("score_",o)]]; y <- d[[paste0("y_",o)]]
    hi_flag <- ave(s, d$quartile, FUN=function(z) z > median(z)) == 1
    q13 <- d$quartile != "Q4_high"
    hi <- which(q13 & hi_flag); lo <- which(q13 & !hi_flag)
    inh <- logical(nrow(d)); inh[hi] <- TRUE; inl <- logical(nrow(d)); inl[lo] <- TRUE
    bt <- vapply(boot_index, function(ib) exp(mean(y[ib[inh[ib]]]) - mean(y[ib[inl[ib]]])), numeric(1))
    ratio_rows[[length(ratio_rows)+1]] <- data.table(version=ver, outcome=o,
      ratio=exp(mean(y[hi])-mean(y[lo])), lo=quantile(bt,.025,names=FALSE), hi=quantile(bt,.975,names=FALSE),
      resp_diff_pp=100*(mean(y[hi]>=RESP)-mean(y[lo]>=RESP)))
  }
}
ratios <- rbindlist(ratio_rows)
cat("\nPooled Q1-Q3 geometric-mean fold-rise ratios (higher vs lower score):\n"); print(ratios)
fwrite(ratios, file.path(OUT, "q13_ratios_by_marker_version.csv"))
saveRDS(allv, file.path(OUT, "scores_by_marker_version.rds"))

# ---------------------------------------------------------------------------
# Supplementary figure
# ---------------------------------------------------------------------------
th <- theme_minimal(base_size=8, base_family='sans') + theme(
  panel.grid.minor=element_blank(), axis.text=element_text(size=6.5, color='grey20'),
  axis.title=element_text(size=7.5), plot.title=element_text(face='bold', size=9),
  plot.tag=element_text(face='bold', size=11), legend.position='bottom', legend.title=element_blank())
pa_d <- w[, .(fold=factor(fold), s30=score_comp_sample_30 - ave(score_comp_sample_30, fold),
              s15=score_comp_sample_15 - ave(score_comp_sample_15, fold))]
r_lab <- sprintf('italic(r)*" = %.2f"', agree[outcome=="comp", r_15_vs_30])
pa <- ggplot(pa_d, aes(s30, s15, color=fold)) +
  geom_abline(slope=1, intercept=0, color='grey70', linewidth=.3, linetype='dashed') +
  geom_point(size=.9, alpha=.75) +
  annotate('text', x=-Inf, y=Inf, hjust=-.15, vjust=1.4, label=r_lab, parse=TRUE, size=2.6) +
  scale_color_brewer(palette='Dark2', name='Fold') +
  labs(title='Composite Score', x='Score from 30 markers (fold-centered)', y='Score from 15 markers (fold-centered)') +
  th + theme(legend.title=element_text(size=7))
vlab <- c(full_30='30 markers, all cells', sample_30='30 markers, cell sample', sample_15='15 markers, cell sample')
olab <- c(h1n1='H1N1', h3n2='H3N2', bvic='B/Victoria', comp='Composite')
pb_d <- copy(ratios)[, `:=`(version=factor(vlab[version], levels=vlab), outcome=factor(olab[outcome], levels=olab))]
pb <- ggplot(pb_d, aes(outcome, ratio, color=version)) +
  geom_hline(yintercept=1, color='grey55', linewidth=.35) +
  geom_linerange(aes(ymin=lo, ymax=hi), position=position_dodge(width=.6), linewidth=.6) +
  geom_point(position=position_dodge(width=.6), size=2) +
  scale_y_continuous(trans='log2', breaks=c(1,1.5,2,3,4,6)) +
  scale_color_manual(values=c('#1B4F72','#5DADE2','#E67E22')) +
  guides(color=guide_legend(ncol=2, byrow=TRUE)) +
  labs(title='Higher vs Lower Score, Q1–Q3', x=NULL, y='Ratio of geometric mean fold-rises') + th
fig <- (pa | pb) + plot_layout(widths=c(1,1.25)) + plot_annotation(tag_levels='A', tag_suffix='.')
ggsave(file.path(FIG_DIR, "E12_score_15_vs_30.png"), fig, width=7.2, height=3.6, dpi=600, bg='white')
cat("Saved", file.path(FIG_DIR, "E12_score_15_vs_30.png"), "\n")
