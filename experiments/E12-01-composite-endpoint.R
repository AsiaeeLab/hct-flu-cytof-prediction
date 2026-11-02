#!/usr/bin/env Rscript
# Composite antibody endpoint and same-baseline score effects.
#
# Composite outcome: mean of the three antigen log fold-rises (H1N1, H3N2,
# B/Victoria), i.e. the log of the geometric-mean fold-rise across antigens.
#
# 1. Four-model ladder on the composite outcome. Each fold reuses the atlas and
#    resolution of the primary analysis and refits single-outcome elastic nets
#    (alpha = 0.9, lambda.min, the same grouped inner folds as the primary fits).
#    Recipient-level paired bootstrap as in E3.
# 2. Composite cellular score: held-out CLR features times the composite M3
#    cluster coefficients (intercept and titer terms excluded), with the same
#    within-quartile median split and recipient bootstrap as E6-01.
# 3. Same-baseline effect of the cellular score (per antigen and composite):
#    log fold-rise regressed on the fold-centred score, the baseline log titer
#    and fold indicators. At a fixed baseline titer, a log fold-rise difference
#    is a log post-vaccination titer difference, so exp(slope x IQR) is the ratio
#    of post-vaccination titers between episodes at the 75th and 25th score
#    percentiles with the same starting titer.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({library(data.table); library(glmnet); library(here);
  library(ggplot2); library(patchwork)})
source(here("methods", "35-per-fold-atlas-cv.R"))

ROOT <- FLU_ROOT
FOLD_ROOT <- file.path(ROOT, "per_fold_atlas_5000")
STRAT_FILE <- file.path(ROOT, "post_hoc_5000/stratified/E2_stratified.RData")
OUT <- file.path(ROOT, "post_hoc_5000/composite"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FIG_DIR <- "/home/amir/repo/flu-cellular-prediction/paper/images"
V <- c("h1n1","h3n2","bvic")
B <- 2000L
RESP <- log(4) - 1e-12

# ---------------------------------------------------------------------------
# 1. Refit the ladder on the composite outcome, fold by fold
# ---------------------------------------------------------------------------
rows <- list(); fold_info <- list()
for (f in 1:5) {
  e <- new.env(); load(file.path(FOLD_ROOT, paste0("fold_",f), sprintf("fold_%d_all.RData", f)), e)
  tr <- e$train_p_ids; te <- e$test_p_ids
  stopifnot(identical(rownames(e$Y_test), te), identical(rownames(e$clr_test), te))
  yb_tr <- rowMeans(e$Y_train); yb_te <- rowMeans(e$Y_test)
  fid <- compute_grouped_foldid(tr, 5L, 7100 + f*1000 + 99)
  pred <- list(); cm3 <- NULL
  for (m in paste0("M",1:4)) {
    Xtr <- e[[paste0("X_",m,"_train")]]; Xte <- e[[paste0("X_",m,"_test")]]
    fit <- cv.glmnet(Xtr, yb_tr, family="gaussian", alpha=0.9, foldid=fid, standardize=TRUE)
    pred[[m]] <- as.vector(predict(fit, newx=Xte, s="lambda.min"))
    if (m == "M3") cm3 <- as.matrix(coef(fit, s="lambda.min"))[,1]
  }
  cl <- colnames(e$clr_test)
  cscore <- as.vector(e$clr_test %*% cm3[cl])
  ascore <- sapply(V, function(v) as.vector(e$clr_test %*% e$coef_M3_en[[v]][cl]))
  rows[[f]] <- data.table(episode_id=te, fold=f,
    y_h1n1=e$Y_test[,"h1n1"], y_h3n2=e$Y_test[,"h3n2"], y_bvic=e$Y_test[,"bvic"], y_comp=yb_te,
    base_h1n1=e$X_M1_test[,"baseline_h1n1"], base_h3n2=e$X_M1_test[,"baseline_h3n2"],
    base_bvic=e$X_M1_test[,"baseline_bvic"], base_comp=rowMeans(e$X_M1_test),
    pred_M1=pred$M1, pred_M2=pred$M2, pred_M3=pred$M3, pred_M4=pred$M4,
    score_comp=cscore, score_h1n1=ascore[,"h1n1"], score_h3n2=ascore[,"h3n2"], score_bvic=ascore[,"bvic"],
    avgmt_M1=rowMeans(e$pred_M1_en_test), avgmt_M3=rowMeans(e$pred_M3_en_test))
  fold_info[[f]] <- data.table(fold=f, k=ncol(e$clr_test),
    nonzero_composite=sum(cm3[cl]!=0), nonzero_multitask=sum(e$coef_M3_en$h3n2[cl]!=0))
}
dt <- rbindlist(rows); fold_info <- rbindlist(fold_info)
print(fold_info)

# Composite scores must equal those of the null-bootstrap observed fit.
obs_file <- file.path(ROOT, "score_null_bootstrap_composite/observed_composite_scores.rds")
if (file.exists(obs_file)) {
  o <- readRDS(obs_file)
  d <- max(abs(o$cscore[dt$episode_id] - dt$score_comp))
  cat("Max abs difference vs null-bootstrap observed composite scores:", d, "\n")
  stopifnot(d < 1e-9)
}

dt[, base_id := sub("_yr[0-9]+$", "", episode_id)]
mse <- sapply(paste0("M",1:4), function(m) mean((dt[[paste0("pred_",m)]] - dt$y_comp)^2))
cat("\nComposite MSE by model:\n"); print(round(mse, 4))
cat("Reference: composite of multitask predictions, M1 vs M3 MSE:",
    round(mean((dt$avgmt_M1-dt$y_comp)^2),4), round(mean((dt$avgmt_M3-dt$y_comp)^2),4), "\n")

# Recipient-level paired bootstrap (E3 convention)
recips <- unique(dt$base_id); rows_by <- split(seq_len(nrow(dt)), dt$base_id)[recips]
cmp <- data.table(comparison=c("M1_to_M2","M1_to_M3","M1_to_M4","M2_to_M3","M2_to_M4","M3_to_M4"),
                  a=c("M1","M1","M1","M2","M2","M3"), b=c("M2","M3","M4","M3","M4","M4"))
boot_res <- rbindlist(lapply(seq_len(nrow(cmp)), function(i) {
  ea <- (dt[[paste0("pred_",cmp$a[i])]]-dt$y_comp)^2; eb <- (dt[[paste0("pred_",cmp$b[i])]]-dt$y_comp)^2
  obs <- (mean(ea)-mean(eb))/mean(ea)*100
  set.seed(7000L + i*1000L + 1L)
  bd <- replicate(B, {ix <- unlist(rows_by[sample.int(length(recips), replace=TRUE)], use.names=FALSE)
                      (mean(ea[ix])-mean(eb[ix]))/mean(ea[ix])*100})
  ci <- quantile(bd, c(.025,.975), names=FALSE)
  data.table(comparison=cmp$comparison[i], delta_pct=obs, ci_low=ci[1], ci_high=ci[2],
             p_two_sided=min(1, 2*min(mean(bd<=0), mean(bd>=0))))
}))
cat("\nComposite ladder comparisons:\n"); print(boot_res)
fwrite(data.table(model=names(mse), mse=mse), file.path(OUT, "composite_ladder_mse.csv"))
fwrite(boot_res, file.path(OUT, "composite_ladder_bootstrap.csv"))

# ---------------------------------------------------------------------------
# 2. Score groups within baseline-titer quartiles (E6-01 convention)
# ---------------------------------------------------------------------------
s <- new.env(); load(STRAT_FILE, s)
qv <- as.character(s$quartile); names(qv) <- names(s$avg_baseline)
dt[, quartile := factor(qv[episode_id], levels=c("Q1_low","Q2","Q3","Q4_high"))]
stopifnot(!anyNA(dt$quartile))
outs <- c(V, "comp")
for (o in outs) {
  sc <- paste0("score_",o)
  dt[, (paste0("group_",o)) := ifelse(get(sc) > median(get(sc)), "High", "Low"), by=quartile]
}
base_levels <- sort(unique(dt$base_id))
rows_by_base <- split(seq_len(nrow(dt)), factor(dt$base_id, levels=base_levels))
set.seed(20260728L)
boot_index <- lapply(seq_len(B), function(b)
  unlist(rows_by_base[sample.int(length(base_levels), length(base_levels), replace=TRUE)], use.names=FALSE))
ci_of <- function(x) {x <- x[is.finite(x)]; quantile(x, c(.025,.975), names=FALSE)}
summ <- function(hi, lo, y, o, stratum) {
  gm <- function(r) exp(mean(y[r])); rr <- function(r) mean(y[r] >= RESP)
  inh <- logical(nrow(dt)); inh[hi] <- TRUE; inl <- logical(nrow(dt)); inl[lo] <- TRUE
  bt <- t(vapply(boot_index, function(ib) { bh <- ib[inh[ib]]; bl <- ib[inl[ib]]
    if (length(bh)<2 || length(bl)<2) return(rep(NA_real_,6))
    c(gm(bh), gm(bl), rr(bh), rr(bl), gm(bh)/gm(bl), rr(bh)-rr(bl)) }, numeric(6)))
  data.table(outcome=o, stratum=stratum, n_high=length(hi), n_low=length(lo),
    gmfr_high=gm(hi), gmfr_high_lo=ci_of(bt[,1])[1], gmfr_high_hi=ci_of(bt[,1])[2],
    gmfr_low=gm(lo), gmfr_low_lo=ci_of(bt[,2])[1], gmfr_low_hi=ci_of(bt[,2])[2],
    resp_high=rr(hi), resp_high_lo=ci_of(bt[,3])[1], resp_high_hi=ci_of(bt[,3])[2],
    resp_low=rr(lo), resp_low_lo=ci_of(bt[,4])[1], resp_low_hi=ci_of(bt[,4])[2],
    ratio=gm(hi)/gm(lo), ratio_lo=ci_of(bt[,5])[1], ratio_hi=ci_of(bt[,5])[2],
    diff_pp=100*(rr(hi)-rr(lo)), diff_lo=100*ci_of(bt[,6])[1], diff_hi=100*ci_of(bt[,6])[2])
}
grp <- rbindlist(lapply(outs, function(o) {
  y <- dt[[paste0("y_",o)]]; g <- dt[[paste0("group_",o)]]
  rbind(rbindlist(lapply(levels(dt$quartile), function(q)
          summ(which(dt$quartile==q & g=="High"), which(dt$quartile==q & g=="Low"), y, o, q))),
        summ(which(dt$quartile!="Q4_high" & g=="High"), which(dt$quartile!="Q4_high" & g=="Low"), y, o, "Q1-Q3"))
}))
cat("\nScore groups (pooled Q1-Q3):\n"); print(grp[stratum=="Q1-Q3", .(outcome, n_high, n_low, gmfr_high, gmfr_low, ratio, ratio_lo, ratio_hi, resp_high, resp_low, diff_pp)])
# The per-antigen ratios must reproduce the published values.
chk <- grp[stratum=="Q1-Q3" & outcome %in% V, setNames(round(ratio,2), outcome)]
cat("Per-antigen Q1-Q3 ratios (manuscript: 1.81, 2.21, 3.22):", chk, "\n")
if (Sys.getenv("FLU_RESULTS_ROOT") == "") stopifnot(isTRUE(all.equal(unname(chk), c(1.81,2.21,3.22), tolerance=0.006)))
fwrite(grp, file.path(OUT, "score_groups_composite_and_antigens.csv"))

# ---------------------------------------------------------------------------
# 3. Same-baseline effect of the cellular score
# ---------------------------------------------------------------------------
cont <- rbindlist(lapply(outs, function(o) {
  d <- data.table(y=dt[[paste0("y_",o)]], s=dt[[paste0("score_",o)]], b=dt[[paste0("base_",o)]], fold=factor(dt$fold))
  d[, sc := s - mean(s), by=fold]
  iqr <- diff(quantile(d$sc, c(.25,.75), names=FALSE))
  est <- function(ix) { m <- lm(y ~ sc + b + fold, data=d[ix]); coef(m)[["sc"]] * iqr }
  pt <- est(seq_len(nrow(d)))
  bt <- vapply(boot_index, function(ib) tryCatch(est(ib), error=function(err) NA_real_), numeric(1))
  ci <- ci_of(bt)
  data.table(outcome=o, iqr_score=iqr, log_effect=pt, titer_ratio=exp(pt),
             ratio_lo=exp(ci[1]), ratio_hi=exp(ci[2]), spearman=cor(d$sc, d$y, method="spearman"))
}))
cat("\nSame-baseline titer ratio, 75th vs 25th score percentile:\n"); print(cont)
fwrite(cont, file.path(OUT, "same_baseline_score_effect.csv"))
saveRDS(dt, file.path(OUT, "composite_episode_table.rds"))

# ---------------------------------------------------------------------------
# 4. Supplementary figure: composite score groups across baseline quartiles
# ---------------------------------------------------------------------------
blue <- '#1B4F72'; orange <- '#E67E22'
th <- theme_minimal(base_size=8, base_family='sans') + theme(
  panel.grid.minor=element_blank(), axis.text=element_text(size=6.5, color='grey20'),
  axis.title=element_text(size=7.5), plot.title=element_text(face='bold', size=9),
  plot.tag=element_text(face='bold', size=11), legend.position='bottom', legend.title=element_blank())
gl <- rbind(grp[outcome=="comp", .(stratum, group="High", n=n_high, gmfr=gmfr_high, lo=gmfr_high_lo, hi=gmfr_high_hi,
                                   r=resp_high, rlo=resp_high_lo, rhi=resp_high_hi)],
            grp[outcome=="comp", .(stratum, group="Low", n=n_low, gmfr=gmfr_low, lo=gmfr_low_lo, hi=gmfr_low_hi,
                                   r=resp_low, rlo=resp_low_lo, rhi=resp_low_hi)])
ql <- c("Q1_low","Q2","Q3","Q4_high","Q1-Q3")
gl[, x := c(1,2,3,4,5.6)[match(stratum, ql)]]
gl[, group := factor(group, levels=c("High","Low"))]
ns <- dcast(gl, stratum ~ group, value.var="n")[match(ql, stratum)]
xl <- paste0(c('Q1\nLowest','Q2\n','Q3\n','Q4\nHighest','Q1–Q3\nCombined'), '\nn=', ns$High, '/', ns$Low)
labs_g <- c('Higher composite score','Lower composite score')
base <- ggplot(gl, aes(x, color=group, shape=group, group=group)) +
  scale_color_manual(values=c(High=blue, Low=orange), labels=labs_g) +
  scale_shape_manual(values=c(High=16, Low=17), labels=labs_g) +
  scale_x_continuous(breaks=c(1,2,3,4,5.6), labels=xl) +
  geom_vline(xintercept=4.8, color='grey80', linewidth=.3) + th +
  theme(panel.grid.major.x=element_blank()) + labs(x='Baseline antibody quartile')
dodge <- position_dodge(width=.45)
p1 <- base + geom_hline(yintercept=4, linetype='dashed', color='grey55', linewidth=.35) +
  geom_linerange(aes(ymin=lo, ymax=hi), position=dodge, linewidth=.6) +
  geom_point(aes(y=gmfr), position=dodge, size=2.3) +
  scale_y_continuous(trans='log2', breaks=c(1,2,4,8,16,32,64)) +
  labs(y='Geometric mean fold-rise across antigens')
p2 <- base + geom_linerange(aes(ymin=100*rlo, ymax=100*rhi), position=dodge, linewidth=.6) +
  geom_point(aes(y=100*r), position=dodge, size=2.3) +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y='Composite fold-rise of at least 4 (%)')
fig <- (p1 | p2) + plot_layout(guides='collect') + plot_annotation(tag_levels='A', tag_suffix='.') &
  theme(legend.position='bottom')
ggsave(file.path(FIG_DIR, "E12_composite_score.png"), fig, width=7.2, height=3.8, dpi=600, bg='white')
cat("\nSaved figure and tables to", OUT, "\n")
