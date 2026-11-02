#!/usr/bin/env Rscript
# Top-down marker panel for the populations that drive the cellular score.
#
# The bottom-up analysis drops markers and re-clusters. This one starts from the
# ten populations that entered the six-population signature in at least one
# training fold (positive 4, 32, 16, 14; negative 3, 1, 29, 6, 8, 22) and asks how
# few markers identify them well enough to reproduce the cellular score.
#
# Reference score. Composite M3 (33 population CLR features plus the three log
# baseline titers) is fitted once on all 343 episodes (alpha 0.9, recipient-
# grouped 5-fold lambda CV, lambda.min). The score of an episode is the CLR
# vector times the population coefficients (titers and intercept excluded). It
# is an in-sample, descriptive score used only to weight populations here.
#
# Cell classifiers. Recipients are split 70/30 (stratified by cohort) into
# training and test sets; every episode of a recipient stays on one side.
# Greedy marker choices are made on a validation quarter of the training
# recipients, so the test cells play no part in choosing markers. Training and
# validation use at most 300 sampled cells per episode; test metrics use every
# sampled cell.
#   Per population: one-vs-rest logistic regression on the markers and their
#   squares, with the probability cutoff that maximizes training F1; greedy
#   forward selection up to 10 markers. Few populations reach F1 0.80, so each
#   population's set is also reported as the smallest one within 0.02 of its
#   best F1 (the plateau set).
#   Joint panel: greedy forward selection up to 15 markers on the macro-F1 of an
#   eleven-class labeling (the ten populations and "other"). For speed the greedy
#   search uses quadratic discriminant analysis over all 33 populations, collapsed
#   to the eleven classes; the chosen order is then scored on the test cells with
#   both QDA and a random forest.
#
# Score reproduction. For panels of the first k joint-greedy markers, a random
# forest assigns test cells to all 33 populations. Each test episode's sampled
# cells are recounted, CLR-transformed and scored with the reference
# coefficients, and compared with the score from the same cells under their
# 30-marker atlas labels. The atlas sample caps the larger populations at 8,000
# cells, so a second version weights each cell by the inverse sampling fraction
# of its atlas population to recover full-data proportions. The 15-marker
# consensus panel of the bottom-up analysis and all 30 markers (classifier
# ceiling) are run through the same pipeline. A second classifier, nearest
# atlas microcluster centroid on the panel markers (centroids from training
# cells) mapped to its population, is reported as a sensitivity check.
# Correlation intervals come from a recipient-level bootstrap of the test
# episodes (2,000 resamples).
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(ranger); library(parallel)
  library(MASS); library(FNN); library(ggplot2); library(patchwork); library(ggrepel); library(here)
})

ROOT <- FLU_ROOT
OUT <- file.path(ROOT, "post_hoc_5000/marker_topdown"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
BIO_FILE <- Sys.getenv("FLU_BIO_FILE", here("results", "final_atlas_production", "biology_production_90203bf3.RData"))
ANN_FILE <- Sys.getenv("FLU_ANN_FILE", here("results", "reports", "cluster_annotations_corrected.csv"))
FIG <- here("paper", "images", "E13_marker_topdown.png")
NCORE <- as.integer(Sys.getenv("N_CORES", "8")); CAP <- 300L; SEED <- 7099L
POS <- c(4L, 32L, 16L, 14L); NEG <- c(3L, 1L, 29L, 6L, 8L, 22L)
# First-season atlas: the populations that entered the six-population signature
# in at least one training fold, and that atlas's consensus 15-marker panel.
if (nzchar(Sys.getenv("FLU_RESULTS_ROOT"))) {
  sg <- fread(file.path(FLU_ROOT, "post_hoc_5000/shared_atlas_signature/signature_selection_frequency.csv"))
  POS <- sg[direction == "positive", as.integer(cluster_id)]; NEG <- sg[direction == "negative", as.integer(cluster_id)]
}
POPS <- c(POS, NEG)
CONSENSUS <- c("CD19","CD20","CD3","CD4","CD8","CD21","CXCR5","CCR6","CCR2","HLA-DR",
               "CD86","CD25","CD38","CD27","CD57")
if (nzchar(Sys.getenv("FLU_RESULTS_ROOT")))
  CONSENSUS <- fread(file.path(FLU_ROOT, "post_hoc_5000/marker_reduction/E2_marker_consensus.csv"))[order(rank)][1:15, marker]
K_POP <- 10L; K_JOINT <- 15L; K_SCORE <- c(4L, 6L, 8L, 10L, 12L, 15L)
RF_TREES <- 100L; RF_NODE <- 10L
t_start <- Sys.time()
msg <- function(...) cat(sprintf("[%6.1f min] ", as.numeric(Sys.time() - t_start, units="mins")), ..., "\n", sep="")

clr <- function(counts) { l <- log(counts + 0.5); l - rowMeans(l) }

# ---- data ----------------------------------------------------------------------
bio <- new.env(); load(BIO_FILE, bio)
ids <- bio$common_p_ids
counts_full <- bio$meta_counts_all; rownames(counts_full) <- ids
NPOP <- ncol(counts_full); pop_cols <- paste0("cluster_", 1:NPOP)
stopifnot(identical(colnames(counts_full), pop_cols))
clr33 <- clr(counts_full)
chk <- as.matrix(bio$clr_features); rownames(chk) <- ids
stopifnot(max(abs(clr33 - chk[, pop_cols])) < 1e-9)              # same CLR as the paper
stopifnot(all(CONSENSUS %in% bio$marker_names))
ann <- fread(ANN_FILE)

dt <- bio$dt_keep
X <- bio$X_scaled_keep
stopifnot(nrow(X) == nrow(dt), identical(colnames(X), bio$marker_names))
# the display sample also holds cells of episodes without complete titers; drop them
keep_cells <- dt$patient_id %in% ids; dt <- dt[keep_cells]; X <- X[keep_cells, , drop = FALSE]
MK <- bio$marker_names
cell_ep <- dt$patient_id
stopifnot(all(cell_ep %in% ids))
cell_pop <- dt$meta

comp <- readRDS(file.path(ROOT, "post_hoc_5000/composite/composite_episode_table.rds"))
comp <- comp[match(ids, episode_id)]
stopifnot(identical(comp$episode_id, ids))
cell_rec <- comp$base_id[match(cell_ep, ids)]

# ---- 1. reference score ----------------------------------------------------------
Xref <- cbind(clr33, as.matrix(comp[, .(base_h1n1, base_h3n2, base_bvic)]))
urec <- unique(comp$base_id)
set.seed(SEED)
rec_fold <- setNames(sample(rep(1:5, length.out=length(urec))), urec)
foldid <- as.integer(rec_fold[comp$base_id])
cvfit <- cv.glmnet(Xref, comp$y_comp, family="gaussian", alpha=0.9, standardize=TRUE, foldid=foldid)
beta_all <- as.matrix(coef(cvfit, s="lambda.min"))[, 1]
beta <- beta_all[pop_cols]
ref_score <- drop(clr33 %*% beta)
msg("reference score: lambda.min=", signif(cvfit$lambda.min, 3), ", ", sum(beta != 0),
    " nonzero population coefficients")
fwrite(data.table(population=1:NPOP, beta=beta,
                  annotation=ann$annotation[match(1:NPOP, ann$cluster_id)]),
       file.path(OUT, "reference_betas.csv"))

# ---- splits ----------------------------------------------------------------------
stopifnot(is.null(names(bio$cohort_vec)) || identical(names(bio$cohort_vec), ids))
rec_cohort <- unique(data.table(base_id=comp$base_id, cohort=bio$cohort_vec))
stopifnot(!anyDuplicated(rec_cohort$base_id))
set.seed(SEED + 1L)
test_rec <- unlist(lapply(split(rec_cohort$base_id, rec_cohort$cohort),
                          function(r) sample(r, round(0.3 * length(r)))))
train_rec <- setdiff(rec_cohort$base_id, test_rec)
tr_cohort <- rec_cohort[base_id %in% train_rec]
val_rec <- unlist(lapply(split(tr_cohort$base_id, tr_cohort$cohort),
                         function(r) sample(r, round(0.25 * length(r)))))
cap_cells <- function(idx) unlist(lapply(split(idx, cell_ep[idx]),
                                         function(i) if (length(i) > CAP) sample(i, CAP) else i),
                                  use.names=FALSE)
set.seed(SEED + 2L)
idx_test <- which(cell_rec %in% test_rec)                               # all sampled cells
idx_train <- cap_cells(which(cell_rec %in% train_rec))                  # final fits
idx_fit <- intersect(idx_train, which(!cell_rec %in% val_rec))          # greedy fits
idx_val <- intersect(idx_train, which(cell_rec %in% val_rec))           # greedy choices
test_ep <- ids[comp$base_id %in% test_rec]
msg("recipients train/test = ", length(train_rec), "/", length(test_rec),
    " (validation ", length(val_rec), "); episodes test = ", length(test_ep),
    "; cells fit/val/train/test = ", length(idx_fit), "/", length(idx_val), "/",
    length(idx_train), "/", length(idx_test))

# ---- classifiers -----------------------------------------------------------------
f1 <- function(y, p) { tp <- sum(y & p); if (tp == 0) 0 else 2 * tp / (sum(y) + sum(p)) }
macro_f1 <- function(yt, yp) mean(vapply(levels(yt), function(l) f1(yt == l, yp == l), 0))
qdesign <- function(idx, cols) { Z <- X[idx, cols, drop=FALSE]; cbind(1, Z, Z^2) }
fit_ovr <- function(cols, y, tr, te) {        # one-vs-rest quadratic logistic, F1 cutoff
  Ztr <- qdesign(tr, cols)
  b <- suppressWarnings(glm.fit(Ztr, y[tr], family=binomial()))$coefficients
  b[is.na(b)] <- 0
  s <- drop(Ztr %*% b); ytr <- y[tr] == 1
  th <- unique(quantile(s, seq(0.40, 0.999, length.out=300)))
  cut <- th[which.max(vapply(th, function(t) f1(ytr, s >= t), 0))]
  f1(y[te] == 1, drop(qdesign(te, cols) %*% b) >= cut)
}
fit_rf <- function(cols, y, tr, trees, threads) {
  ranger(x=X[tr, cols, drop=FALSE], y=y[tr], num.trees=trees, min.node.size=RF_NODE,
         num.threads=threads, seed=SEED, verbose=FALSE)
}
pred_rf <- function(m, cols, te, threads) predict(m, X[te, cols, drop=FALSE], num.threads=threads, seed=SEED)$predictions
greedy <- function(score_fn, kmax) {          # score_fn(cols) -> validation metric
  sel <- character(0); val <- numeric(0)
  for (k in seq_len(kmax)) {
    cand <- setdiff(MK, sel)
    sc <- unlist(mclapply(cand, function(m) score_fn(c(sel, m)), mc.cores=NCORE))
    sel <- c(sel, cand[which.max(sc)]); val <- c(val, max(sc))
  }
  list(markers=sel, val=val)
}

# ---- 2. per-population greedy ----------------------------------------------------
pop_rows <- list()
for (p in POPS) {
  y <- as.integer(cell_pop == p)
  g <- greedy(function(cols) fit_ovr(cols, y, idx_fit, idx_val), K_POP)
  test_f1 <- unlist(mclapply(seq_len(K_POP), function(k)
    fit_ovr(g$markers[seq_len(k)], y, idx_train, idx_test), mc.cores=NCORE))
  pop_rows[[length(pop_rows) + 1]] <- data.table(
    population=p, direction=ifelse(p %in% POS, "positive", "negative"),
    annotation=ann$annotation[match(p, ann$cluster_id)], n_markers=seq_len(K_POP),
    added_marker=g$markers, markers=vapply(seq_len(K_POP), function(k)
      paste(g$markers[seq_len(k)], collapse=";"), ""),
    f1_validation=g$val, f1_test=test_f1)
  msg("population ", p, ": ", paste(g$markers, collapse=" "), " | test F1 ",
      paste(sprintf("%.2f", test_f1), collapse=" "))
}
pop_f1 <- rbindlist(pop_rows)
# all-30-marker reference for each population
all30 <- unlist(mclapply(POPS, function(p) fit_ovr(MK, as.integer(cell_pop == p), idx_train, idx_test),
                         mc.cores=NCORE))
pop_f1[, f1_test_all30 := all30[match(population, POPS)]]
fwrite(pop_f1, file.path(OUT, "per_population_f1.csv"))
min_set <- function(d, thr) {                  # smallest set reaching thr, else the best set
  hit <- d[f1_test >= thr]
  if (nrow(hit)) hit[which.min(n_markers), .(k=n_markers, markers, f1=f1_test, reached=TRUE)]
  else d[which.max(f1_test), .(k=n_markers, markers, f1=f1_test, reached=FALSE)]
}
plateau_set <- function(d) {                   # smallest set within 0.02 of the best of 1-10 markers
  d[f1_test >= max(f1_test) - 0.02][which.min(n_markers), .(k_plateau=n_markers, markers_plateau=markers,
                                                            f1_plateau=f1_test, f1_best=max(d$f1_test))]
}
pop_min <- pop_f1[, c(list(annotation=annotation[1], f1_all30=f1_test_all30[1]),
                      setNames(min_set(.SD, 0.80), paste0(c("k","markers","f1","reached"), "_80")),
                      setNames(min_set(.SD, 0.90), paste0(c("k","markers","f1","reached"), "_90")),
                      plateau_set(.SD)),
                  by=.(population, direction)]
fwrite(pop_min, file.path(OUT, "per_population_min_sets.csv"))
union_plateau <- unique(unlist(strsplit(pop_min$markers_plateau, ";")))
msg("union of per-population plateau sets: ", length(union_plateau), " markers")

# ---- 3. joint eleven-class greedy -------------------------------------------------
y11 <- factor(ifelse(cell_pop %in% POPS, as.character(cell_pop), "other"),
              levels=c(as.character(POPS), "other"))
y33 <- factor(cell_pop, levels=1:NPOP)
pop11 <- function(lab) factor(ifelse(lab %in% POPS, as.character(lab), "other"), levels=levels(y11))
qda_pred <- function(cols, tr, te) {
  m <- MASS::qda(X[tr, cols, drop=FALSE], y33[tr])
  as.integer(as.character(predict(m, X[te, cols, drop=FALSE])$class))
}
gj <- greedy(function(cols) macro_f1(y11[idx_val], pop11(qda_pred(cols, idx_fit, idx_val))), K_JOINT)
joint_test_qda <- unlist(mclapply(seq_len(K_JOINT), function(k)
  macro_f1(y11[idx_test], pop11(qda_pred(gj$markers[seq_len(k)], idx_train, idx_test))), mc.cores=NCORE))
joint_test_rf <- vapply(seq_len(K_JOINT), function(k) {
  cols <- gj$markers[seq_len(k)]
  m <- fit_rf(cols, y11, idx_train, RF_TREES, NCORE)
  macro_f1(y11[idx_test], pred_rf(m, cols, idx_test, NCORE))
}, 0)
joint <- data.table(rank=seq_len(K_JOINT), marker=gj$markers, macro_f1_validation_qda=gj$val,
                    macro_f1_test_qda=joint_test_qda, macro_f1_test_rf=joint_test_rf,
                    in_consensus15=gj$markers %in% CONSENSUS,
                    in_union_plateau=gj$markers %in% union_plateau)
fwrite(joint, file.path(OUT, "joint_greedy_order.csv"))
msg("joint order: ", paste(gj$markers, collapse=" "))

# ---- 4. score reproduction ---------------------------------------------------------
micro2pop <- unique(dt[, .(micro, meta)]); stopifnot(!anyDuplicated(micro2pop$micro))
nearest_micro <- function(cols, tr, te) {
  C <- rowsum(X[tr, cols, drop=FALSE], dt$micro[tr]); C <- C / as.vector(table(dt$micro[tr])[rownames(C)])
  nn <- get.knnx(C, X[te, cols, drop=FALSE], k=1)$nn.index[, 1]
  micro2pop$meta[match(as.integer(rownames(C))[nn], micro2pop$micro)]
}
# Inverse atlas sampling fraction of each cell's 30-marker population, rescaled to
# mean 1 so weighted counts stay on the sample scale (same pseudocount meaning).
samp_frac <- as.numeric(table(y33)) / colSums(counts_full)
w_cell <- 1 / samp_frac[cell_pop]; w_cell <- w_cell / mean(w_cell)
te_ep_f <- factor(cell_ep[idx_test], levels=test_ep)
count_mat <- function(lab, w) {
  m <- as.matrix(xtabs(w ~ te_ep_f + lab, data.frame(te_ep_f, lab=factor(lab, levels=1:NPOP), w)))
  matrix(m, nrow(m), ncol(m), dimnames=list(test_ep, pop_cols))
}
score_from <- function(lab, weighted)
  drop(clr(count_mat(lab, if (weighted) w_cell[idx_test] else rep(1, length(idx_test)))) %*% beta)
orig_s <- score_from(cell_pop[idx_test], FALSE)
orig_w <- score_from(cell_pop[idx_test], TRUE)
y_te <- comp$y_comp[match(test_ep, ids)]
ref_te <- ref_score[test_ep]

# recipient-level bootstrap over test episodes for the correlations
B_BOOT <- 2000L
te_rows_by_rec <- split(seq_along(test_ep), comp$base_id[match(test_ep, ids)])
set.seed(SEED + 3L)
boot_rows <- lapply(seq_len(B_BOOT), function(b)
  unlist(te_rows_by_rec[sample(length(te_rows_by_rec), replace=TRUE)], use.names=FALSE))
boot_cor <- function(a, b) vapply(boot_rows, function(i) cor(a[i], b[i]), 0)

panels <- c(setNames(lapply(K_SCORE, function(k) gj$markers[seq_len(k)]), paste0("greedy_", K_SCORE)),
            list(consensus_15=CONSENSUS, all_30=MK))
rep_rows <- list(); boot_r <- list()
for (nm in names(panels)) for (clf in c("random_forest", "nearest_microcluster")) {
  cols <- panels[[nm]]
  lab <- if (clf == "random_forest") {
    m <- fit_rf(cols, y33, idx_train, RF_TREES, NCORE)
    as.integer(as.character(pred_rf(m, cols, idx_test, NCORE)))
  } else nearest_micro(cols, which(cell_rec %in% train_rec), idx_test)
  s <- score_from(lab, FALSE); sw <- score_from(lab, TRUE)
  key <- paste(nm, clf)
  boot_r[[key]] <- boot_cor(s, orig_s)
  rep_rows[[key]] <- data.table(
    panel=sub("_[0-9]+$", "", nm), classifier=clf, k=length(cols), markers=paste(cols, collapse=";"),
    pearson=cor(s, orig_s), pearson_lo=quantile(boot_r[[key]], 0.025, names=FALSE),
    pearson_hi=quantile(boot_r[[key]], 0.975, names=FALSE),
    spearman=cor(s, orig_s, method="spearman"),
    cor_with_outcome=cor(s, y_te),
    pearson_weighted=cor(sw, orig_w), spearman_weighted=cor(sw, orig_w, method="spearman"),
    cor_with_outcome_weighted=cor(sw, y_te),
    cell_accuracy_33=mean(lab == cell_pop[idx_test]),
    macro_f1_33=macro_f1(y33[idx_test], factor(lab, levels=1:NPOP)),
    macro_f1_10pops=mean(vapply(POPS, function(p) f1(cell_pop[idx_test] == p, lab == p), 0)))
  msg(key, ": r=", sprintf("%.3f", rep_rows[[key]]$pearson), " (weighted ",
      sprintf("%.3f", rep_rows[[key]]$pearson_weighted), "), accuracy ",
      sprintf("%.3f", rep_rows[[key]]$cell_accuracy_33))
}
rep_dt <- rbindlist(rep_rows)
boot_ref <- boot_cor(orig_s, ref_te)
ref_rows <- data.table(
  panel=c("atlas_labels_30", "full_data_reference"), classifier="atlas", k=30L, markers=paste(MK, collapse=";"),
  pearson=c(1, cor(orig_s, ref_te)), pearson_lo=c(1, quantile(boot_ref, 0.025, names=FALSE)),
  pearson_hi=c(1, quantile(boot_ref, 0.975, names=FALSE)),
  spearman=c(1, cor(orig_s, ref_te, method="spearman")),
  cor_with_outcome=c(cor(orig_s, y_te), cor(ref_te, y_te)),
  pearson_weighted=c(1, cor(orig_w, ref_te)), spearman_weighted=c(1, cor(orig_w, ref_te, method="spearman")),
  cor_with_outcome_weighted=c(cor(orig_w, y_te), cor(ref_te, y_te)),
  cell_accuracy_33=NA_real_, macro_f1_33=NA_real_, macro_f1_10pops=NA_real_)
rep_dt <- rbind(rep_dt, ref_rows)
fwrite(rep_dt, file.path(OUT, "score_reproduction.csv"))
# paired bootstrap: greedy panel minus 15-marker consensus (random forest)
diff_dt <- rbindlist(lapply(K_SCORE, function(kk) {
  d <- boot_r[[paste0("greedy_", kk, " random_forest")]] - boot_r[["consensus_15 random_forest"]]
  data.table(k=kk, diff=rep_dt[panel == "greedy" & classifier == "random_forest" & k == kk, pearson] -
               rep_dt[panel == "consensus" & classifier == "random_forest", pearson],
             lo=quantile(d, 0.025, names=FALSE), hi=quantile(d, 0.975, names=FALSE))
}))
fwrite(diff_dt, file.path(OUT, "greedy_minus_consensus_bootstrap.csv"))
fwrite(data.table(episode_id=test_ep, y_comp=y_te, ref_score=ref_te, sample_score=orig_s,
                  sample_score_weighted=orig_w, n_cells=as.integer(table(te_ep_f))),
       file.path(OUT, "test_episode_scores.csv"))

# ---- 5. figure -------------------------------------------------------------------
NAVY <- "#1B4F72"; ORANGE <- "#E67E22"
short <- c(`4`="CD21-low B", `32`="Activated DNT", `16`="EM CD4 T CD25+", `14`="Naive CD4 T CD25\u2212",
           `3`="CD38+ lin\u2212", `1`="Exh. EM CD8 T", `29`="Activated lin\u2212",
           `6`="EM CD4 T CD25\u2212", `8`="CM CD8 T CCR4+", `22`="EMRA CD8 T TIGIT\u2212",
           `27`="Naive CD8 T", `17`="Cycling mono", `11`="Activated APC")
pa_dt <- copy(pop_f1)[, `:=`(lab=paste0(population, " ", short[as.character(population)]),
                             Direction=factor(ifelse(direction == "positive", "Positive weight", "Negative weight"),
                                              levels=c("Positive weight", "Negative weight")))]
end_dt <- pa_dt[n_markers == K_POP]
th <- theme_minimal(base_size=8) +
  theme(panel.grid.minor=element_blank(), axis.text=element_text(size=6.5, color="grey20"),
        plot.tag=element_text(face="bold"))
pA <- ggplot(pa_dt, aes(n_markers, f1_test, group=population, color=Direction)) +
  annotate("segment", x=1, xend=K_POP, y=c(0.8, 0.9), yend=c(0.8, 0.9), linetype="dotted",
           color="grey55", linewidth=0.3) +
  geom_line(linewidth=0.45) + geom_point(size=0.7) +
  geom_text_repel(data=end_dt, aes(label=lab), size=2, direction="y", hjust=0, nudge_x=0.4,
                  xlim=c(K_POP + 0.3, K_POP + 6), segment.size=0.2, segment.color="grey60",
                  min.segment.length=0, box.padding=0.08, seed=1, show.legend=FALSE) +
  scale_color_manual(values=c("Positive weight"=NAVY, "Negative weight"=ORANGE), name=NULL) +
  scale_x_continuous(breaks=1:K_POP, limits=c(1, K_POP + 6), expand=expansion(add=c(0.3, 0))) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  coord_cartesian(clip="off") +
  labs(x="Number of markers", y=expression("Held-out "*F[1])) + th +
  theme(legend.position=c(0.36, 0.12), legend.key.height=unit(7, "pt"),
        legend.background=element_rect(fill="white", color=NA))
gr <- rep_dt[panel == "greedy" & classifier == "random_forest"]
cons <- rep_dt[panel == "consensus" & classifier == "random_forest"]
ceil <- rep_dt[panel == "all" & classifier == "random_forest"]
pb_dt <- rbind(gr[, .(k, r=pearson, lo=pearson_lo, hi=pearson_hi, series="Greedy top-down panel")],
               cons[, .(k=k + 0.35, r=pearson, lo=pearson_lo, hi=pearson_hi, series="15-marker consensus panel")])
pb_dt[, series := factor(series, levels=c("Greedy top-down panel", "15-marker consensus panel"))]
pB <- ggplot(pb_dt, aes(k, r, color=series, shape=series)) +
  geom_hline(yintercept=ceil$pearson, linetype="dashed", color="grey55", linewidth=0.3) +
  annotate("text", x=4, y=ceil$pearson, label="All 30 markers", vjust=-0.6, hjust=0, size=2.1, color="grey35") +
  geom_line(data=pb_dt[series == "Greedy top-down panel"], linewidth=0.5) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.25, linewidth=0.35) +
  geom_point(size=1.5, fill="white", stroke=0.7) +
  scale_color_manual(values=c("Greedy top-down panel"=NAVY, "15-marker consensus panel"=ORANGE), name=NULL) +
  scale_shape_manual(values=c("Greedy top-down panel"=16, "15-marker consensus panel"=23), name=NULL) +
  scale_x_continuous(breaks=K_SCORE) +
  scale_y_continuous(limits=c(min(0.5, min(pb_dt$lo)), 1), breaks=seq(0.5, 1, 0.1)) +
  labs(x="Number of markers", y=expression("Pearson "*italic(r)*" with 30-marker score")) + th +
  theme(legend.position=c(0.62, 0.14), legend.key.height=unit(7, "pt"),
        legend.background=element_rect(fill="white", color=NA))
fig <- (pA | pB) + plot_layout(widths=c(1.25, 1)) + plot_annotation(tag_levels="A", tag_suffix=".")
ggsave(FIG, fig, width=7.2, height=3.2, dpi=600, bg="white")
msg("figure: ", FIG)

# ---- 6. summary ------------------------------------------------------------------
f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x)
ci <- function(d) sprintf("%s (95%% CI %s to %s)", f2(d$pearson), f2(d$pearson_lo), f2(d$pearson_hi))
nz <- which(beta != 0)
pm <- pop_min[order(direction != "positive", -f1_best)]
pop_lines <- pm[, sprintf("- %d %s (%s): %s gives F1 %s (best of 1-10 markers %s; all 30 markers %s).",
                          population, annotation, direction, gsub(";", ", ", markers_plateau),
                          f2(f1_plateau), f2(f1_best), f2(f1_all30))]
reach80 <- pm[reached_80 == TRUE]
rep_tab <- rbind(
  gr[, .(panel=paste0("Greedy, first ", k), r=ci(.SD), rho=f2(spearman), r_w=f2(pearson_weighted),
         acc=f2(cell_accuracy_33), y=f2(cor_with_outcome)), by=k][, k := NULL],
  cons[, .(panel="Consensus 15", r=ci(.SD), rho=f2(spearman), r_w=f2(pearson_weighted),
           acc=f2(cell_accuracy_33), y=f2(cor_with_outcome))],
  ceil[, .(panel="All 30 (random forest)", r=ci(.SD), rho=f2(spearman), r_w=f2(pearson_weighted),
           acc=f2(cell_accuracy_33), y=f2(cor_with_outcome))])
nm30 <- rep_dt[panel == "all" & classifier == "nearest_microcluster"]
nmg <- rep_dt[panel == "greedy" & classifier == "nearest_microcluster"]
g6 <- gr[k == 6]; g15 <- gr[k == 15]; d15 <- diff_dt[k == 15]; d6 <- diff_dt[k == 6]
ref_row <- rep_dt[panel == "full_data_reference"]; smp_row <- rep_dt[panel == "atlas_labels_30"]
sm <- c(
  "# Top-down marker panel for the score-driving populations (E13-05)", "",
  "Script: `experiments/E13-05-marker-topdown.R`. Figure: `paper/images/E13_marker_topdown.png`.", "",
  "**Reference score.** Composite M3 (33 population CLR features plus three log baseline titers) fitted once on",
  sprintf("all 343 episodes (alpha 0.9, recipient-grouped 5-fold CV, lambda.min = %s). %d populations have nonzero",
          signif(cvfit$lambda.min, 3), length(nz)),
  sprintf("coefficients (%s). The score is the CLR vector times these coefficients. It is an in-sample,",
          paste(nz, collapse=", ")),
  "descriptive score used only to weight populations; it is not a validated prediction.", "",
  sprintf("**Design.** %d test recipients (%d episodes, %s sampled cells) and %d training recipients (%d of them",
          length(test_rec), length(test_ep), format(length(idx_test), big.mark=","), length(train_rec), length(val_rec)),
  "used only to choose markers). At most 300 cells per episode for fitting; every sampled test cell for evaluation.",
  "Cells come from the atlas sample (about 500 per episode), not the full data.", "",
  "**Per-population identification** (one-vs-rest logistic regression on markers and their squares, greedy",
  sprintf("forward selection). No population reaches held-out F1 0.90. %s",
          if (nrow(reach80)) paste0("Only ", paste(sprintf("%d (%d markers, F1 %s)", reach80$population, reach80$k_80, f2(reach80$f1_80)), collapse=" and "),
                                    if (nrow(reach80) == 1) " reaches" else " reach", " 0.80.") else "None reaches 0.80."),
  sprintf("With all 30 markers the same classifier reaches F1 %s to %s, so the atlas boundaries are not simple gates.",
          f2(min(pm$f1_all30)), f2(max(pm$f1_all30))),
  "Smallest set within 0.02 of each population's best F1:", pop_lines,
  sprintf("The union of these sets has %d markers.", length(union_plateau)), "",
  sprintf("**Joint panel** (eleven classes, greedy on macro-F1): %s. Held-out macro-F1 (random forest) is %s at 6,",
          paste(gj$markers, collapse=", "), f2(joint_test_rf[6])),
  sprintf("%s at 10 and %s at 15 markers.", f2(joint_test_rf[10]), f2(joint_test_rf[15])), "",
  "**Score reproduction** (random forest on 33 populations, test episodes; r = Pearson with the 30-marker score",
  "from the same sampled cells; rho = Spearman; r_w = design-weighted counts; acc = 33-class cell accuracy;",
  "y = Pearson with the composite outcome):", "",
  "| Panel | r | rho | r_w | acc | y |", "|---|---|---|---|---|---|",
  rep_tab[, sprintf("| %s | %s | %s | %s | %s | %s |", panel, r, rho, r_w, acc, y)], "",
  sprintf("The 30-marker sample score correlates %s with the outcome in these episodes (full-data score %s) and",
          f2(smp_row$cor_with_outcome), f2(ref_row$cor_with_outcome)),
  sprintf("%s with the full-data score, so about 500 sampled cells per episode already add noise of their own.",
          f2(ref_row$pearson)), "",
  sprintf("**Findings.** Six markers (%s) reproduce the 30-marker sample score at r = %s,",
          paste(gj$markers[1:6], collapse=", "), ci(g6)),
  sprintf("although they assign only %d%% of test cells to the right population; four markers give r = %s.",
          round(100 * g6$cell_accuracy_33), f2(gr[k == 4, pearson])),
  sprintf("Gains after six markers are small (r = %s at 10, %s at 15). At 15 markers the greedy panel and the",
          f2(gr[k == 10, pearson]), f2(g15$pearson)),
  sprintf("bottom-up consensus differ by %s (greedy minus consensus, 95%% CI %s to %s). At six markers the",
          f2(d15$diff), f2(d15$lo), f2(d15$hi)),
  if (d6$hi < 0) sprintf("greedy panel trails the 15-marker consensus by %s (95%% CI %s to %s).", f2(-d6$diff), f2(-d6$hi), f2(-d6$lo))
  else sprintf("difference from the 15-marker consensus is %s (95%% CI %s to %s).", f2(d6$diff), f2(d6$lo), f2(d6$hi)),
  sprintf("The random-forest ceiling with all 30 markers is %s; nearest-microcluster assignment gives %s with all",
          f2(ceil$pearson), f2(nm30$pearson)),
  sprintf("30 markers and %s to %s for the greedy panels. The populations themselves are not identifiable with",
          f2(min(nmg$pearson)), f2(max(nmg$pearson))),
  sprintf("five or six markers: held-out F1 at six markers ranges from %s to %s. About six markers approximate",
          f2(min(pop_f1[n_markers == 6, f1_test])), f2(max(pop_f1[n_markers == 6, f1_test]))),
  "the score; they do not identify the populations that define it.", "",
  "**Caveats.** In-sample reference coefficients; sampled cells rather than full counts (rare populations such",
  "as 32 have one or two sampled cells per episode); one 70/30 recipient split; greedy choices use QDA (joint)",
  "or quadratic logistic regression (per population), scoring uses a random forest. The atlas sample caps nine",
  "populations at 8,000 cells; r_w reweights cells by the inverse sampling fraction of their atlas population.")
writeLines(sm, file.path(OUT, "summary.md"))
saveRDS(list(beta=beta, gj=gj, pop_f1=pop_f1, pop_min=pop_min, rep_dt=rep_dt, diff_dt=diff_dt,
             splits=list(test=test_rec, val=val_rec)), file.path(OUT, "marker_topdown.rds"))
msg("done")
