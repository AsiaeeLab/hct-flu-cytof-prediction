#!/usr/bin/env Rscript
# Cellular signature and lineage-negative sensitivity on the shared atlas.
#
# Both analyses use the shared 33-population atlas (the populations named in the
# paper) and the primary recipient-grouped outer folds. The atlas was built from
# all cells without outcomes, so these analyses evaluate prediction with a shared
# atlas, as in the per-cohort analysis.
#
# Signature. Inside each training fold only, M3 (titers plus the 33 population
# features) is refitted on 100 recipient-level bootstrap resamples at the fold's
# cross-validated penalty. Populations are ranked by their bootstrap-averaged
# standardized coefficient, averaged across the three antigens. The three most
# positive and the three most negative populations define an unweighted signature:
# the mean standardized CLR value of the positive set minus that of the negative
# set (standardized with training-fold means and SDs). The held-out fold is then
# scored, and titers plus the signature are compared with titers alone.
#
# Lineage-negative sensitivity. M3 is refitted after removing the four
# lineage-negative populations (clusters 3, 19, 24, 29) and recomputing the
# centered log-ratio on the remaining 29 populations, which is equivalent to
# excluding those cells before computing composition.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({library(data.table); library(glmnet); library(here)})
source(here("methods", "35-per-fold-atlas-cv.R"))

ROOT <- FLU_ROOT
FOLD_ROOT <- file.path(ROOT, "per_fold_atlas_5000")
OUT <- file.path(ROOT, "post_hoc_5000/shared_atlas_signature"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
BIO_FILE <- FLU_BIO_FILE
ANN_FILE <- FLU_ANN_FILE
V <- c("h1n1","h3n2","bvic"); B <- 2000L; B_SIG <- 100L; K_SIG <- 3L
LINEAGE_NEG <- paste0("cluster_", data.table::fread(ANN_FILE)[broad_lineage == "Other", cluster_id])  # lineage-negative populations

bio <- new.env(); load(BIO_FILE, bio)
ids <- bio$common_p_ids
counts <- bio$meta_counts_all
stopifnot(is.null(rownames(counts)) || identical(rownames(counts), ids))
rownames(counts) <- ids
clr33 <- transform_counts(counts, "clr")
chk <- as.matrix(bio$clr_features); rownames(chk) <- ids
stopifnot(max(abs(clr33 - chk[, colnames(clr33)])) < 1e-9)    # same CLR as the paper
keep <- setdiff(colnames(counts), LINEAGE_NEG)
clr29 <- transform_counts(counts[, keep], "clr")
ann <- fread(ANN_FILE)

comp <- readRDS(file.path(ROOT, "post_hoc_5000/composite/composite_episode_table.rds"))
agg <- new.env(); load(file.path(FOLD_ROOT, "aggregated/E2_aggregated.RData"), agg)

fit_pair <- function(Xtr, Ytr, Xte, fid) {
  m <- cv.glmnet(Xtr, Ytr, family="mgaussian", alpha=0.9, foldid=fid, standardize=TRUE)
  c <- cv.glmnet(Xtr, rowMeans(Ytr), family="gaussian", alpha=0.9, foldid=fid, standardize=TRUE)
  list(pred=as_pred_matrix(predict(m, newx=Xte, s="lambda.min"), V),
       cpred=as.vector(predict(c, newx=Xte, s="lambda.min")), lambda=m$lambda.min)
}

rank_populations <- function(X, Y, base, lambda, seed) {
  set.seed(seed)
  uniq <- unique(base); rows_by <- split(seq_len(nrow(X)), base)
  pops <- grep("^cluster_", colnames(X), value=TRUE)
  sdv <- apply(X[, pops, drop=FALSE], 2, sd)
  acc <- matrix(0, length(pops), 3, dimnames=list(pops, V))
  for (b in seq_len(B_SIG)) {
    idx <- unlist(rows_by[sample(uniq, length(uniq), replace=TRUE)], use.names=FALSE)
    fit <- glmnet(X[idx,,drop=FALSE], Y[idx,,drop=FALSE], family="mgaussian", alpha=0.9,
                  lambda=lambda, standardize=TRUE)
    cf <- coef(fit, s=lambda)
    for (v in seq_along(V)) acc[, v] <- acc[, v] + as.matrix(cf[[v]])[pops, 1]
  }
  imp <- rowMeans(sweep(acc / B_SIG, 1, sdv, `*`))
  pos <- names(sort(imp[imp > 0], decreasing=TRUE))[seq_len(min(K_SIG, sum(imp > 0)))]
  neg <- names(sort(imp[imp < 0]))[seq_len(min(K_SIG, sum(imp < 0)))]
  list(imp=imp, pos=pos, neg=neg)
}

signature_value <- function(clr_tr, clr_te, pos, neg) {
  mu <- colMeans(clr_tr); s <- apply(clr_tr, 2, sd)
  z <- function(M) sweep(sweep(M, 2, mu, `-`), 2, s, `/`)
  zt <- z(clr_tr); ze <- z(clr_te)
  f <- function(Z) rowMeans(Z[, pos, drop=FALSE]) - rowMeans(Z[, neg, drop=FALSE])
  list(train=f(zt), test=f(ze))
}

pred <- list(); sel <- list()
for (f in 1:5) {
  e <- new.env(); load(file.path(FOLD_ROOT, paste0("fold_",f), sprintf("fold_%d_all.RData", f)), e)
  tr <- e$train_p_ids; te <- e$test_p_ids
  fid <- compute_grouped_foldid(tr, 5L, 7100 + f*1000 + 99)
  btr <- sub("_yr[0-9]+$", "", tr)
  X33tr <- cbind(clr33[tr,], e$X_M1_train); X33te <- cbind(clr33[te,], e$X_M1_test)
  X29tr <- cbind(clr29[tr,], e$X_M1_train); X29te <- cbind(clr29[te,], e$X_M1_test)
  m33 <- fit_pair(X33tr, e$Y_train, X33te, fid)
  m29 <- fit_pair(X29tr, e$Y_train, X29te, fid)
  rk <- rank_populations(X33tr, e$Y_train, btr, m33$lambda, 20261011L + f)
  sg <- signature_value(clr33[tr,], clr33[te,], rk$pos, rk$neg)
  msg <- fit_pair(cbind(e$X_M1_train, signature=sg$train), e$Y_train,
                  cbind(e$X_M1_test, signature=sg$test), fid)
  pred[[f]] <- data.table(episode_id=te, fold=f,
    m33_h1n1=m33$pred[,1], m33_h3n2=m33$pred[,2], m33_bvic=m33$pred[,3], m33_comp=m33$cpred,
    m29_h1n1=m29$pred[,1], m29_h3n2=m29$pred[,2], m29_bvic=m29$pred[,3], m29_comp=m29$cpred,
    sig_h1n1=msg$pred[,1], sig_h3n2=msg$pred[,2], sig_bvic=msg$pred[,3], sig_comp=msg$cpred,
    signature=sg$test)
  sel[[f]] <- data.table(fold=f, population=c(rk$pos, rk$neg),
                         direction=c(rep("positive", length(rk$pos)), rep("negative", length(rk$neg))),
                         importance=rk$imp[c(rk$pos, rk$neg)])
  cat(sprintf("fold %d: positive %s | negative %s\n", f, paste(rk$pos, collapse=","), paste(rk$neg, collapse=",")))
}
pr <- rbindlist(pred); sel <- rbindlist(sel)
d <- merge(comp[, .(episode_id, base_id, quartile, y_h1n1, y_h3n2, y_bvic, y_comp, pred_M1)], pr, by="episode_id")
m1 <- agg$oof_pred_M1_en; rownames(m1) <- agg$common_p_ids
for (i in seq_along(V)) d[, (paste0("m1_", V[i])) := m1[episode_id, i]]
d[, m1_comp := pred_M1]

# Recipient-level paired bootstrap of percent MSE reduction versus titers alone
recips <- unique(d$base_id); rows_by <- split(seq_len(nrow(d)), d$base_id)[recips]
compare <- function(model, scope, seed) {
  if (scope == "overall") {
    ea <- rowSums(sapply(V, function(v) (d[[paste0("m1_",v)]] - d[[paste0("y_",v)]])^2))
    eb <- rowSums(sapply(V, function(v) (d[[paste0(model,"_",v)]] - d[[paste0("y_",v)]])^2))
  } else {
    ea <- (d$m1_comp - d$y_comp)^2; eb <- (d[[paste0(model,"_comp")]] - d$y_comp)^2
  }
  obs <- (sum(ea) - sum(eb)) / sum(ea) * 100
  set.seed(seed)
  bd <- replicate(B, {ix <- unlist(rows_by[sample.int(length(recips), replace=TRUE)], use.names=FALSE)
                      (sum(ea[ix]) - sum(eb[ix])) / sum(ea[ix]) * 100})
  ci <- quantile(bd, c(.025,.975), names=FALSE)
  data.table(model=model, scope=scope, delta_pct=obs, ci_low=ci[1], ci_high=ci[2],
             p_two_sided=min(1, 2*min(mean(bd<=0), mean(bd>=0))))
}
res <- rbindlist(list(
  compare("m33","overall",101L), compare("m29","overall",102L), compare("sig","overall",103L),
  compare("m33","composite",201L), compare("m29","composite",202L), compare("sig","composite",203L)))
res[, label := c(m33=sprintf("Titers + %d populations (shared atlas)", ncol(clr33)),
                 m29=sprintf("Titers + %d populations (lineage-negative removed)", length(keep)),
                 sig="Titers + 6-population signature")[model]]
cat("\nPercent MSE reduction versus titers alone:\n"); print(res)
fwrite(res, file.path(OUT, "signature_and_lineage_negative_mse.csv"))

# Signature score groups, pooled Q1-Q3 (E6 seed and row order)
d <- d[match(comp$episode_id, episode_id)]
base_levels <- sort(unique(d$base_id))
rows_by_base <- split(seq_len(nrow(d)), factor(d$base_id, levels=base_levels))
set.seed(20260728L)
boot_index <- lapply(seq_len(B), function(b)
  unlist(rows_by_base[sample.int(length(base_levels), length(base_levels), replace=TRUE)], use.names=FALSE))
hi_flag <- ave(d$signature, d$quartile, FUN=function(z) z > median(z)) == 1
q13 <- d$quartile != "Q4_high"
grp <- rbindlist(lapply(c(V, "comp"), function(o) {
  y <- d[[paste0("y_",o)]]; hi <- which(q13 & hi_flag); lo <- which(q13 & !hi_flag)
  inh <- logical(nrow(d)); inh[hi] <- TRUE; inl <- logical(nrow(d)); inl[lo] <- TRUE
  bt <- vapply(boot_index, function(ib) exp(mean(y[ib[inh[ib]]]) - mean(y[ib[inl[ib]]])), numeric(1))
  data.table(outcome=o, ratio=exp(mean(y[hi]) - mean(y[lo])), lo=quantile(bt,.025,names=FALSE),
             hi=quantile(bt,.975,names=FALSE))
}))
cat("\nSignature score groups, pooled Q1-Q3 fold-rise ratio:\n"); print(grp)
fwrite(grp, file.path(OUT, "signature_score_groups.csv"))

# Which populations were chosen, and how often
sel[, cluster_id := as.integer(sub("cluster_", "", population))]
sel <- merge(sel, ann[, .(cluster_id, annotation)], by="cluster_id", all.x=TRUE)
freq <- sel[, .(folds_selected=.N, mean_importance=mean(importance)), by=.(cluster_id, annotation, direction)][order(direction, -folds_selected)]
cat("\nSelection across the five training folds:\n"); print(freq)
fwrite(sel, file.path(OUT, "signature_selection_by_fold.csv"))
fwrite(freq, file.path(OUT, "signature_selection_frequency.csv"))

# Descriptive signature from all episodes (not used for evaluation)
allX <- cbind(clr33[comp$episode_id,], as.matrix(comp[, .(baseline_h1n1=base_h1n1, baseline_h3n2=base_h3n2, baseline_bvic=base_bvic)]))
allY <- as.matrix(comp[, .(h1n1=y_h1n1, h3n2=y_h3n2, bvic=y_bvic)])
fidall <- compute_grouped_foldid(comp$episode_id, 5L, 7000L)
lam_all <- cv.glmnet(allX, allY, family="mgaussian", alpha=0.9, foldid=fidall, standardize=TRUE)$lambda.min
rk_all <- rank_populations(allX, allY, comp$base_id, lam_all, 20261011L)
full <- data.table(population=names(rk_all$imp), importance=rk_all$imp)[order(-importance)]
full[, cluster_id := as.integer(sub("cluster_", "", population))]
full <- merge(full, ann[, .(cluster_id, annotation)], by="cluster_id")[order(-importance)]
full[, in_signature := population %in% c(rk_all$pos, rk_all$neg)]
cat("\nAll-episode ranking (descriptive):\n"); print(full[in_signature == TRUE])
fwrite(full, file.path(OUT, "signature_all_episodes_ranking.csv"))
