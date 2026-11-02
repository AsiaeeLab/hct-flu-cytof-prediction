#!/usr/bin/env Rscript
# How many populations carry the cellular prediction? (October 2026 round)
#
# Shared 33-population atlas and the primary recipient-grouped outer folds, as
# in E12-04. Inside each training fold only, populations are ranked by the
# magnitude of their bootstrap-averaged standardized coefficient (100 recipient
# bootstrap refits of M3 at the fold's cross-validated penalty, averaged across
# antigens). For k = 1, 2, ..., 33, M3 is refitted on the three baseline titers
# plus the k top-ranked populations (CLR computed over all 33), and the held-out
# fold is predicted. Reported: percent reduction in out-of-fold MSE relative to
# titers alone, overall across antigens and for the composite outcome, with
# recipient-level paired bootstrap 95% CIs and two-sided p-values.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({library(data.table); library(glmnet); library(here)})
source(here("methods", "35-per-fold-atlas-cv.R"))

ROOT <- FLU_ROOT
FOLD_ROOT <- file.path(ROOT, "per_fold_atlas_5000")
OUT <- file.path(ROOT, "post_hoc_5000/top_k_populations"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
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

KS <- unique(c(1:8, 10, 12, 15, 20, 25, 30, 35, ncol(clr33)))
pred <- list(); ranks <- list()
for (f in 1:5) {
  e <- new.env(); load(file.path(FOLD_ROOT, paste0("fold_",f), sprintf("fold_%d_all.RData", f)), e)
  tr <- e$train_p_ids; te <- e$test_p_ids
  fid <- compute_grouped_foldid(tr, 5L, 7100 + f*1000 + 99)
  btr <- sub("_yr[0-9]+$", "", tr)
  X33tr <- cbind(clr33[tr,], e$X_M1_train); X33te <- cbind(clr33[te,], e$X_M1_test)
  m33 <- fit_pair(X33tr, e$Y_train, X33te, fid)
  rk <- rank_populations(X33tr, e$Y_train, btr, m33$lambda, 20261011L + f)
  ord <- names(sort(abs(rk$imp), decreasing = TRUE))
  ranks[[f]] <- data.table(fold = f, population = ord, rank = seq_along(ord), importance = rk$imp[ord])
  for (k in KS) {
    keepk <- ord[seq_len(k)]
    mk <- fit_pair(cbind(clr33[tr, keepk, drop=FALSE], e$X_M1_train), e$Y_train,
                   cbind(clr33[te, keepk, drop=FALSE], e$X_M1_test), fid)
    pred[[length(pred) + 1]] <- data.table(episode_id = te, fold = f, k = k,
      p_h1n1 = mk$pred[,1], p_h3n2 = mk$pred[,2], p_bvic = mk$pred[,3], p_comp = mk$cpred)
  }
  cat(sprintf("fold %d done; top 6: %s\n", f, paste(ord[1:6], collapse = ",")))
}
pr <- rbindlist(pred); ranks <- rbindlist(ranks)
fwrite(ranks, file.path(OUT, "population_ranks_by_fold.csv"))
d0 <- comp[, .(episode_id, base_id, y_h1n1, y_h3n2, y_bvic, y_comp, pred_M1)]
m1 <- agg$oof_pred_M1_en; rownames(m1) <- agg$common_p_ids
for (i in seq_along(V)) d0[, (paste0("m1_", V[i])) := m1[episode_id, i]]
d0[, m1_comp := pred_M1]
recips <- unique(d0$base_id)
saveRDS(pr, file.path(OUT, "top_k_predictions.rds"))
res <- rbindlist(lapply(KS, function(kk) {
  k <- kk
  d <- merge(d0, pr[pr$k == kk], by = "episode_id"); rows_by <- split(seq_len(nrow(d)), d$base_id)[recips]
  rbindlist(lapply(c("overall", "composite"), function(scope) {
    if (scope == "overall") {
      ea <- rowSums(sapply(V, function(v) (d[[paste0("m1_",v)]] - d[[paste0("y_",v)]])^2))
      eb <- rowSums(sapply(V, function(v) (d[[paste0("p_",v)]] - d[[paste0("y_",v)]])^2))
    } else { ea <- (d$m1_comp - d$y_comp)^2; eb <- (d$p_comp - d$y_comp)^2 }
    obs <- (sum(ea) - sum(eb)) / sum(ea) * 100
    set.seed(300L + k)
    bd <- replicate(B, {ix <- unlist(rows_by[sample.int(length(recips), replace=TRUE)], use.names=FALSE)
                        (sum(ea[ix]) - sum(eb[ix])) / sum(ea[ix]) * 100})
    data.table(k = k, scope = scope, delta_pct = obs, ci_low = quantile(bd, .025), ci_high = quantile(bd, .975),
               p_two_sided = min(1, 2 * min(mean(bd <= 0), mean(bd >= 0))))
  }))
}))
fwrite(res, file.path(OUT, "top_k_mse.csv")); print(res)
