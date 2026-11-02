#!/usr/bin/env Rscript
# E16-02: Composite-outcome versions of the analyses that were computed only for
# the three-antigen (multitask) models. The composite outcome is the mean of the
# three antigen log fold-rises; every model below is a single-outcome elastic net
# (alpha 0.9, lambda.min) or OLS on that outcome. First-season episodes only.
#   A. Prediction by baseline-titer quartile (primary per-fold models).
#   B. Primary (pooled) models evaluated within each cohort.
#   C. Models refitted within each cohort on the shared atlas (E5 scheme), their
#      quartile gains, and bootstrap selection frequencies and standardized
#      coefficients of M3 in each cohort.
#   D. Training in one cohort and testing in the other (shared atlas).
#   E. Pairwise interactions among titers and clinical factors (per-fold atlases).
#   F. Illustrative assay-noise calculation.
#   G. Fitted responses under the two vaccine-dose models (random forests).
# Recipient-level bootstrap: B = 2,000, seed 20260728, two-sided p.
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(randomForest) })
ROOT <- Sys.getenv("FLU_RESULTS_ROOT"); stopifnot(nzchar(ROOT))
OUT <- file.path(ROOT, "post_hoc_5000/composite_main"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B <- 2000L
boot_gain <- function(ea, eb, seed = 20260728) {
  est <- 100 * (1 - sum(eb) / sum(ea)); set.seed(seed)
  bs <- replicate(B, { i <- sample.int(length(ea), replace = TRUE); 100 * (1 - sum(eb[i]) / sum(ea[i])) })
  data.table(gain = est, lo = quantile(bs, .025), hi = quantile(bs, .975), p = min(1, 2 * min(mean(bs <= 0), mean(bs >= 0))))
}
# ---- Primary per-fold composite predictions ---------------------------------
fold_rows <- rbindlist(lapply(1:5, function(f) {
  e <- new.env(); load(file.path(ROOT, sprintf("per_fold_atlas_5000/fold_%d/fold_%d_all.RData", f, f)), e)
  data.table(fold = f, episode_id = e$test_p_ids, y = rowMeans(e$Y_test), base = rowMeans(e$X_M1_test),
             M1 = e$comp_M1$test, M2 = e$comp_M2$test, M3 = e$comp_M3$test, M4 = e$comp_M4$test)
}))
fold_rows[, cohort := ifelse(grepl("^ASCT", episode_id), "Adults", "Children")]
stopifnot(nrow(fold_rows) %in% c(262L, 260L))
se <- function(m) (fold_rows[[m]] - fold_rows$y)^2
# A. Quartiles of the mean baseline log titer (same cut points as Table 2)
fold_rows[, quartile := cut(base, quantile(base, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))]
qa <- fold_rows[, .(n = .N, titer_lo = exp(min(base)), titer_hi = exp(max(base)),
                    M1 = mean((M1 - y)^2), M2 = mean((M2 - y)^2), M3 = mean((M3 - y)^2), M4 = mean((M4 - y)^2)), by = quartile][order(quartile)]
qa[, `:=`(gain_M1_M3 = 100 * (1 - M3 / M1), gain_M2_M4 = 100 * (1 - M4 / M2))]
qa <- cbind(qa, rbindlist(lapply(levels(fold_rows$quartile), function(q) {
  d <- fold_rows[quartile == q]; g <- boot_gain((d$M1 - d$y)^2, (d$M3 - d$y)^2)
  setnames(g, paste0(names(g), "_M1_M3")) })))
print(qa); fwrite(qa, file.path(OUT, "quartiles.csv"))
# B. Pooled models within each cohort
pb <- rbindlist(lapply(c("All", "Children", "Adults"), function(co) {
  d <- if (co == "All") fold_rows else fold_rows[cohort == co]
  rbind(cbind(cohort = co, comparison = "M1->M3", n = nrow(d), boot_gain((d$M1 - d$y)^2, (d$M3 - d$y)^2)),
        cbind(cohort = co, comparison = "M2->M4", n = nrow(d), boot_gain((d$M2 - d$y)^2, (d$M4 - d$y)^2)))
}))
print(pb); fwrite(pb, file.path(OUT, "pooled_within_cohort.csv"))

# ---- Shared atlas design -------------------------------------------------------
a <- new.env(); load(Sys.getenv("FLU_PROD_ATLAS"), a)
ids <- a$common_p_ids; yc <- rowMeans(a$Y_mat[ids, ]); cohort <- a$cohort_vec[ids]
clr <- as.matrix(a$clr_features[ids, ]); Xb <- as.matrix(a$X_base[ids, ])
titer_cols <- grep("^baseline_", colnames(Xb), value = TRUE)
X1 <- Xb[, titer_cols, drop = FALSE]; X2 <- Xb; X3 <- cbind(clr, X1); X4 <- cbind(clr, X2)
base_id <- sub("_yr[0-9]+$", "", ids)
stopifnot(length(ids) %in% c(262L, 260L), all(cohort %in% c("peds", "adult")))
# C. Per-cohort refits (E5 scheme: 5 recipient folds, seed 7000; 3 inner folds)
make_recipient_folds <- function(br, k = 5, seed = 7000) { set.seed(seed); u <- unique(br)
  fo <- sample(rep(seq_len(k), length.out = length(u))); names(fo) <- u; unname(fo[br]) }
inner_fid <- function(p, k = 3, seed) { set.seed(seed); u <- unique(p); fa <- sample(rep(seq_len(k), length.out = length(u))); names(fa) <- u; unname(fa[p]) }
run_cohort <- function(mask) {
  fo <- make_recipient_folds(base_id[mask])
  preds <- lapply(list(M1 = X1, M2 = X2, M3 = X3, M4 = X4), function(X) {
    Xc <- X[mask, , drop = FALSE]; yy <- yc[mask]; out <- numeric(sum(mask))
    for (f in 1:5) { te <- fo == f; tr <- !te
      fit <- cv.glmnet(Xc[tr, , drop = FALSE], yy[tr], family = "gaussian", alpha = 0.9,
                       foldid = inner_fid(base_id[mask][tr], 3, 7000 + f * 17), standardize = TRUE)
      out[te] <- as.vector(predict(fit, Xc[te, , drop = FALSE], s = "lambda.min")) }
    out })
  list(y = yc[mask], pred = preds, base = rowMeans(X1[mask, , drop = FALSE]))
}
pc <- list(Children = run_cohort(cohort == "peds"), Adults = run_cohort(cohort == "adult"))
pcr <- rbindlist(lapply(names(pc), function(co) { r <- pc[[co]]
  rbind(cbind(cohort = co, comparison = "M1->M3", n = length(r$y), boot_gain((r$pred$M1 - r$y)^2, (r$pred$M3 - r$y)^2)),
        cbind(cohort = co, comparison = "M2->M4", n = length(r$y), boot_gain((r$pred$M2 - r$y)^2, (r$pred$M4 - r$y)^2))) }))
print(pcr); fwrite(pcr, file.path(OUT, "per_cohort_refit.csv"))
pcq <- rbindlist(lapply(names(pc), function(co) { r <- pc[[co]]
  q <- cut(r$base, quantile(r$base, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
  rbindlist(lapply(levels(q), function(qq) { i <- q == qq
    data.table(cohort = co, quartile = qq, n = sum(i), M1 = mean((r$pred$M1[i] - r$y[i])^2), M3 = mean((r$pred$M3[i] - r$y[i])^2)) })) }))
pcq[, improvement_pct := 100 * (1 - M3 / M1)]; print(pcq); fwrite(pcq, file.path(OUT, "per_cohort_quartiles.csv"))
# Selection frequency and standardized coefficients of M3 in each cohort
set.seed(7000)
lam <- cv.glmnet(X3, yc, family = "gaussian", alpha = 0.9, foldid = inner_fid(base_id, 3, 7000), standardize = TRUE)$lambda.min
sel <- rbindlist(lapply(c(peds = "Children", adult = "Adults"), function(co) {
  m <- cohort == names(which(c(peds = "Children", adult = "Adults") == co))
  Xc <- X3[m, , drop = FALSE]; yy <- yc[m]; br <- base_id[m]; u <- unique(br); rows <- split(seq_len(nrow(Xc)), br)
  set.seed(7000); cnt <- cf <- setNames(numeric(ncol(Xc)), colnames(Xc))
  for (b in 1:200) { i <- unlist(rows[sample(u, length(u), replace = TRUE)], use.names = FALSE)
    fit <- glmnet(Xc[i, , drop = FALSE], yy[i], family = "gaussian", alpha = 0.9, lambda = lam, standardize = TRUE)
    bb <- as.matrix(coef(fit))[-1, 1]; cnt <- cnt + (bb != 0); cf <- cf + bb }
  data.table(cohort = co, feature = colnames(Xc), freq = cnt / 200, imp = (cf / 200) * apply(Xc, 2, sd)) }))
sel_w <- dcast(sel, feature ~ cohort, value.var = c("freq", "imp"))
fwrite(sel_w, file.path(OUT, "per_cohort_selection.csv")); cat("lambda for cohort selection:", lam, "\n")
# D. Transfer between cohorts (composite EN on titers + clinical + composition vs OLS on titers + clinical)
tr_res <- rbindlist(lapply(list(c("peds", "adult"), c("adult", "peds")), function(d) {
  tr <- cohort == d[1]; te <- cohort == d[2]
  ols <- lm.fit(cbind(1, X2[tr, ]), yc[tr]); cfo <- ols$coefficients; cfo[is.na(cfo)] <- 0
  pb_ <- drop(cbind(1, X2[te, ]) %*% cfo)
  set.seed(2027 + 500)
  fit <- cv.glmnet(X4[tr, ], yc[tr], family = "gaussian", alpha = 0.9, nfolds = 5)
  pc_ <- as.vector(predict(fit, X4[te, ], s = "lambda.min"))
  cbind(direction = paste(d, collapse = "_to_"), n_train = sum(tr), n_test = sum(te),
        boot_gain((pb_ - yc[te])^2, (pc_ - yc[te])^2))
}))
print(tr_res); fwrite(tr_res, file.path(OUT, "transfer.csv"))
# E. Interactions (per-fold atlases at the selected resolution; E3 seeds)
build_grouped_foldid <- function(p, k = 3, seed = 7000) { p <- sub("_yr[0-9]+$", "", p); set.seed(seed); u <- unique(p)
  fa <- sample(rep(seq_len(k), length.out = length(u))); names(fa) <- u; unname(fa[p]) }
ia <- rbindlist(lapply(1:5, function(f) {
  e <- new.env(); load(file.path(ROOT, sprintf("per_fold_atlas_5000/fold_%d/fold_%d_all.RData", f, f)), e)
  ytr <- rowMeans(e$Y_train); yte <- rowMeans(e$Y_test)
  Xi_tr <- model.matrix(~ .^2 - 1, data = as.data.frame(e$X_M2_train)); Xi_te <- model.matrix(~ .^2 - 1, data = as.data.frame(e$X_M2_test))
  keep <- apply(Xi_tr, 2, sd) > 0; Xi_tr <- Xi_tr[, keep]; Xi_te <- Xi_te[, keep]
  cl <- setdiff(colnames(e$X_M3_train), grep("^baseline_", colnames(e$X_M3_train), value = TRUE))
  set.seed(7000 + f * 17L)
  f2 <- cv.glmnet(Xi_tr, ytr, family = "gaussian", alpha = 0.9, foldid = build_grouped_foldid(e$train_p_ids, 3, 7000 + f * 13L), standardize = TRUE)
  f4 <- cv.glmnet(cbind(Xi_tr, e$X_M3_train[, cl]), ytr, family = "gaussian", alpha = 0.9,
                  foldid = build_grouped_foldid(e$train_p_ids, 3, 7000 + f * 13L), standardize = TRUE)
  data.table(fold = f, episode_id = e$test_p_ids, y = yte,
             M2int = as.vector(predict(f2, Xi_te, s = "lambda.min")),
             M4int = as.vector(predict(f4, cbind(Xi_te, e$X_M3_test[, cl]), s = "lambda.min")))
}))
ia <- merge(ia, fold_rows[, .(episode_id, M1, M2, M3, M4)], by = "episode_id")
mse <- sapply(c("M1", "M2", "M2int", "M3", "M4", "M4int"), function(m) mean((ia[[m]] - ia$y)^2))
itab <- data.table(model = names(mse), mse = mse, reduction_vs_M1 = 100 * (1 - mse / mse["M1"]))
itab <- rbind(itab, data.table(model = c("M2int->M4int", "M2int->M3"), mse = NA_real_,
                               reduction_vs_M1 = c(100 * (1 - mse["M4int"] / mse["M2int"]), 100 * (1 - mse["M3"] / mse["M2int"]))))
print(itab); fwrite(itab, file.path(OUT, "interactions.csv"))
# F. Noise floor (composite)
s2 <- log(2)^2; vy <- var(fold_rows$y); m <- sapply(c("M1", "M2", "M3", "M4"), function(x) mean(se(x)))
nf <- data.table(model = c(names(m), "noise"), mse = c(m, s2), rmse = sqrt(c(m, s2)), r2 = 1 - c(m, s2) / vy)
nf[, pct_ceiling := r2 / (1 - s2 / vy)]; attr(nf, "var_y") <- vy
print(nf); cat("Var(y_comp) =", vy, "; ceiling R2 =", 1 - s2 / vy, "; gap closed M1->M3 =", (m["M1"] - m["M3"]) / (m["M1"] - s2),
              "; with 2*(log2)^2:", 1 - 2 * s2 / vy, (m["M1"] - m["M3"]) / (m["M1"] - 2 * s2), "\n")
fwrite(cbind(nf, var_y = vy), file.path(OUT, "noise_floor.csv"))
# G. Fitted responses under the two dose models (descriptive, partly in-sample)
W <- a$covariates[ids, "dose_HD"]; Xrf <- cbind(clr, X1, a$covariates[ids, colnames(a$covariates) != "dose_HD", drop = FALSE])
set.seed(20270222)
hd <- randomForest(Xrf[W == 1, ], yc[W == 1], ntree = 500); sd_ <- randomForest(Xrf[W == 0, ], yc[W == 0], ntree = 500)
tl <- data.table(episode_id = ids, observed = yc, y_hat_hd = as.numeric(predict(hd, Xrf)), y_hat_sd = as.numeric(predict(sd_, Xrf)))
tl[, subgroup := ifelse(y_hat_hd < quantile(y_hat_hd, .25) & y_hat_sd < quantile(y_hat_sd, .25), "Lower predicted response",
                 ifelse(y_hat_hd > quantile(y_hat_hd, .75) & y_hat_sd > quantile(y_hat_sd, .75), "Higher predicted response", "Intermediate"))]
fwrite(tl, file.path(OUT, "dose_models.csv"))
cat("Done.\n")
