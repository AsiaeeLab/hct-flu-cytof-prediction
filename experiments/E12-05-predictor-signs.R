#!/usr/bin/env Rscript
# Direction of the cellular predictors: pooled model and each cohort.
#
# 1. Pooled signs: the bootstrap-averaged standardized coefficients behind the
#    signed-importance figure (M4, 200 recipient bootstraps, shared atlas).
# 2. Per-cohort signs: repeats the per-cohort bootstrap of E5-per-cohort-refit.R
#    (M3 on the shared atlas, fixed reference penalty, 200 recipient bootstraps,
#    seed 7000) and additionally keeps the coefficients, so each population's
#    bootstrap-averaged standardized coefficient can be compared between
#    children and adults. The selection frequencies are checked against the
#    saved E5 results first.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
FLU_BOOT_FILE <- Sys.getenv("FLU_BOOT_FILE", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/bootstrap_importance_3f0f3ae8.RData")
suppressPackageStartupMessages({library(data.table); library(glmnet)})

ROOT <- FLU_ROOT
OUT <- file.path(ROOT, "post_hoc_5000/predictor_signs"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
BIO_FILE <- FLU_BIO_FILE
FOLD_FILE <- file.path(ROOT, "per_fold_atlas_5000/fold_1/fold_1_all.RData")
IMP_FILE <- FLU_BOOT_FILE
ANN_FILE <- FLU_ANN_FILE
V <- c("h1n1","h3n2","bvic")
ann <- fread(ANN_FILE)

# ---------------------------------------------------------------------------
# 1. Pooled signs from the signed-importance analysis
# ---------------------------------------------------------------------------
imp <- new.env(); load(IMP_FILE, imp)
st <- as.data.table(imp$stability)[grepl("^cluster_", feature)]
st[, cluster_id := as.integer(sub("cluster_", "", feature))]
pooled <- dcast(st, cluster_id ~ virus, value.var=c("mean_coef","selection_freq","sign_consistency"))
pooled <- merge(pooled, ann[, .(cluster_id, annotation)], by="cluster_id")
pooled[, direction := fifelse(mean_coef_h1n1 > 0 & mean_coef_h3n2 > 0 & mean_coef_bvic > 0, "positive in all three",
                      fifelse(mean_coef_h1n1 < 0 & mean_coef_h3n2 < 0 & mean_coef_bvic < 0, "negative in all three", "mixed"))]
cat("Pooled M4 signs for the populations in the Section 3.3 edit:\n")
print(pooled[cluster_id %in% c(4, 32, 16, 29, 1, 8),
             .(cluster_id, annotation, sel=selection_freq_h3n2, direction,
               h1n1=signif(mean_coef_h1n1,3), h3n2=signif(mean_coef_h3n2,3), bvic=signif(mean_coef_bvic,3))])
cat("\nPlasmablast (28) and transitional B (30):\n")
print(pooled[cluster_id %in% c(28, 30), .(cluster_id, annotation, sel=selection_freq_h3n2, direction)])
fwrite(pooled, file.path(OUT, "pooled_signs.csv"))

# ---------------------------------------------------------------------------
# 2. Per-cohort bootstrap with coefficients (E5 data set-up, verbatim)
# ---------------------------------------------------------------------------
e_bio <- new.env(parent=emptyenv()); load(BIO_FILE, envir=e_bio)
e_fold <- new.env(parent=emptyenv()); load(FOLD_FILE, envir=e_fold)
ids <- e_bio$common_p_ids; cohort <- e_bio$cohort_vec
clr_feat <- as.matrix(e_bio$clr_features)
all_Y <- rbind(e_fold$Y_train, e_fold$Y_test)
all_X_M1 <- rbind(e_fold$X_M1_train, e_fold$X_M1_test)
all_X_M3 <- cbind(clr_feat[rownames(all_Y), , drop=FALSE], all_X_M1)
ord <- match(ids, rownames(all_Y))
all_Y <- all_Y[ord, , drop=FALSE]; all_X_M3 <- all_X_M3[ord, , drop=FALSE]
stopifnot(identical(rownames(all_Y), ids))
base_id <- sub("_yr[0-9]+$", "", ids)

src <- readLines("/home/amir/repo/flu-cellular-prediction/experiments/E5-per-cohort-refit.R")
i0 <- grep("^inner_grouped_foldid <- function", src); i1 <- i0 + which(src[(i0+1):length(src)] == "}")[1]
eval(parse(text=src[i0:i1]))           # use the exact E5 helper

boot_with_coefs <- function(X, Y, p_ids, lambda_min, B = 200L, seed = 7000) {
  set.seed(seed)
  uniq <- unique(p_ids); rows_by_recip <- split(seq_len(nrow(X)), p_ids)
  fn <- colnames(X)
  counts <- matrix(0, length(fn), 3, dimnames=list(fn, V)); coefs <- counts
  for (b in seq_len(B)) {
    smpl <- sample(uniq, length(uniq), replace=TRUE)
    idx <- unlist(rows_by_recip[smpl], use.names=FALSE)
    fit <- glmnet(X[idx,,drop=FALSE], Y[idx,,drop=FALSE], family="mgaussian", alpha=0.9,
                  lambda=lambda_min, standardize=TRUE)
    cf <- coef(fit, s=lambda_min)
    for (v in V) {
      cv <- as.matrix(cf[[v]]); cv <- cv[rownames(cv) != "(Intercept)", , drop=FALSE]
      counts[, v] <- counts[, v] + (abs(cv[,1]) > 0)
      coefs[, v] <- coefs[, v] + cv[,1]
    }
  }
  list(freq=counts/B, mean_coef=coefs/B)
}

fit_cv_M3 <- cv.glmnet(all_X_M3, all_Y, family="mgaussian", alpha=0.9,
                       foldid=inner_grouped_foldid(base_id, 3, 7000), standardize=TRUE)
lambda_ref <- fit_cv_M3$lambda.min
saved <- new.env(); load(file.path(ROOT, "post_hoc_5000/per_cohort/E5_per_cohort.RData"), saved)
stopifnot(abs(lambda_ref - saved$lambda_ref) < 1e-12)
res <- list()
for (co in c("peds","adult")) {
  ix <- cohort == co
  r <- boot_with_coefs(all_X_M3[ix,,drop=FALSE], all_Y[ix,,drop=FALSE], base_id[ix], lambda_ref, B=200)
  ref <- if (co == "peds") saved$sel_peds else saved$sel_adult
  stopifnot(max(abs(r$freq[rownames(ref), ] - ref)) < 1e-12)    # reproduces E5 exactly
  sdv <- apply(all_X_M3[ix,,drop=FALSE], 2, sd)
  std <- sweep(r$mean_coef, 1, sdv, `*`)
  res[[co]] <- data.table(feature=rownames(std), cohort=co, freq=rowMeans(r$freq),
                          imp_h1n1=std[,1], imp_h3n2=std[,2], imp_bvic=std[,3], imp_mean=rowMeans(std))
}
cat("Per-cohort selection frequencies reproduce E5 exactly.\n")
pc <- rbindlist(res)[grepl("^cluster_", feature)]
pc[, cluster_id := as.integer(sub("cluster_", "", feature))]
w <- dcast(pc, cluster_id ~ cohort, value.var=c("freq","imp_mean","imp_h1n1","imp_h3n2","imp_bvic"))
w <- merge(w, ann[, .(cluster_id, annotation)], by="cluster_id")
w[, same_sign := sign(imp_mean_peds) == sign(imp_mean_adult)]
focus <- c(4, 32, 21, 16, 29, 1, 8)
cat("\nPopulations named in the Results (shared in both cohorts: 4, 32, 21, 16, 29; Spyros's list adds 1, 8):\n")
print(w[cluster_id %in% focus, .(cluster_id, annotation, freq_peds, freq_adult,
        imp_mean_peds=round(imp_mean_peds,3), imp_mean_adult=round(imp_mean_adult,3), same_sign)])
both70 <- w[freq_peds >= .7 & freq_adult >= .7]
cat("\nPopulations selected in >=70% of bootstraps in both cohorts:", nrow(both70),
    "| same direction in both:", sum(both70$same_sign), "\n")
fwrite(w, file.path(OUT, "per_cohort_signs.csv"))
saveRDS(list(per_cohort=pc, lambda_ref=lambda_ref), file.path(OUT, "per_cohort_signs.rds"))
