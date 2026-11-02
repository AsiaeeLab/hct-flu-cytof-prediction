# ==============================================================================
# Prediction Models
# ==============================================================================
# This file contains functions for different supervised learning models:
# linear, regularized (Lasso/Ridge/ElasticNet), and non-linear (Random Forest).

library(glmnet)      # For Lasso, Ridge, ElasticNet
library(randomForest) # For Random Forest


#' Fit linear regression model
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @return List with train and test predictions
fit_linear_model <- function(X_train, y_train, X_test) {
  # Convert to data frame for lm
  train_df <- as.data.frame(X_train)
  train_df$y <- y_train

  # Fit model
  model <- lm(y ~ ., data = train_df)

  # Predict
  train_pred <- predict(model, newdata = as.data.frame(X_train))
  test_pred <- predict(model, newdata = as.data.frame(X_test))

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    model = model,
    method = "linear"
  ))
}


#' Fit regularized regression model (Lasso/Ridge/ElasticNet)
#'
#' Uses cross-validation to select optimal lambda.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @param alpha Elastic net mixing parameter (0=Ridge, 1=Lasso, 0<alpha<1=ElasticNet)
#' @param nfolds Number of CV folds for lambda selection
#' @return List with train and test predictions
fit_regularized_model <- function(X_train, y_train, X_test,
                                   alpha = 1, nfolds = 10) {
  # Fit with cross-validation to select lambda
  cv_model <- cv.glmnet(X_train, y_train, alpha = alpha,
                        nfolds = nfolds, intercept = TRUE)

  # Predict using lambda.min
  train_pred <- as.numeric(predict(cv_model, newx = X_train, s = "lambda.min"))
  test_pred <- as.numeric(predict(cv_model, newx = X_test, s = "lambda.min"))

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    model = cv_model,
    lambda_min = cv_model$lambda.min,
    alpha = alpha,
    method = ifelse(alpha == 0, "ridge",
                   ifelse(alpha == 1, "lasso", "elasticnet"))
  ))
}


#' Fit regularized model with fixed lambda
#'
#' For grid search where lambda is predetermined.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @param alpha Elastic net mixing parameter
#' @param lambda Regularization parameter
#' @return List with train and test predictions
fit_regularized_fixed_lambda <- function(X_train, y_train, X_test,
                                         alpha = 1, lambda = 0.01) {
  # Fit with fixed lambda
  model <- glmnet(X_train, y_train, alpha = alpha,
                  lambda = lambda, intercept = TRUE)

  # Predict
  train_pred <- as.numeric(predict(model, newx = X_train, s = lambda))
  test_pred <- as.numeric(predict(model, newx = X_test, s = lambda))

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    model = model,
    lambda = lambda,
    alpha = alpha,
    method = ifelse(alpha == 0, "ridge",
                   ifelse(alpha == 1, "lasso", "elasticnet"))
  ))
}


#' Fit Random Forest model
#'
#' Non-linear ensemble method using decision trees.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @param ntree Number of trees
#' @param mtry Number of variables to try at each split (default: p/3)
#' @param nodesize Minimum size of terminal nodes
#' @return List with train and test predictions
fit_random_forest <- function(X_train, y_train, X_test,
                               ntree = 500, mtry = NULL, nodesize = 5) {
  # Default mtry for regression: p/3
  if (is.null(mtry)) {
    mtry <- max(floor(ncol(X_train) / 3), 1)
  }

  # Convert to data frame
  train_df <- as.data.frame(X_train)
  train_df$y <- y_train

  # Fit model
  rf_model <- randomForest(y ~ ., data = train_df,
                           ntree = ntree, mtry = mtry,
                           nodesize = nodesize,
                           importance = TRUE)

  # Predict
  train_pred <- predict(rf_model, newdata = as.data.frame(X_train))
  test_pred <- predict(rf_model, newdata = as.data.frame(X_test))

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    model = rf_model,
    importance = importance(rf_model),
    method = "random_forest"
  ))
}


#' Fit gradient boosting model (XGBoost)
#'
#' Non-linear ensemble method using gradient boosting.
#' Note: Requires xgboost package to be installed.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @param nrounds Number of boosting iterations
#' @param max_depth Maximum tree depth
#' @param eta Learning rate
#' @return List with train and test predictions
fit_xgboost <- function(X_train, y_train, X_test,
                        nrounds = 100, max_depth = 6, eta = 0.3) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost package is required but not installed.")
  }

  library(xgboost)

  # Prepare data
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgb.DMatrix(data = X_test)

  # Fit model
  xgb_model <- xgboost(data = dtrain,
                       nrounds = nrounds,
                       max_depth = max_depth,
                       eta = eta,
                       objective = "reg:squarederror",
                       verbose = 0)

  # Predict
  train_pred <- predict(xgb_model, dtrain)
  test_pred <- predict(xgb_model, dtest)

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    model = xgb_model,
    method = "xgboost"
  ))
}


#' Unified prediction interface
#'
#' Fits a model based on specified method and returns predictions.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param X_test Test feature matrix
#' @param method Prediction method: "linear", "lasso", "ridge", "elasticnet", "random_forest", "xgboost"
#' @param ... Additional parameters passed to specific model function
#' @return List with train and test predictions
fit_prediction_model <- function(X_train, y_train, X_test, method = "lasso", ...) {
  if (method == "linear") {
    return(fit_linear_model(X_train, y_train, X_test))
  } else if (method == "lasso") {
    return(fit_regularized_model(X_train, y_train, X_test, alpha = 1, ...))
  } else if (method == "ridge") {
    return(fit_regularized_model(X_train, y_train, X_test, alpha = 0, ...))
  } else if (method == "elasticnet") {
    return(fit_regularized_model(X_train, y_train, X_test, alpha = 0.5, ...))
  } else if (method == "random_forest") {
    return(fit_random_forest(X_train, y_train, X_test, ...))
  } else if (method == "xgboost") {
    return(fit_xgboost(X_train, y_train, X_test, ...))
  } else {
    stop("Unknown prediction method: ", method)
  }
}


#' Grid search for Elastic Net hyperparameters
#'
#' Tests a grid of lambda values and returns the best one based on CV error.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param alpha Elastic net mixing parameter (fixed at 0.5)
#' @param lambda_grid Vector of lambda values to test
#' @param n_inner_folds Number of CV folds for inner validation
#' @return List with best lambda and CV errors
grid_search_elasticnet <- function(X_train, y_train,
                                   alpha = 0.5,
                                   lambda_grid = 10^seq(-3, 1, length.out = 20),
                                   n_inner_folds = 5) {
  # Ensure y_train matches X_train rows
  if (!is.null(rownames(X_train))) {
    y_train <- y_train[rownames(X_train)]
  }

  # Create inner CV folds on row indices
  set.seed(123)
  n_samples <- nrow(X_train)
  inner_folds <- createFolds(1:n_samples, k = n_inner_folds, list = TRUE)

  cv_errors <- rep(NA, length(lambda_grid))

  # Test each lambda
  for (lambda_idx in seq_along(lambda_grid)) {
    lambda <- lambda_grid[lambda_idx]
    fold_errors <- numeric(n_inner_folds)

    # Inner CV
    for (fold in 1:n_inner_folds) {
      val_idx <- inner_folds[[fold]]
      train_idx <- setdiff(1:n_samples, val_idx)

      X_inner_train <- X_train[train_idx, , drop = FALSE]
      y_inner_train <- y_train[train_idx]
      X_inner_val <- X_train[val_idx, , drop = FALSE]
      y_inner_val <- y_train[val_idx]

      # Fit with fixed lambda
      model <- glmnet(X_inner_train, y_inner_train,
                     alpha = alpha, lambda = lambda, intercept = TRUE)
      pred <- as.numeric(predict(model, newx = X_inner_val, s = lambda))
      fold_errors[fold] <- mean((pred - y_inner_val)^2)
    }

    cv_errors[lambda_idx] <- mean(fold_errors)
  }

  # Select best lambda
  best_idx <- which.min(cv_errors)
  best_lambda <- lambda_grid[best_idx]

  return(list(
    best_lambda = best_lambda,
    cv_errors = cv_errors,
    lambda_grid = lambda_grid,
    best_idx = best_idx
  ))
}


#' Grid search for Random Forest hyperparameters
#'
#' Tests a grid of mtry and nodesize values.
#'
#' @param X_train Training feature matrix
#' @param y_train Training response vector
#' @param mtry_grid Vector of mtry values to test
#' @param nodesize_grid Vector of nodesize values to test
#' @param ntree Number of trees (fixed)
#' @param n_inner_folds Number of CV folds for inner validation
#' @return List with best hyperparameters and CV errors
grid_search_rf <- function(X_train, y_train,
                           mtry_grid = NULL,
                           nodesize_grid = c(1, 5, 10),
                           ntree = 500,
                           n_inner_folds = 3) {  # Fewer folds for RF (slow)

  # Ensure y_train matches X_train rows
  if (!is.null(rownames(X_train))) {
    y_train <- y_train[rownames(X_train)]
  }

  # Default mtry grid
  if (is.null(mtry_grid)) {
    p <- ncol(X_train)
    mtry_grid <- c(max(1, floor(sqrt(p))),
                  max(1, floor(p / 3)),
                  max(1, floor(p / 2)))
    mtry_grid <- unique(mtry_grid)
  }

  # Create inner CV folds on row indices
  set.seed(123)
  n_samples <- nrow(X_train)
  inner_folds <- createFolds(1:n_samples, k = n_inner_folds, list = TRUE)

  # Initialize results
  grid_results <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
  grid_results$cv_error <- NA

  # Test each combination
  for (i in 1:nrow(grid_results)) {
    mtry <- grid_results$mtry[i]
    nodesize <- grid_results$nodesize[i]
    fold_errors <- numeric(n_inner_folds)

    # Inner CV
    for (fold in 1:n_inner_folds) {
      val_idx <- inner_folds[[fold]]
      train_idx <- setdiff(1:n_samples, val_idx)

      X_inner_train <- X_train[train_idx, , drop = FALSE]
      y_inner_train <- y_train[train_idx]
      X_inner_val <- X_train[val_idx, , drop = FALSE]
      y_inner_val <- y_train[val_idx]

      # Fit RF
      train_df <- as.data.frame(X_inner_train)
      train_df$y <- y_inner_train

      rf_model <- randomForest(y ~ ., data = train_df,
                              ntree = ntree, mtry = mtry,
                              nodesize = nodesize)

      pred <- predict(rf_model, newdata = as.data.frame(X_inner_val))
      fold_errors[fold] <- mean((pred - y_inner_val)^2)
    }

    grid_results$cv_error[i] <- mean(fold_errors)
  }

  # Select best
  best_idx <- which.min(grid_results$cv_error)
  best_mtry <- grid_results$mtry[best_idx]
  best_nodesize <- grid_results$nodesize[best_idx]

  return(list(
    best_mtry = best_mtry,
    best_nodesize = best_nodesize,
    best_ntree = ntree,
    cv_errors = grid_results$cv_error,
    grid_results = grid_results,
    best_idx = best_idx
  ))
}


#' Calculate model performance metrics
#'
#' @param y_true True response values
#' @param y_pred Predicted response values
#' @return List with MSE, RMSE, MAE, and R-squared
calculate_metrics <- function(y_true, y_pred) {
  mse <- mean((y_true - y_pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_true - y_pred))

  # R-squared
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  r_squared <- 1 - (ss_res / ss_tot)

  return(list(
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    R2 = r_squared
  ))
}
