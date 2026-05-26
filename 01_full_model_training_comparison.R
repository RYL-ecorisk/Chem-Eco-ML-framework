# ==============================================================================
# 01_full_model_training_comparison.R
# ------------------------------------------------------------------------------
# Cleaned main script for full-model development and comparison.
# Computational settings are kept from the working analysis script:
#   - 10 x 10 repeated cross-validation
#   - Bayesian optimization search settings and bounds unchanged
#   - Candidate-pool and generalization-score rules unchanged
#   - Outputs retain the FINAL run tag
# ==============================================================================

# ==============================================================================
# Train and compare four machine-learning models with repeated cross-validation
# Models: Elastic Net, Support Vector Machine, Random Forest, and XGBoost
# Outputs: model performance, selected model, hyperparameters, and predictions
# ==============================================================================


# 0. Configuration --------------------------------------------------------------

RUN_TAG <- "FINAL"

OUT_FILE <- function(name, ext = "csv") {
  paste0(name, "_", RUN_TAG, ".", ext)
}

CFG <- list(
  cv_folds = 10,
  cv_repeats = 10,
  
  gap_good = 0.15,
  gap_ok = 0.25,
  low_cv_cutoff = 0.60,
  
  enet_penalty_gap = 0.12,
  enet_penalty_sd  = 0.05,
  enet_sd_cutoff   = 0.15,
  
  svm_penalty_gap = 0.18,
  svm_penalty_sd  = 0.06,
  svm_sd_cutoff   = 0.15,
  
  rf_penalty_gap = 0.25,
  rf_penalty_sd  = 0.05,
  rf_sd_cutoff   = 0.12,
  
  xgb_penalty_gap = 0.15,
  xgb_penalty_sd  = 0.05,
  xgb_sd_cutoff   = 0.15,
  xgb_bonus_gap_ok = 0.002,
  
  pool_pre_top_n = 20,
  pool_final_top_n = 10,
  pool_pre_delta_cv = 0.015,
  pool_final_delta_cv = 0.010
)


# 1. Packages -------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  caret,
  dplyr,
  xgboost,
  ranger,
  glmnet,
  kernlab,
  doParallel,
  foreach,
  rBayesianOptimization
)

try(parallel::stopCluster(cl), silent = TRUE)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

cat("Parallel cores for caret models:", n_cores, "\n")
cat("xgboost version:", as.character(packageVersion("xgboost")), "\n")
cat("Run tag:", RUN_TAG, "\n")


# 2. Data -----------------------------------------------------------------------

df <- read.csv("datasets.csv")

features <- c(
  "Solubility", "log.Kow", "HaCount", "MW", "Complexity", "Rings",
  "Taxavalue", "Respiration", "Locomotion", "Feeding", "Size", "HLC"
)

if (!all(features %in% names(df))) {
  missing_cols <- setdiff(features, names(df))
  stop("Missing predictor columns: ", paste(missing_cols, collapse = ", "))
}

if (!"Logtox" %in% names(df)) {
  stop("Response column Logtox is missing.")
}

X <- df[, features]
y <- df$Logtox

if (!all(sapply(X, is.numeric))) {
  non_numeric <- names(X)[!sapply(X, is.numeric)]
  stop("Non-numeric predictor columns: ", paste(non_numeric, collapse = ", "))
}

if (!is.numeric(y)) stop("The response variable must be numeric.")

X_mat <- as.matrix(X)

if (any(is.na(X_mat)) || any(is.infinite(X_mat))) {
  stop("Predictor matrix contains NA or Inf values.")
}

if (any(is.na(y)) || any(is.infinite(y))) {
  stop("Response variable contains NA or Inf values.")
}

cat("N =", nrow(X_mat), "\n")
cat("Number of predictors =", ncol(X_mat), "\n")


# 3. Metric functions -----------------------------------------------------------

calc_r2 <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  pred <- pred[ok]
  obs <- obs[ok]
  
  if (length(obs) < 2 || length(unique(obs)) <= 1 || length(unique(pred)) <= 1) {
    return(NA_real_)
  }
  
  r <- suppressWarnings(cor(pred, obs))
  if (!is.finite(r)) return(NA_real_)
  r^2
}

calc_rmse <- function(pred, obs) {
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}

calc_mae <- function(pred, obs) {
  mean(abs(pred - obs), na.rm = TRUE)
}

calc_nrmse <- function(pred, obs) {
  y_range <- diff(range(obs, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) return(NA_real_)
  calc_rmse(pred, obs) / y_range
}

gap_status <- function(gap) {
  ifelse(
    gap < 0.10, "Excellent",
    ifelse(
      gap < CFG$gap_good, "Good",
      ifelse(gap <= CFG$gap_ok, "Acceptable", "High Variance")
    )
  )
}

qualification_flag <- function(gap, cv_r2) {
  if (is.na(gap) || is.na(cv_r2)) return("Unknown")
  if (cv_r2 < CFG$low_cv_cutoff) return("Lower Predictive Power")
  if (gap < CFG$gap_good) return("Hard Qualified")
  if (gap <= CFG$gap_ok) return("Qualified")
  return("Needs More Tuning")
}

generalization_score <- function(cv_r2, gap_r2, cv_r2_sd,
                                 gap_penalty = 0.15,
                                 sd_penalty = 0.05,
                                 gap_cutoff = 0.20,
                                 sd_cutoff = 0.15,
                                 bonus_if_ok = 0) {
  score <- cv_r2 -
    gap_penalty * max(0, gap_r2 - gap_cutoff) -
    sd_penalty * max(0, cv_r2_sd - sd_cutoff)
  
  if (is.finite(gap_r2) && gap_r2 <= CFG$gap_ok) {
    score <- score + bonus_if_ok
  }
  
  if (!is.finite(score)) score <- -999
  score
}


# 4. Cross-validation setup -----------------------------------------------------

set.seed(123)

cv_index <- caret::createMultiFolds(
  y = y,
  k = CFG$cv_folds,
  times = CFG$cv_repeats
)

ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = CFG$cv_folds,
  repeats = CFG$cv_repeats,
  index = cv_index,
  savePredictions = "final",
  returnResamp = "final",
  summaryFunction = caret::defaultSummary,
  allowParallel = TRUE
)

cat("Number of CV resamples:", length(cv_index), "\n")


# 5. Helper functions -----------------------------------------------------------

append_history <- function(history_df, row_list) {
  rbind(history_df, as.data.frame(row_list, stringsAsFactors = FALSE))
}

dedup_integer_history <- function(history,
                                  int_cols = c("mtry", "min.node.size")) {
  keep <- int_cols[int_cols %in% names(history)]
  if (length(keep) == 0) return(history)
  
  tmp <- history
  for (nm in keep) tmp[[paste0(nm, "_int")]] <- round(tmp[[nm]])
  
  dedup_cols <- paste0(keep, "_int")
  tmp <- tmp[!duplicated(tmp[, dedup_cols, drop = FALSE]), , drop = FALSE]
  tmp[, dedup_cols] <- NULL
  tmp
}

build_candidate_pool <- function(history,
                                 cv_col = "CV_R2",
                                 gap_col = "Gap_R2",
                                 sd_col = "CV_R2_SD",
                                 pre_top_n = CFG$pool_pre_top_n,
                                 final_top_n = CFG$pool_final_top_n,
                                 pre_delta_cv = CFG$pool_pre_delta_cv,
                                 final_delta_cv = CFG$pool_final_delta_cv) {
  history <- history[
    is.finite(history[[cv_col]]) & is.finite(history[[gap_col]]),
    ,
    drop = FALSE
  ]
  
  if (nrow(history) == 0) stop("No valid candidates to build a pool.")
  
  history <- history[order(history[[cv_col]], decreasing = TRUE), , drop = FALSE]
  best_cv <- max(history[[cv_col]], na.rm = TRUE)
  
  pre_top_pool <- head(history, min(pre_top_n, nrow(history)))
  pre_delta_pool <- history[history[[cv_col]] >= (best_cv - pre_delta_cv), , drop = FALSE]
  pre_pool <- unique(rbind(pre_top_pool, pre_delta_pool))
  
  pre_pool <- pre_pool[order(pre_pool[[cv_col]], decreasing = TRUE), , drop = FALSE]
  final_top_pool <- head(pre_pool, min(final_top_n, nrow(pre_pool)))
  final_delta_pool <- pre_pool[pre_pool[[cv_col]] >= (best_cv - final_delta_cv), , drop = FALSE]
  final_pool <- unique(rbind(final_top_pool, final_delta_pool))
  
  pick_best <- function(df) {
    if (nrow(df) == 0) return(NULL)
    
    if (sd_col %in% names(df)) {
      df <- df[
        order(
          df[[cv_col]],
          df[[gap_col]],
          df[[sd_col]],
          decreasing = c(TRUE, FALSE, FALSE)
        ),
        ,
        drop = FALSE
      ]
    } else {
      df <- df[
        order(df[[cv_col]], df[[gap_col]], decreasing = c(TRUE, FALSE)),
        ,
        drop = FALSE
      ]
    }
    
    df[1, , drop = FALSE]
  }
  
  pool_gap15 <- final_pool[final_pool[[gap_col]] < CFG$gap_good, , drop = FALSE]
  pool_gap25 <- final_pool[final_pool[[gap_col]] <= CFG$gap_ok, , drop = FALSE]
  
  chosen <- pick_best(pool_gap15)
  level <- paste0("Gap < ", CFG$gap_good)
  
  if (is.null(chosen)) {
    chosen <- pick_best(pool_gap25)
    level <- paste0("Gap <= ", CFG$gap_ok)
  }
  
  if (is.null(chosen)) {
    chosen <- pick_best(final_pool)
    level <- "Final high-CV pool only"
  }
  
  list(
    pre_pool = pre_pool,
    final_pool = final_pool,
    chosen = chosen,
    level = level
  )
}

get_caret_summary <- function(model, model_name, X_data, y_data) {
  train_pred <- predict(model, X_data)
  
  train_r2 <- calc_r2(train_pred, y_data)
  train_rmse <- calc_rmse(train_pred, y_data)
  train_mae <- calc_mae(train_pred, y_data)
  train_nrmse <- calc_nrmse(train_pred, y_data)
  
  cv_r2 <- mean(model$resample$Rsquared, na.rm = TRUE)
  cv_rmse <- mean(model$resample$RMSE, na.rm = TRUE)
  cv_mae <- mean(model$resample$MAE, na.rm = TRUE)
  cv_nrmse <- cv_rmse / diff(range(y_data, na.rm = TRUE))
  
  gap_r2 <- train_r2 - cv_r2
  
  data.frame(
    Model = model_name,
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse,
    CV_nRMSE = cv_nrmse,
    CV_MAE = cv_mae,
    Train_R2 = train_r2,
    Train_RMSE = train_rmse,
    Train_nRMSE = train_nrmse,
    Train_MAE = train_mae,
    CV_R2_SD = sd(model$resample$Rsquared, na.rm = TRUE),
    CV_RMSE_SD = sd(model$resample$RMSE, na.rm = TRUE),
    CV_MAE_SD = sd(model$resample$MAE, na.rm = TRUE),
    Gap_R2 = gap_r2,
    Status = gap_status(gap_r2),
    Qualification = qualification_flag(gap_r2, cv_r2),
    stringsAsFactors = FALSE
  )
}


# 6. Elastic Net ----------------------------------------------------------------

enet_history <- data.frame()

enet_fit <- function(alpha, lambda) {
  set.seed(123)
  
  alpha <- max(0, min(alpha, 1))
  lambda <- max(0.00005, min(lambda, 0.12))
  
  out <- tryCatch({
    model <- caret::train(
      x = X_mat,
      y = y,
      method = "glmnet",
      trControl = ctrl,
      preProcess = c("center", "scale"),
      tuneGrid = data.frame(alpha = alpha, lambda = lambda),
      metric = "Rsquared"
    )
    
    cv_r2 <- mean(model$resample$Rsquared, na.rm = TRUE)
    cv_rmse <- mean(model$resample$RMSE, na.rm = TRUE)
    cv_mae <- mean(model$resample$MAE, na.rm = TRUE)
    cv_r2_sd <- sd(model$resample$Rsquared, na.rm = TRUE)
    
    train_pred <- predict(model, X_mat)
    train_r2 <- calc_r2(train_pred, y)
    gap_r2 <- train_r2 - cv_r2
    
    score <- generalization_score(
      cv_r2 = cv_r2,
      gap_r2 = gap_r2,
      cv_r2_sd = cv_r2_sd,
      gap_penalty = CFG$enet_penalty_gap,
      sd_penalty = CFG$enet_penalty_sd,
      gap_cutoff = 0.20,
      sd_cutoff = CFG$enet_sd_cutoff
    )
    
    enet_history <<- append_history(enet_history, list(
      alpha = alpha,
      lambda = lambda,
      CV_RMSE = cv_rmse,
      CV_R2 = cv_r2,
      CV_MAE = cv_mae,
      CV_R2_SD = cv_r2_sd,
      Train_R2 = train_r2,
      Gap_R2 = gap_r2,
      Score = score,
      Error = ""
    ))
    
    list(Score = score, Pred = 0)
    
  }, error = function(e) {
    enet_history <<- append_history(enet_history, list(
      alpha = alpha,
      lambda = lambda,
      CV_RMSE = NA,
      CV_R2 = NA,
      CV_MAE = NA,
      CV_R2_SD = NA,
      Train_R2 = NA,
      Gap_R2 = NA,
      Score = -999,
      Error = conditionMessage(e)
    ))
    
    list(Score = -999, Pred = 0)
  })
  
  out
}

set.seed(999)

enet_search <- tryCatch({
  BayesianOptimization(
    FUN = enet_fit,
    bounds = list(
      alpha = c(0, 1),
      lambda = c(0.00005, 0.12)
    ),
    init_points = 10,
    n_iter = 20,
    acq = "ucb",
    kappa = 2.2,
    eps = 1e-4,
    verbose = TRUE
  )
}, error = function(e) {
  message("ElasticNet BayesianOptimization stopped early: ", conditionMessage(e))
  NULL
})

valid_enet_history <- enet_history[
  is.finite(enet_history$CV_R2) & is.finite(enet_history$Gap_R2),
  ,
  drop = FALSE
]

if (nrow(valid_enet_history) == 0) {
  stop("No valid ElasticNet optimization result was found.")
}

valid_enet_history <- valid_enet_history[
  order(valid_enet_history$CV_R2, decreasing = TRUE),
  ,
  drop = FALSE
]

enet_pool_obj <- build_candidate_pool(
  history = valid_enet_history,
  cv_col = "CV_R2",
  gap_col = "Gap_R2",
  sd_col = "CV_R2_SD"
)

best_enet <- enet_pool_obj$chosen

model_enet <- caret::train(
  x = X_mat,
  y = y,
  method = "glmnet",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneGrid = data.frame(
    alpha = best_enet$alpha,
    lambda = best_enet$lambda
  ),
  metric = "Rsquared"
)

write.csv(enet_history, OUT_FILE("ElasticNet_Search_History"), row.names = FALSE)


# 7. Support Vector Machine -----------------------------------------------------

svm_history <- data.frame()

svm_fit <- function(sigma, C) {
  set.seed(123)
  
  sigma <- max(0.0003, min(sigma, 0.03))
  C <- max(0.1, min(C, 8))
  
  out <- tryCatch({
    model <- caret::train(
      x = X_mat,
      y = y,
      method = "svmRadial",
      trControl = ctrl,
      preProcess = c("center", "scale"),
      tuneGrid = data.frame(sigma = sigma, C = C),
      metric = "Rsquared"
    )
    
    cv_r2 <- mean(model$resample$Rsquared, na.rm = TRUE)
    cv_rmse <- mean(model$resample$RMSE, na.rm = TRUE)
    cv_mae <- mean(model$resample$MAE, na.rm = TRUE)
    cv_r2_sd <- sd(model$resample$Rsquared, na.rm = TRUE)
    
    train_pred <- predict(model, X_mat)
    train_r2 <- calc_r2(train_pred, y)
    gap_r2 <- train_r2 - cv_r2
    
    score <- generalization_score(
      cv_r2 = cv_r2,
      gap_r2 = gap_r2,
      cv_r2_sd = cv_r2_sd,
      gap_penalty = CFG$svm_penalty_gap,
      sd_penalty = CFG$svm_penalty_sd,
      gap_cutoff = 0.20,
      sd_cutoff = CFG$svm_sd_cutoff
    )
    
    svm_history <<- append_history(svm_history, list(
      sigma = sigma,
      C = C,
      CV_RMSE = cv_rmse,
      CV_R2 = cv_r2,
      CV_MAE = cv_mae,
      CV_R2_SD = cv_r2_sd,
      Train_R2 = train_r2,
      Gap_R2 = gap_r2,
      Score = score,
      Error = ""
    ))
    
    list(Score = score, Pred = 0)
    
  }, error = function(e) {
    svm_history <<- append_history(svm_history, list(
      sigma = sigma,
      C = C,
      CV_RMSE = NA,
      CV_R2 = NA,
      CV_MAE = NA,
      CV_R2_SD = NA,
      Train_R2 = NA,
      Gap_R2 = NA,
      Score = -999,
      Error = conditionMessage(e)
    ))
    
    list(Score = -999, Pred = 0)
  })
  
  out
}

set.seed(999)

svm_search <- tryCatch({
  BayesianOptimization(
    FUN = svm_fit,
    bounds = list(
      sigma = c(0.0003, 0.03),
      C = c(0.1, 8)
    ),
    init_points = 10,
    n_iter = 20,
    acq = "ucb",
    kappa = 2.2,
    eps = 0.0005,
    verbose = TRUE
  )
}, error = function(e) {
  message("SVM BayesianOptimization stopped early: ", conditionMessage(e))
  NULL
})

valid_svm_history <- svm_history[
  is.finite(svm_history$CV_R2) & is.finite(svm_history$Gap_R2),
  ,
  drop = FALSE
]

if (nrow(valid_svm_history) == 0) {
  stop("No valid SVM optimization result was found.")
}

valid_svm_history <- valid_svm_history[
  order(valid_svm_history$CV_R2, decreasing = TRUE),
  ,
  drop = FALSE
]

svm_pool_obj <- build_candidate_pool(
  history = valid_svm_history,
  cv_col = "CV_R2",
  gap_col = "Gap_R2",
  sd_col = "CV_R2_SD"
)

best_svm <- svm_pool_obj$chosen

model_svm <- caret::train(
  x = X_mat,
  y = y,
  method = "svmRadial",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneGrid = data.frame(
    sigma = best_svm$sigma,
    C = best_svm$C
  ),
  metric = "Rsquared"
)

write.csv(svm_history, OUT_FILE("SVM_Search_History"), row.names = FALSE)


# 8. Random Forest --------------------------------------------------------------

rf_history <- data.frame()

rf_fit <- function(mtry, min_node_size) {
  set.seed(123)
  
  mtry <- as.integer(round(mtry))
  min_node_size <- as.integer(round(min_node_size))
  
  mtry <- max(3, min(mtry, 7))
  min_node_size <- max(3, min(min_node_size, 14))
  
  out <- tryCatch({
    model <- caret::train(
      x = X_mat,
      y = y,
      method = "ranger",
      trControl = ctrl,
      num.trees = 1000,
      importance = "none",
      tuneGrid = data.frame(
        mtry = mtry,
        min.node.size = min_node_size,
        splitrule = "variance"
      ),
      metric = "Rsquared"
    )
    
    cv_r2 <- mean(model$resample$Rsquared, na.rm = TRUE)
    cv_rmse <- mean(model$resample$RMSE, na.rm = TRUE)
    cv_mae <- mean(model$resample$MAE, na.rm = TRUE)
    cv_r2_sd <- sd(model$resample$Rsquared, na.rm = TRUE)
    
    train_pred <- predict(model, X_mat)
    train_r2 <- calc_r2(train_pred, y)
    gap_r2 <- train_r2 - cv_r2
    
    score <- generalization_score(
      cv_r2 = cv_r2,
      gap_r2 = gap_r2,
      cv_r2_sd = cv_r2_sd,
      gap_penalty = CFG$rf_penalty_gap,
      sd_penalty = CFG$rf_penalty_sd,
      gap_cutoff = 0.20,
      sd_cutoff = CFG$rf_sd_cutoff
    )
    
    rf_history <<- append_history(rf_history, list(
      mtry = mtry,
      min.node.size = min_node_size,
      CV_RMSE = cv_rmse,
      CV_R2 = cv_r2,
      CV_MAE = cv_mae,
      CV_R2_SD = cv_r2_sd,
      Train_R2 = train_r2,
      Gap_R2 = gap_r2,
      Score = score,
      Error = ""
    ))
    
    list(Score = score, Pred = 0)
    
  }, error = function(e) {
    rf_history <<- append_history(rf_history, list(
      mtry = mtry,
      min.node.size = min_node_size,
      CV_RMSE = NA,
      CV_R2 = NA,
      CV_MAE = NA,
      CV_R2_SD = NA,
      Train_R2 = NA,
      Gap_R2 = NA,
      Score = -999,
      Error = conditionMessage(e)
    ))
    
    list(Score = -999, Pred = 0)
  })
  
  out
}

set.seed(999)

rf_search <- tryCatch({
  BayesianOptimization(
    FUN = rf_fit,
    bounds = list(
      mtry = c(3, 7),
      min_node_size = c(3.5, 14)
    ),
    init_points = 10,
    n_iter = 20,
    acq = "ucb",
    kappa = 2.2,
    eps = 0.01,
    verbose = TRUE
  )
}, error = function(e) {
  message("RF BayesianOptimization stopped early: ", conditionMessage(e))
  NULL
})

valid_rf_history <- rf_history[
  is.finite(rf_history$CV_R2) & is.finite(rf_history$Gap_R2),
  ,
  drop = FALSE
]

if (nrow(valid_rf_history) == 0) {
  stop("No valid RF optimization result was found.")
}

valid_rf_history <- dedup_integer_history(
  valid_rf_history,
  int_cols = c("mtry", "min.node.size")
)

rf_pool_obj <- build_candidate_pool(
  history = valid_rf_history,
  cv_col = "CV_R2",
  gap_col = "Gap_R2",
  sd_col = "CV_R2_SD",
  pre_top_n = 18,
  final_top_n = 10,
  pre_delta_cv = 0.012,
  final_delta_cv = 0.008
)

best_rf <- rf_pool_obj$chosen
best_rf_mtry <- as.integer(round(best_rf$mtry))
best_rf_min_node_size <- as.integer(round(best_rf$min.node.size))

model_rf <- caret::train(
  x = X_mat,
  y = y,
  method = "ranger",
  trControl = ctrl,
  num.trees = 2200,
  importance = "permutation",
  tuneGrid = data.frame(
    mtry = best_rf_mtry,
    min.node.size = best_rf_min_node_size,
    splitrule = "variance"
  ),
  metric = "Rsquared"
)

write.csv(rf_history, OUT_FILE("RF_Search_History"), row.names = FALSE)

rf_importance <- ranger::importance(model_rf$finalModel)

rf_importance_df <- data.frame(
  Feature = names(rf_importance),
  Importance = as.numeric(rf_importance),
  stringsAsFactors = FALSE
)

rf_importance_df <- rf_importance_df[
  order(rf_importance_df$Importance, decreasing = TRUE),
]

write.csv(
  rf_importance_df,
  OUT_FILE("RF_Permutation_Importance"),
  row.names = FALSE
)


# 9. XGBoost --------------------------------------------------------------------

cat("\n>>> XGBoost optimization: Bayesian search + deterministic final evaluation.\n")
cat(">>> Search space: max_depth = 4-5, min_child_weight = 2-4.\n")

try(parallel::stopCluster(cl), silent = TRUE)
foreach::registerDoSEQ()

xgb_history <- data.frame()

xgb_generalization_score <- function(cv_r2, gap_r2, cv_r2_sd) {
  score <- cv_r2
  
  score <- score - 0.45 * pmax(0, gap_r2 - 0.198)
  score <- score - 1.25 * pmax(0, gap_r2 - 0.200)
  score <- score - 0.04 * pmax(0, cv_r2_sd - 0.15)
  
  score <- score + ifelse(gap_r2 < 0.200, 0.004, 0)
  score <- score + ifelse(gap_r2 >= 0.180 & gap_r2 < 0.200, 0.002, 0)
  
  bad <- !is.finite(cv_r2) | !is.finite(gap_r2) | !is.finite(cv_r2_sd)
  score[bad] <- -999
  
  score
}

make_xgb_params <- function(row_df, nthread_use = 1L, seed_use = 2026L) {
  max_depth_use <- if ("max_depth" %in% names(row_df)) {
    as.integer(round(row_df$max_depth[1]))
  } else {
    5L
  }
  
  min_child_weight_use <- if ("min_child_weight" %in% names(row_df)) {
    as.integer(round(row_df$min_child_weight[1]))
  } else {
    2L
  }
  
  list(
    objective = "reg:squarederror",
    booster = "gbtree",
    max_depth = max_depth_use,
    eta = row_df$eta[1],
    gamma = row_df$gamma[1],
    colsample_bytree = row_df$colsample_bytree[1],
    subsample = row_df$subsample[1],
    min_child_weight = min_child_weight_use,
    lambda = row_df$lambda[1],
    alpha = row_df$alpha[1],
    seed = seed_use,
    seed_per_iteration = FALSE,
    nthread = nthread_use
  )
}

evaluate_xgb_combo <- function(params_row, X_mat, y, cv_index, seed_offset = 0) {
  resample_metrics <- data.frame()
  best_rounds_vec <- c()
  
  params <- make_xgb_params(
    row_df = params_row,
    nthread_use = 1L,
    seed_use = 2026L + seed_offset
  )
  
  for (i in seq_along(cv_index)) {
    nm <- names(cv_index)[i]
    
    set.seed(2026L + seed_offset + i)
    
    train_idx <- cv_index[[i]]
    test_idx <- setdiff(seq_along(y), train_idx)
    
    dtrain_outer <- xgb.DMatrix(
      data = X_mat[train_idx, , drop = FALSE],
      label = y[train_idx]
    )
    
    dtest_outer <- xgb.DMatrix(
      data = X_mat[test_idx, , drop = FALSE],
      label = y[test_idx]
    )
    
    inner_valid_size <- max(5, floor(length(train_idx) * 0.20))
    inner_valid_local <- sample(seq_along(train_idx), size = inner_valid_size)
    inner_valid_idx <- train_idx[inner_valid_local]
    inner_train_idx <- setdiff(train_idx, inner_valid_idx)
    
    dtrain_inner <- xgb.DMatrix(
      data = X_mat[inner_train_idx, , drop = FALSE],
      label = y[inner_train_idx]
    )
    
    dvalid_inner <- xgb.DMatrix(
      data = X_mat[inner_valid_idx, , drop = FALSE],
      label = y[inner_valid_idx]
    )
    
    fold_params <- params
    fold_params$seed <- 2026L + seed_offset + i
    
    model_inner <- tryCatch({
      xgb.train(
        params = fold_params,
        data = dtrain_inner,
        nrounds = 5000,
        watchlist = list(train = dtrain_inner, valid = dvalid_inner),
        early_stopping_rounds = 200,
        maximize = FALSE,
        verbose = 0
      )
    }, error = function(e) NULL)
    
    if (is.null(model_inner)) return(NULL)
    
    best_nrounds <- model_inner$best_iteration
    
    if (is.null(best_nrounds) || is.na(best_nrounds) || best_nrounds <= 0) {
      best_nrounds <- 200L
    }
    
    best_rounds_vec <- c(best_rounds_vec, best_nrounds)
    
    model_outer <- tryCatch({
      xgb.train(
        params = fold_params,
        data = dtrain_outer,
        nrounds = best_nrounds,
        verbose = 0
      )
    }, error = function(e) NULL)
    
    if (is.null(model_outer)) return(NULL)
    
    pred_test <- predict(model_outer, dtest_outer)
    obs_test <- y[test_idx]
    
    resample_metrics <- rbind(
      resample_metrics,
      data.frame(
        Resample = nm,
        RMSE = calc_rmse(pred_test, obs_test),
        Rsquared = calc_r2(pred_test, obs_test),
        MAE = calc_mae(pred_test, obs_test),
        Best_nrounds = best_nrounds,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(resample_metrics) == 0) return(NULL)
  
  cv_r2 <- mean(resample_metrics$Rsquared, na.rm = TRUE)
  cv_rmse <- mean(resample_metrics$RMSE, na.rm = TRUE)
  cv_mae <- mean(resample_metrics$MAE, na.rm = TRUE)
  cv_r2_sd <- sd(resample_metrics$Rsquared, na.rm = TRUE)
  
  final_nrounds <- as.integer(round(median(best_rounds_vec, na.rm = TRUE)))
  
  if (!is.finite(final_nrounds) || final_nrounds <= 0) {
    final_nrounds <- 200L
  }
  
  dall <- xgb.DMatrix(data = X_mat, label = y)
  
  full_params <- params
  full_params$seed <- 2026L + seed_offset
  
  model_full <- tryCatch({
    xgb.train(
      params = full_params,
      data = dall,
      nrounds = final_nrounds,
      verbose = 0
    )
  }, error = function(e) NULL)
  
  if (is.null(model_full)) return(NULL)
  
  pred_full <- predict(model_full, dall)
  train_r2 <- calc_r2(pred_full, y)
  gap_r2 <- train_r2 - cv_r2
  
  score <- xgb_generalization_score(
    cv_r2 = cv_r2,
    gap_r2 = gap_r2,
    cv_r2_sd = cv_r2_sd
  )
  
  data.frame(
    max_depth = params$max_depth,
    min_child_weight = params$min_child_weight,
    gamma = params$gamma,
    eta = params$eta,
    colsample_bytree = params$colsample_bytree,
    subsample = params$subsample,
    lambda = params$lambda,
    alpha = params$alpha,
    Mean_nrounds = mean(best_rounds_vec, na.rm = TRUE),
    Median_nrounds = final_nrounds,
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse,
    CV_nRMSE = cv_rmse / diff(range(y, na.rm = TRUE)),
    CV_MAE = cv_mae,
    Train_R2 = train_r2,
    Gap_R2 = gap_r2,
    CV_R2_SD = cv_r2_sd,
    Score = as.numeric(score[1]),
    Status = gap_status(gap_r2),
    stringsAsFactors = FALSE
  )
}

xgb_bo_objective <- function(max_depth,
                             min_child_weight,
                             gamma,
                             eta,
                             colsample_bytree,
                             subsample,
                             lambda,
                             alpha) {
  p <- data.frame(
    max_depth = as.integer(round(max_depth)),
    min_child_weight = as.integer(round(min_child_weight)),
    gamma = gamma,
    eta = eta,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    lambda = lambda,
    alpha = alpha
  )
  
  p$max_depth <- max(4L, min(p$max_depth, 5L))
  p$min_child_weight <- max(2L, min(p$min_child_weight, 4L))
  p$gamma <- max(0.15, min(p$gamma, 0.55))
  p$eta <- max(0.020, min(p$eta, 0.040))
  p$colsample_bytree <- max(0.75, min(p$colsample_bytree, 0.92))
  p$subsample <- max(0.82, min(p$subsample, 0.96))
  p$lambda <- max(0.05, min(p$lambda, 0.80))
  p$alpha <- max(0.00, min(p$alpha, 0.20))
  
  out <- tryCatch({
    res <- evaluate_xgb_combo(
      params_row = p,
      X_mat = X_mat,
      y = y,
      cv_index = cv_index,
      seed_offset = 0
    )
    
    if (is.null(res)) stop("Evaluation returned NULL.")
    
    xgb_history <<- append_history(xgb_history, as.list(res[1, ]))
    
    cat(sprintf(
      "XGB | depth=%d mcw=%d gamma=%.4f eta=%.4f col=%.4f sub=%.4f lambda=%.4f alpha=%.4f | CV=%.4f Train=%.4f Gap=%.4f SD=%.4f Score=%.4f MedianRounds=%d\n",
      res$max_depth,
      res$min_child_weight,
      res$gamma,
      res$eta,
      res$colsample_bytree,
      res$subsample,
      res$lambda,
      res$alpha,
      res$CV_R2,
      res$Train_R2,
      res$Gap_R2,
      res$CV_R2_SD,
      res$Score,
      res$Median_nrounds
    ))
    
    list(Score = as.numeric(res$Score[1]), Pred = 0)
    
  }, error = function(e) {
    xgb_history <<- append_history(xgb_history, list(
      max_depth = p$max_depth,
      min_child_weight = p$min_child_weight,
      gamma = p$gamma,
      eta = p$eta,
      colsample_bytree = p$colsample_bytree,
      subsample = p$subsample,
      lambda = p$lambda,
      alpha = p$alpha,
      Mean_nrounds = NA,
      Median_nrounds = NA,
      CV_R2 = NA,
      CV_RMSE = NA,
      CV_nRMSE = NA,
      CV_MAE = NA,
      Train_R2 = NA,
      Gap_R2 = NA,
      CV_R2_SD = NA,
      Score = -999,
      Status = "Error"
    ))
    
    list(Score = -999, Pred = 0)
  })
  
  out
}

set.seed(2026)

xgb_search <- tryCatch({
  BayesianOptimization(
    FUN = xgb_bo_objective,
    bounds = list(
      max_depth = c(4L, 5L),
      min_child_weight = c(2L, 4L),
      gamma = c(0.15, 0.55),
      eta = c(0.020, 0.040),
      colsample_bytree = c(0.75, 0.92),
      subsample = c(0.82, 0.96),
      lambda = c(0.05, 0.80),
      alpha = c(0.00, 0.20)
    ),
    init_points = 30,
    n_iter = 60,
    acq = "ucb",
    kappa = 1.9,
    eps = 0.001,
    verbose = TRUE
  )
}, error = function(e) {
  message("XGB BayesianOptimization stopped early: ", conditionMessage(e))
  message("Proceeding with valid parameter sets already recorded in xgb_history.")
  NULL
})

write.csv(xgb_history, OUT_FILE("XGB_Search_History"), row.names = FALSE)

valid_xgb_history <- xgb_history %>%
  dplyr::filter(
    is.finite(CV_R2),
    is.finite(Train_R2),
    is.finite(Gap_R2),
    is.finite(CV_R2_SD),
    is.finite(Score)
  )

if (nrow(valid_xgb_history) == 0) {
  stop("No valid XGB optimization result was found.")
}

candidate_pool <- valid_xgb_history %>%
  dplyr::filter(Gap_R2 <= CFG$gap_ok) %>%
  dplyr::arrange(desc(Score), desc(CV_R2), Gap_R2, CV_R2_SD)

if (nrow(candidate_pool) < 12) {
  candidate_pool <- valid_xgb_history %>%
    dplyr::arrange(desc(Score), desc(CV_R2), Gap_R2, CV_R2_SD)
}

top12_xgb <- head(candidate_pool, min(12, nrow(candidate_pool)))

write.csv(top12_xgb, OUT_FILE("XGB_Top12Candidates"), row.names = FALSE)

verify_one_xgb <- function(row_df, seed_offset = 0) {
  evaluate_xgb_combo(
    params_row = row_df[, c(
      "max_depth",
      "min_child_weight",
      "gamma",
      "eta",
      "colsample_bytree",
      "subsample",
      "lambda",
      "alpha"
    )],
    X_mat = X_mat,
    y = y,
    cv_index = cv_index,
    seed_offset = seed_offset
  )
}

verified_xgb <- data.frame()

for (i in seq_len(nrow(top12_xgb))) {
  res_i <- verify_one_xgb(top12_xgb[i, , drop = FALSE], seed_offset = 0)
  
  if (!is.null(res_i)) {
    res_i$Candidate_ID <- paste0("Top", i)
    verified_xgb <- rbind(verified_xgb, res_i)
  }
}

if (nrow(verified_xgb) == 0) {
  stop("No valid verified XGB candidate was found.")
}

write.csv(verified_xgb, OUT_FILE("XGB_VerifiedTop12"), row.names = FALSE)

verified_xgb_ranked <- verified_xgb %>%
  dplyr::arrange(desc(Score), desc(CV_R2), Gap_R2, CV_R2_SD)

top_for_rounds <- head(verified_xgb_ranked, min(8, nrow(verified_xgb_ranked)))

write.csv(
  top_for_rounds,
  OUT_FILE("XGB_TopCandidates_For_Nrounds"),
  row.names = FALSE
)

evaluate_xgb_fixed_nrounds <- function(row_df,
                                       nrounds_fixed,
                                       X_mat,
                                       y,
                                       cv_index,
                                       seed_base = 2026L,
                                       return_model = TRUE) {
  resample_metrics <- data.frame()
  xgb_oof_pred <- data.frame()
  
  base_params <- make_xgb_params(
    row_df = row_df,
    nthread_use = 1L,
    seed_use = seed_base
  )
  
  for (i in seq_along(cv_index)) {
    set.seed(seed_base + i)
    
    train_idx <- cv_index[[i]]
    test_idx <- setdiff(seq_along(y), train_idx)
    
    dtrain_i <- xgb.DMatrix(
      data = X_mat[train_idx, , drop = FALSE],
      label = y[train_idx]
    )
    
    dtest_i <- xgb.DMatrix(
      data = X_mat[test_idx, , drop = FALSE],
      label = y[test_idx]
    )
    
    fold_params <- base_params
    fold_params$seed <- seed_base + i
    
    model_i <- tryCatch({
      xgb.train(
        params = fold_params,
        data = dtrain_i,
        nrounds = nrounds_fixed,
        verbose = 0
      )
    }, error = function(e) NULL)
    
    if (is.null(model_i)) return(NULL)
    
    pred_i <- predict(model_i, dtest_i)
    obs_i <- y[test_idx]
    
    resample_metrics <- rbind(
      resample_metrics,
      data.frame(
        Resample = names(cv_index)[i],
        RMSE = calc_rmse(pred_i, obs_i),
        Rsquared = calc_r2(pred_i, obs_i),
        MAE = calc_mae(pred_i, obs_i),
        stringsAsFactors = FALSE
      )
    )
    
    xgb_oof_pred <- rbind(
      xgb_oof_pred,
      data.frame(
        rowIndex = test_idx,
        obs = obs_i,
        pred = pred_i,
        Resample = names(cv_index)[i],
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(resample_metrics) == 0) return(NULL)
  
  set.seed(seed_base)
  
  dall <- xgb.DMatrix(data = X_mat, label = y)
  
  final_params <- base_params
  final_params$seed <- seed_base
  
  model_full <- tryCatch({
    xgb.train(
      params = final_params,
      data = dall,
      nrounds = nrounds_fixed,
      verbose = 0
    )
  }, error = function(e) NULL)
  
  if (is.null(model_full)) return(NULL)
  
  xgb_train_pred <- predict(model_full, dall)
  
  cv_r2 <- mean(resample_metrics$Rsquared, na.rm = TRUE)
  cv_rmse <- mean(resample_metrics$RMSE, na.rm = TRUE)
  cv_mae <- mean(resample_metrics$MAE, na.rm = TRUE)
  train_r2 <- calc_r2(xgb_train_pred, y)
  gap_r2 <- train_r2 - cv_r2
  cv_r2_sd <- sd(resample_metrics$Rsquared, na.rm = TRUE)
  
  xgb_metrics <- data.frame(
    Model = "XGB_Bayes",
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse,
    CV_nRMSE = cv_rmse / diff(range(y, na.rm = TRUE)),
    CV_MAE = cv_mae,
    Train_R2 = train_r2,
    Train_RMSE = calc_rmse(xgb_train_pred, y),
    Train_nRMSE = calc_nrmse(xgb_train_pred, y),
    Train_MAE = calc_mae(xgb_train_pred, y),
    CV_R2_SD = cv_r2_sd,
    CV_RMSE_SD = sd(resample_metrics$RMSE, na.rm = TRUE),
    CV_MAE_SD = sd(resample_metrics$MAE, na.rm = TRUE),
    Gap_R2 = gap_r2,
    Status = gap_status(gap_r2),
    Qualification = qualification_flag(gap_r2, cv_r2),
    stringsAsFactors = FALSE
  )
  
  summary_row <- data.frame(
    max_depth = base_params$max_depth,
    min_child_weight = base_params$min_child_weight,
    gamma = base_params$gamma,
    eta = base_params$eta,
    colsample_bytree = base_params$colsample_bytree,
    subsample = base_params$subsample,
    lambda = base_params$lambda,
    alpha = base_params$alpha,
    Mean_nrounds = nrounds_fixed,
    Median_nrounds = nrounds_fixed,
    CV_R2 = xgb_metrics$CV_R2,
    CV_RMSE = xgb_metrics$CV_RMSE,
    CV_nRMSE = xgb_metrics$CV_nRMSE,
    CV_MAE = xgb_metrics$CV_MAE,
    Train_R2 = xgb_metrics$Train_R2,
    Gap_R2 = xgb_metrics$Gap_R2,
    CV_R2_SD = xgb_metrics$CV_R2_SD,
    Score = as.numeric(xgb_generalization_score(
      cv_r2 = xgb_metrics$CV_R2,
      gap_r2 = xgb_metrics$Gap_R2,
      cv_r2_sd = xgb_metrics$CV_R2_SD
    )),
    Status = xgb_metrics$Status,
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary_row,
    metrics = xgb_metrics,
    model = if (return_model) model_full else NULL,
    params = final_params,
    nrounds = nrounds_fixed,
    train_pred = xgb_train_pred,
    resample_metrics = resample_metrics,
    oof_pred = xgb_oof_pred
  )
}

round_grid <- sort(unique(c(
  240L, 260L, 280L, 300L, 320L, 340L,
  350L, 360L, 370L, 375L, 380L, 385L, 390L, 395L, 400L,
  405L, 410L, 415L, 420L, 430L, 440L, 460L, 480L
)))

round_tuning_results <- data.frame()

cat("\n>>> Tuning nrounds for verified XGBoost candidates...\n")

for (i in seq_len(nrow(top_for_rounds))) {
  for (nr in round_grid) {
    cat(sprintf(
      "Nrounds tuning | candidate %d/%d | nrounds = %d\n",
      i,
      nrow(top_for_rounds),
      nr
    ))
    
    tmp <- evaluate_xgb_fixed_nrounds(
      row_df = top_for_rounds[i, , drop = FALSE],
      nrounds_fixed = nr,
      X_mat = X_mat,
      y = y,
      cv_index = cv_index,
      seed_base = 2026L,
      return_model = FALSE
    )
    
    if (!is.null(tmp)) {
      tmp_row <- tmp$summary
      tmp_row$Candidate_ID <- top_for_rounds$Candidate_ID[i]
      round_tuning_results <- rbind(round_tuning_results, tmp_row)
    }
  }
}

if (nrow(round_tuning_results) == 0) {
  stop("No valid XGB nrounds tuning result was found.")
}

write.csv(
  round_tuning_results,
  OUT_FILE("XGB_NroundsTuning_All"),
  row.names = FALSE
)

best_cv_any <- max(round_tuning_results$CV_R2, na.rm = TRUE)

strict_pool <- round_tuning_results %>%
  dplyr::filter(
    is.finite(CV_R2),
    is.finite(Gap_R2),
    is.finite(CV_R2_SD),
    Gap_R2 < 0.200,
    CV_R2 >= best_cv_any - 0.020
  ) %>%
  dplyr::mutate(
    Gap_Target_Distance = abs(Gap_R2 - 0.1985),
    FinalScore = xgb_generalization_score(CV_R2, Gap_R2, CV_R2_SD)
  ) %>%
  dplyr::arrange(
    desc(CV_R2),
    Gap_Target_Distance,
    CV_R2_SD
  )

if (nrow(strict_pool) == 0) {
  strict_pool <- round_tuning_results %>%
    dplyr::filter(
      is.finite(CV_R2),
      is.finite(Gap_R2),
      is.finite(CV_R2_SD),
      Gap_R2 < 0.200
    ) %>%
    dplyr::mutate(
      Gap_Target_Distance = abs(Gap_R2 - 0.1985),
      FinalScore = xgb_generalization_score(CV_R2, Gap_R2, CV_R2_SD)
    ) %>%
    dplyr::arrange(
      desc(CV_R2),
      Gap_Target_Distance,
      CV_R2_SD
    )
}

if (nrow(strict_pool) == 0) {
  closest_xgb <- round_tuning_results %>%
    dplyr::filter(
      is.finite(CV_R2),
      is.finite(Gap_R2),
      is.finite(CV_R2_SD)
    ) %>%
    dplyr::mutate(
      Gap_Distance = abs(Gap_R2 - 0.200)
    ) %>%
    dplyr::arrange(
      Gap_Distance,
      desc(CV_R2),
      CV_R2_SD
    )
  
  closest_xgb <- head(closest_xgb, min(10, nrow(closest_xgb)))
  
  write.csv(
    closest_xgb,
    OUT_FILE("XGB_Closest_To_Gap020_NoStrictCandidate"),
    row.names = FALSE
  )
  
  stop(
    "No XGB candidate with Gap_R2 < 0.20 was found. ",
    "Closest candidates were saved to XGB_Closest_To_Gap020_NoStrictCandidate."
  )
}

best_xgb_verified <- strict_pool %>%
  dplyr::slice(1) %>%
  dplyr::select(-Gap_Target_Distance, -FinalScore)

write.csv(
  best_xgb_verified,
  OUT_FILE("XGB_FinalBest_StrictGap"),
  row.names = FALSE
)

cat("\n>>> Best strict-gap XGBoost candidate:\n")
print(best_xgb_verified, digits = 4)

cat("\n>>> Running final deterministic XGBoost evaluation...\n")

best_xgb_nrounds <- as.integer(round(best_xgb_verified$Median_nrounds))

if (!is.finite(best_xgb_nrounds) || best_xgb_nrounds <= 0) {
  stop("Invalid XGB nrounds.")
}

final_xgb_eval <- evaluate_xgb_fixed_nrounds(
  row_df = best_xgb_verified,
  nrounds_fixed = best_xgb_nrounds,
  X_mat = X_mat,
  y = y,
  cv_index = cv_index,
  seed_base = 2026L,
  return_model = TRUE
)

if (is.null(final_xgb_eval)) {
  stop("Final deterministic XGB evaluation failed.")
}

model_xgb_bayes <- final_xgb_eval$model
final_xgb_params <- final_xgb_eval$params
xgb_train_pred <- final_xgb_eval$train_pred
xgb_resample_metrics <- final_xgb_eval$resample_metrics
xgb_oof_pred <- final_xgb_eval$oof_pred
xgb_metrics <- final_xgb_eval$metrics

best_xgb_verified <- final_xgb_eval$summary

write.csv(
  best_xgb_verified,
  OUT_FILE("XGB_FinalBest_Deterministic"),
  row.names = FALSE
)

write.csv(
  xgb_resample_metrics,
  OUT_FILE("XGB_Resample_Metrics"),
  row.names = FALSE
)

write.csv(
  xgb_oof_pred,
  OUT_FILE("XGB_OOF_Predictions"),
  row.names = FALSE
)

write.csv(
  xgb_metrics,
  OUT_FILE("XGB_Final_Metrics_Deterministic"),
  row.names = FALSE
)

xgb.save(model_xgb_bayes, OUT_FILE("Final_XGB_Bayes", ext = "json"))

cat("\n>>> Deterministic final XGBoost metrics:\n")
print(xgb_metrics, digits = 4)

if (xgb_metrics$Gap_R2 >= 0.20) {
  stop("Final deterministic XGB Gap_R2 is still >= 0.20.")
}


# 10. Model comparison ----------------------------------------------------------

caret_models <- list(
  ElasticNet_Bayes = model_enet,
  SVM_Bayes = model_svm,
  RF_Bayes = model_rf
)

caret_metrics <- do.call(
  rbind,
  lapply(names(caret_models), function(n) {
    get_caret_summary(caret_models[[n]], n, X_mat, y)
  })
)

all_metrics <- rbind(caret_metrics, xgb_metrics)
all_metrics <- all_metrics[order(all_metrics$CV_R2, decreasing = TRUE), ]

cat("\n>>> Final model comparison based on main CV metrics:\n")
print(all_metrics, digits = 4)


# 11. Supplementary pooled CV R2 ------------------------------------------------

calc_true_r2 <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  pred <- pred[ok]
  obs <- obs[ok]
  
  if (length(obs) < 2) return(NA_real_)
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  
  1 - (ss_res / ss_tot)
}

si_true_cvr2 <- data.frame()

for (model_name in names(caret_models)) {
  final_model <- caret_models[[model_name]]
  pred_df <- final_model$pred
  
  tune_cols <- names(final_model$bestTune)
  
  for (nm in tune_cols) {
    pred_df <- pred_df[pred_df[[nm]] == final_model$bestTune[[nm]], , drop = FALSE]
  }
  
  pooled_pred <- pred_df %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      obs = mean(obs, na.rm = TRUE),
      pred = mean(pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  si_true_cvr2 <- rbind(
    si_true_cvr2,
    data.frame(
      Model = model_name,
      CV_R2_Traditional = calc_true_r2(pooled_pred$pred, pooled_pred$obs),
      stringsAsFactors = FALSE
    )
  )
}

xgb_pooled_pred <- xgb_oof_pred %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(
    obs = mean(obs, na.rm = TRUE),
    pred = mean(pred, na.rm = TRUE),
    .groups = "drop"
  )

si_true_cvr2 <- rbind(
  si_true_cvr2,
  data.frame(
    Model = "XGB_Bayes",
    CV_R2_Traditional = calc_true_r2(xgb_pooled_pred$pred, xgb_pooled_pred$obs),
    stringsAsFactors = FALSE
  )
)

all_metrics <- all_metrics %>%
  dplyr::left_join(si_true_cvr2, by = "Model") %>%
  dplyr::select(
    Model,
    CV_R2,
    CV_R2_Traditional,
    CV_RMSE,
    CV_nRMSE,
    CV_MAE,
    Train_R2,
    Train_RMSE,
    Train_nRMSE,
    Train_MAE,
    CV_R2_SD,
    CV_RMSE_SD,
    CV_MAE_SD,
    Gap_R2,
    Status,
    Qualification
  )

write.csv(
  all_metrics,
  OUT_FILE("Model_Performance_Summary"),
  row.names = FALSE
)

write.csv(
  si_true_cvr2,
  OUT_FILE("All_Models_True_CVR2_SI"),
  row.names = FALSE
)

cat("\n>>> Supplementary pooled CV R2:\n")
print(si_true_cvr2, digits = 4)

cat("\n>>> Final model comparison with pooled CV R2:\n")
print(all_metrics, digits = 4)


# 12. Final model selection and export -----------------------------------------

all_params <- dplyr::bind_rows(
  data.frame(Model = "ElasticNet_Bayes", model_enet$bestTune),
  data.frame(Model = "SVM_Bayes", model_svm$bestTune),
  data.frame(Model = "RF_Bayes", model_rf$bestTune),
  data.frame(
    Model = "XGB_Bayes",
    nrounds = best_xgb_nrounds,
    max_depth = final_xgb_params$max_depth,
    eta = final_xgb_params$eta,
    gamma = final_xgb_params$gamma,
    colsample_bytree = final_xgb_params$colsample_bytree,
    min_child_weight = final_xgb_params$min_child_weight,
    subsample = final_xgb_params$subsample,
    lambda = final_xgb_params$lambda,
    alpha = final_xgb_params$alpha
  )
)

write.csv(
  all_params,
  OUT_FILE("Model_Best_Hyperparameters"),
  row.names = FALSE
)

model_selection_table <- all_metrics %>%
  dplyr::mutate(
    Predictive_OK = CV_R2 >= CFG$low_cv_cutoff,
    Gap_Strict_OK = Gap_R2 < 0.20,
    Gap_Acceptable_OK = Gap_R2 <= CFG$gap_ok
  )

primary_pool <- model_selection_table %>%
  dplyr::filter(Predictive_OK, Gap_Strict_OK)

if (nrow(primary_pool) > 0) {
  best_model <- primary_pool %>%
    dplyr::arrange(
      desc(CV_R2),
      CV_RMSE,
      CV_MAE,
      Gap_R2,
      CV_R2_SD
    ) %>%
    dplyr::slice(1)
  
  selection_level <- "Predictive models with Gap_R2 < 0.20"
  
} else {
  secondary_pool <- model_selection_table %>%
    dplyr::filter(Predictive_OK, Gap_Acceptable_OK)
  
  if (nrow(secondary_pool) > 0) {
    best_model <- secondary_pool %>%
      dplyr::arrange(
        desc(CV_R2),
        CV_RMSE,
        CV_MAE,
        Gap_R2,
        CV_R2_SD
      ) %>%
      dplyr::slice(1)
    
    selection_level <- paste0("Predictive models with Gap_R2 <= ", CFG$gap_ok)
    
  } else {
    best_model <- model_selection_table %>%
      dplyr::arrange(
        desc(CV_R2),
        CV_RMSE,
        CV_MAE,
        Gap_R2,
        CV_R2_SD
      ) %>%
      dplyr::slice(1)
    
    selection_level <- "Fallback: highest CV_R2 among all models"
  }
}

best_model_name <- best_model$Model

cat("\nSelected model:", best_model_name, "\n")
cat("Selection level:", selection_level, "\n")

write.csv(
  model_selection_table,
  OUT_FILE("Model_Selection_Table"),
  row.names = FALSE
)

if (best_model_name == "XGB_Bayes") {
  saveRDS(
    list(
      model = model_xgb_bayes,
      params = final_xgb_params,
      nrounds = best_xgb_nrounds,
      metrics = xgb_metrics,
      resample_metrics = xgb_resample_metrics,
      oof_pred = xgb_oof_pred,
      train_pred = xgb_train_pred,
      best_xgb_verified = best_xgb_verified
    ),
    OUT_FILE("Selected_Model_XGB_Bayes", ext = "rds")
  )
} else {
  saveRDS(
    caret_models[[best_model_name]],
    OUT_FILE(paste0("Selected_Model_", best_model_name), ext = "rds")
  )
}

writeLines(
  capture.output(sessionInfo()),
  OUT_FILE("sessionInfo", ext = "txt")
)


# 13. Prediction export ---------------------------------------------------------

cat("\n>>> Exporting unified prediction file...\n")

if (best_model_name == "XGB_Bayes") {
  train_pred_df <- data.frame(
    rowIndex = seq_along(y),
    Observed = y,
    Predicted = xgb_train_pred,
    Set = "Training",
    stringsAsFactors = FALSE
  )
  
  cv_pred_df <- xgb_oof_pred %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      Observed = mean(obs, na.rm = TRUE),
      Predicted = mean(pred, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Set = "CrossValidation") %>%
    dplyr::select(rowIndex, Observed, Predicted, Set)
  
  final_pred_export <- dplyr::bind_rows(train_pred_df, cv_pred_df)
  
} else {
  final_model <- caret_models[[best_model_name]]
  train_pred <- predict(final_model, X_mat)
  
  train_pred_df <- data.frame(
    rowIndex = seq_along(y),
    Observed = y,
    Predicted = train_pred,
    Set = "Training",
    stringsAsFactors = FALSE
  )
  
  pred_df <- final_model$pred
  tune_cols <- names(final_model$bestTune)
  
  for (nm in tune_cols) {
    pred_df <- pred_df[pred_df[[nm]] == final_model$bestTune[[nm]], , drop = FALSE]
  }
  
  cv_pred_df <- pred_df %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      Observed = mean(obs, na.rm = TRUE),
      Predicted = mean(pred, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Set = "CrossValidation") %>%
    dplyr::select(rowIndex, Observed, Predicted, Set)
  
  final_pred_export <- dplyr::bind_rows(train_pred_df, cv_pred_df)
}

write.csv(
  final_pred_export,
  OUT_FILE("Best_Model_Predictions"),
  row.names = FALSE
)

try(parallel::stopCluster(cl), silent = TRUE)
foreach::registerDoSEQ()

cat("Saved:", OUT_FILE("Best_Model_Predictions"), "\n")
cat("Done.\n")





