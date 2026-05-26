# ==============================================================================
# 03_streamlining_m4_m8_m12.R
# ------------------------------------------------------------------------------
# Streamlined model comparison and final 4-variable model bundle export.
# Computational settings are kept unchanged:
#   - M4_Core and M8_Extended feature sets unchanged
#   - 10 x 10 repeated cross-validation
#   - XGB and RF rough-grid/local Bayesian optimization settings unchanged
#   - Final RF/XGB tuning settings unchanged
# ==============================================================================

# ==============================================================================
# Streamlining comparison
# ------------------------------------------------------------------------------
# Purpose:
#   Train only two streamlined candidate models:
#     M4_Core:
#       Solubility + Size + log.Kow + Feeding
#
#     M8_Extended:
#       Solubility + log.Kow + Complexity + Respiration + Locomotion
#       + Feeding + Size + HLC
#
#   Read M12_Full performance directly from the previous full-model training
#   summary, then combine M4 / M8 / M12 into one SI-ready comparison table.
#
# Design:
#   1) Use 10 x 10 repeated CV for M4_Core and M8_Extended
#   2) XGB: rough grid -> local Bayesian optimization -> nrounds grid
#   3) RF : rough grid -> local Bayesian optimization -> final evaluation
#   4) Read M12_Full from Model_Performance_Summary_FINAL.csv
#   5) Save final M4_Core XGB and RF prediction bundles
#   6) Add prediction helper functions and an example prediction script
# ==============================================================================


# 0. Global settings ------------------------------------------------------------

# Output naming is simplified by saving all outputs into one run folder.
# File names are kept short and publication/GitHub friendly.

output_dir <- "outputs_streamlining_m4_m8_m12_10x10"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

OUT_FILE <- function(name, ext = "csv") {
  file.path(output_dir, paste0(name, ".", ext))
}

# Full-precision M12 full-model summary generated from the main training script.
M12_METRICS_FILE <- "Model_Performance_Summary_FINAL.csv"

CFG <- list(
  cv_folds = 10,
  cv_repeats = 10,
  
  gap_good = 0.15,
  gap_ok = 0.25,
  low_cv_cutoff = 0.60,
  
  xgb_penalty_gap = 0.15,
  xgb_penalty_sd  = 0.05,
  xgb_sd_cutoff   = 0.15,
  xgb_bonus_gap_ok = 0.002,
  
  rf_penalty_gap = 0.25,
  rf_penalty_sd  = 0.05,
  rf_sd_cutoff   = 0.12,
  
  # XGB search
  xgb_bo_init_points = 5,
  xgb_bo_n_iter = 8,
  xgb_nrounds_grid = c(150L, 250L, 400L),
  xgb_top_candidates_for_rounds = 2,
  
  # RF search
  rf_num_trees_tuning = 1200,
  rf_num_trees_final  = 2200,
  rf_min_node_candidates = c(3, 5, 8, 12),
  rf_bo_init_points = 5,
  rf_bo_n_iter = 8
)


# 1. Packages -------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  caret,
  dplyr,
  xgboost,
  ranger,
  doParallel,
  foreach,
  rBayesianOptimization
)

try(parallel::stopCluster(cl), silent = TRUE)

n_cores <- min(10, max(1, parallel::detectCores() - 2))
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

cat("Parallel cores for caret models:", n_cores, "\n")
cat("xgboost version:", as.character(packageVersion("xgboost")), "\n")
cat("Output folder:", output_dir, "\n")


# 2. Data -----------------------------------------------------------------------

df <- read.csv("datasets.csv")

all_features <- c(
  "Solubility", "log.Kow", "HaCount", "MW", "Complexity", "Rings",
  "Taxavalue", "Respiration", "Locomotion", "Feeding", "Size", "HLC"
)

if (!all(all_features %in% names(df))) {
  missing_cols <- setdiff(all_features, names(df))
  stop("Missing predictor columns: ", paste(missing_cols, collapse = ", "))
}

if (!"Logtox" %in% names(df)) {
  stop("Response column Logtox is missing.")
}

X_full <- df[, all_features]
y <- df$Logtox

if (!all(sapply(X_full, is.numeric))) {
  non_numeric <- names(X_full)[!sapply(X_full, is.numeric)]
  stop("Non-numeric predictor columns: ", paste(non_numeric, collapse = ", "))
}

if (!is.numeric(y)) {
  stop("The response variable must be numeric.")
}

if (any(is.na(X_full)) || any(is.infinite(as.matrix(X_full)))) {
  stop("Predictor matrix contains NA or Inf values.")
}

if (any(is.na(y)) || any(is.infinite(y))) {
  stop("Response variable contains NA or Inf values.")
}

cat("N =", nrow(X_full), "\n")
cat("Number of available predictors =", ncol(X_full), "\n")


# 3. Streamlining feature sets --------------------------------------------------

feature_sets <- list(
  M4_Core = c(
    "Solubility",
    "Size",
    "log.Kow",
    "Feeding"
  ),
  
  M8_Extended = c(
    "Solubility",
    "log.Kow",
    "Complexity",
    "Respiration",
    "Locomotion",
    "Feeding",
    "Size",
    "HLC"
  )
)

feature_path_table <- dplyr::bind_rows(
  lapply(names(feature_sets), function(nm) {
    data.frame(
      Subset = nm,
      n_vars = length(feature_sets[[nm]]),
      Features = paste(feature_sets[[nm]], collapse = "; "),
      Removed = paste(setdiff(all_features, feature_sets[[nm]]), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
)

m12_path_row <- data.frame(
  Subset = "M12_Full",
  n_vars = length(all_features),
  Features = paste(all_features, collapse = "; "),
  Removed = "",
  stringsAsFactors = FALSE
)

feature_path_table <- dplyr::bind_rows(
  feature_path_table,
  m12_path_row
)

write.csv(
  feature_path_table,
  OUT_FILE("feature_sets_m4_m8_m12"),
  row.names = FALSE
)

cat("\nFeature sets:\n")
print(feature_path_table)


# 4. Metrics --------------------------------------------------------------------

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


# 5. Cross-validation -----------------------------------------------------------

set.seed(123)

cv_index <- caret::createMultiFolds(
  y = y,
  k = CFG$cv_folds,
  times = CFG$cv_repeats
)

make_caret_seeds <- function(n_resamples, n_models = 1, base_seed = 202600) {
  set.seed(base_seed)
  seeds <- vector(mode = "list", length = n_resamples + 1)
  for (i in seq_len(n_resamples)) {
    seeds[[i]] <- sample.int(1000000L, n_models)
  }
  seeds[[n_resamples + 1]] <- sample.int(1000000L, 1)
  seeds
}

ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = CFG$cv_folds,
  repeats = CFG$cv_repeats,
  index = cv_index,
  savePredictions = "final",
  returnResamp = "final",
  summaryFunction = caret::defaultSummary,
  allowParallel = TRUE,
  seeds = make_caret_seeds(length(cv_index), n_models = 1)
)

cat("Number of CV resamples:", length(cv_index), "\n")


# 6. Helper functions -----------------------------------------------------------

summarise_resamples <- function(resample_df,
                                train_pred,
                                y_obs,
                                algo,
                                subset_name,
                                kept_features,
                                removed_features) {
  
  cv_r2 <- mean(resample_df$Rsquared, na.rm = TRUE)
  cv_rmse <- mean(resample_df$RMSE, na.rm = TRUE)
  cv_mae <- mean(resample_df$MAE, na.rm = TRUE)
  
  train_r2 <- calc_r2(train_pred, y_obs)
  train_rmse <- calc_rmse(train_pred, y_obs)
  train_mae <- calc_mae(train_pred, y_obs)
  train_nrmse <- calc_nrmse(train_pred, y_obs)
  
  gap_r2 <- train_r2 - cv_r2
  
  data.frame(
    Algorithm = algo,
    Subset = subset_name,
    n_vars = length(kept_features),
    Removed = paste(removed_features, collapse = "; "),
    Features = paste(kept_features, collapse = "; "),
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse,
    CV_nRMSE = cv_rmse / diff(range(y_obs, na.rm = TRUE)),
    CV_MAE = cv_mae,
    Train_R2 = train_r2,
    Train_RMSE = train_rmse,
    Train_nRMSE = train_nrmse,
    Train_MAE = train_mae,
    CV_R2_SD = sd(resample_df$Rsquared, na.rm = TRUE),
    CV_RMSE_SD = sd(resample_df$RMSE, na.rm = TRUE),
    CV_MAE_SD = sd(resample_df$MAE, na.rm = TRUE),
    Gap_R2 = gap_r2,
    Status = gap_status(gap_r2),
    Qualification = qualification_flag(gap_r2, cv_r2),
    stringsAsFactors = FALSE
  )
}


# 7. Prediction helper functions ------------------------------------------------

predict_m4_core_xgb <- function(bundle, newdata) {
  
  required_predictors <- bundle$input_feature_order
  
  missing_cols <- setdiff(required_predictors, names(newdata))
  if (length(missing_cols) > 0) {
    stop("Missing required predictor columns: ", paste(missing_cols, collapse = ", "))
  }
  
  X_new <- as.data.frame(newdata[, required_predictors, drop = FALSE])
  
  if (!all(sapply(X_new, is.numeric))) {
    non_numeric <- names(X_new)[!sapply(X_new, is.numeric)]
    stop("Non-numeric predictor columns: ", paste(non_numeric, collapse = ", "))
  }
  
  dnew <- xgboost::xgb.DMatrix(data = as.matrix(X_new))
  as.numeric(predict(bundle$model, newdata = dnew))
}

predict_m4_core_rf <- function(bundle, newdata) {
  
  required_predictors <- bundle$input_feature_order
  
  missing_cols <- setdiff(required_predictors, names(newdata))
  if (length(missing_cols) > 0) {
    stop("Missing required predictor columns: ", paste(missing_cols, collapse = ", "))
  }
  
  X_new <- as.data.frame(newdata[, required_predictors, drop = FALSE])
  
  if (!all(sapply(X_new, is.numeric))) {
    non_numeric <- names(X_new)[!sapply(X_new, is.numeric)]
    stop("Non-numeric predictor columns: ", paste(non_numeric, collapse = ", "))
  }
  
  as.numeric(predict(bundle$model, newdata = X_new))
}


# 8. Random Forest functions ----------------------------------------------------

get_rf_rough_grid <- function(p) {
  
  mtry_candidates <- sort(unique(c(
    1,
    max(1, floor(sqrt(p))),
    max(1, ceiling(p / 2)),
    p
  )))
  
  expand.grid(
    mtry = mtry_candidates,
    min.node.size = CFG$rf_min_node_candidates,
    stringsAsFactors = FALSE
  )
}

get_rf_local_bounds <- function(center, p) {
  
  list(
    mtry = c(
      max(1L, as.integer(round(center$mtry)) - 2L),
      min(p, as.integer(round(center$mtry)) + 2L)
    ),
    min_node_size = c(
      max(1L, as.integer(round(center$min.node.size)) - 4L),
      min(20L, as.integer(round(center$min.node.size)) + 4L)
    )
  )
}

evaluate_rf_combo <- function(X_sub,
                              y,
                              ctrl,
                              subset_name,
                              removed_features,
                              mtry,
                              min_node_size,
                              num_trees,
                              stage = "rf_candidate") {
  
  p <- ncol(X_sub)
  
  mtry <- as.integer(round(mtry))
  min_node_size <- as.integer(round(min_node_size))
  
  mtry <- max(1L, min(mtry, p))
  min_node_size <- max(1L, min(min_node_size, 30L))
  
  candidate_seed <- 2026L +
    sum(utf8ToInt(subset_name)) +
    100L * mtry +
    10L * min_node_size +
    ifelse(stage == "final_evaluation", 50000L, 0L)
  
  set.seed(candidate_seed)
  
  model_rf <- tryCatch({
    caret::train(
      x = as.data.frame(X_sub),
      y = y,
      method = "ranger",
      trControl = ctrl,
      num.trees = num_trees,
      num.threads = 1,
      importance = "none",
      tuneGrid = data.frame(
        mtry = mtry,
        min.node.size = min_node_size,
        splitrule = "variance"
      ),
      metric = "Rsquared"
    )
  }, error = function(e) NULL)
  
  if (is.null(model_rf)) return(NULL)
  
  train_pred <- predict(model_rf, newdata = as.data.frame(X_sub))
  
  summary_row <- summarise_resamples(
    resample_df = model_rf$resample,
    train_pred = train_pred,
    y_obs = y,
    algo = "RF",
    subset_name = subset_name,
    kept_features = colnames(X_sub),
    removed_features = removed_features
  )
  
  score <- generalization_score(
    cv_r2 = summary_row$CV_R2,
    gap_r2 = summary_row$Gap_R2,
    cv_r2_sd = summary_row$CV_R2_SD,
    gap_penalty = CFG$rf_penalty_gap,
    sd_penalty = CFG$rf_penalty_sd,
    gap_cutoff = 0.20,
    sd_cutoff = CFG$rf_sd_cutoff
  )
  
  history_row <- data.frame(
    Subset = subset_name,
    n_vars = p,
    stage = stage,
    mtry = mtry,
    min.node.size = min_node_size,
    num.trees = num_trees,
    candidate_seed = candidate_seed,
    CV_R2 = summary_row$CV_R2,
    CV_RMSE = summary_row$CV_RMSE,
    CV_nRMSE = summary_row$CV_nRMSE,
    CV_MAE = summary_row$CV_MAE,
    Train_R2 = summary_row$Train_R2,
    Train_RMSE = summary_row$Train_RMSE,
    Train_nRMSE = summary_row$Train_nRMSE,
    Train_MAE = summary_row$Train_MAE,
    CV_R2_SD = summary_row$CV_R2_SD,
    CV_RMSE_SD = summary_row$CV_RMSE_SD,
    CV_MAE_SD = summary_row$CV_MAE_SD,
    Gap_R2 = summary_row$Gap_R2,
    Score = score,
    Status = summary_row$Status,
    Qualification = summary_row$Qualification,
    stringsAsFactors = FALSE
  )
  
  list(
    model = model_rf,
    summary = summary_row,
    history = history_row,
    bestTune = model_rf$bestTune,
    score = score,
    candidate_seed = candidate_seed
  )
}

run_rf_subset <- function(X_sub,
                          y,
                          ctrl,
                          subset_name,
                          removed_features) {
  
  p <- ncol(X_sub)
  
  rough_grid <- get_rf_rough_grid(p)
  rough_history <- data.frame()
  
  cat("\n>>> RF rough grid:", subset_name, "\n")
  
  for (i in seq_len(nrow(rough_grid))) {
    
    row_i <- rough_grid[i, , drop = FALSE]
    
    cat(sprintf(
      "RF rough grid %s | %d/%d | mtry = %d | min.node.size = %d\n",
      subset_name,
      i,
      nrow(rough_grid),
      row_i$mtry,
      row_i$min.node.size
    ))
    
    res_i <- evaluate_rf_combo(
      X_sub = X_sub,
      y = y,
      ctrl = ctrl,
      subset_name = subset_name,
      removed_features = removed_features,
      mtry = row_i$mtry,
      min_node_size = row_i$min.node.size,
      num_trees = CFG$rf_num_trees_tuning,
      stage = "rough_grid"
    )
    
    if (!is.null(res_i)) {
      rough_history <- dplyr::bind_rows(rough_history, res_i$history)
    }
  }
  
  if (nrow(rough_history) == 0) {
    stop("No valid RF rough-grid result for ", subset_name)
  }
  
  rough_history <- rough_history %>%
    dplyr::arrange(desc(Score), desc(CV_R2), Gap_R2, CV_R2_SD)
  
  rough_best <- rough_history[1, , drop = FALSE]
  
  bounds <- get_rf_local_bounds(
    center = rough_best,
    p = p
  )
  
  bo_history <- data.frame()
  
  rf_bo_objective <- function(mtry, min_node_size) {
    
    out <- tryCatch({
      
      res <- evaluate_rf_combo(
        X_sub = X_sub,
        y = y,
        ctrl = ctrl,
        subset_name = subset_name,
        removed_features = removed_features,
        mtry = mtry,
        min_node_size = min_node_size,
        num_trees = CFG$rf_num_trees_tuning,
        stage = "bayesopt"
      )
      
      if (is.null(res)) stop("RF evaluation returned NULL.")
      
      bo_history <<- dplyr::bind_rows(bo_history, res$history)
      
      list(
        Score = res$score,
        Pred = 0
      )
      
    }, error = function(e) {
      list(
        Score = -999,
        Pred = 0
      )
    })
    
    out
  }
  
  cat("\n>>> RF local Bayesian optimization:", subset_name, "\n")
  
  set.seed(2026)
  
  suppressWarnings(
    tryCatch({
      BayesianOptimization(
        FUN = rf_bo_objective,
        bounds = bounds,
        init_points = CFG$rf_bo_init_points,
        n_iter = CFG$rf_bo_n_iter,
        acq = "ucb",
        kappa = 2.0,
        eps = 0.001,
        verbose = TRUE
      )
    }, error = function(e) {
      message("RF BayesianOptimization stopped for ", subset_name, ": ", conditionMessage(e))
      NULL
    })
  )
  
  combined_history <- dplyr::bind_rows(rough_history, bo_history) %>%
    dplyr::filter(
      is.finite(CV_R2),
      is.finite(Gap_R2),
      is.finite(CV_R2_SD),
      is.finite(Score)
    )
  
  if (nrow(combined_history) == 0) {
    stop("No valid RF combined-history result for ", subset_name)
  }
  
  strict_pool <- combined_history %>%
    dplyr::filter(Gap_R2 < 0.20) %>%
    dplyr::arrange(desc(CV_R2), CV_RMSE, CV_MAE, Gap_R2, CV_R2_SD)
  
  if (nrow(strict_pool) > 0) {
    
    best_rf <- strict_pool %>% dplyr::slice(1)
    selection_level <- "Gap_R2 < 0.20"
    
  } else {
    
    best_rf <- combined_history %>%
      dplyr::arrange(desc(CV_R2), CV_RMSE, CV_MAE, Gap_R2, CV_R2_SD) %>%
      dplyr::slice(1)
    
    selection_level <- "Fallback: highest CV_R2 among tuned RF candidates"
  }
  
  cat("\n>>> RF final evaluation:", subset_name, "\n")
  cat(
    "Selected RF params:",
    "mtry =", best_rf$mtry,
    "| min.node.size =", best_rf$min.node.size,
    "| selection =", selection_level,
    "\n"
  )
  
  final_res <- evaluate_rf_combo(
    X_sub = X_sub,
    y = y,
    ctrl = ctrl,
    subset_name = subset_name,
    removed_features = removed_features,
    mtry = best_rf$mtry,
    min_node_size = best_rf$min.node.size,
    num_trees = CFG$rf_num_trees_final,
    stage = "final_evaluation"
  )
  
  if (is.null(final_res)) {
    stop("Final RF evaluation failed for ", subset_name)
  }
  
  final_summary <- final_res$summary %>%
    dplyr::mutate(
      Selection_Level = selection_level,
      mtry = as.integer(round(best_rf$mtry)),
      min.node.size = as.integer(round(best_rf$min.node.size)),
      num.trees = CFG$rf_num_trees_final,
      candidate_seed = final_res$candidate_seed,
      .after = n_vars
    )
  
  final_history <- dplyr::bind_rows(
    rough_history,
    bo_history,
    final_res$history
  )
  
  list(
    model = final_res$model,
    summary = final_summary,
    history = final_history,
    bestTune = final_res$model$bestTune,
    candidate_seed = final_res$candidate_seed
  )
}


# 9. XGBoost functions ----------------------------------------------------------

get_xgb_rough_grid <- function() {
  expand.grid(
    eta = c(0.02, 0.04),
    gamma = c(0.00, 0.20),
    colsample_bytree = 0.80,
    subsample = c(0.80, 1.00),
    lambda = c(0.00, 1.00),
    alpha = c(0.00, 0.10),
    stringsAsFactors = FALSE
  )
}

get_xgb_local_bounds <- function(center) {
  list(
    eta = c(
      max(0.01, center$eta * 0.8),
      min(0.08, center$eta * 1.2)
    ),
    gamma = c(
      max(0.00, center$gamma - 0.15),
      min(1.50, center$gamma + 0.15)
    ),
    colsample_bytree = c(
      max(0.65, center$colsample_bytree - 0.08),
      min(1.00, center$colsample_bytree + 0.08)
    ),
    subsample = c(
      max(0.70, center$subsample - 0.08),
      min(1.00, center$subsample + 0.08)
    ),
    lambda = c(
      max(0.00, center$lambda * 0.5),
      min(5.00, center$lambda * 1.5 + 1e-8)
    ),
    alpha = c(
      max(0.00, center$alpha - 0.10),
      min(2.00, center$alpha + 0.10)
    )
  )
}

evaluate_xgb_combo <- function(X_mat,
                               y,
                               cv_index,
                               params_row,
                               nrounds_mode = "earlystop",
                               fixed_nrounds = NULL,
                               seed_offset = 0) {
  
  resample_metrics <- data.frame()
  best_rounds_vec <- c()
  
  params <- list(
    objective = "reg:squarederror",
    booster = "gbtree",
    max_depth = 4L,
    min_child_weight = 2L,
    eta = params_row$eta,
    gamma = params_row$gamma,
    colsample_bytree = params_row$colsample_bytree,
    subsample = params_row$subsample,
    lambda = params_row$lambda,
    alpha = params_row$alpha,
    nthread = 1,
    seed = 123,
    seed_per_iteration = FALSE
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
    
    if (nrounds_mode == "earlystop") {
      
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
          nrounds = 1500,
          watchlist = list(train = dtrain_inner, valid = dvalid_inner),
          early_stopping_rounds = 100,
          maximize = FALSE,
          verbose = 0
        )
      }, error = function(e) NULL)
      
      if (is.null(model_inner)) return(NULL)
      
      best_nrounds <- model_inner$best_iteration
      
      if (is.null(best_nrounds) || is.na(best_nrounds) || best_nrounds <= 0) {
        best_nrounds <- 200L
      }
      
    } else {
      
      best_nrounds <- fixed_nrounds
    }
    
    best_rounds_vec <- c(best_rounds_vec, best_nrounds)
    
    fold_params <- params
    fold_params$seed <- 2026L + seed_offset + i
    
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
  train_rmse <- calc_rmse(pred_full, y)
  train_mae <- calc_mae(pred_full, y)
  train_nrmse <- calc_nrmse(pred_full, y)
  gap_r2 <- train_r2 - cv_r2
  
  score <- generalization_score(
    cv_r2 = cv_r2,
    gap_r2 = gap_r2,
    cv_r2_sd = cv_r2_sd,
    gap_penalty = CFG$xgb_penalty_gap,
    sd_penalty = CFG$xgb_penalty_sd,
    gap_cutoff = 0.20,
    sd_cutoff = CFG$xgb_sd_cutoff,
    bonus_if_ok = CFG$xgb_bonus_gap_ok
  )
  
  data.frame(
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse,
    CV_nRMSE = cv_rmse / diff(range(y, na.rm = TRUE)),
    CV_MAE = cv_mae,
    Train_R2 = train_r2,
    Train_RMSE = train_rmse,
    Train_nRMSE = train_nrmse,
    Train_MAE = train_mae,
    CV_R2_SD = cv_r2_sd,
    CV_RMSE_SD = sd(resample_metrics$RMSE, na.rm = TRUE),
    CV_MAE_SD = sd(resample_metrics$MAE, na.rm = TRUE),
    Gap_R2 = gap_r2,
    Mean_nrounds = mean(best_rounds_vec, na.rm = TRUE),
    Median_nrounds = final_nrounds,
    Score = score,
    Status = gap_status(gap_r2),
    Qualification = qualification_flag(gap_r2, cv_r2),
    stringsAsFactors = FALSE
  )
}

run_xgb_subset <- function(X_sub,
                           y,
                           cv_index,
                           subset_name,
                           removed_features) {
  
  X_mat <- as.matrix(X_sub)
  
  rough_grid <- get_xgb_rough_grid()
  rough_history <- data.frame()
  
  cat("\n>>> XGB rough grid:", subset_name, "\n")
  
  for (i in seq_len(nrow(rough_grid))) {
    
    p <- rough_grid[i, , drop = FALSE]
    
    cat(sprintf(
      "XGB rough grid %s | %d/%d\n",
      subset_name,
      i,
      nrow(rough_grid)
    ))
    
    res <- tryCatch({
      evaluate_xgb_combo(
        X_mat = X_mat,
        y = y,
        cv_index = cv_index,
        params_row = p,
        nrounds_mode = "earlystop",
        seed_offset = 0
      )
    }, error = function(e) NULL)
    
    if (!is.null(res)) {
      rough_history <- dplyr::bind_rows(
        rough_history,
        cbind(
          data.frame(
            Subset = subset_name,
            n_vars = ncol(X_sub),
            stage = "rough_grid",
            eta = p$eta,
            gamma = p$gamma,
            colsample_bytree = p$colsample_bytree,
            subsample = p$subsample,
            lambda = p$lambda,
            alpha = p$alpha,
            stringsAsFactors = FALSE
          ),
          res
        )
      )
    }
  }
  
  if (nrow(rough_history) == 0) {
    stop("No valid XGB rough-grid result for ", subset_name)
  }
  
  rough_history <- rough_history %>%
    dplyr::arrange(desc(CV_R2), Gap_R2, CV_R2_SD)
  
  rough_best <- rough_history[1, ]
  
  center <- list(
    eta = rough_best$eta,
    gamma = rough_best$gamma,
    colsample_bytree = rough_best$colsample_bytree,
    subsample = rough_best$subsample,
    lambda = rough_best$lambda,
    alpha = rough_best$alpha
  )
  
  bounds <- get_xgb_local_bounds(center)
  bo_history <- data.frame()
  
  xgb_bo_objective <- function(gamma,
                               eta,
                               colsample_bytree,
                               subsample,
                               lambda,
                               alpha) {
    
    p <- data.frame(
      gamma = gamma,
      eta = eta,
      colsample_bytree = colsample_bytree,
      subsample = subsample,
      lambda = lambda,
      alpha = alpha
    )
    
    out <- tryCatch({
      
      res <- evaluate_xgb_combo(
        X_mat = X_mat,
        y = y,
        cv_index = cv_index,
        params_row = p,
        nrounds_mode = "earlystop",
        seed_offset = 0
      )
      
      if (is.null(res)) stop("Evaluation returned NULL.")
      
      bo_history <<- dplyr::bind_rows(
        bo_history,
        cbind(
          data.frame(
            Subset = subset_name,
            n_vars = ncol(X_sub),
            stage = "bayesopt",
            eta = eta,
            gamma = gamma,
            colsample_bytree = colsample_bytree,
            subsample = subsample,
            lambda = lambda,
            alpha = alpha,
            stringsAsFactors = FALSE
          ),
          res
        )
      )
      
      list(Score = res$Score[1], Pred = 0)
      
    }, error = function(e) {
      list(Score = -999, Pred = 0)
    })
    
    out
  }
  
  cat("\n>>> XGB local Bayesian optimization:", subset_name, "\n")
  
  set.seed(2026)
  
  suppressWarnings(
    tryCatch({
      BayesianOptimization(
        FUN = xgb_bo_objective,
        bounds = bounds,
        init_points = CFG$xgb_bo_init_points,
        n_iter = CFG$xgb_bo_n_iter,
        acq = "ucb",
        kappa = 2.0,
        eps = 0.001,
        verbose = TRUE
      )
    }, error = function(e) {
      message("XGB BayesianOptimization stopped for ", subset_name, ": ", conditionMessage(e))
      NULL
    })
  )
  
  combined_history <- dplyr::bind_rows(rough_history, bo_history) %>%
    dplyr::filter(
      is.finite(CV_R2),
      is.finite(Gap_R2),
      is.finite(CV_R2_SD)
    ) %>%
    dplyr::arrange(desc(CV_R2), Gap_R2, CV_R2_SD)
  
  if (nrow(combined_history) == 0) {
    stop("No valid XGB combined-history result for ", subset_name)
  }
  
  top_candidates <- head(
    combined_history,
    min(CFG$xgb_top_candidates_for_rounds, nrow(combined_history))
  )
  
  round_tuning_results <- data.frame()
  
  cat("\n>>> XGB nrounds tuning:", subset_name, "\n")
  
  for (i in seq_len(nrow(top_candidates))) {
    
    row_df <- top_candidates[i, , drop = FALSE]
    
    params_row <- row_df[, c(
      "eta",
      "gamma",
      "colsample_bytree",
      "subsample",
      "lambda",
      "alpha"
    )]
    
    for (nr in CFG$xgb_nrounds_grid) {
      
      cat(sprintf(
        "XGB nrounds %s | candidate %d/%d | nrounds = %d\n",
        subset_name,
        i,
        nrow(top_candidates),
        nr
      ))
      
      tmp <- tryCatch({
        evaluate_xgb_combo(
          X_mat = X_mat,
          y = y,
          cv_index = cv_index,
          params_row = params_row,
          nrounds_mode = "fixed",
          fixed_nrounds = nr,
          seed_offset = 0
        )
      }, error = function(e) NULL)
      
      if (!is.null(tmp)) {
        round_tuning_results <- dplyr::bind_rows(
          round_tuning_results,
          cbind(
            data.frame(
              Subset = subset_name,
              n_vars = ncol(X_sub),
              stage = "nrounds_grid",
              eta = row_df$eta,
              gamma = row_df$gamma,
              colsample_bytree = row_df$colsample_bytree,
              subsample = row_df$subsample,
              lambda = row_df$lambda,
              alpha = row_df$alpha,
              stringsAsFactors = FALSE
            ),
            tmp
          )
        )
      }
    }
  }
  
  if (nrow(round_tuning_results) == 0) {
    stop("No valid XGB nrounds tuning result for ", subset_name)
  }
  
  strict_pool <- round_tuning_results %>%
    dplyr::filter(
      is.finite(CV_R2),
      is.finite(Gap_R2),
      is.finite(CV_R2_SD),
      Gap_R2 < 0.20
    ) %>%
    dplyr::arrange(desc(CV_R2), CV_RMSE, CV_MAE, Gap_R2, CV_R2_SD)
  
  if (nrow(strict_pool) > 0) {
    
    best_xgb <- strict_pool %>% dplyr::slice(1)
    selection_level <- "Gap_R2 < 0.20"
    
  } else {
    
    best_xgb <- round_tuning_results %>%
      dplyr::filter(
        is.finite(CV_R2),
        is.finite(Gap_R2),
        is.finite(CV_R2_SD)
      ) %>%
      dplyr::arrange(desc(CV_R2), CV_RMSE, CV_MAE, Gap_R2, CV_R2_SD) %>%
      dplyr::slice(1)
    
    selection_level <- "Fallback: highest CV_R2 among tuned XGB candidates"
  }
  
  final_params <- list(
    objective = "reg:squarederror",
    booster = "gbtree",
    max_depth = 4L,
    eta = best_xgb$eta,
    gamma = best_xgb$gamma,
    colsample_bytree = best_xgb$colsample_bytree,
    subsample = best_xgb$subsample,
    min_child_weight = 2L,
    lambda = best_xgb$lambda,
    alpha = best_xgb$alpha,
    nthread = 1,
    seed = 2026,
    seed_per_iteration = FALSE
  )
  
  final_nrounds <- as.integer(best_xgb$Median_nrounds)
  
  dtrain_all <- xgb.DMatrix(data = X_mat, label = y)
  
  model_final <- xgb.train(
    params = final_params,
    data = dtrain_all,
    nrounds = final_nrounds,
    verbose = 0
  )
  
  train_pred <- predict(model_final, dtrain_all)
  
  train_r2 <- calc_r2(train_pred, y)
  train_rmse <- calc_rmse(train_pred, y)
  train_mae <- calc_mae(train_pred, y)
  train_nrmse <- calc_nrmse(train_pred, y)
  gap_r2_final <- train_r2 - best_xgb$CV_R2
  
  final_summary <- data.frame(
    Algorithm = "XGB",
    Subset = subset_name,
    n_vars = ncol(X_sub),
    Removed = paste(removed_features, collapse = "; "),
    Features = paste(colnames(X_sub), collapse = "; "),
    Selection_Level = selection_level,
    CV_R2 = best_xgb$CV_R2,
    CV_RMSE = best_xgb$CV_RMSE,
    CV_nRMSE = best_xgb$CV_nRMSE,
    CV_MAE = best_xgb$CV_MAE,
    Train_R2 = train_r2,
    Train_RMSE = train_rmse,
    Train_nRMSE = train_nrmse,
    Train_MAE = train_mae,
    CV_R2_SD = best_xgb$CV_R2_SD,
    CV_RMSE_SD = best_xgb$CV_RMSE_SD,
    CV_MAE_SD = best_xgb$CV_MAE_SD,
    Gap_R2 = gap_r2_final,
    Mean_nrounds = best_xgb$Mean_nrounds,
    Median_nrounds = best_xgb$Median_nrounds,
    Status = gap_status(gap_r2_final),
    Qualification = qualification_flag(gap_r2_final, best_xgb$CV_R2),
    stringsAsFactors = FALSE
  )
  
  list(
    model = model_final,
    summary = final_summary,
    history = dplyr::bind_rows(rough_history, bo_history, round_tuning_results),
    final_params = final_params,
    final_nrounds = final_nrounds
  )
}


# 10. Run M4_Core and M8_Extended -----------------------------------------------

xgb_path_summary <- data.frame()
rf_path_summary  <- data.frame()
xgb_path_history <- data.frame()
rf_path_history  <- data.frame()

saved_models <- list()
prediction_bundles <- list()

for (subset_name in names(feature_sets)) {
  
  kept_features <- feature_sets[[subset_name]]
  removed_features <- setdiff(all_features, kept_features)
  X_sub <- df[, kept_features, drop = FALSE]
  
  cat("\n============================================================\n")
  cat("Running subset:", subset_name, "| n_vars =", length(kept_features), "\n")
  cat("Features:", paste(kept_features, collapse = ", "), "\n")
  
  cat("\n>>> RF on", subset_name, "\n")
  
  rf_res <- run_rf_subset(
    X_sub = X_sub,
    y = y,
    ctrl = ctrl,
    subset_name = subset_name,
    removed_features = removed_features
  )
  
  rf_path_summary <- dplyr::bind_rows(rf_path_summary, rf_res$summary)
  rf_path_history <- dplyr::bind_rows(rf_path_history, rf_res$history)
  
  rf_model_key <- paste0("RF_", subset_name, "_final_model")
  saved_models[[rf_model_key]] <- rf_res$model
  
  prediction_bundles[[paste0("RF_", subset_name, "_prediction_bundle")]] <- list(
    algorithm = "RF",
    subset = subset_name,
    n_vars = length(kept_features),
    features = kept_features,
    removed = removed_features,
    response_name = "Logtox",
    response_scale = "logtox",
    model = rf_res$model,
    best_tune = rf_res$bestTune,
    candidate_seed = rf_res$candidate_seed,
    input_feature_order = kept_features,
    predict_fun = predict_m4_core_rf,
    notes = "Prediction bundle for streamlined RF model. Predictions are on the Logtox scale."
  )
  
  cat("\n>>> XGB on", subset_name, "\n")
  
  xgb_res <- run_xgb_subset(
    X_sub = X_sub,
    y = y,
    cv_index = cv_index,
    subset_name = subset_name,
    removed_features = removed_features
  )
  
  xgb_path_summary <- dplyr::bind_rows(xgb_path_summary, xgb_res$summary)
  xgb_path_history <- dplyr::bind_rows(xgb_path_history, xgb_res$history)
  
  xgb_model_key <- paste0("XGB_", subset_name, "_final_model")
  saved_models[[xgb_model_key]] <- xgb_res$model
  
  prediction_bundles[[paste0("XGB_", subset_name, "_prediction_bundle")]] <- list(
    algorithm = "XGB",
    subset = subset_name,
    n_vars = length(kept_features),
    features = kept_features,
    removed = removed_features,
    response_name = "Logtox",
    response_scale = "logtox",
    model = xgb_res$model,
    final_params = xgb_res$final_params,
    final_nrounds = xgb_res$final_nrounds,
    input_feature_order = kept_features,
    predict_fun = predict_m4_core_xgb,
    notes = "Prediction bundle for streamlined native XGBoost model. Predictions are on the Logtox scale."
  )
}


# 11. Export M4/M8 summaries and search histories -------------------------------

streamlining_summary_trained <- dplyr::bind_rows(
  xgb_path_summary,
  rf_path_summary
)

write.csv(
  streamlining_summary_trained,
  OUT_FILE("streamlining_summary_m4_m8_trained"),
  row.names = FALSE
)

write.csv(
  xgb_path_summary,
  OUT_FILE("xgb_streamlining_summary_m4_m8"),
  row.names = FALSE
)

write.csv(
  rf_path_summary,
  OUT_FILE("rf_streamlining_summary_m4_m8"),
  row.names = FALSE
)

write.csv(
  xgb_path_history,
  OUT_FILE("xgb_search_history_m4_m8"),
  row.names = FALSE
)

write.csv(
  rf_path_history,
  OUT_FILE("rf_search_history_m4_m8"),
  row.names = FALSE
)


# 12. Save model objects and prediction bundles ---------------------------------

subset_file_prefix <- function(subset_name) {
  switch(
    subset_name,
    "M4_Core" = "m4_core",
    "M8_Extended" = "m8_extended",
    tolower(subset_name)
  )
}

for (subset_name in names(feature_sets)) {
  
  prefix <- subset_file_prefix(subset_name)
  
  xgb_model_key <- paste0("XGB_", subset_name, "_final_model")
  rf_model_key  <- paste0("RF_", subset_name, "_final_model")
  
  xgb_bundle_key <- paste0("XGB_", subset_name, "_prediction_bundle")
  rf_bundle_key  <- paste0("RF_", subset_name, "_prediction_bundle")
  
  if (!is.null(saved_models[[xgb_model_key]])) {
    saveRDS(
      saved_models[[xgb_model_key]],
      OUT_FILE(paste0(prefix, "_xgb_model"), ext = "rds")
    )
  }
  
  if (!is.null(saved_models[[rf_model_key]])) {
    saveRDS(
      saved_models[[rf_model_key]],
      OUT_FILE(paste0(prefix, "_rf_model"), ext = "rds")
    )
  }
  
  if (!is.null(prediction_bundles[[xgb_bundle_key]])) {
    saveRDS(
      prediction_bundles[[xgb_bundle_key]],
      OUT_FILE(paste0(prefix, "_xgb_prediction_bundle"), ext = "rds")
    )
  }
  
  if (!is.null(prediction_bundles[[rf_bundle_key]])) {
    saveRDS(
      prediction_bundles[[rf_bundle_key]],
      OUT_FILE(paste0(prefix, "_rf_prediction_bundle"), ext = "rds")
    )
  }
}

# Dedicated final M4_Core prediction bundles, saved with the same clean names.
saveRDS(
  prediction_bundles[["XGB_M4_Core_prediction_bundle"]],
  OUT_FILE("m4_core_xgb_prediction_bundle", ext = "rds")
)

saveRDS(
  prediction_bundles[["RF_M4_Core_prediction_bundle"]],
  OUT_FILE("m4_core_rf_prediction_bundle", ext = "rds")
)


# 13. Traditional pooled CV R2 for M4/M8 ----------------------------------------

calc_true_r2 <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  pred <- pred[ok]
  obs  <- obs[ok]
  
  if (length(obs) < 2) return(NA_real_)
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - (ss_res / ss_tot)
}

safe_calc_rf_traditional_r2 <- function(model_obj) {
  
  if (is.null(model_obj)) return(NA_real_)
  if (is.null(model_obj$pred)) return(NA_real_)
  if (!is.data.frame(model_obj$pred) || nrow(model_obj$pred) == 0) return(NA_real_)
  if (is.null(model_obj$bestTune)) return(NA_real_)
  
  pred_df <- model_obj$pred
  tune_cols <- names(model_obj$bestTune)
  
  for (nm in tune_cols) {
    if (!nm %in% names(pred_df)) return(NA_real_)
    pred_df <- pred_df[pred_df[[nm]] == model_obj$bestTune[[nm]], , drop = FALSE]
  }
  
  if (nrow(pred_df) == 0) return(NA_real_)
  if (!all(c("rowIndex", "obs", "pred") %in% names(pred_df))) return(NA_real_)
  
  pooled_pred <- pred_df %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      obs = mean(obs, na.rm = TRUE),
      pred = mean(pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(pooled_pred) == 0) return(NA_real_)
  calc_true_r2(pooled_pred$pred, pooled_pred$obs)
}

safe_calc_xgb_traditional_r2 <- function(X_sub,
                                         y,
                                         cv_index,
                                         final_params,
                                         final_nrounds) {
  
  X_mat <- as.matrix(X_sub)
  xgb_oof_pred <- data.frame()
  
  for (i in seq_along(cv_index)) {
    
    train_idx <- cv_index[[i]]
    test_idx  <- setdiff(seq_along(y), train_idx)
    
    dtrain_i <- xgb.DMatrix(
      data = X_mat[train_idx, , drop = FALSE],
      label = y[train_idx]
    )
    
    dtest_i <- xgb.DMatrix(
      data = X_mat[test_idx, , drop = FALSE],
      label = y[test_idx]
    )
    
    fold_params <- final_params
    fold_params$seed <- 2026L + i
    
    model_i <- tryCatch({
      xgb.train(
        params = fold_params,
        data = dtrain_i,
        nrounds = final_nrounds,
        verbose = 0
      )
    }, error = function(e) NULL)
    
    if (is.null(model_i)) next
    
    pred_i <- tryCatch({
      predict(model_i, dtest_i)
    }, error = function(e) NULL)
    
    if (is.null(pred_i)) next
    
    xgb_oof_pred <- rbind(
      xgb_oof_pred,
      data.frame(
        rowIndex = test_idx,
        obs = y[test_idx],
        pred = pred_i,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (!is.data.frame(xgb_oof_pred) || nrow(xgb_oof_pred) == 0) return(NA_real_)
  
  pooled_pred <- xgb_oof_pred %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      obs = mean(obs, na.rm = TRUE),
      pred = mean(pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(pooled_pred) == 0) return(NA_real_)
  calc_true_r2(pooled_pred$pred, pooled_pred$obs)
}

traditional_r2_df <- data.frame()

for (subset_name in names(feature_sets)) {
  
  kept_features <- feature_sets[[subset_name]]
  X_sub <- df[, kept_features, drop = FALSE]
  
  rf_model <- saved_models[[paste0("RF_", subset_name, "_final_model")]]
  rf_trad_r2 <- safe_calc_rf_traditional_r2(rf_model)
  
  traditional_r2_df <- rbind(
    traditional_r2_df,
    data.frame(
      Algorithm = "RF",
      Subset = subset_name,
      CV_R2_Traditional = rf_trad_r2,
      stringsAsFactors = FALSE
    )
  )
  
  xgb_bundle <- prediction_bundles[[paste0("XGB_", subset_name, "_prediction_bundle")]]
  
  xgb_trad_r2 <- safe_calc_xgb_traditional_r2(
    X_sub = X_sub,
    y = y,
    cv_index = cv_index,
    final_params = xgb_bundle$final_params,
    final_nrounds = xgb_bundle$final_nrounds
  )
  
  traditional_r2_df <- rbind(
    traditional_r2_df,
    data.frame(
      Algorithm = "XGB",
      Subset = subset_name,
      CV_R2_Traditional = xgb_trad_r2,
      stringsAsFactors = FALSE
    )
  )
}

xgb_path_summary <- xgb_path_summary %>%
  dplyr::left_join(
    traditional_r2_df %>% dplyr::filter(Algorithm == "XGB"),
    by = c("Algorithm", "Subset")
  )

rf_path_summary <- rf_path_summary %>%
  dplyr::left_join(
    traditional_r2_df %>% dplyr::filter(Algorithm == "RF"),
    by = c("Algorithm", "Subset")
  )

streamlining_summary_trained <- dplyr::bind_rows(
  xgb_path_summary,
  rf_path_summary
)

write.csv(
  streamlining_summary_trained,
  OUT_FILE("streamlining_summary_m4_m8_trained"),
  row.names = FALSE
)

write.csv(
  traditional_r2_df,
  OUT_FILE("traditional_oof_r2_m4_m8"),
  row.names = FALSE
)


# 14. Read M12_Full performance -------------------------------------------------

standard_cols <- c(
  "Algorithm",
  "Subset",
  "n_vars",
  "Selection_Level",
  "mtry",
  "min.node.size",
  "num.trees",
  "Mean_nrounds",
  "Median_nrounds",
  "candidate_seed",
  "Removed",
  "Features",
  "CV_R2",
  "CV_R2_Traditional",
  "CV_RMSE",
  "CV_nRMSE",
  "CV_MAE",
  "Train_R2",
  "Train_RMSE",
  "Train_nRMSE",
  "Train_MAE",
  "CV_R2_SD",
  "CV_RMSE_SD",
  "CV_MAE_SD",
  "Gap_R2",
  "Status",
  "Qualification"
)

align_summary_cols <- function(dat) {
  for (nm in standard_cols) {
    if (!nm %in% names(dat)) {
      dat[[nm]] <- NA
    }
  }
  dat[, standard_cols, drop = FALSE]
}

read_m12_summary <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop(
      "Cannot find M12 full-model summary file: ",
      file_path,
      "\nPlease put Model_Performance_Summary_FINAL.csv in the working directory."
    )
  }
  
  raw <- read.csv(file_path, stringsAsFactors = FALSE)
  
  if (!"Model" %in% names(raw)) {
    stop("M12 summary file must contain a 'Model' column.")
  }
  
  m12_raw <- raw %>%
    dplyr::filter(
      Model %in% c("XGB_Bayes", "RF_Bayes") |
        grepl("XGB", Model, ignore.case = TRUE) |
        grepl("^RF", Model, ignore.case = TRUE)
    ) %>%
    dplyr::mutate(
      Algorithm = dplyr::case_when(
        grepl("XGB", Model, ignore.case = TRUE) ~ "XGB",
        grepl("RF", Model, ignore.case = TRUE) ~ "RF",
        TRUE ~ Model
      ),
      Subset = "M12_Full",
      n_vars = length(all_features),
      Selection_Level = "Read from previous full-model training summary",
      Removed = "",
      Features = paste(all_features, collapse = "; ")
    ) %>%
    dplyr::filter(Algorithm %in% c("XGB", "RF"))
  
  if (nrow(m12_raw) == 0) {
    stop("No XGB/RF rows found in the M12 full-model summary file.")
  }
  
  align_summary_cols(m12_raw)
}

m12_summary <- read_m12_summary(M12_METRICS_FILE)

write.csv(
  m12_summary,
  OUT_FILE("m12_full_readin_summary"),
  row.names = FALSE
)


# 15. Combine M4 / M8 / M12 -----------------------------------------------------

combined_summary <- dplyr::bind_rows(
  align_summary_cols(streamlining_summary_trained),
  m12_summary
)

combined_summary <- combined_summary %>%
  dplyr::mutate(
    Subset = factor(
      Subset,
      levels = c("M4_Core", "M8_Extended", "M12_Full")
    )
  ) %>%
  dplyr::arrange(
    Algorithm,
    Subset
  ) %>%
  dplyr::mutate(
    Subset = as.character(Subset)
  )

comparison_table <- combined_summary %>%
  dplyr::group_by(Algorithm) %>%
  dplyr::mutate(
    M12_CV_R2 = CV_R2[Subset == "M12_Full"][1],
    Delta_CV_R2_vs_M12 = CV_R2 - M12_CV_R2,
    M12_CV_RMSE = CV_RMSE[Subset == "M12_Full"][1],
    Delta_CV_RMSE_vs_M12 = CV_RMSE - M12_CV_RMSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Algorithm,
    Subset,
    n_vars,
    CV_R2,
    CV_R2_Traditional,
    Delta_CV_R2_vs_M12,
    CV_RMSE,
    Delta_CV_RMSE_vs_M12,
    CV_nRMSE,
    CV_MAE,
    Train_R2,
    Gap_R2,
    CV_R2_SD,
    Status,
    Qualification,
    Features,
    Removed
  )

write.csv(
  combined_summary,
  OUT_FILE("streamlining_summary_m4_m8_m12"),
  row.names = FALSE
)

write.csv(
  comparison_table,
  OUT_FILE("streamlining_comparison_m4_m8_m12"),
  row.names = FALSE
)

cat("\n>>> Combined comparison table:\n")
print(comparison_table, digits = 4)


# 16. Recommendation table ------------------------------------------------------

xgb_m4 <- comparison_table %>%
  dplyr::filter(Algorithm == "XGB", Subset == "M4_Core")

xgb_m12 <- comparison_table %>%
  dplyr::filter(Algorithm == "XGB", Subset == "M12_Full")

rf_m4 <- comparison_table %>%
  dplyr::filter(Algorithm == "RF", Subset == "M4_Core")

recommended_subset <- "M4_Core"

recommendation_reason <- paste0(
  "M4_Core was retained as the final mechanistically streamlined model. ",
  "It preserves predictors representing chemical availability, membrane partitioning, ",
  "organismal size, and feeding-related exposure while being directly compared against ",
  "the SHAP-informed M8_Extended model and the read-in M12_Full model. ",
  "Both XGB and RF prediction bundles for M4_Core were saved for downstream prediction."
)

recommendation_df <- data.frame(
  Recommended_Subset = recommended_subset,
  XGB_M4_CV_R2 = ifelse(nrow(xgb_m4) == 1, xgb_m4$CV_R2, NA),
  XGB_M4_CV_R2_Traditional = ifelse(nrow(xgb_m4) == 1, xgb_m4$CV_R2_Traditional, NA),
  XGB_M4_Gap_R2 = ifelse(nrow(xgb_m4) == 1, xgb_m4$Gap_R2, NA),
  XGB_M12_CV_R2 = ifelse(nrow(xgb_m12) == 1, xgb_m12$CV_R2, NA),
  Delta_XGB_M4_vs_M12 = ifelse(nrow(xgb_m4) == 1, xgb_m4$Delta_CV_R2_vs_M12, NA),
  RF_M4_CV_R2 = ifelse(nrow(rf_m4) == 1, rf_m4$CV_R2, NA),
  Reason = recommendation_reason,
  stringsAsFactors = FALSE
)

write.csv(
  recommendation_df,
  OUT_FILE("recommended_streamlined_model"),
  row.names = FALSE
)

cat("\n>>> Recommended streamlined subset:\n")
print(recommendation_df)


# 17. Prediction example script and notes ---------------------------------------

prediction_note <- data.frame(
  Bundle = c(
    "m4_core_xgb_prediction_bundle.rds",
    "m4_core_rf_prediction_bundle.rds"
  ),
  Required_predictors = paste(feature_sets$M4_Core, collapse = "; "),
  Response_scale = "Logtox",
  Note = c(
    "Load the RDS bundle and call bundle$predict_fun(bundle, newdata). Requires package xgboost.",
    "Load the RDS bundle and call bundle$predict_fun(bundle, newdata). Requires package caret/ranger."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  prediction_note,
  OUT_FILE("m4_core_prediction_bundle_notes"),
  row.names = FALSE
)

prediction_example <- c(
  "# ==============================================================================",
  "# Example: predicting Logtox using the final M4_Core prediction bundles",
  "# ============================================================================== ",
  "",
  "library(xgboost)",
  "library(caret)",
  "library(ranger)",
  "",
  "# New data must contain the four M4_Core predictors:",
  "# Solubility, Size, log.Kow, Feeding",
  "",
  "newdata <- read.csv('newdata_for_prediction.csv')",
  "",
  "xgb_bundle <- readRDS('m4_core_xgb_prediction_bundle.rds')",
  "rf_bundle  <- readRDS('m4_core_rf_prediction_bundle.rds')",
  "",
  "pred_logtox_xgb <- xgb_bundle$predict_fun(xgb_bundle, newdata)",
  "pred_logtox_rf  <- rf_bundle$predict_fun(rf_bundle, newdata)",
  "",
  "prediction_output <- data.frame(",
  "  pred_logtox_xgb = pred_logtox_xgb,",
  "  pred_logtox_rf = pred_logtox_rf,",
  "  pred_toxicity_xgb_ugL = 10^pred_logtox_xgb,",
  "  pred_toxicity_rf_ugL = 10^pred_logtox_rf",
  ")",
  "",
  "write.csv(prediction_output, 'm4_core_prediction_output.csv', row.names = FALSE)"
)

writeLines(
  prediction_example,
  OUT_FILE("m4_core_prediction_example", ext = "R")
)


# 18. Session info and cleanup --------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  OUT_FILE("session_info", ext = "txt")
)

try(parallel::stopCluster(cl), silent = TRUE)
foreach::registerDoSEQ()

cat("\nDone.\n")
cat("All outputs saved in folder:\n")
cat(" -", output_dir, "\n\n")

cat("Saved main tables:\n")
cat(" -", OUT_FILE("feature_sets_m4_m8_m12"), "\n")
cat(" -", OUT_FILE("streamlining_summary_m4_m8_trained"), "\n")
cat(" -", OUT_FILE("m12_full_readin_summary"), "\n")
cat(" -", OUT_FILE("streamlining_summary_m4_m8_m12"), "\n")
cat(" -", OUT_FILE("streamlining_comparison_m4_m8_m12"), "\n")
cat(" -", OUT_FILE("recommended_streamlined_model"), "\n")
cat(" -", OUT_FILE("traditional_oof_r2_m4_m8"), "\n")
cat(" -", OUT_FILE("m4_core_prediction_bundle_notes"), "\n")
cat(" -", OUT_FILE("m4_core_prediction_example", ext = "R"), "\n")
cat(" -", OUT_FILE("session_info", ext = "txt"), "\n")

cat("\nSaved model objects and prediction bundles:\n")
for (subset_name in names(feature_sets)) {
  prefix <- subset_file_prefix(subset_name)
  cat(" -", OUT_FILE(paste0(prefix, "_xgb_model"), ext = "rds"), "\n")
  cat(" -", OUT_FILE(paste0(prefix, "_rf_model"), ext = "rds"), "\n")
  cat(" -", OUT_FILE(paste0(prefix, "_xgb_prediction_bundle"), ext = "rds"), "\n")
  cat(" -", OUT_FILE(paste0(prefix, "_rf_prediction_bundle"), ext = "rds"), "\n")
}