# =============================================================================
# 01_full_model_training_comparison.R
# Four-model Bayesian-optimization workflow
# Elastic Net / SVM-RBF / Random Forest / XGBoost
# =============================================================================

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(readr)
  library(glmnet)
  library(kernlab)
  library(ranger)
  library(xgboost)
  library(rBayesianOptimization)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

DATA_FILE <- "datasets.csv"
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

features <- c(
  "Solubility", "log.Kow", "HaCount", "MW", "Complexity", "Rings",
  "Taxavalue", "Respiration", "Locomotion", "Feeding", "Size", "HLC"
)

cfg <- list(
  folds = 10L,
  repeats = 10L,
  gap_limit = 0.20,

  enet_alpha = c(0, 1),
  enet_lambda = c(1e-7, 1e-1),

  svm_sigma = c(2^-12, 2^-1),
  svm_C = c(2^-4, 2^8),

  rf_mtry = c(1, length(features)),
  rf_min_node_size = c(1, 20),
  rf_num_trees = 2200L,

  # Final expanded XGBoost search ranges
  xgb_max_depth = 5L,
  xgb_min_child_weight = 2,
  xgb_eta = c(0.005, 0.20),
  xgb_gamma = c(0, 1),
  xgb_colsample_bytree = c(0.30, 1.00),
  xgb_subsample = c(0.60, 1.00),
  xgb_lambda = c(1e-4, 2.00),
  xgb_alpha = c(0, 0.50),
  xgb_max_rounds = 5000L,
  xgb_early_stopping = 200L,

  bo_init = 20L,
  bo_iter = 30L
)

dat <- read_csv(DATA_FILE, show_col_types = FALSE)
stopifnot(all(c(features, "Logtox") %in% names(dat)))

X <- dat[, features]
y <- dat$Logtox
X_mat <- as.matrix(X)

cv_index <- createMultiFolds(y, k = cfg$folds, times = cfg$repeats)

ctrl <- trainControl(
  method = "repeatedcv",
  number = cfg$folds,
  repeats = cfg$repeats,
  index = cv_index,
  savePredictions = "final",
  allowParallel = TRUE
)

r2 <- function(obs, pred) suppressWarnings(cor(obs, pred, use = "complete.obs")^2)
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

model_summary <- function(model, name) {
  train_pred <- predict(model, X)
  data.frame(
    Model = name,
    CV_R2 = mean(model$resample$Rsquared, na.rm = TRUE),
    CV_RMSE = mean(model$resample$RMSE, na.rm = TRUE),
    CV_MAE = mean(model$resample$MAE, na.rm = TRUE),
    Train_R2 = r2(y, train_pred),
    Train_RMSE = rmse(y, train_pred),
    Train_MAE = mae(y, train_pred)
  ) |>
    mutate(Gap_R2 = Train_R2 - CV_R2)
}

select_stable <- function(history) {
  ok <- history |> filter(is.finite(CV_R2), Gap_R2 <= cfg$gap_limit)
  if (nrow(ok) == 0) ok <- history |> filter(is.finite(CV_R2))
  ok |> arrange(desc(CV_R2), CV_RMSE, Gap_R2) |> slice(1)
}

# -----------------------------------------------------------------------------
# Elastic Net
# -----------------------------------------------------------------------------

enet_history <- tibble()

enet_obj <- function(alpha, log10_lambda) {
  lambda <- 10^log10_lambda
  fit <- train(
    x = X, y = y, method = "glmnet", trControl = ctrl,
    tuneGrid = data.frame(alpha = alpha, lambda = lambda),
    metric = "Rsquared"
  )
  s <- model_summary(fit, "ElasticNet")
  enet_history <<- bind_rows(
    enet_history,
    tibble(alpha, lambda, CV_R2 = s$CV_R2, CV_RMSE = s$CV_RMSE,
           CV_MAE = s$CV_MAE, Train_R2 = s$Train_R2, Gap_R2 = s$Gap_R2)
  )
  list(Score = s$CV_R2, Pred = 0)
}

set.seed(999)
BayesianOptimization(
  enet_obj,
  bounds = list(alpha = cfg$enet_alpha, log10_lambda = log10(cfg$enet_lambda)),
  init_points = 12, n_iter = 20, acq = "ucb", verbose = FALSE
)

best_enet <- select_stable(enet_history)

model_enet <- train(
  x = X, y = y, method = "glmnet", trControl = ctrl,
  tuneGrid = data.frame(alpha = best_enet$alpha, lambda = best_enet$lambda),
  metric = "Rsquared"
)

# -----------------------------------------------------------------------------
# SVM-RBF
# -----------------------------------------------------------------------------

svm_history <- tibble()

svm_obj <- function(log2_sigma, log2_C) {
  sigma <- 2^log2_sigma
  C <- 2^log2_C
  fit <- train(
    x = X, y = y, method = "svmRadial", trControl = ctrl,
    preProcess = c("center", "scale"),
    tuneGrid = data.frame(sigma = sigma, C = C),
    metric = "Rsquared"
  )
  s <- model_summary(fit, "SVM")
  svm_history <<- bind_rows(
    svm_history,
    tibble(sigma, C, CV_R2 = s$CV_R2, CV_RMSE = s$CV_RMSE,
           CV_MAE = s$CV_MAE, Train_R2 = s$Train_R2, Gap_R2 = s$Gap_R2)
  )
  list(Score = s$CV_R2, Pred = 0)
}

set.seed(999)
BayesianOptimization(
  svm_obj,
  bounds = list(log2_sigma = log2(cfg$svm_sigma), log2_C = log2(cfg$svm_C)),
  init_points = 15, n_iter = 25, acq = "ucb", verbose = FALSE
)

best_svm <- select_stable(svm_history)

model_svm <- train(
  x = X, y = y, method = "svmRadial", trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneGrid = data.frame(sigma = best_svm$sigma, C = best_svm$C),
  metric = "Rsquared"
)

# -----------------------------------------------------------------------------
# Random Forest
# -----------------------------------------------------------------------------

rf_history <- tibble()

rf_obj <- function(mtry, min_node_size) {
  mtry <- as.integer(round(mtry))
  min_node_size <- as.integer(round(min_node_size))

  fit <- train(
    x = X, y = y, method = "ranger", trControl = ctrl,
    num.trees = cfg$rf_num_trees,
    importance = "permutation",
    tuneGrid = data.frame(
      mtry = mtry,
      splitrule = "variance",
      min.node.size = min_node_size
    ),
    metric = "Rsquared"
  )

  s <- model_summary(fit, "RandomForest")
  rf_history <<- bind_rows(
    rf_history,
    tibble(mtry, min.node.size = min_node_size,
           CV_R2 = s$CV_R2, CV_RMSE = s$CV_RMSE, CV_MAE = s$CV_MAE,
           Train_R2 = s$Train_R2, Gap_R2 = s$Gap_R2)
  )
  list(Score = s$CV_R2, Pred = 0)
}

set.seed(999)
BayesianOptimization(
  rf_obj,
  bounds = list(mtry = cfg$rf_mtry, min_node_size = cfg$rf_min_node_size),
  init_points = 12, n_iter = 20, acq = "ucb", verbose = FALSE
)

best_rf <- select_stable(rf_history)

model_rf <- train(
  x = X, y = y, method = "ranger", trControl = ctrl,
  num.trees = cfg$rf_num_trees,
  importance = "permutation",
  tuneGrid = data.frame(
    mtry = as.integer(best_rf$mtry),
    splitrule = "variance",
    min.node.size = as.integer(best_rf$min.node.size)
  ),
  metric = "Rsquared"
)

# -----------------------------------------------------------------------------
# XGBoost
# -----------------------------------------------------------------------------

make_xgb_params <- function(gamma, eta, colsample_bytree, subsample, lambda, alpha) {
  list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = cfg$xgb_max_depth,
    min_child_weight = cfg$xgb_min_child_weight,
    gamma = gamma,
    eta = eta,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    lambda = lambda,
    alpha = alpha,
    nthread = 1
  )
}

xgb_eval <- function(params) {
  fold_metrics <- vector("list", length(cv_index))
  best_rounds <- integer(length(cv_index))

  for (i in seq_along(cv_index)) {
    train_idx <- cv_index[[i]]
    test_idx <- setdiff(seq_along(y), train_idx)

    set.seed(5000 + i)
    valid_idx <- sample(train_idx, size = max(5, floor(length(train_idx) * 0.20)))
    inner_train <- setdiff(train_idx, valid_idx)

    dtr <- xgb.DMatrix(X_mat[inner_train, , drop = FALSE], label = y[inner_train])
    dva <- xgb.DMatrix(X_mat[valid_idx, , drop = FALSE], label = y[valid_idx])

    m0 <- xgb.train(
      params = params, data = dtr, nrounds = cfg$xgb_max_rounds,
      evals = list(valid = dva), early_stopping_rounds = cfg$xgb_early_stopping,
      verbose = 0
    )

    br <- as.integer(attr(m0, "best_iteration"))
    if (!is.finite(br)) br <- cfg$xgb_max_rounds
    best_rounds[i] <- br

    dfull <- xgb.DMatrix(X_mat[train_idx, , drop = FALSE], label = y[train_idx])
    m <- xgb.train(params = params, data = dfull, nrounds = br, verbose = 0)
    pr <- predict(m, X_mat[test_idx, , drop = FALSE])

    fold_metrics[[i]] <- data.frame(
      R2 = r2(y[test_idx], pr),
      RMSE = rmse(y[test_idx], pr),
      MAE = mae(y[test_idx], pr)
    )
  }

  fm <- bind_rows(fold_metrics)
  nrounds <- as.integer(round(median(best_rounds)))

  dfull <- xgb.DMatrix(X_mat, label = y)
  model <- xgb.train(params = params, data = dfull, nrounds = nrounds, verbose = 0)
  train_pred <- predict(model, X_mat)

  list(
    CV_R2 = mean(fm$R2, na.rm = TRUE),
    CV_RMSE = mean(fm$RMSE, na.rm = TRUE),
    CV_MAE = mean(fm$MAE, na.rm = TRUE),
    Train_R2 = r2(y, train_pred),
    Train_RMSE = rmse(y, train_pred),
    Train_MAE = mae(y, train_pred),
    Gap_R2 = r2(y, train_pred) - mean(fm$R2, na.rm = TRUE),
    nrounds = nrounds,
    model = model
  )
}

xgb_history <- tibble()

xgb_obj <- function(gamma, log10_eta, colsample_bytree, subsample,
                    log10_lambda, alpha) {
  params <- make_xgb_params(
    gamma = gamma, eta = 10^log10_eta,
    colsample_bytree = colsample_bytree, subsample = subsample,
    lambda = 10^log10_lambda, alpha = alpha
  )
  z <- xgb_eval(params)

  xgb_history <<- bind_rows(
    xgb_history,
    tibble(
      gamma, eta = 10^log10_eta, colsample_bytree, subsample,
      lambda = 10^log10_lambda, alpha,
      CV_R2 = z$CV_R2, CV_RMSE = z$CV_RMSE, CV_MAE = z$CV_MAE,
      Train_R2 = z$Train_R2, Gap_R2 = z$Gap_R2, nrounds = z$nrounds
    )
  )
  list(Score = z$CV_R2, Pred = 0)
}

set.seed(2026)
BayesianOptimization(
  xgb_obj,
  bounds = list(
    gamma = cfg$xgb_gamma,
    log10_eta = log10(cfg$xgb_eta),
    colsample_bytree = cfg$xgb_colsample_bytree,
    subsample = cfg$xgb_subsample,
    log10_lambda = log10(cfg$xgb_lambda),
    alpha = cfg$xgb_alpha
  ),
  init_points = cfg$bo_init, n_iter = cfg$bo_iter,
  acq = "ucb", verbose = FALSE
)

best_xgb <- select_stable(xgb_history)

final_xgb_params <- make_xgb_params(
  gamma = best_xgb$gamma,
  eta = best_xgb$eta,
  colsample_bytree = best_xgb$colsample_bytree,
  subsample = best_xgb$subsample,
  lambda = best_xgb$lambda,
  alpha = best_xgb$alpha
)

dall <- xgb.DMatrix(X_mat, label = y)
model_xgb <- xgb.train(
  params = final_xgb_params, data = dall,
  nrounds = as.integer(best_xgb$nrounds), verbose = 0
)

xgb_metrics <- data.frame(
  Model = "XGBoost",
  CV_R2 = best_xgb$CV_R2,
  CV_RMSE = best_xgb$CV_RMSE,
  CV_MAE = best_xgb$CV_MAE,
  Train_R2 = r2(y, predict(model_xgb, X_mat)),
  Train_RMSE = rmse(y, predict(model_xgb, X_mat)),
  Train_MAE = mae(y, predict(model_xgb, X_mat))
) |>
  mutate(Gap_R2 = Train_R2 - CV_R2)

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

metrics <- bind_rows(
  model_summary(model_enet, "ElasticNet"),
  model_summary(model_svm, "SVM-RBF"),
  model_summary(model_rf, "RandomForest"),
  xgb_metrics
) |>
  arrange(desc(CV_R2), CV_RMSE)

write_csv(metrics, file.path(OUT_DIR, "Model_Performance_Summary_FINAL.csv"))

# Pooled out-of-fold predictions for the three caret models.
pooled_caret <- function(model, name) {
  model$pred |>
    group_by(rowIndex) |>
    summarise(
      Observed_Log10 = mean(obs, na.rm = TRUE),
      Predicted_Log10 = mean(pred, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(Model = name, Set = "Cross-validation")
}

cv_plot <- bind_rows(
  pooled_caret(model_enet, "ElasticNet"),
  pooled_caret(model_svm, "SVM-RBF"),
  pooled_caret(model_rf, "RandomForest")
)

# XGBoost pooled OOF predictions using final parameters and fixed nrounds.
xgb_oof <- lapply(seq_along(cv_index), function(i) {
  tr <- cv_index[[i]]
  te <- setdiff(seq_along(y), tr)
  m <- xgb.train(
    final_xgb_params,
    xgb.DMatrix(X_mat[tr, , drop = FALSE], label = y[tr]),
    nrounds = as.integer(best_xgb$nrounds), verbose = 0
  )
  data.frame(rowIndex = te, obs = y[te], pred = predict(m, X_mat[te, , drop = FALSE]))
}) |>
  bind_rows() |>
  group_by(rowIndex) |>
  summarise(
    Observed_Log10 = mean(obs),
    Predicted_Log10 = mean(pred),
    .groups = "drop"
  ) |>
  mutate(Model = "XGBoost", Set = "Cross-validation")

cv_plot <- bind_rows(cv_plot, xgb_oof)

train_plot <- bind_rows(
  tibble(Model = "ElasticNet", Observed_Log10 = y,
         Predicted_Log10 = predict(model_enet, X), Set = "Training"),
  tibble(Model = "SVM-RBF", Observed_Log10 = y,
         Predicted_Log10 = predict(model_svm, X), Set = "Training"),
  tibble(Model = "RandomForest", Observed_Log10 = y,
         Predicted_Log10 = predict(model_rf, X), Set = "Training"),
  tibble(Model = "XGBoost", Observed_Log10 = y,
         Predicted_Log10 = predict(model_xgb, X_mat), Set = "Training")
)

figure3_data <- bind_rows(train_plot, cv_plot) |>
  mutate(
    Observed_ug_L = 10^Observed_Log10,
    Predicted_ug_L = 10^Predicted_Log10
  )

write_csv(figure3_data, file.path(OUT_DIR, "Figure3_Scatter_Plot_Data_FINAL.csv"))

p3 <- ggplot(
  figure3_data,
  aes(Observed_ug_L, Predicted_ug_L, colour = Set)
) +
  geom_point(alpha = 0.55, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~Model, ncol = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = expression(paste("Observed toxicity (", mu, "g/L)")),
    y = expression(paste("Predicted toxicity (", mu, "g/L)")),
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "Figure3_Main_Model_Comparison.png"),
       p3, width = 8.2, height = 7.2, dpi = 600)

bundle <- list(
  data = dat,
  X = X,
  y = y,
  features = features,
  cv_index = cv_index,
  models = list(
    ElasticNet = model_enet,
    SVM_RBF = model_svm,
    RandomForest = model_rf,
    XGBoost = model_xgb
  ),
  xgb_params = final_xgb_params,
  xgb_nrounds = as.integer(best_xgb$nrounds),
  model_performance = metrics,
  configuration = cfg
)

saveRDS(bundle, file.path(OUT_DIR, "All_Final_Models_Analysis_Bundle_FINAL.rds"))

message("01 complete.")
