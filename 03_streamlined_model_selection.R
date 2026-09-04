# =============================================================================
# 03_streamlined_model_selection.R
# Streamline the full XGBoost model directly to the final Four4 model
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(xgboost)
  library(rBayesianOptimization)
  library(ggplot2)
})

FULL_BUNDLE <- "outputs/All_Final_Models_Analysis_Bundle_FINAL.rds"
OUT_DIR <- "outputs"

FOUR4 <- c("Solubility", "Size", "log.Kow", "Feeding")

cfg <- list(
  max_depth = 5L,
  min_child_weight = 2,
  eta = c(0.005, 0.20),
  gamma = c(0, 1),
  colsample_bytree = c(0.30, 1.00),
  subsample = c(0.60, 1.00),
  lambda = c(1e-4, 2.00),
  alpha = c(0, 0.50),
  max_rounds = 5000L,
  early_stopping = 200L,
  gap_limit = 0.20,
  init_points = 20L,
  n_iter = 30L
)

full <- readRDS(FULL_BUNDLE)
dat <- as.data.frame(full$data)
y <- full$y
cv_index <- full$cv_index

r2 <- function(obs, pred) suppressWarnings(cor(obs, pred, use = "complete.obs")^2)
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

make_params <- function(gamma, eta, colsample_bytree, subsample, lambda, alpha) {
  list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = cfg$max_depth,
    min_child_weight = cfg$min_child_weight,
    gamma = gamma,
    eta = eta,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    lambda = lambda,
    alpha = alpha,
    nthread = 1
  )
}

evaluate <- function(params, fixed_nrounds = NULL) {
  X <- as.matrix(dat[, FOUR4, drop = FALSE])
  fold <- vector("list", length(cv_index))
  rounds <- integer(length(cv_index))

  for (i in seq_along(cv_index)) {
    tr <- cv_index[[i]]
    te <- setdiff(seq_along(y), tr)

    if (is.null(fixed_nrounds)) {
      set.seed(5000 + i)
      va <- sample(tr, size = max(5, floor(length(tr) * 0.20)))
      it <- setdiff(tr, va)

      m0 <- xgb.train(
        params,
        xgb.DMatrix(X[it, , drop = FALSE], label = y[it]),
        nrounds = cfg$max_rounds,
        evals = list(valid = xgb.DMatrix(X[va, , drop = FALSE], label = y[va])),
        early_stopping_rounds = cfg$early_stopping,
        verbose = 0
      )

      nr <- as.integer(attr(m0, "best_iteration"))
      if (!is.finite(nr)) nr <- cfg$max_rounds
    } else {
      nr <- fixed_nrounds
    }

    rounds[i] <- nr

    m <- xgb.train(
      params,
      xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr]),
      nrounds = nr,
      verbose = 0
    )

    pr <- predict(m, X[te, , drop = FALSE])

    fold[[i]] <- data.frame(
      R2 = r2(y[te], pr),
      RMSE = rmse(y[te], pr),
      MAE = mae(y[te], pr)
    )
  }

  fm <- bind_rows(fold)
  nr <- if (is.null(fixed_nrounds)) as.integer(round(median(rounds))) else fixed_nrounds

  mfull <- xgb.train(
    params,
    xgb.DMatrix(X, label = y),
    nrounds = nr,
    verbose = 0
  )

  train_pred <- predict(mfull, X)

  list(
    CV_R2 = mean(fm$R2, na.rm = TRUE),
    CV_RMSE = mean(fm$RMSE, na.rm = TRUE),
    CV_MAE = mean(fm$MAE, na.rm = TRUE),
    Train_R2 = r2(y, train_pred),
    Gap_R2 = r2(y, train_pred) - mean(fm$R2, na.rm = TRUE),
    nrounds = nr,
    model = mfull
  )
}

history <- tibble()

objective <- function(gamma, log10_eta, colsample_bytree, subsample,
                      log10_lambda, alpha) {
  params <- make_params(
    gamma = gamma,
    eta = 10^log10_eta,
    colsample_bytree = colsample_bytree,
    subsample = subsample,
    lambda = 10^log10_lambda,
    alpha = alpha
  )

  z <- evaluate(params)

  history <<- bind_rows(
    history,
    tibble(
      gamma,
      eta = 10^log10_eta,
      colsample_bytree,
      subsample,
      lambda = 10^log10_lambda,
      alpha,
      CV_R2 = z$CV_R2,
      CV_RMSE = z$CV_RMSE,
      CV_MAE = z$CV_MAE,
      Train_R2 = z$Train_R2,
      Gap_R2 = z$Gap_R2,
      nrounds = z$nrounds
    )
  )

  list(Score = z$CV_R2, Pred = 0)
}

set.seed(22026)
BayesianOptimization(
  objective,
  bounds = list(
    gamma = cfg$gamma,
    log10_eta = log10(cfg$eta),
    colsample_bytree = cfg$colsample_bytree,
    subsample = cfg$subsample,
    log10_lambda = log10(cfg$lambda),
    alpha = cfg$alpha
  ),
  init_points = cfg$init_points,
  n_iter = cfg$n_iter,
  acq = "ucb",
  verbose = FALSE
)

stable <- history |> filter(Gap_R2 <= cfg$gap_limit)
if (nrow(stable) == 0) stable <- history

best <- stable |>
  arrange(desc(CV_R2), CV_RMSE, Gap_R2) |>
  slice(1)

final_params <- make_params(
  best$gamma,
  best$eta,
  best$colsample_bytree,
  best$subsample,
  best$lambda,
  best$alpha
)

final <- evaluate(final_params, fixed_nrounds = as.integer(best$nrounds))

bundle <- list(
  model_name = "Four4",
  model = final$model,
  features = FOUR4,
  final_params = final_params,
  final_nrounds = as.integer(final$nrounds),
  CV_RMSE = final$CV_RMSE,
  CV_R2 = final$CV_R2,
  CV_MAE = final$CV_MAE,
  training_data = dat[, c(FOUR4, "Logtox"), drop = FALSE],
  response_name = "Logtox",
  cv_index = cv_index
)

saveRDS(bundle, file.path(OUT_DIR, "Four4_Optimal_Model_Bundle.rds"))

full_row <- full$model_performance |>
  filter(grepl("XGBoost", Model, ignore.case = TRUE)) |>
  slice(1) |>
  transmute(
    Model = "Full12",
    N_predictors = 12,
    CV_R2, CV_RMSE, CV_MAE, Train_R2, Gap_R2
  )

four_row <- tibble(
  Model = "Four4",
  N_predictors = 4,
  CV_R2 = final$CV_R2,
  CV_RMSE = final$CV_RMSE,
  CV_MAE = final$CV_MAE,
  Train_R2 = final$Train_R2,
  Gap_R2 = final$Gap_R2
)

comparison <- bind_rows(full_row, four_row)
write_csv(comparison, file.path(OUT_DIR, "XGB_Full12_Four4_Comparison.csv"))

p <- comparison |>
  select(Model, CV_R2, CV_RMSE, CV_MAE) |>
  tidyr::pivot_longer(-Model, names_to = "Metric", values_to = "Value") |>
  ggplot(aes(Model, Value, fill = Model)) +
  geom_col(width = 0.65) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL)

ggsave(
  file.path(OUT_DIR, "Figure_Streamlined_Full12_Four4.png"),
  p, width = 7.2, height = 3.8, dpi = 600
)

message("03 complete. Final model: Four4.")
