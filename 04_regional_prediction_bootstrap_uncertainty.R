# =============================================================================
# 04_regional_prediction_bootstrap_uncertainty.R
# Regional prediction and bootstrap uncertainty for England and Jiangsu
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(xgboost)
})

BUNDLE_FILE <- "outputs/Four4_Optimal_Model_Bundle.rds"
ENG_FILE <- "England_PAH_input_4vars.csv"
JS_FILE <- "Jiangsu_PAH_input_4vars.csv"
OUT_DIR <- "outputs"

B <- 1000L
SEED <- 2026L

bundle <- readRDS(BUNDLE_FILE)

features <- bundle$features
train <- bundle$training_data
response <- bundle$response_name
cv_rmse <- bundle$CV_RMSE

eng <- read_csv(ENG_FILE, show_col_types = FALSE)
js  <- read_csv(JS_FILE, show_col_types = FALSE)

stopifnot(all(features %in% names(eng)))
stopifnot(all(features %in% names(js)))

X_train <- as.matrix(train[, features, drop = FALSE])
y_train <- train[[response]]
X_eng <- as.matrix(eng[, features, drop = FALSE])
X_js  <- as.matrix(js[, features, drop = FALSE])

eng$Pred_LogTox_Best <- predict(bundle$model, X_eng)
js$Pred_LogTox_Best  <- predict(bundle$model, X_js)

set.seed(SEED)

pred_eng <- matrix(NA_real_, nrow(eng), B)
pred_js  <- matrix(NA_real_, nrow(js), B)

for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(train)), replace = TRUE)

  m <- xgb.train(
    params = bundle$final_params,
    data = xgb.DMatrix(X_train[idx, , drop = FALSE], label = y_train[idx]),
    nrounds = bundle$final_nrounds,
    verbose = 0
  )

  pred_eng[, b] <- predict(m, X_eng)
  pred_js[, b]  <- predict(m, X_js)
}

summarise_boot <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  s_boot <- apply(mat, 1, sd, na.rm = TRUE)

  tibble(
    Pred_LogTox_Mean = mu,
    Sigma_Boot = s_boot,
    Sigma_Total = sqrt(s_boot^2 + cv_rmse^2),
    Pred_LogTox_Q2.5 = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    Pred_LogTox_Q50 = apply(mat, 1, quantile, 0.50, na.rm = TRUE),
    Pred_LogTox_Q97.5 = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
}

eng_out <- bind_cols(eng, summarise_boot(pred_eng))
js_out  <- bind_cols(js, summarise_boot(pred_js))

write_csv(eng_out, file.path(OUT_DIR, "England_PAH_pred_with_uncertainty_Four4.csv"))
write_csv(js_out, file.path(OUT_DIR, "Jiangsu_PAH_pred_with_uncertainty_Four4.csv"))

message("04 complete.")
