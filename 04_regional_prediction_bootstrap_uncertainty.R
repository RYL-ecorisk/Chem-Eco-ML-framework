# ==============================================================================
# 04_regional_prediction_bootstrap_uncertainty.R
# ------------------------------------------------------------------------------
# Final 4-variable XGBoost regional prediction and bootstrap uncertainty.
# Computational settings are kept unchanged:
#   - Final M4_Core XGB bundle
#   - B = 1000 bootstrap refits
#   - SEED_MAIN = 2026
#   - CV_RMSE fallback = 0.707291251
# ==============================================================================

# ==============================================================================
# Final XGB 4-variable prediction + bootstrap uncertainty
# England / Jiangsu
#
# Updated version:
#   1. Use app-compatible model bundle:
#        m4_core_xgb_prediction_bundle.rds
#   2. Read CV_RMSE_4VARS from bundle if available;
#        otherwise use fallback value from final M4_Core result: 0.707291251
#   3. Update ssdtools syntax:
#        ssd_gof(fits, wt = TRUE)
#        ssd_hc(fits, proportion = 0.05, average = FALSE)
# ==============================================================================

library(readr)
library(dplyr)
library(xgboost)

# 0. Working directory ----------------------------------------------------------
# Set the working directory to the project root before running this script.
# 1. File names -----------------------------------------------------------------

bundle_file  <- "m4_core_xgb_prediction_bundle.rds"

train_file   <- "datasets.csv"
eng_file     <- "England_PAH_input_4vars.csv"
js_file      <- "Jiangsu_PAH_input_4vars.csv"

eng_out_file <- "England_PAH_pred_with_uncertainty_XGB4.csv"
js_out_file  <- "Jiangsu_PAH_pred_with_uncertainty_XGB4.csv"

# 2. Settings -------------------------------------------------------------------

# Fallback value only used when CV RMSE is not stored in the bundle
# Updated according to final M4_Core result:
# CV_RMSE = 0.707291251
CV_RMSE_FALLBACK <- 0.707291251

# Bootstrap iterations
B <- 1000

# Random seed
SEED_MAIN <- 2026

# 3. Load final model bundle ----------------------------------------------------

bundle <- readRDS(bundle_file)

# Safety check: this must be the full prediction bundle, not model-only RDS
required_bundle_items <- c("model", "features", "final_params", "final_nrounds")
missing_items <- setdiff(required_bundle_items, names(bundle))

if (length(missing_items) > 0) {
  stop(
    "The prediction bundle is missing required item(s): ",
    paste(missing_items, collapse = ", "),
    "\nPlease check whether this is the full prediction bundle rather than model-only RDS."
  )
}

features <- bundle$features
response_name <- "Logtox"

# 3.1 Read CV RMSE from bundle if available ------------------------------------

get_cv_rmse_from_bundle <- function(bundle, fallback_value) {
  
  candidate_list <- list(
    cv_rmse = bundle[["cv_rmse"]],
    CV_RMSE = bundle[["CV_RMSE"]],
    CV_RMSE_4VARS = bundle[["CV_RMSE_4VARS"]],
    cv_rmse_4vars = bundle[["cv_rmse_4vars"]],
    bundle_cv_rmse = bundle[["bundle.cv_rmse"]],
    
    metrics_CV_RMSE = if (!is.null(bundle[["metrics"]])) {
      bundle[["metrics"]][["CV_RMSE"]]
    } else {
      NULL
    },
    
    performance_CV_RMSE = if (!is.null(bundle[["performance"]])) {
      bundle[["performance"]][["CV_RMSE"]]
    } else {
      NULL
    },
    
    model_performance_CV_RMSE = if (!is.null(bundle[["model_performance"]])) {
      bundle[["model_performance"]][["CV_RMSE"]]
    } else {
      NULL
    },
    
    final_metrics_CV_RMSE = if (!is.null(bundle[["final_metrics"]])) {
      bundle[["final_metrics"]][["CV_RMSE"]]
    } else {
      NULL
    }
  )
  
  for (nm in names(candidate_list)) {
    v <- candidate_list[[nm]]
    
    if (!is.null(v)) {
      v_num <- suppressWarnings(as.numeric(v[1]))
      
      if (length(v_num) == 1 && is.finite(v_num) && v_num > 0) {
        cat("CV RMSE loaded from bundle field:", nm, "\n")
        return(v_num)
      }
    }
  }
  
  warning(
    "CV RMSE was not found in the bundle. ",
    "Using fallback value: ", fallback_value
  )
  
  return(fallback_value)
}

CV_RMSE_4VARS <- get_cv_rmse_from_bundle(
  bundle = bundle,
  fallback_value = CV_RMSE_FALLBACK
)

cat("Final model bundle loaded:\n")
cat(" -", bundle_file, "\n")
cat("Final model features:\n")
print(features)
cat("Final nrounds:", bundle$final_nrounds, "\n")
cat("CV_RMSE_4VARS used for uncertainty propagation:", CV_RMSE_4VARS, "\n")

# 4. Load training data ---------------------------------------------------------

df_train <- read.csv(train_file, stringsAsFactors = FALSE)

if (!all(features %in% names(df_train))) {
  stop(
    "Training data missing required columns: ",
    paste(setdiff(features, names(df_train)), collapse = ", ")
  )
}

if (!response_name %in% names(df_train)) {
  stop("Training data missing response column: ", response_name)
}

X_train <- as.matrix(df_train[, features])
y_train <- as.numeric(df_train[[response_name]])

if (any(!is.finite(X_train))) {
  stop("Training feature matrix contains NA/Inf.")
}

if (any(!is.finite(y_train))) {
  stop("Training response contains NA/Inf.")
}

cat("Training set loaded: n =", nrow(df_train), "\n")
cat("Using CV_RMSE as residual predictive error proxy:", CV_RMSE_4VARS, "\n")

# 5. Load prediction datasets ---------------------------------------------------

eng <- read_csv(eng_file, show_col_types = FALSE)
js  <- read_csv(js_file,  show_col_types = FALSE)

if (!all(features %in% names(eng))) {
  stop(
    "England dataset missing required columns: ",
    paste(setdiff(features, names(eng)), collapse = ", ")
  )
}

if (!all(features %in% names(js))) {
  stop(
    "Jiangsu dataset missing required columns: ",
    paste(setdiff(features, names(js)), collapse = ", ")
  )
}

X_eng <- as.matrix(eng[, features])
X_js  <- as.matrix(js[, features])

if (any(!is.finite(X_eng))) {
  stop("England feature matrix contains NA/Inf.")
}

if (any(!is.finite(X_js))) {
  stop("Jiangsu feature matrix contains NA/Inf.")
}

cat("England rows :", nrow(eng), "\n")
cat("Jiangsu rows :", nrow(js), "\n")

# 6. Best-model point prediction -----------------------------------------------

eng$Pred_LogTox_Best <- predict(bundle$model, X_eng)
js$Pred_LogTox_Best  <- predict(bundle$model, X_js)

eng$Pred_Tox_ugL_Best <- 10^(eng$Pred_LogTox_Best)
js$Pred_Tox_ugL_Best  <- 10^(js$Pred_LogTox_Best)

# 7. Bootstrap refits -----------------------------------------------------------

set.seed(SEED_MAIN)

pred_mat_eng <- matrix(NA_real_, nrow = nrow(eng), ncol = B)
pred_mat_js  <- matrix(NA_real_, nrow = nrow(js),  ncol = B)

cat("\nStarting bootstrap refits (B =", B, ") ...\n")

for (b in seq_len(B)) {
  
  idx <- sample(
    seq_len(nrow(df_train)),
    size = nrow(df_train),
    replace = TRUE
  )
  
  dboot <- xgb.DMatrix(
    data  = as.matrix(df_train[idx, features]),
    label = as.numeric(df_train[idx, response_name])
  )
  
  model_b <- xgb.train(
    params  = bundle$final_params,
    data    = dboot,
    nrounds = bundle$final_nrounds,
    verbose = 0
  )
  
  pred_mat_eng[, b] <- predict(model_b, X_eng)
  pred_mat_js[, b]  <- predict(model_b, X_js)
  
  if (b %% 50 == 0 || b == B) {
    cat("  -> bootstrap", b, "/", B, "done\n")
  }
}

# 8. Summarise bootstrap distributions -----------------------------------------

summarise_boot_preds <- function(pred_mat, cv_rmse) {
  
  pred_mean  <- rowMeans(pred_mat, na.rm = TRUE)
  sigma_boot <- apply(pred_mat, 1, sd, na.rm = TRUE)
  
  out <- data.frame(
    Pred_LogTox_Mean   = pred_mean,
    Sigma_Boot         = sigma_boot,
    Sigma_Total        = sqrt(sigma_boot^2 + cv_rmse^2),
    Pred_LogTox_Q2.5   = apply(pred_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    Pred_LogTox_Q50    = apply(pred_mat, 1, quantile, probs = 0.5,   na.rm = TRUE),
    Pred_LogTox_Q97.5  = apply(pred_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out$Pred_Tox_ugL_Mean  <- 10^(out$Pred_LogTox_Mean)
  out$Pred_Tox_ugL_Q2.5  <- 10^(out$Pred_LogTox_Q2.5)
  out$Pred_Tox_ugL_Q50   <- 10^(out$Pred_LogTox_Q50)
  out$Pred_Tox_ugL_Q97.5 <- 10^(out$Pred_LogTox_Q97.5)
  
  return(out)
}

eng_sum <- summarise_boot_preds(pred_mat_eng, CV_RMSE_4VARS)
js_sum  <- summarise_boot_preds(pred_mat_js,  CV_RMSE_4VARS)

# 9. Combine outputs ------------------------------------------------------------

eng_out <- bind_cols(eng, eng_sum)
js_out  <- bind_cols(js,  js_sum)

# 10. Consistency check ---------------------------------------------------------

r2_eng <- cor(eng_out$Pred_LogTox_Best, eng_out$Pred_LogTox_Mean)^2
r2_js  <- cor(js_out$Pred_LogTox_Best,  js_out$Pred_LogTox_Mean)^2

rmse_eng <- sqrt(mean((eng_out$Pred_LogTox_Best - eng_out$Pred_LogTox_Mean)^2))
rmse_js  <- sqrt(mean((js_out$Pred_LogTox_Best  - js_out$Pred_LogTox_Mean)^2))

cat("\nConsistency check:\n")
cat(
  "England | R2(best vs mean) =", round(r2_eng, 6),
  "| RMSE =", round(rmse_eng, 6), "\n"
)
cat(
  "Jiangsu | R2(best vs mean) =", round(r2_js, 6),
  "| RMSE =", round(rmse_js, 6), "\n"
)

# 11. Save prediction outputs ---------------------------------------------------

write_csv(eng_out, eng_out_file)
write_csv(js_out,  js_out_file)

cat("\nPrediction + uncertainty step done.\n")
cat("Saved files:\n")
cat(" -", eng_out_file, "\n")
cat(" -", js_out_file, "\n")









