# =============================================================================
# 02_model_interpretation_shap.R
# Exact native TreeSHAP interpretation of the final 12-variable XGBoost model
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
})

BUNDLE_FILE <- "outputs/All_Final_Models_Analysis_Bundle_FINAL.rds"
OUT_DIR <- "outputs"
TOP_N <- 6L

bundle <- readRDS(BUNDLE_FILE)

model <- bundle$models$XGBoost
X <- as.data.frame(bundle$X)
features <- bundle$features
X_mat <- as.matrix(X)

# Exact TreeSHAP from XGBoost. The final column is the baseline contribution.
contrib <- predict(
  model,
  xgb.DMatrix(X_mat),
  predcontrib = TRUE,
  approxcontrib = FALSE
) |>
  as.matrix()

shap <- contrib[, seq_along(features), drop = FALSE]
colnames(shap) <- features

importance <- tibble(
  Feature = features,
  MeanAbsSHAP = colMeans(abs(shap))
) |>
  mutate(Percentage = 100 * MeanAbsSHAP / sum(MeanAbsSHAP)) |>
  arrange(desc(MeanAbsSHAP))

write_csv(importance, file.path(OUT_DIR, "XGB_Native_TreeSHAP_Importance.csv"))

shap_long <- as.data.frame(shap) |>
  mutate(Row = row_number()) |>
  pivot_longer(-Row, names_to = "Feature", values_to = "SHAP") |>
  left_join(
    X |>
      mutate(Row = row_number()) |>
      pivot_longer(-Row, names_to = "Feature", values_to = "Value"),
    by = c("Row", "Feature")
  )

write_csv(shap_long, file.path(OUT_DIR, "XGB_Native_TreeSHAP_Long.csv"))

top_features <- head(importance$Feature, TOP_N)

labels <- c(
  Solubility = "Water solubility",
  `log.Kow` = "log KOW",
  HaCount = "Heavy atom count",
  MW = "Molecular weight",
  Complexity = "Molecular complexity",
  Rings = "Ring count",
  Taxavalue = "Taxonomic value",
  Respiration = "Respiration mode",
  Locomotion = "Locomotion mode",
  Feeding = "Feeding group",
  Size = "Maximum body size",
  HLC = "Henry's law constant"
)

label_feature <- function(x) {
  y <- unname(labels[x])
  ifelse(is.na(y), x, y)
}

# Global importance panel.
p_importance <- importance |>
  mutate(
    Display = factor(label_feature(Feature),
                     levels = rev(label_feature(Feature)))
  ) |>
  ggplot(aes(Percentage, Display)) +
  geom_col(width = 0.70) +
  labs(
    x = "Normalized mean |TreeSHAP| (%)",
    y = NULL,
    title = "(a) Global TreeSHAP importance"
  ) +
  theme_classic(base_size = 12)

# Dependence panels.
make_dependence <- function(f) {
  d <- shap_long |> filter(Feature == f)

  # Feeding is an encoded ecological category; show distributions by code.
  if (f == "Feeding") {
    return(
      ggplot(d, aes(factor(Value), SHAP)) +
        geom_boxplot(outlier.shape = NA, width = 0.55) +
        geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.4) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
        labs(x = label_feature(f), y = "TreeSHAP value") +
        theme_classic(base_size = 11)
    )
  }

  ggplot(d, aes(Value, SHAP)) +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_smooth(
      method = "lm",
      formula = y ~ poly(x, 2, raw = TRUE),
      se = TRUE,
      linewidth = 0.8
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    labs(x = label_feature(f), y = "TreeSHAP value") +
    theme_classic(base_size = 11)
}

dep_plots <- lapply(top_features, make_dependence)

figure4 <- wrap_plots(
  c(list(p_importance), dep_plots),
  ncol = 2,
  guides = "collect"
)

ggsave(
  file.path(OUT_DIR, "Figure4_Main_XGBoost_TreeSHAP.png"),
  figure4, width = 9.2, height = 11.5, dpi = 600
)

message("02 complete.")
