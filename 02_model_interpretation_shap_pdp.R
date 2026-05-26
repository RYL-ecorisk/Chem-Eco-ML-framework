# ==============================================================================
# 02_model_interpretation_shap_pdp.R
# ------------------------------------------------------------------------------
# SHAP, PDP, and feeding-group interpretation.
# Run after 01_full_model_training_comparison.R in the same R session so that
# the selected model and training objects are available.
# Computational settings are kept from the working analysis script.
# ==============================================================================

# ==============================================================================
# 15. Export SHAP-ready bundle for independent interpretation
# ==============================================================================

cat("\n>>> Exporting SHAP-ready bundle...\n")

# Keep this section aligned with the training pipeline.
# Required objects from the training section:
# X_mat, features, y, best_model_name, caret_models,
# model_xgb_bayes, final_xgb_params, best_xgb_nrounds.

X_export <- as.data.frame(X_mat)
colnames(X_export) <- features

if (best_model_name == "XGB_Bayes") {
  
  shap_bundle <- list(
    best_model_name = best_model_name,
    model_object = model_xgb_bayes,
    X = X_export,
    y = y,
    features = features,
    xgb_params = final_xgb_params,
    xgb_nrounds = best_xgb_nrounds,
    model_class = "xgboost_native"
  )
  
} else {
  
  shap_bundle <- list(
    best_model_name = best_model_name,
    model_object = caret_models[[best_model_name]],
    X = X_export,
    y = y,
    features = features,
    model_class = "caret"
  )
}

saveRDS(
  shap_bundle,
  OUT_FILE("SHAP_Ready_Bundle", ext = "rds")
)

cat("Saved:", OUT_FILE("SHAP_Ready_Bundle", ext = "rds"), "\n")


# ==============================================================================
# 16. SHAP analysis
# ==============================================================================

cat("\n>>> Running SHAP analysis...\n")

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  dplyr,
  tidyr,
  ggplot2,
  ggbeeswarm,
  fastshap,
  foreach,
  doParallel,
  scales,
  grid,
  xgboost
)

if (!exists("RUN_TAG")) {
  RUN_TAG <- "FINAL"
}

if (!exists("OUT_FILE")) {
  OUT_FILE <- function(name, ext = "csv") {
    paste0(name, "_", RUN_TAG, ".", ext)
  }
}

bundle_file <- OUT_FILE("SHAP_Ready_Bundle", ext = "rds")

if (!file.exists(bundle_file)) {
  stop("Cannot find SHAP-ready bundle: ", bundle_file)
}

bundle <- readRDS(bundle_file)

best_model_name <- bundle$best_model_name
best_model_obj  <- bundle$model_object
X_df            <- as.data.frame(bundle$X)
y               <- bundle$y
features        <- bundle$features
model_class     <- bundle$model_class

cat("Loaded bundle:", bundle_file, "\n")
cat("Best model:", best_model_name, "\n")
cat("Model class:", model_class, "\n")
cat("Samples:", nrow(X_df), "| Features:", ncol(X_df), "\n")

# ------------------------------------------------------------------------------
# Optional SHAP subsampling
# ------------------------------------------------------------------------------

use_subsample <- FALSE
subsample_n <- min(500, nrow(X_df))

if (use_subsample) {
  set.seed(123)
  sample_indices <- sample(seq_len(nrow(X_df)), subsample_n)
  X_shap <- X_df[sample_indices, , drop = FALSE]
} else {
  X_shap <- X_df
}

# ------------------------------------------------------------------------------
# Feature display names
# ------------------------------------------------------------------------------

name_map <- c(
  "Solubility"  = "Water solubility",
  "log.Kow"     = "Log Kow",
  "HaCount"     = "Heavy atom count",
  "MW"          = "Molecular weight",
  "Complexity"  = "Molecular complexity",
  "Rings"       = "Ring count",
  "Taxavalue"   = "Taxonomic classification",
  "Respiration" = "Respiration mode",
  "Locomotion"  = "Locomotion mode",
  "Feeding"     = "Feeding habits",
  "Size"        = "Maximum body size",
  "HLC"         = "Henry's law constant"
)

# ------------------------------------------------------------------------------
# Compute SHAP values
# ------------------------------------------------------------------------------

raw_shap_matrix <- NULL

if (model_class == "xgboost_native") {
  
  cat("Using native XGBoost SHAP...\n")
  
  dshap <- xgboost::xgb.DMatrix(data = as.matrix(X_shap))
  
  shap_pred <- predict(
    best_model_obj,
    newdata = dshap,
    predcontrib = TRUE
  )
  
  shap_pred <- as.data.frame(shap_pred)
  
  bias_col_idx <- which(tolower(colnames(shap_pred)) %in% c("biasterm", "bias"))
  
  if (length(bias_col_idx) > 0) {
    raw_shap_matrix <- shap_pred[, -bias_col_idx, drop = FALSE]
  } else {
    if (ncol(shap_pred) == (ncol(X_shap) + 1)) {
      raw_shap_matrix <- shap_pred[, seq_len(ncol(X_shap)), drop = FALSE]
    } else {
      raw_shap_matrix <- shap_pred
    }
  }
  
  raw_shap_matrix <- raw_shap_matrix[, colnames(X_shap), drop = FALSE]
  
} else if (model_class == "caret") {
  
  cat("Using fastshap for caret model...\n")
  
  cl_shap <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
  doParallel::registerDoParallel(cl_shap)
  
  pred_wrapper <- function(object, newdata) {
    predict(object, newdata = as.data.frame(newdata))
  }
  
  raw_shap_matrix <- fastshap::explain(
    object = best_model_obj,
    X = X_shap,
    pred_wrapper = pred_wrapper,
    nsim = 200
  )
  
  parallel::stopCluster(cl_shap)
  foreach::registerDoSEQ()
  
  raw_shap_matrix <- as.data.frame(raw_shap_matrix)
  raw_shap_matrix <- raw_shap_matrix[, colnames(X_shap), drop = FALSE]
  
} else {
  stop("Unknown model_class: ", model_class)
}

cat("SHAP matrix created:", nrow(raw_shap_matrix), "x", ncol(raw_shap_matrix), "\n")

# ------------------------------------------------------------------------------
# SHAP long-format data
# ------------------------------------------------------------------------------

plot_data <- as.data.frame(raw_shap_matrix)
plot_data$id <- seq_len(nrow(plot_data))

shap_long <- plot_data %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "Feature",
    values_to = "SHAP_Value"
  )

features_long <- X_shap %>%
  dplyr::mutate(id = seq_len(nrow(X_shap))) %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "Feature",
    values_to = "Feature_Value"
  )

plot_data_long <- dplyr::left_join(
  shap_long,
  features_long,
  by = c("id", "Feature")
)

plot_data_renamed <- plot_data_long %>%
  dplyr::mutate(
    Feature = dplyr::recode(Feature, !!!name_map)
  ) %>%
  dplyr::group_by(Feature) %>%
  dplyr::mutate(
    min_val_feat = min(Feature_Value, na.rm = TRUE),
    max_val_feat = max(Feature_Value, na.rm = TRUE),
    normalized_per_feature = dplyr::if_else(
      (max_val_feat - min_val_feat) == 0,
      0.5,
      (Feature_Value - min_val_feat) / (max_val_feat - min_val_feat)
    )
  ) %>%
  dplyr::ungroup()

feature_rank <- plot_data_renamed %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    mean_abs_shap = mean(abs(SHAP_Value), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_abs_shap)

feature_levels <- feature_rank$Feature

plot_data_renamed$Feature <- factor(
  plot_data_renamed$Feature,
  levels = feature_levels
)

imp_labels <- feature_rank %>%
  dplyr::mutate(
    label_text = sprintf("%.3f", mean_abs_shap),
    Feature = factor(Feature, levels = feature_levels)
  )

# ------------------------------------------------------------------------------
# SHAP summary plot
# ------------------------------------------------------------------------------

data_min <- min(plot_data_renamed$SHAP_Value, na.rm = TRUE)
data_max <- max(plot_data_renamed$SHAP_Value, na.rm = TRUE)
range_span <- data_max - data_min

if (!is.finite(range_span) || range_span == 0) {
  range_span <- 1
}

text_position_x <- data_max + (range_span * 0.02)

pub_colors <- c(
  "#4575b4", "#74add1", "#abd9e9",
  "#fee090", "#fdae61", "#f46d43", "#d73027"
)

p_summary_final <- ggplot(
  plot_data_renamed,
  aes(x = SHAP_Value, y = Feature, color = normalized_per_feature)
) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.5) +
  ggbeeswarm::geom_quasirandom(
    groupOnX = FALSE,
    varwidth = TRUE,
    bandwidth = 0.2,
    size = 2.5,
    alpha = 0.8,
    stroke = 0.1
  ) +
  geom_text(
    data = imp_labels,
    aes(x = text_position_x, y = Feature, label = label_text),
    inherit.aes = FALSE,
    color = "black",
    size = 4.5,
    fontface = "bold",
    hjust = 0
  ) +
  annotate(
    "text",
    x = text_position_x - 0.2,
    y = length(levels(plot_data_renamed$Feature)) + 0.8,
    label = "Mean(|SHAP|)",
    size = 4,
    fontface = "bold",
    hjust = 0,
    vjust = 0
  ) +
  scale_color_gradientn(
    colors = pub_colors,
    name = "Feature Value",
    breaks = c(0.05, 0.95),
    labels = c("Low", "High"),
    guide = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      barheight = unit(5.0, "cm"),
      barwidth = unit(0.8, "cm"),
      ticks = FALSE,
      frame.colour = "black",
      frame.linewidth = 0.5,
      title.theme = element_text(angle = 90, size = 16, face = "bold")
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  theme_classic() +
  theme(
    text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", margin = ggplot2::margin(t = 10)),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = c(0.95, 0.5),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 14, face = "bold", color = "black"),
    plot.margin = unit(c(1, 3.8, 0.5, 0.5), "cm"),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      margin = ggplot2::margin(b = 15)
    )
  ) +
  labs(
    x = "SHAP value",
    title = "(a) Feature importance"
  ) +
  coord_cartesian(xlim = c(data_min, data_max), clip = "off")

print(p_summary_final)

# ------------------------------------------------------------------------------
# Save SHAP outputs
# ------------------------------------------------------------------------------

write.csv(
  plot_data_renamed,
  OUT_FILE("SHAP_Plot_Data"),
  row.names = FALSE
)

write.csv(
  imp_labels,
  OUT_FILE("SHAP_Feature_Importance_Rank"),
  row.names = FALSE
)

saveRDS(
  list(
    bundle_file = bundle_file,
    best_model_name = best_model_name,
    plot_data = plot_data_renamed,
    ranking_data = imp_labels,
    shap_matrix = raw_shap_matrix,
    final_plot = p_summary_final
  ),
  OUT_FILE("SHAP_Reproducible_Environment", ext = "rds")
)

ggsave(
  filename = OUT_FILE("SHAP_Summary_Plot", ext = "png"),
  plot = p_summary_final,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = OUT_FILE("SHAP_Summary_Plot", ext = "pdf"),
  plot = p_summary_final,
  width = 10,
  height = 8,
  bg = "white"
)

cat("SHAP analysis finished.\n")


# ==============================================================================
# 17. PDP analysis for the selected XGBoost model
# ==============================================================================

cat("\n>>> Running PDP analysis...\n")

pacman::p_load(
  ggplot2,
  dplyr,
  pdp,
  patchwork,
  cowplot,
  scales,
  xgboost
)

bundle_file <- OUT_FILE("SHAP_Ready_Bundle", ext = "rds")
raw_data_file <- "datasets.csv"

if (!file.exists(bundle_file)) {
  stop("Cannot find SHAP-ready bundle: ", bundle_file)
}

if (!file.exists(raw_data_file)) {
  stop("Cannot find raw dataset file: ", raw_data_file)
}

bundle <- readRDS(bundle_file)
raw_df <- read.csv(raw_data_file, stringsAsFactors = FALSE)

if (bundle$model_class != "xgboost_native") {
  stop(
    "PDP section is currently aligned to the selected native XGBoost model. ",
    "The selected model is: ", bundle$best_model_name
  )
}

target_model <- bundle$model_object
pdp_train_data <- as.data.frame(bundle$X)

if (!inherits(target_model, "xgb.Booster")) {
  stop("The selected model object is not a native xgboost booster.")
}

continuous_features <- c("Solubility", "log.Kow", "Size")
categorical_feature_model <- "Feeding"
categorical_feature_raw   <- "Feeding.habits"

missing_pdp_features <- setdiff(
  c(continuous_features, categorical_feature_model),
  colnames(pdp_train_data)
)

if (length(missing_pdp_features) > 0) {
  stop(
    "Missing PDP feature columns in bundle$X: ",
    paste(missing_pdp_features, collapse = ", ")
  )
}

if (!(categorical_feature_raw %in% colnames(raw_df))) {
  stop("Cannot find raw category column: ", categorical_feature_raw)
}

if (nrow(raw_df) != nrow(pdp_train_data)) {
  stop("Row count mismatch between raw_df and bundle$X.")
}

output_width  <- 13
output_height <- 10
output_dpi    <- 600

feature_labels <- c(
  "Solubility" = "Water solubility (mg/L)",
  "log.Kow"    = "Log Kow (Hydrophobicity)",
  "Size"       = "Maximum body size (mm)",
  "Feeding"    = "Feeding habits"
)

panel_titles <- c(
  "Solubility" = "(b) Environmental fate",
  "Feeding"    = "(c) Behavioral exposure",
  "log.Kow"    = "(d) Membrane transport",
  "Size"       = "(e) Internal metabolism"
)

use_log_x <- c(
  "Solubility" = FALSE,
  "log.Kow"    = FALSE,
  "Size"       = TRUE
)

custom_breaks_list <- list(
  "Solubility" = NULL,
  "log.Kow"    = NULL,
  "Size"       = c(1, 10, 50, 100)
)

custom_labels_list <- list(
  "Solubility" = NULL,
  "log.Kow"    = NULL,
  "Size"       = c("1", "10", "50", "100")
)

x_limits_list <- list(
  "Solubility" = NULL,
  "log.Kow"    = NULL,
  "Size"       = c(0.8, 150)
)

feeding_order <- c("Predator", "Scraper", "Gatherer", "Shredder", "Filterer")

colors_pub <- c(
  Blue   = "#4DBBD5FF",
  Red    = "#E63946",
  Purple = "#8338EC"
)

pub_theme <- cowplot::theme_half_open(font_size = 16, font_family = "sans") +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0,
      margin = ggplot2::margin(b = 10),
      color = "black"
    ),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(color = "black", size = 14, face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = ggplot2::margin(20, 20, 10, 20),
    legend.position = "none"
  )

predict_xgb_numeric <- function(model, newdata_df) {
  new_mat <- as.matrix(newdata_df)
  dnew <- xgboost::xgb.DMatrix(data = new_mat)
  as.numeric(predict(model, newdata = dnew))
}

normalize_feeding_labels <- function(x) {
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_low <- tolower(x_chr)
  
  dplyr::case_when(
    x_low %in% c("predator", "predators") ~ "Predator",
    x_low %in% c("scraper", "scrapers") ~ "Scraper",
    x_low %in% c("gatherer", "gatherers") ~ "Gatherer",
    x_low %in% c("shredder", "shredders") ~ "Shredder",
    x_low %in% c("filterer", "filterers", "filter feeder", "filter-feeder") ~ "Filterer",
    TRUE ~ x_chr
  )
}

make_clean_pdp <- function(model, train_data, feature, log_x = FALSE, grid_res = 100) {
  
  pdp_raw <- pdp::partial(
    object = model,
    pred.var = feature,
    train = train_data,
    pred.fun = function(object, newdata) predict_xgb_numeric(object, newdata),
    ice = FALSE,
    grid.resolution = grid_res
  )
  
  pdp_raw <- as.data.frame(pdp_raw)
  names(pdp_raw)[1:2] <- c("x", "yhat")
  
  pdp_dat <- pdp_raw %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(
      yhat = mean(yhat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(x)
  
  raw_vals <- train_data[[feature]]
  
  if (log_x) {
    valid_vals <- raw_vals[is.finite(raw_vals) & raw_vals > 0]
    if (length(valid_vals) > 1) {
      dens <- density(log10(valid_vals), na.rm = TRUE)
      dens_df <- data.frame(x = 10^(dens$x), y = dens$y)
    } else {
      dens_df <- data.frame(x = numeric(0), y = numeric(0))
    }
  } else {
    valid_vals <- raw_vals[is.finite(raw_vals)]
    if (length(valid_vals) > 1) {
      dens <- density(valid_vals, na.rm = TRUE)
      dens_df <- data.frame(x = dens$x, y = dens$y)
    } else {
      dens_df <- data.frame(x = numeric(0), y = numeric(0))
    }
  }
  
  list(pdp = pdp_dat, dens = dens_df)
}

plot_pdp_smooth <- function(pdp_dat, dens_df, label_x, title, color_line,
                            log_x = FALSE,
                            custom_breaks = NULL,
                            custom_labels = NULL,
                            x_limits = NULL,
                            strip_ratio = 0.035,
                            y_expand_top = 0.05,
                            y_expand_bottom = 0.06,
                            smooth_method = "loess",
                            smooth_span = 0.35) {
  
  y_hat_min <- min(pdp_dat$yhat, na.rm = TRUE)
  y_hat_max <- max(pdp_dat$yhat, na.rm = TRUE)
  y_range_total <- y_hat_max - y_hat_min
  
  if (!is.finite(y_range_total) || y_range_total == 0) {
    y_range_total <- 1
  }
  
  strip_height <- y_range_total * strip_ratio
  strip_baseline <- y_hat_min - y_range_total * y_expand_bottom - strip_height
  y_top <- y_hat_max + y_range_total * y_expand_top
  
  if (nrow(dens_df) > 0) {
    dens_df$y_scaled <- strip_baseline + (dens_df$y / max(dens_df$y)) * strip_height
  }
  
  p <- ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(x = x, ymin = strip_baseline, ymax = y_scaled),
      fill = color_line,
      alpha = 0.22
    ) +
    geom_smooth(
      data = pdp_dat,
      aes(x = x, y = yhat),
      method = smooth_method,
      formula = y ~ x,
      span = smooth_span,
      se = TRUE,
      level = 0.80,
      color = color_line,
      fill = "gray85",
      linewidth = 1.7,
      alpha = 0.45
    ) +
    labs(
      title = title,
      x = label_x,
      y = NULL
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = c(strip_baseline, y_top),
      expand = FALSE
    ) +
    pub_theme
  
  if (log_x) {
    if (!is.null(custom_breaks)) {
      p <- p + scale_x_log10(breaks = custom_breaks, labels = custom_labels)
    } else {
      p <- p + scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
    }
  } else {
    if (!is.null(custom_breaks)) {
      p <- p + scale_x_continuous(
        breaks = custom_breaks,
        labels = custom_labels
      )
    }
  }
  
  p
}

# ------------------------------------------------------------------------------
# PDP data
# ------------------------------------------------------------------------------

pdp_sol  <- make_clean_pdp(target_model, pdp_train_data, "Solubility", log_x = FALSE, grid_res = 100)
pdp_kow  <- make_clean_pdp(target_model, pdp_train_data, "log.Kow",    log_x = FALSE, grid_res = 100)
pdp_size <- make_clean_pdp(target_model, pdp_train_data, "Size",       log_x = TRUE,  grid_res = 100)

all_predictions <- predict_xgb_numeric(target_model, pdp_train_data)

feeding_labels_raw <- normalize_feeding_labels(raw_df[[categorical_feature_raw]])

pdp_feeding_raw_manual <- data.frame(
  yhat = all_predictions,
  Feeding = feeding_labels_raw
)

real_levels <- intersect(
  feeding_order,
  unique(as.character(pdp_feeding_raw_manual$Feeding))
)

if (length(real_levels) == 0) {
  real_levels <- unique(as.character(pdp_feeding_raw_manual$Feeding))
}

pdp_feeding_raw_manual$Feeding <- factor(
  pdp_feeding_raw_manual$Feeding,
  levels = real_levels
)

# Save PDP source data
write.csv(
  pdp_sol$pdp,
  OUT_FILE("PDP_Data_Solubility"),
  row.names = FALSE
)

write.csv(
  pdp_kow$pdp,
  OUT_FILE("PDP_Data_LogKow"),
  row.names = FALSE
)

write.csv(
  pdp_size$pdp,
  OUT_FILE("PDP_Data_Size"),
  row.names = FALSE
)

write.csv(
  pdp_feeding_raw_manual,
  OUT_FILE("PDP_Data_Feeding"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# PDP plots
# ------------------------------------------------------------------------------

p1 <- plot_pdp_smooth(
  pdp_dat = pdp_sol$pdp,
  dens_df = pdp_sol$dens,
  label_x = feature_labels["Solubility"],
  title = panel_titles["Solubility"],
  color_line = colors_pub["Blue"],
  log_x = use_log_x["Solubility"],
  custom_breaks = custom_breaks_list[["Solubility"]],
  custom_labels = custom_labels_list[["Solubility"]],
  x_limits = x_limits_list[["Solubility"]],
  smooth_span = 0.30
)

p2 <- plot_pdp_smooth(
  pdp_dat = pdp_kow$pdp,
  dens_df = pdp_kow$dens,
  label_x = feature_labels["log.Kow"],
  title = panel_titles["log.Kow"],
  color_line = colors_pub["Red"],
  log_x = use_log_x["log.Kow"],
  custom_breaks = custom_breaks_list[["log.Kow"]],
  custom_labels = custom_labels_list[["log.Kow"]],
  x_limits = x_limits_list[["log.Kow"]],
  smooth_span = 0.35
)

p3 <- plot_pdp_smooth(
  pdp_dat = pdp_size$pdp,
  dens_df = pdp_size$dens,
  label_x = feature_labels["Size"],
  title = panel_titles["Size"],
  color_line = colors_pub["Purple"],
  log_x = use_log_x["Size"],
  custom_breaks = custom_breaks_list[["Size"]],
  custom_labels = custom_labels_list[["Size"]],
  x_limits = x_limits_list[["Size"]],
  smooth_span = 0.40
)

p4 <- ggplot(pdp_feeding_raw_manual, aes(x = Feeding, y = yhat)) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.6,
    color = "gray40"
  ) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.6,
    alpha = 0.9,
    outlier.shape = NA,
    color = "gray40",
    fill = "#C7E9C0",
    fatten = 0
  ) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    width = 0.5,
    linewidth = 0.8,
    color = "darkblue"
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    color = "darkgreen",
    alpha = 0.55,
    shape = 16
  ) +
  labs(
    title = panel_titles["Feeding"],
    x = feature_labels["Feeding"],
    y = NULL
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(face = "bold", size = 13)
  )

final_plot_core <- (p1 + p4) / (p2 + p3) +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    theme = theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
  )

final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(
    final_plot_core,
    x = 0.045,
    y = 0,
    width = 0.955,
    height = 1
  ) +
  cowplot::draw_label(
    "Predicted toxicity (Log scale \u03bcg/L)",
    x = 0.05,
    y = 0.5,
    angle = 90,
    fontface = "bold",
    size = 16
  )

print(final_plot)

ggsave(
  filename = OUT_FILE("PDP_4panel_XGB_SmoothCI", ext = "png"),
  plot = final_plot,
  width = output_width,
  height = output_height,
  dpi = output_dpi,
  bg = "white"
)

ggsave(
  filename = OUT_FILE("PDP_4panel_XGB_SmoothCI", ext = "pdf"),
  plot = final_plot,
  width = output_width,
  height = output_height,
  bg = "white"
)

saveRDS(
  list(
    bundle_file = bundle_file,
    pdp_sol = pdp_sol,
    pdp_kow = pdp_kow,
    pdp_size = pdp_size,
    pdp_feeding = pdp_feeding_raw_manual,
    plot = final_plot
  ),
  OUT_FILE("PDP_Reproducible_Environment", ext = "rds")
)

cat("Saved PDP plot:", OUT_FILE("PDP_4panel_XGB_SmoothCI", ext = "png"), "\n")


# ==============================================================================
# 18. Feeding prediction differences
# ==============================================================================

cat("\n>>> Testing feeding-group prediction differences...\n")

pacman::p_load(
  FSA,
  multcompView,
  dplyr
)

feeding_test_df <- pdp_feeding_raw_manual %>%
  dplyr::filter(!is.na(Feeding), !is.na(yhat)) %>%
  dplyr::mutate(
    Feeding = droplevels(as.factor(Feeding))
  )

cat("\n================ Feeding prediction significance test ================\n")
cat("Group sizes:\n")
print(table(feeding_test_df$Feeding))

cat("\nGroup summaries:\n")

feeding_summary <- feeding_test_df %>%
  dplyr::group_by(Feeding) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(yhat, na.rm = TRUE),
    sd = sd(yhat, na.rm = TRUE),
    median = median(yhat, na.rm = TRUE),
    IQR = IQR(yhat, na.rm = TRUE),
    min = min(yhat, na.rm = TRUE),
    max = max(yhat, na.rm = TRUE),
    .groups = "drop"
  )

print(feeding_summary)

kw_res <- kruskal.test(yhat ~ Feeding, data = feeding_test_df)

cat("\nKruskal-Wallis test:\n")
print(kw_res)

dunn_res <- FSA::dunnTest(
  yhat ~ Feeding,
  data = feeding_test_df,
  method = "bh"
)

cat("\nDunn post hoc test, BH-adjusted:\n")
print(dunn_res$res)

pvec <- dunn_res$res$P.adj
names(pvec) <- gsub(" - ", "-", dunn_res$res$Comparison)

letters_obj <- multcompView::multcompLetters(pvec)

letters_df <- data.frame(
  Feeding = names(letters_obj$Letters),
  Letters = letters_obj$Letters,
  stringsAsFactors = FALSE
)

letters_df$Feeding <- factor(
  letters_df$Feeding,
  levels = levels(feeding_test_df$Feeding)
)

letters_df <- letters_df %>%
  dplyr::arrange(Feeding)

cat("\nCompact letter display for manual annotation:\n")
print(letters_df)

write.csv(
  feeding_summary,
  OUT_FILE("Feeding_Group_Summary"),
  row.names = FALSE
)

write.csv(
  dunn_res$res,
  OUT_FILE("Feeding_Dunn_Posthoc_BH"),
  row.names = FALSE
)

write.csv(
  letters_df,
  OUT_FILE("Feeding_Dunn_Letters"),
  row.names = FALSE
)

saveRDS(
  list(
    feeding_test_df = feeding_test_df,
    feeding_summary = feeding_summary,
    kruskal_test = kw_res,
    dunn_test = dunn_res,
    letters = letters_df
  ),
  OUT_FILE("Feeding_Group_Test_Results", ext = "rds")
)

cat("\nSaved files:\n")
cat(" -", OUT_FILE("Feeding_Group_Summary"), "\n")
cat(" -", OUT_FILE("Feeding_Dunn_Posthoc_BH"), "\n")
cat(" -", OUT_FILE("Feeding_Dunn_Letters"), "\n")
cat(" -", OUT_FILE("Feeding_Group_Test_Results", ext = "rds"), "\n")
cat("=====================================================================\n")
cat("All downstream interpretation scripts are aligned with the training pipeline.\n")