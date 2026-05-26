# ==============================================================================
# 06_ssd_figures_and_hc5_tables.R
# ------------------------------------------------------------------------------
# Supporting SSD figures, HC5 formatting, solubility comparison and related tables.
# Computational settings are kept unchanged where simulations are used:
#   - n_mc = 10000 for SSD visualization
#   - target PAHs and ACR/HC5 mapping unchanged
# This script contains plotting/table sections used for manuscript figures and SI.
# ==============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)

# =========================
# 0. 文件与参数
# =========================
pred_file <- "Jiangsu_PAH_pred_with_uncertainty_XGB4.csv"
hc5_file  <- "Jiangsu_PAH_HC5_Summary_10k_Median.csv"
acr_file  <- "ACR_final_16PAHs.csv"

target_pahs <- c("Nap", "Flu", "BaA", "BaP", "BghiP")
n_mc <- 10000
set.seed(123)

# =========================
# 1. 读取数据
# =========================
df_species <- read_csv(pred_file, show_col_types = FALSE)
df_hc5     <- read_csv(hc5_file, show_col_types = FALSE)
df_acr     <- read_csv(acr_file, show_col_types = FALSE)

# =========================
# 2. 统一映射与标签
# =========================
target_labels <- c(
  "Nap"   = "Naphthalene (2-ring)",
  "Flu"   = "Fluorene (3-ring)",
  "BaA"   = "Benzo[a]anthracene (4-ring)",
  "BaP"   = "Benzo[a]pyrene (5-ring)",
  "BghiP" = "Benzo[ghi]perylene (6-ring)"
)

pah_colors <- c(
  "Nap"   = "#3182bd",
  "Flu"   = "#6baed6",
  "BaA"   = "#2ca02c",
  "BaP"   = "#ff7f0e",
  "BghiP" = "#e31a1c"
)

pah_shapes <- c(
  "Nap"   = 21,
  "Flu"   = 22,
  "BaA"   = 24,
  "BaP"   = 23,
  "BghiP" = 25
)

# =========================
# 3. 自动识别列名
# =========================
# species prediction file
pah_col_pred <- if ("PAHs" %in% names(df_species)) "PAHs" else if ("PAH" %in% names(df_species)) "PAH" else stop("No PAH column found in species prediction file.")
mean_col_pred <- if ("Pred_LogTox_Mean" %in% names(df_species)) "Pred_LogTox_Mean" else stop("Pred_LogTox_Mean not found.")
sd_col_pred   <- if ("Sigma_Total" %in% names(df_species)) "Sigma_Total" else stop("Sigma_Total not found.")

# HC5 file
region_col <- if ("Region" %in% names(df_hc5)) "Region" else NULL
pah_col_hc5 <- if ("PAH" %in% names(df_hc5)) "PAH" else if ("PAHs" %in% names(df_hc5)) "PAHs" else stop("No PAH column found in HC5 file.")

acute_col <- if ("HC5_Median_Acute" %in% names(df_hc5)) {
  "HC5_Median_Acute"
} else if ("HC5_Mean" %in% names(df_hc5)) {
  "HC5_Mean"
} else {
  stop("No acute HC5 median column found in HC5 file.")
}

acute_lcl_col <- if ("HC5_LCL" %in% names(df_hc5)) "HC5_LCL" else NULL
acute_ucl_col <- if ("HC5_UCL" %in% names(df_hc5)) "HC5_UCL" else NULL

# ACR file
pah_col_acr <- if ("PAH" %in% names(df_acr)) "PAH" else if ("Abb" %in% names(df_acr)) "Abb" else stop("No PAH column found in ACR file.")
acr_final_col <- if ("ACR_final" %in% names(df_acr)) "ACR_final" else stop("ACR_final not found in ACR file.")

# =========================
# 4. 提取 acute HC5 和 ACR
# =========================
hc5_acute <- df_hc5 %>%
  {
    if (!is.null(region_col)) filter(., .data[[region_col]] == "Jiangsu") else .
  } %>%
  filter(.data[[pah_col_hc5]] %in% target_pahs) %>%
  transmute(
    PAH = .data[[pah_col_hc5]],
    HC5_Median_Acute = .data[[acute_col]],
    HC5_LCL_Acute = if (!is.null(acute_lcl_col)) .data[[acute_lcl_col]] else NA_real_,
    HC5_UCL_Acute = if (!is.null(acute_ucl_col)) .data[[acute_ucl_col]] else NA_real_
  ) %>%
  mutate(PAH = factor(PAH, levels = target_pahs))

acr_map_df <- df_acr %>%
  filter(.data[[pah_col_acr]] %in% target_pahs) %>%
  transmute(
    PAH = .data[[pah_col_acr]],
    ACR_final = .data[[acr_final_col]]
  )

acr_map <- setNames(acr_map_df$ACR_final, acr_map_df$PAH)

hc5_chronic <- hc5_acute %>%
  mutate(
    HC5_Median_Chronic = HC5_Median_Acute / acr_map[as.character(PAH)],
    HC5_LCL_Chronic = ifelse(!is.na(HC5_LCL_Acute), HC5_LCL_Acute / acr_map[as.character(PAH)], NA_real_),
    HC5_UCL_Chronic = ifelse(!is.na(HC5_UCL_Acute), HC5_UCL_Acute / acr_map[as.character(PAH)], NA_real_)
  )

# =========================
# 5. 筛选目标 PAHs 的物种预测数据
# =========================
df_sub <- df_species %>%
  filter(.data[[pah_col_pred]] %in% target_pahs) %>%
  mutate(
    PAHs = factor(.data[[pah_col_pred]], levels = target_pahs),
    Pred_LogTox_Mean = .data[[mean_col_pred]],
    Sigma_Total = .data[[sd_col_pred]]
  )

# =========================
# 6. Monte Carlo：先生成 acute，再平移成 chronic
# =========================
sim_list <- vector("list", n_mc)

for (i in seq_len(n_mc)) {
  sim_list[[i]] <- df_sub %>%
    group_by(PAHs) %>%
    mutate(
      sim_tox_acute = 10^rnorm(n(), mean = Pred_LogTox_Mean, sd = Sigma_Total),
      sim_id = i
    ) %>%
    arrange(PAHs, sim_tox_acute) %>%
    mutate(prob = row_number() / (n() + 1)) %>%
    ungroup()
}

all_sims <- bind_rows(sim_list) %>%
  mutate(
    PAHs = factor(PAHs, levels = target_pahs),
    acr_val = acr_map[as.character(PAHs)],
    sim_tox_chronic = sim_tox_acute / acr_val
  )

# =========================
# 7. 锚点表
# =========================
acute_points <- hc5_acute %>%
  transmute(
    PAHs = factor(PAH, levels = target_pahs),
    HC5 = HC5_Median_Acute
  )

chronic_points <- hc5_chronic %>%
  transmute(
    PAHs = factor(PAH, levels = target_pahs),
    HC5 = HC5_Median_Chronic
  )

# =========================
# 8. 公共主题
# =========================
base_theme <- theme_classic(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(color = "black", face = "bold", size = 16),
    legend.position = c(0.22, 0.62),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    legend.spacing.y = unit(0.35, "cm"),
    legend.key.height = unit(0.8, "cm"),
    plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm")
  )

# =========================
# 9. Acute SSD
# =========================
p_acute <- ggplot() +
  geom_line(
    data = all_sims,
    aes(
      x = sim_tox_acute,
      y = prob,
      group = interaction(PAHs, sim_id),
      color = PAHs
    ),
    alpha = 0.006, linewidth = 0.1
  ) +
  geom_hline(
    yintercept = 0.05,
    color = "grey60",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    data = acute_points,
    aes(x = HC5, y = 0.05, fill = PAHs, shape = PAHs),
    size = 5.5, color = "white", stroke = 1.1
  ) +
  annotate(
    "text",
    x = 0.02, y = 0.98,
    label = "Acute",
    hjust = 0, vjust = 1,
    size = 8, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2e5, y = 0.20,
    label = paste0("MC simulations\n= ", format(n_mc, big.mark = ",")),
    size = 5.3, hjust = 0.5
  ) +
  scale_x_log10(
    limits = c(0.01, 1e6),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    labels = function(x) x * 100,
    breaks = c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0),
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = pah_colors, labels = target_labels) +
  scale_fill_manual(values = pah_colors, labels = target_labels) +
  scale_shape_manual(values = pah_shapes, labels = target_labels) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL, fill = NULL, shape = NULL
  ) +
  guides(
    fill = guide_legend(byrow = TRUE),
    color = guide_legend(byrow = TRUE),
    shape = guide_legend(byrow = TRUE)
  ) +
  base_theme

# =========================
# 10. Chronic SSD（由 acute 平移）
# =========================
p_chronic <- ggplot() +
  geom_line(
    data = all_sims,
    aes(
      x = sim_tox_chronic,
      y = prob,
      group = interaction(PAHs, sim_id),
      color = PAHs
    ),
    alpha = 0.006, linewidth = 0.1
  ) +
  geom_hline(
    yintercept = 0.05,
    color = "grey60",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    data = chronic_points,
    aes(x = HC5, y = 0.05, fill = PAHs, shape = PAHs),
    size = 5.5, color = "white", stroke = 1.1
  ) +
  annotate(
    "text",
    x = 2e-4, y = 0.98,
    label = "Chronic",
    hjust = 0, vjust = 1,
    size = 8, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2e5, y = 0.20,
    label = "n = 360\n(Number of species)",
    size = 5.3, hjust = 0.5
  ) +
  scale_x_log10(
    limits = c(1e-4, 1e6),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    labels = function(x) x * 100,
    breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = pah_colors, labels = target_labels) +
  scale_fill_manual(values = pah_colors, labels = target_labels) +
  scale_shape_manual(values = pah_shapes, labels = target_labels) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL, fill = NULL, shape = NULL
  ) +
  guides(
    fill = guide_legend(byrow = TRUE),
    color = guide_legend(byrow = TRUE),
    shape = guide_legend(byrow = TRUE)
  ) +
  base_theme

# 如果你只想保留一个图例，可把下面这句打开
# p_chronic <- p_chronic + theme(legend.position = "none")

# =========================
# 11. 拼图：共用 x / y 轴标题
# =========================
p_panels <- p_acute / p_chronic +
  plot_layout(heights = c(1, 1))

shared_y <- wrap_elements(
  full = textGrob(
    "Fraction of Species Affected (%)",
    x = 0.65,   # 往右靠一点（默认大致0.5）
    y = 0.5,
    rot = 90,
    gp = gpar(fontsize = 22, fontface = "bold")
  )
)

shared_x <- wrap_elements(
  full = textGrob(
    "Concentration (\u03bcg/L)",
    x = 0.5,
    y = 0.72,   # 往上靠一点
    gp = gpar(fontsize = 22, fontface = "bold")
  )
)

p_body <- shared_y + p_panels +
  plot_layout(widths = c(0.05, 1))   # 左侧再紧一点

p_combined <- p_body / shared_x +
  plot_layout(heights = c(1, 0.05))  # 底部再紧一点

print(p_combined)

ggsave(
  "Jiangsu_Acute_Chronic_SSD_combined_Nap_Flu_BaA_BaP_BghiP.png",
  p_combined,
  width = 11,
  height = 15,
  dpi = 600
)


library(readr)
library(dplyr)

# =========================
# 1. 读取最新数据
# =========================
hc5_all <- read_csv(
  "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv",
  show_col_types = FALSE
)

acr_df <- read_csv(
  "ACR_final_16PAHs.csv",
  show_col_types = FALSE
)

# =========================
# 2. 定义格式化函数
# =========================
fmt_hc5 <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 100 ~ sprintf("%.1f", x),
    x >= 10  ~ sprintf("%.2f", x),
    TRUE     ~ sprintf("%.3f", x)
  )
}

# =========================
# 3. 合并 ACR 信息
#    最新 HC5 文件里是 acute HC5
# =========================
hc5_joined <- hc5_all %>%
  left_join(
    acr_df %>%
      select(
        PAH, Name, Ring,
        ACR_source, ACR_final, ACR_model_error
      ),
    by = "PAH"
  )

# =========================
# 4. 计算 chronic HC5 并格式化
# =========================
final_out <- hc5_joined %>%
  mutate(
    # ---- chronic = acute / ACR_final ----
    HC5_Median_Chronic_calc = HC5_Median_Acute / ACR_final,
    HC5_LCL_Chronic_calc    = HC5_LCL_Acute / ACR_final,
    HC5_UCL_Chronic_calc    = HC5_UCL_Acute / ACR_final,
    
    # ---- acute 格式化 ----
    Acute_HC5_fmt = paste0(
      fmt_hc5(HC5_Median_Acute), " (",
      fmt_hc5(HC5_LCL_Acute), "-",
      fmt_hc5(HC5_UCL_Acute), ")"
    ),
    
    # ---- chronic 格式化 ----
    Chronic_HC5_fmt = paste0(
      fmt_hc5(HC5_Median_Chronic_calc), " (",
      fmt_hc5(HC5_LCL_Chronic_calc), "-",
      fmt_hc5(HC5_UCL_Chronic_calc), ")"
    )
  ) %>%
  select(
    Region, PAH, Name, Ring, N_species,
    ACR_source, ACR_final, ACR_model_error,
    HC5_Median_Acute, HC5_LCL_Acute, HC5_UCL_Acute, Acute_HC5_fmt,
    HC5_Median_Chronic_calc, HC5_LCL_Chronic_calc, HC5_UCL_Chronic_calc,
    Chronic_HC5_fmt,
    Best_Dist, Win_Rate_Pct, Second_Best, Second_Rate_Pct, Success_Rate
  ) %>%
  arrange(
    factor(Region, levels = c("Jiangsu", "England")),
    Ring,
    PAH
  )

print(final_out)

write_csv(final_out, "HC5_Acute_Chronic_Formatted_latest.csv")


# =========================
# 5. 论文表格简洁版
# =========================
table_out <- final_out %>%
  select(
    Region, PAH, Name, Ring,
    Acute_HC5_fmt,
    Chronic_HC5_fmt
  ) %>%
  arrange(
    Ring,
    PAH,
    factor(Region, levels = c("Jiangsu", "England"))
  )

print(table_out)

write_csv(table_out, "HC5_Acute_Chronic_Formatted_Table_latest.csv")




library(readr)
library(dplyr)

# =========================
# 1. 读取最新数据
# =========================
hc5_all <- read_csv(
  "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv",
  show_col_types = FALSE
)

acr_df <- read_csv(
  "ACR_final_16PAHs.csv",
  show_col_types = FALSE
)

# 看一下列名，确认用的是最新表
print(names(hc5_all))
print(names(acr_df))

# =========================
# 2. 定义格式化函数
# =========================
fmt_hc5 <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 100 ~ sprintf("%.1f", x),
    x >= 10  ~ sprintf("%.2f", x),
    TRUE     ~ sprintf("%.3f", x)
  )
}

# =========================
# 3. 合并 ACR，并计算 chronic HC5
#    当前 HC5_Median_Acute 是 acute HC5
# =========================
simple_out <- hc5_all %>%
  left_join(
    acr_df %>%
      select(PAH, Name, Ring, ACR_final),
    by = "PAH"
  ) %>%
  mutate(
    Acute_HC5 = HC5_Median_Acute,
    Chronic_HC5 = HC5_Median_Acute / ACR_final,
    
    Acute_HC5_fmt = fmt_hc5(Acute_HC5),
    Chronic_HC5_fmt = fmt_hc5(Chronic_HC5)
  ) %>%
  arrange(
    factor(Region, levels = c("Jiangsu", "England")),
    Ring,
    PAH
  ) %>%
  select(
    Region,
    PAH,
    Acute_HC5,
    Chronic_HC5,
    Acute_HC5_fmt,
    Chronic_HC5_fmt
  )

print(simple_out)

write_csv(simple_out, "HC5_Acute_Chronic_without_CI_by_Region.csv")

table_simple_out <- simple_out %>%
  select(
    Region,
    PAH,
    Acute_HC5 = Acute_HC5_fmt,
    Chronic_HC5 = Chronic_HC5_fmt
  )

print(table_simple_out)

write_csv(table_simple_out, "HC5_Acute_Chronic_without_CI_by_Region_Table.csv")



library(readr)
library(dplyr)
library(stringr)

# =========================
# 1. 读取数据
# =========================
hc5 <- read_csv(
  "HC5_Acute_Chronic_without_CI_by_Region_Table.csv",
  show_col_types = FALSE
)

pah <- read_csv(
  "PAHs4.20.csv",
  show_col_types = FALSE
)

# =========================
# 2. 检查关键列
# =========================
required_hc5_cols <- c("Region", "PAH", "Acute_HC5", "Chronic_HC5")
required_pah_cols <- c("Abb", "Name", "Ring", "Solubility.mgL")

missing_hc5 <- setdiff(required_hc5_cols, names(hc5))
missing_pah <- setdiff(required_pah_cols, names(pah))

if (length(missing_hc5) > 0) {
  stop("Missing columns in HC5 file: ", paste(missing_hc5, collapse = ", "))
}

if (length(missing_pah) > 0) {
  stop("Missing columns in PAH file: ", paste(missing_pah, collapse = ", "))
}

# =========================
# 3. 合并溶解度并统一单位
#    Solubility: mg/L -> ug/L
# =========================
compare_df <- hc5 %>%
  mutate(
    Acute_HC5 = as.numeric(Acute_HC5),
    Chronic_HC5 = as.numeric(Chronic_HC5)
  ) %>%
  left_join(
    pah %>%
      transmute(
        PAH = Abb,
        Name,
        Ring,
        Solubility_mgL = Solubility.mgL,
        Solubility_ugL = Solubility.mgL * 1000
      ),
    by = "PAH"
  ) %>%
  mutate(
    Acute_ratio_to_solubility = Acute_HC5 / Solubility_ugL,
    Chronic_ratio_to_solubility = Chronic_HC5 / Solubility_ugL,
    
    Acute_exceed_solubility = Acute_HC5 > Solubility_ugL,
    Chronic_exceed_solubility = Chronic_HC5 > Solubility_ugL,
    
    Acute_status = case_when(
      is.na(Solubility_ugL) ~ "Missing solubility",
      Acute_exceed_solubility ~ "Acute HC5 > solubility",
      TRUE ~ "Acute HC5 <= solubility"
    ),
    
    Chronic_status = case_when(
      is.na(Solubility_ugL) ~ "Missing solubility",
      Chronic_exceed_solubility ~ "Chronic HC5 > solubility",
      TRUE ~ "Chronic HC5 <= solubility"
    )
  ) %>%
  arrange(
    factor(Region, levels = c("Jiangsu", "England")),
    Ring,
    PAH
  )

# =========================
# 4. 输出完整比较表
# =========================
compare_out <- compare_df %>%
  select(
    Region, PAH, Name, Ring,
    Acute_HC5, Chronic_HC5,
    Solubility_mgL, Solubility_ugL,
    Acute_ratio_to_solubility,
    Chronic_ratio_to_solubility,
    Acute_status,
    Chronic_status
  )

cat("\n===== Full HC5 vs solubility comparison =====\n")
print(compare_out, n = Inf)

write_csv(compare_out, "HC5_vs_Solubility_Comparison.csv")

# =========================
# 5. 分区域汇总
# =========================
summary_out <- compare_df %>%
  group_by(Region) %>%
  summarise(
    n_total = n(),
    
    n_acute_exceed = sum(Acute_exceed_solubility, na.rm = TRUE),
    acute_exceed_PAHs = paste(PAH[Acute_exceed_solubility], collapse = ", "),
    
    n_chronic_exceed = sum(Chronic_exceed_solubility, na.rm = TRUE),
    chronic_exceed_PAHs = paste(PAH[Chronic_exceed_solubility], collapse = ", "),
    
    .groups = "drop"
  ) %>%
  mutate(
    acute_exceed_PAHs = ifelse(acute_exceed_PAHs == "", "None", acute_exceed_PAHs),
    chronic_exceed_PAHs = ifelse(chronic_exceed_PAHs == "", "None", chronic_exceed_PAHs)
  )

cat("\n===== Summary by region =====\n")
print(summary_out)

write_csv(summary_out, "HC5_vs_Solubility_Summary.csv")

# =========================
# 6. 分区域打印结果
# =========================
for (reg in unique(compare_df$Region)) {
  
  df_reg <- compare_df %>% filter(Region == reg)
  
  cat("\n=============================\n")
  cat("Region:", reg, "\n")
  cat("=============================\n")
  
  cat("\nAcute HC5 > solubility:\n")
  acute_exceed <- df_reg %>%
    filter(Acute_exceed_solubility) %>%
    pull(PAH)
  print(if (length(acute_exceed) == 0) "None" else acute_exceed)
  
  cat("\nChronic HC5 > solubility:\n")
  chronic_exceed <- df_reg %>%
    filter(Chronic_exceed_solubility) %>%
    pull(PAH)
  print(if (length(chronic_exceed) == 0) "None" else chronic_exceed)
}

















# =========================================================
# Jiangsu PAHs SSD figures for SI
# 11 remaining PAHs divided into 3 groups
# Acute SSD = Monte Carlo + log-normal fit
# Chronic SSD = shifted from acute SSD using ACR
# Updated for latest input files
# =========================================================

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)
library(tibble)

# =========================================================
# 1. File names
# =========================================================
pred_file <- "Jiangsu_PAH_pred_with_uncertainty_XGB4.csv"
hc5_file  <- "Jiangsu_PAH_HC5_Summary_10k_Median.csv"
acr_file  <- "ACR_final_16PAHs.csv"

# =========================================================
# 2. Parameters
# =========================================================
n_mc <- 10000
n_plot_lines <- 10000
set.seed(123)

# =========================================================
# 3. Helper: pick column name automatically
# =========================================================
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(
      "Cannot find any of these columns in data: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

# =========================================================
# 4. Read data
# =========================================================
pred_df <- read_csv(pred_file, show_col_types = FALSE)
hc5_df  <- read_csv(hc5_file, show_col_types = FALSE)
acr_df  <- read_csv(acr_file, show_col_types = FALSE)

cat("\nPrediction file columns:\n")
print(names(pred_df))

cat("\nHC5 file columns:\n")
print(names(hc5_df))

cat("\nACR file columns:\n")
print(names(acr_df))

# If Region column exists, keep Jiangsu only
if ("Region" %in% names(hc5_df)) {
  hc5_df <- hc5_df %>% filter(Region == "Jiangsu")
}

# ---- prediction file common column candidates ----
species_col <- pick_col(pred_df, c("Latin", "Species", "species", "Taxon"))
pah_col_pred <- pick_col(pred_df, c("PAHs", "PAH", "Abb"))

pred_mean_col <- pick_col(
  pred_df,
  c(
    "Pred_LogTox_Mean", "Pred_mean_log10", "Pred_Mean_log10",
    "Mean_log10_ugL", "Pred_log10_ugL", "Pred_log10",
    "Predicted_log10", "Pred_mean", "PredMean", "Mean"
  )
)

pred_sd_col <- pick_col(
  pred_df,
  c(
    "Sigma_Total", "Pred_sd_log10", "Pred_SD_log10", "SD_log10",
    "Pred_total_sd_log10", "Total_SD_log10",
    "Prediction_SD_log10", "Pred_SD", "SD"
  )
)

# ---- HC5 file common column candidates ----
# 最新 HC5 文件是 acute HC5；如果列名是 HC5_Median，也能识别
pah_col_hc5 <- pick_col(hc5_df, c("PAH", "PAHs", "Abb"))

acute_hc5_col <- pick_col(
  hc5_df,
  c("HC5_Median_Acute", "HC5_Median", "Acute_HC5", "HC5_Acute")
)

# ---- ACR file common column candidates ----
pah_col_acr <- pick_col(acr_df, c("PAH", "PAHs", "Abb"))

acr_col <- pick_col(
  acr_df,
  c("ACR_final", "ACR", "ACR_geomean", "ACR_Geomean")
)

# =========================================================
# 5. Full names
# =========================================================
pah_full_labels <- c(
  "Nap"   = "Naphthalene",
  "Acy"   = "Acenaphthylene",
  "Ace"   = "Acenaphthene",
  "Flu"   = "Fluorene",
  "Ant"   = "Anthracene",
  "Phe"   = "Phenanthrene",
  "Flt"   = "Fluoranthene",
  "Pyr"   = "Pyrene",
  "BaA"   = "Benzo[a]anthracene",
  "Chry"  = "Chrysene",
  "BaP"   = "Benzo[a]pyrene",
  "BbF"   = "Benzo[b]fluoranthene",
  "BkF"   = "Benzo[k]fluoranthene",
  "DBA"   = "Dibenz[a,h]anthracene",
  "BghiP" = "Benzo[ghi]perylene",
  "InP"   = "Indeno[1,2,3-cd]pyrene"
)

# =========================================================
# 6. Groups: remaining 11 PAHs
#    Main text already used Nap, Flu, BaA, BaP, BghiP
# =========================================================
group_list <- list(
  "S1_3ring"   = c("Acy", "Ace", "Ant", "Phe"),
  "S2_4ring"   = c("Flt", "Pyr", "Chry"),
  "S3_5_6ring" = c("BbF", "BkF", "DBA", "InP")
)

group_titles <- list(
  "S1_3ring"   = "3-ring PAHs",
  "S2_4ring"   = "4-ring PAHs",
  "S3_5_6ring" = "5- and 6-ring PAHs"
)

# =========================================================
# 7. Colors and shapes
# =========================================================
group_color_list <- list(
  "S1_3ring" = c(
    "Acy" = "#0072B2",
    "Ace" = "#D55E00",
    "Ant" = "#009E73",
    "Phe" = "#CC79A7"
  ),
  "S2_4ring" = c(
    "Flt"  = "#E69F00",
    "Pyr"  = "#56B4E9",
    "Chry" = "#009E73"
  ),
  "S3_5_6ring" = c(
    "BbF" = "#D73027",
    "BkF" = "#4575B4",
    "DBA" = "#1A9850",
    "InP" = "#762A83"
  )
)

get_shape_map <- function(target_pahs) {
  shape_pool <- c(16, 17, 15, 18, 8, 3)
  setNames(shape_pool[seq_along(target_pahs)], target_pahs)
}

# =========================================================
# 8. Data cleaning
# =========================================================
pred_use <- pred_df %>%
  transmute(
    Species = .data[[species_col]],
    PAH = .data[[pah_col_pred]],
    Pred_mean_log10 = as.numeric(.data[[pred_mean_col]]),
    Pred_sd_log10 = as.numeric(.data[[pred_sd_col]])
  ) %>%
  filter(
    !is.na(Species),
    !is.na(PAH),
    !is.na(Pred_mean_log10),
    !is.na(Pred_sd_log10)
  )

acr_use <- acr_df %>%
  transmute(
    PAH = .data[[pah_col_acr]],
    ACR_final = as.numeric(.data[[acr_col]])
  ) %>%
  filter(
    !is.na(PAH),
    !is.na(ACR_final),
    ACR_final > 0
  )

hc5_use <- hc5_df %>%
  transmute(
    PAH = .data[[pah_col_hc5]],
    HC5_Median_Acute = as.numeric(.data[[acute_hc5_col]])
  ) %>%
  left_join(acr_use, by = "PAH") %>%
  mutate(
    HC5_Median_Chronic = HC5_Median_Acute / ACR_final
  ) %>%
  filter(
    !is.na(PAH),
    !is.na(HC5_Median_Acute),
    !is.na(HC5_Median_Chronic)
  )

# Check unmatched PAHs
missing_acr <- hc5_df %>%
  transmute(PAH = .data[[pah_col_hc5]]) %>%
  distinct() %>%
  anti_join(acr_use %>% distinct(PAH), by = "PAH")

if (nrow(missing_acr) > 0) {
  warning(
    "Some PAHs in HC5 file have no matched ACR: ",
    paste(missing_acr$PAH, collapse = ", ")
  )
}

cat("\nCleaned prediction data:\n")
print(pred_use %>% count(PAH))

cat("\nCleaned HC5 data with calculated chronic HC5:\n")
print(hc5_use)

# =========================================================
# 9. Helper functions
# =========================================================

# position on log axis
get_xpos <- function(xlim, frac) {
  10^(log10(xlim[1]) + frac * (log10(xlim[2]) - log10(xlim[1])))
}

# build acute/chronic SSD cloud for a given PAH group
build_group_sim <- function(target_pahs, pred_use, hc5_use,
                            n_mc = 10000, n_plot_lines = 1000) {
  
  acute_sims_list <- list()
  chronic_sims_list <- list()
  acute_points_list <- list()
  chronic_points_list <- list()
  
  for (pah in target_pahs) {
    
    dat_pred <- pred_use %>% filter(PAH == pah)
    dat_hc5  <- hc5_use %>% filter(PAH == pah)
    
    if (nrow(dat_pred) == 0) {
      warning("No prediction data for PAH: ", pah)
      next
    }
    
    if (nrow(dat_hc5) == 0) {
      warning("No HC5 data for PAH: ", pah)
      next
    }
    
    hc5_a <- dat_hc5$HC5_Median_Acute[1]
    hc5_c <- dat_hc5$HC5_Median_Chronic[1]
    shift_ratio <- hc5_c / hc5_a
    
    # only plot subset of lines, but still annotate n_mc = 10000
    sim_ids <- seq_len(n_plot_lines)
    
    acute_lines_one_pah <- vector("list", length(sim_ids))
    
    for (i in seq_along(sim_ids)) {
      
      # one Monte Carlo draw per species on log10 scale
      sampled_log10 <- rnorm(
        n = nrow(dat_pred),
        mean = dat_pred$Pred_mean_log10,
        sd   = dat_pred$Pred_sd_log10
      )
      
      sampled_tox <- 10^sampled_log10
      sampled_tox <- sampled_tox[is.finite(sampled_tox) & sampled_tox > 0]
      
      if (length(sampled_tox) < 5) next
      
      # log-normal SSD fit
      meanlog <- mean(log(sampled_tox), na.rm = TRUE)
      sdlog   <- sd(log(sampled_tox), na.rm = TRUE)
      
      # x grid for curve
      x_min <- min(sampled_tox) / 5
      x_max <- max(sampled_tox) * 5
      x_grid <- exp(seq(log(x_min), log(x_max), length.out = 250))
      prob <- plnorm(x_grid, meanlog = meanlog, sdlog = sdlog)
      
      acute_lines_one_pah[[i]] <- tibble(
        PAH = pah,
        sim_id = i,
        sim_tox = x_grid,
        prob = prob
      )
    }
    
    acute_df <- bind_rows(acute_lines_one_pah)
    
    chronic_df <- acute_df %>%
      mutate(sim_tox = sim_tox * shift_ratio)
    
    acute_sims_list[[pah]] <- acute_df
    chronic_sims_list[[pah]] <- chronic_df
    
    acute_points_list[[pah]] <- tibble(
      PAH = pah,
      HC5 = hc5_a
    )
    
    chronic_points_list[[pah]] <- tibble(
      PAH = pah,
      HC5 = hc5_c
    )
  }
  
  list(
    acute_sims = bind_rows(acute_sims_list),
    chronic_sims = bind_rows(chronic_sims_list),
    acute_points = bind_rows(acute_points_list),
    chronic_points = bind_rows(chronic_points_list)
  )
}

# =========================================================
# 10. Main plotting function
# =========================================================
make_two_panel_figure <- function(target_pahs,
                                  group_key,
                                  group_name,
                                  pred_use,
                                  hc5_use,
                                  n_mc = 10000,
                                  n_plot_lines = 1000,
                                  acute_xlim = c(1e-2, 1e6),
                                  chronic_xlim = c(1e-4, 1e6),
                                  out_prefix = "SI_Figure") {
  
  sim_obj <- build_group_sim(
    target_pahs = target_pahs,
    pred_use = pred_use,
    hc5_use = hc5_use,
    n_mc = n_mc,
    n_plot_lines = n_plot_lines
  )
  
  all_sims_acute   <- sim_obj$acute_sims
  all_sims_chronic <- sim_obj$chronic_sims
  acute_points     <- sim_obj$acute_points
  chronic_points   <- sim_obj$chronic_points
  
  # fixed order
  all_sims_acute$PAH   <- factor(all_sims_acute$PAH, levels = target_pahs)
  all_sims_chronic$PAH <- factor(all_sims_chronic$PAH, levels = target_pahs)
  acute_points$PAH     <- factor(acute_points$PAH, levels = target_pahs)
  chronic_points$PAH   <- factor(chronic_points$PAH, levels = target_pahs)
  
  pah_colors <- group_color_list[[group_key]][target_pahs]
  pah_shapes <- get_shape_map(target_pahs)
  legend_labels <- pah_full_labels[target_pahs]
  
  # number of species in this group
  n_species_total <- pred_use %>%
    filter(PAH %in% target_pahs) %>%
    distinct(Species) %>%
    nrow()
  
  # base theme: no outer box, only axes
  base_theme <- theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold", size = 16),
      axis.text  = element_text(color = "black", size = 12),
      axis.line  = element_line(color = "black", linewidth = 0.7),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 10.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  acute_label_x   <- get_xpos(acute_xlim, 0.03)
  acute_note_x    <- get_xpos(acute_xlim, 0.62)
  chronic_label_x <- get_xpos(chronic_xlim, 0.03)
  chronic_note_x  <- get_xpos(chronic_xlim, 0.62)
  
  # -------------------------
  # Acute panel
  # -------------------------
  p_acute <- ggplot() +
    geom_line(
      data = all_sims_acute,
      aes(
        x = sim_tox,
        y = prob,
        group = interaction(PAH, sim_id),
        color = PAH
      ),
      alpha = 0.015,
      linewidth = 0.18
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.6
    ) +
    geom_point(
      data = acute_points,
      aes(x = HC5, y = 0.05, color = PAH, shape = PAH),
      size = 3.8,
      stroke = 0.9
    ) +
    annotate(
      "text",
      x = acute_label_x, y = 0.98,
      label = "Acute",
      hjust = 0, vjust = 1,
      fontface = "bold", size = 5.5
    ) +
    annotate(
      "text",
      x = acute_note_x, y = 0.16,
      label = paste0("MC simulations = ", format(n_mc, big.mark = ",")),
      hjust = 0, size = 4.2
    ) +
    scale_x_log10(
      limits = acute_xlim,
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      limits = c(0, 1.05),
      breaks = c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0),
      labels = function(x) x * 100,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = pah_colors, labels = legend_labels) +
    scale_shape_manual(values = pah_shapes, labels = legend_labels) +
    labs(x = NULL, y = NULL) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 1.2, alpha = 1, size = 3.2)
      ),
      shape = "none"
    ) +
    base_theme +
    theme(
      legend.position = c(0.03, 0.82),
      legend.justification = c(0, 1)
    )
  
  # -------------------------
  # Chronic panel
  # -------------------------
  p_chronic <- ggplot() +
    geom_line(
      data = all_sims_chronic,
      aes(
        x = sim_tox,
        y = prob,
        group = interaction(PAH, sim_id),
        color = PAH
      ),
      alpha = 0.015,
      linewidth = 0.18
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.6
    ) +
    geom_point(
      data = chronic_points,
      aes(x = HC5, y = 0.05, color = PAH, shape = PAH),
      size = 3.8,
      stroke = 0.9
    ) +
    annotate(
      "text",
      x = chronic_label_x, y = 0.98,
      label = "Chronic",
      hjust = 0, vjust = 1,
      fontface = "bold", size = 5.5
    ) +
    annotate(
      "text",
      x = chronic_note_x, y = 0.16,
      label = paste0("n = ", n_species_total, " species"),
      hjust = 0, size = 4.2
    ) +
    scale_x_log10(
      limits = chronic_xlim,
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      limits = c(0, 1.05),
      breaks = c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0),
      labels = function(x) x * 100,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = pah_colors, labels = legend_labels) +
    scale_shape_manual(values = pah_shapes, labels = legend_labels) +
    labs(x = NULL, y = NULL) +
    guides(
      color = "none",
      shape = "none"
    ) +
    base_theme
  
  # -------------------------
  # Combine panels with shared axes
  # -------------------------
  p_panels <- p_acute / p_chronic +
    plot_layout(heights = c(1, 1))
  
  shared_y <- wrap_elements(
    full = textGrob(
      "Fraction of Species Affected (%)",
      rot = 90,
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
  
  shared_x <- wrap_elements(
    full = textGrob(
      "Concentration (\u03bcg/L)",
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
  
  p_body <- shared_y + p_panels +
    plot_layout(widths = c(0.065, 1))
  
  p_combined <- p_body / shared_x +
    plot_layout(heights = c(1, 0.06))
  
  print(p_combined)
  
  file_stub <- paste0(out_prefix, "_", group_key)
  
  ggsave(
    filename = paste0(file_stub, ".png"),
    plot = p_combined,
    width = 10,
    height = 12,
    dpi = 600
  )
  
  ggsave(
    filename = paste0(file_stub, ".pdf"),
    plot = p_combined,
    width = 10,
    height = 12
  )
  
  invisible(p_combined)
}

# =========================================================
# 11. Draw 3 figures
# =========================================================

fig_S1 <- make_two_panel_figure(
  target_pahs = group_list[["S1_3ring"]],
  group_key = "S1_3ring",
  group_name = group_titles[["S1_3ring"]],
  pred_use = pred_use,
  hc5_use = hc5_use,
  n_mc = n_mc,
  n_plot_lines = n_plot_lines,
  acute_xlim = c(1e-2, 1e6),
  chronic_xlim = c(1e-4, 1e6),
  out_prefix = "Fig_S1_remaining_PAHs"
)

fig_S2 <- make_two_panel_figure(
  target_pahs = group_list[["S2_4ring"]],
  group_key = "S2_4ring",
  group_name = group_titles[["S2_4ring"]],
  pred_use = pred_use,
  hc5_use = hc5_use,
  n_mc = n_mc,
  n_plot_lines = n_plot_lines,
  acute_xlim = c(1e-2, 1e6),
  chronic_xlim = c(1e-4, 1e6),
  out_prefix = "Fig_S2_remaining_PAHs"
)

fig_S3 <- make_two_panel_figure(
  target_pahs = group_list[["S3_5_6ring"]],
  group_key = "S3_5_6ring",
  group_name = group_titles[["S3_5_6ring"]],
  pred_use = pred_use,
  hc5_use = hc5_use,
  n_mc = n_mc,
  n_plot_lines = n_plot_lines,
  acute_xlim = c(1e-2, 1e6),
  chronic_xlim = c(1e-4, 1e6),
  out_prefix = "Fig_S3_remaining_PAHs"
)




library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# =========================
# 1. 读取数据
# =========================
dat <- read_csv("datasets.csv", show_col_types = FALSE)

# 查看列名
print(names(dat))

# =========================
# 2. 基本信息
# =========================
n_records <- nrow(dat)
n_species <- dat %>% distinct(Latin) %>% nrow()
n_pahs <- dat %>% distinct(PAHs) %>% nrow()

cat("Number of toxicity records:", n_records, "\n")
cat("Number of freshwater macroinvertebrate species:", n_species, "\n")
cat("Number of PAHs represented:", n_pahs, "\n")

# =========================
# 3. 检查 16 priority PAHs 覆盖情况
# =========================
priority_16 <- c(
  "Nap", "Acy", "Ace", "Flu", "Ant", "Phe", "Flt", "Pyr",
  "BaA", "Chry", "BaP", "BbF", "BkF", "BghiP", "InP", "DBA"
)

present_pahs <- sort(unique(dat$PAHs))
missing_pahs <- setdiff(priority_16, present_pahs)

cat("\nPAHs present in dataset:\n")
print(present_pahs)

cat("\nPriority PAHs with no eligible records:\n")
print(missing_pahs)

# =========================
# 4. LMW / HMW 占比
#    常见定义：2–3 rings = LMW; 4–6 rings = HMW
# =========================
dat_lmw_hmw <- dat %>%
  mutate(
    MW_group = case_when(
      Rings <= 3 ~ "LMW",
      Rings >= 4 ~ "HMW",
      TRUE ~ NA_character_
    )
  )

lmw_hmw_summary <- dat_lmw_hmw %>%
  count(MW_group) %>%
  mutate(
    Percent = 100 * n / sum(n)
  )

cat("\nLMW / HMW summary:\n")
print(lmw_hmw_summary)

# =========================
# 5. 主要类群占比
#    这里默认 Group 列就是大类群
# =========================
group_summary <- dat %>%
  count(Group, sort = TRUE) %>%
  mutate(
    Percent = 100 * n / sum(n)
  )

cat("\nTaxonomic group summary:\n")
print(group_summary)

# =========================
# 6. 毒性范围（原始浓度 ug）
# =========================
tox_summary <- dat %>%
  summarise(
    min_ug = min(ug, na.rm = TRUE),
    max_ug = max(ug, na.rm = TRUE),
    fold_range = max_ug / min_ug,
    min_logtox = min(Logtox, na.rm = TRUE),
    max_logtox = max(Logtox, na.rm = TRUE),
    log_range = max_logtox - min_logtox
  )

cat("\nToxicity range summary:\n")
print(tox_summary)

# =========================
# 7. 输出适合文中描述的汇总表
# =========================
table_s1_summary <- tibble(
  Metric = c(
    "Number of toxicity records",
    "Number of PAHs represented",
    "Number of freshwater macroinvertebrate species",
    "LMW proportion (%)",
    "HMW proportion (%)",
    "Min toxicity (ug/L)",
    "Max toxicity (ug/L)",
    "Fold range"
  ),
  Value = c(
    n_records,
    n_pahs,
    n_species,
    round(lmw_hmw_summary$Percent[lmw_hmw_summary$MW_group == "LMW"], 2),
    round(lmw_hmw_summary$Percent[lmw_hmw_summary$MW_group == "HMW"], 2),
    signif(tox_summary$min_ug, 4),
    signif(tox_summary$max_ug, 4),
    signif(tox_summary$fold_range, 4)
  )
)

cat("\nTable S1 summary info:\n")
print(table_s1_summary)

# =========================
# 8. 输出主要类群表
# =========================
group_summary_out <- group_summary %>%
  mutate(
    Percent = round(Percent, 2)
  )

# =========================
# 9. 导出
# =========================
write_csv(table_s1_summary, "Table_S1_summary_metrics.csv")
write_csv(group_summary_out, "Table_S1_taxonomic_groups.csv")
write_csv(lmw_hmw_summary %>% mutate(Percent = round(Percent, 2)),
          "Table_S1_LMW_HMW_summary.csv")

# =========================
# 10. 自动生成一段可直接放文中的结果文字
# =========================
lmw_pct <- round(lmw_hmw_summary$Percent[lmw_hmw_summary$MW_group == "LMW"], 2)
hmw_pct <- round(lmw_hmw_summary$Percent[lmw_hmw_summary$MW_group == "HMW"], 2)

# 如果你想固定四大类群顺序
group_vec <- group_summary_out %>%
  mutate(txt = paste0(Group, " (", Percent, "%)")) %>%
  pull(txt)

text_out <- paste0(
  "The compiled training dataset comprised ", n_records,
  " toxicity records for ", n_species,
  " freshwater macroinvertebrate species across ", n_pahs,
  " PAHs, with no eligible records available for ",
  paste(missing_pahs, collapse = " or "), ". ",
  "Low-molecular-weight (LMW) and high-molecular-weight (HMW) PAHs were represented in broadly similar proportions, accounting for ",
  lmw_pct, "% and ", hmw_pct, "% of the dataset, respectively. ",
  "The dataset included the following major taxonomic groups: ",
  paste(group_vec, collapse = ", "), ". ",
  "Toxicity values spanned from ",
  signif(tox_summary$min_ug, 4), " to ",
  signif(tox_summary$max_ug, 4), " ug/L, corresponding to a range of >",
  signif(tox_summary$fold_range, 3), "-fold."
)

cat("\nAuto-generated paragraph:\n")
cat(text_out, "\n")




library(readr)
library(dplyr)
library(ggplot2)

# ============================================================
# 0. 参数设置
# ============================================================

input_file <- "streamlining_summary_m4_m8_m12.csv"

# 可选：
# "main" = 只画 Full-feature model vs 4-variable model
# "si"   = 画 Full-feature model + 8-variable model + 4-variable model
plot_mode <- "si"

# ============================================================
# 1. 读取数据
# ============================================================

df <- read_csv(input_file, show_col_types = FALSE)

cat("\n===== Raw data =====\n")
print(df)

cat("\n===== Column names =====\n")
print(names(df))

# ============================================================
# 2. 检查关键列
# ============================================================

required_cols <- c(
  "Algorithm", "Subset",
  "CV_R2", "CV_R2_SD",
  "CV_RMSE", "CV_nRMSE", "CV_MAE",
  "Train_R2", "Gap_R2",
  "Status", "Qualification"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# 3. 统一接口：Algorithm / Subset 命名
# ============================================================

df_clean <- df %>%
  mutate(
    Algorithm = case_when(
      Algorithm %in% c("RF", "Bayes-RF") ~ "Bayes-RF",
      Algorithm %in% c("XGB", "XGBoost", "Bayes-XGBoost") ~ "Bayes-XGBoost",
      TRUE ~ Algorithm
    ),
    
    Version = case_when(
      Subset %in% c("M12_Full", "12vars", "Full", "Full features") ~ "Full-feature model",
      Subset %in% c("M8_Extended", "8vars", "8 variables") ~ "8-variable model",
      Subset %in% c("M4_Core", "4vars", "4 critical variables", "4 variables") ~ "4-variable model",
      TRUE ~ Subset
    ),
    
    Version_short = case_when(
      Version == "Full-feature model" ~ "Full",
      Version == "8-variable model" ~ "8V",
      Version == "4-variable model" ~ "4V",
      TRUE ~ Version
    )
  )

# ============================================================
# 4. 根据主文 / SI 选择展示版本
# ============================================================

if (plot_mode == "main") {
  
  keep_versions <- c("Full-feature model", "4-variable model")
  version_levels <- c("Full-feature model", "4-variable model")
  version_labels <- c("Full", "4V")
  output_name <- "Fig5a_Streamlined_Model_CV_R2_Main_Full_vs_4V.png"
  
} else if (plot_mode == "si") {
  
  keep_versions <- c("Full-feature model", "8-variable model", "4-variable model")
  version_levels <- c("Full-feature model", "8-variable model", "4-variable model")
  version_labels <- c("Full", "8V", "4V")
  output_name <- "FigS_Streamlining_CV_R2_Full_8V_4V.png"
  
} else {
  stop("plot_mode must be either 'main' or 'si'.")
}

df_plot <- df_clean %>%
  filter(Version %in% keep_versions) %>%
  mutate(
    Algorithm = factor(
      Algorithm,
      levels = c("Bayes-RF", "Bayes-XGBoost")
    ),
    Version = factor(
      Version,
      levels = version_levels
    ),
    Version_short = factor(
      Version_short,
      levels = version_labels
    )
  ) %>%
  arrange(Algorithm, Version)

cat("\n===== Data used for plotting =====\n")
print(df_plot)

# ============================================================
# 5. 导出整理后的表格
# ============================================================

table_out <- df_plot %>%
  select(
    Algorithm, Version, n_vars,
    CV_R2, CV_R2_SD,
    CV_RMSE, CV_nRMSE, CV_MAE,
    Train_R2, Gap_R2,
    Status, Qualification
  )

write_csv(table_out, "Streamlining_summary_for_plot.csv")

# ============================================================
# 6. 绘图参数
# ============================================================

pd <- position_dodge(width = 0.72)

# 自动设置 y 轴上限，避免标签被裁掉
y_max <- max(df_plot$CV_R2 + df_plot$CV_R2_SD, na.rm = TRUE)
y_limit <- min(1.05, ceiling((y_max + 0.06) * 10) / 10)

# ============================================================
# 7. 绘图：CV R2 bar plot
# ============================================================

p_metrics_final <- ggplot(
  df_plot,
  aes(x = Algorithm, y = CV_R2, fill = Version)
) +
  
  geom_col(
    position = pd,
    width = 0.62,
    color = "black",
    linewidth = 0.75
  ) +
  
  # 单向误差线：只显示向上 SD
  geom_errorbar(
    aes(
      ymin = CV_R2,
      ymax = CV_R2 + CV_R2_SD
    ),
    position = pd,
    width = 0.10,
    linewidth = 0.70
  ) +
  
  # 柱顶数值
  geom_text(
    aes(
      label = sprintf("%.3f", CV_R2),
      y = CV_R2 + CV_R2_SD
    ),
    position = pd,
    vjust = -0.95,
    fontface = "bold",
    size = 5.2
  ) +
  
  # 配色：Full / 8V / 4V
  scale_fill_manual(
    values = c(
      "Full-feature model" = "#92C5DE",
      "8-variable model"   = "#F4A582",
      "4-variable model"   = "#A6DBA0"
    ),
    labels = c(
      "Full-feature model" = "Full variables",
      "8-variable model"   = "8 critical variables",
      "4-variable model"   = "4 critical variables"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1.0, 0.2),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = expression(bold("Cross-validation " * italic(R)^2)),
    fill = NULL
  ) +
  
  theme_classic(base_size = 20) +
  theme(
    axis.line = element_line(linewidth = 1.15, color = "black"),
    axis.text = element_text(face = "bold", color = "black", size = 17),
    axis.title = element_text(face = "bold", size = 21),
    axis.text.x = element_text(size = 18, face = "bold"),
    
    legend.position = c(0.985, 0.985),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(0.85, "cm"),
    legend.spacing.y = grid::unit(0.18, "cm"),
    legend.text = element_text(face = "bold", size = 17),
    legend.title = element_blank(),
    
    plot.margin = margin(10, 12, 10, 10)
  )

print(p_metrics_final)

ggsave(
  output_name,
  p_metrics_final,
  width = 8,
  height = 7,
  dpi = 600
)

ggsave(
  sub("\\.png$", ".pdf", output_name),
  p_metrics_final,
  width = 8,
  height = 7
)



# -------------------------------------------------------------------------
# Figure 4c / 5b: Williams Plot
# Harmonized style with panels a and b
# Structural leverage based on streamlined 4-variable model
# -------------------------------------------------------------------------

rm(list = ls())

library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(xgboost)
library(grid)

# -------------------------------------------------------------------------
# 1. Read data
# -------------------------------------------------------------------------
train_df <- read_csv("datasets.csv", show_col_types = FALSE)
pred_df  <- read_csv("Jiangsu_PAH_input_4vars.csv", show_col_types = FALSE)
obj      <- readRDS("m4_core_xgb_prediction_bundle.rds")

# -------------------------------------------------------------------------
# 2. Variables
# -------------------------------------------------------------------------
features <- obj$features
response_col <- "Logtox"

cat("\nFeatures used:\n")
print(features)

# Check required columns
missing_train_features <- setdiff(features, names(train_df))
missing_pred_features  <- setdiff(features, names(pred_df))

if (length(missing_train_features) > 0) {
  stop("Missing features in training data: ",
       paste(missing_train_features, collapse = ", "))
}

if (length(missing_pred_features) > 0) {
  stop("Missing features in prediction data: ",
       paste(missing_pred_features, collapse = ", "))
}

if (!(response_col %in% names(train_df))) {
  stop("Response column not found in training data: ", response_col)
}

# -------------------------------------------------------------------------
# 3. Build X and y
# -------------------------------------------------------------------------
train_x <- as.data.frame(train_df[, features, drop = FALSE])
pred_x  <- as.data.frame(pred_df[, features, drop = FALSE])
train_y <- train_df[[response_col]]

# -------------------------------------------------------------------------
# 4. Encode categorical variables consistently
# -------------------------------------------------------------------------
ref_levels <- list()

for (nm in names(train_x)) {
  if (is.character(train_x[[nm]]) || is.factor(train_x[[nm]])) {
    ref_levels[[nm]] <- unique(train_x[[nm]])
  }
}

encode_by_ref <- function(df_in, ref_levels) {
  df_out <- df_in
  
  for (nm in names(df_out)) {
    
    if (is.character(df_out[[nm]]) || is.factor(df_out[[nm]])) {
      
      if (!(nm %in% names(ref_levels))) {
        stop("No reference levels available for categorical variable: ", nm)
      }
      
      df_out[[nm]] <- as.numeric(
        factor(df_out[[nm]], levels = ref_levels[[nm]])
      )
      
    } else {
      df_out[[nm]] <- as.numeric(df_out[[nm]])
    }
  }
  
  df_out
}

train_x_num <- encode_by_ref(train_x, ref_levels)
pred_x_num  <- encode_by_ref(pred_x, ref_levels)

# Check encoded data
if (anyNA(train_x_num)) {
  stop("NA values found in encoded training predictors. Please check training data encoding.")
}

if (anyNA(pred_x_num)) {
  warning(
    "NA values found in encoded prediction predictors. ",
    "Some prediction categories may not exist in the training data. ",
    "Rows with NA will be removed for leverage calculation."
  )
}

# Keep complete prediction rows for leverage calculation
pred_complete_idx <- complete.cases(pred_x_num)
pred_x_num_clean <- pred_x_num[pred_complete_idx, , drop = FALSE]

cat("\nTraining samples:", nrow(train_x_num), "\n")
cat("Prediction samples:", nrow(pred_x_num), "\n")
cat("Prediction samples used for leverage:", nrow(pred_x_num_clean), "\n")

# -------------------------------------------------------------------------
# 5. Training prediction and standardized residuals
# -------------------------------------------------------------------------
dtrain <- xgb.DMatrix(data = as.matrix(train_x_num))
train_pred <- predict(obj$model, dtrain)

res_train <- train_y - train_pred
res_train_std <- as.numeric(scale(res_train))

# -------------------------------------------------------------------------
# 6. Hat values / structural leverage
# -------------------------------------------------------------------------

# Robust matrix inverse using SVD
safe_inv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  positive <- s$d > tol * max(s$d)
  
  if (!any(positive)) {
    stop("Matrix inversion failed: all singular values are near zero.")
  }
  
  s$v[, positive, drop = FALSE] %*%
    diag(1 / s$d[positive], nrow = sum(positive)) %*%
    t(s$u[, positive, drop = FALSE])
}

calc_hat_train <- function(X) {
  X <- as.matrix(X)
  X <- cbind(Intercept = 1, X)
  
  XtX_inv <- safe_inv(t(X) %*% X)
  
  as.numeric(diag(X %*% XtX_inv %*% t(X)))
}

calc_hat_pred <- function(X_train, X_pred) {
  X_train <- as.matrix(X_train)
  X_pred  <- as.matrix(X_pred)
  
  X_train <- cbind(Intercept = 1, X_train)
  X_pred  <- cbind(Intercept = 1, X_pred)
  
  XtX_inv <- safe_inv(t(X_train) %*% X_train)
  
  as.numeric(diag(X_pred %*% XtX_inv %*% t(X_pred)))
}

h_train <- calc_hat_train(train_x_num)
h_pred  <- calc_hat_pred(train_x_num, pred_x_num_clean)

# -------------------------------------------------------------------------
# 7. Critical leverage
# -------------------------------------------------------------------------
p <- ncol(train_x_num)
n <- nrow(train_x_num)

h_star <- 3 * (p + 1) / n
pct_pred_in <- round(mean(h_pred <= h_star, na.rm = TRUE) * 100, 1)

cat("\nCritical leverage h*:", round(h_star, 4), "\n")
cat("Prediction within AD:", pct_pred_in, "%\n")

# -------------------------------------------------------------------------
# 8. Plot data
# -------------------------------------------------------------------------
df_main <- data.frame(
  h   = h_train,
  res = res_train_std
)

df_dist <- bind_rows(
  data.frame(h = h_train, Group = "Training dataset"),
  data.frame(h = h_pred,  Group = "Prediction dataset")
) %>%
  filter(is.finite(h), h > 0)

# Dynamic x-axis range
x_min <- min(c(h_train, h_pred), na.rm = TRUE) * 0.8
x_max <- max(c(h_train, h_pred, h_star), na.rm = TRUE) * 1.25

# Avoid too narrow lower bound
x_min <- max(x_min, 1e-4)

x_breaks_all <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
x_breaks_use <- x_breaks_all[x_breaks_all >= x_min & x_breaks_all <= x_max]

# -------------------------------------------------------------------------
# 9. Harmonized palette
# -------------------------------------------------------------------------
col_train <- "#92C5DE"
col_pred  <- "#A6DBA0"
col_ref   <- "#D95F5F"

# -------------------------------------------------------------------------
# 10. Shared theme
# -------------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 22) +
  theme(
    axis.line  = element_line(color = "black", linewidth = 1.1),
    axis.title = element_text(face = "bold", size = 24),
    axis.text  = element_text(face = "bold", size = 19, color = "black")
  )

# -------------------------------------------------------------------------
# 11. Top density panel
# -------------------------------------------------------------------------
p_top <- ggplot(df_dist, aes(x = h, fill = Group)) +
  geom_density(
    alpha = 0.55,
    color = NA,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = h_star,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1.0
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = x_breaks_use,
    labels = x_breaks_use
  ) +
  scale_fill_manual(
    values = c(
      "Training dataset" = col_train,
      "Prediction dataset" = col_pred
    )
  ) +
  
  # 右上角文字：和图例形成右侧一列
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Prediction within AD: ", pct_pred_in, "%"),
    fontface = "bold",
    size = 6.2,
    hjust = 1.03,
    vjust = 1.35
  ) +
  
  labs(
    y = "Density",
    x = NULL
  ) +
  
  theme_pub +
  theme(
    # 图例靠右，并放在 Prediction within AD 下面
    legend.position = c(0.985, 0.70),
    legend.justification = c(1, 1),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 16),
    legend.background = element_blank(),
    legend.key = element_blank(),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.y  = element_text(face = "bold", size = 16),
    
    plot.margin = margin(b = -8, t = 8, l = 16, r = 20)
  )

# -------------------------------------------------------------------------
# 12. Main Williams plot
# -------------------------------------------------------------------------
p_main <- ggplot(df_main, aes(x = h, y = res)) +
  geom_point(
    fill = col_train,
    color = "black",
    size = 4.6,
    alpha = 0.82,
    shape = 21,
    stroke = 0.7
  ) +
  geom_vline(
    xintercept = h_star,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1.0
  ) +
  geom_hline(
    yintercept = c(-3, 3),
    linetype = "dotted",
    color = "grey55",
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.75
  ) +
  annotate(
    "text",
    x = h_star * 1.05,
    y = 4.05,
    label = paste0("Critical leverage = ", round(h_star, 3)),
    color = col_ref,
    fontface = "bold",
    hjust = 0,
    size = 6.0
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = x_breaks_use,
    labels = x_breaks_use
  ) +
  scale_y_continuous(
    limits = c(-4.5, 4.5),
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = expression(bold("Structural leverage (Hat-value, log scale)")),
    y = expression(bold("Standardized residuals (training set)"))
  ) +
  theme_pub +
  theme(
    plot.margin = margin(t = 6, l = 16, r = 16, b = 14)
  )

# -------------------------------------------------------------------------
# 13. Combine
# -------------------------------------------------------------------------
final_plot <- p_top / p_main +
  plot_layout(heights = c(1, 3.4))

print(final_plot)

# -------------------------------------------------------------------------
# 14. Save
# -------------------------------------------------------------------------
ggsave(
  "Williams_Plot_4vars_harmonized.png",
  final_plot,
  width = 10,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Williams_Plot_4vars_harmonized.pdf",
  final_plot,
  width = 10,
  height = 9,
  bg = "white"
)

# -------------------------------------------------------------------------
# 15. Console output
# -------------------------------------------------------------------------
cat("\n====================================\n")
cat("Williams plot summary\n")
cat("====================================\n")
cat("Number of training samples:", n, "\n")
cat("Number of predictors:", p, "\n")
cat("Critical leverage h*:", round(h_star, 3), "\n")
cat("Prediction within AD:", pct_pred_in, "%\n")
cat("Training samples outside h*:", sum(h_train > h_star, na.rm = TRUE), "\n")
cat("Prediction samples outside h*:", sum(h_pred > h_star, na.rm = TRUE), "\n")
cat("Training samples outside ±3 residual:", sum(abs(res_train_std) > 3, na.rm = TRUE), "\n")





# -------------------------------------------------------------------------
# Figure 4c / 5b: Williams Plot
# Harmonized style with panels a and b
# Structural leverage based on streamlined 4-variable model
# -------------------------------------------------------------------------

rm(list = ls())

library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(xgboost)
library(grid)

# -------------------------------------------------------------------------
# 1. Read data
# -------------------------------------------------------------------------
train_df <- read_csv("datasets.csv", show_col_types = FALSE)
pred_df  <- read_csv("England_PAH_input_4vars.csv", show_col_types = FALSE)
obj      <- readRDS("m4_core_xgb_prediction_bundle.rds")

# -------------------------------------------------------------------------
# 2. Variables
# -------------------------------------------------------------------------
features <- obj$features
response_col <- "Logtox"

cat("\nFeatures used:\n")
print(features)

# Check required columns
missing_train_features <- setdiff(features, names(train_df))
missing_pred_features  <- setdiff(features, names(pred_df))

if (length(missing_train_features) > 0) {
  stop("Missing features in training data: ",
       paste(missing_train_features, collapse = ", "))
}

if (length(missing_pred_features) > 0) {
  stop("Missing features in prediction data: ",
       paste(missing_pred_features, collapse = ", "))
}

if (!(response_col %in% names(train_df))) {
  stop("Response column not found in training data: ", response_col)
}

# -------------------------------------------------------------------------
# 3. Build X and y
# -------------------------------------------------------------------------
train_x <- as.data.frame(train_df[, features, drop = FALSE])
pred_x  <- as.data.frame(pred_df[, features, drop = FALSE])
train_y <- train_df[[response_col]]

# -------------------------------------------------------------------------
# 4. Encode categorical variables consistently
# -------------------------------------------------------------------------
ref_levels <- list()

for (nm in names(train_x)) {
  if (is.character(train_x[[nm]]) || is.factor(train_x[[nm]])) {
    ref_levels[[nm]] <- unique(train_x[[nm]])
  }
}

encode_by_ref <- function(df_in, ref_levels) {
  df_out <- df_in
  
  for (nm in names(df_out)) {
    
    if (is.character(df_out[[nm]]) || is.factor(df_out[[nm]])) {
      
      if (!(nm %in% names(ref_levels))) {
        stop("No reference levels available for categorical variable: ", nm)
      }
      
      df_out[[nm]] <- as.numeric(
        factor(df_out[[nm]], levels = ref_levels[[nm]])
      )
      
    } else {
      df_out[[nm]] <- as.numeric(df_out[[nm]])
    }
  }
  
  df_out
}

train_x_num <- encode_by_ref(train_x, ref_levels)
pred_x_num  <- encode_by_ref(pred_x, ref_levels)

# Check encoded data
if (anyNA(train_x_num)) {
  stop("NA values found in encoded training predictors. Please check training data encoding.")
}

if (anyNA(pred_x_num)) {
  warning(
    "NA values found in encoded prediction predictors. ",
    "Some prediction categories may not exist in the training data. ",
    "Rows with NA will be removed for leverage calculation."
  )
}

# Keep complete prediction rows for leverage calculation
pred_complete_idx <- complete.cases(pred_x_num)
pred_x_num_clean <- pred_x_num[pred_complete_idx, , drop = FALSE]

cat("\nTraining samples:", nrow(train_x_num), "\n")
cat("Prediction samples:", nrow(pred_x_num), "\n")
cat("Prediction samples used for leverage:", nrow(pred_x_num_clean), "\n")

# -------------------------------------------------------------------------
# 5. Training prediction and standardized residuals
# -------------------------------------------------------------------------
dtrain <- xgb.DMatrix(data = as.matrix(train_x_num))
train_pred <- predict(obj$model, dtrain)

res_train <- train_y - train_pred
res_train_std <- as.numeric(scale(res_train))

# -------------------------------------------------------------------------
# 6. Hat values / structural leverage
# -------------------------------------------------------------------------

# Robust matrix inverse using SVD
safe_inv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  positive <- s$d > tol * max(s$d)
  
  if (!any(positive)) {
    stop("Matrix inversion failed: all singular values are near zero.")
  }
  
  s$v[, positive, drop = FALSE] %*%
    diag(1 / s$d[positive], nrow = sum(positive)) %*%
    t(s$u[, positive, drop = FALSE])
}

calc_hat_train <- function(X) {
  X <- as.matrix(X)
  X <- cbind(Intercept = 1, X)
  
  XtX_inv <- safe_inv(t(X) %*% X)
  
  as.numeric(diag(X %*% XtX_inv %*% t(X)))
}

calc_hat_pred <- function(X_train, X_pred) {
  X_train <- as.matrix(X_train)
  X_pred  <- as.matrix(X_pred)
  
  X_train <- cbind(Intercept = 1, X_train)
  X_pred  <- cbind(Intercept = 1, X_pred)
  
  XtX_inv <- safe_inv(t(X_train) %*% X_train)
  
  as.numeric(diag(X_pred %*% XtX_inv %*% t(X_pred)))
}

h_train <- calc_hat_train(train_x_num)
h_pred  <- calc_hat_pred(train_x_num, pred_x_num_clean)

# -------------------------------------------------------------------------
# 7. Critical leverage
# -------------------------------------------------------------------------
p <- ncol(train_x_num)
n <- nrow(train_x_num)

h_star <- 3 * (p + 1) / n
pct_pred_in <- round(mean(h_pred <= h_star, na.rm = TRUE) * 100, 1)

cat("\nCritical leverage h*:", round(h_star, 4), "\n")
cat("Prediction within AD:", pct_pred_in, "%\n")

# -------------------------------------------------------------------------
# 8. Plot data
# -------------------------------------------------------------------------
df_main <- data.frame(
  h   = h_train,
  res = res_train_std
)

df_dist <- bind_rows(
  data.frame(h = h_train, Group = "Training dataset"),
  data.frame(h = h_pred,  Group = "Prediction dataset")
) %>%
  filter(is.finite(h), h > 0)

# Dynamic x-axis range
x_min <- min(c(h_train, h_pred), na.rm = TRUE) * 0.8
x_max <- max(c(h_train, h_pred, h_star), na.rm = TRUE) * 1.25

# Avoid too narrow lower bound
x_min <- max(x_min, 1e-4)

x_breaks_all <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
x_breaks_use <- x_breaks_all[x_breaks_all >= x_min & x_breaks_all <= x_max]

# -------------------------------------------------------------------------
# 9. Harmonized palette
# -------------------------------------------------------------------------
col_train <- "#92C5DE"
col_pred  <- "#A6DBA0"
col_ref   <- "#D95F5F"

# -------------------------------------------------------------------------
# 10. Shared theme
# -------------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 22) +
  theme(
    axis.line  = element_line(color = "black", linewidth = 1.1),
    axis.title = element_text(face = "bold", size = 24),
    axis.text  = element_text(face = "bold", size = 19, color = "black")
  )

# -------------------------------------------------------------------------
# 11. Top density panel
# -------------------------------------------------------------------------
p_top <- ggplot(df_dist, aes(x = h, fill = Group)) +
  geom_density(
    alpha = 0.55,
    color = NA,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = h_star,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1.0
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = x_breaks_use,
    labels = x_breaks_use
  ) +
  scale_fill_manual(
    values = c(
      "Training dataset" = col_train,
      "Prediction dataset" = col_pred
    )
  ) +
  
  # 右上角文字：和图例形成右侧一列
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Prediction within AD: ", pct_pred_in, "%"),
    fontface = "bold",
    size = 6.2,
    hjust = 1.03,
    vjust = 1.35
  ) +
  
  labs(
    y = "Density",
    x = NULL
  ) +
  
  theme_pub +
  theme(
    # 图例靠右，并放在 Prediction within AD 下面
    legend.position = c(0.985, 0.70),
    legend.justification = c(1, 1),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 16),
    legend.background = element_blank(),
    legend.key = element_blank(),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.y  = element_text(face = "bold", size = 16),
    
    plot.margin = margin(b = -8, t = 8, l = 16, r = 20)
  )

# -------------------------------------------------------------------------
# 12. Main Williams plot
# -------------------------------------------------------------------------
p_main <- ggplot(df_main, aes(x = h, y = res)) +
  geom_point(
    fill = col_train,
    color = "black",
    size = 4.6,
    alpha = 0.82,
    shape = 21,
    stroke = 0.7
  ) +
  geom_vline(
    xintercept = h_star,
    linetype = "dashed",
    color = col_ref,
    linewidth = 1.0
  ) +
  geom_hline(
    yintercept = c(-3, 3),
    linetype = "dotted",
    color = "grey55",
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.75
  ) +
  annotate(
    "text",
    x = h_star * 1.05,
    y = 4.05,
    label = paste0("Critical leverage = ", round(h_star, 3)),
    color = col_ref,
    fontface = "bold",
    hjust = 0,
    size = 6.0
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = x_breaks_use,
    labels = x_breaks_use
  ) +
  scale_y_continuous(
    limits = c(-4.5, 4.5),
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = expression(bold("Structural leverage (Hat-value, log scale)")),
    y = expression(bold("Standardized residuals (training set)"))
  ) +
  theme_pub +
  theme(
    plot.margin = margin(t = 6, l = 16, r = 16, b = 14)
  )

# -------------------------------------------------------------------------
# 13. Combine
# -------------------------------------------------------------------------
final_plot <- p_top / p_main +
  plot_layout(heights = c(1, 3.4))

print(final_plot)

# -------------------------------------------------------------------------
# 14. Save
# -------------------------------------------------------------------------
ggsave(
  "Williams_Plot_4vars_harmonized.png",
  final_plot,
  width = 10,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Williams_Plot_4vars_harmonized.pdf",
  final_plot,
  width = 10,
  height = 9,
  bg = "white"
)

# -------------------------------------------------------------------------
# 15. Console output
# -------------------------------------------------------------------------
cat("\n====================================\n")
cat("Williams plot summary\n")
cat("====================================\n")
cat("Number of training samples:", n, "\n")
cat("Number of predictors:", p, "\n")
cat("Critical leverage h*:", round(h_star, 3), "\n")
cat("Prediction within AD:", pct_pred_in, "%\n")
cat("Training samples outside h*:", sum(h_train > h_star, na.rm = TRUE), "\n")
cat("Prediction samples outside h*:", sum(h_pred > h_star, na.rm = TRUE), "\n")
cat("Training samples outside ±3 residual:", sum(abs(res_train_std) > 3, na.rm = TRUE), "\n")



# ==============================================================================
# Regional ecological risk visualization for PAHs
# Input: Risk.csv
# RQ values are already calculated.
#
# Risk levels:
#   Negligible risk: RQ < 0.01
#   Low risk:        0.01 <= RQ < 0.1
#   Moderate risk:   0.1 <= RQ < 1
#   High risk:       RQ >= 1
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(grid)

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

input_file <- "Risk.csv"

# 如果你当前文件名是 Risk(1).csv，就改成：
# input_file <- "Risk(1).csv"

output_site_risk <- "Sites_Mixture_Risk_With_Coordinates.csv"
output_long_risk <- "Sites_AllPAHs_Risk_Long.csv"
output_percent   <- "PAH_Risk_Percentages_Long_New.csv"

output_pdf <- "Figure_PAH_FullRisk_Nature.pdf"
output_png <- "Figure_PAH_FullRisk_Nature.png"

# ------------------------------------------------------------------------------
# 1. Column interface
# ------------------------------------------------------------------------------

id_cols <- c("SITE_ID", "LAT", "LNG", "Note")

pah_codes <- c(
  "Ace", "Acy", "Ant", "BaA", "BaP", "BbF", "BghiP", "BkF",
  "Chry", "DBA", "Flt", "Flu", "InP", "Nap", "Phe", "Pyr"
)

risk_vars <- c(pah_codes, "Mixture")

name_map <- c(
  "Ace"     = "Acenaphthene",
  "Acy"     = "Acenaphthylene",
  "Ant"     = "Anthracene",
  "BaA"     = "Benzo[a]anthracene",
  "BaP"     = "Benzo[a]pyrene",
  "BbF"     = "Benzo[b]fluoranthene",
  "BghiP"   = "Benzo[ghi]perylene",
  "BkF"     = "Benzo[k]fluoranthene",
  "Chry"    = "Chrysene",
  "DBA"     = "Dibenz[a,h]anthracene",
  "Flt"     = "Fluoranthene",
  "Flu"     = "Fluorene",
  "InP"     = "Indeno[1,2,3-cd]pyrene",
  "Nap"     = "Naphthalene",
  "Phe"     = "Phenanthrene",
  "Pyr"     = "Pyrene",
  "Mixture" = "Mixture"
)

risk_order <- c("High risk", "Moderate risk", "Low risk", "Negligible risk")

risk_colors <- c(
  "High risk"       = "#D73027",
  "Moderate risk"   = "#FC8D59",
  "Low risk"        = "#91BFDB",
  "Negligible risk" = "#B8DDAA"
)

# ------------------------------------------------------------------------------
# 2. Helper function
# ------------------------------------------------------------------------------

classify_risk <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.01 ~ "Negligible risk",
    x < 0.1  ~ "Low risk",
    x < 1    ~ "Moderate risk",
    x >= 1   ~ "High risk"
  )
}

# ------------------------------------------------------------------------------
# 3. Read and check data
# ------------------------------------------------------------------------------

raw_risk <- read.csv(
  input_file,
  fileEncoding = "GBK",
  check.names = FALSE
)

# 统一去除列名前后空格，避免 Excel 导出造成接口不齐
names(raw_risk) <- trimws(names(raw_risk))

required_cols <- c(id_cols, risk_vars)
missing_cols <- setdiff(required_cols, names(raw_risk))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from the input file:\n",
    paste(missing_cols, collapse = ", ")
  )
}

extra_cols <- setdiff(names(raw_risk), required_cols)

cat("\nColumns used for risk calculation:\n")
print(risk_vars)

cat("\nExtra columns ignored:\n")
print(extra_cols)

# 只转换真正参与计算的 RQ 列
raw_risk <- raw_risk %>%
  mutate(
    across(
      all_of(risk_vars),
      ~ suppressWarnings(as.numeric(as.character(.x)))
    )
  )

# ------------------------------------------------------------------------------
# 4. Export site-level mixture risk table
# ------------------------------------------------------------------------------

site_risk_table <- raw_risk %>%
  mutate(
    Mixture_Risk_Level = classify_risk(Mixture)
  ) %>%
  select(SITE_ID, LAT, LNG, Note, Mixture, Mixture_Risk_Level)

write_csv(site_risk_table, output_site_risk)

# ------------------------------------------------------------------------------
# 5. Convert to long format
# ------------------------------------------------------------------------------

risk_long <- raw_risk %>%
  select(all_of(id_cols), all_of(risk_vars)) %>%
  pivot_longer(
    cols = all_of(risk_vars),
    names_to = "PAH",
    values_to = "RQ"
  ) %>%
  mutate(
    Risk_Level = classify_risk(RQ)
  )

write_csv(risk_long, output_long_risk)

# ------------------------------------------------------------------------------
# 6. Calculate risk-level percentage for each PAH / Mixture
# ------------------------------------------------------------------------------

risk_percent_long <- risk_long %>%
  filter(!is.na(Risk_Level)) %>%
  group_by(PAH, Risk_Level) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(PAH) %>%
  mutate(
    Percentage = n / sum(n) * 100
  ) %>%
  ungroup()

write_csv(risk_percent_long, output_percent)

# ------------------------------------------------------------------------------
# 7. Complete all risk levels
# ------------------------------------------------------------------------------

processed_data <- risk_percent_long %>%
  mutate(
    PAH_Full = ifelse(PAH %in% names(name_map), name_map[PAH], PAH)
  ) %>%
  complete(
    PAH_Full,
    Risk_Level = risk_order,
    fill = list(Percentage = 0)
  )

# ------------------------------------------------------------------------------
# 8. Identify main PAHs and Other PAHs
# ------------------------------------------------------------------------------

risk_check <- processed_data %>%
  filter(
    Risk_Level != "Negligible risk",
    PAH_Full != "Mixture"
  ) %>%
  group_by(PAH_Full) %>%
  summarise(
    AtRisk_Total = sum(Percentage),
    .groups = "drop"
  ) %>%
  filter(AtRisk_Total > 0) %>%
  arrange(AtRisk_Total)

main_pahs <- risk_check$PAH_Full

other_pahs <- setdiff(
  unique(processed_data$PAH_Full),
  c("Mixture", main_pahs)
)

# ------------------------------------------------------------------------------
# 9. Prepare plotting data
# ------------------------------------------------------------------------------

plot_data <- processed_data %>%
  mutate(
    Variable = ifelse(PAH_Full %in% other_pahs, "Other PAHs", PAH_Full)
  ) %>%
  group_by(Variable, Risk_Level) %>%
  summarise(
    Percentage = mean(Percentage),
    .groups = "drop"
  )

final_order <- c("Other PAHs", main_pahs, "Mixture")

plot_data <- plot_data %>%
  mutate(
    Variable = factor(Variable, levels = final_order),
    Risk_Level = factor(Risk_Level, levels = risk_order)
  )

# ------------------------------------------------------------------------------
# 10. Plot
# ------------------------------------------------------------------------------

p_final <- ggplot(plot_data, aes(x = Variable, y = Percentage, fill = Risk_Level)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.74,
    color = "white",
    linewidth = 0.20
  ) +
  coord_flip() +
  scale_fill_manual(values = risk_colors, drop = FALSE) +
  geom_text(
    aes(label = ifelse(Percentage >= 5, paste0(round(Percentage, 0), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 6.2,              # 原来 4.8，百分比标签明显放大
    family = "serif",
    fontface = "bold",
    color = "black"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100.1),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    x = NULL,
    y = "Percentage of sampling locations (%)"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(family = "serif", color = "black", face = "bold"),
    
    axis.title.x = element_text(
      size = 23,
      face = "bold",
      margin = margin(t = 16)
    ),
    
    axis.text.x = element_text(
      size = 19,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 20,             # PAH 名称标签，重点放大
      face = "bold",
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 1.15,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 1.15,
      color = "black"
    ),
    
    axis.ticks.length = unit(0.22, "cm"),
    
    legend.position = "top",
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 19,             # 图例文字放大
      face = "bold"
    ),
    
    legend.key.width = unit(1.45, "cm"),
    legend.key.height = unit(0.65, "cm"),
    
    plot.margin = margin(t = 14, r = 34, b = 14, l = 16)
  )

# ------------------------------------------------------------------------------
# 11. Export
# ------------------------------------------------------------------------------

ggsave(output_pdf, p_final, width = 13.5, height = 8.2, device = cairo_pdf)
ggsave(output_png, p_final, width = 13.5, height = 8.2, dpi = 600)

print(p_final)













# =========================================================
# PAH uncertainty propagation (final simplified version)
# 1) Acute HC5 total uncertainty (already includes ML + SSD + Monte Carlo)
# 2) ACR uncertainty (Observed = 1; Predicted = fixed residual error)
# 3) Chronic HC5 total uncertainty
# =========================================================

rm(list = ls())

library(dplyr)
library(readr)
library(stringr)

# =========================================================
# 0. File paths
# =========================================================
file_acr <- "ACR_final_16PAHs.csv"
file_hc5 <- "Jiangsu_PAH_HC5_Summary_10k_Median.csv"

# output
out_csv <- "PAH_uncertainty_summary_final.csv"

# =========================================================
# 1. Read data
# =========================================================
acr_df <- read_csv(file_acr, show_col_types = FALSE)
hc5_df <- read_csv(file_hc5, show_col_types = FALSE)

# =========================================================
# 2. Standardize source labels
# =========================================================
acr_df <- acr_df %>%
  mutate(
    ACR_source = case_when(
      str_to_lower(ACR_source) %in% c("observed", "obs")   ~ "Observed",
      str_to_lower(ACR_source) %in% c("predicted", "pred") ~ "Predicted",
      TRUE ~ ACR_source
    )
  )

# =========================================================
# 3. Acute HC5 total uncertainty
# ---------------------------------------------------------
# This already includes:
# - ML prediction error
# - SSD uncertainty
# - Monte Carlo propagation
# =========================================================
acute_uncertainty <- hc5_df %>%
  mutate(
    acute_HC5_uncertainty_fold = HC5_UCL / HC5_LCL,
    acute_halfwidth_log10 = log10(acute_HC5_uncertainty_fold) / 2
  ) %>%
  select(
    PAH, Region, N_species,
    HC5_Median, HC5_LCL, HC5_UCL,
    acute_HC5_uncertainty_fold,
    acute_halfwidth_log10
  )

# =========================================================
# 4. Fit ring-based ACR model (for residual sigma only)
# ---------------------------------------------------------
# log10(ACR_geomean) ~ Ring
# Observed ACR values are treated as fixed (uncertainty = 1)
# Predicted ACR values use a common residual error on log10 scale
# =========================================================
acr_obs_fit <- acr_df %>%
  filter(ACR_source == "Observed", !is.na(ACR_geomean), !is.na(Ring)) %>%
  mutate(log10_ACR_geomean = log10(ACR_geomean))

fit_log_ring <- lm(log10_ACR_geomean ~ Ring, data = acr_obs_fit)

cat("\n==============================\n")
cat("ACR ring model summary\n")
cat("==============================\n")
print(summary(fit_log_ring))

# model residual SE
model_sigma <- summary(fit_log_ring)$sigma
acr_z <- qnorm(0.975)

cat("\nModel residual SE (log10 scale) = ", round(model_sigma, 4), "\n", sep = "")
cat("95% half-width for predicted ACR on log10 scale = ", round(acr_z * model_sigma, 4), "\n", sep = "")
cat("Predicted ACR uncertainty fold = ", round(10^(2 * acr_z * model_sigma), 3), "\n", sep = "")

# =========================================================
# 5. ACR uncertainty
# ---------------------------------------------------------
# Observed ACR: no extra uncertainty
# Predicted ACR: common residual uncertainty based on sigma
# =========================================================
acr_uncertainty <- acr_df %>%
  mutate(
    ACR_used = ACR_final,
    ACR_uncertainty_fold = case_when(
      ACR_source == "Observed"  ~ 1,
      ACR_source == "Predicted" ~ 10^(2 * acr_z * model_sigma),
      TRUE ~ NA_real_
    ),
    ACR_halfwidth_log10 = case_when(
      ACR_source == "Observed"  ~ 0,
      ACR_source == "Predicted" ~ acr_z * model_sigma,
      TRUE ~ NA_real_
    )
  ) %>%
  select(
    PAH, Name, Ring, ACR_source, ACR_used,
    ACR_uncertainty_fold, ACR_halfwidth_log10
  )

# =========================================================
# 6. Chronic HC5 total uncertainty
# ---------------------------------------------------------
# chronic HC5 = acute HC5 / ACR
#
# uncertainty propagated on log10 scale:
# h_chronic = sqrt(h_acute^2 + h_ACR^2)
#
# chronic_HC5_uncertainty_fold = 10^(2 * h_chronic)
# =========================================================
final_df <- acute_uncertainty %>%
  left_join(acr_uncertainty, by = "PAH") %>%
  mutate(
    chronic_HC5 = HC5_Median / ACR_used,
    
    chronic_halfwidth_log10 = sqrt(
      acute_halfwidth_log10^2 +
        ACR_halfwidth_log10^2
    ),
    
    chronic_HC5_uncertainty_fold = 10^(2 * chronic_halfwidth_log10),
    
    chronic_HC5_LCL = chronic_HC5 / (10^chronic_halfwidth_log10),
    chronic_HC5_UCL = chronic_HC5 * (10^chronic_halfwidth_log10)
  ) %>%
  arrange(Ring, PAH)

# =========================================================
# 7. Final summary table
# =========================================================
summary_table <- final_df %>%
  transmute(
    PAH,
    Ring,
    ACR_source,
    
    # acute HC5
    Acute_HC5 = HC5_Median,
    Acute_HC5_LCL = HC5_LCL,
    Acute_HC5_UCL = HC5_UCL,
    Acute_fold = acute_HC5_uncertainty_fold,
    
    # ACR
    ACR_used,
    ACR_fold = ACR_uncertainty_fold,
    
    # chronic HC5
    Chronic_HC5 = chronic_HC5,
    Chronic_HC5_LCL = chronic_HC5_LCL,
    Chronic_HC5_UCL = chronic_HC5_UCL,
    Chronic_fold = chronic_HC5_uncertainty_fold
  )

# =========================================================
# 8. Save results
# =========================================================
write_csv(summary_table, out_csv)

cat("\n=====================================\n")
cat("Saved final uncertainty summary to:\n")
cat(out_csv, "\n")
cat("=====================================\n\n")

# =========================================================
# 9. Print full table
# =========================================================
print(summary_table, n = Inf)

# =========================================================
# 10. Compact manuscript-style table
# =========================================================
manuscript_table <- summary_table %>%
  mutate(
    Acute_fold = round(Acute_fold, 2),
    ACR_fold = round(ACR_fold, 2),
    Chronic_fold = round(Chronic_fold, 2),
    Acute_HC5 = signif(Acute_HC5, 4),
    Chronic_HC5 = signif(Chronic_HC5, 4)
  )

cat("\nConcise manuscript table:\n")
print(manuscript_table, n = Inf)




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
  "Feeding"     = "Feeding group",
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
  "Feeding"    = "Feeding group"
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






library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(grid)

# =========================
# 1. 读取数据
# =========================
df <- read_csv("dddd.csv", show_col_types = FALSE)

cat("\n===== Column names =====\n")
print(names(df))

# =========================
# 2. 固定顺序
#    这里写你想显示的纵向顺序
# =========================
taxa_levels    <- c("Molluscs", "Worms", "Insects", "Crustaceans")
feeding_levels <- c("Predator", "Shredder", "Scraper", "Gatherer", "Filterer")
size_levels    <- c("Very Large", "Large", "Medium", "Small")

# =========================
# 3. 数据整理
# =========================
df_final <- df %>%
  dplyr::mutate(
    Phylum = stringr::str_to_lower(Phylum),
    Class  = stringr::str_to_lower(Class),
    
    Taxa = dplyr::case_when(
      Phylum == "mollusca" ~ "Molluscs",
      Phylum == "annelida" ~ "Worms",
      Phylum == "arthropoda" & Class == "insecta" ~ "Insects",
      Phylum == "arthropoda" & Class %in% c("malacostraca", "branchiopoda") ~ "Crustaceans",
      TRUE ~ NA_character_
    ),
    
    Feeding = stringr::str_to_title(`Feeding.habits`),
    
    Size = dplyr::case_when(
      stringr::str_to_lower(Size_class) == "small" ~ "Small",
      stringr::str_to_lower(Size_class) == "medium" ~ "Medium",
      stringr::str_to_lower(Size_class) == "large" ~ "Large",
      stringr::str_to_lower(Size_class) == "very large" ~ "Very Large",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Taxa), !is.na(Feeding), !is.na(Size)) %>%
  dplyr::mutate(
    Taxa    = factor(as.character(Taxa), levels = taxa_levels),
    Feeding = factor(as.character(Feeding), levels = feeding_levels),
    Size    = factor(as.character(Size), levels = size_levels)
  )

cat("\n===== Cleaned data preview =====\n")
print(head(df_final))

# =========================
# 4. 检查未分类行
# =========================
check_unclassified <- df %>%
  dplyr::mutate(
    Phylum = stringr::str_to_lower(Phylum),
    Class  = stringr::str_to_lower(Class),
    Taxa_check = dplyr::case_when(
      Phylum == "mollusca" ~ "Molluscs",
      Phylum == "annelida" ~ "Worms",
      Phylum == "arthropoda" & Class == "insecta" ~ "Insects",
      Phylum == "arthropoda" & Class %in% c("malacostraca", "branchiopoda") ~ "Crustaceans",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(is.na(Taxa_check))

cat("\n===== Unclassified rows =====\n")
print(check_unclassified)

# =========================
# 5. 汇总 alluvial 数据
# =========================
df_plot <- df_final %>%
  dplyr::count(Taxa, Feeding, Size, name = "Freq")

cat("\n===== Summarised alluvial data =====\n")
print(df_plot)

# =========================
# 6. 配色
# =========================
matched_palette <- c(
  "Crustaceans" = "#FFA94D",
  "Insects"     = "#4DBBD5",
  "Molluscs"    = "#7AD3A8",
  "Worms"       = "#C6A0F6"
)

# =========================
# 7. 先画一个不带标签的底图
#    关键：decreasing = NA, reverse = FALSE
# =========================
p_base <- ggplot(
  df_plot,
  aes(y = Freq, axis1 = Taxa, axis2 = Feeding, axis3 = Size)
) +
  geom_alluvium(
    aes(fill = Taxa),
    width = 0.16,
    alpha = 0.65,
    knot.pos = 0.35,
    color = "white",
    linewidth = 0.18,
    decreasing = NA,
    reverse = FALSE
  ) +
  geom_stratum(
    width = 0.24,
    fill = "grey98",
    color = "grey35",
    linewidth = 0.4,
    decreasing = NA,
    reverse = FALSE
  ) +
  scale_fill_manual(values = matched_palette) +
  scale_x_discrete(
    limits = c("Taxonomic group", "Feeding group", "Size class"),
    expand = c(0.08, 0.08)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      size = 22,
      color = "black",
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(12, 14, 12, 14)
  )

# =========================
# 8. 从 ggplot_build 里提取 stratum 真实位置
# =========================
gb <- ggplot_build(p_base)
stratum_data <- gb$data[[2]]

label_df <- stratum_data %>%
  dplyr::transmute(
    x = x,
    y = (ymin + ymax) / 2,
    stratum = as.character(stratum)
  )

# =========================
# 9. 计算百分比
# =========================
taxa_pct <- df_final %>%
  dplyr::count(Taxa, name = "Freq") %>%
  dplyr::mutate(
    stratum = as.character(Taxa),
    Percent = 100 * Freq / sum(Freq)
  ) %>%
  dplyr::select(stratum, Percent)

feeding_pct <- df_final %>%
  dplyr::count(Feeding, name = "Freq") %>%
  dplyr::mutate(
    stratum = as.character(Feeding),
    Percent = 100 * Freq / sum(Freq)
  ) %>%
  dplyr::select(stratum, Percent)

size_pct <- df_final %>%
  dplyr::count(Size, name = "Freq") %>%
  dplyr::mutate(
    stratum = as.character(Size),
    Percent = 100 * Freq / sum(Freq)
  ) %>%
  dplyr::select(stratum, Percent)

pct_df <- dplyr::bind_rows(taxa_pct, feeding_pct, size_pct)

label_df <- label_df %>%
  dplyr::left_join(pct_df, by = "stratum") %>%
  dplyr::mutate(
    label = paste0(stratum, "\n", sprintf("%.1f%%", Percent))
  )

cat("\n===== Final label data (matched to real strata positions) =====\n")
print(label_df)

# =========================
# 10. 最终图：加标签
# =========================
p_alluvial_final <- p_base +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 6.2,
    fontface = "plain",
    fill = "white",
    color = "black",
    label.size = 0.22,
    label.padding = unit(0.26, "lines"),
    lineheight = 1.00
  )

# =========================
# 11. 输出
# =========================
print(p_alluvial_final)

ggsave(
  "Figure3A_Alluvial_order_fixed_bigfont_v2.pdf",
  p_alluvial_final,
  width = 10,
  height = 7.6
)

ggsave(
  "Figure3A_Alluvial_order_fixed_bigfont_v2.png",
  p_alluvial_final,
  width = 10,
  height = 7.6,
  dpi = 600
)



# ==============================================================================
# Regional ecological risk visualization for PAHs
# Input: Risk.csv
# RQ values are already calculated.
#
# Risk levels:
#   Negligible risk: RQ < 0.01
#   Low risk:        0.01 <= RQ < 0.1
#   Moderate risk:   0.1 <= RQ < 1
#   High risk:       RQ >= 1
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(grid)

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

input_file <- "Risk.csv"

# 如果你当前文件名是 Risk(1).csv，就改成：
# input_file <- "Risk(1).csv"

output_site_risk <- "Sites_Mixture_Risk_With_Coordinates.csv"
output_long_risk <- "Sites_AllPAHs_Risk_Long.csv"
output_percent   <- "PAH_Risk_Percentages_Long_New.csv"

output_pdf <- "Figure_PAH_FullRisk_Nature.pdf"
output_png <- "Figure_PAH_FullRisk_Nature.png"

# ------------------------------------------------------------------------------
# 1. Column interface
# ------------------------------------------------------------------------------

id_cols <- c("SITE_ID", "LAT", "LNG", "Note")

pah_codes <- c(
  "Ace", "Acy", "Ant", "BaA", "BaP", "BbF", "BghiP", "BkF",
  "Chry", "DBA", "Flt", "Flu", "InP", "Nap", "Phe", "Pyr"
)

risk_vars <- c(pah_codes, "Mixture")

name_map <- c(
  "Ace"     = "Acenaphthene",
  "Acy"     = "Acenaphthylene",
  "Ant"     = "Anthracene",
  "BaA"     = "Benzo[a]anthracene",
  "BaP"     = "Benzo[a]pyrene",
  "BbF"     = "Benzo[b]fluoranthene",
  "BghiP"   = "Benzo[ghi]perylene",
  "BkF"     = "Benzo[k]fluoranthene",
  "Chry"    = "Chrysene",
  "DBA"     = "Dibenz[a,h]anthracene",
  "Flt"     = "Fluoranthene",
  "Flu"     = "Fluorene",
  "InP"     = "Indeno[1,2,3-cd]pyrene",
  "Nap"     = "Naphthalene",
  "Phe"     = "Phenanthrene",
  "Pyr"     = "Pyrene",
  "Mixture" = "Mixture"
)

risk_order <- c("High risk", "Moderate risk", "Low risk", "Negligible risk")

risk_colors <- c(
  "High risk"       = "#D73027",
  "Moderate risk"   = "#FC8D59",
  "Low risk"        = "#91BFDB",
  "Negligible risk" = "#B8DDAA"
)

# ------------------------------------------------------------------------------
# 2. Helper function
# ------------------------------------------------------------------------------

classify_risk <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.01 ~ "Negligible risk",
    x < 0.1  ~ "Low risk",
    x < 1    ~ "Moderate risk",
    x >= 1   ~ "High risk"
  )
}

# ------------------------------------------------------------------------------
# 3. Read and check data
# ------------------------------------------------------------------------------

raw_risk <- read.csv(
  input_file,
  fileEncoding = "GBK",
  check.names = FALSE
)

names(raw_risk) <- trimws(names(raw_risk))

required_cols <- c(id_cols, risk_vars)
missing_cols <- setdiff(required_cols, names(raw_risk))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from the input file:\n",
    paste(missing_cols, collapse = ", ")
  )
}

extra_cols <- setdiff(names(raw_risk), required_cols)

cat("\nColumns used for risk calculation:\n")
print(risk_vars)

cat("\nExtra columns ignored:\n")
print(extra_cols)

raw_risk <- raw_risk %>%
  mutate(
    across(
      all_of(risk_vars),
      ~ suppressWarnings(as.numeric(as.character(.x)))
    )
  )

# ------------------------------------------------------------------------------
# 4. Export site-level mixture risk table
# ------------------------------------------------------------------------------

site_risk_table <- raw_risk %>%
  mutate(
    Mixture_Risk_Level = classify_risk(Mixture)
  ) %>%
  select(SITE_ID, LAT, LNG, Note, Mixture, Mixture_Risk_Level)

write_csv(site_risk_table, output_site_risk)

# ------------------------------------------------------------------------------
# 5. Convert to long format
# ------------------------------------------------------------------------------

risk_long <- raw_risk %>%
  select(all_of(id_cols), all_of(risk_vars)) %>%
  pivot_longer(
    cols = all_of(risk_vars),
    names_to = "PAH",
    values_to = "RQ"
  ) %>%
  mutate(
    Risk_Level = classify_risk(RQ)
  )

write_csv(risk_long, output_long_risk)

# ------------------------------------------------------------------------------
# 6. Calculate risk-level percentage for each PAH / Mixture
# ------------------------------------------------------------------------------

risk_percent_long <- risk_long %>%
  filter(!is.na(Risk_Level)) %>%
  group_by(PAH, Risk_Level) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(PAH) %>%
  mutate(
    Percentage = n / sum(n) * 100
  ) %>%
  ungroup()

write_csv(risk_percent_long, output_percent)

# ------------------------------------------------------------------------------
# 7. Complete all risk levels
# ------------------------------------------------------------------------------

processed_data <- risk_percent_long %>%
  mutate(
    PAH_Full = ifelse(PAH %in% names(name_map), name_map[PAH], PAH)
  ) %>%
  complete(
    PAH_Full,
    Risk_Level = risk_order,
    fill = list(Percentage = 0)
  )

# ------------------------------------------------------------------------------
# 8. Identify main PAHs and Other PAHs
# ------------------------------------------------------------------------------

risk_check <- processed_data %>%
  filter(
    Risk_Level != "Negligible risk",
    PAH_Full != "Mixture"
  ) %>%
  group_by(PAH_Full) %>%
  summarise(
    AtRisk_Total = sum(Percentage),
    .groups = "drop"
  ) %>%
  filter(AtRisk_Total > 0) %>%
  arrange(AtRisk_Total)

main_pahs <- risk_check$PAH_Full

other_pahs <- setdiff(
  unique(processed_data$PAH_Full),
  c("Mixture", main_pahs)
)

# ------------------------------------------------------------------------------
# 9. Prepare plotting data
# ------------------------------------------------------------------------------

plot_data <- processed_data %>%
  mutate(
    Variable = ifelse(PAH_Full %in% other_pahs, "Other PAHs", PAH_Full)
  ) %>%
  group_by(Variable, Risk_Level) %>%
  summarise(
    Percentage = mean(Percentage),
    .groups = "drop"
  )

final_order <- c("Other PAHs", main_pahs, "Mixture")

plot_data <- plot_data %>%
  mutate(
    Variable = factor(Variable, levels = final_order),
    Risk_Level = factor(Risk_Level, levels = risk_order)
  )

# ------------------------------------------------------------------------------
# 10. Plot
# ------------------------------------------------------------------------------

p_final <- ggplot(plot_data, aes(x = Variable, y = Percentage, fill = Risk_Level)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.74,
    color = "white",
    linewidth = 0.22
  ) +
  coord_flip() +
  scale_fill_manual(values = risk_colors, drop = FALSE) +
  
  geom_text(
    aes(label = ifelse(Percentage >= 5, paste0(round(Percentage, 0), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 5.8,
    family = "serif",
    fontface = "bold",
    color = "black"
  ) +
  
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100.1),
    breaks = seq(0, 100, 20),
    labels = function(x) as.character(x)
  ) +
  labs(
    x = NULL,
    y = "Proportion of sampling locations (%)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 22,
      face = "bold",
      margin = margin(t = 15)
    ),
    
    axis.text.x = element_text(
      size = 18,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 1.15,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 1.15,
      color = "black"
    ),
    
    axis.ticks.length = unit(0.22, "cm"),
    
    legend.position = "top",
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 19,
      face = "bold",
      color = "black"
    ),
    
    legend.key.width = unit(1.35, "cm"),
    legend.key.height = unit(0.60, "cm"),
    
    legend.spacing.x = unit(0.30, "cm"),
    
    plot.margin = margin(
      t = 12,
      r = 28,
      b = 14,
      l = 14
    )
  )
# ------------------------------------------------------------------------------
# 11. Export
# ------------------------------------------------------------------------------

ggsave(
  output_pdf,
  p_final,
  width = 12.2,
  height = 7.5,
  device = cairo_pdf
)

ggsave(
  output_png,
  p_final,
  width = 12.2,
  height = 7.5,
  dpi = 600
)

print(p_final)