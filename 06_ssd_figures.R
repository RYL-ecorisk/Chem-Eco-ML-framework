# =============================================================================
# 06_ssd_figures.R
# Main-text and supplementary Jiangsu SSD figures
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

OUT_DIR <- "outputs"
N_MC <- 10000L
SEED <- 123L

PRED_FILE <- file.path(
  OUT_DIR,
  "Jiangsu_PAH_pred_with_uncertainty_Four4.csv"
)

HC5_FILE <- file.path(
  OUT_DIR,
  "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv"
)

ACR_FILE <- file.path(
  OUT_DIR,
  "Final_ACR_all16.csv"
)

if (!file.exists(ACR_FILE)) {
  stop("Final_ACR_all16.csv not found. Run 07_acr_model.R first.")
}

pred <- read_csv(PRED_FILE, show_col_types = FALSE)
hc5 <- read_csv(HC5_FILE, show_col_types = FALSE)
acr <- read_csv(ACR_FILE, show_col_types = FALSE)

set.seed(SEED)

pred_use <- pred |>
  transmute(
    Species = Latin,
    PAH = PAHs,
    Mean = Pred_LogTox_Mean,
    SD = Sigma_Total
  ) |>
  filter(
    is.finite(Mean),
    is.finite(SD),
    SD >= 0
  )

acr_use <- acr |>
  transmute(
    PAH = Abb,
    Final_ACR
  )

hc5_js <- hc5 |>
  filter(Region == "Jiangsu") |>
  left_join(acr_use, by = "PAH") |>
  mutate(
    Acute_HC5 = HC5_Median,
    Chronic_HC5 = HC5_Median / Final_ACR
  )

pah_names <- c(
  Nap = "Naphthalene",
  Acy = "Acenaphthylene",
  Ace = "Acenaphthene",
  Flu = "Fluorene",
  Ant = "Anthracene",
  Phe = "Phenanthrene",
  Flt = "Fluoranthene",
  Pyr = "Pyrene",
  BaA = "Benzo[a]anthracene",
  Chry = "Chrysene",
  BaP = "Benzo[a]pyrene",
  BbF = "Benzo[b]fluoranthene",
  BkF = "Benzo[k]fluoranthene",
  DBA = "Dibenz[a,h]anthracene",
  BghiP = "Benzo[ghi]perylene",
  InP = "Indeno[1,2,3-cd]pyrene"
)

simulate_curves <- function(target_pahs) {
  acute_all <- list()
  chronic_all <- list()

  for (pah in target_pahs) {
    d <- pred_use |> filter(PAH == pah)
    a <- acr_use |> filter(PAH == pah)

    if (nrow(d) < 5 || nrow(a) != 1) next

    curves <- vector("list", N_MC)

    for (i in seq_len(N_MC)) {
      tox <- 10^rnorm(
        nrow(d),
        mean = d$Mean,
        sd = d$SD
      )

      tox <- sort(
        tox[is.finite(tox) & tox > 0]
      )

      if (length(tox) < 5) next

      curves[[i]] <- tibble(
        PAH = pah,
        sim_id = i,
        tox = tox,
        prob = seq_along(tox) / (length(tox) + 1)
      )
    }

    acute <- bind_rows(curves)

    acute_all[[pah]] <- acute
    chronic_all[[pah]] <- acute |>
      mutate(
        tox = tox / a$Final_ACR
      )
  }

  list(
    acute = bind_rows(acute_all),
    chronic = bind_rows(chronic_all)
  )
}

make_ssd_figure <- function(
    target_pahs,
    file_stub,
    colors
) {
  sim <- simulate_curves(target_pahs)

  points_a <- hc5_js |>
    filter(PAH %in% target_pahs) |>
    transmute(PAH, HC5 = Acute_HC5)

  points_c <- hc5_js |>
    filter(PAH %in% target_pahs) |>
    transmute(PAH, HC5 = Chronic_HC5)

  base_theme <- theme_classic(base_size = 15) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )

  p_acute <- ggplot(
    sim$acute,
    aes(
      tox,
      prob,
      group = interaction(PAH, sim_id),
      color = PAH
    )
  ) +
    geom_line(
      alpha = 0.006,
      linewidth = 0.10
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_point(
      data = points_a,
      aes(HC5, 0.05, fill = PAH),
      shape = 21,
      size = 4,
      color = "white",
      inherit.aes = FALSE
    ) +
    scale_x_log10() +
    scale_y_continuous(
      labels = function(x) x * 100,
      limits = c(0, 1.02)
    ) +
    scale_color_manual(
      values = colors,
      labels = pah_names[target_pahs]
    ) +
    scale_fill_manual(
      values = colors,
      labels = pah_names[target_pahs]
    ) +
    labs(
      x = NULL,
      y = "Fraction of species affected (%)",
      title = "Acute"
    ) +
    base_theme

  p_chronic <- ggplot(
    sim$chronic,
    aes(
      tox,
      prob,
      group = interaction(PAH, sim_id),
      color = PAH
    )
  ) +
    geom_line(
      alpha = 0.006,
      linewidth = 0.10
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_point(
      data = points_c,
      aes(HC5, 0.05, fill = PAH),
      shape = 21,
      size = 4,
      color = "white",
      inherit.aes = FALSE
    ) +
    scale_x_log10() +
    scale_y_continuous(
      labels = function(x) x * 100,
      limits = c(0, 1.02)
    ) +
    scale_color_manual(
      values = colors,
      labels = pah_names[target_pahs]
    ) +
    scale_fill_manual(
      values = colors,
      labels = pah_names[target_pahs]
    ) +
    labs(
      x = expression(paste("Concentration (", mu, "g/L)")),
      y = "Fraction of species affected (%)",
      title = "Chronic"
    ) +
    base_theme +
    theme(
      legend.position = "none"
    )

  fig <- p_acute / p_chronic

  ggsave(
    file.path(
      OUT_DIR,
      paste0(file_stub, ".png")
    ),
    fig,
    width = 10.5,
    height = 12.5,
    dpi = 600
  )
}

# Main-text figure
make_ssd_figure(
  c("Nap", "Flu", "BaA", "BaP", "BghiP"),
  "Figure_Main_Jiangsu_Acute_Chronic_SSD",
  c(
    Nap = "#3182bd",
    Flu = "#6baed6",
    BaA = "#2ca02c",
    BaP = "#ff7f0e",
    BghiP = "#e31a1c"
  )
)

# Supplementary figures
make_ssd_figure(
  c("Acy", "Ace", "Ant", "Phe"),
  "Figure_S5_3ring_PAHs",
  c(
    Acy = "#0072B2",
    Ace = "#D55E00",
    Ant = "#009E73",
    Phe = "#CC79A7"
  )
)

make_ssd_figure(
  c("Flt", "Pyr", "Chry"),
  "Figure_S6_4ring_PAHs",
  c(
    Flt = "#E69F00",
    Pyr = "#56B4E9",
    Chry = "#009E73"
  )
)

make_ssd_figure(
  c("BbF", "BkF", "DBA", "InP"),
  "Figure_S7_5_6ring_PAHs",
  c(
    BbF = "#D73027",
    BkF = "#4575B4",
    DBA = "#1A9850",
    InP = "#762A83"
  )
)

message("06 complete.")
