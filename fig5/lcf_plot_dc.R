rm(list=ls())

library(ggplot2)
library(spatstat)
library(lcfstat)
library(dplyr)
library(tidyr)
library(jsonlite)
library(RColorBrewer)
library(ggmagnify)

source("../utils_lcf.R")
source("../settings.R")

rmax <- 500
min_cell_num <- 40

prediction_path <- "data/cell"

output_folder <- "img/tmas"
dir.create(output_folder, showWarnings = FALSE, recursive = FALSE)

dataset <- "2022-05-19_Iris_DCpanelTMA9"
slide <- "TMA9_2A_Scan1"
patient <- "T99-11009_IF"
tma_num <- "1"

purple_cols <- c("#800066", "#b3008f",  "#e600b8", "#ff1ad1")
col_func <- colorRampPalette(purple_cols)
b_cell_palette <- col_func(7)
b_cell_palette

# Helper T palette
red_palette <- brewer.pal(9, "Reds")
col_func <- colorRampPalette(red_palette[8:4])
cdc_palette <- col_func(7)
cdc_palette

# Regulatory T color
pdc_palette <- c("#1b6031", "#20733a", "#258644", "#2a994e",  "#30ac57", "#49c571", "#5dcc81")
pdc_palette

# NK palette
yellow_cols <- c("#807400", "#b3a300", "#e6d200", "#ffeb1a")
col_func <- colorRampPalette(yellow_cols)
myeloid_palette <- col_func(7)
myeloid_palette

b_col <- "#CC00A3"
cdc_col <- "#EF3B2C"
pdc_col <- "#2a994e"
myeloid_col <- "#CCBA00"

# Load polygon from disk
tma_roi <- load_tma_pol(dataset, slide, patient, tma_num, prediction_path)

roi_coords <- get_roi_coords(tma_roi, origin_tranform=FALSE)
roi_coords[nrow(roi_coords), ] <- roi_coords[1,]

roi_coords_anticlock <- if (anticlockwise(roi_coords)) roi_coords else roi_coords[nrow(roi_coords):1, ]
tma_window <- window_from_roi(roi_coords_anticlock)

cells <- load_cells(dataset, slide, patient, tma_num, NULL, NULL,  prediction_path=prediction_path)

cells <- cells %>%
  filter(!(type %in% c("cDC1", "Double pos cDC2"))) %>%
  mutate(t=case_when(type == "B cell" ~ "B",
                     type == "cDC2" ~ "cDC2",
                     type == "pDC" ~ "pDC",
                     type == "Myeloid cell" ~ "M"))

# TODO: extract plotting code for this figure

cells_stat <- cells %>%
  group_by(t) %>%
  summarise(n = n()) %>%
  ungroup()

cells_stat <- cells_stat %>%
  filter(n >= min_cell_num)

types <- as.character(cells_stat$t)

lcf_cts <- NULL
for (ctype in types) {
  print(ctype)

  ct_coords <- cells %>%
    filter(t == ctype) %>%
    dplyr::select(X, Y)

  ct_pp <- ppp(ct_coords$X, ct_coords$Y, window = tma_window)
  lcf_ct <- LCFest(ct_pp, rmax=rmax, dim_lims = c(10, 50))

  lcf_ct <- as.data.frame(lcf_ct)
  lcf_ct[["type"]] <- ctype

  lcf_cts <- rbind(lcf_cts, lcf_ct)
}

lcf_cts <- lcf_cts %>%
  mutate(type = factor(type, levels=c("cDC2", "pDC", "B", "M")))

line_width <- 1 / ggp2_magic_number

lcf_plot <- ggplot(lcf_cts) +
  geom_hline(yintercept = 0, linewidth=line_width, color="gray") +
  geom_line(aes(x=r, y=iso, col=type), linewidth=line_width) +
  scale_y_continuous(expand = c(0, 0), limits=c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_x_continuous(expand = c(0, 0), limits=c(0, rmax)) +
  scale_colour_manual(values=c(cdc_col, pdc_col, b_col, myeloid_col)) +
  ylab("LCF") +
  xlab(expression("r," ~ mu * m)) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  theme(plot.margin = margin(5, 7.5, 0, 2),
        strip.background = element_blank(),
        strip.placement = "inside",
        axis.line = element_line(colour = "black", linewidth=line_width),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.title = element_text(size = default_pointsize, family = default_font),
        legend.position="none") +
  geom_magnify(from = c(0, 40, 0, 0.25), to = c(100, 212, -0.85, -0.15), shape = "rect",
               expand = 0, proj.linetype=2,
               colour="black", linewidth=line_width / 2)


lcf_plot

lcf_pdf_path <- file.path(output_folder, sprintf("%s_%s_%s_all_LCF_DC.pdf", slide, patient, tma_num))

pdf_out(lcf_pdf_path, width=2.1, height=2)
print(lcf_plot)
dev.off()

