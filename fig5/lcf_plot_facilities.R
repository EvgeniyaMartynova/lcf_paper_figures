rm(list=ls())

library(ggplot2)
library(lcfstat)
library(dplyr)
library(tidyr)
library(jsonlite)
library(RColorBrewer)

source("../settings.R")

data_folder <- "data"

air_stat_fp <- file.path(data_folder, "lcf_stat_ohio_air.csv")
waste_stat_fp <- file.path(data_folder, "lcf_stat_ohio_waste.csv")

output_folder <- "img/facilities"
dir.create(output_folder, showWarnings = FALSE)

air_stat_df <- read.csv(air_stat_fp)
air_stat_df <- air_stat_df %>%
  select(r, val)
air_stat_df$ftype <- "Air"

waste_stat_df <- read.csv(waste_stat_fp)
waste_stat_df <- waste_stat_df %>%
  select(r, val)
waste_stat_df$ftype <- "Waste"

stat_df <- rbind(air_stat_df, waste_stat_df)
stat_df <- stat_df %>%
  mutate(ftype = factor(ftype)) %>%
  mutate(r = r / 1000)

air_color <- "#A63603" 
waste_color <- "#08519C"

rmax <- 75

line_width <- 1 / ggp2_magic_number

lcf_plot <- ggplot(stat_df) +
  geom_hline(yintercept = 0, linewidth=line_width, color="gray") +
  geom_line(aes(x=r, y=val, col=ftype), linewidth=line_width) +
  scale_y_continuous(expand = c(0, 0), limits=c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_x_continuous(expand = c(0, 0), limits=c(0, rmax), breaks = c(0, 25, 50, 75)) +
  scale_colour_manual(values=c(air_color, waste_color)) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  ylab("LCF") +
  xlab("r, km") +
  theme(plot.margin = margin(5, 7.5, 0, 2),
        strip.background = element_blank(),
        strip.placement = "inside",
        axis.line = element_line(colour = "black", linewidth=line_width),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.title = element_text(size = default_pointsize, family = default_font),
        legend.position="none")

lcf_plot

lcf_pdf_path <- file.path(output_folder, "all_LCF_facilities.pdf")

pdf_out(lcf_pdf_path, width=2.1, height=2)
print(lcf_plot)
dev.off()
