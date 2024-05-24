rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggbeeswarm)
library(lme4)

source("../utils_lcf.R")
source("../settings.R")

stat_file <- "data/lcf_stat_tmas.csv"
output_folder <- "img/tmas"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Colors
cdc2_color <- "#EF3B2C"
pdc_color <- "#2a994e"
b_color <- "#CC00A3"
myeloid_color <-"#CCBA00"

ctype_colors <- c(cdc2_color, pdc_color, b_color, myeloid_color)

lcf_stat_df <- read.csv(stat_file)

lcf_stat_df <- lcf_stat_df %>%
  filter(!(ctype %in% c("Double pos cDC2", "cDC1"))) %>%
  dplyr::select(ctype, cell_num, lcf_integral)

lcf_stat_df <- lcf_stat_df %>%
  mutate(ctype = case_when(ctype == "B cell" ~ "B",
                            ctype == "cDC2" ~ "cDC2",
                            ctype == "pDC" ~ "pDC",
                            ctype == "Myeloid cell" ~ "M"
  )) %>%
  mutate(ctype = factor(ctype, levels=c("cDC2", "pDC", "B", "M")))

lcf_summ <- lcf_stat_df %>% 
  group_by(ctype) %>% 
  summarise(med = median(lcf_integral), 
            min = min(lcf_integral), 
            max = max(lcf_integral), 
            n = n(),
            lcf_cn_corr = cor(lcf_integral, log10(cell_num))) %>%
  mutate(x_tick = sprintf("n=%d", n))
  #mutate(x_tick = ctype)

y_lim <- c(min(lcf_stat_df$lcf_integral) - 0.05, max(lcf_stat_df$lcf_integral) + 0.05)

line_width <- 1 / ggp2_magic_number 

# Beeswarm plot
lcf_ctype_beeswarm_p <- ggplot(lcf_stat_df, aes(x=ctype, y=lcf_integral, col=ctype)) +
  geom_beeswarm(alpha=0.8, corral="random", cex=2.5, size=1.2) + 
  geom_boxplot(color="black", alpha=0.6, width=0.6, outlier.shape=NA, linewidth=line_width*0.8) +
  ylab("AUC") +
  scale_y_continuous(limits=y_lim) +
  scale_x_discrete(labels=lcf_summ$x_tick) +
  scale_color_manual(values=ctype_colors) +
  geom_text(data = lcf_summ, aes(x = ctype, y = med, label = sprintf("%.2f", med)), 
            vjust = -0.22, colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = "Helvetica Neue", fontface="bold") +
  geom_text(data = lcf_summ, aes(x = ctype, y = min, label = sprintf("%.2f", min)), 
            vjust = 1.5, colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = "Helvetica Neue", fontface="bold") +
  geom_text(data = lcf_summ, aes(x = ctype, y = max, label = sprintf("%.2f", max)), 
            vjust = -0.5, colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = "Helvetica Neue", fontface="bold") +
  theme(plot.margin = margin(0, 0, 0, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth=line_width),
        axis.text = element_text(colour = "black", size = default_pointsize, family = default_font),
        axis.title = element_text(colour = "black", size = default_pointsize, family = default_font),
        text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.title.x = element_blank(),
        legend.position="none") 

lcf_ctype_beeswarm_p

lcf_stat_file_path <- file.path(output_folder, "LFC_AUC_beeswarm.pdf")

pdf_out(lcf_stat_file_path, width=2.1, height=2)
print(lcf_ctype_beeswarm_p)
dev.off()

facet_labels <- lcf_summ$x_tick
names(facet_labels) <- lcf_summ$ctype

# show pearson correlation
lcf_vs_pn_p <- ggplot(lcf_stat_df) +
  geom_point(aes(x=cell_num, y=lcf_integral, col=ctype), size=1.2) +
  geom_smooth(aes(x=cell_num, y=lcf_integral), method="loess", se=FALSE, 
              color="black", linewidth=line_width) +
  scale_color_manual(values=ctype_colors) +
  facet_wrap(~ctype, nrow=2, dir = "h", labeller = labeller(ctype = facet_labels)) +
  scale_x_log10(breaks = c(10^2, 10^3, 10^4),
                labels = c(expression(10^2), expression(10^3), expression(10^4))) +
  scale_y_continuous(expand = c(0, 0), limits=y_lim) +
  ylab("AUC") +
  xlab("Number of cells") +
  geom_text(data = lcf_summ, aes(x = 4500, y = y_lim[2] - 0.1, label = sprintf("r = %.2f", lcf_cn_corr)), 
            colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = default_font) +
  theme(plot.margin = margin(0, 0, 0, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.spacing = unit(0.3, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth=line_width),
        axis.text = element_text(colour = "black", size = default_pointsize, family = default_font),
        axis.title = element_text(colour = "black", size = default_pointsize, family = default_font),
        text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        legend.position="none")

lcf_vs_pn_p

lcf_stat_file_path <- file.path(output_folder, "LCF_AUC_vs_pn.pdf")

pdf_out(lcf_stat_file_path, width=2.1, height=2)
print(lcf_vs_pn_p)
dev.off()
   