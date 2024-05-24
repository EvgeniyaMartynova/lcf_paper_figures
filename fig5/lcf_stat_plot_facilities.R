rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggbeeswarm)
library(lme4)

source("../settings.R")

# File names
data_folder <- "data"
stat_file_name <- "lcf_stat_air.csv"
stat_file_path <- file.path(data_folder, stat_file_name)

output_folder <- "img/facilities"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

air_color <- "#A63603" 
waste_color <- "#08519C" 

# Facility type
ftype_colors <- c(air_color, waste_color)

lcf_stat_df <- read.csv(stat_file_path)

lcf_stat_df <- lcf_stat_df %>%
  mutate(ftype = as.factor(ftype))

lcf_stat_df <- lcf_stat_df %>%
  mutate(ftype = case_when(ftype == "air" ~ "Air",
                            ftype == "waste" ~ "Waste")) %>%
  mutate(ftype = factor(ftype, levels=c("Air", "Waste")))

lcf_summ <- lcf_stat_df %>% 
  group_by(ftype) %>% 
  summarise(med = median(lcf_integral), 
            min = min(lcf_integral), 
            max = max(lcf_integral), 
            n = n(),
            lcf_cn_corr = cor(lcf_integral, log10(facility_num))) %>%
  mutate(x_tick = sprintf("n=%d", n))

y_lim <- c(min(lcf_stat_df$lcf_integral) - 0.05, max(lcf_stat_df$lcf_integral) + 0.05)

line_width <- 1 / ggp2_magic_number

# Beeswarm plot
lcf_ctype_beeswarm_p <- ggplot(lcf_stat_df, aes(x=ftype, y=lcf_integral, col=ftype)) +
  geom_beeswarm(alpha=0.8, corral="random", cex=2.5, size=1.2) + 
  geom_boxplot(color="black", alpha=0.6, width=0.5, outlier.shape=NA, linewidth=line_width*0.8) +
  ylab("AUC") +
  scale_y_continuous(limits=y_lim) +
  scale_x_discrete(labels=lcf_summ$x_tick) +
  scale_color_manual(values=ftype_colors) +
  geom_text(data = lcf_summ, aes(x = ftype, y = med, label = sprintf("%.2f", med)), 
            vjust = -0.22, colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = "Helvetica Neue", fontface="bold") +
  geom_text(data = lcf_summ, aes(x = ftype, y = min, label = sprintf("%.2f", min)), 
            vjust = 1.5, colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = "Helvetica Neue", fontface="bold") +
  geom_text(data = lcf_summ, aes(x = ftype, y = max, label = sprintf("%.2f", max)), 
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

lcf_stat_file_path <- file.path(output_folder, "LCF_AUC_beeswarm.pdf")

pdf_out(lcf_stat_file_path, width=1.5, height=2)
print(lcf_ctype_beeswarm_p)
dev.off()

facet_labels <- lcf_summ$x_tick
names(facet_labels) <- lcf_summ$ftype

# show pearson correlation
lcf_vs_pn_p <- ggplot(lcf_stat_df) +
  geom_point(aes(x=facility_num, y=lcf_integral, col=ftype), size=1.2) +
  geom_smooth(aes(x=facility_num, y=lcf_integral), method='loess', se=FALSE, 
              color="black", linewidth=line_width) +
  scale_color_manual(values=ftype_colors) +
  facet_wrap(~ftype, nrow=2, dir = "h", labeller = labeller(ftype = facet_labels)) +
  scale_x_log10(breaks = c(10^2, 10^3),
                labels = c(expression(10^2), expression(10^3))) +
  scale_y_continuous(expand = c(0, 0), limits=y_lim, breaks = c(0.2, 0.4, 0.6)) +
  ylab("AUC") +
  xlab("Number of facilities") +
  geom_text(data = lcf_summ, aes(x = 1500, y = y_lim[2] - 0.1, label = sprintf("r = %.2f", lcf_cn_corr)), 
            colour="black", size = (default_pointsize - 2) * gtext_magic_number, 
            family = default_font) +
  theme(plot.margin = margin(0, 5, 0, 2),
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

pdf_out(lcf_stat_file_path, width=1.5, height=2)
print(lcf_vs_pn_p)
dev.off()


