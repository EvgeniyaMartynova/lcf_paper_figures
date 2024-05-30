rm(list=ls())

library(dplyr)
library(ggplot2)

source("../settings.R")

h_color <- "#D6984D"
pcf_color <- "#008F0C"
lcf_color <- "#4A16DB"

file_path_template <- "data/Calibration_S%s_PN%s.csv"

square_size <- c(141, 200, 282)

parameters <- list(list(square_size=141, point_nums=c(25, 50, 100)),
                   list(square_size=200, point_nums=c(25, 50, 100)),
                   list(square_size=282, point_nums=c(25, 50, 100)))


breaks <- c(10^(-3), 10^(-2), 10^(-1), 10^0)
axes_labels <- log10(breaks)

calibration_df_glob <- NULL

for (param in parameters) {
  square_size <- param$square_size
  point_nums <- param$point_nums

  for (i in 1:length(point_nums)) {
    point_num <- point_nums[i]
    file_path <- sprintf(file_path_template, square_size, point_num)

    calibration_df <- read.csv(file_path)
    calibration_df[["s"]] <- square_size
    calibration_df[["pn"]] <- point_num
    calibration_df_glob <- rbind(calibration_df_glob, calibration_df)
  }
}

calibration_df_glob <- calibration_df_glob %>%
  mutate(stat = case_when(stat == "h" ~ "H",
                          stat == "pcf" ~ "PCF",
                          stat == "lcf" ~ "LCF")) %>%
  mutate(stat = factor(stat, levels=c("H", "PCF", "LCF"))) %>%
  mutate(area_base = case_when(s == 141 ~ 2,
                          s == 200 ~ 4,
                          s == 282 ~ 8)) %>%
  mutate(area = case_when(s == 141 ~ "|W|",
                          s == 200 ~ "2|W|",
                          s == 282 ~ "4|W|"))


labels_df <- calibration_df_glob %>%
  dplyr::select(area, pn, area_base) %>%
  distinct(area, pn, .keep_all = TRUE) %>%
  mutate(area_pow = 4) %>%
  mutate(dens = pn / area_base / 10)

labels_df$x <- 0.0015
labels_df$y1 <- 0.4

my_exp <- as.character(expression(gamma ~ "=" ~ 0.25))

labels_df <- labels_df %>%
  mutate(label= as.character(dens))

x_lim <- c(0.001, 1)
y_lim <- c(0.001, 1)

grid_labeller <- label_bquote(cols = .(pn)~"points", rows= .(area))
line_width <- 1 / ggp2_magic_number

expansion_xy <- expansion(mult = 0.01, add = 0)

calibration_p <- ggplot(calibration_df_glob) +
  geom_line(aes(x=sl, y=mean, col=stat), linewidth=line_width) +
  geom_ribbon(aes(x=sl, ymin = low, ymax = high, fill=stat), alpha = 0.2) +
  geom_text(aes(x=x, y=y1, label = label, hjust = 0), data = labels_df, parse = T,
            size=(default_pointsize - 1) * gtext_magic_number, family=default_font ) +
  scale_colour_manual(values = c(h_color, pcf_color, lcf_color)) +
  scale_fill_manual(values = c(h_color, pcf_color, lcf_color)) +
  geom_abline(intercept=0, slope=1, colour="black", linewidth=line_width)  +
  facet_grid(area ~ pn, labeller=grid_labeller, scales = "free_x") +
  scale_x_log10(breaks = breaks, labels = axes_labels, expand = expansion_xy, limits=x_lim) +
  scale_y_log10(breaks = breaks, labels = axes_labels,  expand = expansion_xy,  limits=y_lim) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  ylab(expression(log[10]~"Type I error")) +
  xlab(expression(log[10]~"Significance level")) +
  theme(plot.margin = margin(0, 0, 0.5, 0.5),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        panel.spacing.x = unit(0.7 , "lines"),
        panel.spacing.y = unit(0.7, "lines"),
        strip.background = element_blank(),
        strip.placement = "inside",
        text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        strip.text.x = element_text(size = default_pointsize, family = default_font, colour = "black"),
        strip.text.y = element_text(angle = 0, vjust=0.5, size = default_pointsize, family = default_font, colour = "black"),
        axis.text.x = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.text.y = element_text(size = default_pointsize, family = default_font, colour = "black"),
        legend.title=element_blank(),
        axis.line = element_line(colour = "black", linewidth=line_width),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position="top")

calibration_p

pdf_out("figure3.pdf", width=3.3, height=3)
print(calibration_p)
dev.off()
