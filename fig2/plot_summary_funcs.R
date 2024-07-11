rm(list=ls())

library(readr)
library(ggplot2)
library(egg)
library(ggmagnify)

source("../settings.R")

plot_sf <- function(sf_df,
                    fun_name,
                    breaks_y = waiver(),
                    labels_y = waiver(),
                    ylim=NULL,
                    add_xaxis = FALSE,
                    add_facet_labels = FALSE,
                    h_lines = NULL,
                    y_lab_vjust = NULL) {

  if (is.null(ylim)) {
    ylim <- c(min(sf_df$sf), max(sf_df$sf))
  }

  breaks_x <- c(0, 0.5, 1)
  labels_x <- c("0", "0.5", "1")

  if (add_xaxis) {
    x_lab <- "r"
    text_x <- element_text(size = default_pointsize, family = default_font, colour = "black")
    ticks_x <- element_line(colour = "black", linewidth=line_width)
    line_x <- element_line(colour = "black", linewidth=line_width)
  } else {
    x_lab <- ""
    text_x <- element_blank()
    ticks_x <- element_blank()
    line_x <- element_blank()
  }

  if (add_facet_labels) {
    strip_text <- element_text(size = default_pointsize, family = default_font, colour = "black")
  } else {
    strip_text <- element_blank()
  }

  square_labeller <- label_bquote(cols = "|W|="~.(as.numeric(levels(square))[square]))

  sf_p <- ggplot(sf_df) +
    facet_wrap(~square, scales = "free_x", strip.position = "top", labeller = square_labeller, nrow=1) +
    geom_line(aes(x=r, y=sf, col=pattern), linewidth=line_width) +
    scale_colour_manual(values = c(clustered_color, dispersed_color, random_color)) +
    scale_x_continuous(breaks = breaks_x, labels=labels_x, expand = stand_expansion) +
    scale_y_continuous(breaks = breaks_y, labels=labels_y, limits=ylim, expand = stand_expansion) +
    theme_bw(base_size = default_pointsize, base_family = default_font) +
    ylab(fun_name) +
    xlab(x_lab) +
    theme(plot.margin = margin(0, 2, 0, 4),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_blank(),
          strip.text.x = strip_text,
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_text(vjust = y_lab_vjust),
          axis.text.y = element_text(size = default_pointsize, family = default_font,
                                     colour = "black"),
          axis.text.x = text_x,
          axis.line.y = element_line(colour = "black", linewidth=line_width),
          axis.line.x = line_x,
          axis.ticks.y = element_line(colour = "black", linewidth=line_width),
          axis.ticks.x = ticks_x,
          legend.position="none")

  if (!is.null(h_lines)) {
    for (h_line in h_lines) {
      sf_p <- sf_p +
        geom_hline(yintercept = h_line, linewidth=line_width * 0.75, col="gray")
    }
  }

  sf_p
}

img_folder <- "img"
data_folder <- "data"

# ggplot2 parameters
line_width <- 1 / ggp2_magic_number
stand_expansion <- expansion(mult = 0.01, add = 0)

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

# Plot K
k_file <- file.path(data_folder, "k.csv")
k_df <- read_csv(k_file) %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))

# Make y limits large enough to fit the inset
k_ylim <- c(-0.02, 4)
k_p <- plot_sf(k_df, "K", breaks_y = c(0, 1, 2, 4), labels_y = c("0.0", "1.0", "2.0", "4.0"),
               ylim = k_ylim, add_facet_labels=TRUE)

# Add squared inset to the plot with the smallest observation window
magnify_data <- data.frame(square = factor(c("1"), levels = c("1", "2", "4")))

# Make an area to magnify on the plot (bottom lest square)
inset_min <- -0.015
inset_max <- 0.408
magnify_data$from <- list(c(inset_min, inset_max, inset_min, inset_max))

# The ratio of the inset size on the plot and ymax
range_share <- 0.575

# Make an area where a magnified plot is shown
# We want the inset to be in the middle of x scale and at the top of y scale
ymax <- 4
ymin <- ymax - ymax * range_share

# X is depends on the observation window side length, 1 for the smallest one
rmax <- sqrt(1 / pi)
xmin <- rmax * (1 - range_share) / 2
xmax <- xmin + rmax * range_share

magnify_data$to <- list(c(xmin, xmax, ymin, ymax))

# Add an inset to the plot
k_p_ins <- k_p +
  scale_colour_manual(values = c(clustered_color, dispersed_color, random_color)) +
  geom_magnify(mapping=aes(from = from, to = to), shape = "rect",
               data=magnify_data, expand = 0, proj.linetype=2,
               colour="black", linewidth=line_width / 2)
k_p_ins

# Plot normalized H
h1_file <- file.path(data_folder, "h1.csv")
h1_df <- read_csv(h1_file) %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))

ylim_h1 <- c(min(h1_df$sf), 1)

h1_p <- plot_sf(h1_df, "H1", breaks_y = c(0, 0.5, 1), labels_y = c("0.0", "0.5", "1.0"),
                ylim = ylim_h1)
h1_p

# Plot PCF
pcf_file <- file.path(data_folder, "pcf.csv")
pcf_df <- read_csv(pcf_file) %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))

# For maximally clustered pattern, PCF's estimate is huge at small distances
# We truncate PCF that is greased than 1.3 * maximum value for the maximally dispersed pattern
max_md <- max((pcf_df %>% filter(pattern=="md"))$sf)

pcf_df_plot <- pcf_df %>%
  mutate(plot = case_when(pattern=="mc" & sf > max_md * 1.3 ~ FALSE,
                          TRUE ~ TRUE)) %>%
  filter(plot) %>%
  dplyr::select(-plot)

pcf_p <- plot_sf(pcf_df_plot, "PCF", breaks_y = c(0, 2, 4, 6), labels_y = c("0.0", "2.0", "4.0", "6.0"))
pcf_p

# Plot LCF
lcf_file <- file.path(data_folder, "lcf.csv")
lcf_df <- read_csv(lcf_file) %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))

lcf_p <- plot_sf(lcf_df, "LCF", breaks_y = c(-1, 0, 1), labels_y = c("-1.0", "0.0", "1.0"),
                 add_xaxis = TRUE, y_lab_vjust = -0.5)
lcf_p

# Plot together
# Inches
pdf_width <- 3.3
facet_height <- 0.9
facet_space <- 0.02
labels_height <- 0.03
pdf_height <- 4 * facet_height + 4 * facet_space + labels_height

pdf_out(file.path(img_folder, "sfs_plot.pdf"), width=pdf_width, height=pdf_height)
ggarrange(k_p_ins, h1_p, pcf_p, lcf_p,
          ncol = 1, newpage = FALSE)
dev.off()
