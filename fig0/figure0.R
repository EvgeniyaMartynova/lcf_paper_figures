rm(list=ls())

library(RColorBrewer)
library(dplyr)
library(lcfstat)
library(spatstat)
library(ggplot2)
library(ggmagnify)

source("../settings.R")
source("../utils_lcf.R")


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

  square_labeller <- label_bquote(cols = "|W|="~.(as.numeric(levels(square_s))[square_s]))

  sf_p <- ggplot(sf_df) +
    facet_wrap(~square_s, scales = "free_x", strip.position = "top", labeller = square_labeller, nrow=1) +
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

set.seed(2906)

img_folder <- "img"
dir.create(img_folder, showWarnings = FALSE)

data_folder <- "data"
dir.create(data_folder, showWarnings = FALSE)

square_sides <- c(1, sqrt(2), 2)

clust_radius <- 0.05
disp_distance <- 0.15

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

width <- 3.3
height <- 0.9

line_width <- 1 / ggp2_magic_number
stand_expansion <- expansion(mult = 0.01, add = 0)

expansion_x <- expansion(mult = 0.01, add = 0)
expansion_y <- expansion(mult = 0.02, add = 0)

# Point patterns parameters
# Maximum clustering
num_points <- 200
# Generate random pattern in a cluster domain
circle <- disc(radius = clust_radius)
circle_a <- spatstat.geom::area(circle)
intensity_mc <- num_points / circle_a
pp_loc <- rpoispp(intensity_mc, win=circle)

# Random
intensity_rand <- 200

# Summary functions' parameters
correction <- "Ripley"
dim_lims <- c(4, 50)

k_df <- NULL
pcf_df <- NULL
lcf_df <- NULL

steps_in <- 513

lcf_mc_dims <- c(30, 40, 50)
lcf_md_dims <- c(30, 45, 60)

for (i in 1:length(square_sides)) {
  square_side <- square_sides[i]
  # Increase the number of bins for r
  steps_r <- round(steps_in * square_side)

  square_df <- make_square_df(square_side)
  square_window <- window_from_roi(square_df)
  square_area <- square_side^2

  # Change window to be the target domain
  pp_mc <- ppp(pp_loc$x, pp_loc$y, window = square_window)

  # Random
  pp_random <- rpoispp(intensity_rand, win=square_window)

  # Maximum dispersion
  column_distance <- sqrt(disp_distance^2 - (disp_distance / 2)^2)
  grid_pp <- max_disp_pp(square_window, disp_distance)

  rmax <- sqrt(square_area/pi)

  # K
  k_mc <- Kest(pp_mc, correction = correction, r = seq(0, rmax, length.out=steps_r))
  k_rand <- Kest(pp_random, correction = correction, r = seq(0, rmax, length.out=steps_r))
  k_md <- Kest(grid_pp, correction = correction, r = seq(0, rmax, length.out=steps_r))

  k_df_ss <- data.frame(r=k_mc$r, sf = c(k_mc$iso, k_rand$iso, k_md$iso),
                        pattern = c(rep("mc", nrow(k_mc)), rep("rand", nrow(k_rand)), rep("md", nrow(k_md))),
                        square=square_side^2)

  k_df <- rbind(k_df, k_df_ss)

  # PCF
  pcf_mc <- pcf(pp_mc, correction = correction, r = seq(0, rmax, length.out=steps_r))
  pcf_rand <- pcf(pp_random, correction = correction, r = seq(0, rmax, length.out=steps_r))
  pcf_md <- pcf(grid_pp, correction = correction, r = seq(0, rmax, length.out=steps_r))

  # Remove infinity for plotting
  pcf_mc_iso <- c(pcf_mc$iso[2], pcf_mc$iso[-1])
  pcf_rand_iso <- c(pcf_rand$iso[2], pcf_rand$iso[-1])
  pcf_md_iso <- c(pcf_md$iso[2], pcf_md$iso[-1])

  pcf_df_ss <- data.frame(r=pcf_mc$r, sf = c(pcf_mc_iso, pcf_rand_iso, pcf_md_iso),
                          pattern = c(rep("mc", nrow(pcf_mc)), rep("rand", nrow(pcf_rand)), rep("md", nrow(pcf_md))),
                          square=square_side^2)

  pcf_df <- rbind(pcf_df, pcf_df_ss)

  # LCF
  dim_mc <- lcf_mc_dims[i]
  lcf_mc <- LCFest(pp_mc, correction = correction, r = seq(0, rmax, length.out=steps_r), dim=dim_mc)
  lcf_rand <- LCFest(pp_random, correction = correction, r = seq(0, rmax, length.out=steps_r))
  dim_md <- lcf_md_dims[i]
  lcf_md <- LCFest(grid_pp, correction = correction, r = seq(0, rmax, length.out=steps_r), dim=dim_md)

  lcf_df_ss <- data.frame(r=lcf_mc$r, sf = c(lcf_mc$iso, lcf_rand$iso, lcf_md$iso),
                          pattern = c(rep("mc", nrow(lcf_mc)), rep("rand", nrow(lcf_rand)), rep("md", nrow(lcf_md))),
                          square=square_side^2)

  lcf_df <- rbind(lcf_df, lcf_df_ss)

  # Visualize the point patterns
  if (i == 1) {
    pdf_side <- 1.05
    # Change window to be the target domain
    pdf_out(file.path(img_folder, "Max_clust.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(pp_mc, clustered_color, square_df, scale=0.1,
                   scale_offset=0.05, scale_lwd = 3, text_height=0.1, cex=0.5)
    dev.off()

    # Random
    pdf_out(file.path(img_folder, "Rand.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(pp_random, random_color, square_df, cex=0.5)
    dev.off()

    # Maximum dispersion
    pdf_out(file.path(img_folder, "Max_disp.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(grid_pp, dispersed_color, square_df, cex=0.5)
    dev.off()
  }
}

# Make plot of K-function
k_df <- k_df %>%
  mutate(pattern = as.factor(pattern),
         square_s = factor(square, levels=c("1", "2", "4")))

write.csv(k_df, file.path(data_folder, "k.csv"))

k_ylim <- c(-0.02, 4)

k_p <- plot_sf(k_df, "K", breaks_y = c(0, 1, 2, 4), labels_y = c("0.0", "1.0", "2.0", "4.0"),
               ylim = k_ylim, add_facet_labels=TRUE)

# Add squared inset
magnify_data <- data.frame(square_s = factor(c("1"), levels = c("1", "2", "4")))

# Magnify the same area for each observation window
inset_min <- -0.015
inset_max <- 0.408

magnify_data$from <- list(c(inset_min, inset_max, inset_min, inset_max))

# Make approximately squared target area, same y, x depends on the rmax
range_share <- 0.575

# We want the inset to be in the middle of x scale and at the bottom of y scale
# Y is the same for all observation windows, manually picked ymin that looks good
ymax <- 4
ymin <- ymax - 4 * range_share

# X is depends on the observation window size
rmax <- square_sides[1] * sqrt(1 / pi)
xmin <- rmax * (1 - range_share) / 2
xmax <- xmin + rmax * range_share

magnify_data$to <- list(c(xmin, xmax, ymin, ymax))

k_p_ins <- k_p +
  scale_colour_manual(values = c(clustered_color, dispersed_color, random_color)) +
  geom_magnify(mapping=aes(from = from, to = to), shape = "rect",
               data=magnify_data, expand = 0, proj.linetype=2,
               colour="black", linewidth=line_width / 2)
k_p_ins

# Normalized H
h1_df <- k_df %>%
  mutate(sf = sqrt(sf / pi) / sqrt(square / pi) - r / sqrt(square / pi))

write.csv(h1_df, file.path(data_folder, "h1.csv"))

ylim_h1 <- c(min(h1_df$sf), 1)

h1_p <- plot_sf(h1_df, "H1", breaks_y = c(0, 0.5, 1), labels_y = c("0.0", "0.5", "1.0"),
                ylim = ylim_h1)
h1_p

# PCF
pcf_df <- pcf_df %>%
  mutate(pattern = as.factor(pattern),
         square_s = as.factor(square))

write.csv(pcf_df, file.path(data_folder, "pcf.csv"))

# For maximally clustered pattern, PCF's estimate is huge
# we do not want to show it on the plot
max_md <- max((pcf_df %>% filter(pattern=="md"))$sf)

pcf_df_plot <- pcf_df %>%
  mutate(plot = case_when(pattern=="mc" & sf > max_md * 1.3 ~ FALSE,
                          TRUE ~ TRUE)) %>%
  filter(plot) %>%
  dplyr::select(-plot)

pcf_p <- plot_sf(pcf_df_plot, "PCF", breaks_y = c(0, 2, 4, 6), labels_y = c("0.0", "2.0", "4.0", "6.0"))
pcf_p

# LCF
lcf_df <- lcf_df %>%
  mutate(pattern = as.factor(pattern),
         square_s = as.factor(square))

write.csv(lcf_df, file.path(data_folder, "lcf.csv"))

lcf_p <- plot_sf(lcf_df, "LCF", breaks_y = c(-1, 0, 1), labels_y = c("-1.0", "0.0", "1.0"),
                 add_xaxis = TRUE, y_lab_vjust = -0.5)
lcf_p

# Plot together
# Inches
facet_height <- 0.9
facet_space <- 0.02
labels_height <- 0.03
height <- 4 * facet_height + 4 * facet_space + labels_height

pdf_out(file.path(img_folder, "sfs_plot.pdf"), width=width, height=height)
ggarrange(k_p_ins, h1_p, pcf_p, lcf_p,
          ncol = 1, newpage = FALSE)
dev.off()
