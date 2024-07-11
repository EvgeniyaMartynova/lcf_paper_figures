rm(list=ls())

library(lcfstat)
library(spatstat)
library(dplyr)
library(MASS)
library(egg)
library(RColorBrewer)
library(RANN)
library(grid)

source("../settings.R")
source("../utils_lcf.R")

seed <- 1901
set.seed(seed)

# Output folder
output_folder <- "img"
dir.create(output_folder, showWarnings = FALSE)

clustered_color <- "#DB165C"
clustered_color_tr <- alpha(clustered_color, alpha=0.3)

# Window - square
square_side <- 1000
square_window <- make_square_win(square_side)
square_a <- square_side^2

# Average number of points in all point patterns
num_points <- 2000

# Domain radius, same for all patterns
rad <- 25
circle_a <- pi * rad^2

rmax_plot <- square_side / 4

# Display window parameters
# Empirically found based on the pattern with 3 clusters
hw_ratio <- 0.73
pdf_width <- 1.5
pdf_height <- pdf_width * hw_ratio

# ggplot2 parameters
expansion_y <- expansion(mult = 0.02, add = 0)
expansion_x <- expansion(mult = 0.005, add = 0)
line_width <- 1 / ggp2_magic_number

ppunif_from_clust_centers <- function(centers, clust_rad, clust_int, square_window) {

  # Make clusters
  pps <- vector(mode = "list", length = 0)
  for (i in 1:nrow(centers)) {
    xo <- centers$x[i]
    yo <- centers$y[i]

    # Points are distributed uniformly inside clusters
    circle_df <- make_circle_df(xo, yo, clust_rad)
    cluster_window <- window_from_roi(circle_df)
    pp <- rpoispp(clust_int, win = cluster_window, nsim=1)

    pps[[i]] <- pp
  }

  pps[["W"]] <- square_window

  pp <- do.call("superimpose", pps)
  pp
}

get_max_dens_point <- function(pp, x_lim=NULL, y_lim=NULL, n=50) {

  kde <- kde2d(pp$x, pp$y, n = n)

  if (!is.null(x_lim)) {
    inds_x_vis <- which((kde$x > x_lim[1]) & (kde$x < x_lim[2]))
  } else {
    inds_x_vis <- 1:length(kde$x)
  }

  if (!is.null(y_lim)) {
    inds_y_vis <- which((kde$y > y_lim[1]) & (kde$y < y_lim[2]))
  } else {
    inds_y_vis <- 1:length(kde$y)
  }

  kde_z_vis <- kde$z[inds_x_vis, inds_y_vis]

  max_inds_vis <- which(kde_z_vis == max(kde_z_vis), arr.ind = TRUE)

  argmax_x_vis <- as.numeric(max_inds_vis[1, "row"])
  argmax_y_vis <- as.numeric(max_inds_vis[1, "col"])

  argmax_x <- inds_x_vis[argmax_x_vis]
  argmax_y <- inds_y_vis[argmax_y_vis]

  x_max <- kde$x[argmax_x]
  y_max <- kde$y[argmax_y]
  list(x=x_max, y=y_max)
}

plot_lcf <- function(df, domain_rad, rmax, col,
                     lwd, y_lab,
                     int_dists=NULL,
                     x_breaks=waiver(),
                     x_labels=waiver(),
                     y_random=0,
                     y_lim=NULL,
                     y_breaks=waiver(),
                     show_x_axis=FALSE,
                     show_y_axis=FALSE) {

  rmax_plot <- max(df$r)

  if (is.null(x_breaks) && show_x_axis) {
    x_breaks <- c(0, domain_rad * 2, rmax/2, rmax)
  }

  if (is.null(x_labels) && show_x_axis) {
    x_labels <- as.integer(round(x_breaks))
  }

  lcf_p <- ggplot(df) +
    geom_hline(yintercept = y_random, col = "gray", linewidth=lwd) +
    geom_line(aes(x=r, y=iso, col=col), linewidth=lwd) +
    scale_colour_manual(values = c(col, "gray")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits=c(0, rmax_plot), expand = expansion_x) +
    scale_y_continuous(breaks = y_breaks, limits=y_lim, expand = expansion_y) +
    theme_bw(base_size = default_pointsize, base_family = default_font)

  if (show_x_axis) {
    line_x <- element_line(colour = "black", linewidth=lwd)
    ticks_x <- element_line(colour = "black", linewidth=lwd)
    text_x <- element_text(size = default_pointsize, family = default_font, colour = "black")
  } else {
    line_x <- element_blank()
    ticks_x <- element_blank()
    text_x <- element_blank()
  }

  if (show_y_axis) {
    lcf_p <- lcf_p +
      ylab(y_lab)

    line_y <- element_line(colour = "black", linewidth=lwd)
    ticks_y <- element_line(colour = "black", linewidth=lwd)
    title_y <- element_text(size = default_pointsize, family = default_font, colour = "black")
    text_y <- element_text(size = default_pointsize, family = default_font, colour = "black")
  } else {
    line_y <- element_blank()
    ticks_y <- element_blank()
    title_y <- element_blank()
    text_y <- element_blank()
  }

  if (!is.null(int_dists)) {
    for (dist in int_dists) {
      lcf_p <- lcf_p +
        geom_vline(xintercept = dist, linewidth=lwd, col="black", linetype="dashed")
    }
  }

  lcf_p <- lcf_p +
    theme(plot.margin = margin(8, 4, 0, 2),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.y = text_y,
          axis.title.y = title_y,
          axis.line.y = line_y,
          axis.ticks.y = ticks_y,
          axis.text.x = text_x,
          axis.title.x = element_blank(),
          axis.line.x = line_x,
          axis.ticks.x = ticks_x,
          legend.position="none")

  lcf_p
}

# Three clusters
cluster_dist <- 120
column_distance <- sqrt(cluster_dist^2 - (cluster_dist / 2)^2)

# Make a point pattern with three clusters
cluster_centers <- data.frame(x=c(-cluster_dist / 2, 0, cluster_dist / 2) + square_side / 2,
                              y=c(-column_distance/2, column_distance/2, -column_distance/2) + square_side / 2)
clust_intensity <- num_points / 3 / circle_a

three_clust_pp <- ppunif_from_clust_centers(cluster_centers, rad, clust_intensity, square_window)

# Compute the display limits
padding <- 75
disp_rect_width <- cluster_dist + 2 * rad + 2 * padding
disp_rect_height <- disp_rect_width * hw_ratio

ox <- (square_side - disp_rect_width) / 2
disp_x_lim <- c(0, disp_rect_width) + ox

oy <- (square_side - disp_rect_height) / 2
disp_y_lim <- c(0, disp_rect_height) + oy

# Save to file
three_clust_pp_file <- file.path(output_folder, "Three clusters.pdf")
pdf_out(three_clust_pp_file, width = pdf_width, height = pdf_height)
save_pp_as_pdf(three_clust_pp, clustered_color, disp_x_lim, disp_y_lim,
               scale=50, scale_offset=12.5, cex=0.5)
dev.off()

# Matern cluster
parent_num <- 75
parent_intensity <- parent_num / square_a
clust_pn <- num_points / parent_num

mat_clust_pp <- rMatClust(parent_intensity,
                          rad,
                          clust_pn,
                          win = square_window)

# Find the location with the highest point density to center the display window at
# Restrict the search to the central area of the window in such a way that
# any selected central point gives the display window that fully lies in the main window
x_min <- square_window$xrange[1] + disp_rect_width / 2
x_max <- square_window$xrange[2] - disp_rect_width / 2
y_min <- square_window$yrange[1] + disp_rect_height / 2
y_max <- square_window$yrange[2] - disp_rect_height / 2

vis_point <- get_max_dens_point(mat_clust_pp, c(x_min, x_max), c(y_min, y_max))

# Make display rectangle
ox <- vis_point$x - disp_rect_width / 2
disp_x_lim <- c(0, disp_rect_width) + ox

oy <- vis_point$y - disp_rect_height / 2
disp_y_lim <- c(0, disp_rect_height) + oy

# Get a subset of a point pattern for visualization
mat_clust_pp_sub <- subset(mat_clust_pp, (x >= disp_x_lim[1] & x <= disp_x_lim[2]
                                          & y >= disp_y_lim[1] & y <= disp_y_lim[2]))

# Save to file
mat_clust_pp_file <- file.path(output_folder, "Mattern Cluster.pdf")
pdf_out(mat_clust_pp_file, width = pdf_width, height = pdf_height)
save_pp_as_pdf(mat_clust_pp_sub, clustered_color, disp_x_lim, disp_y_lim, cex=0.5)
dev.off()

# Compute H, PCF and LCF
# Triangle
h_df_three_clust <- h_func(three_clust_pp, correction = "Ripley", rmax=rmax_plot)
pcf_df_three_clust <- pc_func(three_clust_pp, correction = "Ripley", rmax=rmax_plot)
lcf_df_three_clust <- LCFest(three_clust_pp, rmax=rmax_plot, dim_lims = c(10, 50))

# Matern cluster process
h_df_mat_clust <- h_func(mat_clust_pp, correction = "Ripley", rmax=rmax_plot)
pcf_df_mat_clust <- pc_func(mat_clust_pp, correction = "Ripley", rmax=rmax_plot)
lcf_df_mat_clust <- LCFest(mat_clust_pp, rmax=rmax_plot, dim_lims = c(10, 50))

# Plots for the pattern with 3 clusters
# LCF
three_clust_dist <- c(50, 70, 120, 170)
x_labels <- c("50  ", "  70", "120", "170")

lcf_three_clust_p <- plot_lcf(lcf_df_three_clust, rad, rmax, clustered_color,
                              line_width, "LCF",
                              int_dists = three_clust_dist,
                              y_lim = c(-1, 1),
                              y_breaks = c(-1, 0, 1),
                              show_y_axis=TRUE)

# A weird hack to fix the saving of the funcs.pdf plot
pdf_out(file.path(output_folder, "LCF_three.pdf"), width=1.6, height=1.6)
print(lcf_three_clust_p)
dev.off()

# H
h_y_lim_min <- min(h_df_three_clust$iso, h_df_mat_clust$iso)
h_y_lim_max <- max(h_df_three_clust$iso, h_df_mat_clust$iso)
h_y_lim <- c(h_y_lim_min * 1.5, h_y_lim_max + 3)

h_three_clust_p <- plot_lcf(h_df_three_clust, rad, rmax, clustered_color,
                            line_width, "H",
                            int_dists = three_clust_dist,
                            y_lim = h_y_lim,
                            show_y_axis=TRUE)


# PCF
pcf_y_lim_min <- min(pcf_df_three_clust$iso, pcf_df_mat_clust$iso)
pcf_y_lim_max <- max(pcf_df_three_clust$iso, pcf_df_mat_clust$iso)
pcf_y_lim <- c(-5, pcf_y_lim_max + 1)

pcf_three_clust_p <- plot_lcf(pcf_df_three_clust, rad, rmax, clustered_color,
                              line_width, "PCF",
                              y_random = 1,
                              int_dists = three_clust_dist,
                              x_breaks =  three_clust_dist,
                              x_labels = x_labels,
                              y_lim = pcf_y_lim,
                              show_x_axis=TRUE,
                              show_y_axis=TRUE)


# Plots for Mater cluster pattern
mattern_dist <- c(50)

lcf_matern_p <- plot_lcf(lcf_df_mat_clust, rad, rmax, clustered_color,
                        line_width, "LCF",
                        int_dists = mattern_dist,
                        y_lim = c(-1, 1))

h_matern_p <- plot_lcf(h_df_mat_clust, rad, rmax, clustered_color,
                       line_width, "H",
                       int_dists = mattern_dist,
                       y_lim = h_y_lim)

pcf_matern_p <- plot_lcf(pcf_df_mat_clust, rad, rmax, clustered_color,
                         line_width, "PCF",
                         y_random = 1,
                         int_dists = mattern_dist,
                         x_breaks =  mattern_dist,
                         y_lim = pcf_y_lim,
                         show_x_axis=TRUE)

# Combine all plots
row_height <- 0.92
row_space <- 0.02
height <- row_height * 3 + 2 * row_space

pdf_out(file.path(output_folder, "funcs.pdf"), width=3.3, height=height)
ggarrange(lcf_three_clust_p, lcf_matern_p,
          h_three_clust_p, h_matern_p,
          pcf_three_clust_p, pcf_matern_p,
          ncol = 2,
          bottom=textGrob("r", gp = gpar(fontsize=default_pointsize, fontfamily=default_font)),
          newpage = FALSE)
dev.off()


