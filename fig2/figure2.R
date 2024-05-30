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

# Clustered patterns
num_points <- 2000

clustered_color <- "#DB165C"
clustered_color_tr <- alpha(clustered_color, alpha=0.3)

# Window - square
square_side <- 1000
# TODO: get rid of the square_df and only use square_window?
square_df <- make_square_df(square_side)
square_window <- window_from_roi(square_df)
square_a <- square_side^2

# Domain radius, same for all patterns
rad <- 25
circle_a <- pi * rad^2

rmax_plot <- square_side / 4

# Display window params
# Based on the pattern with 3 clusters
wh_ratio <- 2.1 / sqrt(3)

expansion_y <- expansion(mult = 0.02, add = 0)
expansion_x <- expansion(mult = 0.005, add = 0)

line_width <- 1 / ggp2_magic_number

ppunif_from_clust_centers <- function(centers, clust_rad, clust_int, square_window) {

  # Points are distributed uniformly inside clusters
  xs <- NULL
  ys <- NULL

  for (i in 1:nrow(centers)) {
    xo <- centers$x[i]
    yo <- centers$y[i]

    circle_df <- make_circle_df(xo, yo, clust_rad)
    cluster_window <- window_from_roi(circle_df)
    pp <- rpoispp(clust_int, win = cluster_window, nsim=1)

    x <- pp[["x"]]
    xs <- c(xs, x)

    y <- pp[["y"]]
    ys <- c(ys, y)
  }

  pp <- ppp(xs, ys, window = square_window)
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
    title_y <- element_text(angle=0, vjust=0.5, size = default_pointsize, family = default_font, colour = "black")
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
    theme(plot.margin = margin(7.5, 12, 0, 2.5),
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


# Theree clusters
cluster_dist <- 120
column_distance <- sqrt(cluster_dist^2 - (cluster_dist / 2)^2)

disp_rect_width <- 3 * column_distance + 2 * rad
disp_rect_height <- disp_rect_width / wh_ratio

# Add more space around border
disp_rect_width <- disp_rect_width + 4

cluster_centers <- data.frame(x=c(-cluster_dist / 2, 0, cluster_dist / 2),
                              y=c(-column_distance/2, column_distance/2, -column_distance/2))
clust_intensity <- num_points / 3 / circle_a

three_clust_pp <- ppunif_from_clust_centers(cluster_centers, rad, clust_intensity, square_window)

# Make display rectangle
three_clust_disp_df <- make_rect_df(disp_rect_width, disp_rect_height)

# Save to file
three_clust_pp_file <- file.path(output_folder, "Three clusters.pdf")
x_lim <- c(min(three_clust_disp_df$X), max(three_clust_disp_df$X))
y_lim <- c(min(three_clust_disp_df$Y), max(three_clust_disp_df$Y))
wh_ratio <- (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1])
pdf_width <- 1.5
pdf_out(three_clust_pp_file, width = pdf_width, height = pdf_width * wh_ratio)
save_pp_as_pdf(three_clust_pp, clustered_color, three_clust_disp_df, scale=50,
               scale_offset=20, cex=0.5)
dev.off()

# Matern cluster
parent_num <- 75
parent_intensity <- parent_num / square_a
clust_pn <- num_points / parent_num

mat_clust_pp <- rMatClust(parent_intensity,
                          rad,
                          clust_pn,
                          win = square_window)

x_min <- min(square_df$X) + disp_rect_width / 2
x_max <- max(square_df$X) - disp_rect_width / 2
y_min <- min(square_df$Y) + disp_rect_height / 2
y_max <- max(square_df$Y) - disp_rect_height / 2

vis_point <- get_max_dens_point(mat_clust_pp, c(x_min, x_max), c(y_min, y_max))

# Make display rectangle
mat_clust_disp_df <- make_rect_df(disp_rect_width, disp_rect_height, origin = c(vis_point$x, vis_point$y))

# Get subset of pp in the display window
mat_clust_window <- window_from_roi(mat_clust_disp_df)

# Get subset of pp in the display window
mat_clust_pp_subset <- ppp(mat_clust_pp$x, mat_clust_pp$y, window = mat_clust_window)

# Save to file
mat_clust_pp_file <- file.path(output_folder, "Mattern Cluster.pdf")

pdf_out(mat_clust_pp_file, width = pdf_width, height = pdf_width * wh_ratio)
save_pp_as_pdf(mat_clust_pp_subset, clustered_color, mat_clust_disp_df, cex=0.5)
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

# Plot results
# Plot F and H
# Y axis limits H
h_y_lim_min <- min(h_df_three_clust$iso, h_df_mat_clust$iso)
h_y_lim_max <- max(h_df_three_clust$iso, h_df_mat_clust$iso)
h_y_lim <- c(h_y_lim_min * 1.5, h_y_lim_max + 3)

# Y axis limits PCF
pcf_y_lim_min <- min(pcf_df_three_clust$iso, pcf_df_mat_clust$iso)
pcf_y_lim_max <- max(pcf_df_three_clust$iso, pcf_df_mat_clust$iso)
pcf_y_lim <- c(-5, pcf_y_lim_max + 1)

# Plots for the pattern with 3 clusters
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

h_three_clust_p <- plot_lcf(h_df_three_clust, rad, rmax, clustered_color,
                            line_width, "H",
                            int_dists = three_clust_dist,
                            y_lim = h_y_lim,
                            show_y_axis=TRUE)


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


pdf_out(file.path(output_folder, "funcs.pdf"), width=3.3, height=3.3)
print(ggarrange(lcf_three_clust_p, lcf_matern_p,
                h_three_clust_p, h_matern_p,
                pcf_three_clust_p, pcf_matern_p,
                ncol = 2,
                bottom=textGrob("r", gp = gpar(fontsize=default_pointsize, fontfamily=default_font)),
                newpage = FALSE))
dev.off()


