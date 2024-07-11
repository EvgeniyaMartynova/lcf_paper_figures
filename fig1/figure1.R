rm(list=ls())

library(RColorBrewer)
library(dplyr)
library(RANN)
library(plotrix)
library(lcfstat)
library(spatstat)
library(ggplot2)
library(egg)

source("../settings.R")
source("../utils_lcf.R")

get_dispaly_win_origin <- function(point_vis, disp_square_side, large_radius) {
  xo <- point_vis$X - disp_square_side / 2
  yo <- min(point_vis$Y + disp_square_side /3 , point_vis$Y + 1.1 * large_radius) - disp_square_side / 2

  list(xo=xo, yo=yo)
}

make_vis_target_segms <- function(point, size) {
  segments_x0 <- c(point$X - size, point$X)
  segments_y0 <- c(point$Y, point$Y - size)
  segments_x1 <- c(point$X + size, point$X)
  segments_y1 <- c(point$Y, point$Y + size)

  return(list(x0=segments_x0, x1=segments_x1, y0=segments_y0, y1=segments_y1))
}

make_inset_segms <- function(angles, point_vis, rad_vis, point_vis_inset, rad_vis_inset) {
  get_x <- function (rad, angle, xo) rad * cos(angle) + xo
  get_y <- function (rad, angle, xo) rad * sin(angle) + yo

  xo <- point_vis$X[1]
  yo <- point_vis$Y[1]

  x0 <- NULL
  y0 <- NULL
  for (i in 1:2) {
    angle <- angles[i]
    x0 <- c(x0, get_x(rad_vis, angle, xo))
    y0 <- c(y0, get_y(rad_vis, angle, yo))
  }

  xo <- point_vis_inset$X[1]
  yo <- point_vis_inset$Y[1]

  x1 <- NULL
  y1 <- NULL
  for (i in 3:4) {
    angle <- angles[i]
    x1 <- c(x1, get_x(rad_vis_inset, angle, xo))
    y1 <- c(y1, get_y(rad_vis_inset, angle, yo))
  }

  return(list(x0=x0, x1=x1, y0=y0, y1=y1))
}

plot_radius <- function(point, rad, text,
                        text_scale, text_x_sift, text_y_sift,
                        rangle = 0) {
  rad_angle <- 0
  rad_seg_x1 <- point$X + cos(rangle) * rad
  rad_seg_y1 <- point$Y + sin(rangle) * rad

  segments(point$X, point$Y, x1 = rad_seg_x1, y1 = rad_seg_y1,
           col = "black")

  r_label_x <- point$X + cos(rangle) * text_scale + text_x_sift
  r_label_y <- point$Y + sin(rangle) * text_scale + text_y_sift
  text(r_label_x, r_label_y, text, cex=0.53, srt=rad2deg(rangle))
}

plot_full_and_zoomed_pp <- function(points,
                                    col,
                                    window,
                                    point_vis,
                                    small_radius,
                                    large_radius,
                                    zoom_ratio=0.75,
                                    angles=c(0, pi, pi / 20, 19 * pi / 20),
                                    plot_radii=FALSE,
                                    add_scale=FALSE) {

  square_side <- window$xrange[2] - window$xrange[1]
  # Scale a large radius to display in the inset
  large_radius_sc <- zoom_ratio * square_side / 2

  x_lim <- window$xrange
  # Fit the display window of the point pattern and the inset
  y_lim <- c(window$yrange[1] - large_radius_sc * 2 - 10, window$yrange[2])

  # Plot the point pattern
  par(mfrow = c(1,1), mar = c(0,0,0,0), lwd=1)
  plot(Y ~ X, data=points, pch=16, col=col, asp=1, axes=FALSE, xlab="", ylab="",
       xlim=x_lim,  ylim=y_lim, cex=0.5)
  plot(window, add = TRUE)

  # Make the point at the center of the inset larger
  points(point_vis$X, point_vis$Y, pch=16, cex=0.65, col=col)

  # Draw a target at this point
  target_size <- 2.6
  target <- make_vis_target_segms(point_vis, target_size)

  segments(target$x0, target$y0, x1 = target$x1, y1 = target$y1,
            col = "black", lwd=1)

  # Draw circles to indicate the inset
  draw.circle(point_vis$X, point_vis$Y, small_radius, lw=1)
  draw.circle(point_vis$X, point_vis$Y, large_radius, lw=1)

  # Plot the magnified inset
  # Extract the points indise the inset
  encl_circle_win <- disc(large_radius, c(point_vis$X, point_vis$Y))
  pp_encl_subset <- ppp(points$X, points$Y, window = encl_circle_win)
  points_sub <- data.frame(X=pp_encl_subset$x, Y=pp_encl_subset$y)

  # Scale a small radius
  scalar <- large_radius_sc / large_radius
  small_radius_sc <- small_radius * scalar

  # The coordinates of the central point in the magnified inset
  point_vis_inset <- data.frame(X=point_vis$X[1], Y = window$yrange[1] - large_radius_sc - 10)

  # Compute the sift of the coordinates of other points in the inset
  x_sh <- point_vis_inset$X[1] - scalar * point_vis$X[1]
  y_sh <- point_vis_inset$Y[1] - scalar * point_vis$Y[1]

  # Scale and shift the points
  points_sub_sc <- points_sub %>%
    mutate(X = scalar * X + x_sh,
           Y = scalar * Y + y_sh)

  # Plots the points in the magnified inset
  points(Y ~ X, data=points_sub_sc, pch=16, cex=0.8, col=col)

  # Enlarge the central point
  points(Y ~ X, data=point_vis_inset, pch=16, cex=1, col=col)

  # Draw a target at the central point
  target_size <- target_size * scalar
  target <- make_vis_target_segms(point_vis_inset, target_size)

  segments(target$x0, target$y0, x1 = target$x1, y1 = target$y1,
           col = "black", lwd=1)

  # Also draw radii
  if (plot_radii) {
    # Small radius
    plot_radius(point_vis_inset, small_radius_sc, "r",
                small_radius_sc / 2,
                cos(pi/3) * 5, sin(pi/3) * 5)

    # Large radius
    plot_radius(point_vis_inset, large_radius_sc, "hr",
                (small_radius_sc + large_radius_sc) / 2,
                cos(pi/4) * 2, sin(pi/4) * 7,
                - pi / 9)
  }

  # Draw magnified circles
  draw.circle(point_vis_inset$X, point_vis_inset$Y, small_radius_sc, lw=1)
  draw.circle(point_vis_inset$X, point_vis_inset$Y, large_radius_sc, lw=1)

  # Draw lines to connect
  inset_segms <- make_inset_segms(angles, point_vis, large_radius, point_vis_inset, large_radius_sc)
  segments(inset_segms$x0, inset_segms$y0, x1 = inset_segms$x1, y1 = inset_segms$y1)

  # Add scale
  if (add_scale) {
    scale_offset <- 5
    text_width <- 25
    text_height <- 9
    plot_scale(small_radius, scale_offset, x_lim, y_lim, 3, text_width, text_height, TRUE)
  }
}

save_pp_with_zoom <- function(pp, point_vis, pdf_name, point_color, plot_radii=FALSE, add_scale=FALSE) {

  # Build the square around this point
  win_origin <- get_dispaly_win_origin(point_vis, disp_square_side, large_radius)
  disp_square_window <- make_square_win(disp_square_side, c(win_origin$xo, win_origin$yo))

  # Extract points in the plotting window
  pp_plot_subset <- ppp(pp$x, pp$y, window = disp_square_window)
  points_plot <- data.frame(X=pp_plot_subset$x, Y=pp_plot_subset$y)

  # Full visualization of point pattern subset and F-func ROI
  pp_file <- file.path(output_folder, pdf_name)
  pdf_out(pp_file, width = pdf_width, height = pdf_height)
  plot_full_and_zoomed_pp(points_plot, point_color, disp_square_window, point_vis,
                          small_radius, large_radius, zoom_ratio, plot_radii=plot_radii, add_scale=add_scale)
  dev.off()
}

plot_summary_functon <- function(sf_df, x_breaks, x_labels, y_lim, sf_name) {

  rmax <- max(sf_df$r)

  sf_p <- ggplot(sf_df) +
    geom_segment(aes(x = 0, y = 0, xend = rmax, yend = 0, col = "black"),
                 linewidth=line_width) +
    geom_line(aes(x=r, y=iso, col=type), linewidth=line_width) +
    geom_vline(xintercept = small_radius, linewidth=line_width, col="black", linetype="dashed") +
    geom_vline(xintercept = small_radius * 2, linewidth=line_width, col="black", linetype="dashed") +
    scale_colour_manual(values = c(dispersed_color, random_color, clustered_color, "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits=c(0, rmax), expand = stand_expansion) +
    scale_y_continuous(limits=y_lim, expand = stand_expansion) +
    theme_bw(base_size = default_pointsize, base_family = default_font) +
    ylab(sf_name) +
    xlab("r") +
    theme(plot.margin = margin(3, 8, 0, 0),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
          axis.line.y = element_line(colour = "black", linewidth=line_width),
          axis.line.x = element_line(colour = "black", linewidth=line_width),
          axis.ticks.y = element_line(colour = "black", linewidth=line_width),
          axis.ticks.x = element_line(colour = "black", linewidth=line_width),
          legend.position="none")

  sf_p
}

seed <- 1645
set.seed(seed)

output_folder <- "img"
dir.create(output_folder, showWarnings = FALSE)

sim_square_side <- 400
sim_square_window <- make_square_win(sim_square_side)
sim_square_area <- sim_square_side^2

# Parameters of a window to zoom in
disp_square_side <- 100

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

# Pick visualization circle radii
# For visualization purposes only
h <- 2

# Radii to visualise
small_radius <- 15
large_radius <- small_radius * sqrt(h)
# Size of the inset area relative to the display window's size
zoom_ratio <- 0.65

# Compute the width and height of PDFs with saved plots
height_plot <- disp_square_side + disp_square_side * zoom_ratio + 10
hw_ratio <- height_plot / disp_square_side
pdf_width <- 1.043307
pdf_height <- pdf_width * hw_ratio

# Pick an anchor point for the choice of a magnified region 1/3 y and 1/2 x
xt <- sim_square_side / 2
yt <- sim_square_side / 2 - disp_square_side / 6

# Generate point patterns
# Random
point_num <- 6000
intensity <- point_num / sim_square_area
rand_pp <- rpoispp(intensity, win = sim_square_window)

# Pick a point to emphasize
rand_points <- data.frame(X=rand_pp$x, Y=rand_pp$y)
point_vis <- closest_point(rand_points, xt, yt)

save_pp_with_zoom(rand_pp, point_vis, "Rand_full.pdf", random_color)

# Clustered
avg_parent_num <- 50
parent_int <- avg_parent_num / sim_square_area
avg_child_num <- point_num / avg_parent_num
clust_rad <- 7.5

clust_pp <- rThomas(parent_int,
                    clust_rad,
                    avg_child_num,
                    win = sim_square_window,
                    saveparents = TRUE)

# Pick a point to emphasize
parents <- attr(clust_pp, "parents")
parents_df <- data.frame(X=parents$x, Y=parents$y)
point_vis <- closest_point(parents_df, xt + small_radius, yt)

save_pp_with_zoom(clust_pp, point_vis, "Clust_full.pdf", clustered_color, add_scale=TRUE)

# Dispersed
strauss_pp <- sim_strauss_mh(point_num, 0.4, small_radius, sim_square_window, 1)

# Pick a point to emphasize
points_dispersed <- data.frame(X=strauss_pp$x, Y=strauss_pp$y)
point_vis <- closest_point(points_dispersed, xt, yt)

save_pp_with_zoom(strauss_pp, point_vis, "Disp_full.pdf", dispersed_color, plot_radii=TRUE)

 # Plot summary functions
rmax <- sim_square_side / 4
x_breaks <- c(0, small_radius, 2 * small_radius, rmax / 2, rmax)
x_labels <- as.integer(round(x_breaks))

# ggplot2 parameters
line_width <- 1 / ggp2_magic_number
stand_expansion <- expansion(mult = 0.005, add = 0)

# Calculate LCF for point patterns
dim_lims <- c(10, 50)

random_lcf <- LCFest(rand_pp, rmax=rmax, dim_lims = dim_lims)
random_lcf_df <- as.data.frame(random_lcf) %>%
  dplyr::select(-theo) %>%
  mutate(type="random")

clustered_lcf <- LCFest(clust_pp, rmax=rmax, dim_lims = dim_lims)
clustered_lcf_df <- as.data.frame(clustered_lcf) %>%
  dplyr::select(-theo) %>%
  mutate(type="clustered")

dispersed_lcf <- LCFest(strauss_pp, rmax=rmax, dim_lims = dim_lims)
dispersed_lcf_df <- as.data.frame(dispersed_lcf) %>%
  dplyr::select(-theo) %>%
  mutate(type="dispersed")

lcf_df <- rbind(random_lcf_df, clustered_lcf_df, dispersed_lcf_df) %>%
  mutate(type = factor(type, levels=c("dispersed", "random", "clustered")))

# Plot LCF
lcf_ylim <- c(-1, 1)
lcf_p <- plot_summary_functon(lcf_df, x_breaks, x_labels, lcf_ylim, "LCF")

pdf_out(file.path(output_folder, "LCF_plot.pdf"), width=1.6, height=1.15)
print(lcf_p)
dev.off()

# Calculate H function for point patterns
random_h_df <- h_func(rand_pp, rmax=rmax, correction="Ripley")
random_h_df <- random_h_df %>%
  dplyr::select(-theo) %>%
  mutate(type="random")

clustered_h_df <- h_func(clust_pp, rmax=rmax, correction="Ripley")
clustered_h_df <- clustered_h_df %>%
  dplyr::select(-theo) %>%
  mutate(type="clustered")

dispersed_h_df <- h_func(strauss_pp, rmax=rmax, correction="Ripley")
dispersed_h_df <- dispersed_h_df %>%
  dplyr::select(-theo) %>%
  mutate(type="dispersed")

h_df <- rbind(random_h_df, clustered_h_df, dispersed_h_df) %>%
  mutate(type = factor(type, levels=c("dispersed", "random", "clustered")))

# Plot H
h_ylim <- c(min(dispersed_h_df$iso), max(clustered_h_df$iso))
h_p <- plot_summary_functon(h_df, x_breaks, x_labels, h_ylim, "H")

pdf_out(file.path(output_folder, "H-func.pdf"), width=1.6, height=1.15)
print(h_p)
dev.off()

# Plot together
pdf_out(file.path(output_folder, "funcs.pdf"), width=3.3, height=1.15)
ggarrange(h_p, lcf_p, ncol = 2, nrow=1, newpage = FALSE)
dev.off()

