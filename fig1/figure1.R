rm(list=ls())

library(RColorBrewer)
library(dplyr)
library(RANN)
library(egg)
library(plotrix)
library(lcfstat)
library(spatstat)
library(ggplot2)

source("../settings.R")
source("../utils_lcf.R")

seed <- 1645
set.seed(seed)

output_folder <- "img"
dir.create(output_folder, showWarnings = FALSE)

sim_square_side <- 400
sim_square_df <- make_square_df(sim_square_side)
sim_square_window <- window_from_roi(sim_square_df)
sim_square_area <- sim_square_side^2

# Parameters of a window to zoom in
disp_square_side <- 100
disp_square_area <- disp_square_side^2

# Make a display window to zoom in
disp_square_df <- make_square_df(disp_square_side)
disp_square_window <- window_from_roi(disp_square_df)

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

# Pick visualization circle radii
# For visualization purposes only
k <- 2

small_radius <- 15
large_radius <- small_radius * sqrt(k)
zoom_ratio <- 0.65

# Pick a  1/3 y and 1/2 x
xt <- 0
yt <- disp_square_window$yrange[1] + (disp_square_window$yrange[2] - disp_square_window$yrange[1]) /3

# Functions for visualisation of the point patterns
plot_full_pp <- function(points, col, window_df, point_vis, small_radius, large_radius) {
  x_lim <- c(min(window_df$X), max(window_df$X))
  y_lim <- c(min(window_df$Y), max(window_df$Y))

  par(mfrow = c(1,1), mar = c(0,0,0,0))
  plot(points$X, points$Y, pch=16, col=col, asp=1, axes=FALSE, xlab="", ylab="",
       xlim=x_lim,  ylim=y_lim)
  polygon(window_df$X, window_df$Y)
  draw.circle(point_vis$X, point_vis$Y, small_radius)
  draw.circle(point_vis$X, point_vis$Y, large_radius)
}

plot_zoomed_pp <- function(points, col, point_vis, small_radius, large_radius) {
  x_lim <- c(point_vis$X - large_radius, point_vis$X + large_radius)
  y_lim <- c(point_vis$Y - large_radius, point_vis$Y + large_radius)

  par(mfrow = c(1,1), mar = c(0,0,0,0))
  plot(points$X, points$Y, pch=16, cex=2, col=col, asp=1, axes=FALSE, xlab="", ylab="",
       xlim=x_lim,  ylim=y_lim)
  draw.circle(point_vis$X, point_vis$Y, small_radius)
  draw.circle(point_vis$X, point_vis$Y, large_radius)
}

make_vis_target_segms <- function(point, offset) {
  segments_x0 <- c(point$X - offset, point$X)
  segments_y0 <- c(point$Y, point$Y - offset)
  segments_x1 <- c(point$X + offset, point$X)
  segments_y1 <- c(point$Y, point$Y + offset)

  return(list(x0=segments_x0, x1=segments_x1, y0=segments_y0, y1=segments_y1))
}

plot_full_and_zoomed_pp <- function(points,
                                    col,
                                    window_df,
                                    point_vis,
                                    small_radius,
                                    large_radius,
                                    zoom_ratio=0.75,
                                    angles=c(0, pi, pi / 20, 19 * pi / 20),
                                    plot_radii=FALSE,
                                    add_scale=FALSE) {

  large_radius_sc <- zoom_ratio * disp_square_side / 2

  height <- max(window_df$Y) - min(window_df$Y)
  x_lim <- c(min(window_df$X), max(window_df$X))
  y_lim <- c(min(window_df$Y) - large_radius_sc * 2 - 10, max(window_df$Y))

  par(mfrow = c(1,1), mar = c(0,0,0,0), lwd=1)
  plot(points$X, points$Y, pch=16, col=col, asp=1, axes=FALSE, xlab="", ylab="",
       xlim=x_lim,  ylim=y_lim, cex=0.5)
  polygon(window_df$X, window_df$Y)

  points(point_vis$X, point_vis$Y, pch=16, cex=0.65, col=col)

  point_vis_offset <- 2.6
  segments_vis <- make_vis_target_segms(point_vis, point_vis_offset)

  segments(segments_vis$x0, segments_vis$y0, x1 = segments_vis$x1, y1 = segments_vis$y1,
            col = "black", lwd=1)

  draw.circle(point_vis$X, point_vis$Y, small_radius, lw=1)
  draw.circle(point_vis$X, point_vis$Y, large_radius, lw=1)

  # Add scale
  if (add_scale) {
    scale_offset <- 5
    text_height <- 9
    seg_x0 <- x_lim[2] - scale_offset - small_radius
    seg_x1 <- seg_x0 + small_radius
    # For Y axis offset should be negative
    seg_y0 <- y_lim[2] - scale_offset
    seg_y1 <- seg_y0

    segments(seg_x0, seg_y0, x1 = seg_x1, y1 = seg_y1,
             col = "black", lwd=3, lend=2)

    label_x <- seg_x0 - 25
    label_y <- seg_y0 - text_height
    text(label_x, label_y, as.expression(bquote(bold(.(as.character(small_radius)))~bold(units))), cex=0.675, adj = 0)
  }

  # Vis ROI
  encl_circle_win <- disc(large_radius, c(point_vis$X, point_vis$Y))
  pp_encl_subset <- ppp(points$X, points$Y, window = encl_circle_win)
  points_sub <- data.frame(X=pp_encl_subset$x, Y=pp_encl_subset$y)

  scalar <- large_radius_sc / large_radius
  small_radius_sc <- small_radius * scalar

  point_vis_new <- data.frame(X=point_vis$X[1], Y=min(window_df$Y) - large_radius_sc - 10)

  x_sh <- point_vis_new$X[1] - scalar * point_vis$X[1]
  y_sh <- point_vis_new$Y[1] - scalar * point_vis$Y[1]

  points_sub_sc <- points_sub %>%
    mutate(X=scalar * X + x_sh, Y=scalar * Y + y_sh)

  points(points_sub_sc$X, points_sub_sc$Y, pch=16, cex=0.8, col=col)

  points(point_vis_new$X, point_vis_new$Y, pch=16, cex=1, col=col)

  point_vis_offset <- point_vis_offset * scalar
  segments_vis_new <- make_vis_target_segms(point_vis_new, point_vis_offset)

  segments(segments_vis_new$x0, segments_vis_new$y0, x1 = segments_vis_new$x1, y1 = segments_vis_new$y1,
           col = "black", lwd=1)

  # Also draw radii
  if (plot_radii) {
    # Small radius
    rad_angle <- 0
    rad_seg_x1 <- point_vis_new$X + cos(rad_angle) * small_radius_sc
    rad_seg_y1 <- point_vis_new$Y + sin(rad_angle) * small_radius_sc

    segments(point_vis_new$X, point_vis_new$Y, x1 = rad_seg_x1, y1 = rad_seg_y1,
             col = "black")

    r_label_x <- point_vis_new$X + cos(rad_angle) * small_radius_sc / 2 + cos(pi/3) * 5
    r_label_y <- point_vis_new$Y + sin(rad_angle) * small_radius_sc / 2 + sin(pi/3) * 5
    text(r_label_x, r_label_y, "r", cex=0.53, srt=0)

    # Large radius
    rad_angle <- - pi / 9
    rad_seg_x1 <- point_vis_new$X + cos(rad_angle) * large_radius_sc
    rad_seg_y1 <- point_vis_new$Y + sin(rad_angle) * large_radius_sc

    segments(point_vis_new$X, point_vis_new$Y, x1 = rad_seg_x1, y1 = rad_seg_y1,
             col = "black")

    r_label_x <- point_vis_new$X + cos(rad_angle) * (small_radius_sc + large_radius_sc) / 2 + cos(pi/4) *2
    r_label_y <- point_vis_new$Y + sin(rad_angle) * (small_radius_sc + large_radius_sc) / 2 + sin(pi/4) * 7
    text(r_label_x, r_label_y, "hr", cex=0.53, srt=-15)
  }

  # Draw circles
  draw.circle(point_vis_new$X, point_vis_new$Y, small_radius_sc, lw=1)
  draw.circle(point_vis_new$X, point_vis_new$Y, large_radius_sc, lw=1)

  # Draw lines to connect
  inner_left_angle <- angles[1]
  inner_left_x <- large_radius * cos(inner_left_angle) + point_vis$X[1]
  inner_left_y <- large_radius * sin(inner_left_angle) + point_vis$Y[1]

  inner_right_angle <- angles[2]
  inner_right_x <- large_radius * cos(inner_right_angle) + point_vis$X[1]
  inner_right_y <- large_radius * sin(inner_right_angle) + point_vis$Y[1]

  outer_left_angle <- angles[3]
  outer_left_x <- large_radius_sc * cos(outer_left_angle) + point_vis_new$X[1]
  outer_left_y <- large_radius_sc * sin(outer_left_angle) + point_vis_new$Y[1]

  outer_right_angle <- angles[4]
  outer_right_x <- large_radius_sc * cos(outer_right_angle) + point_vis_new$X[1]
  outer_right_y <- large_radius_sc * sin(outer_right_angle) + point_vis_new$Y[1]

  segments(c(inner_left_x, inner_right_x), c(inner_left_y, inner_right_y),
           c(outer_left_x, outer_right_x), c(outer_left_y, outer_right_y))
}


# Random
point_num <- 6000
intensity <- point_num / sim_square_area

rand_pp <- rpoispp(intensity, win = sim_square_window)
points_rand <- data.frame(X=rand_pp$x, Y=rand_pp$y)

# Pick a point to emphasize
point_vis <- closest_point(points_rand, xt, yt)
disp_square_xo <- point_vis$X
disp_square_yo <- min(point_vis$Y + (disp_square_window$yrange[2] - disp_square_window$yrange[1]) /3 , point_vis$Y + 1.1*large_radius)

# Build the square around this point
disp_square_df <- make_square_df(disp_square_side, c(disp_square_xo, disp_square_yo))
disp_square_window <- window_from_roi(disp_square_df)

# Extract points in the plotting window
pp_plot_subset <- ppp(rand_pp$x, rand_pp$y, window = disp_square_window)
rand_points_plot <- data.frame(X=pp_plot_subset$x, Y=pp_plot_subset$y)

# Full visualisation of point pattern subset and F-func ROI
# Compute this only once, maybe move the code above
height_win <- max(sim_square_df$Y) - min(sim_square_df$Y)
width_win <- max(sim_square_df$X) - min(sim_square_df$X)

height_plot <- height_win + height_win * zoom_ratio + 20
hw_ratio <- height_plot / width_win

pdf_width <- 1.043307
pdf_height <- pdf_width * hw_ratio

rand_pp_file <- file.path(output_folder, "Rand_full.pdf")
pdf_out(rand_pp_file, width = pdf_width, height = pdf_height)
plot_full_and_zoomed_pp(rand_points_plot, random_color, disp_square_df, point_vis, small_radius, large_radius, zoom_ratio)
dev.off()

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

parents <- attr(clust_pp, "parents")
parents_df <- data.frame(X=parents$x, Y=parents$y)

# Pick a point to emphasize
point_vis <- closest_point(parents_df, xt + small_radius, yt)
disp_square_xo <- point_vis$X
disp_square_yo <- min(point_vis$Y + (disp_square_window$yrange[2] - disp_square_window$yrange[1]) /3, point_vis$Y + 1.1 * large_radius)

# Build the square around this point
disp_square_df <- make_square_df(disp_square_side, c(disp_square_xo, disp_square_yo))
disp_square_window <- window_from_roi(disp_square_df)

# Extract points in the plotting window
pp_plot_subset <- ppp(clust_pp$x, clust_pp$y, window = disp_square_window)
clust_points_plot <- data.frame(X=pp_plot_subset$x, Y=pp_plot_subset$y)

# Full visualisation of point pattern subset and F-func ROI
clust_pp_file <- file.path(output_folder, "Clust_full.pdf")
pdf_out(clust_pp_file, width = pdf_width, height = pdf_height)
plot_full_and_zoomed_pp(clust_points_plot, clustered_color, disp_square_df, point_vis,
                        small_radius, large_radius, zoom_ratio, add_scale=TRUE)
dev.off()

# Dispersed
strauss_pp <- sim_strauss_mh(point_num, 0.4, small_radius, sim_square_window, 1)
strauss_pp$n

points_dispersed <- data.frame(X=strauss_pp$x, Y=strauss_pp$y)

# Pick a point to emphasize
point_vis <- closest_point(points_dispersed, xt, yt)
disp_square_xo <- point_vis$X
disp_square_yo <- min(point_vis$Y + disp_square_side / 3, point_vis$Y + 1.1*large_radius)

# Make display window
disp_square_df <- make_square_df(disp_square_side, c(disp_square_xo, disp_square_yo))
disp_square_window <- window_from_roi(disp_square_df)

# Extract points in the plotting window
pp_plot_subset <- ppp(points_dispersed$X, points_dispersed$Y, window = disp_square_window)
disp_points_plot <- data.frame(X=pp_plot_subset$x, Y=pp_plot_subset$y)

rand_pp_file <- file.path(output_folder, "Disp_full.pdf")
pdf_out(rand_pp_file, width = pdf_width, height = pdf_height, family="Helvetica Neue Medium")
plot_full_and_zoomed_pp(disp_points_plot, dispersed_color, disp_square_df, point_vis, small_radius, large_radius, zoom_ratio, plot_radii=TRUE)
dev.off()

 # Plot everything together
rmax <- sim_square_side / 4
dim_lims <- c(10, 50)
x_breaks <- c(0, small_radius, 2*small_radius, rmax / 2, rmax)
x_labels <- as.integer(round(x_breaks))

# Calculate F function for point patterns
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

spec_lines_df <- data.frame(x1=0, x2=rmax, y_rand=0)

line_width <- 1 / ggp2_magic_number
stand_expansion <- expansion(mult = 0.005, add = 0)

lcf_p <- ggplot(lcf_df) +
  geom_segment(aes(x = x1, y = y_rand, xend = x2, yend = y_rand, col = "black"),
               data = spec_lines_df, linewidth=line_width) +
  geom_line(aes(x=r, y=iso, col=type), linewidth=line_width) +
  geom_vline(xintercept = small_radius, linewidth=line_width, col="black", linetype="dashed") +
  geom_vline(xintercept = small_radius * 2, linewidth=line_width, col="black", linetype="dashed") +
  scale_colour_manual(values = c(dispersed_color, random_color, clustered_color, "black")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits=c(0, rmax), expand = stand_expansion) +
  scale_y_continuous(limits=c(-1, 1), expand = stand_expansion) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  ylab("LCF") +
  xlab("r") +
  theme(plot.margin = margin(5, 10, 5, 2.5),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.line.y = element_line(colour = "black", linewidth=line_width),
        axis.line.x = element_line(colour = "black", linewidth=line_width),
        axis.title.y = element_text(angle=0, vjust=0.5),
        axis.ticks.y = element_line(colour = "black", linewidth=line_width),
        axis.ticks.x = element_line(colour = "black", linewidth=line_width),
        legend.position="none")

pdf_out(file.path(output_folder, "LCF_plot.pdf"), width=1.6, height=1.6)
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

h_ylim <- c(min(dispersed_h_df$iso), max(clustered_h_df$iso))

spec_lines_df <- data.frame(x1=0, x2=rmax, y_rand=0)

h_p <- ggplot(h_df) +
  geom_segment(aes(x = x1, y = y_rand, xend = x2, yend = y_rand, col = "black"),
               data = spec_lines_df, linewidth=line_width) +
  geom_line(aes(x=r, y=iso, col=type), linewidth=line_width) +
  geom_vline(xintercept = small_radius, linewidth=line_width, col="black", linetype="dashed") +
  geom_vline(xintercept = small_radius * 2, linewidth=line_width, col="black", linetype="dashed") +
  scale_colour_manual(values = c(dispersed_color, random_color, clustered_color, "black")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits=c(0, rmax), expand = stand_expansion) +
  scale_y_continuous(limits=h_ylim, expand = stand_expansion) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  ylab("H") +
  xlab("r") +
  theme(plot.margin = margin(5, 10, 5, 2.5),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.line.y = element_line(colour = "black", linewidth=line_width),
        axis.line.x = element_line(colour = "black", linewidth=line_width),
        axis.title.y = element_text(angle=0, vjust=0.5),
        axis.ticks.y = element_line(colour = "black", linewidth=line_width),
        legend.position="none")


pdf_out(file.path(output_folder, "H-func.pdf"), width=1.6, height=1.6)
print(h_p)
dev.off()


pdf_out(file.path(output_folder, "funcs.pdf"), width=3.3, height=1.6)
print(ggarrange(h_p, lcf_p, ncol = 2, nrow=1, newpage = FALSE))
dev.off()

