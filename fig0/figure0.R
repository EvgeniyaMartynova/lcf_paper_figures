rm(list=ls())

library(RColorBrewer)
library(dplyr)
library(lcfstat)
library(spatstat)
library(ggplot2)
library(ggmagnify)
library(grid)

source("../settings.R")
source("../utils_lcf.R")


plot_sf <- function(sf_df, fun_name) {
  
  sf_ylim <- c(min(sf_df$sf), max(sf_df$sf))
  
  breaks_x <- c(0, 0.5, 1)
  labels_x <- c("0", "0.5", "1")
  
  sf_p <- ggplot(sf_df) +
    facet_wrap(~square, scales = "free_x", strip.position = "top", nrow=1) +
    geom_line(aes(x=r, y=sf, col=pattern), linewidth=line_width) +
    scale_colour_manual(values = c(clustered_color, dispersed_color, random_color)) +
    scale_x_continuous(breaks = breaks_x, labels=labels_x, expand = stand_expansion) +
    scale_y_continuous(limits=sf_ylim, expand = stand_expansion) +
    theme_bw(base_size = default_pointsize, base_family = default_font) +
    ylab(fun_name) +
    xlab("r") +
    theme(plot.margin = margin(5, 5, 5, 2.5),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
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
  
  sf_p
}

set.seed(2706)

img_folder <- "img"
dir.create(img_folder, showWarnings = FALSE)

square_sides <- c(1, sqrt(2), 2)

clust_radius <- 0.05
disp_distance <- 0.15

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

width <- 3.3
height <- 1.2

line_width <- 1 / ggp2_magic_number
stand_expansion <- expansion(mult = 0.01, add = 0)

expansion_x <- expansion(mult = 0.01, add = 0)
expansion_y <- expansion(mult = 0.02, add = 0)

# Point patterns parameters
# Maximum clustering
num_points <- 500
# Generate random pattern in a cluster domain
circle <- disc(radius = clust_radius)
circle_a <- spatstat.geom::area(circle)
intensity_mc <- num_points / circle_a
pp_loc <- rpoispp(intensity_mc, win=circle)

# Random 
intensity_rand <- 500

# Summary functions' parameters
correction <- "Ripley"
dim_lims <- c(4, 50)

k_df <- NULL
pcf_df <- NULL
lcf_df <- NULL

steps_in <- 513

lcf_mc_dims <- c(30, 40, 50)
#lcf_md_dims <- c(12, 18, 25)
lcf_md_dims <- c(30, 45, 60)
i <- 1

for (square_side in square_sides) {
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
  #lcf_md <- LCFest(grid_pp, correction = correction, r = seq(0, rmax, length.out=steps_r))
  
  lcf_df_ss <- data.frame(r=lcf_mc$r, sf = c(lcf_mc$iso, lcf_rand$iso, lcf_md$iso), 
                          pattern = c(rep("mc", nrow(lcf_mc)), rep("rand", nrow(lcf_rand)), rep("md", nrow(lcf_md))), 
                          square=square_side^2)
  
  lcf_df <- rbind(lcf_df, lcf_df_ss)
  
  # Update parameters
  i <- i + 1
}


if (FALSE) {
  plot(pcf_mc)
  plot(pcf_rand)
  
  plot(pcf_md)
  
  plot(lcf_md)
  abline(v=0.15, col=2)
  abline(v=2*column_distance, col=2)
  abline(v=0.3, col=2)
  abline(v=0.398, col=2)
  abline(v=0.45, col=2)
  abline(v=4*column_distance, col=2)
  abline(v=0.54, col=2)
}

if (FALSE) {
# Visualization of the point patterns
square_side <- square_sides[1]

square_df <- make_square_df(square_side)
square_window <- window_from_roi(square_df)
square_area <- square_side^2

pdf_width <- 1.05
# Change window to be the target domain
pp_mc <- ppp(pp_loc$x, pp_loc$y, window = square_window)
pdf_out(file.path(img_folder, "Max_clust.pdf"), width = pdf_width, height = pdf_width)
save_pp_as_pdf(pp_mc, clustered_color, square_df, scale=0.1, scale_offset=0.075, text_height=0.1, cex=0.5)
dev.off()

# Random
pp_random <- rpoispp(intensity_rand, win=square_window)
pdf_out(file.path(img_folder, "Rand.pdf"), width = pdf_width, height = pdf_width)
save_pp_as_pdf(pp_random, random_color, square_df, cex=0.5)
dev.off()

# Maximum dispersion
column_distance <- sqrt(disp_distance^2 - (disp_distance / 2)^2)
grid_pp <- max_disp_pp(square_window, disp_distance)

pdf_out(file.path(img_folder, "Max_disp.pdf"), width = pdf_width, height = pdf_width)
save_pp_as_pdf(grid_pp, dispersed_color, square_df, cex=0.5)
dev.off()
}

# White
w_df <- k_df %>%
  mutate(sf = sqrt(sf / pi) / sqrt(square / pi) - r / sqrt(square / pi)) %>%
  mutate(pattern = as.factor(pattern),
         square = factor(square, levels=c("1", "2", "4")))

# K
k_df <- k_df %>%
  mutate(pattern = as.factor(pattern),
         square = factor(square, levels=c("1", "2", "4")))


k_p <- plot_sf(k_df, "K")
k_p

if (FALSE) {
  pdf_out(file.path(img_folder, "K_plot.pdf"), width=width, height=height)
  print(k_p)
  dev.off()
}

# H
h_df <- k_df %>%
  mutate(sf = sqrt(sf / pi) - r)


h_p <- plot_sf(h_df, "H")
h_p

if (FALSE) {
  pdf_out(file.path(img_folder, "H_plot.pdf"), width=width, height=height)
  print(h_p)
  dev.off()
}

# White
w_p <- plot_sf(w_df, "W")
w_p

if (FALSE) {
  pdf_out(file.path(img_folder, "White_plot.pdf"), width=width, height=height)
  print(w_p)
  dev.off()
}


# PCF
pcf_df <- pcf_df %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))

max_md <- max((pcf_df %>% filter(pattern=="md"))$sf)

pcf_df1 <- pcf_df %>%
  mutate(plot = case_when(pattern=="mc" & sf > max_md * 1.3 ~ FALSE,
                          TRUE ~ TRUE)) %>%
  filter(plot) %>%
  select(-plot)

pcf_p <- plot_sf(pcf_df, "PCF")
pcf_p

if (FALSE) {
  pdf_out(file.path(img_folder, "PCF_plot.pdf"), width=width, height=height)
  print(pcf_p)
  dev.off()
}

# LCF
lcf_df <- lcf_df %>%
  mutate(pattern = as.factor(pattern),
         square = as.factor(square))


lcf_p <- plot_sf(lcf_df, "LCF")
lcf_p

if (FALSE) {
  pdf_out(file.path(img_folder, "LCF_plot.pdf"), width=width, height=height)
  print(lcf_p)
  dev.off()
}

# Plot together 

k_df["func"] <- "K"
w_df["func"] <- "H*"
pcf_df1["func"] <- "PCF"
lcf_df["func"] <- "LCF"

sfs_df <- rbind(k_df, w_df, pcf_df1, lcf_df)
sfs_df <- sfs_df %>%
  mutate(func = factor(func, levels = c("K", "H*", "PCF", "LCF")))

# Save the dataframe as .csv
write.csv(sfs_df, "sfs_df.csv")

grid_labeller <- label_bquote(cols = "|W|"~"="~.(as.character(square)), rows= .(func))

breaks_x <- c(0, 0.5, 1)
labels_x <- c("0", "0.5", "1")

dummy_h <- data.frame(r = rep(0.1, 3), sf = rep(1, 3),
                      pattern = rep("mc", 3), square= factor(c(1,2,4)), 
                      func=factor("H*", levels = c("K", "H*", "PCF", "LCF")))

dummy_k <- data.frame(r = rep(0.1, 3), sf = rep(-3.5, 3),
                      pattern = rep("mc", 3), square= factor(c(1,2,4)), 
                      func=factor("K", levels = c("K", "H*", "PCF", "LCF")))

dummy <- rbind(dummy_h, dummy_k)

magnify_data <- data.frame(func = factor(c("K", "K", "K"), levels = c("K", "H*", "PCF", "LCF")),
                           square = factor(c("1", "2", "4"), levels = c("1", "2", "4")))

magnify_data$from <- list(c(0, 0.405, -0.01, 0.405), c(0, 0.405, -0.01, 0.405), c(0, 0.405, -0.01, 0.405))
magnify_data$to <- list(c(0.106, 0.458, -0.75, -3.25), c(0.15, 0.648, -0.75, -3.25), c(0.216, 0.917, -0.75, -3.25))

sf_p <- ggplot(sfs_df) +
  facet_grid(func~square, labeller=grid_labeller, scales = "free", switch = "y") +
  geom_line(aes(x=r, y=sf, col=pattern), linewidth=line_width) +
  geom_point(aes(x=r, y=sf, col="white"), data=dummy) +
  scale_colour_manual(values = c(clustered_color, dispersed_color, random_color, "white")) +
  scale_x_continuous(breaks = breaks_x, labels=labels_x, expand = expansion_x) +
  scale_y_continuous(expand = expansion_y) +
  theme_bw(base_size = default_pointsize, base_family = default_font) +
  ylab("") +
  xlab("r") +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.7, "lines"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = default_pointsize, family = default_font, colour = "black"),
        strip.text.y = element_text(angle = 45, vjust=0, size = default_pointsize, family = default_font, colour = "black"),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = default_pointsize, family = default_font, colour = "black"),
        axis.line.y = element_line(colour = "black", linewidth=line_width),
        axis.line.x = element_line(colour = "black", linewidth=line_width),
        axis.title.y = element_text(angle=0, vjust=0.5),
        axis.ticks.y = element_line(colour = "black", linewidth=line_width),
        legend.position="none") +
  geom_magnify(mapping=aes(from = from, to = to), shape = "rect",
               data=magnify_data, expand = 0, proj.linetype=2,
               colour="black", proj = "corresponding", linewidth=line_width / 2) 

sf_p

coef <- (4 + 3.5) / 4

gt <- ggplot_gtable(ggplot_build(sf_p))
gt$heights[10] <- coef * gt$heights[10]
grid.draw(gt)

facet_height <- 0.9
facet_space <- 0.02
height <- facet_height * coef + 3 * facet_height + 4 * facet_space

pdf_out(file.path(img_folder, "sfs_plot_adj.pdf"), width=width, height=height)
grid.draw(gt)
dev.off()

library(gtable)
gtable_show_layout(gt)


pdf_out(file.path(img_folder, "sfs_plot.pdf"), width=width, height=4.175)
print(sf_p)
dev.off()

 # Try combining plots with


if (FALSE) {
  disp_distance <- 0.15
  
  square_side <- 2
  square_df <- make_square_df(square_side)
  square_window <- window_from_roi(square_df)
  square_area <- square_side^2
  
  grid_pp <- max_disp_pp(square_window, disp_distance)
  plot(grid_pp)
  grid_pp$n
  
  rmax <- sqrt(square_area/pi)
  steps_r <- round(steps_in * square_side)
  
  lcf_md <- LCFest(grid_pp, correction = correction, r = seq(0, rmax, length.out=steps_r), dim=18)
  plot(lcf_md)
}

# 13, 18, 25