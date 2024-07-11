rm(list=ls())

library(RColorBrewer)
library(readr)
library(dplyr)
library(lcfstat)
library(spatstat)

source("../settings.R")
source("../utils_lcf.R")

set.seed(2906)

img_folder <- "img"
dir.create(img_folder, showWarnings = FALSE)

data_folder <- "data"
dir.create(data_folder, showWarnings = FALSE)

# Plotting parameters
random_color <- "#4EACD7"
clustered_color <- "#DB165C"
dispersed_color <- "#DBCB16"

# Point patterns parameters
# Maximum clustering
num_points <- 200
# Generate random pattern in a cluster domain
clust_radius <- 0.05
circle <- disc(radius = clust_radius)
circle_a <- spatstat.geom::area(circle)
intensity_mc <- num_points / circle_a
pp_loc <- rpoispp(intensity_mc, win=circle)

# Random
intensity_rand <- 200

# Maximally dispersed
disp_distance <- 0.15

# Summary functions' parameters
correction <- "Ripley"
dim_lims <- c(4, 50)

steps_in <- 513

lcf_mc_dims <- c(30, 40, 50)
lcf_md_dims <- c(30, 45, 60)

square_sides <- c(1, sqrt(2), 2)

k_df <- NULL
pcf_df <- NULL
lcf_df <- NULL

for (i in seq_along(square_sides)) {
  square_side <- square_sides[i]
  # Increase the number of bins for r
  breaks_r <- round(steps_in * square_side)

  square_window <- make_square_win(square_side, c(-square_side/2, -square_side/2))
  square_area <- square_side^2

  # Maximal clustering
  pp_mc <- ppp(pp_loc$x, pp_loc$y, window = square_window)

  # Random
  pp_random <- rpoispp(intensity_rand, win=square_window)

  # Maximal dispersion
  column_distance <- sqrt(disp_distance^2 - (disp_distance / 2)^2)
  grid_pp <- max_disp_pp(square_window, disp_distance)

  # Compute rmax based on the H1 properties
  rmax <- sqrt(square_area/pi)

  # Compute summary functions
  # K
  k_mc <- Kest(pp_mc, correction = correction, r = seq(0, rmax, length.out=breaks_r))
  k_rand <- Kest(pp_random, correction = correction, r = seq(0, rmax, length.out=breaks_r))
  k_md <- Kest(grid_pp, correction = correction, r = seq(0, rmax, length.out=breaks_r))

  k_df_ss <- data.frame(r=k_mc$r, sf = c(k_mc$iso, k_rand$iso, k_md$iso),
                        pattern = c(rep("mc", nrow(k_mc)), rep("rand", nrow(k_rand)), rep("md", nrow(k_md))),
                        square=square_side^2)

  k_df <- rbind(k_df, k_df_ss)

  # PCF
  pcf_mc <- pcf(pp_mc, correction = correction, r = seq(0, rmax, length.out=breaks_r))
  pcf_rand <- pcf(pp_random, correction = correction, r = seq(0, rmax, length.out=breaks_r))
  pcf_md <- pcf(grid_pp, correction = correction, r = seq(0, rmax, length.out=breaks_r))

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
  lcf_mc <- LCFest(pp_mc, correction = correction, r = seq(0, rmax, length.out=breaks_r), rmax=rmax, dim=dim_mc)
  lcf_rand <- LCFest(pp_random, correction = correction, r = seq(0, rmax, length.out=breaks_r), rmax=rmax)
  dim_md <- lcf_md_dims[i]
  lcf_md <- LCFest(grid_pp, correction = correction, r = seq(0, rmax, length.out=breaks_r), rmax=rmax, dim=dim_md)

  lcf_df_ss <- data.frame(r=lcf_mc$r, sf = c(lcf_mc$iso, lcf_rand$iso, lcf_md$iso),
                          pattern = c(rep("mc", nrow(lcf_mc)), rep("rand", nrow(lcf_rand)), rep("md", nrow(lcf_md))),
                          square=square_side^2)

  lcf_df <- rbind(lcf_df, lcf_df_ss)

  # Visualize the point patterns
  if (i == 1) {
    pdf_side <- 1.05
    # Change window to be the target domain
    pdf_out(file.path(img_folder, "Max_clust.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(pp_mc, clustered_color, square_window, scale=0.1,
                   scale_offset=0.05, scale_lwd = 3, scale_right=TRUE,
                   text_height=0.1, text_width=0.32, cex=0.5)
    dev.off()

    # Random
    pdf_out(file.path(img_folder, "Rand.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(pp_random, random_color, square_window, cex=0.5)
    dev.off()

    # Maximum dispersion
    pdf_out(file.path(img_folder, "Max_disp.pdf"), width = pdf_side, height = pdf_side)
    save_pp_as_pdf(grid_pp, dispersed_color, square_window, cex=0.5)
    dev.off()
  }
}

# Save K-function
write_csv(k_df, file.path(data_folder, "k.csv"))

# Save normalized H
h1_df <- k_df %>%
  mutate(sf = sqrt(sf / pi) / sqrt(square / pi) - r / sqrt(square / pi))

write_csv(h1_df, file.path(data_folder, "h1.csv"))

# Save PCF
write_csv(pcf_df, file.path(data_folder, "pcf.csv"))

# Save LCF
write_csv(lcf_df, file.path(data_folder, "lcf.csv"))
