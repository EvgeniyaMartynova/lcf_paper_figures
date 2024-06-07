rm(list=ls())

library(sf)
library(spatstat)
library(USA.state.boundaries)
library(dplyr)
library(aistat)
library(dbscan)
library(scales)
library(RColorBrewer)
library(plotrix)
library(ggplot2)

source("../utils_lcf.R")

data_folder <- "data_to_plot/cell"
output_folder <- "img/tmas"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Palettes
# Tumor color
blue_palette <- brewer.pal(9, "Blues")
tumor_color <- rgb(200 / 255, 209 / 255, 230 / 255, alpha=0.75)

# B cell palette
purple_cols <- c("#800066", "#b3008f",  "#e600b8", "#ff1ad1")
col_func <- colorRampPalette(purple_cols)
b_cell_palette <- col_func(7)

# cDC2 palette
red_palette <- brewer.pal(9, "Reds")
col_func <- colorRampPalette(red_palette[8:4])
cdc2_palette <- col_func(7)

# pDC
pdc_palette <- c("#1b6031", "#20733a", "#258644", "#2a994e",  "#30ac57", "#49c571", "#5dcc81")

# Myeloid palette
yellow_cols <- c("#807400", "#b3a300", "#e6d200", "#ffeb1a")
col_func <- colorRampPalette(yellow_cols)
myeloid_palette <- col_func(7)

palettes <- list("B"=b_cell_palette,
                 "cDC2"=cdc2_palette,
                 "pDC"=pdc_palette,
                 "M"=myeloid_palette)

rotate <- function(a) matrix(c(-cos(a), sin(a), sin(a), cos(a)), 2, 2)

clust_color_intens <- function(points_num, intesity, pallete) {
  if (points_num == 1) {
    return(pallete[1])
  }

  color <-  case_when(intesity < 0.005 ~ pallete[2],
                      intesity > 0.005 && intesity < 0.0065 ~ pallete[3],
                      intesity > 0.0065 && intesity < 0.008 ~ pallete[4],
                      intesity > 0.008 && intesity < 0.010 ~ pallete[5],
                      intesity > 0.010 && intesity < 0.012 ~ pallete[6],
                      intesity > 0.012 ~ pallete[7])

  return(color)
}


cluster_cells <- function(coords, dist_th=8) {
  cell_dist <- dist(ct_coords)
  clust <- dbscan(cell_dist, dist_th, minPts = 2)
  return(clust$cluster)
}

plot_clusters <- function(ct_coords, cluster_ind, palette, rad) {

  single_points <- ct_coords[cluster_ind == 0,]
  for (i in 1:nrow(single_points)) {
    draw.circle(single_points$X[i], single_points$Y[i], rad,
                col=palette[1], border = palette[1])
  }

  sf_cells <- apply(as.matrix(ct_coords), 1, st_point, simplify = FALSE)
  sfc_cells <- st_sfc(sf_cells)

  num_clusters <- max(cluster_ind)
  clust_dens <- NULL
  clust_lens <- NULL
  clust_buffs <- vector(mode = "list", length = num_clusters)
  clust_colors <- NULL

  if (num_clusters > 0) {
    for (i in 1:num_clusters) {
      cluster <- ct_coords[cluster_ind == i,]

      cluster_mp <- st_multipoint(as.matrix(cluster))
      cluster_buffer <- st_buffer(cluster_mp, dist = rad)
      clust_buffs[[i]] <- cluster_buffer

      cluster_buffer_dens <- st_buffer(cluster_mp, dist = rad * 1.2)
      cluster_ch <- st_convex_hull(cluster_buffer_dens)
      inCluster <- st_intersects(sfc_cells, cluster_ch, sparse = FALSE)
      cells_in_ch <- ct_coords[inCluster, ]

      dens <- nrow(cluster) / spatstat.geom::area(cluster_ch)
      clust_dens <- c(clust_dens, dens)

      clust_lens <- c(clust_lens, nrow(cluster))

      color <- clust_color_intens(nrow(cluster_mp), dens, palette)
      color <- alpha(color, 1)
      clust_colors <- c(clust_colors, color)
    }

    clust_num_df <- st_sf(ind=1:num_clusters,
                          num=clust_lens,
                          buffer=clust_buffs,
                          color=clust_colors) %>%
      arrange(desc(num))

    for (i in 1:num_clusters) {
      color <- clust_num_df$color[i]
      plot(clust_num_df$buffer[i],
           add=TRUE,
           col=color,
           border=color)
    }
  }
}

dataset <- "2022-05-19_Iris_DCpanelTMA9"
slide <- "TMA9_2A_Scan1"
patient <- "T99-11009_IF"
tma_num <- "1"

# Load tumor and stroma segmentation
tumor_path <- str_c(data_folder, dataset, slide, "tumor.wkt", sep = "/")
tumor_pols <- sf_polygons_from_wkt(tumor_path)

stroma_path <- str_c(data_folder, dataset, slide, "stroma.wkt", sep = "/")
stroma_pols <- sf_polygons_from_wkt(stroma_path)

# Load polygon from disk
tma_roi <- load_tma_pol(dataset, slide, patient, tma_num, data_folder)

roi_coords <- get_roi_coords(tma_roi, origin_tranform=FALSE)
roi_coords[nrow(roi_coords), ] <- roi_coords[1,]

roi_tumor <- pols_in_roi(roi_coords, tumor_pols)
roi_stroma <- pols_in_roi(roi_coords, stroma_pols)

cells <- load_cells(dataset, slide, patient, tma_num, NULL, NULL,  prediction_path=data_folder)

cells <- cells %>%
  filter(!(type %in% c("cDC1", "Double pos cDC2"))) %>%
  mutate(t=case_when(type == "B cell" ~ "B",
                     type == "cDC2" ~ "cDC2",
                     type == "pDC" ~ "pDC",
                     type == "Myeloid cell" ~ "M"))


roi_mtpol <- st_multipolygon(c(roi_tumor, roi_stroma))
roi_contours <- get_roi_contours(roi_mtpol)

for (ctype in c("B", "cDC2", "pDC", "M")) {

  ct_coords <- cells %>%
    filter(t == ctype) %>%
    dplyr::select(X, Y)

  palette <- palettes[[ctype]]

  cells_pdf_path <- file.path(output_folder, sprintf("%s_vis.pdf", ctype))

  pdf(cells_pdf_path, width = 1, height = 1)
  par(mfrow=c(1,1), mar = c(0,0,0,0), lwd=1)

  plot(roi_tumor, col=tumor_color, border=tumor_color)
  #plot(roi_contours, lwd=0.7, add=TRUE)

  cluster_ind <- cluster_cells(ct_coords)
  plot_clusters(ct_coords, cluster_ind, palette, 6)
  dev.off()
}
