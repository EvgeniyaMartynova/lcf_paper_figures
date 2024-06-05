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

data_folder <- "data/facilities"

output_folder <- "img/facilities"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

orange_palette <- brewer.pal(9, "Oranges")
orange_palette <- orange_palette[8:4]
col_func <- colorRampPalette(orange_palette)
air_palette <- col_func(7) 

blue_palette <- brewer.pal(9, "Blues")
blue_palette <- blue_palette[8:4]
col_func <- colorRampPalette(blue_palette)
waste_palette <- col_func(7) 

target_crs <- "ESRI:102010"

extract_sf_points <- function(df, x_col, y_col, map_df) {
  
  x <- df[[x_col]]
  y <- df[[y_col]]
  coords_df <- data.frame(x=x, y=y)
  coords_df <- coords_df[!duplicated(coords_df), ]
  
  points <- apply(as.matrix(coords_df), 1, sf::st_point, simplify = FALSE)
  points <- sf::st_sfc(points)
  st_crs(points) <- st_crs(map_df)
  points
}

extract_points_in_state <- function(points, state) {
  inROI <- sf::st_intersects(points, state, sparse = FALSE)
  state_points <- st_transform(points[inROI], crs = target_crs)
  state_points
}

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


cluster_points <- function(coords, dist_th=5000) {
  points_dist <- dist(coords)
  clust <- dbscan(points_dist, dist_th, minPts = 2)
  return(clust$cluster)
}

plot_clusters <- function(ct_coords, cluster_ind, palette, rad=8) {
  
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

# State boundary
usa <- state_boundaries_wgs84
usa <- usa[usa$TYPE == "Land", ]

state_name <- "Ohio"
state_bdry <- usa$Shape[usa$NAME == state_name]

# Facilities
# Air polluters
air_polluters <- read.csv(file.path(data_folder, "Air_Majors.csv"))

# Waste generators
load(file.path(data_folder, "tsdfs.RData")) 
waste_generators <- x
rm(x)

inds <- waste_generators$lqg==1 & waste_generators$airmajor==0
waste_generators <- waste_generators[inds, ]

# Transform facilities coordinates to SF points
air_pol_points <- extract_sf_points(air_polluters, "LONGITUDE", "LATITUDE", usa)
waste_gen_points <- extract_sf_points(waste_generators, "longitude", "latitude", usa)

# Get air pollution and waste generation facilities in the state
air_polluters_state <- extract_points_in_state(air_pol_points, state_bdry)
waste_generators_state <- extract_points_in_state(waste_gen_points, state_bdry)

# Project state border coordinates to Cartesian
state_bdry <- st_transform(state_bdry, crs = target_crs)

air_coords <- do.call(rbind, st_geometry(air_polluters_state)) %>%
  as_tibble() %>% setNames(c("X","Y"))

waste_coords <- do.call(rbind, st_geometry(waste_generators_state)) %>%
  as_tibble() %>% setNames(c("X","Y"))

state_bb <- st_bbox(state_bdry)
wh_ratio <- as.numeric((state_bb$xmax - state_bb$xmin) / 
  (state_bb$ymax - state_bb$ymin))

rad <- sqrt(area(state_bdry)) / 100

# Make air polluters visualization
width <- 0.827
height <- width / wh_ratio

air_vis_pdf_path <- file.path(output_folder, "air_vis_ohio.pdf")
pdf(air_vis_pdf_path, width = width, height = height)
par(mfrow=c(1,1), mar = c(0,0,0,0))

plot(state_bdry)
cluster_ind <- cluster_points(air_coords, dist_th = rad)
plot_clusters(air_coords, cluster_ind, air_palette, rad)

dev.off()

# Make waste generators visualization
width <- 0.827
height <- width / wh_ratio

waste_vis_pdf_path <- file.path(output_folder, "waste_vis_ohio.pdf")
pdf(waste_vis_pdf_path, width = width, height = height)
par(mfrow=c(1,1), mar = c(0,0,0,0))

plot(state_bdry)
cluster_ind <- cluster_points(waste_coords, dist_th = rad)
plot_clusters(waste_coords, cluster_ind, waste_palette, rad)

dev.off()

