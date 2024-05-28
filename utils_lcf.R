## Utils for the code for paper figures only
library(jsonlite)
library(stringr)
library(sf)

make_rect_df <- function(width, height, origin=c(0, 0)) {
  xo <- origin[1]
  yo <- origin[2]

  square_df <- data.frame(X=c(xo + width / 2, xo - width / 2, xo - width / 2, xo + width / 2),
                          Y=c(yo + height / 2, yo + height / 2, yo - height / 2, yo - height / 2))
  return(square_df)
}

make_square_df <- function(side, origin=c(0, 0)) {
  square_df <- make_rect_df(side, side, origin)
  return(square_df)
}

make_circle_df <- function(xo, yo, radius) {
  angles <- seq(0, 2*pi, length.out=100)
  xs <- cos(angles) * radius + xo
  ys <- sin(angles) * radius + yo
  circle_df <- data.frame(x=xs, y=ys)
  return(circle_df)
}

# Generates point pattern maximally dispersed at a given distance
max_disp_pp <- function(window, dist, show_plot=FALSE) {
  
  window_df <- as.data.frame(cbind(window$bdry[[1]]$x, window$bdry[[1]]$y))
  window_df[nrow(window_df) + 1,] <- window_df[1,]
  
  coords_list <- list(as.matrix(window_df))
  polygon <- st_polygon(coords_list)
  
  side <- dist * sin(pi/3) *2
  grid <- st_make_grid(polygon, c(side, side), square=FALSE)
  
  if (show_plot) {
    plot(grid)
    plot(polygon, add = TRUE)
  }
  
  tess_verts <- NULL
  for (hex in grid) {
    vert <- data.frame(hex[[1]][1:6, ])
    middle_y <- (max(vert$y) + min(vert$y)) / 2
    middle_x <- (max(vert$x) + min(vert$x)) / 2
    vert[nrow(vert) +1, ] <- c(middle_y, middle_x)
    
    tess_verts <- rbind(tess_verts, vert)
  }
  tess_verts_unique <- distinct(tess_verts)
  
  disp_pp <- ppp(tess_verts_unique$y, tess_verts_unique$x, window = window)
  return(disp_pp)
}

# Creates an owin (spatstat) object based on the roi coordinates
# given in anticlockwise order
window_from_roi <- function(roi_coords_anticlock) {

  roi_coords_anticlock <- as.data.frame(roi_coords_anticlock)
  roi_coords_anticlock[nrow(roi_coords_anticlock) + 1,] <- roi_coords_anticlock[1,]
  colnames(roi_coords_anticlock) <- c("x", "y")

  max_x_tr <- max(roi_coords_anticlock[,1])
  max_y_tr <- max(roi_coords_anticlock[,2])
  min_x_tr <- min(roi_coords_anticlock[,1])
  min_y_tr <- min(roi_coords_anticlock[,2])

  tma_window <- owin(xrange = c(min_x_tr, max_x_tr), yrange = c(min_y_tr, max_y_tr), poly=roi_coords_anticlock)
  return(tma_window)
}

closest_point <- function(points, xt, yt) {

  nn <- nn2(points, query = data.frame(X=xt, Y=yt), radius = 100, searchtype = c("radius"), k = 1)
  point_vis <- points[nn$nn.idx, ]
  return(point_vis)
}

# Data loading
# Loads polygon coordinates
load_tma_pol <- function(dataset, slide, patient, tma_num, prediction_path) {
  roi_path <- str_c(prediction_path, dataset, slide, patient, tma_num, "roi.json", sep = "/")
  roi <- as.data.frame(fromJSON(roi_path))
  
  colnames(roi) <- c("X", "Y")
  
  roi
}

# Extracts ROI coordinates and min X and max Y 
# from the object parsed from a JSON request of ROI polygon
# by default transforms the coordinates to cartesian coordinates
# with (minX, minY) = (0, 0)
get_roi_coords <- function(roi_pol, origin_tranform=TRUE) {
  roi_xs <- roi_pol$X
  roi_ys <- roi_pol$Y
  
  roi_min_x <- min(roi_xs)
  roi_max_y <- max(roi_ys)
  
  if (origin_tranform) {
    roi_xs <- roi_xs - roi_min_x 
    roi_ys <- roi_max_y - roi_ys
  } else {
    roi_ys <- -roi_ys
  }
  
  roi_coords <- data.frame(X=roi_xs, Y=roi_ys)
  roi_coords
}

# Creates an owin (spatstat) object based on the roi coordinates
# given in anticlockwise order
window_from_roi <- function(roi_coords_anticlock) {
  
  roi_coords_anticlock <- as.data.frame(roi_coords_anticlock)
  roi_coords_anticlock[nrow(roi_coords_anticlock) + 1,] <- roi_coords_anticlock[1,]
  colnames(roi_coords_anticlock) <- c("x", "y")
  
  max_x_tr <- max(roi_coords_anticlock[,1])
  max_y_tr <- max(roi_coords_anticlock[,2])
  min_x_tr <- min(roi_coords_anticlock[,1])
  min_y_tr <- min(roi_coords_anticlock[,2])
  
  tma_window <- owin(xrange = c(min_x_tr, max_x_tr), yrange = c(min_y_tr, max_y_tr), poly=roi_coords_anticlock)
  return(tma_window)
}

# Loads sf polygons from WKT files
sf_polygons_from_wkt <- function(wkt_path, rotate=TRUE) {
  wkt_file <- file(wkt_path, open="r")
  polygon_strs <-readLines(wkt_file)
  close(wkt_file)
  
  polygons <- NULL
  for (polygon_str in polygon_strs) {
    sf_pol <- st_as_sfc(polygon_str)
    if (rotate) {
      sf_pol <- sf_pol * rotate(pi)
    }
    polygons <- c(polygons, sf_pol)
  }
  
  return(polygons)
}

# Check which of polygons from the given list are inside the ROI
# roi_coords - data frame
# polygons - a vector of sf polygons
pols_in_roi <- function(roi_coords, polygons) {
  # Convert roi_coords to an sf object
  roi_st_pol <- st_polygon(list(as.matrix(roi_coords)))
  roi_sf_pol <- st_sf(st_sfc(roi_st_pol))
  
  # Convert polygons to an sf object
  sf_polys <- st_sfc(polygons)
  # Find indices of polygons that intersect ROI
  res <- st_intersects(roi_sf_pol, sf_polys, sparse=TRUE)
  inds <- res[[1]]
  
  pols_inside <- st_multipolygon(polygons[inds])
  return(pols_inside)
}

# Computes external contours of a multypolygon
get_roi_contours <- function(roi_polygon, dist=5) {
  buff <- st_buffer(roi_polygon, dist)
  
  if (is(buff,"MULTIPOLYGON")) {
    extern_pols <- vector(mode = "list", length = length(buff))
    for (i in 1:length(buff)) {
      extern_pol <- st_polygon(list(buff[[i]][[1]]))
      extern_pols[[i]] <- extern_pol
    }
  } else {
    extern_pols <- vector(mode = "list", length = 1)
    extern_pols[[1]] <- st_polygon(list(buff[[1]]))
  }
  
  contours <- st_multipolygon(extern_pols)
  return(contours)
}

# Determines whether polygon coordinates are in anticlockwise order
anticlockwise <- function(coords) {
  
  x.coords <- c(coords[,1], coords[1,1])
  y.coords <- c(coords[,2], coords[1,2])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area < 0
} 

# Loads predicted cells from disk, filters out the invalid cell types, makes
# if roi_max_y and roi_min_x are given, transforms coordinates to cartesian with (minX, minY) = (0, 0)
load_cells <- function(dataset, slide, patient, tma_num, 
                       roi_max_y=NULL, roi_min_x=NULL,
                       col_names = c("Y", "X", "BDCA1", "BDCA2", "XCR", "CD14", "CD19", "type"),
                       types_to_exclude = c("Other cell", "Invalid"),
                       prediction_path = "prediction") {
  
  cells_path <- str_c(prediction_path, dataset, slide, patient, tma_num, "prediction.txt", sep = "/")
  cells <- read.csv(cells_path, sep="\t", header = FALSE, col.names = col_names)
  
  cells <- cells %>%
    filter(!(type %in% types_to_exclude)) %>% 
    mutate(type = as.factor(type))
  
  if (!is.null(roi_max_y)) {
    cells <- cells %>%
      mutate(Y = roi_max_y - Y)
  } else {
    cells <- cells %>%
      mutate(Y = -Y)
  }
  
  if (!is.null(roi_min_x)) {
    cells <- cells %>%
      mutate(X = X - roi_min_x)
  }
  
  cells
}

# Simulate strauss process for efficiency
sim_strauss_mh <- function(num_points, gamma, inh_distance, window, nsim, nrep=1e6) {

  window_area <- spatstat.geom::area(window)
  intensity <- num_points / window_area

  strauss_model <- list(cif="strauss",
                        par=list(beta=intensity, gamma=gamma, r=inh_distance),
                        w=window)

  hc_pps <- rmh(model=strauss_model,
                start=list(n.start=800),
                control=list(nrep=nrep, p=1),
                nsim=nsim)

  return(hc_pps)
}

h_func <- function(pp, ...) {
  l_fun <- Lest(pp, ...)
  emp_name <- attr(l_fun, "names")[3]

  h_fun_df <- as.data.frame(l_fun) %>%
    mutate(theo = theo - r)

  h_fun_df[[emp_name]] <- l_fun[[emp_name]] - l_fun$r

  return(h_fun_df)
}

pc_func <- function(pp, ...) {
  pc_fun <- pcf(pp, ...)
  emp_name <- attr(pc_fun, "names")[3]

  pc_fun_df <- as.data.frame(pc_fun) %>%
    mutate(val_exp = theo, val = pc_fun[[emp_name]])

  return(pc_fun_df[2:nrow(pc_fun_df), ])
}

save_pp_as_pdf <- function(pp, col, disp_window_df, clust_rad=NULL, clust_o_x=0, clust_o_y=0, cex=1, 
                           scale=NULL, scale_offset=NULL, text_height=30) {
  points <- data.frame(X=pp$x, Y=pp$y)
  x_lim <- c(min(disp_window_df$X), max(disp_window_df$X))
  y_lim <- c(min(disp_window_df$Y), max(disp_window_df$Y))

  par(mfrow = c(1,1), mar = c(0,0,0,0), lwd=1)
  plot(1, 1, col = "white", asp=1, axes=FALSE, xlab="", ylab="",
       xlim=x_lim,  ylim=y_lim)
  polygon(disp_window_df$X, disp_window_df$Y)

  if (!is.null(clust_rad)) {
    draw_rad <- clust_rad + 5
    draw.circle(clust_o_x, clust_o_y, draw_rad, col="gray", border="black")

    # Plot domain size
    angle1 <- 3 * pi / 4
    angle2 <- - pi / 4

    domain_seg_x1 <- clust_o_x + cos(angle1) * draw_rad
    domain_seg_y1 <- clust_o_y + sin(angle1) * draw_rad
    domain_seg_x2 <- clust_o_x + cos(angle2) * draw_rad
    domain_seg_y2 <- clust_o_y + sin(angle2) * draw_rad
  }
  if (!is.null(scale) && !is.null(scale_offset)) {
    seg_x0 <- x_lim[1] + scale_offset
    seg_x1 <- seg_x0 + scale
    # For Y axis offset should be negative
    seg_y0 <- y_lim[1] + scale_offset * 1.3 + text_height
    seg_y1 <- seg_y0

    segments(seg_x0, seg_y0, x1 = seg_x1, y1 = seg_y1,
             col = "black", lwd=4, lend=2)

    label_x <- seg_x0
    label_y <- seg_y0 - text_height
    text(label_x, label_y, as.expression(bquote(bold(.(as.character(scale)))~bold(units))), cex=0.675, adj = 0)
  }
  points(points$X, points$Y, pch=16, col=col, cex=cex)
}


h_at_r <- function(pp, r, ...) {
  l_fun <- Lest(pp, r=r, ...)
  emp_name <- attr(l_fun, "names")[3]
  
  l_fun <- as.data.frame(l_fun)
  h_fun_df <- l_fun %>%
    mutate(theo = theo - r) %>%
    dplyr::select(r, theo)
  
  h_fun_df[[emp_name]] <- l_fun[[emp_name]] - l_fun[["r"]]
  
  return(h_fun_df)
}

pcf_at_r <- function(pp, r, ...) {
  pc_fun <- pcf(pp, ...)
  emp_name <- attr(pc_fun, "names")[3]
  
  pcf_at_r <- approx(pc_fun$r, pc_fun[[emp_name]], r)
  pc_fun_df <- data.frame(r=r, theo=rep(1, length(r)))
  pc_fun_df[[emp_name]] <- pcf_at_r$y
  
  return(pc_fun_df)
}

# Power - make necessary changes
sim_pp_pnum <- function(num_points, pp_fun, num_sim, intensity_gen, ...) {
  
  p_processes <- pp_fun(intensity_gen, nsim=round(num_sim*1.1), ...)
  
  if (num_sim == 1) {
    p_processes <- list(p_processes)
  }
  
  processes_trimmed <- vector(mode="list", length = 0)
  i <- 1
  j <- 1
  for (p_processe in p_processes) {
    if (p_processe$n >= num_points && i <= num_sim) {
      sample_inds <- sample(p_processe$n, num_points)
      trimmed_pp <- p_processe[sample_inds]
      processes_trimmed[[i]] <- trimmed_pp
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  
  if (length(processes_trimmed) < num_sim) {
    sim_left <- num_sim - length(processes_trimmed)
    for (m in 1:sim_left) {
      gen_num_points <- 0
      while (gen_num_points < num_points) {
        p_process <- pp_fun(intensity_gen, ...)
        gen_num_points <- p_process$n
      }
      
      sample_inds <- sample(p_process$n, num_points)
      trimmed_pp <- p_process[sample_inds]
      processes_trimmed[[i]] <- trimmed_pp
      i <- i + 1
    }
  }
  
  if (length(processes_trimmed) > num_sim) {
    processes_trimmed <- processes_trimmed[1:num_sim]
  }
  
  stopifnot(length(processes_trimmed)==num_sim)
  return(processes_trimmed)
}

summ_stat <- function(p_processes, r_check, stat_func, ...) {
  # Calculate test statistics for K-function, PCF, F-Function
  stats <- NULL
  
  for (p_process in p_processes) {
    stat_val <- stat_func(p_process, r=c(0, r_check), ...)
    stat <- stat_val[[3]][2]
    stats <- c(stats, stat)
  }
  
  return(stats)
}


sim_hardcore_mh <- function(num_points, inh_distance, window, nsim, nrep=1e6) {
  # Hard core process:
  hardcore_model <- list(cif="hardcore",
                         par=list(beta=num_points/square_a, hc=inh_distance),
                         w=window)
  
  hc_pps <- rmh(model=hardcore_model,
                start=list(n.start=num_points),
                control=list(p=1, nrep=nrep),
                nsim=nsim)
  
  return(hc_pps)
}


sampl_distr_plot <- function(null_test_stat_df, xlab) {
  null_test_stat_p <- ggplot(null_test_stat_df, aes(test_stat)) +
    geom_histogram() +
    ylab("Count") +
    xlab(xlab) +
    theme_bw()
  
  return(null_test_stat_p)
}


process_pps_lcf <- function(trial_pps,
                            pp_fun,
                            num_points,
                            intensity,
                            r_check,
                            dim_lims,
                            correction,
                            ...) {
  
  trial_stat_lcf <- NULL
  for (i in seq_along(trial_pps)) {
    lcf <- NULL
    
    while (is.null(lcf)) {
      p_process <- trial_pps[[i]]
      
      output <- tryCatch(
        error = function(cnd) {
          print(paste0("LCF computation failed, reason ", conditionMessage(cnd)))
          num_LCFest_failed <<- num_LCFest_failed + 1
          new_pp <- sim_pp_pnum(num_points, pp_fun, 1, intensity, ...)[[1]]
          return(new_pp)
        }, {
          lcf <- LCFest(p_process, r=c(0, r_check), dim_lims = dim_lims, correction=correction)
        }
      )
      
      if (inherits(output, "ppp")) {
        trial_pps[[i]] <- output
      }
    }
    
    stat <- lcf[[3]][2]
    trial_stat_lcf <- c(trial_stat_lcf, stat)
  }
  
  list(lcf_stat=trial_stat_lcf, updated_pps=trial_pps)
}

