rm(list=ls())

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(egg)

source("../settings.R")

output_folder <- "data"

h_color <- "#D6984D"
pclcf_color <- "#008F0C"
lcf_color <- "#4A16DB"

expansion_x <- expansion(mult = 0.01, add = 0)
expansion_y <- expansion(mult = 0.02, add = 0)
line_width <- 1 / ggp2_magic_number

num_to_str <- function(num) {
  return(sprintf("N = %d", num) )
}

read_power_sample_file <- function(power_csv_path) {
  power_df <- read.csv(power_csv_path)
  
  power_df <- power_df %>% 
    mutate(stat = case_when(stat == "h" ~ "H(r)",
                            stat == "pcf" ~ "PCF(r)",
                            stat == "lcf" ~ "LCF(r)")) %>%
    mutate(stat = factor(stat, levels=c("H(r)", "PCF(r)", "LCF(r)")))
  
  return(power_df)
}

build_power_sample_plot <- function(power_df, facet_labeller, breaks_x=waiver(), uncertainty=FALSE, vj=0.5, show_y_title=FALSE) {
  
  breaks_y <- c(0, 0.5, 1)
  
  x_lim <- c(min(power_df$num_points), max(power_df$num_points))
  y_lim <- c(0, 1)
  
  power_df <- power_df %>%
    mutate(title=factor("title"))
  
  y_lab <- if (show_y_title) "Power" else ""
  
  power_p <- ggplot(power_df) +
    geom_line(aes(x=num_points, y=mean, col=stat), linewidth=line_width) +
    scale_colour_manual(values = c(h_color, pclcf_color, lcf_color)) +
    scale_fill_manual(values = c(h_color, pclcf_color, lcf_color)) +
    facet_wrap(~title, labeller = facet_labeller, strip.position = "top", nrow=1) +
    scale_x_continuous(breaks = breaks_x, expand = expansion_x, limits=x_lim) +
    scale_y_continuous(breaks = breaks_y, expand = expansion_y, limits=y_lim) +
    theme_bw(base_size = default_pointsize, base_family = default_font) + 
    ylab(y_lab) +
    xlab("N") +
    theme(plot.margin = margin(0, 7, 0, 2),
          strip.background = element_blank(),
          strip.clip = "off",
          strip.text.x = element_text(vjust = vj, size = default_pointsize, family = default_font, colour = "black"),
          legend.title=element_blank(),
          axis.line = element_line(colour = "black", linewidth=line_width),
          axis.text = element_text(size = default_pointsize - 1, family = default_font, colour = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position="none")
  
  if (uncertainty) {
    power_p <- power_p +
      geom_ribbon(aes(x=num_points, ymin = low, ymax = high, fill=stat), alpha = 0.2)
  }
  
  return(power_p)
}

read_power_effect_data <- function(file_path_template, point_nums) {
  
  power_df_glob <- NULL
  
  for (i in 1:length(point_nums)) {
    
    point_num <- point_nums[i]
    file_path <- sprintf(file_path_template, point_num)
    
    power_df <- read.csv(file_path)
    power_df[["point_num"]] <- point_num
    power_df[["ind"]] <- sprintf("N = %d", point_num) 
    power_df_glob <- rbind(power_df_glob, power_df)
  }
  
  ind_levels <- sapply(point_nums, num_to_str)
  
  power_df_glob <- power_df_glob %>% 
    mutate(stat = case_when(stat == "h" ~ "H(r)",
                            stat == "pcf" ~ "PCF(r)",
                            stat == "lcf" ~ "LCF(r)")) %>%
    mutate(stat = factor(stat, levels=c("H(r)", "PCF(r)", "LCF(r)")),
           ind = factor(ind, levels=ind_levels))
  
  return(power_df_glob)
}

build_power_effect_plot <- function(power_df, x_col, x_name, x_lim,
                                    breaks_x, breaks_y, 
                                    labels_x=waiver(),
                                    test_x_val=NULL,
                                    show_legend = FALSE,
                                    vj=0.5) {

  if (show_legend) {
    legend_text = element_text(size = default_pointsize, family = default_font, colour = "black")
    legend_position = "top"
  } else {
    legend_text = element_blank()
    legend_position = "none"
  }
   
  power_p <- ggplot(power_df) +
    geom_line(aes(x=!!x_col, y=mean, col=stat), linewidth=line_width) +
    geom_ribbon(aes(x=!!x_col, ymin = low, ymax = high, fill=stat), alpha = 0.2) +
    scale_colour_manual(values = c(h_color, pclcf_color, lcf_color)) +
    scale_fill_manual(values = c(h_color, pclcf_color, lcf_color)) +
    facet_wrap(~ind, scales = "free_x", strip.position = "top", nrow=1) +
    scale_x_continuous(breaks = breaks_x, labels=labels_x, expand = expansion_x, limits=x_lim) +
    scale_y_continuous(breaks = breaks_y, expand = expansion_y, limits=c(0, 1)) +
    theme_bw(base_size = default_pointsize, base_family = default_font) + 
    xlab(x_name) +
    theme(plot.margin = margin(0, 3, 0, 0),
          panel.spacing = unit(0.3, "lines"),
          strip.background = element_blank(),
          strip.text.x = element_text(vjust = vj, size = default_pointsize, family = default_font, colour = "black"),
          axis.text.x = element_text(size = default_pointsize - 1, family = default_font, colour = "black"),
          axis.line.x = element_line(colour = "black", linewidth=line_width),
          axis.line.y =  element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          legend.text = legend_text,
          legend.position = legend_position)
  
  if (!is.null(test_x_val)) {
    power_p <- power_p +
      geom_segment(aes(x = test_x_val, y = 0, xend = test_x_val, yend = 1), linewidth=line_width, color="gray", linetype="dashed") 
  }
  
  return(power_p)
}

# Hardcore
# Sample size
data_folder <- "data/hardcore"
power_csv_path <- file.path(data_folder, "power_sample.csv")

power_sample_df <- read_power_sample_file(power_csv_path)
hardcore_labeller <-label_bquote(cols = "R = 5")

power_sample_hardcore_p <- build_power_sample_plot(power_sample_df, hardcore_labeller, breaks_x = c(50, 60, 70), uncertainty = FALSE)
power_sample_hardcore_p

pdf_out(file.path(data_folder, "Power_hardcore_sample.pdf"), width=1, height=0.9)
print(power_sample_hardcore_p)
dev.off()

# Effect size
file_path_template <- file.path(data_folder, "power_pn%s.csv")
point_nums <- c(60, 120, 240)
  
power_effect_df <- read_power_effect_data(file_path_template, point_nums)
  
breaks_x <- c(0, 1, 3, 5)
breaks_y <- c(0, 0.5, 1)
  
column <- sym("ihd")
test_x_val <- 5
x_name <- "R"
x_lim <- c(min(power_effect_df$ihd), max(power_effect_df$ihd))
  
power_effect_hardcore_p <- build_power_effect_plot(power_effect_df, column,
                                                   x_name, x_lim, breaks_x, breaks_y)

power_effect_hardcore_p

pdf_out(file.path(data_folder, "Power_effect.pdf"), width=2.5, height=0.9)
print(power_effect_hardcore_p)
dev.off()


funcs_hardcore_plot <- ggarrange(power_sample_hardcore_p, 
                                 power_effect_hardcore_p, 
                                 nrow = 1,
                                 widths = c(1, 3))

funcs_hardcore_plot

pdf_out(file.path(data_folder, "Power_hardcore.pdf"), width=3.3, height=1.2)
print(funcs_hardcore_plot)
dev.off()


# Strauss
# Sample size
data_folder <- "data/strauss"
power_csv_path <- file.path(data_folder, "power_sample.csv")

power_sample_df <- read_power_sample_file(power_csv_path)
strauss_labeller <-label_bquote(cols = "R = 5," ~ gamma ~ "=" ~ 0.25)

power_sample_strauss_p <- build_power_sample_plot(power_sample_df, strauss_labeller, show_y_title=TRUE,
                                                  breaks_x = c(40, 120, 200), uncertainty = TRUE,
                                                  vj=0)
power_sample_strauss_p

pdf_out(file.path(data_folder, "Power_strauss_sample.pdf"), width=1, height=0.9)
print(power_sample_strauss_p)
dev.off()

# Effect size
file_path_template <- file.path(data_folder, "power_pn%s.csv")
point_nums <- c(60, 120, 240)

power_effect_df <- read_power_effect_data(file_path_template, point_nums)
power_effect_df <- power_effect_df %>%
  mutate(gamma_inv = 1 - gamma)

breaks_x <- c(0, 0.5, 1)
breaks_y <- c(0, 0.5, 1)

column <- sym("gamma_inv")
x_name <- expression(1 ~ "-" ~ gamma)
x_lim <- c(min(power_effect_df$gamma), max(power_effect_df$gamma))

power_effect_strauss_p <- build_power_effect_plot(power_effect_df, column,
                                                   x_name, x_lim, breaks_x, breaks_y,
                                                  c("0", "0.5", "1"), vj=1)

power_effect_strauss_p

pdf_out(file.path(data_folder, "Power_effect.pdf"), width=2.5, height=0.9)
print(power_effect_strauss_p)
dev.off()

funcs_strauss_plot <- ggarrange(power_sample_strauss_p, 
                                power_effect_strauss_p, 
                                nrow = 1,
                                widths = c(1, 3))

funcs_strauss_plot

pdf_out(file.path(data_folder, "Power_strauss.pdf"), width=3.3, height=1.2)
print(funcs_strauss_plot)
dev.off()

# Matern cluster
# Sample size
data_folder <- "data/matern_cluster"
power_csv_path <- file.path(data_folder, "power_sample.csv")

power_sample_df <- read_power_sample_file(power_csv_path)
matern_labeller <-label_bquote(cols = R[C] ~ "=" ~ 5 ~ "," ~  N[p] ~ "=" ~ 50)
 
power_sample_mat_clust_p <- build_power_sample_plot(power_sample_df, matern_labeller, 
                                                    uncertainty=TRUE,
                                                    vj= -0.5)
power_sample_mat_clust_p

pdf_out(file.path(data_folder, "Power_matern_sample.pdf"), width=1.2, height=1.2)
print(power_sample_mat_clust_p)
dev.off()

point_nums <- c(25, 50, 100)

# Effect size, radius
file_path_template <- file.path(data_folder, "power_pn%s.csv")

power_effect_df <- read_power_effect_data(file_path_template, point_nums)

breaks_x <- c(10, 50, 100)
breaks_y <- c(0, 0.5, 1)

column <- sym("cluster_rad")
x_name <- expression(R[C])
x_lim <- c(min(power_effect_df$cluster_rad), max(power_effect_df$cluster_rad))

test_x_val <- 10
power_effect_matern_radius_plot <- build_power_effect_plot(power_effect_df, column,
                                                         x_name, x_lim, breaks_x, breaks_y,
                                                         vj=1)

power_effect_matern_radius_plot <- power_effect_matern_radius_plot + 
  scale_x_reverse()

power_effect_matern_radius_plot

pdf_out(file.path(data_folder, "Power_effect_radius.pdf"), width=2.5, height=1.2)
print(power_effect_matern_radius_plot)
dev.off()

funcs_matern_plot <- ggarrange(power_sample_mat_clust_p, 
                               power_effect_matern_radius_plot, 
                               nrow = 1,
                               widths = c(1, 3))

funcs_matern_plot

pdf_out(file.path(data_folder, "Power_matern.pdf"), width=3.3, height=1.2)
print(funcs_matern_plot)
dev.off()



