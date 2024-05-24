rm(list=ls())

library(dplyr)
library(spatstat)
library(tidyr)
library(purrr)
library(stats)
library(tictoc)
library(ggplot2)
library(lcfstat)
library(boot)

# Script for running the power analysis

source("../utils_lcf.R")

output_folder <- "Power strauss sample"
dir.create(output_folder, showWarnings = FALSE)

square_side <- 200
square_window <- owin(xrange=c(0, square_side), yrange=c(0, square_side))
square_a <- square_side^2

# Methods parameters
correction <- "Ripley"
dim_lims <- c(4, 50)

# Parameters to change
inh_distance <- 5
nums_points <-c(50, 100, 200)  #seq(40, 200, 20)
h1_pprocess <- "strauss"
gamma <- 0.25

# Num of patterns generated with H1
num_h1 <- 100#1000
num_h0 <- 50#5000
num_boot <- 2000
significance_level <- 0.05

# Stat to track
num_LCFest_failed <- 0

output_folder <- file.path(output_folder, paste0("H1_", num_h1, "_H0_", num_h0, "_B", num_boot, "_S", square_side, "_G", gamma, "_", h1_pprocess))
dir.create(output_folder, showWarnings = FALSE)

reject_stat <- function(trial_stats, i, stat_fun_data) {
  trial_stats <- trial_stats[i]
  num_trials <- length(trial_stats)
  sample_h0_stat <- stat_fun_data$h0
  sample_size <- length(sample_h0_stat)
  rejected_sl <- NULL

  pvals <- NULL
  for (trial_stat in trial_stats) {
    pval <- sum(sample_h0_stat <= trial_stat) / sample_size
    pvals <- c(pvals, pval)
  }

  rejected <- sum(pvals <= significance_level) / num_trials
  return(rejected)
}

pval_boot_bca <- function(tiral_h0_stat, stat_fun, stat_fun_data=NULL, conf.int=0.95, B=2000) {
  b <- boot(tiral_h0_stat, statistic = stat_fun, R = B, stat_fun_data=stat_fun_data)

  num_cases <- length(b$t0)
  boot_stats <- NULL
  for (ind in 1:num_cases) {
    result <- tryCatch(
      expr = {
        bootci <- boot.ci(b, conf = conf.int, "bca", index = ind)
        if (!is.null(bootci)) {
          mean <- bootci$t0
          low <- bootci$bca[4]
          high <- bootci$bca[5]
          boot_stat <- data.frame(mean=mean, low=low, high=high)
        } else {
          boot_stat <- data.frame(mean=b$t0[ind], low=b$t0[ind], high=b$t0[ind])
        }
      },
      error = function(e){
        print(e)
        boot_stat <- data.frame(mean=b$t0[ind], low=b$t0[ind], high=b$t0[ind])
        return(boot_stat)
      }
    )
    boot_stats <- rbind(boot_stats, result)
  }

  return(boot_stats)
}

set.seed(2711)

power_h <- NULL
power_pcf <- NULL
power_lcf <- NULL

for (num_points in nums_points) {
  print(paste("Number of points", num_points))

  tic()

  experiment_name <- paste0("IHD", inh_distance, "_PN", num_points)
  experiment_folder <- file.path(output_folder, experiment_name)
  dir.create(file.path(output_folder, experiment_name), showWarnings = FALSE)
  
  # Approximate null distribution of summary functions
  intensity_pois_gen <- num_points * 1.1 / square_a
  poisson_processes <- sim_pp_pnum(num_points, rpoispp, num_h0, intensity_pois_gen, win=square_window)
  
  lcf_out <- process_pps_lcf(poisson_processes, rpoispp, num_points, intensity_pois_gen,
                            inh_distance, dim_lims, correction,
                            win=square_window)
  
  poisson_processes <- lcf_out$updated_pps
  # Vector of length num_h1
  lcf_h0_stats <- lcf_out$lcf_stat
  lcf_h0_data <- list(h0=lcf_h0_stats)
  
  h_func_h0_stats <- summ_stat(poisson_processes, inh_distance, h_at_r, correction=correction)
  h_h0_data <- list(h0=h_func_h0_stats)
  
  pcf_h0_stats <- summ_stat(poisson_processes, inh_distance, pcf_at_r, correction=correction)
  pcf_h0_data <- list(h0=pcf_h0_stats)
  
  # Calculate summary statics for H1
  intensity_gen <- num_points * 1.5 / square_a
  hc_processes <- sim_pp_pnum(num_points, rStrauss, num_h1, intensity_gen, gamma = gamma, R=inh_distance, W=square_window)

  lcf_out <- process_pps_lcf(hc_processes, rStrauss, num_points, intensity_gen,
                             inh_distance, dim_lims, correction,
                             W=square_window, R=inh_distance, gamma = gamma)

  hc_processes <- lcf_out$updated_pps
  # Vector of length num_h1
  lcf_h1_stats <- lcf_out$lcf_stat

  # Calculate test statistics for K-function, PCF, F-Function
  # Vector of length num_h1
  h_func_h1_stats <- summ_stat(hc_processes, inh_distance, h_at_r, correction=correction)
  # Vector of length num_h1
  pcf_h1_stats <- summ_stat(hc_processes, inh_distance, pcf_at_r, correction=correction)

  # Get power estimate with bootstrap
  h_boot_rej <- pval_boot_bca(h_func_h1_stats, reject_stat, h_h0_data, B=num_boot)
  pcf_boot_rej <- pval_boot_bca(pcf_h1_stats, reject_stat, pcf_h0_data, B=num_boot)
  lcf_boot_rej <- pval_boot_bca(lcf_h1_stats, reject_stat, lcf_h0_data, B=num_boot)

  # Visualize the sampling distribution of test stats under H0
  h_xlab <- sprintf("H(%.2f)", inh_distance)
  h_sampl_distr_df <- data.frame(test_stat=h_func_h0_stats)
  write.csv(h_sampl_distr_df, file.path(experiment_folder, sprintf("%s_H0_sampl_distr.csv", h_xlab)))
  h_sampl_distr_p <- sampl_distr_plot(h_sampl_distr_df, h_xlab)
  ggsave(file.path(experiment_folder, sprintf("%s_H0_sampl_distr.pdf", h_xlab)), h_sampl_distr_p, width=5, height=4)

  pcf_xlab <- sprintf("PCF(%.2f)", inh_distance)
  pcf_sampl_distr_df <- data.frame(test_stat=pcf_h0_stats)
  write.csv(pcf_sampl_distr_df, file.path(experiment_folder, sprintf("%s_H0_sampl_distr.csv", pcf_xlab)))
  pcf_sampl_distr_p <- sampl_distr_plot(pcf_sampl_distr_df, pcf_xlab)
  ggsave(file.path(experiment_folder, sprintf("%s_H0_sampl_distr.pdf", pcf_xlab)), pcf_sampl_distr_p, width=5, height=4)

  lcf_xlab <- sprintf("LCF(%.2f)", inh_distance)
  lcf_sampl_distr_df <- data.frame(test_stat=lcf_h0_stats)
  write.csv(lcf_sampl_distr_df, file.path(experiment_folder, sprintf("%s_H0_sampl_distr.csv", lcf_xlab)))
  lcf_sampl_distr_p <- sampl_distr_plot(lcf_sampl_distr_df, lcf_xlab)
  ggsave(file.path(experiment_folder, sprintf("%s_H0_sampl_distr.pdf", lcf_xlab)), lcf_sampl_distr_p, width=5, height=4)

  time_elapsed <- toc()

  power_file_path <- file.path(experiment_folder, "power.txt")
  file_conn <- file(power_file_path)
  lines <- c(
    paste("H", h_boot_rej$mean, h_boot_rej$low, h_boot_rej$high, sep = "\t"),
    paste("PCF", pcf_boot_rej$mean, pcf_boot_rej$low, pcf_boot_rej$high, sep = "\t"),
    paste("LCF", lcf_boot_rej$mean, lcf_boot_rej$low, lcf_boot_rej$high, sep = "\t"),
    paste("time", as.numeric(time_elapsed$toc - time_elapsed$tic), sep = "\t"),
    paste("num_LCFest_failed", as.numeric(num_LCFest_failed), sep = "\t")
  )
  writeLines(lines, file_conn)
  close(file_conn)

  power_h <- rbind(power_h, h_boot_rej)
  power_pcf <- rbind(power_pcf, pcf_boot_rej)
  power_lcf <- rbind(power_lcf, lcf_boot_rej)
}

power_h$stat <- "h"
power_h$num_points <- nums_points

power_pcf$stat <- "pcf"
power_pcf$num_points <- nums_points

power_lcf$stat <- "lcf"
power_lcf$num_points <- nums_points

power_df <- rbind(power_h, power_pcf, power_lcf) %>%
  mutate(stat=factor(stat, levels=c("h", "pcf", "lcf")))

write.csv(power_df, file.path(output_folder, "power_sample.csv"))
