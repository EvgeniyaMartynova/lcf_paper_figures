rm(list=ls())

library(dplyr)
library(spatstat)
library(tidyr)
library(purrr)
library(ggplot2)
library(lcfstat)
library(tictoc)
library(boot)

source("../utils_lcf.R")

calibration_folder <- "Calibration"
dir.create(calibration_folder, showWarnings = FALSE)

h_color <- "#D6984D"
pcf_color <- "#008F0C"
lcf_color <- "#4A16DB"

# Keep intensity the same based on quantities for

# Window square
square_sides <- c(200 / sqrt(2), 200, 200 * sqrt(2))

# Methods params
correction <- "Ripley"
dim_lims <- c(4, 50)

# Params to change
inh_distance <- 5
nums_points <- c(25, 50, 100)

# Monte Carlo test parameters
num_h0_sampling <- 5000
num_h0_trials <- 1000
num_boot <- 2000

# Stat to track
num_LCFest_failed <- 0


sampl_distr_plot <- function(null_test_stat_df, xlab) {
  null_test_stat_p <- ggplot(null_test_stat_df, aes(test_stat)) +
    geom_histogram() +
    ylab("Count") +
    xlab(xlab) +
    theme_bw()

  return(null_test_stat_p)
}

build_calibration_tb <- function(significance_levels,
                                 stat_rej_by_sign_h,
                                 stat_rej_by_sign_pcf,
                                 stat_rej_by_sign_f) {

  stat_rej_by_sign_h <- stat_rej_by_sign_h %>%
    mutate(stat="h", sl=significance_levels)

  stat_rej_by_sign_pcf <- stat_rej_by_sign_pcf %>%
    mutate(stat="pcf", sl=significance_levels)

  stat_rej_by_sign_f <- stat_rej_by_sign_f %>%
    mutate(stat="lcf", sl=significance_levels)

  calibration_df <- rbind(stat_rej_by_sign_h, stat_rej_by_sign_pcf, stat_rej_by_sign_f) %>%
    mutate(stat=factor(stat, levels=c("h", "pcf", "lcf")))

  return(calibration_df)
}

calibration_plot <- function(calibration_df) {

  breaks <- c(10^(-3), 10^(-2), 10^(-1), 10^0)

  calibration_p <- ggplot(calibration_df) +
    geom_line(aes(x=sl, y=mean, col=stat), linewidth=1.1) +
    geom_ribbon(aes(x=sl, ymin = low, ymax = high, fill=stat), alpha = 0.2) +
    scale_colour_manual(values = c(h_color, pcf_color, lcf_color)) +
    scale_fill_manual(values = c(h_color, pcf_color, lcf_color)) +
    geom_abline(intercept=0, slope=1, color="black", linetype = "dashed") +
    ylab("Type I error") +
    xlab("Significance level") +
    scale_x_log10(breaks = breaks, labels = breaks) +
    scale_y_log10(breaks = breaks, labels = breaks) +
    coord_fixed() +
    theme_bw()

  return(calibration_p)
}

qqplot_pvals <- function(pvals_stat_df, title) {
  pvals_theo <- seq(0, 1, length.out=nrow(pvals_stat_df))

  pvals_sorted <- pvals_stat_df %>%
    arrange(mean) %>%
    mutate(theo=pvals_theo) %>%
    mutate(low=if_else(low == 0, 1 / length(pvals_theo), low),
           mean=if_else(mean == 0, 1 / length(pvals_theo), mean),
           high=if_else(high == 0, 1 / length(pvals_theo), high))

  breaks <- c(10^(-3), 10^(-2), 10^(-1), 10^0)

  pvals_qq_p <- ggplot(pvals_sorted) +
    geom_line(aes(x=mean, y=theo), linewidth=1.1, col="red") +
    geom_ribbon(aes(y=theo, xmin = low, xmax = high, fill="red"), alpha = 0.5) +
    geom_abline(intercept=0, slope=1, color="black", linetype = "dashed") +
    ylab("P-values theorethical") +
    xlab("P-values empirical") +
    ggtitle(title) +
    coord_fixed() +
    theme_bw() +
    scale_x_log10(breaks = breaks, labels = breaks) +
    scale_y_log10(breaks = breaks, labels = breaks)

  return(pvals_qq_p)
}

build_pvals_tb <- function(stat_pvals_h, stat_pvals_pcf, stat_pvals_lcf) {
  pvals_theo <- seq(0, 1, length.out=length(stat_pvals_h))

  pvals_sorted_h <- data.frame(observed=sort(stat_pvals_h)) %>%
    mutate(stat="h", theo=pvals_theo)

  pvals_sorted_pcf <- data.frame(observed=sort(stat_pvals_pcf)) %>%
    mutate(stat="pcf", theo=pvals_theo)

  pvals_sorted_f <- data.frame(observed=sort(stat_pvals_lcf)) %>%
    mutate(stat="lcf", theo=pvals_theo)

  pvals_qq_df <- rbind(pvals_sorted_h, pvals_sorted_pcf, pvals_sorted_f) %>%
    mutate(stat=factor(stat, levels=c("h", "pcf", "lcf")))

  return(pvals_qq_df)
}

qqplot_pvals_together <- function(pvals_qq_df) {

  breaks <- c(10^(-3), 10^(-2), 10^(-1), 10^0)

  pvals_qq_p <- ggplot(pvals_qq_df) +
    geom_line(aes(x=observed, y=theo, col=stat), linewidth=1.1) +
    scale_colour_manual(values = c(h_color, pcf_color, lcf_color)) +
    scale_fill_manual(values = c(h_color, pcf_color, lcf_color)) +
    geom_abline(intercept=0, slope=1, color="black", linetype = "dashed") +
    ylab("P-values theorethical") +
    xlab("P-values empirical") +
    coord_fixed() +
    theme_bw() +
    scale_x_log10(breaks = breaks, labels = breaks) +
    scale_y_log10(breaks = breaks, labels = breaks)

  return(pvals_qq_p)
}

get_pvals <- function(trial_stats, sample_h0_stat, expected_h0_stat=0) {
  num_trials <- length(trial_stats)
  sample_size <- length(sample_h0_stat)
  rejected_sl <- NULL

  pvals <- NULL
  for (trial_stat in trial_stats) {
    if (trial_stat < expected_h0_stat) {
      symm_stat <- expected_h0_stat + abs(expected_h0_stat - trial_stat)
      pval <- (sum(sample_h0_stat <= trial_stat) + sum(sample_h0_stat >= symm_stat)) / sample_size
    } else {
      symm_stat <- expected_h0_stat - abs(expected_h0_stat - trial_stat)
      pval <- (sum(sample_h0_stat <= symm_stat) + sum(sample_h0_stat >= trial_stat)) / sample_size
    }

    pvals <- c(pvals, pval)
  }

  return(pvals)
}


reject_stat <- function(trial_stats, i, stat_fun_data) {
  trial_stats <- trial_stats[i]
  num_trials <- length(trial_stats)
  sample_h0_stat <- stat_fun_data$h0
  expected_h0_stat <- stat_fun_data$exp_val
  sample_size <- length(sample_h0_stat)
  rejected_sl <- NULL

  pvals <- NULL
  for (trial_stat in trial_stats) {
    if (trial_stat < expected_h0_stat) {
      symm_stat <- expected_h0_stat + abs(expected_h0_stat - trial_stat)
      pval <- (sum(sample_h0_stat <= trial_stat) + sum(sample_h0_stat >= symm_stat)) / sample_size
    } else {
      symm_stat <- expected_h0_stat - abs(expected_h0_stat - trial_stat)
      pval <- (sum(sample_h0_stat <= symm_stat) + sum(sample_h0_stat >= trial_stat)) / sample_size
    }

    pvals <- c(pvals, pval)
  }

  significance_levels <- stat_fun_data$sl
  for (sl in significance_levels) {
    rejected <- sum(pvals <= sl) / num_trials
    rejected_sl <- c(rejected_sl, rejected)
  }

  return(rejected_sl)
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
significance_levels <- c(0.001, 0.005, seq(0.01, 1, by=0.1))

for (square_side in square_sides) {

  square_window <- owin(xrange=c(0, square_side), yrange=c(0, square_side))
  square_a <- square_side^2

  samples_info <- paste0("H0s", num_h0_sampling, "_H0t", num_h0_trials, "_B", num_boot, "_S", sprintf("%.2f", square_side))
  window_folder <- file.path(calibration_folder, samples_info)

  dir.create(window_folder, showWarnings = FALSE)

  for (num_points in nums_points) {
    tic()
    print(sprintf("PN %s", num_points))

    experiment_name <- paste0("IHD", inh_distance, "_PN", num_points)
    experiment_folder <- file.path(window_folder, experiment_name)

    dir.create(experiment_folder, showWarnings = FALSE)

    # Start analysis
    # Sample trial pps
    intensity <- num_points * 1.1 / square_a
    trial_pps <- sim_pp_pnum(num_points, rpoispp, num_h0_trials, intensity, win=square_window)

    output <- process_pps_lcf(trial_pps, num_points, intensity, square_window, inh_distance,
                             dim_lims, correction)

    trial_pps <- output$updated_pps
    # Vector of length num_h1
    trial_stat_lcf <- output$lcf_stat

    # Calculate test statistics for K-function, PCF, F-Function
    # Vector of length num_h1
    trial_stat_h <- summ_stat(trial_pps, inh_distance, h_at_r, correction=correction)
    # Vector of length num_h1
    trial_stat_pcf <- summ_stat(trial_pps, inh_distance, pcf_at_r, correction=correction)

    sample_h0 <- sim_pp_pnum(num_points, rpoispp, num_h0_sampling, intensity, win=square_window)

    output <- process_pps_lcf(sample_h0, num_points, intensity, square_window, inh_distance,
                             dim_lims, correction)

    sample_h0 <- output$updated_pps
    lcf_func_h0_stats <- output$lcf_stat

    h_func_h0_stats <- summ_stat(sample_h0, inh_distance, h_at_r, correction=correction)
    pcf_h0_stats <- summ_stat(sample_h0, inh_distance, pcf_at_r, correction=correction)

    h_data <- list(h0=h_func_h0_stats, sl=significance_levels, exp_val=0)
    h_boot_rej <- pval_boot_bca(trial_stat_h, reject_stat, h_data, B=num_boot)

    pcf_data <- list(h0=pcf_h0_stats, sl=significance_levels, exp_val=1)
    pcf_boot_rej <- pval_boot_bca(trial_stat_pcf, reject_stat, pcf_data, B=num_boot)

    lcf_data <- list(h0=lcf_func_h0_stats, sl=significance_levels, exp_val=0)
    lcf_boot_rej <- pval_boot_bca(trial_stat_lcf, reject_stat, lcf_data, B=num_boot)

    pvals_stat_h <- get_pvals(trial_stat_h, h_func_h0_stats)
    pvals_stat_pcf <- get_pvals(trial_stat_pcf, pcf_h0_stats, expected_h0_stat = 1)
    pvals_stat_lcf <- get_pvals(trial_stat_lcf, lcf_func_h0_stats)

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
    lcf_sampl_distr_df <- data.frame(test_stat=lcf_func_h0_stats)
    write.csv(lcf_sampl_distr_df, file.path(experiment_folder, sprintf("%s_H0_sampl_distr.csv", lcf_xlab)))
    lcf_sampl_distr_p <- sampl_distr_plot(lcf_sampl_distr_df, lcf_xlab)
    ggsave(file.path(experiment_folder, sprintf("%s_H0_sampl_distr.pdf", lcf_xlab)), lcf_sampl_distr_p, width=5, height=4)


    calibration_tb <- build_calibration_tb(significance_levels,
                                           h_boot_rej,
                                           pcf_boot_rej,
                                           lcf_boot_rej)
    write.csv(calibration_tb, file.path(experiment_folder, sprintf("Calibration_S%d_PN%d.csv", floor(square_side), num_points)))

    calibration_p <- calibration_plot(calibration_tb)

    ggsave(file.path(experiment_folder, "Calibration.pdf"), calibration_p, width=5, height=4)


    # Show p-values qq plots together
    pvals_tb <- build_pvals_tb(pvals_stat_h,
                               pvals_stat_pcf,
                               pvals_stat_lcf)
    write.csv(pvals_tb, file.path(experiment_folder, "Pvals_qq.csv"))

    pvals_p <- qqplot_pvals_together(pvals_tb)

    ggsave(file.path(experiment_folder, "Pvals_qq.pdf"), pvals_p, width=5, height=4)

    time_elapsed <- toc()
    stat_file_path <- file.path(experiment_folder, "stat.txt")
    file_conn <- file(stat_file_path)
    lines <- c(
      paste("time", as.numeric(time_elapsed$toc - time_elapsed$tic), sep = "\t"),
      paste("num_LCFest_failed", as.numeric(num_LCFest_failed), sep = "\t")
    )
    writeLines(lines, file_conn)
    close(file_conn)
  }
}
