rm(list=ls())

library(dplyr)
library(spatstat)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(lcfstat)
library(tictoc)
library(boot)

source("../utils_lcf.R")

sampl_distr_plot <- function(null_test_stat_df, xlab) {
  null_test_stat_p <- ggplot(null_test_stat_df, aes(test_stat)) +
    geom_histogram() +
    ylab("Count") +
    xlab(xlab) +
    theme_classic()

  null_test_stat_p
}

build_calibration_tb <- function(significance_levels,
                                 stat_rej_h,
                                 stat_rej_pcf,
                                 stat_rej_lcf) {

  stat_col <- factor(c(rep("h", nrow(stat_rej_h)),
                       rep("pcf", nrow(stat_rej_h)),
                       rep("lcf", nrow(stat_rej_h))), levels=c("h", "pcf", "lcf"))

  calibration_df <- rbind(stat_rej_h, stat_rej_pcf, stat_rej_lcf) %>%
    mutate(stat=stat_col, sl=rep(significance_levels, 3))

  calibration_df
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
    theme_classic()

  calibration_p
}

build_pvals_tb <- function(stat_pvals_h, stat_pvals_pcf, stat_pvals_lcf) {
  pvals_observed <- c(sort(stat_pvals_h),
                      sort(stat_pvals_pcf),
                      sort(stat_pvals_lcf))

  pvals_theo <- rep(seq(0, 1, length.out=length(stat_pvals_h)), 3)

  stat_col <- factor(c(rep("h", length(stat_pvals_h)),
                       rep("pcf", length(stat_pvals_pcf)),
                       rep("lcf", length(stat_pvals_lcf))), levels=c("h", "pcf", "lcf"))

  pvals_qq_df <- tibble(observed = pvals_observed, theo=pvals_theo, stat=stat_col)

  pvals_qq_df
}

qqplot_pvals <- function(pvals_qq_df) {

  breaks <- c(10^(-3), 10^(-2), 10^(-1), 10^0)

  pvals_qq_p <- ggplot(pvals_qq_df) +
    geom_line(aes(x=observed, y=theo, col=stat), linewidth=1.1) +
    scale_colour_manual(values = c(h_color, pcf_color, lcf_color)) +
    geom_abline(intercept=0, slope=1, color="black", linetype = "dashed") +
    ylab("P-values theorethical") +
    xlab("P-values empirical") +
    coord_fixed() +
    theme_classic() +
    scale_x_log10(breaks = breaks, labels = breaks) +
    scale_y_log10(breaks = breaks, labels = breaks)

  pvals_qq_p
}

get_pvals <- function(trial_stats, h0_stat, expected_h0_stat=0) {
  num_trials <- length(trial_stats)
  sample_size <- length(h0_stat)

  pvals <- NULL
  for (trial_stat in trial_stats) {
    if (trial_stat < expected_h0_stat) {
      symm_stat <- expected_h0_stat + abs(expected_h0_stat - trial_stat)
      pval <- (sum(h0_stat <= trial_stat) + sum(h0_stat >= symm_stat)) / sample_size
    } else {
      symm_stat <- expected_h0_stat - abs(expected_h0_stat - trial_stat)
      pval <- (sum(h0_stat <= symm_stat) + sum(h0_stat >= trial_stat)) / sample_size
    }

    pvals <- c(pvals, pval)
  }

  pvals
}


reject_stat <- function(trial_stats, i, stat_fun_data) {
  trial_stats <- trial_stats[i]
  num_trials <- length(trial_stats)
  sample_h0_stat <- stat_fun_data$h0
  expected_h0_stat <- stat_fun_data$exp_val
  sample_size <- length(sample_h0_stat)

  pvals <- get_pvals(trial_stats, sample_h0_stat, expected_h0_stat)

  get_rej_rate <- function(x) sum(pvals <= x) / num_trials

  significance_levels <- stat_fun_data$sl
  rejected_sl <- map_dbl(significance_levels, get_rej_rate)

  rejected_sl
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
        boot_stat
      }
    )
    boot_stats <- rbind(boot_stats, result)
  }

  boot_stats
}

set.seed(2711)

calibration_folder <- "Calibration"
dir.create(calibration_folder, showWarnings = FALSE)

h_color <- "#D6984D"
pcf_color <- "#008F0C"
lcf_color <- "#4A16DB"

# Window square
square_sides <- c(200 / sqrt(2), 200, 200 * sqrt(2))

# Methods params
correction <- "Ripley"
dim_lims <- c(4, 50)

# Params to change
inh_distance <- 5
nums_points <- c(25, 50, 100)

# Monte Carlo test parameters
num_h0_sampling <- 50
num_h0_trials <- 50
num_boot <- 500

# Track how many times LCF computation has failed
num_LCFest_failed <- 0

significance_levels <- c(0.001, 0.005, seq(0.01, 1, by=0.1))

for (square_side in square_sides) {
  print(sprintf("Square %dx%d", round(square_side), round(square_side)))

  square_window <- make_square_win(square_side)
  square_a <- square_side^2

  samples_info <- paste0("H0s", num_h0_sampling, "_H0t", num_h0_trials, "_B", num_boot, "_S", sprintf("%.2f", square_side))
  window_folder <- file.path(calibration_folder, samples_info)
  dir.create(window_folder, showWarnings = FALSE)

  for (num_points in nums_points) {
    tic()
    print(sprintf("PN %s", num_points))

    test_name <- paste0("IHD", inh_distance, "_PN", num_points)
    test_folder <- file.path(window_folder, test_name)
    dir.create(test_folder, showWarnings = FALSE)

    # Start analysis
    # Simulate H0 distribution
    sample_h0 <- runifpoint(num_points, win=square_window, nsim=num_h0_sampling)

    # Calculate test statistics for H-function, PCF and LCF
    output <- process_pps_lcf(sample_h0, runifpoint, inh_distance, dim_lims, correction,
                              n = num_points,  win = square_window)

    sample_h0 <- output$updated_pps
    lcf_func_h0_stats <- output$lcf_stat

    h_func_h0_stats <- summ_stat(sample_h0, inh_distance, h_at_r, correction=correction)
    pcf_h0_stats <- summ_stat(sample_h0, inh_distance, pcf_at_r, correction=correction)

    # Sample trial point patterns
    trial_pps <- runifpoint(num_points, win=square_window, nsim=num_h0_trials)

    # Calculate test statistics for H-function, PCF and LCF
    output <- process_pps_lcf(trial_pps, runifpoint, inh_distance, dim_lims, correction,
                              n = num_points,  win = square_window)

    trial_pps <- output$updated_pps
    trial_stat_lcf <- output$lcf_stat

    trial_stat_h <- summ_stat(trial_pps, inh_distance, h_at_r, correction=correction)
    trial_stat_pcf <- summ_stat(trial_pps, inh_distance, pcf_at_r, correction=correction)

    # Get probability of H0 rejection for each summary function
    h_data <- list(h0=h_func_h0_stats, sl=significance_levels, exp_val=0)
    h_boot_rej <- pval_boot_bca(trial_stat_h, reject_stat, h_data, B=num_boot)

    pcf_data <- list(h0=pcf_h0_stats, sl=significance_levels, exp_val=1)
    pcf_boot_rej <- pval_boot_bca(trial_stat_pcf, reject_stat, pcf_data, B=num_boot)

    lcf_data <- list(h0=lcf_func_h0_stats, sl=significance_levels, exp_val=0)
    lcf_boot_rej <- pval_boot_bca(trial_stat_lcf, reject_stat, lcf_data, B=num_boot)

    # Get p-values
    pvals_stat_h <- get_pvals(trial_stat_h, h_func_h0_stats)
    pvals_stat_pcf <- get_pvals(trial_stat_pcf, pcf_h0_stats, expected_h0_stat = 1)
    pvals_stat_lcf <- get_pvals(trial_stat_lcf, lcf_func_h0_stats)

    # Visualize the sampling distribution of test statistics under H0
    h_xlab <- sprintf("H(%.2f)", inh_distance)
    h_sampl_distr_df <- data.frame(test_stat=h_func_h0_stats)
    write_csv(h_sampl_distr_df, file.path(test_folder, sprintf("%s_H0_sampl_distr.csv", h_xlab)))
    h_sampl_distr_p <- sampl_distr_plot(h_sampl_distr_df, h_xlab)
    ggsave(file.path(test_folder, sprintf("%s_H0_sampl_distr.pdf", h_xlab)), h_sampl_distr_p, width=5, height=4)

    pcf_xlab <- sprintf("PCF(%.2f)", inh_distance)
    pcf_sampl_distr_df <- data.frame(test_stat=pcf_h0_stats)
    write_csv(pcf_sampl_distr_df, file.path(test_folder, sprintf("%s_H0_sampl_distr.csv", pcf_xlab)))
    pcf_sampl_distr_p <- sampl_distr_plot(pcf_sampl_distr_df, pcf_xlab)
    ggsave(file.path(test_folder, sprintf("%s_H0_sampl_distr.pdf", pcf_xlab)), pcf_sampl_distr_p, width=5, height=4)

    lcf_xlab <- sprintf("LCF(%.2f)", inh_distance)
    lcf_sampl_distr_df <- data.frame(test_stat=lcf_func_h0_stats)
    write_csv(lcf_sampl_distr_df, file.path(test_folder, sprintf("%s_H0_sampl_distr.csv", lcf_xlab)))
    lcf_sampl_distr_p <- sampl_distr_plot(lcf_sampl_distr_df, lcf_xlab)
    ggsave(file.path(test_folder, sprintf("%s_H0_sampl_distr.pdf", lcf_xlab)), lcf_sampl_distr_p, width=5, height=4)


    # Save the results of the calibration analysis as a .csv and as a plot
    calibration_tb <- build_calibration_tb(significance_levels,
                                           h_boot_rej,
                                           pcf_boot_rej,
                                           lcf_boot_rej)
    write_csv(calibration_tb, file.path(test_folder, sprintf("Calibration_S%d_PN%d.csv", floor(square_side), num_points)))

    calibration_p <- calibration_plot(calibration_tb)

    ggsave(file.path(test_folder, "Calibration.pdf"), calibration_p, width=5, height=4)


    # Save p-values in a .csv and as a plot
    pvals_tb <- build_pvals_tb(pvals_stat_h,
                               pvals_stat_pcf,
                               pvals_stat_lcf)
    write_csv(pvals_tb, file.path(test_folder, "Pvals_qq.csv"))

    pvals_p <- qqplot_pvals(pvals_tb)

    ggsave(file.path(test_folder, "Pvals_qq.pdf"), pvals_p, width=5, height=4)

    time_elapsed <- toc()
    # Write a file with running time
    stat_file_path <- file.path(test_folder, "stat.txt")
    file_conn <- file(stat_file_path)
    lines <- c(
      paste("time", as.numeric(time_elapsed$toc - time_elapsed$tic), sep = "\t"),
      paste("num_LCFest_failed", as.numeric(num_LCFest_failed), sep = "\t")
    )
    writeLines(lines, file_conn)
    close(file_conn)
  }
}
