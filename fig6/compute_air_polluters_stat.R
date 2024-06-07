rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggbeeswarm)
library(lme4)

# File names
results_folder <- "lcf_states"
stat_file_name <- "state_stat.csv"
lcf_stat_prefix <- "lcf_stat"

output_folder <- "data_to_plot"
dir.create(output_folder, showWarnings = FALSE, recursive = FALSE)

min_sample_size <- 40


# Get list of states folders
states <- list.dirs(results_folder, full.names = FALSE, recursive = FALSE)

get_lcf_integral <- function(lcf_curve, r_first, r_last) {

  norm_c <- lcf_curve$r[nrow(lcf_curve)]

  r_first <- r_first / norm_c
  r_last <- r_last / norm_c

  lcf_curve <- lcf_curve %>%
    mutate(r = r / norm_c) %>%
    filter(r >= r_first & r <= r_last)

  N <- nrow(lcf_curve) - 1
  lcf_integral <- 0
  for (i in 1:N) {
    a <- lcf_curve$r[i]
    b <- lcf_curve$r[i+1]

    f_a <- lcf_curve$val[i]
    f_b <- lcf_curve$val[i+1]

    lcf_integral <- lcf_integral + (b - a) * (f_a + f_b) / 2
  }

  return(lcf_integral)
}

# Stat parameters
r_first <- 7500
r_last <- 75000

states_df <- NULL
ftypes_df <- NULL
facility_nums_df <- NULL

lcf_integrals_df <- NULL

# Get list of Scan's folders
states <- list.dirs(results_folder, full.names = FALSE, recursive = FALSE)

for (state in states) {

  state_folder <- file.path(results_folder, state)
  stat_file_path <- file.path(state_folder, stat_file_name)

  state_stat <- read.csv(stat_file_path)
  state_stat_filt <- state_stat %>%
    filter(n >= min_sample_size)

  facility_types_num <- nrow(state_stat_filt)
  if (facility_types_num > 0) {
    for (ind in 1:facility_types_num) {
      facility_type <- state_stat_filt$type[ind]
      facility_num <- state_stat_filt$n[ind]

      # To get clustering metric, average from 100 to 500
      ftype_lcf_curve_file <- sprintf("%s_%s.csv", lcf_stat_prefix, facility_type)
      ftype_lcf_curve_path <- file.path(state_folder, ftype_lcf_curve_file)

      if (!file.exists(ftype_lcf_curve_path)) {
        print(sprintf("%s does not exist!", ftype_lcf_curve_path))
        next
      }

      ftype_lcf_curve <- read.csv(ftype_lcf_curve_path)

      lcf_integral <- get_lcf_integral(ftype_lcf_curve, r_first, r_last)

      states_df <- c(states_df, state)
      ftypes_df <- c(ftypes_df, facility_type)
      facility_nums_df <- c(facility_nums_df, facility_num)

      lcf_integrals_df <- c(lcf_integrals_df, lcf_integral)
    }
  }
}


lcf_stat_df <- data.frame(state=states_df,
                         ftype=ftypes_df,
                         facility_num=facility_nums_df,
                         lcf_integral=lcf_integrals_df)

lcf_stat_file_path <- file.path(output_folder, "lcf_stat_facilities.csv")
write.csv(lcf_stat_df, lcf_stat_file_path)
