rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggbeeswarm)
library(lme4)

# File names
results_folder <- "lcf_tmas"
stat_file_name <- "cell_stat.csv"
lcf_stat_prefix <- "ai_stat"

# Constants
min_cell_num <- 40

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
r_first <- 50
r_last <- 500

slides_df <- NULL
patients_df <- NULL
tmas_df <- NULL
ctypes_df <- NULL
cell_nums_df <- NULL

lcf_integrals_df <- NULL

# Get list of Scan's folders
slides <- list.dirs(results_folder, full.names = FALSE, recursive = FALSE)

for (slide in slides) {

  slide_folder <- file.path(results_folder, slide)
  patients <- list.dirs(slide_folder, full.names = FALSE, recursive = FALSE)
  for (patient in patients) {

    patient_folder <- file.path(slide_folder, patient)
    tma_nums <- list.dirs(patient_folder, full.names = FALSE, recursive = FALSE)
    for (tma_num in tma_nums) {

      tma_folder <- file.path(patient_folder, tma_num)
      stat_file_path <- file.path(tma_folder, stat_file_name)

      tma_stat <- read.csv(stat_file_path)
      tma_stat_filt <- tma_stat %>%
        filter(n >= min_cell_num)

      cell_types_num <- nrow(tma_stat_filt)
      if (cell_types_num > 0) {
        for (ind in 1:cell_types_num) {
          cell_type <- tma_stat_filt$type[ind]
          cell_num <- tma_stat_filt$n[ind]

          # To get clustering metric, average from 100 to 500
          ctype_lcf_curve_file <- sprintf("%s_%s.csv", lcf_stat_prefix, cell_type)
          ctype_lcf_curve_path <- file.path(tma_folder, ctype_lcf_curve_file)

          if (!file.exists(ctype_lcf_curve_path)) {
            print(sprintf("%s does not exist!", ctype_lcf_curve_path))
            next
          }

          ctype_lcf_curve <- read.csv(ctype_lcf_curve_path)

          lcf_integral <- get_lcf_integral(ctype_lcf_curve, r_first, r_last )

          slides_df <- c(slides_df, slide)
          patients_df <- c(patients_df, patient)
          tmas_df <- c(tmas_df, tma_num)
          ctypes_df <- c(ctypes_df, cell_type)
          cell_nums_df <- c(cell_nums_df, cell_num)

          lcf_integrals_df <- c(lcf_integrals_df, lcf_integral)
        }
      }
    }
  }
}


lcf_stat_df <- data.frame(slide=slides_df,
                         patient=patients_df,
                         tma=tmas_df,
                         ctype=ctypes_df,
                         cell_num=cell_nums_df,
                         lcf_integral=lcf_integrals_df)


lcf_stat_file_path <- "lcf_stat_tmas.csv"
write.csv(lcf_stat_df, lcf_stat_file_path)
