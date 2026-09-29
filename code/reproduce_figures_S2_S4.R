# Run from package root: Rscript --vanilla code/reproduce_figures_S2_S4.R
source("code/00_setup.R")
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
lines <- readLines("code/supplementary_information.R")
run_block <- function(start, end) {
  a <- grep(start, lines); b <- grep(end, lines)
  stopifnot(length(a) == 1L, length(b) == 1L, b > a)
  eval(parse(text = lines[a:(b - 1)]), envir = .GlobalEnv)
}
run_block("^#figure s2", "^#figure s3")
run_block("^#figure s4", "^#figure s5")
