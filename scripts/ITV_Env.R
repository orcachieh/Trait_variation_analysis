#' ---
#'  title: "Seagrass intraspecific trait variation and evironemtal variables"
#'  author: "Chieh"
#'  output: 
#'    html_document:
#'      toc: true
#'      toc_float: true
#'      toc_depth: 3 
#'      theme: lumen
#'      highlight: tango
#' ---

#+  Read library, include = FALSE
#### Read Libaray ####
library(tidyverse)
library(rethinking)
library(brms)
library(lme4)

#+ Source functions from Functions.R 
#### Input functions  ####

source("./Functions.R")


#### Read all sample localities ####
#' # All sample localities
# read and cleanup sample point csv

trait_dat <- read.csv("../trait_raw_data.csv", header = TRUE)
env_var <- read.csv("../environment_variables.csv", header = TRUE)

trait_dat %>%
  mutate(Belowground_dry_mass1 = as.numeric(Belowground_dry_mass1)) %>%
  select(!ends_with("photo.") & !contains("note") &
           !Pack & 
           !Aboveground_dry_mass3 &
           !Belowground_dry_mass3 &
           !Flower_seed) %>%
  pivot_longer(cols = Canopyheight_1:Internodelength_10,
               names_to = c(".value", "trials"),
               names_sep = "_" ) -> trait_long_dat
head(trait_long_dat)