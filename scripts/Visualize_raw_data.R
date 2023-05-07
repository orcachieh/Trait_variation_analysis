#+ Read library, include = FALSE
library(tidyverse)
library(ggplot2)

#+ Read data and Pivot to plotting format
trait_dat <- read.csv("../Raw_data.csv", header = TRUE)

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

ggplot(data = trait_long_dat[trait_long_dat$abb=="HUN", ]) +
  geom_density(aes(x=Leafwidth))
