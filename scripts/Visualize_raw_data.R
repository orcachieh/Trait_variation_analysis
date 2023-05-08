#+ Read library, include = FALSE
library(tidyverse)
library(ggplot2)
library(ggpubr)

#' Read trait measured data and rearanging it to plotting format.

#+ Read data and Pivot to plotting format, warning = FALSE

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
head(trait_long_dat)

#' plot the potential intraspecific trait variation of each species.

#+ Plotting, echo = FALSE, warning = FALSE
HU <- trait_long_dat[trait_long_dat$abb=="HUN"|trait_long_dat$abb=="HUW", ]
HO <- trait_long_dat[trait_long_dat$abb=="HO", ]
ZC <- trait_long_dat[trait_long_dat$abb=="ZC", ]
CS <- trait_long_dat[trait_long_dat$abb=="CS", ]
HS <- trait_long_dat[trait_long_dat$abb=="HS", ]
TH <- trait_long_dat[trait_long_dat$abb=="TH", ]

trait_plot = function(sp){
  SN = ggplot(data = sp) +
    geom_boxplot(aes(y = Shoot_number), width = 0.7) +
    xlim(-1, 1) +
    ylab("Shoot number per core") +
    xlab("") +
    scale_x_continuous(labels = NULL, breaks = NULL) +
    theme(legend.position = "none")
    
  sp %>% 
    mutate(MAB = rowMeans(cbind(Aboveground_dry_mass1, Aboveground_dry_mass2)),
           MBB = rowMeans(cbind(Belowground_dry_mass1, Belowground_dry_mass2))) %>%
    pivot_longer(cols = c(MAB, MBB), names_to = "DryM", values_to = "Weight")-> sp_DM
  
  DM = ggplot(data = sp_DM) +
    geom_boxplot(aes(y = Weight, fill = DryM)) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +
    ylab("Dry weight (mg)") +
    xlab("") +
    scale_x_continuous(labels = c("Aboveground", "Belowground"), 
                       breaks = c(-0.2, 0.2)) +
    theme(legend.position = "none")
  
  LH = ggplot(data = sp) +
    geom_histogram(aes(x = Leaflength), alpha = 0.5) +
    labs(x = "Leaf length (mm)")
  
  LW = ggplot(data = sp) +
    geom_histogram(aes(x = Leafwidth), alpha = 0.5) +
    labs(x = "Leaf width (mm)")
  
  CH = ggplot(data = sp) +
    geom_histogram(aes(x = Canopyheight), alpha = 0.5) +
    labs(x = "Canopy height (mm)")
  
  RL = ggplot(data = sp) +
    geom_histogram(aes(x = Rootlength), alpha = 0.5) +
    labs(x = "Root length (mm)")
  
  RD = ggplot(data = sp) +
    geom_histogram(aes(x = Rhizomediameter), alpha = 0.5) +
    labs(x = "Rhizome diameter (mm)")
  
  IL = ggplot(data = sp) +
    geom_histogram(aes(x = Internodelength), alpha = 0.5) +
    labs(x = "Internode length (mm)")
  
  p_all = ggarrange(SN, DM, LH, LW, CH, RL, RD, IL,
                    ncol = 2, nrow = 4,
                    common.legend = TRUE, legend = "none")
  return(p_all)
}

#' Holudule uninervis traits
#+ HU, warning = FALSE, message = FALSE
trait_plot(HU)

#' Halophila ovalis traits
#+ HO, warning = FALSE, message = FALSE
trait_plot(HO)

#' Zostera muelleri subsp. Capricorni traits
#+ ZC, warning = FALSE, message = FALSE
trait_plot(ZC)

#' Cymodocea serrulata traits
#+ CS, warning = FALSE, message = FALSE
trait_plot(CS)

#' Halophila spinulosa traits
#+ HS, warning = FALSE, message = FALSE
trait_plot(HS)

#' Thalassia hemprichii traits
#+ TH, warning = FALSE, message = FALSE
trait_plot(TH)