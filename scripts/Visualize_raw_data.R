#' ---
#'  title: "Visualize intrasepecific trait variation measurements"
#'  author: "Chieh"
#'  output: 
#'    html_document:
#'      toc: true
#'      toc_float: true
#'      toc_depth: 3 
#'      theme: lumen
#'      highlight: tango
#' ---


#+ Read library, include = FALSE
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(plotly)
library(grid)

#' This document aims to show trait variation within a seagrass species.
#'The raw data were obtained from 12 distinct locations in Queensland, Australia. The measurements primarily took place within a controlled wet laboratory environment, employing appropriate scientific instruments. However, a subset of finer-scale measurements, typically measuring less than 5 mm, were conducted digitally using the ImageJ software on a computer
#'
#' To illustrate the distribution of the traits, two graphical methods were employed: box plots and histograms. The box plots were utilized to depict the shoot number per core and the Belowground/Aboveground ratio. On the other hand, histograms were employed to display the distribution of six additional traits, namely Leaf Length, Leaf Width, Canopy Height, Root Length, Rhizome Diameter, and Internode Distance.
#' 


#' # Data
#' Read trait measured data and rearranging it to plotting format.

#+ Read data and Pivot to plotting format, warning = FALSE

trait_dat <- read.csv("../trait_raw_data.csv", header = TRUE)

trait_dat %>%
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

#+ Plotting function for all species, echo = FALSE, warning = FALSE
HU <- trait_long_dat[trait_long_dat$abb=="HUN"|trait_long_dat$abb=="HUW", ]
HO <- trait_long_dat[trait_long_dat$abb=="HO", ]
ZC <- trait_long_dat[trait_long_dat$abb=="ZC", ]
CS <- trait_long_dat[trait_long_dat$abb=="CS", ]
HS <- trait_long_dat[trait_long_dat$abb=="HS", ]
TH <- trait_long_dat[trait_long_dat$abb=="TH", ]

trait_plot = function(sp, LHbin, LWbin, CHbin, RLbin, RDbin, ILbin){
  SN = ggplot(data = sp) +
    geom_boxplot(aes(y = Shoot_number), width = 0.7) +
    xlim(-1, 1) +
    ylab("Shoot #/core") +
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
    geom_histogram(aes(x = Leaflength), 
                   binwidth = LHbin, alpha = 0.5) +
    labs(x = "Leaf length (mm)")
  
  LW = ggplot(data = sp) +
    geom_histogram(aes(x = Leafwidth), 
                   binwidth = LWbin, alpha = 0.5) +
    labs(x = "Leaf width (mm)")
  
  CH = ggplot(data = sp) +
    geom_histogram(aes(x = Canopyheight), 
                   binwidth = CHbin, alpha = 0.5) +
    labs(x = "Canopy height (mm)")
  
  RL = ggplot(data = sp) +
    geom_histogram(aes(x = Rootlength), 
                   binwidth = RLbin, alpha = 0.5) +
    labs(x = "Root length (mm)")
  
  RD = ggplot(data = sp) +
    geom_histogram(aes(x = Rhizomediameter),
                   binwidth = RDbin, alpha = 0.5) +
    labs(x = "Rhizome diameter (mm)")
  
  IL = ggplot(data = sp) +
    geom_histogram(aes(x = Internodelength),
                   binwidth = ILbin, alpha = 0.5) +
    labs(x = "Internode length (mm)")
  
  p_all = ggarrange(SN, DM, LH, LW, CH, RL, RD, IL,
                    ncol = 2, nrow = 4,
                    common.legend = TRUE, legend = "none")
  return(p_all)
}

#' # Visualize trait distribution

#' ## Holudule uninervis traits
#+ HU, warning = FALSE, message = FALSE
trait_plot(HU, 4, 0.05, 4, 4, 0.04, 1)

#' ## Halophila ovalis traits
#+ HO, warning = FALSE, message = FALSE
trait_plot(HO, 0.5, 0.15, 0.75, 2, 0.02, 1)

#' ## Zostera muelleri subsp. Capricorni traits
#+ ZC, warning = FALSE, message = FALSE
trait_plot(ZC, 4, 0.2, 4, 3, 0.02, 0.5)

#' ## Cymodocea serrulata traits
#+ CS, warning = FALSE, message = FALSE
trait_plot(CS, 3, 0.1, 3, 3, 0.02, 1.2)

#' ## Halophila spinulosa traits
#+ HS, warning = FALSE, message = FALSE
trait_plot(HS, 0.25, 0.025, 2, 1.5, 0.02, 0.5)

#' ## Thalassia hemprichii traits
#+ TH, warning = FALSE, message = FALSE
trait_plot(TH, 1.5, 0.08, 1.5, 1, 0.05, 0.7)

#+ Plotting function for HUN vs HUW, echo = FALSE, warning = FALSE
trait_plot_HU = function(sp){
  SN = ggplot(data = sp) +
    geom_boxplot(aes(x = abb, y = Shoot_number, fill = abb), width = 0.7) +
    ylab("Shoot #/core") +
    xlab("") +
    theme(legend.position = "none")
  
  sp %>% 
    mutate(MAB = rowMeans(cbind(Aboveground_dry_mass1, Aboveground_dry_mass2)),
           MBB = rowMeans(cbind(Belowground_dry_mass1, Belowground_dry_mass2))) %>%
    pivot_longer(cols = c(MAB, MBB), names_to = "DryM", values_to = "Weight")-> sp_DM
  
  DM = ggplot(data = sp_DM) +
    geom_boxplot(aes(x = DryM, y = Weight, fill = abb)) +
    ylab("Dry weight (mg)") +
    xlab("") +
    scale_x_discrete(labels = c("Aboveground", "Belowground")) +
    theme(legend.position = "none")
  
  LH = ggplot(data = sp) +
    geom_histogram(aes(x = Leaflength, fill = abb),
                   binwidth = 4, alpha = 0.5) +
    labs(x = "Leaf length (mm)")
  
  LW = ggplot(data = sp) +
    geom_histogram(aes(x = Leafwidth, fill = abb),
                   binwidth = 0.05, alpha = 0.5) +
    labs(x = "Leaf width (mm)")
  
  CH = ggplot(data = sp) +
    geom_histogram(aes(x = Canopyheight,fill = abb),
                   binwidth = 4, alpha = 0.5) +
    labs(x = "Canopy height (mm)")
  
  RL = ggplot(data = sp) +
    geom_histogram(aes(x = Rootlength, fill = abb), 
                   binwidth = 4, alpha = 0.5) +
    labs(x = "Root length (mm)")
  
  RD = ggplot(data = sp) +
    geom_histogram(aes(x = Rhizomediameter, fill = abb), 
                   binwidth = 0.04, alpha = 0.5) +
    labs(x = "Rhizome diameter (mm)")
  
  IL = ggplot(data = sp) +
    geom_histogram(aes(x = Internodelength, fill = abb), 
                   binwidth = 1, alpha = 0.5) +
    labs(x = "Internode length (mm)")
  
  p_all = ggarrange(SN, DM, LH, LW, CH, RL, RD, IL,
                    ncol = 2, nrow = 4,
                    common.legend = TRUE, legend = "none")
  return(p_all)
}

#' # HUN vs HUW
#' ## Trait of narrow versus wide form of Haludule uninervis
#+ HUN vs HUW, warning = FALSE, message = FALSE, echo = FALSE
trait_plot_HU(HU)
#' The clear difference is at Leaf width. The rhizome diameter is also quite distinct for two form but with a overlap zone.
#' Leaf width indicates there might be three groups instead of two.
#' 

#' ## Location specific leaf width of HUN and HUW.
#' These graph indicate potential groups in HU.
#+ HU Location, warning = FALSE, message = FALSE, echo = FALSE

ggplot(data = HU[HU$abb == "HUN", ], aes(x = Location, y = Leafwidth, fill = Location)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "plasma") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x=element_text(size = 8, angle=45, hjust=1)
  ) +
  ggtitle("HU narrow leaf width at each sample location") +
  xlab("")

ggplot(data = HU[HU$abb == "HUW", ], aes(x = Location, y = Leafwidth, fill = Location)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x=element_text(size = 8, angle=45, hjust=1)
  ) +
  ggtitle("HU wide leaf width at each sample location") +
  xlab("")

ggplot(data = HU, aes(x = Location, y = Leafwidth, fill = Location)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.05, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x=element_text(size = 8, angle=45, hjust=1)
  ) +
  ggtitle("HU leaf width at each sample location") +
  xlab("")

#' Plot HU trait along the latitude gradient
#+ HU latitude gradient, warning = FALSE, message = FALSE, echo = FALSE

ggplot(data = HU, aes(x = Latitude, y = Leafwidth)) +
  geom_point()

ggplot(data = HU, aes(y = Latitude, x = Longitude, size = Leafwidth, color = Leafwidth)) +
  geom_point()

ggplot(data = HU[HU$abb == "HUN", ], aes(x = Latitude, y = Leafwidth)) +
  geom_point()

ggplot(data = HU[HU$abb == "HUN", ], aes(x = Longitude, y = Leafwidth)) +
  geom_point()

plot_ly(data = HU[HU$abb == "HUN", ], 
        x = ~Longitude, y = ~Latitude, z = ~Leafwidth, 
        type = "scatter3d", mode = "markers")


#' HO leaf surface area distribution.
#+ HO leaf surface area, warning = FALSE, message = FALSE, echo = FALSE

HO %>%
  mutate(LSA = Leaflength * Leafwidth) %>%
  ggplot(aes(x = LSA)) +
  geom_histogram(binwidth = 5) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("HO leaf surface area distribution") +
  xlab("Leaf surface area (mm^2)")

#' Tried to see is there any distinct groups of leaf size of HO, because I felt there a big leaf individuals and small leaf ones. However, it seems to be a continous distibution. 

#' HU triat plot for supplementary information
#+ HU triat plot for supplementary figure
HU %>%
  mutate(ABMs = rowMeans(select(.,Aboveground_dry_mass1, Aboveground_dry_mass2), na.rm = TRUE)/Shoot_number,
         BBMs = rowMeans(select(.,Belowground_dry_mass1, Belowground_dry_mass2))/Shoot_number, .keep = "unused") -> HU_splot


box_gf = function(trait, ytitle){
  HU_splot %>% 
    select(SampleID, abb, trait) %>%
    group_by(SampleID, abb) %>%
    dplyr::summarise(nn = mean(get(trait), na.rm = TRUE),
                     .groups = "drop") -> t_data
  
  
  pp <- ggplot(data = t_data) +
    geom_boxplot(aes_string(y = "abb" , x = colnames(t_data)[3], fill = "abb")) +
    scale_fill_manual(values = c(HUN = "#009E73", HUW = "#E69F00")) +
    theme_minimal() +
    labs(fill = "Growth forms",
         y = "",
         x = ytitle) +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.position = "right",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10))
}


Abiomass <- box_gf("ABMs", "Above ground biomass per shoot (mg)")
Bbiomass <- box_gf("BBMs", "Below ground biomass per shoot (mg)")
CH <- box_gf("Canopyheight", "Canopy Height (mm)")
LL <- box_gf("Leaflength", "Leaf Length (mm)")
LW <- box_gf("Leafwidth", "Leaf Width (mm)")
RL <- box_gf("Rootlength", "Root Length (mm)")
RL <- box_gf("Rootlength", "Root Length (mm)")
RD <- box_gf("Rhizomediameter", "Rhizome Diameter (mm)")
ID <- box_gf("Internodelength", "Internode Distance (mm)")


ggarrange(Abiomass, Bbiomass, LL, LW, CH, RL, RD, ID,
          ncol = 2, nrow = 4,
          common.legend = TRUE, legend = "none") + 
  annotation_custom(grid.polygon(x = c(0.5, 0.5, 0, 1, 0, 1, 0, 1),
                                 y = c(0, 1, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75),
                                 id = c(1, 1, 2, 2, 3, 3, 4, 4))) -> Splot
ggexport(Splot, filename = "../plots/Splot.png",
         width =4500, height = 5500, res = 600)

