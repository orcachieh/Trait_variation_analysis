#' ---
#'  title: "Seagrass intraspecific trait variation and environemtal variables"
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
library(rstan)
library(MASS)
library(rethinking)
library(mvtnorm) # produce multinomial model simulation
library(tidybayes)
library(bayesplot) # MCMC diagnostics
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(pcaMethods) # Bayesian PCA, enable missing value impute
# Plotting
library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(gridExtra) # arranging plots
library(psych) # pairwise scatter plot
library(ggrepel)
library(ggpp)


#+ Source functions from Functions.R 
#### Input functions  ####

source("./scripts/Functions.R")


#' # Trait and Environmental Variables 

#+ read and process data
#### Read and match trait and env. variables ####

trait_dat <- read.csv("./trait_raw_data.csv", header = TRUE)
env <- read.csv("./environment_variables.csv", header = TRUE)

# Multiple 6 location environmental data to 3 rows to match the trait_dat format (ZC, CV, HER, KAR, GLD, WP)
env %>%
  add_row(env[33:38, ]) %>%
  add_row(env[33:38, ]) -> env_var

# change location name to match smapleID in trait_dat
env_var$location[13] = "PAL5_09102022"
env_var$location[17] = "PALSUB1_24092022"
env_var$location[18] = "PALSUB2_24092022"
env_var$location[33:38] = c("ZC1_08102022", "CV1_08102022", "HER1_26102022", "KAR1_01112022", "GLD1_11102022", "WP1_27102022")
env_var$location[39:44] = c("ZC2_08102022", "CV2_08102022", "HER2_26102022", "KAR2_01112022", "GLD2_11102022", "WP2_27102022")
env_var$location[45:50] = c("ZC3_08102022", "CV3_08102022", "HER3_26102022", "KAR3_01112022", "GLD3_11102022", "WP3_27102022")

# Create a long table with long format for repeated measurement variables (e.g. Leaf width, root length...)
trait_dat %>%
  # match trait data with environmental variables based on sample location 
  left_join(env_var, by = join_by("SampleID" == "location")) %>%
  # remove descriptive variables
  dplyr::select(!ends_with("photo.") & !contains("note") &
           !Pack & 
           !Aboveground_dry_mass3 &
           !Belowground_dry_mass3 &
           !Flower_seed) %>%
  # change variables with multiple measurement to long tables
  pivot_longer(cols = Canopyheight_1:Internodelength_10,
               names_to = c(".value", "trials"),
               names_sep = "_" ) %>%
  # Average 2 dry mass measurement
  mutate(ADM = rowMeans(cbind(Aboveground_dry_mass1,
                              Aboveground_dry_mass2), na.rm = TRUE),
         BDM = rowMeans(cbind(Belowground_dry_mass1,
                              Belowground_dry_mass2), na.rm = TRUE),
         .keep = "unused") %>%
  # Rename trait abbreviation
  # Reformat ADM and BDM to per shoot to overcome sampling bias at some location (i.e. no focal samples only contain less than 10 shoots)
  mutate(ADM_S = ADM/Shoot_number,
         BDM_S = BDM/Shoot_number, 
         CanopyH = Canopyheight,
         LLength = Leaflength,
         LWidth = Leafwidth,
         RLength = Rootlength,
         Rdia = Rhizomediameter,
         ILength = Internodelength, .keep = "unused") %>%
  # Filter trait NA entry
  filter_at(vars(-Gravel, -Sand, -Silt, -Clay, -sed_mean), all_vars(!is.na(.))) %>%
  # arrange column
  dplyr::select(SampleID:abb, ADM_S, BDM_S,
         CanopyH:ILength, trials,
         total_REI:perci) -> trait_env_long
head(trait_env_long)
colnames(trait_env_long)

# Create a table with averaging repeated measurement variables (e.g. Leaf width, root length...)
trait_dat %>%
  # match trait data with environmental variables based on sample location 
  left_join(env_var, by = join_by("SampleID" == "location")) %>%
  # remove descriptive variables
  dplyr::select(!ends_with("photo.") & !contains("note") &
           !Pack & 
           !Aboveground_dry_mass3 &
           !Belowground_dry_mass3 &
           !Flower_seed) %>%
  # Average 2 dry mass measurements, and other repeated measurements
  mutate(ADM = rowMeans(cbind(Aboveground_dry_mass1,
                              Aboveground_dry_mass2), na.rm = TRUE),
         BDM = rowMeans(cbind(Belowground_dry_mass1,
                              Belowground_dry_mass2), na.rm = TRUE),
         CanopyH = rowMeans(cbind(.[grepl("Canopy", colnames(.))]), na.rm = TRUE),
         LLength = rowMeans(cbind(.[grepl("Leafle", colnames(.))]), na.rm = TRUE),
         LWidth = rowMeans(cbind(.[grepl("Leafwi", colnames(.))]), na.rm = TRUE),
         RLength = rowMeans(cbind(.[grepl("Root", colnames(.))]), na.rm = TRUE),
         Rdia = rowMeans(cbind(.[grepl("Rhizome", colnames(.))]), na.rm = TRUE),
         ILength = rowMeans(cbind(.[grepl("Internode", colnames(.))]), na.rm = TRUE),
         .keep = "unused") %>%
  # Reformat ADM and BDM to per shoot to overcome sampling bias at some location (i.e. no focal samples only contain less than 10 shoots)
  mutate(ADM_S = ADM/Shoot_number,
         BDM_S = BDM/Shoot_number, .keep = "unused") %>%
  # arrange column
  dplyr::select(SampleID:abb, ADM_S, BDM_S, CanopyH:ILength,
         total_REI:perci) -> trait_env_ave
head(trait_env_ave)
colnames(trait_env_ave)

#' # Halodule uninervis

#+ Create HU data
#+ #### HU ####
# Long table
trait_env_long[trait_env_long$abb=="HUN"|trait_env_long$abb=="HUW", ] -> HU_long
# Average table
trait_env_ave[trait_env_ave$abb=="HUN"|trait_env_ave$abb=="HUW", ] -> HU_ave


#' ## HU two peaks distribution and coefficient of variation
#+ HU two peaks and CV analysis
#+ #### HU two peak plots and CV analysis ####

# Histogram
LW_plot = ggplot(data = HU_long) +
  geom_histogram(aes(x = LWidth, fill = abb),
                 alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_classic() +
  labs(x = "Leaf width (mm)", fill = "Growth form")
  
RD_plot = ggplot(data = HU_long) +
  geom_histogram(aes(x = Rdia, fill = abb),
                 alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  theme_classic() +
  labs(x = "Rhizome diameter (mm)", fill = "Growth form")

ggarrange(LW_plot, RD_plot, legend = "right", common.legend = TRUE) %>% 
  ggexport(filename = "../plots/LW_RD_hist.png",
           width = 4000, height = 1600, res= 500)

#+ HU PCA
#### HU PCA ####
# Explanation material: https://personal.utdallas.edu/~herve/abdi-awPCA2010.pdf

# Create data for environmental covariates pva
HU_ave %>%
  dplyr::select(SampleID, abb, total_REI:ave_temp, sed_mean:perci) -> HU_env

# Environmental covariates Bayesian pca with missing value
HU_env_pca <- pcaMethods::pca(HU_env[, -c(1, 2)], method = "bpca", scale = "uv", center = TRUE, nPcs = 2)

summary(HU_env_pca)

# loading for each covariates
loadings <- as.data.frame(loadings(HU_env_pca)) %>%
  mutate(x0 = rep(0, nrow(.)),
         y0 = rep(0, nrow(.)))
# location for each sample
scores <- as.data.frame(pcaMethods::scores(HU_env_pca)) %>% 
  mutate(abb = HU_ave$abb)
## Get the estimated complete observations
cObs <- completeObs(HU_env_pca)

# Create a env_bpca plot
ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2, fill = abb), size = 2.5,
             shape = 21, stroke = 0.5) +
  geom_segment(data = loadings, aes(x = x0, y = y0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, 
                                       label = c("Total REI",
                                                 "Depth",
                                                 "Air exposure",
                                                 "Average Temp.",
                                                 "Mean Sed. size",
                                                 "Runoff",
                                                 "Precipitation")),
                  color = "black", size = 3.5,
                  position = position_nudge_center(0.2, 0.1, 0, 0)) +
  scale_fill_manual(values = c(HUN = "#009E73", HUW = "#E69F00")) + 
  labs(x = "Dim1 (33.3%)", y = "Dim2 (24.6%)", 
       title = "Enivornmental bpca",
       fill = "Clusters") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(1, "lines"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14),
    panel.grid.minor = element_blank()
  ) -> p
ggexport(p, filename = "../plots/Env_bpca.png",
         width = 3000, height = 3000, res = 600)


# Create data for trait pca analysis
HU_ave %>%
  dplyr::select(SampleID, abb, ADM_S:ILength) %>%
  na.omit() -> HU_ave1

# PCA analysis with scale and centered observations (samples) of each variables/columns (traits)
HU_ave_pca <- prcomp(HU_ave1[, -c(1, 2)], 
                 center =  TRUE,
                 scale = TRUE)
summary(HU_ave_pca)

# Plots to show contribution of variables (traits) and characteristic of components

# Store all pca explanation terms in var
# Show the value below accompany with plots
var <- get_pca_var(HU_ave_pca)

# Percentage of explained variance by components
# The value can be found from above, summary(HU_ave_pca): proportion of Variance
fviz_eig(HU_ave_pca, addlabels = TRUE)

# Correlation between variable and component or contribution of variables to component
# This correlation value determine the information shared by variables and components.
# In this way, this value determine the postion of variables in the biplot
# Smaller values in the first 2 components means the variables is less important for the first 2 components

# Real value
var$coord

# Position of the variables on a biplot
fviz_pca_var(HU_ave_pca, col.var = "black")

# Variables contribution to PC1, PC2 (another way of visualization)
p1 <- fviz_contrib(HU_ave_pca, choice = "var", axes = 1)
p2 <- fviz_contrib(HU_ave_pca, choice = "var", axes = 2)
p3 <- fviz_contrib(HU_ave_pca, choice = "var", axes = c(1, 2))
grid.arrange(p1, p2, p3, nrow = 1)


# Quality of representation of components to variables or observation
# This is denote by Cos2 (cosine square)
# Cos2 indicates how well does the component explain variables or observation
# In general the value of variable Cos2 has similar order to correlation between components and variables

# value of component representation on each variable
var$cos2

# Plot the 1st and 2nd component representation on each variable
fviz_cos2(HU_ave_pca, choice = "var", axes = c(1,2))

# Plot the 1st and 2nd component representation on each observation
fviz_cos2(HU_ave_pca, choice = "ind", axes = c(1,2))


# Plot the pca biplot, use customize function based on fviz_pca_biplot
# Save the PCA plot colored by growth form.
ITV_pca_biplot(HU_ave_pca, HU_ave$abb, 2) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11)) -> p
ggexport(p, filename = "../plots/HU_pca_plot.png",
         width = 3000, height = 3000, res = 600)

#' ## HU GLMM
#' Work on GLMM. Perform following steps:\
#' 1. Check environmental variables correlation
#' 2. Identify family\
#' 3. Run model\
#' 4. Compared model prediction with original data\

#+ HU data exploration
#### HU Mixture model ####

# create a data for the Mixture model and GLM
# rename and scale the environmental variables
# Air exposure are levels (0: non expose, 1:10-20% exposure time ... , 9: 80-100% exposure time)
# The levels are actually ordered and continous. For now, treat them as normal variable first

HU_ave %>%
  mutate(REI = scale(total_REI)[, 1],
         Depth =  scale(depth)[, 1],
         Ave_temp = scale(ave_temp)[, 1],
         Max_temp = scale(mean_max)[, 1],
         Min_temp = scale(mean_min)[, 1],
         Var_temp = scale(mean_var)[, 1],
         Runoff = scale(runoff)[, 1],
         Perci = scale(perci)[, 1],
         Air = scale(air_exposure)[, 1],
         Sed = (sed_mean - mean(sed_mean, na.rm = TRUE))/ sd(sed_mean, na.rm = TRUE), .keep = "unused") %>%
         mutate(ID = match(Location, unique(Location))) -> HU_GLM

# Visualize variables distribution
HU_GLM %>%
  dplyr::select(ADM_S:ILength) -> HU_box

HU_GLM %>%
  dplyr::select(REI:Sed) -> HU_box_env

ggplot(data = pivot_longer(HU_box, everything())) +
  geom_boxplot(aes(x = name, y = value))
# The range of traits are quite different, Below two plots exclude the traits that are too tiny to visualize in this plot

HU_box %>%
  dplyr::select(-c(CanopyH, LLength, RLength)) %>%
  pivot_longer(everything()) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value))

HU_box %>%
  dplyr::select(LWidth, Rdia) %>%
  pivot_longer(everything()) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value))


ggplot(data = pivot_longer(HU_box_env, everything())) +
  geom_boxplot(aes(x = name, y = value))
# Air, Depth, REI, Runoff, and Sed are the variables have higher resolution
# Warning for 10 NA in the Sed

pairs.panels(HU_box,
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)

# Check explained variable (x1...xn) correlation
pairs.panels(HU_box_env,
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)


# After processes above, I identifies that a few explained variables need to be log before scale, and all will need to be scaled except categorical variables or proportional variables\
# total_REI: skew distribution -> scale(log(total_REI))

# Replace scale(REI) with scale(log(REI))
hist(HU_GLM$REI)
HU_GLM$REI = scale(log(HU_ave$total_REI))[, 1]
hist(HU_GLM$REI)

#+ Multivariate model
#### Multivariate mixture and Multivariate ####

# After a few combination test. Only Leaf width and Rhizome diameter can be identify by model as two mixture component.
# Below model use these two traits and identify the best mixture proportion for both of them.

# For accuracy, I use raw Leaf width and Rhizome diameter instead of sample average.
# There is no concern about the resolution of the environmental factors in this model
# Provide more data would obtain better model result.

dat_list <- list(
  N = nrow(HU_long),
  K = 2,
  NY = ncol(HU_long[, c(13, 15)]),
  Y = as.matrix(HU_long[, c(13, 15)]) # 13: LWidth, 15: Rdia
)
MUL_mix_mod <- cmdstan_model("../stan_file/Multivariate_mixture.stan")

# Run MCMC using the 'sample' method from cmdstan
MUL_mix <- MUL_mix_mod$sample(
  data = dat_list,
  seed = 568,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 6000,
  iter_sampling = 2000,
  adapt_delta = 0.97
)

# This stan code can be easily change to Mixture multivariate lognormal regression by changing log(y[N]) in the likelihood section and changing positive_ordered to ordered (log(mean) can be negative).
# This is helpful to keep simulation become all positive.
# However, this will lead to the estimation become a log normal distribution which is generally skew.

# saveRDS(MUL_mix, file = "../models/Multivariate_mixture.Rds")
MUL_mix = readRDS(file = "./models/Multivariate_mixture.Rds")

# Check divergence and tree depth
MUL_mix$diagnostic_summary()
# trace plot
MUL_mix$draws() %>% mcmc_trace()
# precis from rethinking is still wotking to get summary table
precis(MUL_mix, depth = 3)

precis(MUL_mix, depth = 3) %>%
  as.data.frame() -> mix_coef

# Here overlap observed data and model prediction as density plot
# Leaf width

sigmas_k1 <- c(mix_coef[15,1], mix_coef[17,1]) # standard deviations
Rho_k1 <- matrix(c(1, mix_coef[9,1], mix_coef[11,1], 1) , nrow=2) # correlation matrix
Sigma_k1 <- diag(sigmas_k1) %*% Rho_k1 %*% diag(sigmas_k1)

sigmas_k2 <- c(mix_coef[16,1], mix_coef[18,1]) # standard deviations
Rho_k2 <- matrix( c(1, mix_coef[10,1], mix_coef[12,1], 1) , nrow=2 ) # correlation matrix
Sigma_k2 <- diag(sigmas_k2) %*% Rho_k2 %*% diag(sigmas_k2)

data.frame(NLW = rmvnorm(10000, c(mix_coef[3, 1], mix_coef[4, 1]),
                         Sigma_k1)[, 1],
           WLN = rmvnorm(10000, c(mix_coef[5, 1], mix_coef[6, 1]),
                         Sigma_k2)[, 1]) %>%
  ggplot() +
  geom_density(aes(x = NLW), fill = "grey", color = "white") +
  geom_density(aes(x = WLN),fill = "grey", color = "white") +
  geom_density(data = HU_long, 
               aes(x = LWidth, fill = abb, col = abb), 
               alpha = 0.7) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#E69F00"), guide = "none") +
  theme_classic() +
  labs(x = "Leaf Width (mm)", y = "Density", fill = "Growth form") +
  xlim(0, 5.2) -> LW_est
LW_est
ggsave(LW_est, filename = "../plots/LW_estimated_density.png",
       width = 20, height = 15, unit = "cm")

# Same for Rhizome diameter
data.frame(NLW = rmvnorm(10000, c(mix_coef[3, 1], mix_coef[4, 1]),
                         Sigma_k1)[, 2],
           WLN = rmvnorm(10000, c(mix_coef[5, 1], mix_coef[6, 1]),
                         Sigma_k2)[, 2])%>%
  ggplot() +
  geom_density(aes(x = NLW), fill = "grey", color = "white") +
  geom_density(aes(x = WLN),fill = "grey", color = "white") +
  geom_density(data = HU_long, 
               aes(x = Rdia, fill = abb, col = abb), 
               alpha = 0.7) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#E69F00"), guide = "none") +
  theme_classic() +
  labs(x = "Rhizome Diameter (mm)", y = "Density", fill = "Growth form") +
  xlim(0, 3) -> RD_est
RD_est
ggsave(RD_est, filename = "../plots/RD_estimated_density.png",
       width = 20, height = 15, unit = "cm")


ggarrange(LW_est, RD_est, common.legend = TRUE, legend = "right") %>%
  ggexport(filename = "../plots/estimated_density.png",
           width = 6000, height = 3000, res = 800)


# Create plot to compare estimated mean
HU_long %>%
  filter(abb == "HUN") %>%
  dplyr::select(LWidth, Rdia) -> HUN_plot

apply(HUN_plot, 2, mean) -> aN

HU_long %>%
  filter(abb == "HUW") %>%
  dplyr::select(LWidth, Rdia) -> HUW_plot

apply(HUW_plot, 2, mean) -> aW

data.frame(Name = fct_relevel(factor(c("Leaf Width Mean",
                           "Rhizome diameter Mean")), 
                           "Rhizome diameter Mean",
                           "Leaf Width Mean"),
           HMean = c(aN[1], aN[2]),
           HMMean = c(mix_coef[3, 1], mix_coef[4, 1]),
           HMsd = c(mix_coef[3, 2], mix_coef[4, 2]),
           WMean = c(aW[1], aW[2]),
           WMMean = c(mix_coef[5, 1], mix_coef[6, 1]),
           WMsd = c(mix_coef[5, 2], mix_coef[6, 2])) %>%
  ggplot() +
  geom_point(aes(x = HMean, y = Name),  color = "#009E73", size = 4) +
  #geom_segment(aes(x = HMMean - HMsd, xend = HMMean + HMsd, 
  #                 y = Name, yend = Name), linewidth = 2, 
  #              color = "gray30") +
  geom_segment(aes(x = HMMean - 2*HMsd, xend = HMMean + 2*HMsd, 
                   y = Name, yend = Name), linewidth = 1, 
               color = "gray30") +
  geom_point(aes(x = HMMean, y = Name), color = "gray30", size = 3) +
    geom_point(aes(x = WMean, y = Name),  color = "#E69F00", size = 4) +  
  # geom_segment(aes(x = WMMean - WMsd, xend = WMMean + WMsd, y
  #                  = Name, yend = Name), linewidth = 2, 
  #              color = "gray30") +
  geom_segment(aes(x = WMMean - 2*WMsd, xend = WMMean + 2*WMsd, 
                   y = Name, yend = Name), linewidth = 1, 
               color = "gray30") +
  geom_point(aes(x = WMMean, y = Name), color = "gray30", size = 3) +
  labs(x = "Length (mm)", y = "") +
  theme_classic() -> HUcom
HUcom
ggsave(HUcom, filename = "../plots/Model_estimation.png",
       width = 20, height = 8, unit = "cm")


#' Separate data frame to two group and run analysis for environmental factors
#' The multivariate mixture model shows the covariance of 2 traits is low. Below, I separate model according to each traits and growth form.
#' Although, it's possible to put all of them together but I think it's not reasonable approach. Reasons are below.
#' 1. The traits covariance are low means that they don't affect others reaction to environmental factors. It's unnecessary to use multivariate respond variables model.
#' 2. Although with multilevel model, random slope for each growth form is a potential approach. This model will assume the traits of each growth are from the same distribution, which is contradicted to previous mixture model results.
#' Therefore, I decided to separate model to each growth form and each traits as below section.

#+ Two growth form GLM
#+ #### Two growth form GLM ####

# Create one data set for each growth form
HU_GLM %>%
  filter(abb == "HUN") -> HUN

pairs.panels(HUN[, 21:30],
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)
# Use r = 0.80 as threshold, remove highly correlated explained variables
# Remove 1. Max_temp, 2. Min temp, 3. Var_temp


HU_GLM %>%
  filter(abb == "HUW") -> HUW
pairs.panels(HUW[, 21:30],
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)
# Due to unbalance sample location, only REI, Depth, Air exposure, and Sediment varied across samples.

# First set is for Leaf width

# HUN leaf width explained by all not correlated variables

dat <- list(
  LW = HUN$LWidth,
  Rei = HUN$REI,
  De = HUN$Depth,
  AveT = HUN$Ave_temp,
  Run = HUN$Runoff,
  Per = HUN$Perci,
  Air = HUN$Air,
  Sed = HUN$Sed,
  ID = match(HUN$ID, unique(HUN$ID)))

m.nlw1 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ dnorm( 0 , 1 ),
    c(abar, bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_a, sigma_S) ~ normal(0, 0.5),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nlw1)
precis(m.nlw1, depth = 2)

# HUN leaf width explained by higher resolution variables
dat <- list(
  LW = HUN$LWidth,
  Rei = HUN$REI,
  De = HUN$Depth,
  Run = HUN$Runoff,
  Air = HUN$Air,
  Sed = HUN$Sed,
  ID = match(HUN$ID, unique(HUN$ID)))

m.nlw2 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar, bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),
  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nlw2)
precis(m.nlw2, depth = 2)

rethinking::WAIC(m.nlw1)
rethinking::WAIC(m.nlw2)
rethinking::compare(m.nlw1, m.nlw2)
# WAIC are same for both models (less than 0.01 difference)
# Choose the simpler model
# m.nlw2 (Model 2) is preferred 

# Residual check for m.nlw1
m.nlw2.post <- extract.samples(m.nlw2)
colMeans(m.nlw2.post$mu) -> m.nlw2.pred
m.nlw2.pred - HUN$LWidth -> m.nlw2.resd
plot(m.nlw2.pred, m.nlw2.resd)
qqplot(m.nlw2.pred, HUN$LWidth)
abline(0, 1, col = 2)

# HUW leaf width explained by all not correlated variables
dat <- list(
  LW = HUW$LWidth,
  Rei = HUW$REI,
  De = HUW$Depth,
  Air = HUW$Air,
  Sed = HUW$Sed,
  ID = match(HUW$ID, unique(HUW$ID)))

m.wlw1 <- ulam(
  alist(
    LW ~ normal( mu, sigma),
    mu <-  abar + z[ID]*sigma_a +  bS * Sed + bRei * Rei  + bD * De + bAir * Air,
    Sed ~ dnorm(nu, sigma_s),
    z[ID] ~ dnorm(0, 1),
    c(abar, bS, bRei, bD, bAir, nu) ~ normal(0, 0.5),
    c(sigma, sigma_a, sigma_s) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99), log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)

dashboard(m.wlw1)
precis(m.wlw1, depth = 2)


# HUW leaf width explained by higher resolution variables
dat <- list(
  LW = HUW$LWidth,
  Rei = HUW$REI,
  De = HUW$Depth,
  Sed = HUW$Sed,
  ID = match(HUW$ID, unique(HUW$ID)))

m.wlw2 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bD * De + bRei * Rei + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(abar, sigma_a),
    c(abar, bRei, bD, bS, nu) ~ normal(0, 0.5),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wlw2)
precis(m.wlw2, depth = 2)


rethinking::WAIC(m.wlw1)
rethinking::WAIC(m.wlw2)
rethinking::compare(m.wlw1, m.wlw2)
# m.wlw2 (Model 2) is preferred 

m.wlw2.post <- extract.samples(m.wlw2)
colMeans(m.wlw2.post$mu) -> m.wlw2.pred
m.wlw2.pred - HUW$LWidth -> m.wlw2.resd
plot(m.wlw2.pred, m.wlw2.resd)
qqplot(m.wlw2.pred, HUW$LWidth)

# extrapolated Leaf width based on model

# Double check the structure of selected model
# m.nlw2: abar + z[ID] * sigma_a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed
# m.wlw2: abar + z[ID] * sigma_a + bD * De + bRei * Rei + bS * Sed

#  Control all variables except mean sediment size
S_seq = seq(-1, 2.5, length.out = 50)
S_seq * sd(HU_ave$sed_mean,na.rm = TRUE) + mean(HU_ave$sed_mean, na.rm = TRUE) -> S_org

post <- extract.samples(m.nlw2)
p_link <- function(seq , ID) {
  mm <- with( post ,
                   a[, ID] + bRei * mean(HUN$REI) + bD * mean(HUN$Depth) + bRun * mean(HUN$Runoff) + bAir * mean(HUN$Air) + bS * S_seq[seq])
  # Randomly sample a row index from the post list
  row_index <- sample(1:8000, 1000)
  
  pred <- with(post, sapply(row_index, function(i) rnorm(1, mm[i], sigma[i])))
  return(pred)
}

N_lw_sim <- lapply(1:ncol(post$a), function(j) sapply(1:50, function(i) p_link(i, ID = j)))

N_lw <- do.call(rbind, N_lw_sim[1:ncol(post$a)])


post <- extract.samples(m.wlw2)
p_link <- function(seq , ID) {
  mm <- with( post ,
              a[, ID] + bRei * mean(HUW$REI) + bD * mean(HUW$Depth) + bS * S_seq[seq])
  # Randomly sample a row index from the post list
  row_index <- sample(1:8000, 1000)
  
  pred <- with(post, sapply(row_index, function(i) rnorm(1, mm[i], sigma[i])))
  return(pred)
}

W_lw_sim <- lapply(1:ncol(post$a), function(j) sapply(1:50, function(i) p_link(i, ID = j)))

W_lw <- do.call(rbind, W_lw_sim[1:ncol(post$a)])

data.frame(LW = c(colMeans(N_lw), colMeans(W_lw)),
           Low = c(apply(N_lw, 2, PI)[1, ], 
                   apply(W_lw, 2, PI)[1, ]),
           High = c(apply(N_lw, 2, PI)[2, ],
                    apply(W_lw, 2, PI)[2, ]),
           Sed = rep(S_org, times = 2),
           Abb = c(rep("HUN", 50), rep("HUW", 50))) %>%
  ggplot() +
  geom_line(aes(x = Sed, y = LW, color = Abb), linewidth = 2) +
  geom_ribbon(aes(x = Sed, ymin = Low, ymax = High, color = Abb), , linetype = 2, alpha = 0.1) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  labs(y = "Simulated leaf width (mm)",
       x = expression(paste("Mean sediment size (", mu, "m)")),
       color = "Growth form") + 
  theme_classic() -> LW_sed
LW_sed

ggsave(LW_sed, filename = "../plots/LW_sed_effect.png",
       width = 20, height = 15, unit = "cm")

#Second set is for Rhizome diameter

# HUN Rhizome diameter explained by all not correlated variables
dat <- list(
  Rd = HUN$Rdia,
  Rei = HUN$REI,
  De = HUN$Depth,
  AveT = HUN$Ave_temp,
  Run = HUN$Runoff,
  Per = HUN$Perci,
  Air = HUN$Air,
  Sed = HUN$Sed,
  ID = match(HUN$ID, unique(HUN$ID)))

m.nrd1 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ dnorm( 0 , 1 ),
    c(abar, bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_a, sigma_S) ~ normal(0, 0.5),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nrd1)
precis(m.nrd1, depth = 2)

# HUN Rhizome diameter explained by higher resolution variables
dat <- list(
  Rd = HUN$Rdia,
  Rei = HUN$REI,
  De = HUN$Depth,
  Run = HUN$Runoff,
  Air = HUN$Air,
  Sed = HUN$Sed,
  ID = match(HUN$ID, unique(HUN$ID)))

m.nrd2 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ dnorm( 0 , 1 ),
    c(abar, bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nrd2)
precis(m.nrd2, depth = 2)

rethinking::WAIC(m.nrd1)
rethinking::WAIC(m.nrd2)
rethinking::compare(m.nrd1, m.nrd2)
# m.nrd2 (Model 2) is preferred. The difference is tiny, but m.nrd2 is still the better and simpler model

m.nrd2.post <- extract.samples(m.nrd2)
colMeans(m.nrd2.post$mu) -> m.nrd2.pred
m.nrd2.pred - HUN$Rdia -> m.nrd2.resd
plot(m.nrd2.pred, m.nrd2.resd)
qqplot(m.nrd2.pred, HUN$Rdia)

# HUW Rhizome diameter explained by all not correlated variables
dat <- list(
  Rd = HUW$Rdia,
  Rei = HUW$REI,
  De = HUW$Depth,
  Air = HUW$Air,
  Sed = HUW$Sed,
  ID = match(HUW$ID, unique(HUW$ID)))

m.wrd1 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar, bRei, bD, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wrd1)
precis(m.wrd1, depth = 2)

# HUW Rhizome diameter explained by higher resolution variables
dat <- list(
  Rd = HUW$Rdia,
  Rei = HUW$REI,
  De = HUW$Depth,
  Sed = HUW$Sed, 
  ID = match(HUW$ID, unique(HUW$ID)))

m.wrd2 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar, bRei, bD, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wrd2)
precis(m.wrd2, depth = 2)

rethinking::WAIC(m.wrd1)
rethinking::WAIC(m.wrd2)
rethinking::compare(m.wrd1, m.wrd2)
# There is 0.9 difference in WAIC, chose the simpler model
# m.wrd2 (Model 2) is preferred 

m.wrd2.post <- extract.samples(m.wrd2)
colMeans(m.wrd2.post$mu) -> m.wrd2.pred
m.wrd2.pred - HUW$Rdia -> m.wrd2.resd
plot(m.wrd2.pred, m.wrd2.resd)
qqplot(m.wrd2.pred, HUW$Rdia)


# extrapolated Rhizome diameter based on model

# Double check the structure of selected model
# m.nrd2: abar + z[ID] * sigma_a + +bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed
# m.wrd2: abar + z[ID] * sigma_a + bRei * Rei + bD * De + bS * Sed

#  Control all variables except mean sediment size
S_seq = seq(-1, 2.5, length.out = 50)
S_seq * sd(HU_ave$sed_mean,na.rm = TRUE) + mean(HU_ave$sed_mean, na.rm = TRUE) -> S_org

post <- extract.samples(m.nrd2)
p_link <- function(seq , ID) {
  mm <- with( post ,
              a[, ID] + bRei * mean(HUN$REI) + bD * mean(HUN$Depth) + bRun * mean(HUN$Runoff) + bAir * mean(HUN$Air) + bS * S_seq[seq])
  # Randomly sample a row index from the post list
  row_index <- sample(1:8000, 1000)
  
  pred <- with(post, sapply(row_index, function(i) rnorm(1, mm[i], sigma[i])))
  return(pred)
}

N_rd_sim <- lapply(1:ncol(post$a), function(j) sapply(1:50, function(i) p_link(i, ID = j)))

N_rd <- do.call(rbind, N_rd_sim[1:ncol(post$a)])

post <- extract.samples(m.wrd2)
p_link <- function(seq , ID) {
  mm <- with( post ,
              a[, ID] + bRei * mean(HUW$REI) + bD * mean(HUW$Depth) + bS * S_seq[seq])
  # Randomly sample a row index from the post list
  row_index <- sample(1:8000, 1000)
  
  pred <- with(post, sapply(row_index, function(i) rnorm(1, mm[i], sigma[i])))
  return(pred)
}

W_rd_sim <- lapply(1:ncol(post$a), function(j) sapply(1:50, function(i) p_link(i, ID = j)))

W_rd <- do.call(rbind, W_rd_sim[1:ncol(post$a)])

data.frame(Rd = c(colMeans(N_rd),colMeans(W_rd)),
           Low = c(apply(N_rd, 2, PI)[1, ], 
                   apply(W_rd, 2, PI)[1, ]),
           High = c(apply(N_rd, 2, PI)[2, ],
                    apply(W_rd, 2, PI)[2, ]),
           Sed = rep(S_org, times = 2),
           Abb = c(rep("HUN", 50), rep("HUW", 50))) %>%
  ggplot() +
  geom_line(aes(x = Sed, y = Rd, color = Abb), linewidth = 2) +
  geom_ribbon(aes(x = Sed, ymin = Low, ymax = High, color = Abb), , linetype = 2, alpha = 0.1) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  labs(y = "Simulated Rhizome diameter (mm)",
       x = expression(paste("Mean sediment size (", mu, "m)")),
       color = "Growth form") + 
  theme_classic() -> RD_sed
RD_sed

ggsave(RD_sed, filename = "../plots/RD_sed.png",
       width = 20, height = 15, unit = "cm")

ggarrange(LW_sed, RD_sed, common.legend = TRUE, legend = "right") %>%
  ggexport(filename = "../plots/Sed_effect.png",
           width = 6000, height = 3000, res = 800)

#' ADM/S and BDM/S GLM regardless growth form
#### ADM/S, BDM/S GLM ####

# Firstly, scan through variable correlation
pairs.panels(HU_GLM[, 21:30],
            gap = 0,
            density = TRUE,  # show density plots
            hist.col = "#00AFBB",
            ellipses = TRUE) # show correlation ellipses)
# Exclude highly correlated (0.9) variables: 1. Max_temp, 2. Min_temp, 3. Var temp
# Only use higher resolution variables

# HU Above ground biomass per shoot explained by all not correlated variables
dat <- list(
  Adm = log(HU_GLM$ADM_S),
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  AveT = HU_GLM$Ave_temp,
  Run = HU_GLM$Runoff,
  Per = HU_GLM$Perci,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed,
  ID = HU_GLM$ID)

m.adm1 <- ulam(
  alist(
    Adm ~ normal(mu, sigma),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar,bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.adm1)
precis(m.adm1, depth = 2)

# HU Above ground biomass per shoot explained by higher resolution variables
dat <- list(
  Adm = log(HU_GLM$ADM_S),
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  Run = HU_GLM$Runoff,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed,
  ID = HU_GLM$ID)

m.adm2 <- ulam(
  alist(
    Adm ~ normal(mu, sigma),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar, bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.adm2)
precis(m.adm2, depth = 2)

rethinking::WAIC(m.adm1)
rethinking::WAIC(m.adm2)
rethinking::compare(m.adm1, m.adm2)
# There is only 0.2 difference in WAIC value
# Choose the simpler model, m.adm2 (Model 2)

m.adm2.post <- extract.samples(m.adm2)
colMeans(m.adm2.post$mu) -> m.adm2.pred
m.adm2.pred - log(HU_GLM$ADM_S) -> m.adm2.resd
plot(m.adm2.pred, m.adm2.resd)
qqplot(m.adm2.pred, log(HU_GLM$ADM_S))
abline(a = 0, b = 1)

# HU Below ground biomass per shoot explained by all not correlated variables
dat <- list(
  Bdm = log(HU_GLM$BDM_S),
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  AveT = HU_GLM$Ave_temp,
  Run = HU_GLM$Runoff,
  Per = HU_GLM$Perci,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed,
  ID = HU_GLM$ID)

m.bdm1 <- ulam(
  alist(
    Bdm ~ normal( mu, sigma ),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar,bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.bdm1)
precis(m.bdm1, depth = 2)

# HU Below ground biomass per shoot explained by higher resolution variables
dat <- list(
  Bdm = log(HU_GLM$BDM_S),
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  Run = HU_GLM$Runoff,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed,
  ID = HU_GLM$ID)

m.bdm2 <- ulam(
  alist(
    Bdm ~ normal(mu, sigma),
    mu <- abar + z[ID]*sigma_a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    z[ID] ~ normal(0, 1),
    c(abar, bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S, sigma_a) ~ exponential(1),
    gq> vector[ID]:a <<- abar + z*sigma_a
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.bdm2)
precis(m.bdm2, depth = 2)

rethinking::WAIC(m.bdm1)
rethinking::WAIC(m.bdm2)
rethinking::compare(m.bdm1, m.bdm2)
# The difference of WAIC is only 0.1, chose less complex model
# m.bdm2 (Model 2) is preferred 

m.bdm2.post <- extract.samples(m.bdm2)
colMeans(m.bdm2.post$mu) -> m.bdm2.pred
m.bdm2.pred - log(HU_GLM$BDM_S) -> m.bdm2.resd
plot(m.bdm2.pred, m.bdm2.resd)
qqplot(m.bdm2.pred, log(HU_GLM$BDM_S))
abline(a = 0, b = 1, col = "red")


#' Fish biomass index using fish PCA combination
#+ #### Fish biomass ####
# Followed the fish abundance and meadow structure analysis approach (Jones et.al, 2021, https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.640528/full).
# Put four traits in the PCA to see PC1 composition

HU_ave1 %>%
  dplyr::select(LLength, LWidth, CanopyH, ADM_S) -> HU_fish

HU_fish_pca <- prcomp(HU_fish, 
                      center =  TRUE,
                      scale = TRUE)
summary(HU_fish_pca)

ITV_pca_biplot(HU_fish_pca, HU_ave$abb, 2) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11)) -> p
ggexport(p, filename = "../plots/HU_fish_plot.png",
         width = 3000, height = 3000, res = 600)

get_pca_var(HU_fish_pca)$coord

# Use PC1 value as input
HU_ave %>%
  mutate(PC1 = HU_fish_pca[["x"]][, 1]) -> HU_glm_pca

HU_glm_pca$abb[HU_glm_pca$abb == "HUN"] = 0
HU_glm_pca$abb[HU_glm_pca$abb == "HUW"] = 1

dat <- list(
  PC1 = HU_glm_pca$PC1,
  Abb = as.numeric(HU_glm_pca$abb)
)

fit_PC1 <- ulam(
  alist(
    PC1 ~ normal(mu, sigma),
    mu <- a + bAbb * Abb,
    a ~ normal(1, 1),
    bAbb ~ normal(0, 1),
    sigma ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.95),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(fit_PC1)
precis(fit_PC1)


extract.samples(fit_PC1)-> fit_PC1_draws

PC1_sim <- as.data.frame(matrix(nrow = 1000, ncol = 5))
for(i in 1:1000){
  PC1_sim[i, 1] = rnorm(1, mean = fit_PC1_draws$a[i], sd = fit_PC1_draws$sigma[i])
  PC1_sim[i, 5] = rnorm(1, mean = fit_PC1_draws$a[i] + fit_PC1_draws$bAbb[i], sd = fit_PC1_draws$sigma[i])
}

PC1_sim[, 2] = c(rep(1, 300), rep(2, 700))
PC1_sim[, 3] = c(rep(1, 500), rep(2, 500))
PC1_sim[, 4] = c(rep(1, 700), rep(2, 300))

for(i in 1:1000){
  for(j in 2:4){
    if(PC1_sim[i, j] == 2){
      PC1_sim[i, j] = rnorm(1, mean = fit_PC1_draws$a[i], sd = fit_PC1_draws$sigma[i])
      } else {
        PC1_sim[i, j] = rnorm(1, mean = fit_PC1_draws$a[i] + fit_PC1_draws$bAbb[i], sd = fit_PC1_draws$sigma[i])
        }
  }
}

colnames(PC1_sim) = c("100% HUN", "70% HUN + 30% HUW", "50% HUN + 50% HUW",
                      "30% HUN + 70% HUW", "100% HUW")
pivot_longer(PC1_sim, 1:5, values_to = "index", names_to = "type") %>%
  mutate(type = fct_reorder(type, index)) -> PC1_sim_plot

colorRampPalette(c("#009E73", "#E69F00")) -> colgra

ggplot(PC1_sim_plot) +
  geom_violin(aes(y = type, x = index, fill = type, color = type)) +
  scale_color_manual(values = colgra(5)) +
  scale_fill_manual(values = colgra(5)) + 
  theme_classic() +
  labs(x = "Fish biomass index",
       y = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12)) -> aa
aa
aa %>% ggexport(filename = "../plots/Biomassindex.png",
           width = 3000, height = 2500, res = 600)

#' Blue carbon represented by three traits
#### Blue carbon ####

HU_ave %>%
  dplyr::select(SampleID:abb,ADM_S, BDM_S, LWidth, CanopyH) %>%
  mutate(Biomass = ADM_S + BDM_S,
         CanopyC = CanopyH * LWidth, .keep = "unused") -> HU_car

HU_car$abb[HU_car$abb == "HUN"] = 0
HU_car$abb[HU_car$abb == "HUW"] = 1

dat <- list(
  Bio = log(HU_car$Biomass),
  Cc = log(HU_car$CanopyC),
  Abb = as.numeric(HU_car$abb)
)

fit_car <- ulam(
  alist(
    c(Bio, Cc) ~  multi_normal(c(muB, muC), Rho, sigma),
    muB <- aB + bB * Abb,
    muC <- aC + bC * Abb,
    c(aB, aC) ~ normal(0, 1),
    c(bB, bC) ~ normal(0, 1),
    Rho ~ lkj_corr(3),
    sigma ~ exponential( 1 )
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.95),
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(fit_car)
precis(fit_car, depth = 3)
precis(fit_car, depth = 3) -> mix_coef
extract.samples(fit_car)

sigmas <- c(mix_coef[9,1], mix_coef[10,1]) # standard deviations
Rho<- matrix(c(1, mix_coef[6,1], mix_coef[7,1], 1), nrow=2) # correlation matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

data.frame(N = rmvnorm(10000, c(mix_coef[1, 1],
                                mix_coef[2, 1]), Sigma),
           W = rmvnorm(10000, c(mix_coef[1, 1] + mix_coef[3, 1],
                                mix_coef[2, 1] + mix_coef[4, 1]), Sigma)) -> car_sim

pick_3 <- sample(1:10000, 3000, replace = FALSE)
pick_5 <- sample(1:10000, 5000, replace = FALSE)
pick_7 <- sample(1:10000, 7000, replace = FALSE)

colnames(car_sim) = c("NCc", "NBi", "WCc", "WBi")

car_meadow <- as.data.frame(matrix(nrow = 10000, ncol = 5))
car_meadow[, 5] = 0.3 * car_sim$WBi + 0.7 * car_sim$WCc 
car_meadow[, 4] = c(0.3 * car_sim$WBi[pick_7] + 
                      0.7 * car_sim$WCc[pick_7],
                    0.3 * car_sim$NBi[pick_3] + 
                      0.7 * car_sim$NCc[pick_3])
car_meadow[, 3] = c(0.3 * car_sim$WBi[pick_5] + 
                      0.7 * car_sim$WCc[pick_5],
                    0.3 * car_sim$NBi[pick_5] + 
                      0.7 * car_sim$NCc[pick_5])
car_meadow[, 2] = c(0.3 * car_sim$WBi[pick_3] + 
                      0.7 * car_sim$WCc[pick_3],
                    0.3 * car_sim$NBi[pick_7] + 
                      0.7 * car_sim$NCc[pick_7])
car_meadow[, 1] = 0.3 * car_sim$NBi + 0.7 * car_sim$NCc 

colnames(car_meadow) = c("100% HUN", "70% HUN + 30% HUW", 
                         "50% HUN + 50% HUW",
                         "30% HUN + 70% HUW", "100% HUW")
pivot_longer(car_meadow, 1:5, values_to = "BC", names_to = "type") %>%
  mutate(type = fct_reorder(type, BC)) -> car_meadow_plot

colorRampPalette(c("#009E73", "#E69F00")) -> colgra

ggplot(car_meadow_plot) +
  geom_violin(aes(y = type, x = BC, fill = type, color = type)) +
  scale_color_manual(values = colgra(5)) +
  scale_fill_manual(values = colgra(5)) + 
  theme_classic() +
  labs(x = "Carbon storage (log transformed)",
       y = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12)) -> p
p %>% ggexport(filename = "../plots/Bluecarbon.png",
                width = 3000, height = 2500, res = 600)