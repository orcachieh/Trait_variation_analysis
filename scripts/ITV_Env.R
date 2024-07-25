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
library(rstan)
library(rethinking)
library(lme4)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra) # arranging plots
library(psych) # pairwise scatter plot
library(DHARMa)     #for residual diagnostics
library(emmeans)
library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(MASS)
library(tidyverse)
library(mvtnorm)
library(tidybayes)
library(bayesplot)

#+ Source functions from Functions.R 
#### Input functions  ####

source("./Functions.R")


#' # Trait and Environmental Variables 

#+ read and process data
#### Read and match trait and env. variables ####

trait_dat <- read.csv("../trait_raw_data.csv", header = TRUE)
env <- read.csv("../environment_variables.csv", header = TRUE)

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
  select(!ends_with("photo.") & !contains("note") &
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
  select(SampleID:abb, ADM_S, BDM_S,
         CanopyH:ILength, trials,
         total_REI:perci) -> trait_env_long
head(trait_env_long)
colnames(trait_env_long)

# Create a table with averaging repeated measurement variables (e.g. Leaf width, root length...)
trait_dat %>%
  # match trait data with environmental variables based on sample location 
  left_join(env_var, by = join_by("SampleID" == "location")) %>%
  # remove descriptive variables
  select(!ends_with("photo.") & !contains("note") &
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
  select(SampleID:abb, ADM_S, BDM_S, CanopyH:ILength,
         total_REI:perci) -> trait_env_ave
head(trait_env_ave)
colnames(trait_env_ave)

#' # Halodule uninervise

#+ Create HU data
#+ 
#### HU ####
# Long table
trait_env_long[trait_env_long$abb=="HUN"|trait_env_long$abb=="HUW", ] -> HU_long
# Average table
trait_env_ave[trait_env_ave$abb=="HUN"|trait_env_ave$abb=="HUW", ] -> HU_ave


#' ## HU two peaks distribution and coefficient of variation
#+ HU two peaks and CV analysis
#### HU two peak plots and CV analysis ####

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

# CV
HU_long %>% select(11:16) %>% colMeans(., na.rm = TRUE) -> mean1
HU_long %>% select(11:16) %>% apply(2, sd, na.rm = TRUE) -> sd1

HU_ave %>% select(9:10) %>% colMeans(., na.rm = TRUE) -> mean2
HU_ave %>% select(9:10) %>% apply(2, sd, na.rm = TRUE) -> sd2

HU_CV <- t((cbind(rbind(mean1, sd1),
                  rbind(mean2, sd2)))) %>%
  as.data.frame() %>%
  mutate(CV = sd1/mean1, Trait = row.names(.))

CV_plot <- ggplot(HU_CV) +
  geom_point(aes(x = reorder(Trait, CV), y = CV)) +
  scale_x_discrete(labels =c("RD", "IL", "CH", "RL", "LL", "LW", "BDM/S", "ADM/S")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.3))
CV_plot
ggsave(filename = "../plots/CV_plot.png", plot = CV_plot, 
       width = 10, height = 6, units = "cm")

# Create HO CV plot
trait_env_long[trait_env_long$abb=="HO", ] -> HO_long
# Average table
trait_env_ave[trait_env_ave$abb=="HO", ] -> HO_ave


trait_env_long[trait_env_long$abb=="HO", ] -> HO_long
# Average table
trait_env_ave[trait_env_ave$abb=="HO", ] -> HO_ave

HO_long %>% select(11:16) %>% colMeans(., na.rm = TRUE) -> mean1
HO_long %>% select(11:16) %>% apply(2, sd, na.rm = TRUE) -> sd1

HO_ave %>% select(9:10) %>% colMeans(., na.rm = TRUE) -> mean2
HO_ave %>% select(9:10) %>% apply(2, sd, na.rm = TRUE) -> sd2

HO_CV <- t((cbind(rbind(mean1, sd1),
                  rbind(mean2, sd2)))) %>%
  as.data.frame() %>%
  mutate(CV = sd1/mean1, Trait = row.names(.))


HO_CV_plot <- ggplot(HO_CV) +
  geom_point(aes(x = reorder(Trait, CV), y = CV)) +
  scale_x_discrete(labels =c("LW", "LL", "RH", "CH", "IL", "RL", "BDM/S", "ADM/S")) +
  theme_classic() +
  theme(axis.title.x = element_blank())
HO_CV_plot
ggsave(filename = "../plots/HO_CV_plot.png", plot = CV_plot, 
       width = 10, height = 6, units = "cm")


#' ## HU Clustering

#+ HU Clustering analysis
#### HU Clustering####

# Use average repeated measurements data to create a cluster data for HU and remove NA
HU_ave %>%
  select(SampleID, abb, ADM_S:ILength) %>%
  na.omit() -> HU_clu_data

# Calculate distance matrix
# 1. Standardize (Normalize) each column (x-mean)/sd
# 2. Calculate dissimilarity matrix based on Euclidean distance

# remove first two columns (ID columns) then scale
HU_clu_data_std<- scale(HU_clu_data[, -c(1, 2)],center=TRUE,scale=TRUE)
head(HU_clu_data_std)

# Calculate dissimilarity distance matrix and visualize it
distance = get_dist(HU_clu_data_std, method = "euclidean")
fviz_dist(distance, 
          gradient = list(low = "#00AFBB", 
                          mid = "white", 
                          high = "#FC4E07"))

#' Hierarchical clustering analysis

#+ H-clustering
#### HU hierarchical clustering ####

# hclust: every sample as single cluster (leaves) and then join two closest clusters at each step
HU_hclu_comp = hclust(distance, "complete")
HU_hclu_WD2 = hclust(distance, "ward.D2")

plot(HU_hclu_comp, hang = -1,
     labels=HU_clu_data$abb,main="method = complete")
plot(HU_hclu_WD2, hang = -1,
     labels=HU_clu_data$abb,main="method = Ward.D2")

#' Partitioning Clustering\
#' K-means clustering is used here

#+ K-means
#### HU K-means ####
set.seed(123)


# 1. Scree plot, sum the wss of each number of clusters
# wss: total within-cluster sum of square. K-means clustering minimize wss in a reasonable number of clusters.
# general people consider bent part as optimal cluster

fviz_nbclust(HU_clu_data_std, kmeans, method = "wss")
# To me the bend seems to be happened at 4 clusters

# 2. Average Silhouette Method
# This method measures the quality of a clustering by how well the object lies within a cluster, indicates by silhouette width.
# High average Silhouette width indicate good cluster number 

fviz_nbclust(HU_clu_data_std, kmeans, method = "silhouette")
# Indicate 2 clusters are the optimal

# 3. Gap statistic
# The gap statistic compares the total intracluster variation for different values of k with their expected values under null reference distribution of the data (i.e. a distribution with no obvious clustering)
# The reference distribution is generate from simulated data
# The large "gap" is a indication of optimal cluster number

gap_stat <- clusGap(HU_clu_data_std, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Print the result
print(gap_stat, method = "firstmax")

fviz_gap_stat(gap_stat)
# The gap statistic result is a bit funny since gap goes up according to increasing of cluster number.


# Three different method indicates 2, 3, 4 could be optimal. Peform K-means clustering for k = 2|3|4 and visualize them
k2 <- kmeans(HU_clu_data_std, centers = 2, nstart = 25)
k3 <- kmeans(HU_clu_data_std, centers = 3, nstart = 25)
k4 <- kmeans(HU_clu_data_std, centers = 4, nstart = 25)

p1 <- fviz_cluster(k2, geom = "point", data = HU_clu_data_std) + 
  ggtitle("K-means, k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = HU_clu_data_std) +
  ggtitle("K-means, k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = HU_clu_data_std) +
  ggtitle("K-means, k = 4")
grid.arrange(p1, p2, p3, nrow = 2)

# Here we extract the result of K = 3 clustering result and append it back to dataset
# If K = 2, the cluster is exactly same as HUN and HUW

HU_ave %>%
  mutate(Cluster = k2$cluster) %>%
  select(SampleID, abb, Cluster) %>%
  head(20)

# Check for outliers and grouping condition of 3 different cluster numbers (2, 3, 4)
# Observations with a large silhouhette Si (almost 1) are very well clustered.
# A small Si (around 0) means that the observation lies between two clusters.
# Observations with a negative Si are probably placed in the wrong cluster.

# Use eclus() to do k-means cluster in a way can be accept by fviz_silhouette()
k2_sil <- eclust(HU_clu_data_std, graph = FALSE,
                 FUNcluster = "kmeans", k = 2)
k3_sil <- eclust(HU_clu_data_std, graph = FALSE,
                 FUNcluster = "kmeans", k = 3)
k4_sil <- eclust(HU_clu_data_std, graph = FALSE,
                 FUNcluster = "kmeans", k = 4)

p1 <- fviz_silhouette(k2_sil)
p2 <- fviz_silhouette(k3_sil)
p3 <- fviz_silhouette(k4_sil)
grid.arrange(p1, p2, p3, nrow = 2)

# Negative Si indicates outliers
# Outliers (potential wrong classification) appear
# Try PAM Partitioning Around Medoids after finishing K-means

#' ## HU PAM
#' Since there are some outliers identified, try Partitioning Around Medoids (PAM). PAM is said to be less sensitive to ouliers and provided a robust alternative to K-means

#+ HU PAM
#### HU PAM ####

# The input is the original daraframe. pam() will take it and calculate dissimilarity distance based on specified method
# Perfome pam on k = 2|3|4 again
pam2 <- pam(HU_clu_data_std, k = 2, diss = FALSE, metric = "euclidean")
pam3 <- pam(HU_clu_data_std, k = 3, diss = FALSE, metric = "euclidean")
pam4 <- pam(HU_clu_data_std, k = 4, diss = FALSE, metric = "euclidean")

p1 <- fviz_cluster(pam2, geom = "point", data = HU_clu_data_std) + 
  ggtitle("PAM, k = 2")
p2 <- fviz_cluster(pam3, geom = "point", data = HU_clu_data_std) + 
  ggtitle("PAM, k = 3")
p3 <- fviz_cluster(pam4, geom = "point", data = HU_clu_data_std) + 
  ggtitle("PAM, k = 4")
grid.arrange(p1, p2, p3, nrow = 2)

# Check the quality of clustering again with Silhouette plot
p1 <- fviz_silhouette(pam2)
p2 <- fviz_silhouette(pam3)
p3 <- fviz_silhouette(pam4)
grid.arrange(p1, p2, p3, nrow = 2)

#' Comparing PAM and K-means\
#' The same cluster number actually produce different grouping (See table below)
#' K-means k=2 provides grouping exactly same as leaf width observation
#' K-means generally closer to observational classes (same for k3 and k4)
#' PAM groups samples  more consistently among different k 
HU_ave %>%
  select(SampleID, Location, abb) %>%
  mutate(k2 = k2$cluster,
         k3 = k3$cluster,
         k4 = k4$cluster,
         pam2 = pam2$clustering,
         pam3 = pam3$clustering,
         pam4 = pam4$clustering) -> HU_kmean_pam
HU_kmean_pam

#' ## HU PCA
#' Start with PCA to identify potential grouping

#+ HU PCA
#### HU PCA ####
# Explanation material: https://personal.utdallas.edu/~herve/abdi-awPCA2010.pdf

# Create data for pca analysis
HU_ave %>%
  select(SampleID, abb, ADM_S:ILength) %>%
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

# Input clusters from Kmean or PAM clustering
# Decide how to present variable arrows
# color_arrow = 1: color the variable arrows according to their contribution
# color_arrow = 2: same color for all variable arrows

# Show the 2 clusters from Kmean and PAM
p1 <- ITV_pca_biplot(HU_ave_pca, HU_kmean_pam$k2, 2)
p2 <- ITV_pca_biplot(HU_ave_pca, HU_kmean_pam$pam2, 2)
grid.arrange(p1, p2, nrow = 1)

# Show 3 cluster with variables contribution
p1 <- ITV_pca_biplot(HU_ave_pca, HU_kmean_pam$k3, 1)
p2 <- ITV_pca_biplot(HU_ave_pca, HU_kmean_pam$pam3, 1)
grid.arrange(p1, p2, nrow = 1) # This code might cause error if the plots window is too small. Just increase it if needed.

# Save the PCA plot colored by growthform. It's the same as k=2 kmeans clustering.
ITV_pca_biplot(HU_ave_pca, HU_kmean_pam$abb, 2) %>%
  ggexport(filename = "../plots/HU_pca_plot.png",
           width = 2000, height = 1500, res = 300)

#' Try the Jones fish approach. Put four traits in the PCA to see PC1 composition

HU_ave1 %>%
  select(LLength, LWidth, CanopyH, ADM_S) -> HU_fish

HU_fish_pca <- prcomp(HU_fish, 
                       center =  TRUE,
                       scale = TRUE)
summary(HU_fish_pca)

ITV_pca_biplot(HU_fish_pca, HU_kmean_pam$abb, 2) %>%
  ggexport(filename = "../plots/HU_fish_plot.png",
         width = 2000, height = 1500, res = 300)

get_pca_var(HU_fish_pca)$coord

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
         Sed = (sed_mean - mean(sed_mean, na.rm = TRUE))/ sd(sed_mean, na.rm = TRUE), .keep = "unused") -> HU_GLM

# Sediment are set to three categories + "0" for NA subtidal site
#for(i in 1:nrow(HU_GLM)){
#  if(HU_GLM$Gravel[i] > 20){
#    HU_GLM$Sed[i] = 1
#  } else if(HU_GLM$Sand[i] > 55){
#    HU_GLM$Sed[i] = 2
#  } else if(HU_GLM$Silt[i] > 50){
#    HU_GLM$Sed[i] = 3
#  } else if(HU_GLM$Clay[i] == 0){
#    HU_GLM$Sed[i] = 0
#  } else{
#    HU_GLM$Sed[i] = 4
#  }
#}
#HU_GLM$Sed = as.factor(HU_GLM$Sed)


# 

HU_GLM %>%
  select(ADM_S:ILength) -> HU_box

HU_GLM %>%
  select(REI:Sed) -> HU_box_env

ggplot(data = pivot_longer(HU_box, everything())) +
  geom_boxplot(aes(x = name, y = value))
# The range of traits are quite different, Below two plots exclude the traits that are too tiny to visualize in this plot

HU_box %>%
  select(-c(CanopyH, LLength, RLength)) %>%
  pivot_longer(everything()) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value))

HU_box %>%
  select(LWidth, Rdia) %>%
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

#' After testing with the simulated data (Can be find in the end of this r script), started to work on real data using the similar model structure.\
#' The stan model is adopted from Bayesian Data Analysis for Cognitive Science Chapter 19(https://vasishth.github.io/bayescogsci/book/ch-mixture.html).
#' The main modification is changing the distribution family to two normal distribution.
#' Non-exchangeable priors for intercepts are used to break the labeling degeneracy.

# Real data, no explained variables
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$LWidth
)

LW_mix_mu <- stan( file = "../stan_file/mixture_regression_mu.stan",
                   data=dat_list, chains=4, cores=4,
                   iter = 4000,
                   warmup = 3500,
                   control=list(adapt_delta=0.99,
                                max_treedepth = 10))
#+ Pure mu output for leaf width
dashboard(LW_mix_mu)
precis(LW_mix_mu, depth = 2)

# real data with rei variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$LWidth,
  rei = HU_GLM$REI)

# The model apply here uses the trick of non-exchangeable prior to break the labeling degeneracy
LW_mix_rei <- stan( file = "../stan_file/mixture_regression_rei.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 2000,
                    warmup = 1500,
                    control=list(adapt_delta=0.99,
                                 max_treedepth = 10))
#+ Leaf width output with one explained variable, relative exposure index (rei)
dashboard(LW_mix_rei)
precis(LW_mix_rei, depth = 2)[1:14,]

# real data with depth variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$LWidth,
  rei = HU_GLM$Dummy
)

# The model apply here uses the trick of non-exchangeable prior to break the labeling degeneracy
LW_mix_de <- stan( file = "../stan_file/mixture_regression_rei.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 2000,
                    warmup = 1500,
                    control=list(adapt_delta=0.99,
                                 max_treedepth = 10))

#+ Leaf width output with one explained variable, depth
dashboard(LW_mix_de)
precis(LW_mix_de, depth = 2)


HU_GLM$Dummy = rnorm(nrow(HU_GLM), mean = 6, sd = 2)
HU_GLM$Dummy2 = rnorm(nrow(HU_GLM), mean = 0, sd = 1)
for(i in 1:nrow(HU_GLM)){
  if(HU_GLM$abb[i] == "HUN"){
    HU_GLM$Dummy[i] = rnorm(1, mean = 0.5, sd = 0.2)
    HU_GLM$Dummy2[i] = rnorm(1, mean = 0.7, sd = 0.2)
  } else{
    HU_GLM$Dummy[i] = rnorm(1, mean = -0.5, sd = 0.2)
    HU_GLM$Dummy2[i] = rnorm(1, mean = -0.7, sd = 0.2)
  }
}


# real data with de, rei variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$LWidth,
  de = HU_GLM$Dummy,
  rei = HU_GLM$Dummy2
)

LW_mix_2var <- stan( file = "../stan_file/mixture_regression_2var.stan",
                         data=dat_list, chains=4, cores=4,
                         iter = 2000,
                         warmup = 1500,
                         control=list(adapt_delta=0.99,
                                      max_treedepth = 15))

dashboard(LW_mix_2var)
precis(LW_mix_2var, depth = 2)
pairs(LW_mix_2var, pars = c("brei1", "brei2", "bde1", "bde2"))


#' In general, the mixture model with explain variable is really hard for the MCMC to walk around the posterior space
#' Change to pure mu model to find the mixing proportion first
#' Here I want to put all traits in the model and see how is the clustering proportion that the model will estimate.
#' The model is modify from this one: https://discourse.mc-stan.org/t/mixture-models/17721/2
#' The correlation matrix (quad_form_diag) is taken from rethinking book for easier understanding.
#' The original approach is also retain in lines with //

#+ Multivariate model without explained variables
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

# the rstan::stan suddenly have problem with running this model, go for cmdstan at above.
#MUL_mix <- stan( file = "../stan_file/Multivariate_mixture.stan",
#                         data=dat_list, chains=4, cores=4,
#                        # iter = 8000,
#                         #warmup = 6000,
#                         control=list(adapt_delta=0.95,
#                                      max_treedepth = 10))
# This stan code can be easily change to Mixture multivariate lognormal regression by changing log(y[N]) in the likelihood section and changing positive_ordered to ordered (log(mean) can be negative).
# This is helpful to keep simulation become all positive.
# However, this will lead to the estimation become a log normal distribution which is generally skew.

# saveRDS(MUL_mix, file = "../models/Multivariate_mixture.Rds")
MUL_mix = readRDS(file = "../models/Multivariate_mixture.Rds")

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
  labs(x = "Leaf Width (mm)", fill = "Growth form") +
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
  labs(x = "Rhizome Diameter (mm)", fill = "Growth form") +
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
#### Two growth form GLM ####

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
  Sed = HUN$Sed)

m.nlw1 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
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
  Sed = HUN$Sed)

m.nlw2 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),
  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nlw2)
precis(m.nlw2, depth = 2)

WAIC(m.nlw1)
WAIC(m.nlw2)
# m.nlw1 (Model 1) is preferred 

# Residual check for m.nlw1
m.nlw1.post <- extract.samples(m.nlw1)
colMeans(m.nlw1.post$mu) -> m.nlw1.pred
m.nlw1.pred - HUN$LWidth -> m.nlw1.resd
plot(m.nlw1.pred, m.nlw1.resd)
qqplot(m.nlw1.pred, HUN$LWidth)
abline(0, 1, col = 2)

# HUW leaf width explained by all not correlated variables
dat <- list(
  LW = HUW$LWidth,
  Rei = HUW$REI,
  De = HUW$Depth,
  Air = HUW$Air,
  Sed = HUW$Sed)

m.wlw1 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wlw1)
precis(m.wlw1, depth = 2)

# HUW leaf width explained by higher resolution variables
dat <- list(
  LW = HUW$LWidth,
  Rei = HUW$REI,
  De = HUW$Depth,
  Sed = HUW$Sed)

m.wlw2 <- ulam(
  alist(
    LW ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wlw2)
precis(m.wlw2, depth = 2)

WAIC(m.wlw1)
WAIC(m.wlw2)
# m.wlw2 (Model 2) is preferred 

m.wlw2.post <- extract.samples(m.wlw2)
colMeans(m.wlw2.post$mu) -> m.wlw2.pred
m.wlw2.pred - HUW$LWidth -> m.wlw2.resd
plot(m.wlw2.pred, m.wlw2.resd)
qqplot(m.wlw2.pred, HUW$LWidth)

# extrapolated Leaf width based on model

# Double check the structure of selected model
# m.nlw1: a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed
# m.wlw2: a + bRei * Rei + bD * De + bS * Sed

#  Control all variables except mean sediment size
S_seq = seq(-1, 2.5, length.out = 50)
S_seq * sd(HU_ave$sed_mean,na.rm = TRUE) + mean(HU_ave$sed_mean, na.rm = TRUE) -> S_org

post <- extract.samples(m.nlw1)
N_lw_sim <- with( post , sapply( 1:50 ,
                              function(i) rnorm( 1e3 , a + bRei * mean(HUN$REI) + bD * mean(HUN$Depth) + bAT * mean(HUN$Ave_temp) + bRun * mean(HUN$Runoff) + bP * mean(HUN$Perci) + bAir * mean(HUN$Air) + bS * S_seq[i], sigma )))

post <- extract.samples(m.wlw2)
W_lw_sim <- with( post , sapply( 1:50 ,
                              function(i) rnorm( 1e3 , a + bRei * mean(HUW$REI) + bD * mean(HUW$Depth) + bS * S_seq[i], sigma)))

data.frame(LW = c(colMeans(N_lw_sim),colMeans(W_lw_sim)),
           Low = c(apply(N_lw_sim, 2, PI)[1, ], 
                   apply(W_lw_sim, 2, PI)[1, ]),
           High = c(apply(N_lw_sim, 2, PI)[2, ],
                    apply(W_lw_sim, 2, PI)[2, ]),
           Sed = rep(S_org, times = 2),
           Abb = c(rep("HUN", 50), rep("HUW", 50))) %>%
  ggplot() +
  geom_line(aes(x = Sed, y = LW, color = Abb), linewidth = 2) +
  geom_ribbon(aes(x = Sed, ymin = Low, ymax = High, color = Abb), , linetype = 2, alpha = 0.1) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  labs(y = "Simulated leaf width (mm)",
       x = "Mean sediment size (mm)",
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
  Sed = HUN$Sed)

m.nrd1 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
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
  Sed = HUN$Sed)

m.nrd2 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.nrd2)
precis(m.nrd2, depth = 2)

WAIC(m.nrd1)
WAIC(m.nrd2)
# m.nrd1 (Model 1) is preferred 

m.nrd1.post <- extract.samples(m.nrd1)
colMeans(m.nrd1.post$mu) -> m.nrd1.pred
m.nrd1.pred - HUN$Rdia -> m.nrd1.resd
plot(m.nrd1.pred, m.nrd1.resd)
qqplot(m.nrd1.pred, HUN$Rdia)

# HUW Rhizome diameter explained by all not correlated variables
dat <- list(
  Rd = HUW$Rdia,
  Rei = HUW$REI,
  De = HUW$Depth,
  Air = HUW$Air,
  Sed = HUW$Sed)

m.wrd1 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
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
  Sed = HUW$Sed)

m.wrd2 <- ulam(
  alist(
    Rd ~ normal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.wrd2)
precis(m.wrd2, depth = 2)

WAIC(m.wrd1)
WAIC(m.wrd2)
# There is only 0.5 difference in WAIC, chose the simpler model
# m.wrd2 (Model 2) is preferred 

m.wrd2.post <- extract.samples(m.wrd2)
colMeans(m.wrd2.post$mu) -> m.wrd2.pred
m.wrd2.pred - HUW$Rdia -> m.wrd2.resd
plot(m.wrd2.pred, m.wrd2.resd)
qqplot(m.wrd2.pred, HUW$Rdia)


# extrapolated Rhizome diameter based on model

# Double check the structure of selected model
# m.nrd1: a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed
# m.wrd2: a + bRei * Rei + bD * De + bS * Sed

#  Control all variables except mean sediment size
S_seq = seq(-1, 2.5, length.out = 50)
S_seq * sd(HU_ave$sed_mean,na.rm = TRUE) + mean(HU_ave$sed_mean, na.rm = TRUE) -> S_org

post <- extract.samples(m.nrd1)
N_rd_sim <- with( post , sapply( 1:50 ,
                                 function(i) rnorm( 1e3 , a + bRei * mean(HUN$REI) + bD * mean(HUN$Depth) + bAT * mean(HUN$Ave_temp) + bRun * mean(HUN$Runoff) + bP * mean(HUN$Perci) + bAir * mean(HUN$Air) + bS * S_seq[i], sigma )))

post <- extract.samples(m.wrd2)
W_rd_sim <- with( post , sapply( 1:50 ,
                                 function(i) rnorm( 1e3 , a + bRei * mean(HUW$REI) + bD * mean(HUW$Depth) + bS * S_seq[i], sigma)))

data.frame(Rd = c(colMeans(N_rd_sim),colMeans(W_rd_sim)),
           Low = c(apply(N_rd_sim, 2, PI)[1, ], 
                   apply(W_rd_sim, 2, PI)[1, ]),
           High = c(apply(N_rd_sim, 2, PI)[2, ],
                    apply(W_rd_sim, 2, PI)[2, ]),
           Sed = rep(S_org, times = 2),
           Abb = c(rep("HUN", 50), rep("HUW", 50))) %>%
  ggplot() +
  geom_line(aes(x = Sed, y = Rd, color = Abb), linewidth = 2) +
  geom_ribbon(aes(x = Sed, ymin = Low, ymax = High, color = Abb), , linetype = 2, alpha = 0.1) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  labs(y = "Simulated Rhizome diameter (mm)",
       x = "Mean sediment size (mm)",
       color = "Growth form") + 
  theme_classic() -> RD_sed
RD_sed

ggsave(Rd_sed, filename = "../plots/RD_sed_effect.png",
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
  Adm = HU_GLM$ADM_S,
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  AveT = HU_GLM$Ave_temp,
  Run = HU_GLM$Runoff,
  Per = HU_GLM$Perci,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed)

m.adm1 <- ulam(
  alist(
    Adm ~ lognormal(mu, sigma),
    mu <- a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.adm1)
precis(m.adm1, depth = 2)

# HU Above ground biomass per shoot explained by higher resolution variables
dat <- list(
  Adm = HU_GLM$ADM_S,
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  Run = HU_GLM$Runoff,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed)

m.adm2 <- ulam(
  alist(
    Adm ~ lognormal(mu, sigma),
    mu <- a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.adm2)
precis(m.adm2, depth = 2)

WAIC(m.adm1)
WAIC(m.adm2)
# m.adm2 (Model 2) is preferred 

m.adm2.post <- extract.samples(m.adm2)
exp(colMeans(m.adm2.post$mu)) -> m.adm2.pred
m.adm2.pred - HU_GLM$ADM_S -> m.adm2.resd
plot(m.adm2.pred, m.adm2.resd)
qqplot(m.adm2.pred, HU_GLM$ADM_S)
abline(a = 0, b = 1)

# HU Below ground biomass per shoot explained by all not correlated variables
dat <- list(
  Bdm = HU_GLM$BDM_S,
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  AveT = HU_GLM$Ave_temp,
  Run = HU_GLM$Runoff,
  Per = HU_GLM$Perci,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed)

m.bdm1 <- ulam(
  alist(
    Bdm ~ lognormal( mu, sigma ),
    mu <- a + bRei * Rei + bD * De + bAT * AveT + bRun * Run + bP * Per + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bAT, bRun, bP, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.bdm1)
precis(m.bdm1, depth = 2)

# HU Below ground biomass per shoot explained by higher resolution variables
dat <- list(
  Bdm = HU_GLM$BDM_S,
  Rei = HU_GLM$REI,
  De = HU_GLM$Depth,
  Run = HU_GLM$Runoff,
  Air = HU_GLM$Air,
  Sed = HU_GLM$Sed)

m.bdm2 <- ulam(
  alist(
    Bdm ~ lognormal(mu, sigma),
    mu <- a + bRei * Rei + bD * De + bRun * Run + bAir * Air + bS * Sed,
    Sed ~ dnorm(nu, sigma_S),
    a ~ normal(1, 1),
    c(bRei, bD, bRun, bAir, bS, nu) ~ normal(0, 1),
    c(sigma, sigma_S) ~ exponential(1)
  ), data = dat, chains = 4, core = 4, 
  control = list(adapt_delta = 0.99),  log_lik = TRUE,
  iter = 8000, warmup = 6000, cmdstan = TRUE
)
dashboard(m.bdm2)
precis(m.bdm2, depth = 2)


WAIC(m.bdm1)
WAIC(m.bdm2)
# The difference of WAIC is only 0.2, chose less complex model
# m.bdm2 (Model 2) is preferred 

m.bdm2.post <- extract.samples(m.bdm2)
exp(colMeans(m.bdm2.post$mu)) -> m.bdm2.pred
m.bdm2.pred - HU_GLM$BDM_S -> m.bdm2.resd
plot(m.bdm2.pred, m.bdm2.resd)
qqplot(m.bdm2.pred, HU_GLM$BDM_S)
abline(a = 0, b = 1, col = "red")

#' Fish biomass index using fish PCA combination
#### Fish biomass ####
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
  theme(legend.position = "none") -> aa
aa %>% ggexport(filename = "../plots/Biomassindex.png",
           width = 4000, height = 3500, res = 600)

#' Blue carbon represented by three traits
#### Blue carbon ####

HU_ave %>%
  select(SampleID:abb,ADM_S, BDM_S, LWidth, CanopyH) %>%
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

sigmas <- c(mix_coef[9,1], mix_coef[10,1]) # standard deviations
Rho<- matrix(c(1, mix_coef[6,1], mix_coef[7,1], 1), nrow=2) # correlation matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

data.frame(N = rmvnorm(10000, c(mix_coef[1, 1],
                                mix_coef[2, 1]), Sigma),
           W = rmvnorm(10000, c(mix_coef[1, 1] + mix_coef[3, 1],
                                mix_coef[2, 1] + mix_coef[4, 1]), Sigma)) -> car_sim

pick_5 <- sample(1:10000, 5000, replace = FALSE)
pick_3 <- sample(1:10000, 3000, replace = FALSE)
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
  theme(legend.position = "none") -> bb
bb
bb %>% ggexport(filename = "../plots/Bluecarbon.png",
                width = 4000, height = 3500, res = 600)

#' For the other six traits. Multivariate models from brms are used to estimate the distribution of the traits.
#' Below there are four models for comparison
#' 1. All seven traits ~ all environment covariates in Gaussian family
#' 2. All seven traits ~ environment covariates (REI, Depth, and Sediments) with more variability in Gaussian family
#' 3. All seven traits ~ all environment covariates in log-normal family
#' 4. All seven traits ~ environment covariates (REI, Depth, and Sediments) with more variability in log-normal family
#' Then compare 4 models using loo function and select the best one

HU_ave %>%
  select(!c(mean_max, mean_min)) %>%
  mutate(REI = scale(total_REI)[, 1],
         Depth =  scale(depth)[, 1],
         Ave_temp = scale(ave_temp)[, 1],
         Var_temp = scale(mean_var)[, 1],
         Runoff = scale(runoff)[, 1],
         Perci = scale(perci)[, 1],
         Air = scale(air_exposure), .keep = "unused") %>%
  mutate(Sed = Gravel) %>%
  mutate_all(~replace(., is.na(.), 0)) -> HU_mul

for(i in 1:nrow(HU_mul)){
  if(HU_mul$Gravel[i] > 20){
    HU_mul$Sed[i] = 1
  } else if(HU_mul$Sand[i] > 55){
    HU_mul$Sed[i] = 2
  } else if(HU_mul$Silt[i] > 50){
    HU_mul$Sed[i] = 3
  } else if(HU_mul$Clay[i] == 0){
    HU_mul$Sed[i] = 0
  } else{
    HU_mul$Sed[i] = 4
  }
}
HU_mul$Sed = as.factor(HU_mul$Sed)

#' Model 1

priors <- prior(normal(0 ,10), class = "Intercept") +
  prior(normal(0, 1), class = "b")


#prior(gamma(2, 1), class = "sigma")

bform1 <- 
  bf(mvbind(ADM_S, BDM_S, CanopyH, LLength, RLength, ILength) ~ (1|Location) + Sed) +
  set_rescor(TRUE)

fit1 <- brm(bform1,
            data = HU_mul, 
            #prior = priors,
            #family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4, 
            cores = 4,
            iter = 2000,
            warmup = 1000)
saveRDS(fit1, file = "../models/Multivariate1.Rds")
fit1 <- readRDS(file = "../models/Multivariate1.Rds")

# fit1 <- add_criterion(fit1, "loo", moment_match = TRUE)

summary(fit1)
conditional_effects(fit1)

#' Model 2

bform2 <- 
  bf(mvbind(Shoot_number, ADM, BDM, CanopyH, LLength, RLength, ILength) ~ REI + Depth + as.factor(Sed)) +
  set_rescor(TRUE)

fit2 <- brm(bform2,
            data = HU_mul, 
            #prior = priors,
            #family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4, 
            cores = 4,
            iter = 2000,
            warmup = 1000)
saveRDS(fit2, file = "../models/Multivariate2.Rds")
fit2 <- readRDS(file = "../models/Multivariate2.Rds")


# fit2 <- add_criterion(fit2, "loo", moment_match = TRUE)
summary(fit2)

#' Model 3

bform3 <- 
  bf(mvbind(Shoot_number, ADM, BDM, CanopyH, LLength, RLength, ILength) ~ REI + Depth + Ave_temp + Var_temp + Runoff + Perci + Sed) +
  set_rescor(FALSE)

fit3 <- brm(bform3,
            data = HU_mul, 
            #prior = priors,
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4, 
            cores = 4,
            iter = 2000,
            warmup = 1000)
saveRDS(fit3, file = "../models/Multivariate3.Rds")
fit3 <- readRDS(file = "../models/Multivariate3.Rds")

# fit3 <- add_criterion(fit3, "loo", moment_match = TRUE)
summary(fit3)

#' Model 4
bform4 <- 
  bf(mvbind(Shoot_number, ADM, BDM, CanopyH, LLength, RLength, ILength) ~ REI + Depth + as.factor(Sed)) +
  set_rescor(FALSE)

fit4 <- brm(bform4,
            data = HU_mul, 
            #prior = priors,
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4, 
            cores = 4,
            iter = 2000,
            warmup = 1000)
saveRDS(fit4, file = "../models/Multivariate4.Rds")
fit4 <- readRDS(file = "../models/Multivariate4.Rds")

# fit4 <- add_criterion(fit4, "loo", moment_match = TRUE)
summary(fit4)

#' Model 5
bform5 <- 
  bf(mvbind(Shoot_number, ADM, BDM, CanopyH, LLength, RLength, ILength) ~ 1) +
  set_rescor(FALSE)

fit5 <- brm(bform5,
            data = HU_mul, 
            #prior = priors,
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4, 
            cores = 4,
            iter = 2000,
            warmup = 1000)
saveRDS(fit5, file = "../models/Multivariate5.Rds")
fit5 <- readRDS(file = "../models/Multivariate5.Rds")

summary(fit5)

loo(fit1, fit2, fit3, fit4, fit5)
# loo stands for "leave one out" estimate how the model look like with one observed data is kick out
# "pareto k diagnostic" indicates if any observed data is hard for loo to kick out and simulate. If there is any data seat at "bad", "very bad", the aic might not be trustable -> refine model
# loo_compare() put the best model at the top row (base on deviance(elpd_diff), the higher the better)

# fit3 is the best model of the 5
fit3$fit %>% stan_trace()
fit3$fit %>% stan_ac()

conditions <- data.frame(Sed = c(1, 2, 3, 4))
fit3 %>% conditional_effects( resp = "ADM",
                             conditions = conditions) 

preds <- fit3 %>% posterior_predict(nasamples = 250, summary = FALSE)
fit3.resids <- createDHARMa(simulatedResponse = t(preds[,,"ADM"]),
                            observedResponse = HU_mul$Shoot_number,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
fit.resids %>% plot()






#+ Exploring mixture model with simulated data and Extra example not run, eval = FALSE, include = FALSE
#### Not run simulation and example ####


#' Before run the model with real data, start with simulated data set.

#+ Simulated data, comments = "", warnings = FALSE

set.seed(1429850)
N = 100
# Simulate two types of depth
de <- sample(c(4, 5, 3), N, replace = TRUE)
# Mean leaf width is depends on the depth
z <- rbern(n = N, prob = 0.4)

lw <- if_else(z == 1,
              rnorm(N,
                    mean = 5 + 4 * de,
                    sd = 1),
              rnorm(N, 
                    mean = 1 + 0.8 *de,
                    sd = 2))

# Leaf width data show 2 set of values according to their means. The proportion is 4:6

#+ Pure mu model, results = "hide", warning = FALSE, message = FALSE, comment = ""
# First run a model without explained variable (de)
dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw
)

LWM <- cstan( file = "../stan_file/mixture_regression_mu.stan",
              data=dat_list, chains=4, cores=4,
              iter = 4000,
              warmup = 3000,
              control=list(adapt_delta=0.95,
                           max_treedepth = 10))
#+ Pure mu output
dashboard(LWM)
precis(LWM, depth = 2)
# "theta" is mixing proportion. The estimation is really close to true value (4:6)
# "sigma" is standard deviation of each group of the data. Less precise than theta estimation. (group[1] is set to 2. group[2]is set to 1). If increase the standard deviation in the simulating process, the model estimation will deviate from the real data more.

#+ Explained variables included, results = "hide", warning = FALSE, message = FALSE, comment = ""
# This model include the explained variable (de)
dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw,
  de = de
)

LW_mix_sim_test <- cstan( file = "../stan_file/mixture_regression_simple_test.stan",
                          data=dat_list, chains=4, cores=4,
                          iter = 10000,
                          warmup = 9000,
                          control=list(adapt_delta=0.95,
                                       max_treedepth = 10))
#+ Explained variables included output
dashboard(LW_mix_sim_test)
precis(LW_mix_sim_test, depth = 2)

# Theta, sigma, and mu are still doing well in this model.
# The coefficient (b = 2 or 5) used to determine mu in the simulation did not pick up by the model well enough. Need a bit more work on model structure,

# Simulate Leaf width depends on 2 explained variables
N <- 200

de <- c(rnorm(200*0.3, 1, 1), rnorm(200*0.7, 5, 1))
rei <- c(rnorm(200*0.3, 5, 1), rnorm(200*0.7, 2, 1))

# Parameters true values:
alpha <- 1
bde1 <- 0.05
brei1 <- 0.3
sigma1 <- 0.4
sigma2 <- 0.5
gamma <- 5.8
bde2 <- 0.4
brei2 <- 0.8

z <- rbern(n = N, prob = 0.3)

lw <- c(rnorm(N*0.3, mean = alpha + bde1 * mean(de[1:60]) + brei1 * mean(rei[1:60]), sd = sigma1),
        rnorm(N*0.7, mean = gamma + bde2 * mean(de[61:200]) + brei2 * mean(rei[61:200]), sd = sigma2))

dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw
)

LW_mix_mu_test <- cstan( file = "../stan_file/mixture_regression_mu.stan",
                         data=dat_list, chains=4, cores=4,
                         iter = 4000,
                         warmup = 3500,
                         control=list(adapt_delta=0.95,
                                      max_treedepth = 10))
#+ Pure mu output
dashboard(LW_mix_mu_test)
precis(LW_mix_mu_test, depth = 2)

# This model include 2 explained variable (de and rei)
dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw,
  de = de,
  rei = rei
)

LW_mix_2_test <- cstan( file = "../stan_file/mixture_regression_2var.stan",
                        data=dat_list, chains=4, cores=4,
                        iter = 2000,
                        warmup = 1500,
                        control=list(adapt_delta=0.95,
                                     max_treedepth = 10))
#+ Explained variables included output
dashboard(LW_mix_2_test)
precis(LW_mix_2_test, depth = 2)

# Simulate Leaf width depends on 3 explained variables and one of them is ordered categorical variable
N <- 200

de <- sample(c(-1, -0.5, 0, 1, 0.5), N, replace = TRUE)
rei <- sample(c(1, 3, 5, 7), N, replace = TRUE)
air <- sample(c(1, 2, 3), N, replace = TRUE, prob = c(0.5, 0.2, 0.3))

# Parameters true values:
alpha <- 5.8
bde1 <- 0.1
brei1 <- 0.3
bair1 <- 0.2
sigma1 <- 0.4
sigma2 <- 0.5
gamma <- 3
bde2 <- 0.4
brei2 <- 0.8
bair2 <- 0.5

delta_j1 <- c(0, 0.2, 0.3, 0.5)
delta_j2 <- c(0, 0.4, 0.2, 0.4)

z <- rbern(n = N, prob = 0.3)

lw = c(1:N)
for(i in 1:N){
  if(z[i] == 1 & air[i] == 1){
    lw[i] = rnorm(1, mean = alpha + bde1 * de[i] + brei1 * rei[i] + bair1 *sum(delta_j1[1:2]), sd = sigma1)} else if(z[i] == 1 & air[i] == 2){
      lw[i] = rnorm(1, mean = alpha + bde1 * de[i] + brei1 * rei[i] + bair1 *sum(delta_j1[1:3]), sd = sigma1)} else if(z[i] == 1 & air[i] == 3){
        lw[i] = rnorm(1, mean = alpha + bde1 * de[i] + brei1 * rei[i] + bair1 *sum(delta_j1[1:4]), sd = sigma1)} else if(z[i] == 0 & air[i] == 1){
          lw[i] = rnorm(1, mean = gamma + bde2 * de[i] + brei2 * rei[i] + bair2 *sum(delta_j2[1:2]), sd = sigma2)} else if(z[i] == 0 & air[i] == 1){
            lw[i] = rnorm(1, mean = gamma + bde2 * de[i] + brei2 * rei[i] + bair2 *sum(delta_j2[1:2]), sd = sigma2)} else if(z[i] == 0 & air[i] == 2){
              lw[i] = rnorm(1, mean = gamma + bde2 * de[i] + brei2 * rei[i] + bair2 *sum(delta_j2[1:3]), sd = sigma2)} else{
                lw[i] = rnorm(1, mean = gamma + bde2 * de[i] + brei2 * rei[i] + bair2 *sum(delta_j2[1:4]), sd = sigma2)}
}

mean(lw[z == 1])
mean(lw[z == 0])

dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw
)

LW_mix_mu_test <- cstan( file = "../stan_file/mixture_regression_mu.stan",
                         data=dat_list, chains=4, cores=4,
                         iter = 4000,
                         warmup = 3000,
                         control=list(adapt_delta=0.95,
                                      max_treedepth = 10))
#+ Pure mu output
dashboard(LW_mix_mu_test)
precis(LW_mix_mu_test, depth = 2)

# This model include 3 explained variable (de, rei, and air)
dat_list <- list(
  K = 2, # number of mixture component
  N = length(lw),
  LW = lw,
  de = de,
  rei = rei,
  air = air,
  pair1 = rep(2, 3), # Prior for ordered categories process
  pair2 = rep(2, 3) # Prior for ordered categories process
)

LW_mix_cat_test <- cstan( file = "../stan_file/mixture_regression_cat_test.stan",
                          data=dat_list, chains=4, cores=4,
                          iter = 8000,
                          warmup = 7500,
                          control=list(adapt_delta=0.95,
                                       max_treedepth = 10))
#+ 3 Explained variables included output
dashboard(LW_mix_cat_test)
precis(LW_mix_cat_test, depth = 2)



# brms on simulated data. surprisingly not doing too good
priors <- c(
  prior(normal(0, 7), Intercept, dpar = mu1),
  prior(normal(4, 7), Intercept, dpar = mu2)
)
mix <- mixture(gaussian, gaussian)
test_dat = data.frame(lw, de = rnorm(100))
fit1_test <- brm(bf(lw ~ de), 
                 test_dat,
                 family = mix,
                 prior = priors,
                 iter = 1000,
                 warmup = 500,
                 chains = 3,
                 cores = 3,
                 control = list(adapt_delta = 0.95))

update(fit1_test)
summary(fit1_test)
brms::stancode(fit1_test)

preds <- fit1_test %>% posterior_predict(nasamples = 250, summary = FALSE)
fit1_test.resids <- createDHARMa(simulatedResponse = t(preds),
                                 observedResponse = lw,
                                 fittedPredictedResponse = apply(preds, 2, median),
                                 integerResponse = FALSE)
fit1_test.resids %>% plot()
fit1_test.resids %>% testDispersion()

# example of ordered categorical expained variables from rethinking package

data(Trolley)
d <- Trolley

edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
d$edu_new <- edu_levels[ d$edu ]

dat <- list(
  R = d$response ,
  action = d$action,
  intention = d$intention,
  contact = d$contact,
  E = as.integer( d$edu_new ), # edu_new as an index
  alpha = rep( 2 , 7 ) ) # delta prior

m12.6 <- ulam(
  alist(
    R ~ ordered_logistic( phi , kappa ),
    phi <- bE*sum( delta_j[1:E] ) + bA*action + bI*intention + bC*contact,
    kappa ~ normal( 0 , 1.5 ),
    c(bA,bI,bC,bE) ~ normal( 0 , 1 ),
    vector[8]: delta_j <<- append_row( 0 , delta ),
    simplex[7]: delta ~ dirichlet( alpha )
  ), data=dat , chains=4 , cores=4 )

precis(m12.6, depth = 2)
rethinking::stancode(m12.6)


dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leafwidth,
  D = HU_GLM$depth,
  REI = as.numeric(scale(HU_GLM$total_REI)),
  AIR = as.integer(HU_GLM$air_exposure+1),
  alpha = rep(2, 4), # Prior for ordered categories process
  GA = HU_GLM$Gravel,
  SA = HU_GLM$Sand,
  SI = HU_GLM$Silt,
  CL = HU_GLM$Clay
)

LW_mix_cat <- stan( file = "../stan_file/mixture_regression.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 5000,
                    warmup = 4000,
                    control=list(adapt_delta=0.9,
                                  max_treedepth = 14) )
dashboard(LW_mix_cat)
precis(LW_mix_cat, depth = 2)




a <- 3.5 # average morning wait time
bgood <- (-1) # average difference afternoon wait time for good cafe
bbad <- (1) # average difference afternoon wait time for bad cafe
sigma_a <- 1 # std dev in intercepts
sigma_b <- 0.5 # std dev in slopes
rho <- (-0.7) # correlation between intercepts and slopes

Mugood <- c( a , bgood )
Mubad <- c(a, bbad)
cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )

sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

N_cafes <- 20


set.seed(5) # used to replicate example
vary_effectsgood <- mvrnorm( N_cafes , Mugood , Sigma )
vary_effectsbad <- mvrnorm( N_cafes , Mubad , Sigma )

a_cafe <- c(vary_effectsgood[1:10,1], vary_effectsbad[1:10,1])
b_cafegood <- vary_effectsgood[1:10,2]
b_cafebad <- vary_effectsbad[1:10, 2]

set.seed(22)
N_visits <- 10
goodbad <- rep(c(5, -5), each = 100)
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep( 1:N_cafes , each=N_visits )
mu <- c((a_cafe[cafe_id] + b_cafegood[cafe_id]*afternoon)[1:100],  (a_cafe[cafe_id] + b_cafebad[cafe_id]*afternoon)[1:100])
sigma <- 0.5 # std dev within cafes
wait <- rnorm( N_visits*N_cafes , mu , sigma )
d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )



set.seed(867530)
m14.1 <- ulam(
  alist(
    wait ~ normal( mu , sigma ),
    mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
    c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ),
    a ~ normal(5,2),
    b ~ normal(-1,0.5),
    sigma_cafe ~ exponential(1),
    sigma ~ exponential(1),
    Rho ~ lkj_corr(2)
  ) , data=d , chains=4 , cores=4 )
precis(m14.1, depth = 3)
rethinking::stancode(m14.1)

set.seed(73)
N <- 500
U_sim <- rnorm( N )
Q_sim <- sample( 1:4 , size=N , replace=TRUE )
E_sim <- rnorm( N , U_sim + Q_sim )
W_sim <- rnorm( N , U_sim + 0*E_sim )
dat_sim <- list(
  W=scale(W_sim) ,
  E=scale(E_sim) ,
  Q=scale(Q_sim))

m14.6 <- ulam(
  alist(
    c(W,E) ~ multi_normal( c(muW,muE) , Rho , Sigma ),
    muW <- aW + bEW*E,
    muE <- aE + bQE*Q,
    c(aW,aE) ~ normal( 0 , 0.2 ),
    c(bEW,bQE) ~ normal( 0 , 0.5 ),
    Rho ~ lkj_corr( 2 ),
    Sigma ~ exponential( 1 )
  ), data=dat_sim , chains=4 , cores=4 )

precis( m14.6 , depth=3 )
rethinking::stancode(m14.6)
