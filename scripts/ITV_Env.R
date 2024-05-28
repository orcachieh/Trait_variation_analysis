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
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra) # arranging plots
library(psych) # pairwise scatter plot
library(DHARMa)     #for residual diagnostics

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
  # arrange column
  select(SampleID:Shoot_number, ADM:BDM,
         Canopyheight:Internodelength, trials,
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
  # arrange column
  select(SampleID:Shoot_number, ADM:ILength,
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

#' ## HU coefficient of variation
#+ HU CV analysis
#### HU CV analysis ####

HU_long %>% select(12:17) %>% colMeans(., na.rm = TRUE) -> mean1
HU_long %>% select(12:17) %>% apply(2, sd, na.rm = TRUE) -> sd1

HU_ave %>% select(9:11) %>% colMeans(., na.rm = TRUE) -> mean2
HU_ave %>% select(9:11) %>% apply(2, sd, na.rm = TRUE) -> sd2

HU_CV <- t((cbind(rbind(mean1, sd1),
                  rbind(mean2, sd2)))) %>%
  as.data.frame() %>%
  mutate(CV = sd1/mean1, Trait = row.names(.))

ggplot(HU_CV) +
  geom_point(aes(x = reorder(Trait, CV), y = CV)) +
  theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#' ## HU Clusteting

#+ HU Clustering analysis
#### HU Clustering####

# Use average repeated measurements data to create a cluster data for HU and remove NA
HU_ave %>%
  select(SampleID, abb, Shoot_number:ILength) %>%
  na.omit() -> HU_clu_data

# Calculate distance matrix
# 1. Standardize (Normalize) each column (x-mean)/sd
# 2. Calculate dissimilarity matrix based on Euclidean distance

# remove first two columns (ID columns) then scale
HU_clu_data_std<- scale(HU_clu_data[, -c(1, 2)],center=TRUE,scale=TRUE)
head(HU_clu_data_std)

# Caluluate dissimilarity distance matrix and visulize it
distance = get_dist(HU_clu_data_std, method = "euclidean")
fviz_dist(distance, 
          gradient = list(low = "#00AFBB", 
                          mid = "white", 
                          high = "#FC4E07"))

#' Hierarchical clustering analysis

#+ H-clustering
#### HU hiraechical clustering ####

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
# Indicate 3 clusters are the optimal

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
# Nevertheless, firstmax method indicates 2 cluster is optimal


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
  mutate(Cluster = k3$cluster) %>%
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
# Explaination material: https://personal.utdallas.edu/~herve/abdi-awPCA2010.pdf

# Create data for pca analysis
HU_ave %>%
  select(SampleID, abb, Shoot_number:ILength) %>%
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

# Variables contribution to PC1, PC2 (another way of visulization)
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
p1 <- ITV_pva_biplot(HU_ave_pca, HU_kmean_pam$k2, 1)
p2 <- ITV_pva_biplot(HU_ave_pca, HU_kmean_pam$pam2, 2)
grid.arrange(p1, p2, nrow = 1)

# Show 3 cluster with variables contribution
p1 <- ITV_pva_biplot(HU_ave_pca, HU_kmean_pam$k3, 1)
p2 <- ITV_pva_biplot(HU_ave_pca, HU_kmean_pam$pam3, 1)
grid.arrange(p1, p2, nrow = 1) # This code might cause error if the plots window is too small. Just increase it if needed.

#' ## HU GLMM
#' Work on GLMM. Perform following steps:\
#' 1. Check environmental variables correlation
#' 2. Identify family\
#' 3. Run model\
#' 4. Compared model prediction with original data\

#+ HU data exploration
#### HU Mixture model ####

HU_long %>%
  na.omit() -> HU_GLM

HU_GLM_box <- HU_GLM
HU_GLM_box[c(19, 20, 22:25, 30, 31)] <- lapply(HU_GLM_box[c(19, 20, 22:25, 30, 31)] , function(x) c(scale(x)))

# Separate data to with repeated measurement data frame and non repeated measuremtn data frame for plotting
HU_GLM_box %>%
  select(Canopyheight:Internodelength) %>%
  distinct() -> HU_GLM_re

HU_GLM_box %>%
  select(Shoot_number:BDM) %>%
  distinct() -> HU_GLM_nonre

HU_GLM_box %>%
  select(total_REI, depth, ave_temp:mean_var, runoff, perci) %>%
  distinct() -> HU_GLM_env

ggplot(data = pivot_longer(HU_GLM_re, everything())) +
  geom_boxplot(aes(x = name, y = value))

ggplot(data = pivot_longer(HU_GLM_nonre, everything())) +
  geom_boxplot(aes(x = name, y = value))

ggplot(data = pivot_longer(HU_GLM_env, everything())) +
  geom_boxplot(aes(x = name, y = value))


pairs.panels(HU_GLM_re,
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)
pairs.panels(HU_GLM_nonre,
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE)
# Check explained variable (x1...xn) correlation
pairs.panels(HU_GLM_env,
             gap = 0,
             density = TRUE,  # show density plots
             hist.col = "#00AFBB",
             ellipses = TRUE) # show correlation ellipses)


#' After processes above, I identifies that a few explained variables need to be log before scale, and all will need to be scaled except categorical variables or proportional variables\
#' 1. total_REI: skew distribution -> scale(log(total_REI))\
#' 2. Air_exposure is ordered categorical data. Need to be treated differently in the model\
#' 3. Sediment need to sum up to 100. Not sure how to deal with this in model. If can't figure it out, maybe just separate to 4 types depends on the dominant sediment types.\


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

de <- sample(c(-1, -0.5, 0, 1, 0.5), N, replace = TRUE)
rei <- sample(c(1, 3, 5, 7), N, replace = TRUE)

# Parameters true values:
alpha <- 5.8
bde1 <- 0.05
brei1 <- 0.3
sigma1 <- 0.4
sigma2 <- 0.5
gamma <- 3
bde2 <- 0.4
brei2 <- 0.8

z <- rbern(n = N, prob = 0.3)

lw <- if_else(z == 1,
              rnorm(N,
                    mean = alpha + bde1 * de + brei1 *rei,
                    sd = sigma1),
              rnorm(N, 
                    mean = gamma + bde2 * de + brei2 *rei,
                    sd = sigma2))

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

LW_mix_2_test <- cstan( file = "../stan_file/mixture_regression_2_test.stan",
                        data=dat_list, chains=4, cores=4,
                        iter = 4000,
                        warmup = 3500,
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

#' After testing with the simulated data, started to work on real data using the similar model structure.\
#' The stan model is adopted from Bayesian Data Analysis for Cognitive Science Chapter 19(https://vasishth.github.io/bayescogsci/book/ch-mixture.html).
#' The main modification is changing the distribution family to two normal distribution.
#' Non-exchangeable priors for intercepts are used to break the labeling degeneracy.

# Real data, no explained variables
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leafwidth
)

LW_mix_mu <- stan( file = "../stan_file/mixture_regression_mu.stan",
                   data=dat_list, chains=4, cores=4,
                   iter = 4000,
                   warmup = 3500,
                   control=list(adapt_delta=0.95,
                                max_treedepth = 10))
#+ Pure mu output for leaf width
dashboard(LW_mix_mu)
precis(LW_mix_mu, depth = 2)

# real data with rei variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leafwidth,
  rei = as.numeric(scale(HU_GLM$total_REI))
)
saveRDS(dat_list, "LW_mix_rei.RData")

# The model apply here uses the trick of non-exchangeable prior to break the labeling degeneracy
LW_mix_rei <- stan( file = "../stan_file/mixture_regression_rei.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 4000,
                    warmup = 3500,
                    control=list(adapt_delta=0.99,
                                 max_treedepth = 10))
#+ Leaf width output with one explained variable, relative exposure index (rei)
dashboard(LW_mix_rei)
precis(LW_mix_rei, depth = 2)

# real data with depth variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leafwidth,
  rei = HU_GLM$depth
)

# The model apply here uses the trick of non-exchangeable prior to break the labeling degeneracy
LW_mix_de <- stan( file = "../stan_file/mixture_regression_rei.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 4000,
                    warmup = 3500,
                    control=list(adapt_delta=0.99,
                                 max_treedepth = 10))

#+ Leaf width output with one explained variable, depth
dashboard(LW_mix_de)
precis(LW_mix_de, depth = 2)

extract(LW_mix_rei) -> aa

# real data with de, rei variable
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leafwidth,
  de = HU_GLM$depth,
  rei = as.numeric(scale(HU_GLM$total_REI))
)
saveRDS(dat_list, "LW_mix_2var.RData")

LW_mix_2var <- stan( file = "../stan_file/mixture_regression_2variables.stan",
                    data=dat_list, chains=4, cores=4,
                    iter = 2000,
                    warmup = 1500,
                    control=list(adapt_delta=0.99,
                                 max_treedepth = 10))

dashboard(LW_mix_2var)
precis(LW_mix_2var, depth = 2)




# Leaf length
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Leaflength
)

LL <- cstan( file = "../stan_file/mixture_regression_mu.stan",
                         data=dat_list, chains=4, cores=4,
                         iter = 4000,
                         warmup = 3000,
                         control=list(adapt_delta=0.95,
                                      max_treedepth = 10))
#+ Pure mu output
dashboard(LL)
precis(LL, depth = 2)


# Above Biomass
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$ADM
)

ADM <- cstan( file = "../stan_file/mixture_regression_mu.stan",
             data=dat_list, chains=4, cores=4,
             iter = 4000,
             warmup = 3000,
             control=list(adapt_delta=0.95,
                          max_treedepth = 10))
#+ Pure mu output
dashboard(ADM)
precis(ADM, depth = 2)

# Below Biomass
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$BDM
)

BDM <- cstan( file = "../stan_file/mixture_regression_mu.stan",
              data=dat_list, chains=4, cores=4,
              iter = 4000,
              warmup = 3000,
              control=list(adapt_delta=0.95,
                           max_treedepth = 10))
#+ Pure mu output
dashboard(BDM)
precis(BDM, depth = 2)

# Shoot number
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Shoot_number
)

SN <- cstan( file = "../stan_file/mixture_regression_mu.stan",
              data=dat_list, chains=4, cores=4,
              iter = 4000,
              warmup = 3000,
              control=list(adapt_delta=0.95,
                           max_treedepth = 10))
#+ Pure mu output
dashboard(SN)
precis(SN, depth = 2)

# Canopy height
dat_list <- list(
  K = 2, # number of mixture component
  N = nrow(HU_GLM),
  LW = HU_GLM$Canopyheight
)

CN <- cstan( file = "../stan_file/mixture_regression_mu.stan",
             data=dat_list, chains=4, cores=4,
             iter = 4000,
             warmup = 3000,
             control=list(adapt_delta=0.95,
                          max_treedepth = 10))
#+ Pure mu output
dashboard(CN)
precis(CN, depth = 2)


precis(LWM, depth = 2)
precis(LL, depth = 2)
precis(ADM, depth = 2)
precis(BDM, depth = 2)
precis(SN, depth = 2)
precis(CN, depth = 2)

WLF_LW_m = 28.3
WLF_LW_var = 0.55
NLF_LW_m = 0.59
NLF_LW_var = 0.23

WLF_LL_m = 59.16
WLF_LL_var = 30.26
NLF_LL_m = 35.32
NLF_LL_var = 23.09

WLF_ADM_m = 59.40
WLF_ADM_var = 220.51
NLF_ADM_m = 8.74
NLF_ADM_var = 86.93

WLF_BDM_m = 47.31
WLF_BDM_var = 948.94
NLF_BDM_m = 36.03
NLF_BDM_var = 251.50

WLF_SN_m = 35.40
WLF_SN_var = 18.39
NLF_SN_m = 98.57
NLF_SN_var = 18.39

WLF_CN_m = 71.92
WLF_CN_var = 53.02
NLF_CN_m = 29.57
NLF_CN_var = 71.72


# Pure WLF
P_WLF <- matrix(c(abs(rnorm(1000, WLF_LW_m, WLF_LW_var)),
               abs(rnorm(1000, WLF_LL_m, WLF_LL_var)),
               abs(rnorm(1000, WLF_SN_m, WLF_SN_var))),
               nrow = 1000, ncol = 3, byrow = FALSE) %>%
  as.data.frame() %>%
  mutate(P_WLF_Fish = V1 * V2 * V3)

# 70% WLF, 30% NLF
WNLF_7030 <- matrix(c(c(abs(rnorm(700, WLF_LW_m, WLF_LW_var)),
                        abs(rnorm(300, NLF_LW_m, NLF_LW_var))),
                    c(abs(rnorm(700, WLF_LL_m, WLF_LL_var)),
                      abs(rnorm(300, NLF_LL_m, NLF_LL_var))),
                    c(abs(rnorm(700, WLF_SN_m, WLF_SN_var)),
                      abs(rnorm(300, NLF_SN_m, NLF_SN_var)))),
                nrow = 1000, ncol = 3, byrow = FALSE) %>%
  as.data.frame() %>%
  mutate(WNLF_7030_Fish = V1 * V2 * V3)

# 50% WLF, 50% NLF
WNLF_5050 <- matrix(c(c(abs(rnorm(500, WLF_LW_m, WLF_LW_var)),
                        abs(rnorm(500, NLF_LW_m, NLF_LW_var))),
                      c(abs(rnorm(500, WLF_LL_m, WLF_LL_var)),
                        abs(rnorm(500, NLF_LL_m, NLF_LL_var))),
                      c(abs(rnorm(500, WLF_SN_m, WLF_SN_var)),
                        abs(rnorm(500, NLF_SN_m, NLF_SN_var)))),
                    nrow = 1000, ncol = 3, byrow = FALSE) %>%
  as.data.frame() %>%
  mutate(WNLF_5050_Fish = V1 * V2 * V3)


# 30% WLF, 70% NLF
WNLF_3070 <- matrix(c(c(abs(rnorm(300, WLF_LW_m, WLF_LW_var)),
                        abs(rnorm(700, NLF_LW_m, NLF_LW_var))),
                      c(abs(rnorm(300, WLF_LL_m, WLF_LL_var)),
                        abs(rnorm(700, NLF_LL_m, NLF_LL_var))),
                      c(abs(rnorm(300, WLF_SN_m, WLF_SN_var)),
                        abs(rnorm(700, NLF_SN_m, NLF_SN_var)))),
                    nrow = 1000, ncol = 3, byrow = FALSE) %>%
  as.data.frame() %>%
  mutate(WNLF_3070_Fish = V1 * V2 * V3)

# Pure NLF
P_NLF <- matrix(c(abs(rnorm(1000, NLF_LW_m, NLF_LW_var)),
                  abs(rnorm(1000, NLF_LL_m, NLF_LL_var)),
                  abs(rnorm(1000, NLF_SN_m, NLF_SN_var))),
                nrow = 1000, ncol = 3, byrow = FALSE) %>%
  as.data.frame() %>%
  mutate(P_NLF_Fish = V1 * V2 * V3)


data_plot <- as.data.frame(cbind(P_WLF$P_WLF_Fish, 
                   WNLF_7030$WNLF_7030_Fish, 
                   WNLF_5050$WNLF_5050_Fish,
                   WNLF_3070$WNLF_3070_Fish,
                   P_NLF$P_NLF_Fish))
colnames(data_plot) <- c("Pure WLF", "WNLF 70%30%",
                         "WNLF 50%50%", "WNLF 30%70%",
                         "Pure NLF")
pivot_longer(data_plot, cols = 1:5, names_to = "type", values_to = "values") -> d_plot

d_plot$values <- d_plot$values/10000


ggplot(data = d_plot) +
  geom_density(aes(x = values, col = type)) 




#+ Extra example not run, eval = FALSE, include = FALSE

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


