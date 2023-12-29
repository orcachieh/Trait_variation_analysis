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
# The result is a bit funny since gap goes up according to increasing of cluster number.
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
  mutate(Cluster = k3$cluster)

# Check for outliers and grouping condition of 3 different cluster numbers (2, 3, 4)
# Observations with a large silhouhette Si (almost 1) are very well clustered.
# A small Si (around 0) means that the observation lies between two clusters.
# Observations with a negative Si are probably placed in the wrong cluster.

# Use eclus() to do k-means cluster in a way can be accept by fviz_silhouette()
k2_sil <- eclust(HU_clu_data_std, FUNcluster = "kmeans", k = 2)
k3_sil <- eclust(HU_clu_data_std, FUNcluster = "kmeans", k = 3)
k4_sil <- eclust(HU_clu_data_std, FUNcluster = "kmeans", k = 4)

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
         pam4 = pam4$clustering)


#' ## HU PCA
#' Start with PCA to identify potential grouping

#+ HU PCA
#### HU PCA ####

HU_ave %>%
  select(SampleID, abb, Shoot_number:ILength) %>%
  na.omit() -> HU_pca_ave

HU_long

scaled_HU_pca_ave <- apply(HU_pca_ave[, -c(1, 2)], 2, scale)
head(scaled_HU_pca_ave)

set.seed(111)
ind <- sample(2, nrow(HU_pca_data),
              replace = TRUE,
              prob = c(0.8, 0.2))
training <- na.omit(HU_pca_data[ind==1,])
testing <- na.omit(HU_pca_data[ind==2,])

library(psych)

pairs.panels(training[,-1],
             gap = 0,
             #bg = c("red", "yellow")[training$abb],
             pch=21)
pairs(training[, -1])

HU_pca <- prcomp(training[, -1],
                 center = TRUE,
                 scale. = TRUE)
attributes(HU_pca)
print(HU_pca)

summary(HU_pca)

pairs.panels(HU_pca$x,
             gap=0,
             bg = c("red", "yellow")[training$abb],
             pch=21)

library(ggbiplot)

g <- ggbiplot(HU_pca,
              obs.scale = 1,
              var.scale = 1,
              groups = training$Location,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)

trg <- predict(HU_pca, training)
trg <- data.frame(trg, training[1])
tst <- predict(HU_pca, testing)
tst <- data.frame(tst, testing[1])

library(nnet)

trg$abb <- relevel(trg$abb, ref = "HUN")
mymodel <- multinom(Location ~ PC1+PC2, data = trg)
summary(mymodel)

p <- predict(mymodel, trg)
tab <- table(p, trg$Location)
tab

1 - sum(diag(tab))/sum(tab)

p1 <- predict(mymodel, tst)
tab1 <- table(p1, tst$Location)
tab1

1 - sum(diag(tab1))/sum(tab1)


#' ## HU GLMM
#' Work on GLMM. Perform following steps:\
#' 1. Check environmental variables correlation
#' 2. Identify family\
#' 3. run model\
#' 4. Compared model prediction with original data\

#+ HU GLMM


lm(Leafwidth ~ total_REI + depth + air_exposure + ave_temp, data = trait_env_HU) -> aa
summary(aa)

