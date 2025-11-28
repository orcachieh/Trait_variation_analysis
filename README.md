# Intraspecifc Trait Variation in *Halodule uninervis*

This repository contains the required data and analyses for intraspecific trait variation of *Halodule uninervis* in Tropical Queensland, Australia.


## Data
- trait_raw_data.csv contains the trait measurements of field collected samples. For more detail, see Seagrass TraitDB.
- Sediment Data contains the sediment grain size distribution collected from the field.
- environment_variables.csv contains the sorted environmental variables from Climate Data Store.
- Climate Data contains the original data downloaded from Climate Data Store.
- Seagrass_sample_points.csv contains the coordinates of samples.

## Scripts
- ITV_env.R contains the main analysis code.
- Visualize_raw_data.R provides the first steps of visualize traits distribution of all collected sample and produce the plot for supplementary figure S1.
- Environmental_variables.R is the process of extracting environmental factors base on sample location.
- Functions.R contains customize utility functions.

## Stan_file and models
- This folders contain the stan code for mixture regression and output model to avoid the extra work required for setting up cmddstan and C++ in R

## Analysis structure
The analyses includes following steps
1. *H. uninervis* trait variation subgrouping (Visualizing distribution and PCA)
2. The definition of each growth form (mixture model)
3. Links between environmental variables of trait variation (GLMMs)
4. Estimating effect of ITV on ecosystem services provided by seagrass meadows (Simulation)