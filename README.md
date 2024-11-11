# Intraspecifc Trait Variation in *Halodule uninervis*

This repository contains the required data and analyses for intraspecific trait variation of *Halodule uninervis* in Tropical Queensland, Australia.
The corresponded publication can be found though this link: 

## Data
- trait_raw_data.csv contains the trait measurements of field collected samples. For more detail, see Seagrass TraitDB.
- Sediment Data contains the sediment grain size distribution collected from the field.
- environment_variables.csv contains the sorted environmental variables from Climate Data Store.
- Climate Data contains the original data downloaded from Climate Data Store.
- Seagrass_sample_points.csv contains the coordinates of samples.

## Scripts
- ITV_env.R contains the main analysis code.
- Visualize_raw_data.R provides the first steps of visulize traits distribtuion of all collected sample and porduce the plot for supplementary fiugre S1.
- Environmental_variables.R is the process of extracting environmetal factors base on sample location.
- Functions.R contains customize utility functions.

## Analysis structure
The analyses includes following steps
1. *H. uninervis* trait variation subgrouping
2. The definition of each growth form
3. Links between environmental varialbes of trait variation
4. Estimating effect of ITV on ecosystem services provided by seagrass meadows