#' ---
#'  title: "Environmental variables"
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

library(marmap)
library(tidyverse)
library(shape)
library(waver) # Calculate fetch
library(rWind) # Calculate wind
library(mapdata)
library(maps)
library(ncdf4) # package for netcdf manipulation
library(terra)
library(ggplot2)
library(doBy)
library(gstat) # For interpolate depth from contour lines

#+ Source functions from Functions.R 
#### Import functions  ####

source("./Functions.R")

#' This document present the all environmental data used in chapter 2


#' # All sample points
#+ Read sample points
#### Read all sample localities ####

# read and cleanup sample point csv
loc <- read.csv("../../Sample points/Sample_points_include_Tropwater.csv")

loc %>%
  select(lon, lat, name) -> loc
loc

# create a sf object and assign correct coordinate reference system (CRS) -> WGS84 (it's a records from GPS. Generally it should be WGS84
loc_point_all <- st_as_sf(x = loc,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#' # Environmental variables\
#' Required data are listed below\
#' 1. Relative physical exposure\
#' 2. Depth\
#' 3. Air exposure\
#' 4. Average temperature\
#' 5. Average of minimum temperature\
#' 6. Average of maximum temperature\
#' 7. Temperature variation\
#' 8. Sediment\
#' 9. Runoff\
#' 10. Precipitation\
#' 

#' ## Import climate data from climate store ####

#' The import data are monthly average of "ERA5 hourly data on single levels from 1940 to present"
#' The download data ranged from 2014 Dec to 2023 Nov. According to sample period (Last sample was collected on 2023 Jan 24), I subset data from 2015 Feb to 2023 Jan (8 years, 96 months in total)

#+ Read climate shop data
#### Environmental variables ####
eastward_10m_wind <- sort_nc("eastward_10m_wind_2015_2023")
northward_10m_wind <- sort_nc("northward_10m_wind_2015_2023")
temp_2m <- sort_nc("2m_temperature_2015_2023")
sea_temp <- sort_nc("Sea_surface_temperature_2015_2023")
runoff <- sort_nc("Runoff_2015_2023")
perci <- sort_nc("Total_precipitation_2015_2023")


#' ##  1. Relative physical exposure
#' There are three availalbe methods to reflect Relative physical exposure \
#' 1. Effective fetch (EF), which present fetch from several angles and one fetch is an averagre of nearby fetch\
#' 2. Relative exposure Index (REI) multiply the EF by wind speed (Masom et. al, 2018).\
#' 3. Wave energy (package "waver") is calculated by three component: fetch, dominant wind speed, and depth.\

#' Starting with fetch.\ 
#' In order to accomplish the problem of some intertidal point are higher than sea level and is "on land" if using normal coastal line, I use the intertidal elevation layer (National Intertidal Digital Elevation Model 25m 1.0.0, Bishop-Taylor et al, 2019) from Geoscience Australia and specify shore line to as close to 0 as possible (the contour line they created does not have exactly 0 m for all the subset map)\
#' Using the consistent method for all sample points should be alright for calculate relative exposure among the sites.\

#' Due  to large geographic range. The shore line are read in different to for potential operation later\
#' to: Pallarenda + Maggi + Cleveland\ 
#' cl: Clairview\
#' gl: Gladstone\
#' he: Hervey bay\
#' ka: karumba\
#' we: Weipa\

#+ Create fetch
#### Relative physical exposure ####

# read intertidal elevation contour layer (National Intertidal Digital Elevation Model 25m 1.0.0)
to1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_24_147.05_-18.95.shp") # contain sample points
to2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_195_146.51_-18.95.shp") # contain sample points

cl1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_191_149.62_-21.95.shp") # contain sample points
cl2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_143_149.79_-22.14.shp")
cl3 <- vect("../Coastal_data/NIDEM/NIDEM_contours_234_149.84_-22.32.shp")
cl4 <- vect("../Coastal_data/NIDEM/NIDEM_contours_39_149.58_-21.77.shp")

gl1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_269_151.42_-23.71.shp") # contain sample points
gl2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_177_151.69_-23.18.shp")
gl3 <- vect("../Coastal_data/NIDEM/NIDEM_contours_218_151.04_-23.36.shp")
gl4 <- vect("../Coastal_data/NIDEM/NIDEM_contours_253_152.38_-24.58.shp")
gl5 <- vect("../Coastal_data/NIDEM/NIDEM_contours_97_152.53_-23.99.shp")

he1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_300_152.96_-25.35.shp") # contain sample points
he2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_168_153.32_-25.48.shp")
he3 <- vect("../Coastal_data/NIDEM/NIDEM_contours_249_152.90_-25.01.shp")

we1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_244_141.39_-12.88.shp") # contain sample points
we2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_99_141.57_-12.19.shp")

ka1 <- vect("../Coastal_data/NIDEM/NIDEM_contours_136_140.66_-17.14.shp") # contain sample points
ka2 <- vect("../Coastal_data/NIDEM/NIDEM_contours_281_140.26_-17.52.shp")

# merge all subset maps

all <- rbind(to1, to2,
             cl1, cl2, cl3, cl4,
             gl1, gl2, gl3, gl4, gl5,
             he1, he2, he3,
             we1, we2,
             ka1, ka2)

# transfer terra object to sf to compatible with "waver" package
allsf <- st_as_sf(all)

# each subset map contain 9 contour lines. Grouping them by 9 and find the minimum positive number in each groups will give us a afetch_plotroximated coastline for fetch calculation
allsf %>%
  mutate(op = rep(1:(length(elev_m)/9), each = 9)) %>%
  group_by(op) %>%
  filter(elev_m > 0) %>%
  slice(which.min(elev_m)) -> coast_elev

# remove 1. sediment sample point
#        2. Seagrass boundary at Nelly bay
#        3. Samples I did not measure in the lab
loc_point_all %>%
  mutate(name = ifelse(row_number() == 26, "PAL9_26122022", name)) %>%
  filter(!grepl("SED", name),
         !grepl("SEAG", name),
         !grepl("PIC*_06112022", name),
         !name %in% c("PAL6_09102022", "PAL10_22012023", 
                      "PAL11_22012023", "PAL12_22012023",
                      "PIC3_10092022",
                      "GEO4_24012023", "GEO5_24012023", "GEO6_24012023")) -> loc_point

# project sample points as coastline data, GDA. This process is required for fecth calculation by "waver" package
loc_point_gda <- st_transform(loc_point, crs(allsf))

# calculate fetch for sample points. 72 fetch lines (360 degree each with 5 degree separation) are calculated for each point.
# maximum fetch is 20000 m if no coastline present at that angle
# shoreline is define as the contour line closest to 0 m elevation in each subset map

all_fetch <- waver::fetch_len_multi(loc_point_gda, bearing = c(seq(0, 330, 5)), shoreline = coast_elev, dmax = 20000)

# Rearrange fetch matrix for geom_spoke plotting
# First just keep the fetch at each 30 degree for visualization
# There is a bit of work due to angles are seen differently from two function.
# The fetch were calculated from north and clockwise, but geom_spoke plot from east (0 degree of x axis) and counter clockwise.

all_fetch %>%
  as_tibble() %>%
  select(seq(1, 67, by = 6)) %>%
  pivot_longer(everything(), names_to = "angle", values_to = "fetch") %>%
  mutate(angle = (abs(as.numeric(angle)-360)+90)*pi/180,
         # Re-assign angle and transform to radian
         X = rep(st_coordinates(loc_point_gda)[, 1], each = 12),
         Y = rep(st_coordinates(loc_point_gda)[, 2], each = 12)) -> all_spoke

# plotting fetch and coastline

#+ Fetch plot, echo = FALSE
# For all site
ggplot() +
  geom_sf(data = coast_elev) +
  geom_spoke(data = all_spoke, 
             aes(x = X, y = Y, angle = angle, radius = fetch)) +
  geom_sf(data = loc_point_gda, color = "red") +
  coord_sf(xlim = c(min(st_coordinates(loc_point_gda)[, 1]) - 20000,
                    max(st_coordinates(loc_point_gda)[, 1]) + 20000),
           ylim = c(min(st_coordinates(loc_point_gda)[, 2]) - 20000,
                    max(st_coordinates(loc_point_gda)[, 2]) + 20000)) -> fetch_plot

# Zoom to different location
# Townsville 1:31
zoom(fetch_plot, 1:31)

# Clairview (32)
zoom(fetch_plot, 32)

# Hervey Bay (33)
zoom(fetch_plot, 33)

# Karumba (34)
zoom(fetch_plot, 34)

# Gladstone (35) 
zoom(fetch_plot, 35)

# Weipa (36)
zoom(fetch_plot, 36)

#' Calculate wind intensity and direction through u and v component of 2m wind
#+ Calculate wind intensity and direction

u10 <- extract(eastward_10m_wind, loc_point)
u10[, -1] %>%
  pivot_longer(everything(), names_to = "time", values_to = "u10") %>%
  mutate(u10 = u10, .keep = "used") %>%
  unlist() -> u10_math

v10 <- extract(northward_10m_wind, loc_point)
v10[, -1] %>%
  pivot_longer(everything(), names_to = "time", values_to = "v10") %>%
  mutate(v10 = v10, .keep = "used") %>%
  unlist() -> v10_math

# Calculate wind speed and direction of each location and each time point
rWind::uv2ds(u10_math, v10_math) -> wind_dir_speed

# the "dir" from uv2ds is angle start from north and going clockwise.
# Same as fecth calculate from waver

# Calculate monthly average wind of each location
wind_10m_month <- as_tibble(wind_dir_speed) %>%
  mutate(location = as.factor(rep(1:36, each = 96)),
         each12 = as.factor(rep(rep((1:12), each = 8), times = 36))) %>%
  group_by(location, each12) %>%
  summarise(mean_dir = min(dir), mean_speed = mean(speed), n = n()) %>%
  ungroup() # This create a long table with prevalence wind speed and direction of each month at each sample location.

# Then multiply the wind to closest directional fetch (every 5 degree) and sum up to get relative exposure

# use match fetch to match the closest fetch to the monthly prevalence wind
wind_fetch <- match_fetch(wind_10m_month, all_fetch)

colnames(wind_fetch)

# Multiply the wind speed with fetch length and sum up 12 month to receive total REI of a year at each location
# a funny rearrange at the end to make this variable have the same format as others

wind_fetch %>%
  select(fetch, mean_dir, mean_speed) %>%
  mutate(REI = fetch * mean_speed, 
         location = rep(loc_point$name, each = 12)) %>%
  group_by(location) %>%
  summarise(total_REI = sum(REI)) %>%
  select(!location) %>%
  cbind(location = loc_point$name) -> total_REI

#' ## 2. Depth
#+ Depth
#### Depth ####

# Create Digital Elevation Models (DEMs) from the contour lines shape files above to receive depth of each sample points.

# Townsville region. Need to merge 2 contour lines file first to cover all points
to <- rbind(to1, to2)
tod <- sample_depth(to)

# Clairview
cld <- sample_depth(cl1)

# Hervey Bay
hed <- sample_depth(he1)

# Karumba
kad <- sample_depth(ka1)

# Gladstone
gld <- sample_depth(gl1)

# Weipa
wed <- sample_depth(we1)

# merge all data frames together and create final elevation dataframe
d_list <- list(tod, cld, hed, kad, gld, wed)

d_list %>% 
  reduce(full_join, by='ID') %>%
  select(-ID) %>%
  mutate(depth = na.omit(unlist(.)), .keep = "none") %>%
  cbind(location = loc_point$name) -> loc_elev

#' ## 3. Air exposure
#' Air exposure read the data from Intertidal Extent Model. It provides the percentage time that the location is expose to air
#' 0: no exposure\
#' 1: at least 0 to 10% time is expose, ...\
#' 9: highest, 80 - 100% time is expose)\

#+ Air exposure
#### Air exposure ####

to1_a <- rast("../Coastal_data/ITEM/ITEM_REL_24_147.05_-18.95.tif")
to2_a <- rast("../Coastal_data/ITEM/ITEM_REL_195_146.51_-18.95.tif")

cl1_a <- rast("../Coastal_data/ITEM/ITEM_REL_191_149.62_-21.95.tif")

gl1_a <- rast("../Coastal_data/ITEM/ITEM_REL_269_151.42_-23.71.tif")

he1_a<- rast("../Coastal_data/ITEM/ITEM_REL_300_152.96_-25.35.tif")

we1_a<- rast("../Coastal_data/ITEM/ITEM_REL_244_141.39_-12.88.tif")

ka1_a <- rast("../Coastal_data/ITEM/ITEM_REL_136_140.66_-17.14.tif")

# Townsville region. Need to merge 2 contour lines file first to cover all points
to_a <- merge(to1_a, to2_a)
to_air <-sample_air(to_a)

# Clairview
cl_air <-sample_air(cl1_a)

# Hervey Bay
he_air <-sample_air(he1_a)

# Karumba
ka_air <-sample_air(ka1_a)

# Gladstone
gl_air <-sample_air(gl1_a)

# Weipa
we_air <-sample_air(we1_a)

# merge all data frames together and create final elevation dataframe
a_list <- list(to_air, cl_air, he_air, ka_air, gl_air, we_air)

a_list %>% 
  reduce(full_join, by='ID') %>%
  select(-ID) %>%
  mutate(air_exposure = na.omit(unlist(.)), .keep = "none") %>%
  cbind(location = loc_point$name) -> air_exposure

#' ##  4. Average temperature

#+ Average temperature
#### Average temperature ####

temp <- extract(temp_2m, loc_point)
temp[,-1] %>%
  rowMeans() %>%  
  as_tibble() %>%
  mutate(ave_temp = value, .keep = "none") %>%
  cbind(location = loc_point$name) -> ave_temp

#' ##  5. Average minimum temperature

#+ Average minimum temperature
#### Average minimum temperature ####

temp[,-1] %>%
  pivot_longer(everything(), names_to = "time", values_to = "temp") %>%
  mutate(location = as.factor(rep(1:36, each = 96)),
         each12 = as.factor(rep(rep((1:12), each = 8), times = 36))) %>%
  group_by(location, each12) %>%
  summarise(min = min(temp), n = n()) %>%
  summarise(mean_min = mean(min), n = n()) %>%
  ungroup() %>%
  mutate(mean_min = mean_min, .keep = "used") %>%
  cbind(location = loc_point$name) -> ave_min_temp


#' ## 6. Average of maximum temperature

#+ Average maximum temperature
#### Average maximum temperature ####

temp[,-1] %>%
  pivot_longer(everything(), names_to = "time", values_to = "temp") %>%
  mutate(location = as.factor(rep(1:36, each = 96)),
         each12 = as.factor(rep(rep((1:12), each = 8), times = 36))) %>%
  group_by(location, each12) %>%
  summarise(max = max(temp), n = n()) %>%
  summarise(mean_max = mean(max), n = n()) %>%
  ungroup() %>%
  mutate(mean_max = mean_max, .keep = "used") %>%
  cbind(location = loc_point$name) -> ave_max_temp

#' ## 7. Average temperature variation
#' Temperature variation is define as the max - min temprature of each 12 month

#+ Average temperature variation
#### Average temperature variation ####

temp[,-1] %>%
  pivot_longer(everything(), names_to = "time", values_to = "temp") %>%
  mutate(location = as.factor(rep(1:36, each = 96)),
         each12 = as.factor(rep(rep((1:12), each = 8), times = 36))) %>%
  group_by(location, each12) %>%
  summarise(max = max(temp), min = min(temp), n = n()) %>%
  mutate(var = max - min) %>%
  summarise(mean_var = mean(var), n = n()) %>%
  ungroup() %>%
  mutate(mean_var = mean_var, .keep = "used") %>%
  cbind(location = loc_point$name) -> ave_temp_var

#' ## 8. Sediment

#+ Sediment
#### Sediment ####

sed <- read.csv("../Sediment_data/sed_four_categories.csv")
colnames(sed) <- c("name", "Gravel", "Sand", "Silt", "Clay")

loc %>%
  filter(grepl("SED", name)) %>%
  mutate(loc = str_split_i(name, "_", 1),
         date = str_split_i(name, "_", 2),
         loc = gsub("\\s|\\d", "", loc),
         date = str_c(str_sub(date, 1, 4), str_sub(date, 7, 8), sep = ""),
         name = str_c(loc, date, sep = "_")) %>%
  select(lat, lon, name)-> loc_sed

sed %>%
  left_join(loc_sed, by = "name") -> sed_partial

# The missing locality of sediment sample are
# 1. ZC -> use ZC081022 locality directly
# 2. GEOSED_061122 -> use GEO3_06112022 locality
# 3. NELSED_240123 -> use NEL3_24012023 locality
# 4. TropWATER samples (*5) -> use seagrass sample localities directly

# 1.
sed_partial[5, c(6, 7)] = loc[47, c(2, 1)]
#2
sed_partial[7, c(6, 7)] = loc[6, c(2, 1)]
#3
sed_partial[9, c(6, 7)] = loc[14, c(2, 1)]
#4
sed_partial[c(12:16), c(6, 7)] = loc[c(52, 49, 50, 48, 51), c(2, 1)]
# A section of ugly manual code. If reoreder any of the 2 tables. This section need to be changed accordingly.

# rename to state the file is complete
sed_all <- sed_partial

# Matching the sediment data to sample location by nearest neighborhood

# seagrass input: loc_point, sediment input: sed_all
loc_sed <- sed_match(loc_point, sed_all)

#' ## 9. Runoff

#+ Runoff
#### Runoff ####

roff <- extract(runoff, loc_point)
roff[,-1] %>%
  rowMeans() %>%
  as_tibble() %>%
  mutate(runoff = value, .keep = "none") %>%
  cbind(location = loc_point$name) -> ave_runoff

#' ## 10. Precipitation
#' The download data from climate store was resampled to get monthly sum precipitation.
#' Here we finalized it as average monthly precipitation

#+ Precipitation
#### Precipitation ####

perc <- extract(perci, loc_point)
perc[,-1] %>%
  rowMeans() %>%
  as_tibble() %>%
  mutate(perci = value*1000, .keep = "none") %>% # multiply 1000 to change the unit to mm
  cbind(location = loc_point$name) -> ave_perci

#' # Combine 10 variables

#+ Combine 10 variables
#### Combine 10 variables ####

var <- list(total_REI, loc_elev, air_exposure,
            ave_temp, ave_min_temp, ave_max_temp,
            ave_temp_var, loc_sed, ave_runoff, ave_perci)
# Sort and clean up to produce finally table of variables
var %>%
  reduce(full_join, by="location") %>%
  select(location, total_REI:perci) -> fin_var

fin_var

# write.csv(fin_var, "../environment_variables.csv", row.names = FALSE)
