#' ---
#'  title: "Environmental covariates"
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
library(marmap)
library(tidyverse)
library(shape)
library(waver)
library(mapdata)
library(maps)
library(terra)
library(ggplot2)

#' This document present the 

#' # All sample localities
# read and cleanup sample point csv
loc <- read.csv("../../Sample points/Sample_points_include_Tropwater.csv")

loc %>%
  select(lon, lat, name, ele) -> loc

# create a sf object and assign correct coordinate reference system (CRS) -> WGS84 (it's a records from GPS. Generally it should be WGS84
loc_point <- st_as_sf(x = loc,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#' # Environmental covariates

#' ## Relative physical exposure (fetch)
#' There are three availalbe methods to reflect Relative physical exposure 
#' 1. Effective fetch (EF), which present fetch from several angles and one fetch is an averagre of nearby fetch
#' 2. Relative exposure Index (REI) multiply the EF by wind speed (Masom et. al, 2018).
#' 3. Wave energy (package "waver") is calculated by three component: fetch, dominant wind speed, and depth.

#' Starting with fetch. 
#' In order to accomplish the problem of some intertidal point are higher than sea level and is "on land" if using normal coastal line, I use the intertidal elevation layer (National Intertidal Digital Elevation Model 25m 1.0.0, Bishop-Taylor et al, 2019) from Geoscience Australia and specify shore line to as close to 0 as possible (the contour line they created does not have exactly 0 m for all the subset map)
#' Using the consistent method for all sample points should be alright for calculate relative exposure among the sites.

#' Due  to large geographic range. The shore line was separate to following location for easier fetch calculation - 
#' to: Pallarenda + Maggi + Cleveland 
#' cl: Clairview
#' gl: Gladstone
#' he: Hervey bay
#' ka: karumba
#' we: Weipa

# read intertidal elevation contour layer (National Intertidal Digital Elevation Model 25m 1.0.0)
to1 <- vect("../Coastal_data/NIDEM_contours_24_147.05_-18.95.shp")
to2 <- vect("../Coastal_data/NIDEM_contours_195_146.51_-18.95.shp")

cl1 <- vect("../Coastal_data/NIDEM_contours_191_149.62_-21.95.shp")
cl2 <- vect("../Coastal_data/NIDEM_contours_143_149.79_-22.14.shp")
cl3 <- vect("../Coastal_data/NIDEM_contours_234_149.84_-22.32.shp")
cl4 <- vect("../Coastal_data/NIDEM_contours_39_149.58_-21.77.shp")



# merge several subset maps to obtain proper range if required
to <- rbind(to1, to2)
cl <- rbind(cl1, cl2, cl3, cl4)

# transfer terra object to sf to compatible with "waver" package
tosf <- st_as_sf(to)
clsf <- st_as_sf(cl)

# remove sediment sample
# project sample points as coastline data, GDA. This process is required for fecth calculation by "waver" package
loc_point_gda <- st_transform(loc_point, crs(tosf))

# calculate fetch for sample points. 12 fetch lines (360 degree each with 30 degree separation) are calculated for each point.
# maximum fetch is 20000 m if no coastline present at that angle
# shoreline is define as the contour line closest to 0 m elevation in each subset map

to_fetch <- fetch_len_multi(loc_point_gda[1:47, ], bearing = c(seq(0, 330, 30)), shoreline = tosf[tosf$elev_m %in% c(0.04, 0.09) ,], dmax = 20000)

cl_fetch <- fetch_len(loc_point_gda[48, ], bearing = c(seq(0, 330, 30)), shoreline = clsf[clsf$elev_m %in% c(0.43, 0.39, 0.43) ,], dmax = 20000)

# Rearrange fetch matrix for geom_spoke plotting
# There is a bit of work due to angles are seen differently from two function.
# The fetch were calculated from north and clockwise, but geom_spoke plot from east (0 degree of x axis) and counter clockwise.

to_fetch %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "angle", values_to = "fetch") %>%
  mutate(angle = (abs(as.numeric(angle)-360)+90)*pi/180,
         # Re-assign angle and transform to radian
         X = rep(st_coordinates(loc_point_gda[1:47, ])[, 1], each = 12),
         Y = rep(st_coordinates(loc_point_gda[1:47, ])[, 2], each = 12)) -> to_spoke

# Single point fetch goes a bit differnet way to sort to plotting format
cl_fetch %>%
  as_tibble() %>%
  mutate(angle = (abs(seq(0, 330, 30)-360)+90)*pi/180,
         fetch = value,
         # Re-assign angle and transform to radian
         X = rep(st_coordinates(loc_point_gda[48, ])[, 1], each = 12),
         Y = rep(st_coordinates(loc_point_gda[48, ])[, 2], each = 12), .keep = "none") -> cl_spoke

# plotting fetch and coastline

ggplot() +
  geom_sf(data = tosf[tosf$elev_m %in% c(0.04, 0.09) ,]) +
  geom_spoke(data = to_spoke, 
             aes(x = X, y = Y, angle = angle, radius = fetch)) +
  geom_sf(data = loc_point_gda[1:47, ], color = "red") +
  coord_sf(xlim = c(min(st_coordinates(loc_point_gda[1:47, ])[, 1]) - 20000,
                    max(st_coordinates(loc_point_gda[1:47, ])[, 1]) + 20000),
           ylim = c(min(st_coordinates(loc_point_gda[1:47, ])[, 2]) - 20000,
                    max(st_coordinates(loc_point_gda[1:47, ])[, 2]) + 20000))

ggplot() +
  geom_sf(data = clsf[clsf$elev_m %in% c(0.43, 0.39, 0.43) ,]) +
  geom_spoke(data = cl_spoke, 
             aes(x = X, y = Y, angle = angle, radius = fetch)) +
  geom_sf(data = loc_point_gda[48, ], color = "red") +
  coord_sf(xlim = c(min(st_coordinates(loc_point_gda[48, ])[, 1]) - 20000,
                    max(st_coordinates(loc_point_gda[48, ])[, 1]) + 20000),
           ylim = c(min(st_coordinates(loc_point_gda[48, ])[, 2]) - 20000,
                    max(st_coordinates(loc_point_gda[48, ])[, 2]) + 20000))

# Bathymetry data obtain from NOAA (highest resolution 0.5 degree)
# Problem intertidal points are actually above higher than sea level

qldt <- getNOAA.bathy(lon1 = 146.5, lon2 = 148, lat1 = -19, lat2 = -19.5, resolution = 0.5)


qlds <- getNOAA.bathy(lon1 = 137, lon2 = 156, lat1 = -20, lat2 = -27, resolution = 0.5)
summary(qld)

blues <- colorRampPalette(c("purple","blue",
                            "cadetblue1","lightsteelblue1"))
plot(qldt, image = TRUE, land = TRUE,
     bpal = list(c(0, max(qldt), "grey"),
                 c(min(qldt), 0, blues(100))))
plot(qldt, deep = 0, shallow = 0, step = 0,
     lwd = 1.5, add = TRUE)
scaleBathy(qldt, deg = 1, x = "bottomleft", inset = 5)

#points(loc$lon, loc$lat, pch = 19, cex = 3, asp = 1)

loc %>%
  filter(lat < -19, lat > -20) -> loct

sp <- get.depth(qldt, loct[, 1:2], locator = FALSE)
sp

mx <- abs(min(sp$depth, na.rm = TRUE))
col.points <- femmecol(round(mx))

points(sp[,1:2], col = "black", bg = col.points[abs(sp$depth)],
       pch = 21, cex = 1.5)
colorlegend(zlim = c(mx, 0), col = rev(col.points),
            main = "depth (m)", posx = c(0.85, 0.88))

