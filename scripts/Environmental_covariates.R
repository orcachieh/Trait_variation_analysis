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
#' In order to accomplish the problem of some intertidal point are higher than sea level and is "on land" if using normal coastal line, I use the intertidal elevation layer (Intertidal Extents Model 25 m v 2.0.0) from Geoscience Australia and specify shore line as level 7 instead of level 0.This mean the shore line is actually 7 meter elevation (model value) to calculate fetch.
#' Using the consistent method for all sample points should be alright for calculate relative exposure among the sites.

#' Due  to large geographic range. The shore line was separate to following location for easier fetch calculation - 
#' to: Pallarenda + Maggi + Cleveland 
#' cl: Clairview
#' gl: Gladstone
#' he: Hervey bay
#' ka: karumba
#' we: Weipa

# read intertidal elevation layer (Intertidal Extents Model 25 m v 2.0.0)
t1 <- rast("../Coastal_data/ITEM_REL_24_147.05_-18.95.tif")
t2 <- rast("../Coastal_data/ITEM_REL_195_146.51_-18.95.tif")

# merge several raster to obtain proper range if required
merge(t1, t2) -> to

# create contour line from raster file
as.contour(to) -> toc

# transfer terra object to sf to compatible with "waver" package
st_as_sf(toc)-> tosf

# remove sediment sample
# project sample points as coastline data, GDA. This process is required for fecth calculation by "waver" package
loc_point_gda <- st_transform(loc_point, crs(tosf))

# calculate fetch for sample points. 12 fetch lines (360 degree each with 30 degree separation) are calculated for each point.
# maximum fetch is 20000 m if no coastline present at that angle

fetch_len_multi(loc_point_gda[1:47, ], bearing = c(seq(0, 330, 30)), shoreline =tosf[tosf$level == 7,], dmax = 20000) -> to_fetch

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
# plotting fetch and coastline

ggplot() +
  geom_sf(data = tosf[tosf$level == 7,]) +
  geom_spoke(data = to_spoke, 
             aes(x = X, y = Y, angle = angle, radius = fetch)) +
  geom_sf(data = loc_point_gda[1:47, ], color = "red") +
  coord_sf(xlim = c(min(st_coordinates(loc_point_gda[1:47, ])[, 1]) - 20000,
                    max(st_coordinates(loc_point_gda[1:47, ])[, 1]) + 20000),
           ylim = c(min(st_coordinates(loc_point_gda[1:47, ])[, 2]) - 20000,
                    max(st_coordinates(loc_point_gda[1:47, ])[, 2]) + 20000))


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

