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

#' This document present the 

#' # All sample localities

loc <- read.csv("../../Sample points/Sample_points_include_Tropwater.csv")

loc %>%
  select(lon, lat, name) -> loc

#' # Environmental covariates

#' ## Relative physical exposure (fetch)

# read coastline data

t1 <- rast("../Coastal_data/ITEM_REL_24_147.05_-18.95.tif")
t2 <- rast("../Coastal_data/ITEM_REL_195_146.51_-18.95.tif")

as.contour(t2) -> t2c

st_as_sf(t2c) -> t2sf
# get world coastline polygon

maps::map("worldHires", plot = FALSE, fill = TRUE) %>%
  st_as_sf(crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") -> whr
sf_aus <- ozmap("states")
st_crs(sf_aus) <- st_crs(sf_point)

sf_point <- st_as_sf(x = loc,
                     coords = c("lon", "lat"),
                     crs = st_crs(t1))

fetch_len(sf_point[31, ], bearing = c(seq(0, 360, 30)), shoreline = t2sf, dmax = 40000)


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

