#### sort_nc ####
# Sort the NetCF data from climate store and store the data in "terra SpatRaster" for future application. 
# Plot the first 2 and last 2 time stamp as example.

sort_nc <- function(File_name){
  nc_data <- nc_open(paste("../Climate_data/", File_name, ".nc",sep = ""))
  # Save the print(nc) dump to a text file
  {
  sink(paste("../Climate_data/", File_name, ".txt",sep = ""))
  print(nc_data)
  sink()
  }
  
  nc_close(nc_data)
  
  rast_terra <- rast(paste("../Climate_data/", File_name, ".nc",sep = ""))
  # According to sample period (Last sample was collected on 2023 Jan 24), I subset data from 2015 Feb to 2023 Jan (8 years, 96 months in total)
  rast_terra <- subset(rast_terra, 3:98)
  
  # Plot first 2 and last 2 layers for reference
  plot(subset(rast_terra, c(1:2, 95:96)))
  
  return(rast_terra)
}

#### match_fetch ####
# This function find the fetch length that match with the direction of prevalence wind each month at each location
# Return a tibble contains prevalence wind speed, direction, and matched fetch length of each month at each location
# the wind_data need to have a column called mean_dir and contain the direction of the (prevalence) wind

match_fetch <- function(wind_data, fetch_data){
  # list the available fetch direction (every 5 degree in this case)
  fetch_direction <- (as.numeric(colnames(fetch_data)))
  
  # Find the closest available fetch angle to the prevalence wind and store it
  sapply(wind_data$mean_dir, 
         function(i) fetch_direction[which.min(abs(i - fetch_direction))]) %>%
    matrix(nrow = nrow(fetch_data), ncol = 12, byrow = TRUE) -> closest_angle
  
  # Then use the closest_anlge to select the fetch for each month at each location
  select_fetch <- matrix(nrow = nrow(fetch_data), ncol = 12)
  for(i in 1:nrow(fetch_data)){
    for(j in 1:12){
      select_fetch[i, j] <- all_fetch[i, colnames(all_fetch) == closest_angle[i, j]]
    }
  }
  
  # combind the selected fetch length and wind data
  select_fetch %>%
    as_tibble() %>%
    pivot_longer(everything(), values_to = "fetch") %>%
    cbind(wind_data) -> wind_fetch
  
  return(wind_fetch)
}

#### Zoom in fetch plot ####
zoom <- function(plot, sites){
  plot + 
    coord_sf(xlim = c(min(st_coordinates(loc_point_gda[sites, ])[, 1]) - 20000,
                      max(st_coordinates(loc_point_gda[sites, ])[, 1]) + 20000),
             ylim = c(min(st_coordinates(loc_point_gda[sites, ])[, 2]) - 20000,
                      max(st_coordinates(loc_point_gda[sites, ])[, 2]) + 20000))
}

#### Depth ####
# This function use the contour lines from NIDEM data set to extract data of each sample points.
# In order to reduce computing pressure, the depth is extract for each region.
# Import "loc_point" for sample points and "temp_2m" for projection as default input

sample_depth <- function(contour_data, ref = temp_2m, loc = loc_point){
  
  # Create an empty raster layer to store elevation value
  conrast <- rast(xmin = xmin(contour_data),
                  xmax = xmax(contour_data),
                  ymin = ymin(contour_data),
                  ymax = ymax(contour_data),
                  crs = crs(contour_data))
  # Create a point layer with elevation information of contour line
  conp <- as.data.frame(as.points(contour_data), geom = "XY")
  
  # Use "gstat" to create a model for interpolate
  gmodel <- gstat(id="elev_m", formula = elev_m~1, 
             locations = ~x+y, data=conp,
             nmax=7, set=list(idp = .5))
  
  # Interpolate the elevation for each cell in the raster layer
  z <- interpolate(conrast, gmodel, debug.level = 0, index = 1)
  
  # Project the elevation layer to temp_2m for extract purpose.
  # And plot the layer as reference
  # method = "near" use nearest neighborhood makes sure the raster value is not average across nearby cells
  z_proj <- project(z, crs(ref), method = "near")
  plot(z_proj)
  
  # Extract depth of all samples points. The points not cover in the raster extent will receive NaN
  sample_d <- extract(z_proj, loc)
  
  return(sample_d)
}

#### Air exposure ####
# This function collect the relative exposure index for each sample site at each subset location (eaier for computer)
# Import "loc_point" for sample points and "temp_2m" for projection as default input

sample_air <- function(air_exposure, ref = temp_2m, loc = loc_point){
  
  # project air exposure layer to temp_2m for extract purpose.
  # method = "near" use nearest neighborhood makes sure the raster value is not average across nearby cells
  project(air_exposure, crs(ref), method = "near") -> air_exposure_proj
  
  # Extract depth of all samples points. The points not cover in the raster extent will receive NaN
  sample_a <- extract(air_exposure_proj, loc_point)
  
  return(sample_a)
}

#### Sediment ####
# This function matches the sediment data to sample location by nearest neighborhood
# Then add the sediment data to the seagrass sample point (loc_point)
# Finally, assign NA to subtidal sample locations since their are no sediment samples available for them

sed_match <- function(seagrass, sediment){
  # To use terra::nearest, both layer need to be SpatVector type
  loc_sed_near <- vect(seagrass)
  sed_near <- vect(sediment, crs = crs(loc_sed_near))
  
  # Run nearest to find cloesest matched seagrass and sediment points
  nearest(loc_sed_near, sed_near) -> sea_sed
  # create a tibble to store the matched sediment sample ID
  tibble(`ID` = values(sea_sed)$to_id) -> to_id
  
  # Create a ID column in sediment for joint purpose
  sediment %>%
    mutate(ID = 1:16) -> sed_ID
  
  # Use the matched sediment ID to create the sediment data for each sample points
  to_id %>% 
    left_join(sed_ID, by = "ID") %>%
    select(Gravel:Clay) %>%
    cbind(location = loc_point$name) -> sed_point
  
  # Assign NA to intertidal points. This need to be caution
  sed_point[grep("*SUB*", sed_point$location), 1:4] = NA
  
  return(sed_point)
}


sed_mean_match <- function(seagrass, sediment){
  # To use terra::nearest, both layer need to be SpatVector type
  loc_sed_near <- vect(seagrass)
  sed_near <- vect(sediment, crs = crs(loc_sed_near))
  
  # Run nearest to find cloesest matched seagrass and sediment points
  nearest(loc_sed_near, sed_near) -> sea_sed
  # create a tibble to store the matched sediment sample ID
  tibble(`ID` = values(sea_sed)$to_id) -> to_id
  
  # Create a ID column in sediment for joint purpose
  sediment %>%
    mutate(ID = 1:16) -> sed_ID
  
  # Use the matched sediment ID to create the sediment data for each sample points
  to_id %>% 
    left_join(sed_ID, by = "ID") %>%
    select(sed_mean) %>%
    cbind(location = loc_point$name) -> sed_point
  
  # Assign NA to intertidal points. This need to be caution
  sed_point[grep("*SUB*", sed_point$location), 1] = NA
  
  return(sed_point)
}


#### PCA biplot ####
ITV_pca_biplot <- function(ITV_pca, cluster, color_arrow){
  if(color_arrow == 1){
    fviz_pca_biplot(ITV_pca, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 2.5,
                    fill.ind = as.factor(cluster),
                    col.ind = "black",
                    # Color variable by groups
                    col.var = "contrib",
                    gradient.cols = "RdYlBu",
                    legend.title = list(fill = "Clusters", color = "Contrib."),
                    repel = TRUE,        # Avoid label overplotting
                    title = paste("Biplot", str_sub(deparse(substitute(cluster)), start = 14), "clusters", sep = " ")
    ) +
      ggpubr::fill_palette("jco")   # Individual fill color
  } else {
    fviz_pca_biplot(ITV_pca, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 2.5,
                    fill.ind = as.factor(cluster),
                    palette = c("#009E73", "#E69F00"), # Individual fill color palette
                    col.ind = "gray7",
                    col.var = "gray20",
                    legend.title = "Clusters",
                    repel = TRUE,        # Avoid label overplotting
                    title = paste("Biplot", str_sub(deparse(substitute(cluster)), start = 14), "clusters", sep = " ")
    )
  }
}

