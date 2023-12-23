# Sort the NetCF data from climate store and store the data in 3-dimensional array for future application. 
# Plot the first time stamp as example.

sort_nc <- function(File_name){
  nc_data <- nc_open(paste("../Climate_data/", File_name, ".nc",sep = ""))
  # Save the print(nc) dump to a text file
  {
  sink(paste("../Climate_data/", File_name, ".txt",sep = ""))
  print(nc_data)
  sink()
  }
  
  ndvi.array <- ncvar_get(nc_data) # store the data in a 3-dimensional array
  
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat", verbose = F)
  t <- ncvar_get(nc_data, "time")
  
  nc_close(nc_data)
  
  # Create a simple plot for first time stamp
  ndvi.slice <- ndvi.array[ , , 1]
  r <- rast(t(ndvi.slice),
            #xmin = min(lon), xmax = max(lon),
            #ymin = min(lat), ymax = max(lat),
            crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
  r <- flip(r, direction = "vertical")
  plot(r)
  
  return(ndvi.array)
}
