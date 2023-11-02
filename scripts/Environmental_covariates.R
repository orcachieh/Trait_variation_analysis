library(marmap)

qld <- getNOAA.bathy(lon1 = 137, lon2 = 156, lat1 = -10, lat2 = -27, resolution = 5)

summary(qld)

plot(qld, land = TRUE)
plot(qld, deep = 0, shallow = 0, step = 0,
     lwd = 0.4, add = TRUE)

qld[1]

values(qld)
