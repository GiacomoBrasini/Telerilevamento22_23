# LiDAR: 3D in R

# LiDAR is a method for determining ranges by targeting an object or a surface
# with a laser and measuring the time for the reflected light to return to the receiver.
# It may scan multiple directions -> 3D

# Load the required packages
library(raster)  # Geographic Data Analysis and Modeling
library(rgdal)  # Geospatial Data Abstraction Library
library(viridis)
library(ggplot2)
library(patchwork)

# Set working directory in Windows
setwd("C:/lab/data")


# Load data: Digital Surface Model 2013
dsm_2013 <- raster("2013Elevation_DigitalElevationModel-0.5m.tif")
dsm_2013  # info of the raster object

# Plot the DSM 2013
plot(dsm_2013, main = "Lidar Digital Surface Model San Genesio/Jenesien")

# Load data: Digital Terrain Model 2013
dtm_2013 <- raster("2013Elevation_DigitalTerrainModel-0.5m.tif")
dtm_2013

# Plot the DTM 2013
plot(dtm_2013, main = "Lidar Digital Terrain Model San Genesio/Jenesien")


# Create CHM 2013 as difference between DEM and DTM
chm_2013 <- dsm_2013 - dtm_2013
chm_2013  # view CHM attributes

# Coerce the CHM 2013 into a dataframe
chm_2013d <- as.data.frame(chm_2013, xy = T)

# Plot CHM 2013 using ggplot
ggplot() +
  geom_raster(chm_2013d,
              mapping = aes(x = x, y=y, fill = layer)) +
  scale_fill_viridis() +
  ggtitle("CHM 2013 San Genesio/Jenesien")
# dark blue represents the meadows around the village,
# while the lighter colors represent the heights of the plants

# Save the CHM on computer
writeRaster(chm_2013,"chm_2013_San_genesio.tif")


# Load data: Digital Surface Model 2004
dsm_2004 <- raster("2004Elevation_DigitalElevationModel-2.5m.tif")
dsm_2004

# Plot the DSM 2004
plot(dsm_2004, main="Lidar Digital Surface Model San Genesio/Jenesien")

# Load data: Digital Terrain Model 2004
dtm_2004 <- raster("2004Elevation_DigitalTerrainModel-2.5m.tif")

# Plot the DTM 2004
plot(dtm_2004, main = "Lidar Digital Terrain Model San Genesio/Jenesien")


# Create CHM 2004 as difference between DSM and DTM
chm_2004 <- dsm_2004 - dtm_2004
chm_2004

# Coerce into a dataframe
chm_2004d <- as.data.frame(chm_2004, xy=T)

# Plot CHM 2004
ggplot() +
  geom_raster(chm_2004d, mapping =aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis() +
  ggtitle("CHM 2004 San Genesio/Jenesien")

# Save CHM on computer
writeRaster(chm_2004,"chm_2004_San_genesio.tif")


#error: rasters have different resolution
difference_chm<-chm_2013-chm_2004

#reseample 2013 to 2004 @2.5m
chm_2013_reseampled<-resample(chm_2013, chm_2004)

#calculate difference in CHM
difference<-chm_2013_reseampled-chm_2004

# Coerce into a dataframe
difference_d <- as.data.frame(difference, xy=T)

#plot the difference
ggplot() +
  geom_raster(difference_d, mapping =aes(x=x, y=y, fill=layer)) +
  scale_fill_viridis() +
  ggtitle("difference CHM San Genesio/Jenesien")
# darker areas indicate that there has been logging over the years,
# lighter areas indicate an increase  


#save the rasters
writeRaster(chm_2013_reseampled,"chm_2013_reseampled_San_genesio.tif")
writeRaster(difference,"difference chm San_genesio.tif")


## point cloud
#install.packages("lidR")
library(lidR)

#load point cloud
point_cloud<-readLAS("C:/lab/point_cloud.laz")

#plot r3 point cloud
plot(point_cloud)
