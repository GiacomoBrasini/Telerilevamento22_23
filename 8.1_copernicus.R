# COPERNICUS: downloading and visualizing data of Earth Observation

#install.packages("ncdf4")
# Load the required packages
library(ncdf4)
library(raster)
library(ggplot2)
library(viridis)

# Set the working directory in Windows
setwd("C:/lab")

# register and download data from:
# https://land.copernicus.vgt.vito.be/PDF/portal/Application.html


# Upload the data as a raster object
soil_moist <- raster("c_gls_SSM1km_202305090000_CEURO_S1CSAR_V1.2.1.nc")
soil_moist
ncell(soil_moist)  # the total size
plot(soil_moist)  # plot

# Coerce data into a dataframe
soil_moist_d <- as.data.frame(soil_moist, xy = T)

# Columns' names of the dataframe
colnames(soil_moist_d)

# Plot with ggplot
ggplot() +
  geom_raster(soil_moist_d,
              mapping = aes(x = x, y = y, fill = Surface.Soil.Moisture)) +
  ggtitle("Soil Moisture from Copernicus")
# fill = third column of the dataframe

# Let's crop the image by coordinates
ext <- c(23, 30, 62, 68) #declare the crop coordinates: x min, x max, y min and y max
soil_moist_crop <- crop(soil_moist, ext)
soil_moist_crop

# Coerce the cropped image into a dataframe
soil_moist_crop_d <- as.data.frame(soil_moist_crop, xy= T)

# Plot the cropped image with ggplot, adding viridis scale color
ggplot() +
  geom_raster(soil_moist_crop_d,
              mapping = aes(x = x, y = y, fill = Surface.Soil.Moisture)) +
  ggtitle("Cropped soil moisture from Copernicus") +
  scale_fill_viridis()


# Exercise: the same path but with an other image chosen by me
forest_cover <- raster("c_gls_FCOVER_202006030000_GLOBE_PROBAV_V1.5.1.nc")
forest_cover
ncell(forest_cover)
plot(forest_cover)

# Too big, let's crop
ext <- c(166, 179, -48, -34)  # New Zealand coordinates
forest_nz <- crop(forest_cover, ext)
forest_nz_d <- as.data.frame(forest_nz, xy = T)
colnames(forest_nz_d)

ggplot() +
  geom_raster(forest_nz_d,
              mapping = aes(x = x, y = y,
                            fill = Fraction.of.green.Vegetation.Cover.1km))
# fill = third column of the dataframe

# Let's focus on Egmnot National Park
ext <- c(173, 176, -40, -38 )
forest_eg <- crop(forest_nz, ext)
plot(forest_eg)
forest_eg_d <- as.data.frame(forest_eg, xy = T)

# Plot using ggplot
ggplot() +
  geom_raster(forest_eg_d,
              mapping = aes(x = x, y = y,
                            fill = Fraction.of.green.Vegetation.Cover.1km)) +
  ggtitle("Cropped forest cover from Copernicus") +
  scale_fill_viridis()

