# Time series analysis

# Greenland increase of temperature ----
# Data and code from Emanuela Cosma

library(raster)
library(rasterVis)
library(rgdal)

setwd("C:/lab/greenland")

# raster() function to import tif images with single layer
lst_2000 <- raster("lst_2000.tif")
lst_2005 <- raster("lst_2005.tif")
lst_2010 <- raster("lst_2010.tif")
lst_2015 <- raster("lst_2015.tif")

# rast() function for using the terra package

ls()  # list of objects imported and loaded

# Plot in a multiframe
par(mfrow = c(2,2))
plot(lst_2000)
plot(lst_2005)
plot(lst_2010)
plot(lst_2015)
# white corresponds to the presence of snow or ice

#Import the whole set altogheter
# List of files:
rlist <- list.files(pattern = "lst")
rlist

# Apply a function over a list or vector
import <- lapply(rlist, raster)  # to apply a function to many files of a list
import

# Stack vectors from a dataframe or list
TGr <- stack(import)  # or terra::c within terra package
TGr  # 4 layers
# stacking vectors concatenates multiple vectors into a single vector

plot(TGr)  # now we have all the images in a single element

cl <- colorRampPalette(c("blue","lightblue","pink","red"))(100)
plot(TGr, col = cl)

dev.off()

plot(TGr$lst_2000, col = cl)
#or 
#plot(TGr[[1]], col = cl)

plotRGB(TGr, 1, 2, 3, stretch = "Lin")
plotRGB(TGr, 2, 3, 4, stretch = "Lin")

dev.off()

# Difference between 2005 and 2000:
dift = TGr[[2]] - TGr[[1]]
cl <- colorRampPalette(c("blue","lightblue","pink","red"))(100)
plot(dift, col = cl)


#NO2 decrease during the lockdown period ----

setwd("C:/lab/en")

# Import the first image file
en01 <- raster("EN_0001.png")
en01
plot(en01) # it's January, before the lockdown

cl <- colorRampPalette(c('red','orange','yellow'))(100) #
plot(en01, col = cl)

en13 <- raster("EN_0013.png")
plot(en13, col = cl) # it's March, after lockdown began

# Let's import the whole set
rlist_2 <- list.files(pattern = "EN")
rlist_2

# lapply(X, FUN)
rimp <- lapply(rlist_2, raster)
rimp

# stack
en <- stack(rimp)
en

# plot everything
plot(en, col = cl)


# Plot first and last images
par(mfrow = c(1,2))
plot(en[[1]], col = cl)
plot(en[[13]], col = cl)
# or by using stack()
en113 <- stack(en[[1]], en[[13]])
plot(en113, col = cl)

# Let's see the difference:
difen <-  en[[1]] - en[[13]]
cl_dif <- colorRampPalette(c('blue','white','red'))(100) 
plot(difen, col = cl_dif)

# plotRGB
par(mfrow = c(2,1))
plotRGB(en, 1, 7, 13, stretch="lin")
plotRGB(en, 1, 7, 13, stretch="hist")
