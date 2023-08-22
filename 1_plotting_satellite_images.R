# Firts lesson of Remote Sensing in GitHub

# Let's install the raster package
library(raster)

# Set the working directory
setwd("C:/lab")

# Brick function to import the data
image_2011 <- brick("p224r63_2011_masked.grd")
image_2011  # ncell = total pixels
# a RasterBrick is a multi-layer raster object

# beware of the size of plots pane, if it's too small it won't work in RStudio

# Plot data in a raw way
plot(image_2011) # default legend

colors <- colorRampPalette(c("red", "orange", "yellow")) (100) # n° of shades
colors_3 <- colorRampPalette(c("red", "orange", "yellow")) (3)

# Plot with the colours palette we chose
plot(image_2011, col = colors)
plot(image_2011, col = colors_3) # less details = worse image quality

# Plot just one element of the multi-layer image
plot(image_2011[[4]], col = colors)

nir <- image_2011$B4_sre # another way to plot it, by assignment
plot(nir, col = colors)
plot(nir, col = colors_3) # less shades mean less details

# Change of palette
new_palette <- colorRampPalette(c("aquamarine", "cadetblue", "cyan")) (100)
plot(image_2011, col = new_palette)
plot(image_2011[[4]], col = new_palette) # the fourth layer

# Exercise: plot the NIR band
# b1 = blue
# b2 = green
# b3 = red
# b4 = NIR

# New colour palette
palette <- colorRampPalette(c("pink", "violet", "darkorchid4")) (100)

# Plot the NIR band of the image in different ways
plot(image_2011$B4_sre, col = new_palette)
plot(image_2011[[4]], col = palette) 

dev.off()  # to clean plots

# Function to export graphs
pdf("my_1st_graph.pdf")
plot(image_2011[[4]], col = palette) 
dev.off()

png("first.png")
plot(image_2011[[4]], col = palette) 
dev.off()


# Function to plot 2 or more graphs in a multiframe
par(mfrow = c(2, 1))  # a multiframe with 2 rows and 1 column
plot(image_2011[[3]], col = palette) #red
plot(image_2011[[4]], col = palette) #NIR

dev.off() # to remove the multiframe

# Here we can see the differences between 2 different shaded palette
par(mfrow = c(2, 1))
plot(nir, col = colors) # 100-shaded palette, more details
plot(nir, col = colors_3) # 3-shaded palette, less details

dev.off()


# Let's plot all bands/layers
par(mfrow = c(2, 2))

#blue
col_blue <- colorRampPalette(c("blue4", "blue2", "lightblue"))(100)
plot(image_2011[[1]], col = col_blue)

#green
col_green <- colorRampPalette(c("darkgreen", "green", "lightgreen"))(100)
plot(image_2011[[2]], col = col_green)

#aquamarine
col_aquamarine <- colorRampPalette(c("aquamarine4", "aquamarine2", "aquamarine"))(100)
plot(image_2011[[3]], col = col_aquamarine)

#pink
col_pink <- colorRampPalette(c("pink", "violet", "darkorchid4"))(100)
plot(image_2011[[4]], col = col_pink)


# RGB plotting, plot of multi-layered raster object
 # 3 bands are combined such that they represent the red, green and blue channel
 # This function can be used to make 'true (or false) color images
 # from Landsat and other multi-band satellite images.

# Blue = 1
# Green = 2
# Red = 3

plotRGB(image_2011, r = 3, g = 2, b = 1, stretch = "Lin") # natural colours
plotRGB(image_2011, r = 4, g = 3, b = 2, stretch = "Lin") # vegetation is red
plotRGB(image_2011, r = 3, g = 2, b = 4, stretch = "Lin") # vegetation is blue
plotRGB(image_2011, r = 3, g = 4, b = 2, stretch = "Lin") # vegetation is green
# the fourth band is the NIR band, plants strongly reflect NIR

# Multiframe with natural and false colours
par(mfrow = c(2, 1))
plotRGB(image_2011, r = 3, g = 2, b = 1, stretch = "Lin")  # natural colours
plotRGB(image_2011, r = 4, g = 3, b = 2, stretch = "Lin")

# Histogram stretching
plotRGB(image_2011, r = 3, g = 2, b = 1, stretch = "Hist")  # natural colours
plotRGB(image_2011, r = 4, g = 3, b = 2, stretch = "Hist")

# Mixed stretching, differences between the 2 types of stretch
plotRGB(image_2011, r = 4, g = 3, b = 2, stretch = "Lin")
plotRGB(image_2011, r = 4, g = 3, b = 2, stretch = "Hist")

dev.off()

# Exercise: plot the NIR band
plot(image_2011[[4]])

# Import the 1988 image
image_1988 <- brick("p224r63_1988_masked.grd")
plot(image_1988)

# Plot in RGB space (natural colours)
plotRGB(image_1988, r = 3, g = 2, b = 1, stretch = "Lin")
plotRGB(image_1988, r = 4, g = 3, b = 2, stretch = "Lin")

plotRGB(image_1988, 4, 3, 2, stretch = "Lin")  # faster way

# Multiframe to see the differences between 1988 and 2011
par(mfrow = c(2,1))
plotRGB(image_1988, 4, 3, 2, stretch = "Lin")
plotRGB(image_2011, 4, 3, 2, stretch = "Lin")

dev.off()
plotRGB(image_1988, 4, 3, 2, stretch = "Hist")  # histogram stretch

# Multiframe with 4 images
par(mfrow=c(2,2))
plotRGB(image_1988, 4, 3, 2, stretch="Lin")
plotRGB(image_2011, 4, 3, 2, stretch="Lin")
plotRGB(image_1988, 4, 3, 2, stretch="Hist")
plotRGB(image_2011, 4, 3, 2, stretch="Hist")
