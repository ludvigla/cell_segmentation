#load the package
library(wholebrain)
library(rgl)
library(magick)
library(EBImage)

# Set this system variable to avoid RStudio from crashing
Sys.setenv(LIBGL_ALWAYS_SOFTWARE=1)
setwd("~/Michael_Ratz/cell_segmentation/hematoxylin_and_eosin_test/")
source("glassbrain_custom.R")

# coordinate AP for the image we are going to analyze
coord <- -1.65

# read images
# Ludvig: 
#   The dapi_adited image was processed in PS by adding two layers of brightness set to 150 and then I 
#   erased the background of the tissue (including the left hemisphere). Note that I exported the 
#   image without any rotations or scaling, I decided to do this from R instead. I just think
#   it's easier to keep track of whatever transformations were applied so that I can apply the same
#   for the spots later.
imdapi <- image_read("images/dapi_edited_extended.tif")

# Ludvig: 
#   Here I have rescaled the edited dapi tif image to 50% and rotated it by 90deg anti-clockwise
imdapi <- image_scale(imdapi, paste0(image_info(imdapi)$width*0.5)) %>% image_rotate(degrees = -90)
image_write(imdapi, path = "images/dapi_edited_extended_small_canvassize.tif")

# this is the filter we used it is a list object with the following parameters.
# Ludvig: 
#   I guess we can play around with these filters, but maybe they don't matter so much now that we just want to 
#   register the Visium spots.
myfilter <- list(alim = c(1, 50), #area limits for what to consider as cell bodies c(minimum, maximum)
               threshold.range = c(10000, 60000), #threshold range in the minimum and maximum in fluorescent intensity where the algorithm should search for neurons.
               eccentricity = 1000, #eccentricity (elongation) of contours sets how round you want cell bodies to be. Default is 1000 and smaller values equal to more round.
               Max = 60000, #Maximum value to display in the 8-bit rendered (sets sort of brightness contrast)
               Min = 0, #Minimum value to display in the 8-bit rendered (sets sort of brightness contrast)
               brain.threshold = 2800, #the exact value where you want to start segmeting the brain outline in autofluorescence.
               resize = 0.28, #resize parameter to match the atlas to your pixel resolution, should be between 0.03 and 0.2 for most applications.
               blur = 13, #blur parameter that sets the smoothness of the tissue outline, if magaded or jagged edges increase. Using a value fo 4 is usually recommended.
               downsample = 0.5 #downsample, default is set to 0.25 and images with a size of 15000 x 8000 pixels can then usually be run smoothly
)

# segmentation
# Ludvig: 
#   segmentation looks aweful, but not sure if it matters now.
imdapi_mod <- "images/dapi_edited_extended_small_canvassize.tif"
quartz()
seg <- segment(imdapi_mod, filter = myfilter, display = TRUE)

# to get the filter from a segmentation list into code.
# edit(seg$filter)

#get the fluorescent intensity of each DAPI nuclei as a gray-scale value
fluorescent.intensity <- scales::rescale(seg$soma$intensity)

#plot the segmentation
plot(seg$soma, pch = 16, cex = 0.25, asp = 1, 
     ylim = rev(range(seg$soma$y)), col = gray(fluorescent.intensity), 
     axes = F, ylab = '', xlab = '')


# run first pass at rgeistration
# Ludvig: 
#   I'm super confused by the resize parameter, but it seems to work with 0.4 :-/
quartz()
regi <- registration(imdapi_mod, filter = myfilter, coordinate = coord)

# Ludvig: 
#   these are coordinates I made, but they are not nearly as good as Daniels :-p
#   I hope you will be able to run this so that you can do it more properly!
regi$correspondance <- readRDS("correspondance.points.extended")

# Ludvig: 
#   This part if what I ran to defines the correspondance points loaded above
# ----------------------------------------
# remove all default 32 corr points
# regi <- remove.corrpoints(regi, 1:32)
# here I manually ad my own points to get it as accurate as I want. Close the command when done by rigt click
# regi <- add.corrpoints(regi)
# ----------------------------------------

# update the registration
regi <- registration(imdapi_mod, filter = myfilter, coordinate = coord, correspondance = regi)

# Ludvig: 
#   You can export the correspondance points as a data.frame. It's probably good to keep!
# Save correspondance points 
# saveRDS(regi$correspondance, file = "correspondance.points.extended")

# dataset with all the cells registered.
dataset <- get.cell.ids(regi, seg, forward.warp = TRUE)
dataset$animal <- 'mouse001'

# make some plots
jpeg(filename = "registration_results_extended.jpeg", width = 1200, height = 400)
par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))
plot.registration(regi, border = FALSE)
plot.registration(regi, border = 'orange', draw.trans.grid = TRUE)
plot.registration(regi, border = 'orange')
points(dataset$x, dataset$y, col = dataset$color, pch = 16, cex = 0.5)
dev.off()
#quartz.save('registration_results_extended.pdf', type = 'pdf', width = 12, height = 4)
pdf(file = "registration_results_extended.pdf", width = dim(im)[2]/100, height = dim(im)[1]/100)
plot.registration(regi, border = 'orange')
dev.off()


#stereotactic plot
schematic.plot(dataset, cex = 0.35, col = F)
quartz.save('cells_in_atlas_extended.png', type='png')

## ST functionality
# --------------------------------------------------------

# read cells pixel coordinates form external file
# cells <- read.table('cell_pixel_centroids.csv', sep=',', header =TRUE)
# get from dataset revious segmentation dataset here instead
st.object <- readRDS("../R_objects/st.object")
spots <- setNames(st.object@meta.data[st.object@meta.data$sample == "1", c("pixel_x", "pixel_y")], nm = c("x", "y"))

# Here we need to apply the same transformations that were used for the DAPI images
# First we need to define a couple of tranformation functions:

# rigid translation allows for shifts along x/y axes
rigid.transl <- function(h = 0, k = 0) {tr <-  matrix(c(1, 0, 0, 0, 1, 0, h, k, 1), nrow = 3); return(tr)}

# rigid transformation allows for rotations around a specified center
rigid.transf <- function(h = 0, k = 0, alpha = 0) {tr <- matrix(c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, h, k, 1), nrow = 3); return(tr)}

# The rotate function centers the x/y coords at origo, applies the rotation and then transforms back the x/y 
# coords to their original center
rotate <- function(angle, center.cur) {
  alpha <- 2*pi*(angle/360)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.transf(center.cur[1], center.cur[2], alpha)%*%tr
  return(tr)
}

# Takes a 3x3 transformation matrix and returns a mapping function
generate.map.affine <- function (
  tr, 
  forward = FALSE
) {
  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(solve(tr)%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      xy <- t(tr%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
}

# Apply rotation of 90 deg anti-clockwise around image center
imd <- readImage(imdapi_mod)
tr <- rotate(-90, center.cur = c(2232, 2232)) # The images were downscaled by a factor 0f 0.5 so the new image center is 4464/2, 4464/2
map.forward <- generate.map.affine(tr, forward = T)
xy <- as.data.frame(do.call(cbind, map.forward(x = spots$x, y = spots$y)))

# Check transformation results
# Now we can see that the coordinates have been transformed -90 deg
par(mfrow = c(1, 2))
plot(spots$x, 4464 - spots$y, xlim = c(0, 4464), ylim = c(0, 4464))
title(main = "original spots")
plot(xy[, 1], 4464 - xy[, 2], xlim = c(0, 4464), ylim = c(0, 4464))
title(main = "transformed spots")

# Now we just need to apply the same scaling factor that was used to downscale the 
# DAPI image earlier when creating "dapi_edited_small_canvassize.tif"
xy <- xy*0.5
xy$x <- xy$x + dim(imd)[1] - 2232 # Here we need to add the height of the added canvas that was added in PS

# We can check if the transformation was correct, i.e. the spots have 
# rotated -90deg and scaled by a factor of 0.5
par(mfrow = c(1, 1))
dapi_mod <- readImage(files = imdapi_mod)
EBImage::display(dapi_mod, method = "raster")
points(xy[, 1], xy[, 2], col = "red")

# Define pixel size and radius
image.size.micron <- 8705 
image.size.pixel <- 4464/2
pixels.per.um <- (image.size.pixel/image.size.micron)
pixelsize <- 1/pixels.per.um # micrometer
spot.radius <- 27.5 # micrometer

# Get x/y coordinates from transformed data.frame
x <- xy$x
y <- xy$y

# function to get polygon for spot
circle.perimeter <- function ( 
  x, y, 
  id, 
  spot.radius, 
  pixelsize
){
  points.per.circle <- 40
  theta = seq(0, 2*pi,length = points.per.circle)
  x = x + (spot.radius*cos(theta))/pixelsize
  y = y + (spot.radius*sin(theta))/pixelsize
  return(data.frame(x, y, id))
}

# plotting polygons
plot.polygon <- function (
  dataframe, 
  col = rgb(0.9, 0, 0.1, 0.7), 
  border = 'black'
){
  lapply(unique(dataframe$id), function(x){
    polygon(dataframe$x[dataframe$id == x], dataframe$y[dataframe$id == x], col = col, border = border)
  })
}

# function for counting cells inside spot
count.spot.inside <- function (
  seg, 
  dataset
){
  cell.count <- integer()
  for(i in unique(seg$soma$contour.ID)){
    cell.count <- c(cell.count, 
                    sum(point.in.polygon(dataset$x, dataset$y, 
                                         seg$soma$contour.x[seg$soma$contour.ID == i],
                                         seg$soma$contour.y[seg$soma$contour.ID == i])))
  }
  return(cell.count)
}

# make polygon dataframe
spots.contours <- do.call("rbind", lapply(seq_along(x), function(i){circle.perimeter(x[i], y[i], i, spot.radius, pixelsize)}))

#normalization function for normalizing cell counts
normalize <- function(x){
  x <- (x - min(x))
  x <- x/max(x)
  return(x)
}

#create a soma list object to contain the spot info in
spotSoma <- list(x = x,
                 y = y,
                 intensity <- rep(100, length(x)),
                 area = rep(pi*(spot.radius^2)*pixelsize, length(x)),
                 contour.x = spots.contours$x,
                 contour.y = spots.contours$y,
                 contour.ID = spots.contours$id)

# create segmentation output list object             
seg.Spots <- list(filter = myfilter, soma = spotSoma)

# plot the previous registration
quartz()
plot.registration(regi, border = 'orange')
# plot exact polygons from contours
plot.polygon(spots.contours, border = "red")

# get identity of features in atlas and in both coordinate systems
# datasetSpots <- get.cell.ids(regi, seg.Spots, forward.warp = TRUE) #use this for no plotting and faster
quartz()
datasetSpots <- inspect.registration(regi, seg.Spots, forward.warps  = TRUE)
datasetSpots$cell.count <- count.spot.inside(seg.Spots, xy)

#plot registration
quartz()
par(mfrow = c(1,3))
plot.registration(regi, border = 'orange', draw.trans.grid =TRUE)
# plot points with gray level intensity according to how many cells in each
points(datasetSpots[, c("x", "y")], pch = 21, col = 'pink', bg = gray(datasetSpots$cell.count))#gray(normalize(datasetSpots$cell.count)))
schematic.plot(datasetSpots,  mm.grid = FALSE, scale.bar = TRUE, device = FALSE)
#overwrite rgeion color with cell count gray scale
datasetSpots$color <- gray(datasetSpots$cell.count)#gray(normalize(datasetSpots$cell.count))
schematic.plot(datasetSpots,  mm.grid = FALSE, scale.bar = TRUE, device = FALSE)

#revert back
datasetSpots$color <- color.from.acronym(datasetSpots$acronym)

# make webmap
# Ludvig:
#   couldn't make this work ...
# makewebmap(imdapi_mod, myfilter, regi, dataset = datasetSpots, scale = 3)
save.image()


# Add info to st.object meta.data
gg <- subset(st.object@meta.data, sample == "1")
gg$acronym <- datasetSpots$acronym
gg$color <- datasetSpots$color
gg$color <- ifelse(is.na(gg$color), "white", gg$color)

library(magrittr)
library(dplyr)
gg %>% group_by(acronym) %>% summarize(mx = median(warped_x), my = median(warped_y))

# Plot results
ggplot(gg, aes(warped_x, warped_y)) +
  geom_point(fill = gg$color, shape = 21, size = 2, stroke = 0.2) +
  coord_flip() +
  theme_void() 
