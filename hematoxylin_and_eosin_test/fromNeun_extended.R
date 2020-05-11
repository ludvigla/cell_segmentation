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
imneun <- image_read("images/Neun_new.tif")

# Ludvig: 
#   Here I have rescaled the edited dapi tif image to 50% and rotated it by 90deg anti-clockwise
imneun <- image_scale(imneun, paste0(image_info(imneun)$width*0.5)) #%>% image_rotate(degrees = -90)
image_write(imneun, path = "images/Neun_new_small_canvassize.tif")

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
               resize = 0.5, #0.28, resize parameter to match the atlas to your pixel resolution, should be between 0.03 and 0.2 for most applications.
               blur = 13, #blur parameter that sets the smoothness of the tissue outline, if magaded or jagged edges increase. Using a value fo 4 is usually recommended.
               downsample = 0.5 #downsample, default is set to 0.25 and images with a size of 15000 x 8000 pixels can then usually be run smoothly
)

# segmentation
# Ludvig: 
#   segmentation looks aweful, but not sure if it matters now.
imneun_mod <- "images/Neun_new_small_canvassize.tif"
quartz()
seg <- segment(imneun_mod, filter = myfilter, display = TRUE)

# to get the filter from a segmentation list into code.
# edit(seg$filter)

#get the fluorescent intensity of each DAPI nuclei as a gray-scale value
fluorescent.intensity <- scales::rescale(seg$soma$intensity)

#plot the segmentation
plot(seg$soma, pch = 16, cex = 0.25, asp = 1, 
     ylim = rev(range(seg$soma$y)), col = gray(fluorescent.intensity), 
     axes = F, ylab = '', xlab = '')


# run first pass at rgeistration
quartz()
regi <- registration(imneun_mod, filter = myfilter, coordinate = coord)

# Ludvig: 
#   these are coordinates I made, but they are not nearly as good as Daniels :-p
#   I hope you will be able to run this so that you can do it more properly!
#regi$correspondance <- readRDS("correspondance.points.extended")

# Ludvig: 
#   This part if what I ran to defines the correspondance points loaded above
# ----------------------------------------
# remove all default 32 corr points
# regi <- remove.corrpoints(regi, 1:32)
# here I manually ad my own points to get it as accurate as I want. Close the command when done by rigt click
# regi <- add.corrpoints(regi)
# ----------------------------------------

# update the registration
#regi <- registration(imneun_mod, filter = myfilter, coordinate = coord, correspondance = regi)

# Ludvig: 
#   You can export the correspondance points as a data.frame. It's probably good to keep!
# Save correspondance points 
# saveRDS(regi$correspondance, file = "correspondance.points.extended")

# Save registration
saveRDS(object = regi, file = "../R_objects/registration_Neun_V9")


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
st.object.v13tov16.rotated <- readRDS("../R_objects/st.object.v13tov16.rotated")

# Select sample
# Now we should be able to select any of the V9-V12 samples (here called 1-4) since they are alreasdy aligned
spots.list <- lapply(1:4, function(s) {
  setNames(st.object.v13tov16.rotated@meta.data[st.object.v13tov16.rotated@meta.data$sample == paste0(s), c("warped_x", "warped_y")], nm = c("x", "y"))
})

# Apply rotation of 90 deg anti-clockwise around image center
ima <- magick::image_read("images/Neun_raw.tif") %>% magick::image_rotate(-90) %>% imager::magick2cimg()
imd <- magick::image_read(imneun_mod) %>% imager::magick2cimg()
#tr <- rotate(-90, center.cur = c(2232, 2232)) # The images were downscaled by a factor 0f 0.5 so the new image center is 4464/2, 4464/2
#map.forward <- generate.map.affine(tr, forward = T)
#xy <- as.data.frame(do.call(cbind, map.forward(x = spots$x, y = spots$y)))

xmax <- dim(ima)[2]
ymax <- dim(ima)[1]
xy.list <- lapply(spots.list, function(xy) {
  rs <- ((dim(imd)[2]/0.5)/xmax)
  xy <- xy*rs
  xy$x <- xy$x + dim(imd)[1]/0.5 - xmax*rs
  return(xy)
})

# Define pixel size and radius
image.size.micron <- 8705 
image.size.pixel <- dim(ima)[1]/2
pixels.per.um <- (image.size.pixel/image.size.micron)
pixelsize <- 1/pixels.per.um # micrometer
spot.radius <- 27.5 # micrometer

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
spots.contours.list <- lapply(xy.list, function(xy) {
  x <- xy$x*0.5
  y <- xy$y*0.5
  do.call("rbind", lapply(seq_along(x), function(i){circle.perimeter(x[i], y[i], i, spot.radius, pixelsize)}))
})

#normalization function for normalizing cell counts
normalize <- function(x){
  x <- (x - min(x))
  x <- x/max(x)
  return(x)
}

#create a soma list object to contain the spot info in
spotSoma.list <- lapply(seq_along(spots.contours.list), function(i) {
  spots.contours <- spots.contours.list[[i]]
  xy <- xy.list[[i]]
  x <- xy$x*0.5
  y <- xy$y*0.5
  list(x = x,
       y = y,
       intensity <- rep(100, length(x)),
       area = rep(pi*(spot.radius^2)*pixelsize, length(x)),
       contour.x = spots.contours$x,
       contour.y = spots.contours$y,
       contour.ID = spots.contours$id)
})

# create segmentation output list object             
seg.Spots.list <- lapply(spotSoma.list, function(spotSoma) {
  list(filter = myfilter, soma = spotSoma)
})

# plot the previous registration
quartz()
par(mfrow=c(1, 4))
plot.registration(regi, border = 'orange')
plot.polygon(spots.contours.list[[1]], border = "red")
plot.registration(regi, border = 'orange')
plot.polygon(spots.contours.list[[2]], border = "red")
plot.registration(regi, border = 'orange')
plot.polygon(spots.contours.list[[3]], border = "red")
plot.registration(regi, border = 'orange')
plot.polygon(spots.contours.list[[4]], border = "red")

# get identity of features in atlas and in both coordinate systems
# datasetSpots <- get.cell.ids(regi, seg.Spots, forward.warp = TRUE) #use this for no plotting and faster
quartz()
datasetSpots.list <- lapply(seq_along(seg.Spots.list), function(i) {
  seg.Spots <- seg.Spots.list[[i]] 
  xy <- xy.list[[i]]*0.5
  datasetSpots <- inspect.registration(regi, seg.Spots, forward.warps = TRUE)
  datasetSpots$cell.count <- count.spot.inside(seg.Spots, xy)
  return(datasetSpots)
})

#plot registration
quartz()
par(mfrow = c(1,3))
plot.registration(regi, border = 'orange', draw.trans.grid =TRUE)
# plot points with gray level intensity according to how many cells in each
points(datasetSpots.list[[1]][, c("x", "y")], pch = 21, col = 'pink', bg = gray(datasetSpots.list[[1]]$cell.count))#gray(normalize(datasetSpots$cell.count)))
schematic.plot(datasetSpots.list[[1]],  mm.grid = FALSE, scale.bar = TRUE, device = FALSE)
#overwrite rgeion color with cell count gray scale
datasetSpots.list[[1]]$color <- gray(datasetSpots.list[[1]]$cell.count)#gray(normalize(datasetSpots$cell.count))
schematic.plot(datasetSpots.list[[1]],  mm.grid = FALSE, scale.bar = TRUE, device = FALSE)

#revert back
datasetSpots.list <- lapply(datasetSpots.list, function(datasetSpots) {
  datasetSpots$color <- color.from.acronym(datasetSpots$acronym)
  return(datasetSpots)
})

# make webmap
# Ludvig:
#   couldn't make this work ...
# makewebmap(imneun_mod, myfilter, regi, dataset = datasetSpots, scale = 3)
save.image()


# Add info to st.object meta.data
gg <- st.object.v13tov16.rotated@meta.data
datasetSpots <- do.call(rbind, datasetSpots.list)
gg$acronym <- datasetSpots$acronym
gg$color <- datasetSpots$color
gg$color <- ifelse(is.na(gg$color), "white", gg$color)

library(magrittr)
library(dplyr)
gg %>% group_by(acronym) %>% summarize(mx = median(warped_x), my = median(warped_y))

# Plot results
ggplot(gg, aes(warped_x, 2000 - warped_y)) +
  geom_point(fill = gg$color, shape = 21, size = 2, stroke = 0.2) +
  theme_void() +
  facet_wrap(~sample)


# Plot outlined reference
plot.registration.only <- function (
  registration, 
  border = rgb(154, 73, 109, maxColorValue = 255)
) {
  scale.factor <- mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height, registration$transformationgrid$width))
  
  img <- paste(registration$outputfile, "_undistorted.png", sep = "")
  img <- readPNG(img)
  par(mar = c(0, 0, 0, 0))
  plot(NULL, ylim = c(dim(img)[1], 0), 
       xlim = c(0, dim(img)[2]), asp = 1, axes = F, xlab = "", ylab = "", 
       col = 0)
  
  if (!is.null(border)) {
    lapply(1:registration$atlas$numRegions, function(x) {
      polygon(registration$atlas$outlines[[x]]$xrT/scale.factor, 
              registration$atlas$outlines[[x]]$yrT/scale.factor, 
              border = border)
    })
    lapply(1:registration$atlas$numRegions, function(x) {
      polygon(registration$atlas$outlines[[x]]$xlT/scale.factor, 
              registration$atlas$outlines[[x]]$ylT/scale.factor, 
              border = border)
    })
  }
}


plot.registration.only(regi, border = "black")

