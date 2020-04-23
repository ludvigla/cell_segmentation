#load the package
library(wholebrain)
library(rgl)
Sys.setenv(LIBGL_ALWAYS_SOFTWARE=1)
setwd("~/Michael_Ratz/hematoxylin_and_eosin/")
source("glassbrain_custom.R")

#coordinate AP for the image we are going to analyze
coord <- -1.65

#system path for folder where images can be found
folder <- 'images'

#get only image files with full paths in the folder (this is used as input for registration and segmentation)
images <- get.images(folder)

#this is the filter we used it is a list object with the following parameters.
myfilter <- list(alim = c(1, 50), #area limits for what to consider as cell bodies c(minimum, maximum)
               threshold.range = c(20000, 65536), #threshold range in the minimum and maximum in fluorescent intensity where the algorithm should search for neurons.
               eccentricity = 1000, #eccentricity (elongation) of contours sets how round you want cell bodies to be. Default is 1000 and smaller values equal to more round.
               Max = 60000, #Maximum value to display in the 8-bit rendered (sets sort of brightness contrast)
               Min = 0, #Minimum value to display in the 8-bit rendered (sets sort of brightness contrast)
               brain.threshold = 2800, #the exact value where you want to start segmeting the brain outline in autofluorescence.
               resize = 0.28, #resize parameter to match the atlas to your pixel resolution, should be between 0.03 and 0.2 for most applications.
               blur = 13, #blur parameter that sets the smoothness of the tissue outline, if magaded or jagged edges increase. Using a value fo 4 is usually recommended.
               downsample = 1.5 #downsample, default is set to 0.25 and images with a size of 15000 x 8000 pixels can then usually be run smoothly
)

#segmentation
quartz()
seg <- segment(images[1], filter = myfilter, display = TRUE)

#to get the filter from a segmentation list into code.
edit(seg$filter)

#normalization function for fluorescent intensity
normalize <- function(x){
  x <- (x - min(x))
  x <- x/max(x)
  return(x)
}

#get the fluorescent intensity of each DAPI nuclei as a gray-scale value
fluorescent.intensity <- normalize(seg$soma$intensity)

#plot the segmentation
plot(seg$soma, pch = 16, cex = 0.25, asp = 1, 
     ylim = rev(range(seg$soma$y)), col = gray(fluorescent.intensity), 
     axes = F, ylab = '', xlab = '')


#run first pass at rgeistration. For cryosections with ST the brain is normally quite deformed so usually I just remove original corr points
quartz()
regi <- registration(images[1], filter = myfilter, coordinate = coord, right.hemisphere =  TRUE)

#these are coordinates I made.
regi$correspondance <- data.frame(
  targetP.x = c(
    106.8,
    102.35,
    341.26,
    700.68,
    496.76,
    540.46,
    554.61,
    583.51,
    222.42,
    349.76,
    496.54,
    361.26,
    199.73,
    1065.95,
    627.4,
    805.82,
    296.98,
    371.86,
    112.46,
    655.46,
    1009.41,
    823.18,
    868.2,
    864.21,
    477.32,
    114.07,
    253.16,
    514.98,
    744.62,
    1066.33,
    397.12,
    394.25,
    456.73,
    257.66,
    149.37,
    151.1,
    883.42,
    934.19,
    953.2,
    1016.81,
    1045.61,
    981.81,
    195.31,
    102.65,
    102.21,
    178.12,
    484.13,
    1179.55,
    1157.53,
    947.12,
    815.54,
    767.49,
    784.47,
    743.18,
    735.92
  ),
  targetP.y = c(
    509.8,
    1474.06,
    1585.34,
    1657.16,
    749.11,
    777.19,
    733,
    546.57,
    623.67,
    595.68,
    522.75,
    435.06,
    567.27,
    1335.16,
    598.88,
    781.46,
    1495.49,
    392.68,
    429.36,
    1572.26,
    1308.18,
    1501.35,
    1558.39,
    1205.96,
    631.55,
    184.16,
    86.82,
    452.15,
    230.49,
    584.3,
    122.58,
    168.33,
    248.07,
    365.4,
    403.63,
    227.2,
    695.79,
    662.9,
    871.08,
    856.98,
    1118.1,
    1107.39,
    1501.56,
    1191.66,
    636.15,
    710.2,
    1663.53,
    1031.07,
    822.19,
    368.97,
    472.27,
    539.38,
    1625.22,
    1540.33,
    1337.34
  ),
  referenceP.x = c(
    197.83,
    192.95,
    492.44,
    678.39,
    562.17,
    601.84,
    622.02,
    598.54,
    300.47,
    409.71,
    508.05,
    394.04,
    250.64,
    921.76,
    647.79,
    754.34,
    435.74,
    392.54,
    199.86,
    692.89,
    886.76,
    800.61,
    838.17,
    792.09,
    519.47,
    195.22,
    293.54,
    518.81,
    685.95,
    933.74,
    383.93,
    386.74,
    444.11,
    298.81,
    226.36,
    230.37,
    809.17,
    849.5,
    850.98,
    905.44,
    924.96,
    864.92,
    310.52,
    195.86,
    200.24,
    270.35,
    583.67,
    1012.96,
    1001.55,
    824.76,
    757.57,
    726.05,
    785.04,
    759.16,
    639.79
  ),
  referenceP.y = c(
    508.36,
    1356.74,
    1267.28,
    1386.12,
    677.58,
    687.6,
    655.57,
    529.23,
    639.86,
    599.42,
    530.36,
    459.32,
    567.8,
    1031.72,
    559.12,
    663.1,
    1202.28,
    413.18,
    437.26,
    1341.18,
    1016.63,
    1200.94,
    1220.28,
    914.82,
    615.31,
    273.26,
    217.9,
    465.13,
    280.54,
    498.39,
    214.23,
    234.58,
    276.72,
    383.37,
    408.79,
    305.02,
    592.53,
    566.57,
    699.3,
    679.56,
    859.79,
    850.54,
    1349.15,
    1040.68,
    662.83,
    706.82,
    1341.74,
    791.51,
    650.21,
    365.33,
    444.49,
    497.94,
    1316.73,
    1286.22,
    1104.61
  ),
  shape = c(
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4,
    4
  )
)

#remove all default 32 corr points
#regi <- remove.corrpoints(regi, 1:32)
#here I manually ad my own points to get it as accurate as I want. Close the command when done by rigt click
#regi <- add.corrpoints(regi)

#update the registration
regi <- registration(images[1], filter = myfilter, coordinate = coord, right.hemisphere =  TRUE, correspondance = regi)

#dataset with all the cells registered.
dataset <- get.cell.ids(regi, seg, forward.warp = TRUE)
dataset$animal <- 'mouse001'

#make some plots
par(mfrow = c(1, 3))
plot.registration(regi, border=FALSE)
plot.registration(regi, border='orange', draw.trans.grid = TRUE)
plot.registration(regi, border='orange')
points(dataset$x, dataset$y, col = dataset$color, pch=16, cex=0.2)
quartz.save('registration_results.pdf', type='pdf')

#stereotactic plot
schematic.plot(dataset, cex=0.15, col = F)
quartz.save('cells_in_atlas.pdf', type='pdf')

getwd()
save.image()

#plot dataset in 3D
glassbrain_custom(dataset)
#drag around to a point of view you like and save the perspective: pp<-par3d(no.readonly = TRUE)

#set persepective on the plot.
pp <- structure(list(FOV = 30, ignoreExtent = FALSE, listeners = 1L, 
                   mouseMode = structure(c("trackball", "zoom", "fov", "pull"
                   ), .Names = c("left", "right", "middle", "wheel")), skipRedraw = FALSE, 
                   userMatrix = structure(c(0.603139936923981, 0.43978950381279, 
                                            -0.665437698364258, 0, -0.0221172589808702, -0.824720859527588, 
                                            -0.565106868743896, 0, -0.79732871055603, 0.355556219816208, 
                                            -0.487695276737213, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), scale = c(1, 
                                                                                                             1, 1), viewport = structure(c(0L, 0L, 1280L, 720L), .Names = c("x", 
                                                                                                                                                                            "y", "width", "height")), zoom = 1, windowRect = c(0L, 45L, 
                                                                                                                                                                                                                               1280L, 765L), family = "sans", font = 1L, cex = 1, useFreeType = TRUE), .Names = c("FOV", 
                                                                                                                                                                                                                                                                                                                  "ignoreExtent", "listeners", "mouseMode", "skipRedraw", "userMatrix", 
                                                                                                                                                                                                                                                                                                                  "scale", "viewport", "zoom", "windowRect", "family", "font", 
                                                                                                                                                                                                                                                                                                                  "cex", "useFreeType"))
#set perspective
par3d(pp)

#update plot but with high resolution 3D mesh
glassbrain(dataset, device = FALSE, high.res = TRUE)

#save 3D plot as PNG
rgl.snapshot('3Dbrain.png')


## ST functionality
# --------------------------------------------------------
# Yes original H&E and Cy3 would be good that way you will also get automatic cell counts per spot/brain region etc.

# All ST like functionality is here:
#   https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fgithub.com%2Ftractatus%2Fwholebrain%2Fblob%2Fmaster%2FR%2Fstaoseq.R&amp;data=02%7C01%7Cmichael.ratz%40ki.se%7C7b990f8d2fcc472ca49708d7d26aea0a%7Cbff7eef1cf4b4f32be3da1dda043c05d%7C0%7C0%7C637209228956089411&amp;sdata=gyPGjc0lquE2maa%2FQ3G62voL0ylEOUo25YmzVIedhN4%3D&amp;reserved=0

# Basically you have rg2gray() which will convert the original H&E and Cy3 for you (for H&E you want to invert it).

# You need to set the hNe.filter and spot.filter once, I can help you with that if you have original images.

# For example this is how beginning of a script could look like:
  
source("staoseq.R")

#get all H&E images and convert them.
hNe <- rgb2gray(images)

#get all Cy3 images and convert them.
images <- get.images(getwd(), type = 'jpg')
spots <- rgb2gray(images, invert = FALSE)

#segment the spots
seg.spots <- segment(spots[1], filter = spot.filter, get.contour = TRUE)
seg.nuclei <- segment(hNe[1], filter = hNe.filter)

#register to H&E only:
regi <- registration(hNe[1], coordinate = -1, filter = hNe.filter, right.hemisphere = TRUE)
regi <- remove.corrpoints(regi, 1:32)
regi <- add.corrpoints(regi)
regi <- registration(hNe[1], coordinate -1, filter= hNe.filter, right.hemisphere = TRUE, correspondance = regi)

#combine TSV file with count matrix. Here corner is the variable that tells you where the small square of dots are 1 is top left, 2 is top right etc….

stdata <- "genecountmatrix.tsv"
dataset <- combine.st.data(seg.spots, regi, stdata, seg.nuclei, corner = 1)


#plot 3D brain

glassbrain(dataset$spots)
