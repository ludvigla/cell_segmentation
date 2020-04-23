#read cells pixel coordinates form external file
#cells <- read.table('cell_pixel_centroids.csv', sep=',', header =TRUE)
#get from dataset revious segmentation dataset here instead
cells <- dataset

#spots <- read.table('spots_pixel_centroids.csv', sep=',', header =TRUE)
#simulate spots instead
width <- 37
height <- 35
spots <- expand.grid(x = seq(60, 1140, length.out = width),
                     y = seq(160, 1200, length.out = height))

pixelsize <- 3 #micrometer
spot.radius <- 25 #micrometer

x <- spots$x
y <- spots$y

#function to get polygon for spot
circle.perimeter <- function(x, y, id, spot.radius, pixelsize){
  points.per.circle <- 40
  theta = seq(0,2*pi,length = points.per.circle)
  x = x + (spot.radius*cos(theta))/pixelsize
  y = y + (spot.radius*sin(theta))/pixelsize
  return(data.frame(x,y,id))
}

#plotting polygons
plot.polygon<-function(dataframe, col = rgb(0.9,0,0.1,0.7), border = 'black'){
  lapply(unique(dataframe$id), function(x){polygon(dataframe$x[dataframe$id==x], dataframe$y[dataframe$id==x], col=col, border = border)})
}

#function for counting cells inside spot
count.spot.inside <- function(seg, dataset){
  cell.count<-integer()
  for(i in unique(seg$soma$contour.ID)){
    cell.count<-c(cell.count, sum(point.in.polygon(dataset$x, 
                                  dataset$y, 
                                  seg$soma$contour.x[seg$soma$contour.ID == i],
                                  seg$soma$contour.y[seg$soma$contour.ID == i]))
    )
  }
  return(cell.count)
}

#make polydon dataframe
spots.contours<-do.call("rbind", 
                        lapply(seq_along(x), 
                               function(i){
                                 circle.perimeter(x[i], y[i], i, spot.radius, pixelsize)
                                 })
                        )

#normalization function for normalizing cell counts
normalize<-function(x){
  x<-(x-min(x))
  x<-x/max(x)
  return(x)
}


#create a soma list object to contain the spot info in
spotSoma <- list(x = x,
             y = y,
             intensity <- rep(100, length(x)),
             area = rep(pi*spot.radius^2*pixelsize, length(x)),
             contour.x = spots.contours$x,
             contour.y = spots.contours$y,
             contour.ID = spots.contours$id)

#create segmentation output list object             
seg.Spots<-list(filter = myfilter, soma = spotSoma)

#plot the previous registration
plot.registration(regi, border='orange')
# plot exact polygons from contours
plot.polygon(spots.contours)

#get identity of features in atlas and in both coordinate systems
#datasetSpots<-get.cell.ids(regi, seg.Spots, forward.warp = TRUE) #use this for no plotting and faster
datasetSpots<-inspect.registration(regi, seg.Spots, forward.warps  = TRUE)

datasetSpots$cell.count <- count.spot.inside(seg.Spots, cells)

#plot registration
par(mfrow=c(1,3))
plot.registration(regi, border='orange', draw.trans.grid =TRUE)
# plot points with gray level intensity according to how many cells in each
points(datasetSpots[,c("x", "y")], pch=21, col='pink', bg=gray(normalize(datasetSpots$cell.count)))
schematic.plot(datasetSpots,  mm.grid = FALSE, scale.bar = TRUE, device =FALSE)
#overwrite rgeion color with cell count gray scale
datasetSpots$color <- gray(normalize(datasetSpots$cell.count))
schematic.plot(datasetSpots,  mm.grid = FALSE, scale.bar = TRUE, device =FALSE)

#revert back
datasetSpots$color <- color.from.acronym(datasetSpots$acronym)

#make webmap
makewebmap(images[1], myfilter, regi, dataset = datasetSpots, scale = 3)
