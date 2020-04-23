data(glassbrain, envir = environment())
library(rgl)
Sys.setenv(LIBGL_ALWAYS_SOFTWARE = 1)

glassbrain_custom <- function (
  dataset, 
  high.res = FALSE, 
  dim = c(720, 1080), 
  device = TRUE, 
  col = "region", 
  cex = 0.5, 
  hemisphere = "right", 
  spheres = FALSE, 
  alpha = 1, 
  laterality = TRUE, 
  plane = "coronal"
) {
  if (sum(dataset$color == "#000000") > 0) {
    dataset <- dataset[-which(dataset$color == "#000000"), ]
  }
  if (device) {
    open3d(windowRect = c(0, 0, 1280, 720))
  }
  par3d(persp)
  if (high.res) {
    drawScene.rgl(list(VOLUME))
  } else {
    drawScene.rgl(list(VOLUMESMALL))
  }
  if (col == "region") {
    color <- as.character(dataset$color)
  } else {
    color <- col
  }
  if (hemisphere == "right") {
    hemisphere <- 2
  } else {
    hemisphere <- 1
  }
  if (plane == "coronal") {
    smp.AP <- rnorm(length(dataset$AP), 0, (320/9.75) * 0.2)
    smp.ML <- 0
  } else {
    smp.ML <- rnorm(length(dataset$AP), 0, (320/9.75) * 0.2)
    smp.AP <- 0
  }
  if (laterality) {
    if (length(unique(dataset$AP)) > 1) {
      laterality <- table(dataset$AP, dataset$right.hemisphere)
      for (i in 1:nrow(laterality)) {
        if (hemisphere == "right") {
          if (laterality[i, 1] > laterality[i, 2]) {
            dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])] <- -dataset$ML[which(dataset$AP == 
                                                                                                         as.numeric(row.names(laterality))[i])]
          }
        } else {
          if (laterality[i, 2] > laterality[i, 1]) {
            dataset$ML[which(dataset$AP == as.numeric(row.names(laterality))[i])] <- -dataset$ML[which(dataset$AP == 
                                                                                                         as.numeric(row.names(laterality))[i])]
          }
        }
      }
    }
  }
  if (spheres) {
    spheres3d(paxTOallen(dataset$AP) - 530/2 + smp.AP, -dataset$DV * 
                1000/25 - 320/2, dataset$ML * 1000/25 + smp.ML, col = color, 
              radius = cex, alpha = alpha)
  } else {
    points3d(paxTOallen(dataset$AP) - 530/2 + smp.AP, (-dataset$DV * 
                                                         1000/25 * 0.95) - 320/2, dataset$ML * 1000/25 + smp.ML, 
             col = color, size = cex, alpha = alpha)
  }
}