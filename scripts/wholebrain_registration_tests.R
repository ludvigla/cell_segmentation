library(wholebrain)
library(shiny)
library(ggiraph)
library(ggplot2)
library(raster)
library(EBImage)

filename <- "/Users/ludviglarsson/Michael_Ratz/cell_segmentation/data/test/he_wholebrain_modified.tif"

im <- readImage(filename)
im <- (1 - im) ^2
im <- normalize(im)
#im[im < 0.11] <- 0
#im[im > 0.6] <- 0.6
im <- EBImage::transpose(x = im) #%>% EBImage::rotate(angle = 90)
display(im, method = "raster")
tiff::writeTIFF(im, "/Users/ludviglarsson/Michael_Ratz/cell_segmentation/data/test/he_wholebrain_modified_transposed.tif", bits.per.sample = 16)


MaskImage <- function (
  d,
  image
) {

  # Craete a grid of points to select from
  x <- y <- c()
  for (i in seq(1, ncol(image), 50)) {
    for (j in seq(1, nrow(image), 50)) {
      x <- c(x, i); y <- c(y, j)
    }
  }
  d <- data.frame(x, y, label = 0)
  d$id <- paste0(1:nrow(d))

  # Convert image to raster to use as background
  g <- grid::rasterGrob(
    as.raster(image),
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )
  annotation <- ggplot2::annotation_custom(g, -Inf, Inf, -Inf, Inf)

  # Create UI for shiny app
  ui <- basicPage(
    girafeOutput("Plot1", width = 1000, height = 1000),
    actionButton(inputId = "confirm", label="Confirm selection"),
    actionButton(inputId = "stopApp", label="Quit annotation tool")
  )

  # Create server side
  server <- function(input, output, session) {
    output$Plot1 <- renderGirafe({

      df <- reactiveValues(label = d$label, id = d$id)

      output$Plot1 <- ggiraph::renderGirafe({

        gg <- ggplot(d, aes(x, y, data_id = id)) +
          annotation +
          ggiraph::geom_point_interactive(color = "orange", size = 0.5, alpha = 0.3) +
          theme_void() +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0))
        x <- ggiraph::girafe(ggobj = gg)
        x <- ggiraph::girafe_options(x,
                                     ggiraph::opts_zoom(max = 5),
                                     ggiraph::opts_selection(type = "multiple"))
        x
      })

      observeEvent(input$confirm, {
        ids.selected <- as.numeric(input$Plot1_selected)
        df$label[which(df$id %in% ids.selected)] <- 1
        session$sendCustomMessage(type = 'Plot1_set', message = character(0))
      })

      observe({
        if(input$stopApp > 0){
          print("Stopped")
          d$label <- df$label
          stopApp(returnValue = d)
        }
      })
    })
  }

  m <- runApp(list(ui = ui, server = server))

  return(m)
}

m <- MaskImage(d, im)
m <- subset(m, label == 1)
ps <- chull(m[, 1:2])
p <- c(ps, ps[1])
pl <- spPolygons(apply(m[p, 1:2], 2, as.integer))
r <- raster(vals = 0, ncols = ncol(im), nrows = nrow(im), xmn = 0, xmx = ncol(im), ymn = 0, ymx = nrow(im))
r <- rasterize(pl, r, fun = sum)
mat <- t(as.matrix(r))
im.test <- im
im.test[is.na(mat)] <- 0
display(im.test, method = "raster")

# Now we can export the modified image and rerun the registration
tiff::writeTIFF(im.test, "/Users/ludviglarsson/Michael_Ratz/cell_segmentation/data/test/Neun_clean.tiff", bits.per.sample = 16)

grDevices::quartz()
seg <- segment("/Users/ludviglarsson/Michael_Ratz/cell_segmentation/data/test/he_wholebrain_modified_transposed.tif", filter = NULL, downsample = 0.25)
seg$filter$resize <- 0.08
seg$filter$resize <- 0.16

# registration
regi <- registration(input = "~/Michael_Ratz/cell_segmentation/data/test/he_wholebrain_modified_transposed.tif",
                     coordinate = -1.5, filter = seg$filter, plane = "coronal", right.hemisphere = T)

# add corr.points
regi <- add.corrpoints(registration = regi, n.points = 10)

# rerun registration
regi <- registration(input = "~/Michael_Ratz/cell_segmentation/data/test/Neun_clean.tiff",
                     coordinate = -1.5, correspondance = regi,
                     filter = seg$filter, plane = "coronal", right.hemisphere = T)
