setwd("~/Michael_Ratz/cell_segmentation/")

segmented <- readRDS("R_objects/segmented")
sizes <- c(0.5, 3, 0.5)

targets <- c("Olig2", "Egfp", "Neun")
fts <- do.call(rbind, lapply(seq_along(segmented), function(i) {
  seg <- segmented[[i]]
  df <- computeFeatures.moment(seg)
  data.frame(x = df[, "m.cx"], y = df[, "m.cy"], target = targets[i], z = i, size = sizes[i])
}))

library(plotly)

plot_ly(fts,
        x = ~max(fts$x) - x, y = ~y, z = ~z, color = ~target, size = ~size, sizes = c(0.5, 3), colors = c("red", "green", "blue"),
        marker = list(showscale = FALSE,
                      opacity = 1)) %>%
  add_markers() %>%
  layout(paper_bgcolor = 'rgb(0, 0, 0)',
         scene = list(zaxis = list(title = '', range = c(-2, 3 + 2), showticks = FALSE, showticklabels = FALSE),
                      xaxis = list(title = '', showticks = FALSE, showticklabels = FALSE),
                      yaxis = list(title = '', showticks = FALSE, showticklabels = FALSE)))


plot_ly(fts,
        x = ~max(fts$x) - x, y = ~y, z = 1, color = ~target, size = ~size, sizes = c(0.5, 3), colors = c("red", "green", "blue"),
        marker = list(showscale = FALSE,
                      opacity = 1)) %>%
  add_markers() %>%
  layout(paper_bgcolor = 'rgb(0, 0, 0)',
         scene = list(zaxis = list(title = '', range = c(-2, 3 + 2), showticks = FALSE, showticklabels = FALSE),
                      xaxis = list(title = '', showticks = FALSE, showticklabels = FALSE),
                      yaxis = list(title = '', showticks = FALSE, showticklabels = FALSE)))


# Or just to show a dummy example of multiple sections
fts_subset <- subset(fts, z == 1)
fts_large <- do.call(rbind, lapply(1:6, function(i) {
  d <- fts_subset
  d$z <- i
  return(d)
}))

plot_ly(fts_large,
        x = ~max(fts$x) - x, y = ~y, z = ~z, color = ~target, colors = "red",
        marker = list(showscale = FALSE,
                      size = 0.5,
                      opacity = 1)) %>%
  add_markers() %>%
  layout(paper_bgcolor = 'rgb(0, 0, 0)',
         scene = list(zaxis = list(title = '', range = c(0, 6), showticks = FALSE, showticklabels = FALSE),
                      xaxis = list(title = '', showticks = FALSE, showticklabels = FALSE),
                      yaxis = list(title = '', showticks = FALSE, showticklabels = FALSE)))
