se.list <- readRDS("~/10x/R_objects/se.list")

se.bc <- se.list[[4]]

se.bc <- ManualAnnotation(se.bc)

png("~/Desktop/selection.png", height = 1000, width = 2200, res = 150)
FeatureOverlay(se.bc, ncols.samples = 2, features = "labels", sampleids = 1:2, dark.theme = T, type = "raw", pt.alpha = 0.6) %>% print()
dev.off()

se.bc <- SetIdent(se.bc, value = "labels")
se.bc <- RegionNeighbours(se.bc, id = "cancer", keep.within.id = T)

png("~/Desktop/selection_neighbors_keep_within.png", height = 1000, width = 2200, res = 150)
FeatureOverlay(se.bc, ncols.samples = 2, features = "nbs_cancer", sampleids = 1:2, dark.theme = T, type = "raw", pt.alpha = 0.6) %>% print()
dev.off()

png("~/Desktop/IHC.png", height = 3801/2, width = 3671/2, res = 150)
display(EBImage::rotate(im, -90), method = "raster")
dev.off()

cells <- readImage(olig2)
cells <- normalize(cells)
cells_crop = cells[2000:2300, 2000:2300]

png("~/Michael_Ratz/cell_segmentation/images/example_raw_data.png", height = 501, width = 501, res = 150)
display(cells_crop, method = "raster")
dev.off()

png("~/Desktop/example_thr.png", height = 501, width = 501, res = 150)
display(cells_crop_th, method = "raster")
dev.off()

png("~/Michael_Ratz/cell_segmentation/images/example_raw_vs_filtered.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop, method = "raster")
text(x = 10, y = 20, label = "raw", adj = c(0,1), col = "orange", cex = 1.5)
display(cells_crop_filtered, method = "raster")
text(x = 10, y = 20, label = "filtered", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/example_raw_vs_corrected.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop, method = "raster")
text(x = 10, y = 20, label = "raw", adj = c(0,1), col = "orange", cex = 1.5)
display(cells_crop_corrected, method = "raster")
text(x = 10, y = 20, label = "corrected", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()

png("~/Michael_Ratz/cell_segmentation/images/example_raw_vs_thresholded.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop, method = "raster")
text(x = 10, y = 20, label = "raw", adj = c(0,1), col = "orange", cex = 1.5)
display(cells_crop_th, method = "raster")
text(x = 10, y = 20, label = "thresholded", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()

png("~/Michael_Ratz/cell_segmentation/images/example_raw_vs_highlighted.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop, method = "raster")
text(x = 10, y = 20, label = "raw", adj = c(0,1), col = "orange", cex = 1.5)
display(cells_crop_colored, method = "raster")
text(x = 10, y = 20, label = "highlighted cells", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/example_highlighted_vs_clean.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop_colored, method = "raster")
text(x = 10, y = 20, label = "highlighted cells ", adj = c(0,1), col = "orange", cex = 1.5)
display(cells_crop_clean_colored, method = "raster")
text(x = 10, y = 20, label = "highlighted clean cells", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/example_raw_vs_watershed.png", height = 501, width = 1002, res = 150)
par(mfrow=c(1, 2))
display(cells_crop, method = "raster")
text(x = 10, y = 20, label = "raw", adj = c(0,1), col = "orange", cex = 1.5)
display(colorLabels(cells_crop_th_watershed), method = "raster")
text(x = 10, y = 20, label = "labelled cells", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()



png("~/Michael_Ratz/cell_segmentation/images/segmented_cells.png", height = 501, width = 501, res = 150)
display(im, method = "raster")
text(x = 10, y = 20, label = "Oligodendrocytes [red], Egfp [green] \n Neurons [blue]", adj = c(0,1), col = "orange", cex = 1)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/segmented_cells.png", height = 501, width = 501, res = 150)
display(im, method = "raster")
text(x = 10, y = 20, label = "Oligodendrocytes [red], Egfp [green] \n Neurons [blue]", adj = c(0,1), col = "orange", cex = 1)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/overlaps.png", height = 1002, width = 1503, res = 150)
par(mfrow = c(2, 3))
display(im, method = "raster")
text(x = 10, y = 20, label = "Oligodendrocytes [red], Egfp [green] \n Neurons [blue]", adj = c(0,1), col = "orange", cex = 1)
display(oligodendrocyte, method = "raster")
text(x = 10, y = 20, label = "oligodendrocyte + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(neuron, method = "raster")
text(x = 10, y = 20, label = "neurons + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(unknown, method = "raster")
text(x = 10, y = 20, label = "unknown + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
oligneu <- rgbImage(red = oligodendrocyte, green = neuron, blue = unknown)
display(oligneu, method = "raster")
text(x = 10, y = 20, label = "(unknown [blue],  olgiodendrocytes[red] \n and neurons [green]) + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(olig2_neun, method = "raster")
text(x = 10, y = 20, label = "overlapping Neun and Olig2", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()


png("~/Michael_Ratz/cell_segmentation/images/overlaps_full.png", height = 3801/4, width = 3671, res = 150)
par(mfrow = c(1, 4))
display(EBImage::rotate(dilate(oligodendrocyte), angle = -90), method = "raster")
text(x = 10, y = 20, label = "oligodendrocyte + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(EBImage::rotate(dilate(neuron), angle = -90),, method = "raster")
text(x = 10, y = 20, label = "neurons + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(EBImage::rotate(dilate(unknown), angle = -90),, method = "raster")
text(x = 10, y = 20, label = "unknown + EGFP", adj = c(0,1), col = "orange", cex = 1.5)
display(EBImage::rotate(dilate(olig2_neun), angle = -90),, method = "raster")
text(x = 10, y = 20, label = "overlapping Neun and Olig2", adj = c(0,1), col = "orange", cex = 1.5)
dev.off()


p1 <- ggplot() +
  geom_point(data = spots, aes(pixel_x, 2000 - pixel_y), size = 1.2, alpha = 0.6, color = "white") +
  #geom_point(data = nuclei_data, aes(nuclei_x, 2000 - nuclei_y), color = "orange", size = 0.5) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_void() +
  ggtitle("spots") +
  theme(plot.background = element_rect(fill = "black", color = "black"), plot.title = element_text(color = "white"))

p2 <- ggplot() +
  #geom_point(data = spots, aes(pixel_x, 2000 - pixel_y), size = 1.2, alpha = 0.6, color = "white") +
  geom_point(data = nuclei_data, aes(nuclei_x, 2000 - nuclei_y), color = "orange", size = 0.3) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_void() +
  ggtitle("nuclei") +
  theme(plot.background = element_rect(fill = "black", color = "black"), plot.title = element_text(color = "white"))

png("~/Michael_Ratz/cell_segmentation/images/spots_vs_nuclei.png", height = 1000, width = 2000, res = 150)
cowplot::plot_grid(p1, p2) %>% print()
dev.off()

p1 <- ggplot() + geom_point(data = spots, aes(pixel_x, 2000 - pixel_y), size = 1.2, alpha = 0.6, color = "white") +
  geom_point(data = nuclei_data, aes(nuclei_x, 2000 - nuclei_y), color = "orange", size = 0.3) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_void() +
  ggtitle("spots + nuclei") +
  theme(plot.background = element_rect(fill = "black", color = "black"), plot.title = element_text(color = "white"))

png("~/Michael_Ratz/cell_segmentation/images/spots_vs_nuclei_combined.png", height = 1000, width = 1000, res = 150)
p1 %>% print()
dev.off()

p1 <- ggplot() + geom_point(data = spots, aes(pixel_x, 2000 - pixel_y, color = label), size = 1.2, alpha = 0.6) +
  geom_point(data = nuclei_data, aes(nuclei_x, 2000 - nuclei_y), color = "orange", size = 0.3) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_void() +
  scale_color_manual(values = RColorBrewer::brewer.pal(name = "Set1", n = 3)[1:2] %>% rev()) +
  ggtitle("spots + nuclei") +
  theme(plot.background = element_rect(fill = "black", color = "black"),
        plot.title = element_text(color = "white"),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))

png("~/Michael_Ratz/cell_segmentation/images/spots_vs_nuclei_combined_labelled.png", height = 1000, width = 1000, res = 150)
p1 %>% print()
dev.off()

p1 <- ggplot() + geom_point(data = subset(spots, label == "include"), aes(pixel_x, 2000 - pixel_y, color = label), size = 1.2, alpha = 0.6) +
  geom_point(data = nuclei_data, aes(nuclei_x, 2000 - nuclei_y), color = "orange", size = 0.3) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  theme_void() +
  scale_color_manual(values = RColorBrewer::brewer.pal(name = "Set1", n = 3)[1]) +
  ggtitle("spots + nuclei") +
  theme(plot.background = element_rect(fill = "black", color = "black"),
        plot.title = element_text(color = "white"),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))

png("~/Michael_Ratz/cell_segmentation/images/spots_vs_nuclei_combined_hits_only.png", height = 1000, width = 1000, res = 150)
p1 %>% print()
dev.off()
