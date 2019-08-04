
library("EBImage")


img.sample <- readImage(
  "/Users/chengt/Documents/OCT_Scan/Edge_Detection/Day_32_04.04_Edited_0333.tif"
)
img.sample <- readImage(
  "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/10.tif"
)
img.sample <- readImage(
  "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/130.tif"
)
plot(img.sample)
img.median <- medianFilter(img.sample, 1)
plot(img.median)

# Write again
{}
# Best ever seen
{
  nmask = thresh(img.median, w=1, h=19, offset=0.01)
  nmask = opening(nmask, makeBrush(7, shape='disc'))
  nmask = opening(nmask, makeBrush(3, shape='line', angle = 0))
  
  nmask = closing(nmask, makeBrush(5, shape='line', angle = 0))
  nmask = closing(nmask, makeBrush(3, shape='diamond'))
  nmask = fillHull(nmask)
  plot(nmask, all = T)
}



{
  # img.median <- medianFilter(img.sample, 1)
  # img.median <- EBImage::resize(img.median, w = 1000)
  nmask = thresh(img.median, w=10, h=10, offset=0.01)
  # nmask = opening(nmask, makeBrush(5, shape='disc'))
  # nmask = opening(nmask, makeBrush(5, shape='box'))
  nmask = opening(nmask, makeBrush(5, shape='diamond'))
  nmask = fillHull(nmask)
  nmask = bwlabel(nmask)
  display(nmask, all = T)
  
  ctmask = opening(img.median, makeBrush(5, shape='disc'))
  cmask = propagate(img.median, seeds=nmask, mask=ctmask)
  display(ctmask, all=TRUE)
  display(cmask)
}

plot(colMeans(img.median))
plot(colMeans(nmask))
dim(nmask)
dim(img.median)
plot(colMeans(img.median[1:100,]))

