library(EBImage)



## Not run:
data(camphora)
data(cryptomeria)
cryptomeria <- rgb2gray(cryptomeria)
img.c1 <- crop(camphora,200,200)
img.c2 <- crop(cryptomeria,300,300)
par(mfrow=c(2,2))
image(rot90c(img.c1), 
      col=gray(c(0:255)/255), main="Raw.img.c1", useRaster=TRUE, axes=FALSE, asp=1)


image(rot90c(edge.detect(img.c1,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
                         method="Canny")),
      col=gray(c(0:255)/255), main="Canny", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c1,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
                         method="Sobel")),col=gray(c(0:255)/255), main="Sobel", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c2,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
                         method="Canny")),col=gray(c(0:255)/255), main="Canny", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c2,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
                         method="Sobel")),col=gray(c(0:255)/255), main="Sobel", useRaster=TRUE, axes=FALSE, asp=1)
## End(Not run)

edge.detect(x, thresh1=1, thresh2=15, noise="gaussian", noise.s=3, method="Canny")


library("EBImage")

img.sample <- readImage(
  "/Users/chengt/Documents/OCT_Scan/14.tif"
)


img.sample <- readImage(
  "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/10.tif"
)


img.sample <- readImage(
  "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/130.tif"
)
display(img.sample)
img.median <- medianFilter(img.sample, 1)
display(img.median)

fhi <- matrix(1, nrow = 3, ncol = 3)
fhi[2, 2] <- -7.5
display(filter2(x = img.median, filter = fhi))
image(rot90c(edge.detect(
  img.median,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
  method="Sobel")
), col=gray(c(0:255)/255), main="Sobel", useRaster=TRUE, axes=FALSE, asp=1)

# img.median <- img.median[]
attr(img.median, 'bits.per.sample') <- 32
attr(img.median, 'samples.per.pixel') <- 1
img.median <- resize(img.median, w = 1000)
image(edge.detect(img.median))
image(edge.detect(img.median, method = 'Sobel'))
ans <- edge.detect(img.median, method = 'Sobel')

dim(ans) # 664, 156

{
  nmask = thresh(img.median, w=10, h=10, offset=0.01)
  # nmask = opening(nmask, makeBrush(5, shape='disc'))
  # nmask = opening(nmask, makeBrush(5, shape='box'))
  nmask = opening(nmask, makeBrush(5, shape='diamond'))
  nmask = fillHull(nmask)
  nmask = bwlabel(nmask)
  display(nmask, all = T)
}
ctmask = opening(img.median>.7, makeBrush(5, shape='disc'))
cmask = propagate(img.median, seeds=nmask, mask=ctmask)

display(ctmask, all=TRUE)
segmented = paintObjects(cmask, img.median, col='#ff00ff')
segmented = paintObjects(nmask, segmented, col='#ffff00')

display(segmented, all=TRUE)




sum(ans>quantile(ans, .99))







nuc = readImage(system.file("images", "nuclei.tif", package="EBImage"))
display(nuc, method = "raster", all = TRUE)
otsu(nuc)
threshold = otsu(nuc)
threshold
nuc_th <- combine(
  mapply(function(frame, th) frame > th, getFrames(nuc), threshold, SIMPLIFY=FALSE) 
)
display(nuc_th, all=TRUE)



# build a filter from 1750 to 666, take average for linear purpose.
plot(c(ans[,1]))
plot(c(img.sample[,1]))

dim(img.sample) # 1750 , 725
lines(lowess(img.sample[,1], f = .3))

try_local <- img.sample[1:3,]
plot(colMeans(try_local))
lines(lowess(colMeans(try_local), f = .3))


plot(colMeans(nmask[,]))
plot(colMeans(nmask[90:100,]))

plot(nmask)
image(nmask)
image(ctmask)
