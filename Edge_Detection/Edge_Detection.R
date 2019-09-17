library(EBImage)
library(plyr)
library(dplyr)
library(doParallel)
# detach('package:EBImage', unload = T)

# can read them directly; must be in the local disk
function() {
    
}

list_all <- EBImage::readImage('~/Documents/OCT_Scan/Flowcell/Resized_Raw/Resized_26_.tif')

img.sample <- list_all[,,1] %>% resize(w = 1000, h = 400)

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

img.median <- medianFilter(x = img.sample, size = 2)
plot(img.median)

Get_BioFilmBoundary <- function(
    image = NULL,
    fname = '/Users/chengt/Documents/OCT_Scan/Edge_Detection/Day_32_04.04_Edited_0333.tif',
    size_median_filter = 2,
    thresh_offset = .03,
    size_movingaveragemask_filter = 50,
    frac_lowess_top_min = 1/10,
    frac_lowess_bottom_max = 1/6,
    flag_plot = F
) {
    require(EBImage)
    if (is.null(image)) {
        img_raw <- readImage(fname)
    } else {img_raw <- image}
    img_temp <- img_raw
    img_temp <- waveslim::denoise.modwt.2d(img_temp)
    img_temp <- medianFilter(x = img_temp, size = size_median_filter)
    img_temp <- waveslim::denoise.modwt.2d(img_temp)
    img_temp <- medianFilter(x = img_temp, size = size_median_filter)
    img_temp <- waveslim::denoise.modwt.2d(img_temp)
    img_temp <- medianFilter(x = img_temp, size = size_median_filter)
    # img_temp <- medianFilter(x = img_temp, size = 10)
    img_temp <- as.Image(img_temp)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    {
        nmask <- thresh(img_temp, w=1, h=19, offset = thresh_offset)
        nmask <- opening(nmask, makeBrush(7, shape='disc'))
        nmask <- opening(nmask, makeBrush(7, shape='disc'))
        nmask <- opening(nmask, makeBrush(3, shape='line', angle = 0))
        nmask <- closing(nmask, makeBrush(5, shape='line', angle = 0))
        nmask <- closing(nmask, makeBrush(3, shape='diamond'))
        nmask <- fillHull(nmask)
    }
    # Regional Edge Detection for the Biofilm, min as top & max as bottom
    index_last_initial <- NROW(nmask) - size_movingaveragemask_filter
    seq_top_min <- seq_bottom_max <- c()
    for (iter_row in 1:index_last_initial) {
        iter_df <- nmask[iter_row:(iter_row + size_movingaveragemask_filter), ]
        iter_df <- colMeans(iter_df)
        iter_peaks <- quantmod::findPeaks(iter_df)
        iter_length <- max(length(iter_peaks)*5/9, 1)
        seq_top_min <- c(seq_top_min, 
                         median(tail(head(iter_peaks, iter_length), 2), 
                                iter_peaks, na.rm = T)
        )
        # 
        iter_length <- max(length(iter_peaks)*1/11, 1)
        seq_bottom_max <- c(seq_bottom_max, 
                            median(head(tail(iter_peaks, iter_length), 2), 
                                   iter_peaks, na.rm = T)
        )
    }
    # 
    seq_top_min <- imputeTS::na_interpolation(seq_top_min, option = 'stine')
    seq_bottom_max <- imputeTS::na_interpolation(seq_bottom_max, option = 'stine')
    seq_top_min <- lowess(seq_top_min, f = frac_lowess_top_min)
    seq_bottom_max <- lowess(seq_bottom_max, f = frac_lowess_bottom_max)
    seq_bottom_max$x <- seq_top_min$x <- 1:NROW(img_raw)
    seq_bottom_max$y <- c(seq_bottom_max$y, rep(NA, size_movingaveragemask_filter))
    seq_top_min$y <- c(seq_top_min$y, rep(NA, size_movingaveragemask_filter))
    seq_top_min <- lowess(seq_top_min, f = 1/10)
    seq_bottom_max <- lowess(seq_bottom_max, f = 1/6)
    seq_top_min$y <- imputeTS::na_interpolation(seq_top_min$y, option = 'stine')
    seq_bottom_max$y <- imputeTS::na_interpolation(seq_bottom_max$y, option = 'stine')
    # 
    if (sum(seq_bottom_max$y > seq_top_min$y)) {
        index_temp <- seq_bottom_max$y > seq_top_min$y
        temp <- seq_bottom_max$y
        seq_bottom_max$y[index_temp] <- seq_top_min$y[index_temp]
        seq_top_min$y[index_temp] <- temp[index_temp]
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    if (flag_plot) {
        par(mfrow = c(2,1))
        # 'TOP' is actually 'BOTTOM' in the image, and vice versa
        # 
        plot(nmask)
        lines(seq_top_min, col = "cyan", lwd = 1.5)
        lines(seq_bottom_max, col = 'deeppink', lwd = 1.5)
        legend('topright', legend = c('top', 'bottom'), col = c('deeppink', 'cyan'), 
               text.col = 'black', lwd = 2, box.lwd = 0)
        # 
        plot(img_temp)
        lines(seq_top_min, col = "cyan", lwd = 1.5)
        lines(seq_bottom_max, col = 'deeppink', lwd = 1.5)
        legend('topright', legend = c('top', 'bottom'), col = c('deeppink', 'cyan'), 
               text.col = 'black', lwd = 2, box.lwd = 0)
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    return(list(
        seq_top_min = seq_top_min,
        seq_bottom_max = seq_bottom_max
    ))
}

Get_BioFilmBoundary(
    fname = "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/10.tif",
    flag_plot = T, thresh_offset = .04, size_movingaveragemask_filter = 100
)
Get_BioFilmBoundary(
    fname = "/Users/chengt/Documents/OCT_Scan/Edge_Detection/Day_32_04.04_Edited_0333.tif",
    flag_plot = T, thresh_offset = .03, size_movingaveragemask_filter = 100
)
Get_BioFilmBoundary(
    fname = "/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/MBR_1_2D/2D for terry/130.tif",
    flag_plot = T, thresh_offset = .04, size_movingaveragemask_filter = 100
)



result <- Get_BioFilmBoundary(
    fname = '~/Documents/OCT_Scan/Flowcell/seq_tiff_/Resized_1_0000.tif', 
    flag_plot = T, thresh_offset = .04, size_movingaveragemask_filter = 60
)




quantmod::findPeaks(c())
min (peaks)

plot(img_temp[1:100,] %>% colMeans)
plot(nmask[1:100,] %>% colMeans)
plot(nmask)
temp <- colMeans(nmask[1:10,])
plot(temp)

valleys <- quantmod::findValleys(temp)
peaks <- quantmod::findPeaks(temp)


points(x = peaks, y = rep(1, length(peaks)), col = 'red')
points(x = valleys, y = rep(1, length(valleys)), col = 'blue')
object.size(nmask)

EBImage::combine(nmask, flip(nmask), flop(nmask)) %>% object.size()
abind(nmask, flip(nmask), flop(nmask), along = 2) %>% plot


# Write again
{}
# Best ever seen



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

# can read them directly; must be in the local disk
list_all <- EBImage::readImage('~/Documents/OCT_Scan/Flowcell/Resized_Raw/Resized_26_.tif')