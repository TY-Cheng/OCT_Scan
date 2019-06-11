rm(list = ls())
setwd("/Users/chengt/Documents/OCT_Scan")
load("Thickness_3D_Raw_.RData")
getwd()
doParallel::registerDoParallel(8)

library(scales)
library(gstat)
library(sp)
# library(gapfill)
library(ggplot2)
library(colorRamps)

# Fiji_Macro_ -------------------------------------------------------------
# 
# setwd('/Volumes/Seagate_Backup/OCT_Scan/OCT_2D_for_Terry_/2D for terry/')
# # Check only tif files left
# fname <- list.files()
# 
# for (iter_i in fname) {
#     cat(
#         paste0(
#             'open("/Volumes/Seagate_Backup/OCT_Scan/OCT_2D_for_Terry_/2D for terry/',
#             iter_i,
#             '");'
#         ), file = 'Fiji_Macro_.txt', sep = '\n', append = T)
#     cat(
#         paste0(
#             'selectWindow("',
#             iter_i,
#             '");'
#         ), file = 'Fiji_Macro_.txt', sep = '\n', append = T)
#     cat(
#         'run("Size...", "width=1970 height=414 depth=1 average interpolation=Bilinear");
#         close();'
#         , file = 'Fiji_Macro_.txt', sep = '\n', append = T)
# }


# Copy&Paste --------------------------------------------------------------
if (0) {
    # assign & read by name in string
    # assign('x', 1)
    # eval(as.symbol('x')) 
    # copy&paste
    table <- clipr::read_clip_tbl(sep = ',')
    table
    clipr::write_clip(table)
}

# Read & Impute -----------------------------------------------------------
if(0) {    
    library(Amelia)
    # 
    df <- clipr::read_clip_tbl()
    rownames(df) <- (1:nrow(df)) + 1
    df$index <- seq_len(nrow(df_p))
    DF <- df_p <- df[,c('index', 'mean', 'sd', 
                        'roughness_Ra', 'roughness_Ra1')]
    # 
    my_cols <- 1:length(df_p)
    pairs(x = df_p,
          upper.panel = panel.cor,
          diag.panel = panel.hist,
          lower.panel = panel.smooth,
          cex.labels = .8)
    if(1){
        DF[df$image=='109.tif',] <- NA
        DF[df$image=='126.tif',] <- NA
        DF[df$image=='127.tif',] <- NA
        DF$index <- seq_len(nrow(DF))
        # 
        set.seed(2019)
        niter <- 100
        a.out <- amelia(DF, ts = 'index',
                        m = niter, p2s = 2, splinetime = 3,
                        parallel = 'multicore'
        )
        missmap(a.out)
        colSums(a.out$missMatrix)
        # combine results
        a <- matrix(0, nrow = nrow(DF), ncol = ncol(DF))
        for(i in 1:niter){
            a <- a + 
                a.out$imputations[[i]][,
                                       !names(DF)%in%c("Date","Day")]
        }
        a <- a/niter
        clipr::write_clip(a)
        rm(niter, i,DF,a,a.out)
        rm(df, df_p, my_cols)
    }
}


# Import Matrix From MATfiles --------------------------------------------------------
Import_Matrix_From_MATfiles <- function(Working_Directory = getwd()){
    # 
    print(Working_Directory)
    setwd(Working_Directory)
    # Initialize the 3D Thickness matrix
    Img_3D <- matrix(data = NA, nrow = 666, ncol = 1970)
    # 
    seq_MATfname <- list.files()
    for (iter_MATfname in seq_MATfname) {
        # Load the mat file, sequentially
        Img_2D <- R.matlab::readMat(iter_MATfname)
        Img_2D <- Img_2D$imageLayer[[3]]
        # numerically, the result should be 'isos' - 'ilm'
        if (Img_2D[,1,1]$name =='isos') {
            isos_y <- as.numeric(Img_2D[,1,1]$pathY)
            isos_x <- as.numeric(Img_2D[,1,1]$pathX)
        }
        if (Img_2D[,1,2]$name =='isos') {
            isos_y <- as.numeric(Img_2D[,1,2]$pathY)
            isos_x <- as.numeric(Img_2D[,1,2]$pathX)
        }
        if (Img_2D[,1,1]$name =='ilm') {
            ilm_y <- as.numeric(Img_2D[,1,1]$pathY)
            ilm_x <- as.numeric(Img_2D[,1,1]$pathX)
        }
        if (Img_2D[,1,2]$name =='ilm') {
            ilm_y <- as.numeric(Img_2D[,1,2]$pathY)
            ilm_x <- as.numeric(Img_2D[,1,2]$pathX)
        }
        # if isos & ilm length wont match, raise alert
        dnt_match <- sum(unique(isos_y)-unique(ilm_y))
        if (dnt_match) {
            cat('\n', iter_MATfname)
            cat('\nisos & ilm match?\t', !dnt_match)
        }
        # RowNumber: take the last number
        row_num <- stringr::str_extract_all(string = iter_MATfname, 
                                            pattern = '\\d+')[[1]]
        row_num <- as.numeric(tail(row_num, n = 1)) +1
        # ColumnNumber & extract the secant (thickness of the row)
        thickness_row <- c()
        col_range <- unique(isos_y)
        for (col_num in col_range) {
            # print(col_num)
            if (col_num/max(col_range)<=.5) {
                thickness_row[col_num] <- 
                    max(isos_x[isos_y==col_num]) - max(ilm_x[ilm_y==col_num])
            }else{
                thickness_row[col_num] <- 
                    min(isos_x[isos_y==col_num]) - min(ilm_x[ilm_y==col_num])
            }
        }
        Img_3D[row_num, col_range] <- thickness_row
    }
    return(Img_3D)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Import MAT file ---------------------------------------------------------
if(1){
    Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan/MATfiles_/'
    setwd(Rootfolder)
    # list.files()
    seq_MATfolder <- list.files()
    # 
    for (iter_MATfolder in seq_MATfolder) {
        print(iter_MATfolder)
        setwd(paste0(Rootfolder, iter_MATfolder))
        # print(getwd())
        assign(x = paste0(iter_MATfolder, 'Img_3D_'), 
               value = Import_Matrix_From_MATfiles(
                   Working_Directory = paste0(Rootfolder, iter_MATfolder)
               )
        )
        cat(iter_MATfolder, '\t\tfinished.\n\n')
    }
    # setwd('/Volumes/Seagate_Backup/OCT_Scan/')
    # save.image("/Volumes/Seagate_Backup/OCT_Scan/Thickness_3D_Raw_.RData")
}


# Crop_Denoise_Image_from_matrix ------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Crop_Denoise_Image_from_matrix <- function(Img_3D, 
                                           flag_col_number_check = 1,
                                           flag_ex_extreme_value = 1,
                                           flag_plot = 1,
                                           flag_save_plot = 1,
                                           denoise_type = 1,
                                           multiplier_pixel2micron = 2.1,
                                           quantiles = c(.01, .99),
                                           scale_range = 60,
                                           save_folder = ''){
    # call by var_name (string)
    # deparse(substitute(Img_3D))
    cat('\n', Img_3D, '\n', 'Denoising...\n')
    Thickness <- eval(as.symbol(Img_3D))
    # From pixels to microns
    Thickness <- Thickness * multiplier_pixel2micron
    # raise alarm if col_number is not '987' for all rows!
    if(flag_col_number_check){
        num_notna <- apply(Thickness, MARGIN = 1, function(x)sum(!is.na(x)))
        dnt_match <- sum((num_notna == 0 )|(num_notna == 987)) != 666
        if (dnt_match) {
            cat('\n\n', Img_3D, '\n col_number is not 987! \n')
            print(unique(num_notna))
        }
        rm(num_notna)
    }
    # print(dim(Thickness))
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Cropping, Resizing, Denoising
    # impute first, resize second
    if (denoise_type==1) {
        Thickness <- na.omit(Thickness[,1:987])
        Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
        # 
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        rm(Thickness_max)
    }
    if (denoise_type==2) {
        Thickness <- na.omit(Thickness[,1:987])
        threshold_row <- lowess(rowMeans(Thickness, na.rm = T), f = 1/2)$y - 
            rowMeans(Thickness, na.rm = T) / 7
        Thickness[rowMeans(Thickness, na.rm = T) < threshold_row, ] <- NA
        Thickness <- na.omit(Thickness)
        # 
        threshold_row <- lowess(rowMeans(Thickness, na.rm = T), f = 1/3)$y - 
            sd(rowMeans(Thickness, na.rm = T)) * 2
        Thickness[rowMeans(Thickness, na.rm = T) < threshold_row, ] <- NA
        Thickness <- na.omit(Thickness)
        # 
        Thickness <- rbind(
            EBImage::resize(
                Thickness[1:(NROW(Thickness)/3), ], 
                w = NROW(Thickness)*2/5, h = 666),
            EBImage::resize(
                Thickness[(NROW(Thickness)/3+1):(NROW(Thickness)*2/3-1), ], 
                w = 666 - NROW(Thickness)/2, h = 666),
            EBImage::resize(
                Thickness[(NROW(Thickness)*2/3):NROW(Thickness), ], 
                w = NROW(Thickness)*2/5, h = 666)
        )
        Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
        #
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        rm(Thickness_max)
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if (denoise_type==3) {
        set.seed(2019)
        Thickness <- Thickness[,1:987]
        threshold_row <- lowess(rowMeans(Thickness, na.rm = T), f = 1/2)$y - 
            rowMeans(Thickness, na.rm = T) / 7
        Thickness[rowMeans(Thickness, na.rm = T) < threshold_row, ] <- NA
        Thickness <- na.omit(Thickness)
        threshold_row <- lowess(rowMeans(Thickness, na.rm = T), f = 1/3)$y - 
            sd(rowMeans(Thickness, na.rm = T)) * 2
        Thickness[rowMeans(Thickness, na.rm = T) < threshold_row, ] <- NA
        # 
        if (sum(is.na(Thickness))>0) {
            Img_3D_grid <- expand.grid(
                X = seq_len(NCOL(Thickness)),
                Y = seq_len(NROW(Thickness))
            )
            Img_3D_grid$Z <- c(t(Thickness)) # Transposed!
            coordinates(Img_3D_grid) <- ~X+Y
            Img_3D_grid_NAfree <- Img_3D_grid[!is.na(Img_3D_grid$Z), ]
            # Variogram
            print(Sys.time())
            tictoc::tic('Iso_Variogram Calculation...')
            Img_variogram <- variogram(
                Z ~ 1, 
                Img_3D_grid_NAfree[sample(NROW(Img_3D_grid_NAfree), 16180),],
                cutoff = 300,
                width = 10,
                verbose = TRUE
            )
            tictoc::toc()
            # Kriging
            tictoc::tic('Kriging...')
            Img_variogram_fit <- fit.variogram(
                Img_variogram, 
                model = vgm(model = 'Ste')
            )
            # print(Img_variogram_fit)
            Img_3D_grid_kriged <- krige(
                formula = Z ~ 1, 
                locations = Img_3D_grid_NAfree[sample(1314), ],
                newdata = Img_3D_grid[is.na(Img_3D_grid$Z), ],
                model = Img_variogram_fit
            )
            print(summary(Img_3D_grid_kriged))
            # 
            Img_3D_grid_kriged <- data.frame(
                coordinates(Img_3D_grid_kriged), 
                Z = Img_3D_grid_kriged$var1.pred
            )
            Thickness_kriged <- t(Thickness) # transposed at first
            apply(
                X = Img_3D_grid_kriged, 
                MARGIN = 1, 
                FUN = function(vec){
                    Thickness_kriged[vec[1], vec[2]] <<- vec[3]
                    return(NULL)
                }
            )
            Thickness <- Thickness_kriged
            #
            tictoc::toc()
            print(Sys.time())
            rm(Img_3D_grid, Img_3D_grid_NAfree, Img_3D_grid_kriged, 
               Thickness_kriged)
        }
        # 
        Thickness <- rbind(
            EBImage::resize(
                Thickness[1:(NROW(Thickness)/3), ], 
                w = NROW(Thickness)*2/5, h = 666),
            EBImage::resize(
                Thickness[(NROW(Thickness)/3+1):(NROW(Thickness)*2/3-1), ], 
                w = 666 - NROW(Thickness)/2, h = 666),
            EBImage::resize(
                Thickness[(NROW(Thickness)*2/3):NROW(Thickness), ], 
                w = NROW(Thickness)*2/5, h = 666)
        )
        # 
        Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        rm(Thickness_max)
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 
    if (flag_ex_extreme_value) {
        # To exclude extreme values
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = quantiles[1]), 
            0)
        maximum <- quantile(Thickness, probs = quantiles[2])
        Thickness_wo_Extreme[Thickness_wo_Extreme<minimum] <- minimum
        Thickness_wo_Extreme[Thickness_wo_Extreme>maximum] <- maximum
        Thickness <- Thickness_wo_Extreme
        rm(Thickness_wo_Extreme, maximum, minimum)
    }
    #
    # Expand the Grid to plot
    Img_3D_grid <- expand.grid(
        X = seq_len(NCOL(Thickness)),
        Y = seq_len(NROW(Thickness))
    )
    Img_3D_grid$Z <- c(t(Thickness))
    # 
    # Levelplot with ggplot2
    p <- ggplot(Img_3D_grid, aes(x = X, y = Y, z = Z)) +
        geom_raster(aes(fill = Z)) +
        coord_fixed() +
        labs(title = paste0(Img_3D, 'Denoised_'),
             x = '', y = '',
             fill = expression(paste('Thickness (',mu,'m)'))
        ) +
        theme_bw()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # scale_range
    # = NULL:       default limits
    # = one number: fixed scale width
    # = two number: fixed scale range
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if (is.null(scale_range)) {
        p <- p + 
            scale_fill_gradientn(
                colours = matlab.like(67)
            )
    }else if (length(scale_range) == 2) {
        p <- p +
            scale_fill_gradientn(
                limits = c(scale_range[1], scale_range[2]),
                colours = matlab.like(67)
            )
    }else if (length(scale_range) == 1) {
        p <- p +
            scale_fill_gradientn(
                limits = c(
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0),
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0) + 
                        scale_range
                ),
                colours = matlab.like(67)
            )
    }
    if (flag_plot) {print(p)}
    # 
    if (flag_save_plot) {
        cat('Saving to', paste0(save_folder, Img_3D, 'Denoised_.png'), '...\n')
        ggsave(filename = paste0(save_folder, Img_3D, 'Denoised_.png'), 
               plot = p, width = 9)
    }
    # 
    return(Thickness)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Cropping, Denoising, Imaging -----------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Manual
    rm(list = ls())
    setwd("/Users/chengt/Documents/OCT_Scan")
    load("Thickness_3D_Raw_.RData")
    df_Img_3D <- data.frame(seq_Img_3D = c("Day_01_04.03_Resized_Img_3D_", "Day_06_09.03_Resized_Img_3D_", "Day_11_14.03_Resized_Img_3D_", "Day_16_19.03_Resized_Img_3D_", "Day_21_24.03_Resized_Img_3D_", "Day_23_26.03_Resized_Img_3D_", "Day_24_27.03_Resized_Img_3D_", "Day_25_28.03_Resized_Img_3D_", "Day_26_29.03_Resized_Img_3D_", "Day_27_30.03_Resized_Img_3D_", "Day_28_31.03_Resized_Img_3D_", "Day_29_01.04_Resized_Img_3D_", "Day_30_02.04_Resized_Img_3D_", "Day_32_04.04_Resized_Img_3D_", "Day_35_07.04_Resized_Img_3D_", "Day_37_09.04_Resized_Img_3D_", "Day_40_12.04_Resized_Img_3D_", "Day_42_14.04_Resized_Img_3D_", "Day_44_16.04_Resized_Img_3D_", "Day_46_18.04_Resized_Img_3D_", "Day_49_21.04_Resized_Img_3D_", "Day_52_24.04_Resized_Img_3D_", "Day_53_25.04_Resized_Img_3D_", "Day_56_28.04_Resized_Img_3D_", "Day_60_02.05_1_Resized_Img_3D_", "Day_60_02.05_10_Resized_Img_3D_", "Day_60_02.05_11_Resized_Img_3D_", "Day_60_02.05_12_Resized_Img_3D_", "Day_60_02.05_13_Resized_Img_3D_", "Day_60_02.05_14_Resized_Img_3D_", "Day_60_02.05_15_Resized_Img_3D_", "Day_60_02.05_2_Resized_Img_3D_", "Day_60_02.05_3_Resized_Img_3D_", "Day_60_02.05_4_Resized_Img_3D_", "Day_60_02.05_5_Resized_Img_3D_", "Day_60_02.05_6_Resized_Img_3D_", "Day_60_02.05_7_Resized_Img_3D_", "Day_60_02.05_8_Resized_Img_3D_", "Day_60_02.05_9_Resized_Img_3D_"), 
                            seq_denoise_type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                            # seq_denoise_type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                            stringsAsFactors = F)
    # 
    # df_Img_3D <- df_Img_3D[df_Img_3D$seq_denoise_type==3,]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Load function here
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if(1){
        Rootfolder <- '/Users/chengt/Documents/OCT_Scan/Img/'
        # 1 Default
        scale_range <- NULL
        for (iter_Img_3D in 1:NROW(df_Img_3D)) {
            Img_3D <- df_Img_3D[iter_Img_3D, 'seq_Img_3D']
            denoise_type <- df_Img_3D[iter_Img_3D, 'seq_denoise_type']
            print(Img_3D)
            assign(
                paste0(Img_3D,'Denoised_'),
                Crop_Denoise_Image_from_matrix(
                    Img_3D = Img_3D, 
                    flag_plot = F,
                    denoise_type = denoise_type, 
                    scale_range = scale_range, 
                    save_folder = paste0(Rootfolder, 'Default/'))
            )
        }
        # 2 Fixed Scale Width
        scale_range <- 60
        for (iter_Img_3D in 1:NROW(df_Img_3D)) {
            Img_3D <- df_Img_3D[iter_Img_3D, 'seq_Img_3D']
            denoise_type <- df_Img_3D[iter_Img_3D, 'seq_denoise_type']
            print(Img_3D)
            assign(
                paste0(Img_3D,'Denoised_'),
                Crop_Denoise_Image_from_matrix(
                    Img_3D = Img_3D, 
                    flag_plot = F,
                    denoise_type = denoise_type, 
                    scale_range = scale_range, 
                    save_folder = paste0(Rootfolder, 'Fixed_Scale_Width/'))
            )
        }
        # 3 Fixed Scale Range
        scale_range <- c(0,151)
        for (iter_Img_3D in 1:NROW(df_Img_3D)) {
            Img_3D <- df_Img_3D[iter_Img_3D, 'seq_Img_3D']
            denoise_type <- df_Img_3D[iter_Img_3D, 'seq_denoise_type']
            print(Img_3D)
            assign(
                paste0(Img_3D,'Denoised_'),
                Crop_Denoise_Image_from_matrix(
                    Img_3D = Img_3D, 
                    flag_plot = F,
                    denoise_type = denoise_type, 
                    scale_range = scale_range, 
                    save_folder = paste0(Rootfolder, 'Fixed_Scale_Range/'))
            )
        }
    }
}


# Test for colors
if(0){
    # 
    ggplot(faithfuld, aes(waiting, eruptions, z = density)) + 
        geom_tile(aes(fill = density)) +
        geom_contour(bins = 10) + # or binwidth
        scale_fill_gradient(low = 'deepskyblue',
                            #'cyan'
                            high = 'firebrick1'
                            #'deeppink'
        ) + 
        # coord_fixed() +
        theme_bw()
}
