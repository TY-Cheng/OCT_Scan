library(doParallel)
library(plyr)

library(scales)
library(gstat)
library(sp)
# library(gapfill)
library(ggplot2)
library(colorRamps)


# Get_Thickness_2D_From_MAT -----------------------------------------------
# Import secant from one single mat file
Get_Thickness_2D_From_MAT <- function(MATfname) {
    # load the mat file
    Img_2D <- R.matlab::readMat(MATfname)
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
        cat('\n', MATfname)
        cat('\nisos & ilm match?\t', !dnt_match)
    }
    # RowNumber: take the last number
    row_num <- stringr::str_extract_all(string = MATfname, 
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
    # assign(x = paste0('row_', row_num), value = thickness_row)
    # return(eval(as.symbol(paste0('row_', row_num))))
    return(list(row_num = row_num, thickness_row = thickness_row))
}


# Get_Thickness_3D_From_MAT -----------------------------------------------
# from the folder, import all secant MAT files, 
# splice 2D thickness_row into one 3D matrix 'Img_3D'
Get_Thickness_3D_From_MAT <- function(Working_Directory, cl) {
    tictoc::tic(Working_Directory)
    old_wd <- getwd()
    print(Working_Directory)
    setwd(Working_Directory)
    # get all mat file names inside the given location
    # numbered well from 0 to 665, no worry
    seq_MATfname <- sort(list.files())
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    Img_2D_list <- parLapplyLB(
        cl = para_socket_cl,
        X = seq_MATfname, 
        fun = Get_Thickness_2D_From_MAT
    )
    stopCluster(para_socket_cl)
    Img_3D <- matrix(NA, nrow = 666, ncol = 1970)
    invisible(lapply(
        X = Img_2D_list,
        FUN = function(rownum_thicknessrow){
            Img_3D[rownum_thicknessrow$row_num,
                   1:length(rownum_thicknessrow$thickness_row)] <<- 
                rownum_thicknessrow$thickness_row}
    ))
    Img_3D <- list(Img_3D[, !colSums(is.na(Img_3D))==666])
    names(Img_3D) <- Working_Directory
    setwd(old_wd)
    tictoc::toc()
    return(Img_3D)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Working_Directory <- list.files()[30]

# Import MAT file ---------------------------------------------------------
if(1){
    Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan/MBR_1_3D/MAT_files_/'
    setwd(Rootfolder)
    seq_MATfolder <- sort(list.files())
    Img_list_MBR1 <- parLapply(seq_MATfolder, Get_Thickness_3D_From_MAT)
    names(Img_list_MBR1) <- seq_MATfolder
    # 
    Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan/MBR_2_3D/MAT_files_/'
    setwd(Rootfolder)
    seq_MATfolder <- sort(list.files())
    Img_list_MBR2 <- lapply(X = seq_MATfolder, FUN = Get_Thickness_3D_From_MAT)
    names(Img_list_MBR2) <- seq_MATfolder
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
                                           quantiles = c(.027, .97),
                                           scale_range = 60,
                                           Fig_Title,
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
        Thickness <- rbind(
            EBImage::resize(
                na.omit(Thickness[
                    1:(NROW(Thickness)*1/3), 
                    ]), 
                w = 300, h = 666),
            EBImage::resize(
                na.omit(Thickness[
                    (NROW(Thickness)*1/3):(NROW(Thickness)*2/3), 
                    ]), 
                w = 300, h = 666),
            EBImage::resize(
                na.omit(Thickness[
                    (NROW(Thickness)*2/3):(NROW(Thickness)), 
                    ]), 
                w = 300, h = 666)
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
        # 
        # Thickness <- na.omit(Thickness)
        # Thickness <- rbind(
        #     EBImage::resize(
        #         Thickness[1:(NROW(Thickness)/3), ], 
        #         w = NROW(Thickness)*2/5, h = 666),
        #     EBImage::resize(
        #         Thickness[(NROW(Thickness)/3+1):(NROW(Thickness)*2/3-1), ], 
        #         w = 666 - NROW(Thickness)/2, h = 666),
        #     EBImage::resize(
        #         Thickness[(NROW(Thickness)*2/3):NROW(Thickness), ], 
        #         w = NROW(Thickness)*2/5, h = 666)
        # )
        # 
        Thickness <- rbind(
            EBImage::resize(
                na.omit(Thickness[
                    1:(NROW(Thickness)*1/3), 
                    ]), 
                w = 300, h = 666),
            EBImage::resize(
                na.omit(Thickness[
                    (NROW(Thickness)*1/3):(NROW(Thickness)*2/3), 
                    ]), 
                w = 300, h = 666),
            EBImage::resize(
                na.omit(Thickness[
                    (NROW(Thickness)*2/3):(NROW(Thickness)), 
                    ]), 
                w = 300, h = 666)
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
    Img_3D_grid$X <- Img_3D_grid$X / 666 * 4
    Img_3D_grid$Y <- Img_3D_grid$Y / 666 * 4
    Img_3D_grid$Z <- c(t(Thickness))
    # 
    # Levelplot with ggplot2
    p <- ggplot(Img_3D_grid, aes(x = X, y = Y, z = Z)) +
        geom_raster(aes(fill = Z)) +
        coord_fixed() +
        labs(
            # title = Fig_Title,
            x = expression(paste('X (mm)')),
            y = expression(paste('Y (mm)')),
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
                colours = matlab.like(300)
            )
    }else if (length(scale_range) == 2) {
        p <- p +
            scale_fill_gradientn(
                limits = c(scale_range[1], scale_range[2]),
                colours = matlab.like(300)
            )
    }else if (length(scale_range) == 1) {
        p <- p +
            scale_fill_gradientn(
                limits = c(
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0),
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0) + 
                        scale_range
                ),
                colours = matlab.like(300)
            )
    }
    if (flag_plot) {print(p)}
    # 
    if (flag_save_plot) {
        cat('Saving to', paste0(save_folder, Fig_Title, '.png'), '...\n')
        ggsave(filename = paste0(save_folder, Fig_Title, '.png'), 
               plot = p, width = 7)
    }
    # 
    return(Thickness)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# Imputation --------------------------------------------------------------
# Impute based on realized features, 
# for further prediction/regression/classification
if (0) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    set.seed(2019)
    tictoc::tic('Imputation...')
    ts_imputed <- data.frame(Date = 1:NROW(ts_raw), ts_raw)
    ts_amelia <- amelia(
        x = ts_imputed, 
        m = 19,
        p2s = 1,
        ts = 'Date',
        polytime = 3,
        splinetime = 3,
        # lags = setdiff(colnames(ts_imputed), 'Date'),
        # leads = setdiff(colnames(ts_imputed), 'Date'),
        logs = 
            colnames(ts_imputed)[
                as.logical(stri_count_fixed(colnames(ts_imputed), 'Close'))
                ],
        max.resample = Inf,
        parallel = 'multicore', 
        ncpus = 8)
    # 
    ts_imputed <- 
        Reduce('+', ts_amelia$imputations) / length(ts_amelia$imputations)
    ts_imputed <- ts_imputed[setdiff(colnames(ts_imputed), 'Date')]
    ts_imputed <- xts(ts_imputed,
                      order.by = as.Date(rownames(ts_imputed)))
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    rm(ts_amelia)
}



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
