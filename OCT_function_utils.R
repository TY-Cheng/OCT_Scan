library(doParallel)
library(plyr)
library(mgcv)
library(EBImage)
# library(scales)
# library(gstat)
# library(sp)
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
Get_Thickness_3D_From_MAT <- function(Working_Directory) {
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
    Img_3D <- Img_3D[, !colSums(is.na(Img_3D))==666]
    setwd(old_wd)
    tictoc::toc()
    return(Img_3D)
}


# Process matrix to image -------------------------------------------------
Impute_Filter_Image_from_Matrix <- function(
    Img_3D, 
    denoise_type = 1,
    multiplier_pixel2micron = 2.1,
    flag_ex_extreme_value = TRUE,
    quantiles = c(.027, .97)
) {
    Thickness <- Img_3D
    # From pixels to microns
    Thickness <- Thickness * multiplier_pixel2micron
    rm(Img_3D)
    
    # Cropping for denoise_type 2
    if (denoise_type == 2) {
        fit_ref <- rowMeans(Thickness)
        fit_loess <- loess(y~x, 
                           data.frame(x = 1:666, y = fit_ref),
                           na.action = na.omit,
                           span = 1/2)
        threshold_row <- predict(fit_loess, 1:666) - fit_ref/7
        Thickness[fit_ref < threshold_row, ] <- NA
        # 
        fit_loess <- loess(y~x, 
                           data.frame(x = 1:666, y = fit_ref),
                           na.action = na.omit,
                           span = 1/3)
        threshold_row <- predict(fit_loess, 1:666) - sd(fit_ref, na.rm = T) * 2.236
        Thickness[fit_ref < threshold_row, ] <- NA
    }
    
    # Imputation, rescale xy ranges as 4 mm
    if (1) {
        Img_3D_grid <- expand.grid(
            X = seq_len(NCOL(Thickness)),
            Y = seq_len(NROW(Thickness))
        )
        Img_3D_grid$Z <- c(t(Thickness))
        Img_3D_grid$X <- Img_3D_grid$X
        Img_3D_grid$Y <- Img_3D_grid$Y
        # # loess
        # fit_loess <- loess(Z~., Img_3D_grid, span = .1, na.action = na.omit)
        # 
        # GAM
        fit_gam <- mgcv::gam(
            Z ~ ti(X, bs = 'ts') + ti(Y, bs = 'ts') + t2(X, Y, bs = c('tp', 'tp')),
            data = Img_3D_grid, na.action = na.omit
        )
        # vis.gam(fit_gam)
        index_gam <- is.na(Img_3D_grid$Z)
        Img_3D_grid[index_gam, 'Z'] <- predict(fit_gam, Img_3D_grid[index_gam, ])
        invisible(
            apply(Img_3D_grid, MARGIN = 1,
                  function(row_vec){Thickness[row_vec[2], row_vec[1]] <<- row_vec[3]})
        )
        # Resize to 666*666
        Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
        rm(Img_3D_grid)
    }
    # Filter
    if (denoise_type %in% 1:2) {
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        rm(Thickness_max)
    }
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
    return(Thickness)
}


# Imaging from Thickness Matrix -------------------------------------------
Plot_Thickness <- function(Thickness, 
                           scale_range = 60,
                           Fig_Title,
                           flag_save_plot= T,
                           save_folder = '') {
    # ReExpand the Grid to plot 666*666 as 4 mm * 4 mm
    Img_3D_grid <- expand.grid(
        X = seq_len(NCOL(Thickness)),
        Y = seq_len(NROW(Thickness))
    )
    Img_3D_grid$Z <- c(t(Thickness))
    Img_3D_grid$X <- Img_3D_grid$X / 666 * 4
    Img_3D_grid$Y <- Img_3D_grid$Y / 666 * 4
    # Levelplot with ggplot2
    p <- ggplot(Img_3D_grid, aes(x = X, y = Y, z = Z)) +
        geom_raster(aes(fill = Z)) +
        coord_fixed() +
        labs(
            # title = Fig_Title,
            x = 'X (mm)',
            y = 'Y (mm)',
            fill = paste('Thickness\n(\u03BCm)', sep = '')
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
    # if (flag_plot) {print(p)}
    # 
    if (flag_save_plot) {
        cat('Saving to', paste0(save_folder, Fig_Title, '.png'), '...\n')
        ggsave(filename = paste0(save_folder, Fig_Title, '.png'), 
               plot = p, width = 7)
    }
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
