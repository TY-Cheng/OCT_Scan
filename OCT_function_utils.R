library(doParallel)
library(plyr)
library(dplyr)
library(mgcv)
library(EBImage)
library(ggplot2)
library(colorRamps)
library(stringr)



# Hyperparameters ---------------------------------------------------------
if (0) {
    df_MBR1 <- data.frame(
        seq_Fig_Index = names(Img_list_MBR1),
        seq_Fig_Title = paste0(
            names(Img_list_MBR1), 
            c(rep('Stable_Flux', 5), 
              rep('Relaxation_1', 9), 
              rep('Relaxation_2', 4), 
              rep('Air_Scouring', 3),
              rep('Relaxation_Air_Scouring', 3),
              rep('Autopsy', 15))
        ),
        seq_Denoise_Type = c(rep(1, 13), rep(2, 8), rep(1, 18)),
        seq_Loess_Span = c(rep(.02, 18), .02, rep(.02, 20)),
        seq_Quantile_Min = c(rep(.01, 14), .01, 
                             rep(.01, 3), 
                             .01, # 44
                             .01, # 46
                             .01, # 49
                             rep(.01, 18)),
        seq_Quantile_Max = c(
            rep(.999, 13), 
            rep(.999, 5),
            rep(.995, 3),
            rep(.999, 18)
        ),
        seq_mean_aim = c(
            48.41, 73.01, 94.16, 127.06, 131.38, 
            133.6, 139.16, 141.47, 129.87, 137.8, 144.39, 144.39, 164.19, 171.25, 
            239.37, 220.3, 230.05, 230.23, 
            176.18, 131.54, 107.82,
            48.39, 39.67, 32.08,
            41.85, 33.35, 38.77, 49.99, 54.45, 
            47.88, 33.69, 31.8, 37.57, 43.06, 
            40.99, 35.85, 31.02, 33.57, 39.21
        ),
        seq_sd_aim = c(
            5.14, 9.68, 8.05, 10.53, 11.71,
            12.08, 11.99, 15.49, 13.31, 14.39, 17.77, 17.77, 14.52, 16.79,
            13.7, 8.73, 9.86, 11.89,
            11.57, 10.46, 13.12,
            11.81, 13.31, 13.4,
            15.15, 10.27, 15.69, 13.2, 15.43, 
            14.79, 18.48, 14.23, 13.76, 16.5, 
            15.54, 19.99, 11.57, 13.72, 12.87
        ),
        stringsAsFactors = F
    )
    # df_MBR2 needs sorting and then add mean&sd aim.
    df_MBR2 <- data.frame(
        seq_Fig_Index = names(Img_list_MBR2),
        seq_Fig_Title = paste0(
            c(
                rep('Reactor_1_Relaxation_', 18),
                rep('Reactor_2_Air_Scouring_', 11),
                rep('Reactor_3_Relaxation_Air_Scouring_', 11)
            ),
            c(
                rep('Day_', 15), rep('Day_40_Autopsy_', 3),
                rep('Day_', 8), rep('Day_40_Autopsy_', 3),
                rep('Day_', 8), rep('Day_40_Autopsy_', 3)
            ),
            stringr::str_pad(
                stringi::stri_extract_last_regex(names(Img_list_MBR2), '(\\d)+'), 
                width = 2, pad = '0')
        ),
        seq_Denoise_Type = c(
            rep(2,4), # 4
            rep(3,8), # Reactor_1_Relaxation_Day_5 6 7 8 9 12 14 20
            2,2, # Reactor_1_Relaxation_Day_26, 34
            3,3,
            2,3,   # Reactor_1_Relaxation_Autopsy_2, 3
            rep(2,3), # Reactor_2_Air_Scouring_Day_1, 4, 8
            rep(3,4), # Reactor_2_Air_Scouring_Day_12, 17, 26, 34
            2,    # Reactor_2_Air_Scouring_Day_39, 
            3,    # Reactor_2_Air_Scouring_Autopsy_1
            2,    # Reactor_2_Air_Scouring_Autopsy_2
            3,    # Reactor_2_Air_Scouring_Autopsy_3
            rep(2,5), # Reactor_3_Relaxation_Air_Scouring_Day_1, 4, 8, 12, 17
            3,    # Reactor_3_Relaxation_Air_Scouring_Day_26
            2,2,  # Reactor_3_Relaxation_Air_Scouring_Day_34, 39
            2,2,3 # Reactor_3_Relaxation_Air_Scouring_Autopsy_1, 2, 3
        ),
        seq_Loess_Span = c(rep(.02, 40)),
        seq_Quantile_Min = c(rep(.005, 40)),
        seq_Quantile_Max = c(rep(.998, 40)),
        stringsAsFactors = F
    )
    df_MBR2 <- df_MBR2[order(df_MBR2$seq_Fig_Title),] # sort according to fig title
    df_MBR2$seq_mean_aim <- c(
        17.45,48.93,58.04,70.69,80.1,86.88,92.81,98.5,105.59,130.98,139.27,191.94,226.54,280.36,281.91,281.91,278.47,290.83,16.25,41.26,47.04,38.28,37.49,47.15,51.44,46.87,46.87,54.47,66.58,27.31,19.85,33.26,30.93,56.72,53.05,60.67,40.62,40.62,53.37,58.95
    )
    df_MBR2$seq_sd_aim <- c(
        5.1,8.22,8.08,7.49,7.36,7.75,8.59,7.59,7.07,7.42,7.97,9.21,9,8.44,6.82,6.82,7.62,9.23,7.15,9.8,6.58,7.02,12.45,8.03,7.37,8.31,8.31,7.07,11.62,11.77,10.76,17.1,18.54,7.63,16.67,7.3,19.76,19.76,8.62,6.5
    )
    # 
    
}

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
Get_Thickness_3D_From_MAT <- function(
    Working_Directory,
    Rootfolder = Rootfolder,
    n_row = 666,
    n_col = 1970
) {
    print(Working_Directory)
    setwd(Rootfolder)
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
    Img_3D <- matrix(NA, nrow = n_row, ncol = n_col)
    invisible(lapply(
        X = Img_2D_list,
        FUN = function(rownum_thicknessrow){
            Img_3D[rownum_thicknessrow$row_num,
                   1:length(rownum_thicknessrow$thickness_row)] <<- 
                rownum_thicknessrow$thickness_row
        }
    ))
    Img_3D <- Img_3D[, !colSums(is.na(Img_3D))==n_row]
    return(Img_3D)
}


# Remove consecutive NAs --------------------------------------------------
Remove_Consecutive_NA <- function(Thickness) {
    fit_ref <- rowMeans(Thickness)
    index_na <- rle(is.na(fit_ref))
    index_na_rm = which(index_na$values == TRUE & index_na$lengths > 2)
    if (any(index_na_rm)) {
        index_na_rm_end <- cumsum(index_na$lengths)[index_na_rm]
        index_na_rm_start <- cumsum(index_na$lengths)[ifelse(
            index_na_rm > 1, index_na_rm - 1, 0
        )] + 1
        if (0 %in% ifelse(index_na_rm > 1, index_na_rm - 1, 0)) {
            index_na_rm_start <-  c(1, index_na_rm_start)
        }
        # 
        index_na_rm <- c()
        for (iter_i in seq_along(index_na_rm_start)) {
            index_na_rm <- c(index_na_rm, 
                             index_na_rm_start[iter_i]:(index_na_rm_end[iter_i]-1))
        }
        Thickness <- Thickness[!seq_len(NROW(Thickness))%in%index_na_rm, ]
    }
    return(Thickness)
}

# Process matrix to image -------------------------------------------------
# from 987 to n_row (666)
Impute_Filter_Image_from_Matrix <- function(
    Img_3D, 
    Fig_Title,
    denoise_type = 1,
    multiplier_pixel2micron = 2.1,
    flag_ex_extreme_value = TRUE,
    loess_span = .02,
    quantiles = c(.027, .97),
    Remove_Consecutive_NA = Remove_Consecutive_NA,
    n_row = 666
) {
    # print(paste0('Processing ', Fig_Title, '...'))
    Thickness <- Img_3D
    # From pixels to microns
    Thickness <- Thickness * multiplier_pixel2micron
    # rm(Img_3D)
    # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Cropping for denoise_type 2
    if (denoise_type %in% 2:3) {
        fit_ref <- rowMeans(Thickness)
        fit_loess <- loess(y~x, 
                           data.frame(x = 1:n_row, y = fit_ref),
                           na.action = na.omit,
                           span = .6)
        threshold_row <- 
            mean(Thickness, na.rm = T) + predict(fit_loess, 1:n_row) - 2 * fit_ref
        Thickness[fit_ref < threshold_row, ] <- NA
        # threshold_row <- predict(fit_loess, 1:n_row) + fit_ref
        # Thickness[fit_ref > threshold_row, ] <- NA
        # 
        fit_ref <- rowMeans(Thickness)
        fit_loess <- loess(y~x, 
                           data.frame(x = 1:n_row, y = fit_ref),
                           na.action = na.omit,
                           span = 1/2)
        threshold_row <- 
            .6 * mean(Thickness, na.rm = T) + 
            .4 * predict(fit_loess, 1:n_row) - 
            sd(fit_ref, na.rm = T) * 1.732
        Thickness[fit_ref < threshold_row, ] <- NA
        threshold_row <- predict(fit_loess, 1:n_row) + sd(fit_ref, na.rm = T) * sqrt(5)
        Thickness[fit_ref > threshold_row, ] <- NA
    }
    # rm too much NA
    Thickness <- Remove_Consecutive_NA(Thickness = Thickness)
    # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    if (flag_ex_extreme_value) {
        # To exclude extreme values
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^4, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-1/10^4, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
        # rm(Thickness_wo_Extreme, maximum, minimum)
    }
    # Specials
    # quantile(Img_list_MBR1$Day_49_21.04_, na.rm = T, 
    # probs = seq(from = 0 , to = 1, by = .01))
    if (Fig_Title %in% 
        c(
            'Day_32_04.04_Relaxation_1', 'Day_35_07.04_Relaxation_2',
            'Day_37_09.04_Relaxation_2', 'Day_40_12.04_Relaxation_2',
            'Day_42_14.04_Relaxation_2'
        )
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^2, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-1/10^4, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    if (Fig_Title %in% c('Day_44_16.04_Air_Scouring', 'Day_46_18.04_Air_Scouring')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^2, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-5/10^3, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    if (Fig_Title %in% c('Day_49_21.04_Air_Scouring')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^2, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-1/10^5, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    if (Fig_Title %in% c('Reactor_1_Relaxation_Day_3',
                         'Reactor_1_Relaxation_Day_4')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 1/10^2, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-2/10^3, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    if (Fig_Title %in% c('Reactor_1_Relaxation_Day_34')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 3/10^2, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-5/10^2, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
        Thickness <- Remove_Consecutive_NA(Thickness)
    }
    if (Fig_Title %in% c('Reactor_2_Air_Scouring_Day_4',
                         'Reactor_2_Air_Scouring_Day_12')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^3, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-1/10^2, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    if (Fig_Title %in% c('Reactor_2_Air_Scouring_Day_34', 
                         'Reactor_2_Air_Scouring_Autopsy_2')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^3, na.rm = T), 
            0)
        maximum <- quantile(Thickness, probs = 1-3/10^2, na.rm = T)
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
    }
    if (Fig_Title %in% c('Reactor_3_Relaxation_Air_Scouring_Autopsy_3')
    ) {
        Thickness_wo_Extreme <- Thickness
        minimum <- max(
            quantile(Thickness, probs = 5/10^3, na.rm = T), 
            0)
        maximum <- 40
        Thickness_wo_Extreme[Thickness_wo_Extreme < minimum] <- NA
        Thickness_wo_Extreme[Thickness_wo_Extreme > maximum] <- NA
        Thickness <- Thickness_wo_Extreme
        Thickness <- Remove_Consecutive_NA(Thickness)
    }
    
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Imputation
    if (sum(is.na(Thickness))>0) {
        Img_3D_grid <- expand.grid(
            X = seq_len(NCOL(Thickness)),
            Y = seq_len(NROW(Thickness))
        )
        Img_3D_grid$Z <- c(t(Thickness))
        Img_3D_grid$X <- Img_3D_grid$X
        Img_3D_grid$Y <- Img_3D_grid$Y
        # loess
        fit_loess <- loess(
            Z~., Img_3D_grid, 
            span = loess_span,
            na.action = na.omit,
            degree = 2, family = 'symmetric',
            control = loess.control(
                trace.hat = 'approximate'
            )
        )
        # 
        # GAM
        # fit_gam <- mgcv::gam(
        #     Z ~ te(X, bs = 'cr', k = 30, m = 30) + 
        #         te(Y, bs = 'cr', k = 30, m = 30) + 
        #         t2(X, Y, bs = c('ds', 'ds')),
        #     # Z ~ te(X, bs = 'cc') + te(Y, bs = 'cc') + s(X, Y, bs = 'sos'),
        #     data = Img_3D_grid, na.action = na.omit,
        #     # method = 'REML',
        #     control = list(nthreads = parallel::detectCores()), 
        #     family = gaussian()
        # )
        # vis.gam(fit_gam, n.grid = 600, plot.type = 'contour')
        # summary(fit_gam)
        # 
        index_na <- is.na(Thickness)
        fit_loess <- t(predict(fit_loess, newdata = Img_3D_grid))
        Thickness[index_na] <- fit_loess[index_na]
        # 
        # Img_3D_grid[index_na, 'Z'] <- predict(fit_gam, Img_3D_grid[index_na, ])
        # Thickness <- xtabs(Z ~ X + Y, data = Img_3D_grid)
        # 
        # Resize to 666*666
        # Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
        # rm(Img_3D_grid)
    }
    
    # a final assurance
    Thickness <- na.omit(Thickness)
    # Resize to 666*666
    Thickness <- EBImage::resize(Thickness, w = 666, h = 666)
    # Filter
    if (denoise_type %in% 1:2) {
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        # rm(Thickness_max)
    }
    if (denoise_type %in% 3) {
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness <- waveslim::denoise.modwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = 1) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(
            Thickness/Thickness_max, size = pi) * Thickness_max
        
    }
    if (flag_ex_extreme_value) {
        # To exclude extreme values
        Thickness_wo_Extreme <- Thickness
        minimum <- max(quantile(Thickness, probs = quantiles[1]), 
                       0)
        maximum <- quantile(Thickness, probs = quantiles[2])
        Thickness_wo_Extreme[Thickness_wo_Extreme<minimum] <- minimum
        Thickness_wo_Extreme[Thickness_wo_Extreme>maximum] <- maximum
        Thickness <- Thickness_wo_Extreme
        # rm(Thickness_wo_Extreme, maximum, minimum)
    }
    return(Thickness)
}


# Imaging from Thickness Matrix -------------------------------------------
Plot_Thickness <- function(Thickness, 
                           scale_range = 60,
                           Fig_Title,
                           flag_plot = F,
                           x_in_mm = 4,
                           y_in_mm = 4,
                           flag_save_plot = T,
                           save_folder = '') {
    print(paste0('Processing ', Fig_Title, '...'))
    # 16 "Look Up Table" or LUT (a.k.a. palette, color table)
    lut_16 <- matrix(
        data = c(0,0,0,
                 1,1,171,
                 1,1,224,
                 0,110,255,
                 1,171,254,
                 1,224,254,
                 1,254,1,
                 190,255,0,
                 255,255,0,
                 255,224,0,
                 255,141,0,
                 250,94,0,
                 245,0,0,
                 245,0,172,
                 222,180,222),
        ncol = 3, byrow = T
    )
    lut_16 <- rgb(lut_16, maxColorValue = 255)
    # ReExpand the Grid to plot 666*666 as 4 mm * 4 mm
    Img_3D_grid <- expand.grid(
        X = seq_len(NCOL(Thickness)),
        Y = seq_len(NROW(Thickness))
    )
    Img_3D_grid$Z <- c(Thickness)
    Img_3D_grid$X <- Img_3D_grid$X / 666 * x_in_mm
    Img_3D_grid$Y <- Img_3D_grid$Y / 666 * y_in_mm
    # Levelplot with ggplot2
    p <- ggplot2::ggplot(Img_3D_grid, ggplot2::aes(x = X, y = Y, z = Z)) +
        ggplot2::geom_raster(ggplot2::aes(fill = Z)) +
        # stat_contour(bins = 15, colour = 'white') +
        ggplot2::coord_fixed() +
        ggplot2::labs(
            title = Fig_Title,
            x = 'X (mm)',
            y = 'Y (mm)',
            fill = paste('Thickness\n(\u03BCm)', sep = '')
        ) +
        ggplot2::theme_bw()
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
            ggplot2::scale_fill_gradientn(
                colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                )
            )
    }else if (length(scale_range) == 2) {
        p <- p +
            ggplot2::scale_fill_gradientn(
                limits = c(scale_range[1], scale_range[2]),
                colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                )
            )
    }else if (length(scale_range) == 1) {
        p <- p +
            ggplot2::scale_fill_gradientn(
                limits = c(
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0),
                    max(round(min(Img_3D_grid$Z) - 20, digits = -1), 0) + 
                        scale_range
                ),
                colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                )
            )
    }
    # 
    if (flag_plot) {print(p)}
    # 
    if (flag_save_plot) {
        cat('Saving', paste0(save_folder, Fig_Title, '.png'), '...\n')
        ggplot2::ggsave(filename = paste0(save_folder, Fig_Title, '.png'), 
                        plot = p, width = 7, height = 6.5)
    }
    return()
}


# Calibration of Thickness ------------------------------------------------
Calibrate_Thickness <- function(Thickness_raw, mean_aim, sd_aim) {
    Thickness <- Thickness_raw
    Thickness <- (Thickness-mean(Thickness))/sd(Thickness)
    Thickness <- Thickness * sd_aim
    Thickness <- Thickness + mean_aim
    return(Thickness)
}

