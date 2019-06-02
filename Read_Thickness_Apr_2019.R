
rm(list = ls())
getwd()


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


# Copy&Paste by 'clipr' ---------------------------------------------------
if (1) {
    # assign & read by name in string
    # assign('x', 1)
    # eval(as.symbol('x')) 
    # copy&paste
    table <- clipr::read_clip_tbl(sep = ',')
    table
    clipr::write_clip(table)
}

# Read & Impute -----------------------------------------------------------
if(1) {    
    library(scales)
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
        library(Amelia)
        
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



# Stable Flux -------------------------------------------------------------
if(0){
    setwd('/Volumes/Seagate_Backup/OCT_Scan/MATfiles_Stable_Flux_')
    fname_1 <- list.files()
    # Day_01_04.03_Resized_, Day_06_09.03_Resized_, 
    # Day_11_14.03_Resized_, Day_16_19.03_Resized_, 
    # Day_21_24.03_Resized_
    for (iter_i in fname_1) {
        cat(iter_i)
        setwd('/Volumes/Seagate_Backup/OCT_Scan/MATfiles_Stable_Flux_')
        setwd(iter_i)
        fname_2 <- list.files()
        # Initialize the Thickness matrix
        Thickness <- matrix(data = NA, nrow = 666, ncol = 1970)
        # List of mat files in specific situation
        for (iter_j in fname_2) {
            # Load the mat file, sequentially
            Img_2D <- R.matlab::readMat(iter_j)
            Img_2D <- Img_2D$imageLayer[[3]]
            {
                # numerically, the result should be 'isos' - 'ilm'
                if (Img_2D[,1,1]$name =='isos') {
                    isos_y <- as.numeric(Img_2D[,1,1]$pathY)
                    isos_x <- as.numeric(Img_2D[,1,1]$pathX)
                }
                if (Img_2D[,1,2]$name =='ilm') {
                    ilm_y <- as.numeric(Img_2D[,1,2]$pathY)
                    ilm_x <- as.numeric(Img_2D[,1,2]$pathX)
                }
                # 
                # RowNumber: take the fourth 4th number for 'stable flux'
                row_num <- as.numeric(stringr::str_extract_all(
                    string = iter_j, 
                    pattern = '\\d+')[[1]][4]) + 1
                # 
                # ColumnNumber & extract the secant thickness
                cat('\n', iter_j)
                cat('\nisos & ilm match?\t', !sum(unique(isos_y)-unique(ilm_y)))
                
                thickness_row <- c()
                col_range <- unique(isos_y)
                for (col_num in col_range) {
                    # print(col_num)
                    if (col_num/max(col_range)<=.5) {
                        # print('less')
                        thickness_row[col_num] <- max(isos_x[isos_y==col_num]) - 
                            max(ilm_x[ilm_y==col_num])
                    }else{
                        # print('more')
                        thickness_row[col_num] <- min(isos_x[isos_y==col_num]) - 
                            min(ilm_x[ilm_y==col_num])
                    }
                }
                Thickness[row_num, col_range] <- thickness_row
            }
        }
        assign(x = paste0('Img_3D_', iter_i), value = Thickness)
    }
    
}


# Relaxation_1 ------------------------------------------------------------
if(0){
    setwd('/Volumes/Seagate_Backup/OCT_Scan/MATfiles_Relaxation_1_')
    # fname_1 <- "Day_32_04.04_Resized_"
    fname_1 <- list.files()
    # Day_23_26.03_Resized_, Day_24_27.03_Resized_, Day_25_28.03_Resized_, 
    # Day_26_29.03_Resized_, Day_27_30.03_Resized_, Day_28_31.03_Resized_, 
    # Day_29_01.04_Resized_, Day_30_02.04_Resized_, Day_32_04.04_Resized_
    for (iter_i in fname_1) {
        cat(iter_i)
        setwd('/Volumes/Seagate_Backup/OCT_Scan/MATfiles_Relaxation_1_')
        setwd(iter_i)
        fname_2 <- list.files()
        # Initialize the Thickness matrix
        Thickness <- matrix(data = NA, nrow = 666, ncol = 1970)
        # List of mat files in specific situation
        for (iter_j in fname_2) {
            # Load the mat file, sequentially
            Img_2D <- R.matlab::readMat(iter_j)
            Img_2D <- Img_2D$imageLayer[[3]]
            {
                # numerically, the result should be 'rpe' - 'isos'
                if (Img_2D[,1,1]$name =='isos') {
                    isos_y <- as.numeric(Img_2D[,1,1]$pathY)
                    isos_x <- as.numeric(Img_2D[,1,1]$pathX)
                }
                if (Img_2D[,1,2]$name =='ilm') {
                    ilm_y <- as.numeric(Img_2D[,1,2]$pathY)
                    ilm_x <- as.numeric(Img_2D[,1,2]$pathX)
                }
                # 
                # RowNumber: take the fourth 4th number for 'relaxation 1'
                row_num <- as.numeric(stringr::str_extract_all(
                    string = iter_j, 
                    pattern = '\\d+')[[1]][4]) + 1
                # 
                # ColumnNumber & extract the secant thickness
                cat('\n', iter_j)
                cat('\nisos & ilm match?\t', !sum(unique(isos_y)-unique(ilm_y)))
                
                thickness_row <- c()
                col_range <- unique(isos_y)
                for (col_num in col_range) {
                    # print(col_num)
                    if (col_num/max(col_range)<=.5) {
                        # print('less')
                        thickness_row[col_num] <- max(isos_x[isos_y==col_num]) - 
                            max(ilm_x[ilm_y==col_num])
                    }else{
                        # print('more')
                        thickness_row[col_num] <- min(isos_x[isos_y==col_num]) - 
                            min(ilm_x[ilm_y==col_num])
                    }
                }
                Thickness[row_num, col_range] <- thickness_row
            }
        }
        assign(x = paste0('Img_3D_', iter_i), value = Thickness)
    }
    
}

# save.image("/Volumes/Seagate_Backup/OCT_Scan/Thickness_3D_Raw_.RData")

# Cropping & Denoising for Img_3D -----------------------------------------
if (1) {
    rm(list = ls())
    library(ggplot2)
    library(colorRamps)
    
    setwd('/Volumes/Seagate_Backup/OCT_Scan')
    load("Thickness_3D_Raw_.RData")
    fname <- c('Img_3D_Day_01_04.03_Resized_', 'Img_3D_Day_06_09.03_Resized_', 
               'Img_3D_Day_11_14.03_Resized_', 'Img_3D_Day_16_19.03_Resized_', 
               'Img_3D_Day_21_24.03_Resized_', 
               'Img_3D_Day_23_26.03_Resized_', 'Img_3D_Day_24_27.03_Resized_', 
               'Img_3D_Day_25_28.03_Resized_', 'Img_3D_Day_26_29.03_Resized_', 
               'Img_3D_Day_27_30.03_Resized_', 'Img_3D_Day_28_31.03_Resized_', 
               'Img_3D_Day_29_01.04_Resized_', 'Img_3D_Day_30_02.04_Resized_', 
               'Img_3D_Day_32_04.04_Resized_')
    # Proved that col_number is '987' for all rows!
    if(1){
        for (iter_i in fname) {
            print(iter_i)
            Thickness <- eval(as.symbol(iter_i))
            num_notna <- apply(Thickness, MARGIN = 1, function(x)sum(!is.na(x)))
            print(sum((num_notna == 0 )|(num_notna == 987)))
            rm(num_notna)
        }
    }
    
    # All except 'Img_3D_Day_32_04.04_Resized_'
    for (iter_i in fname[-14]) {
        cat('\n', iter_i, '\n', 'Denoising...\n')
        Thickness <- na.omit(eval(as.symbol(iter_i))[,1:987])
        print(dim(Thickness))
        Thickness <- EBImage::resize(Thickness, 
                                     w = 666, 
                                     h = 666)
        # times 2.1, from pixels to microns
        Thickness <- Thickness * 2.1
        # Denoising
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(Thickness/Thickness_max,
                                           size = exp(1)) * Thickness_max
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(Thickness/Thickness_max,
                                           size = pi) * Thickness_max
        # To constrain extreme values
        Thickness_wo_Extreme <- Thickness
        minimum <- quantile(Thickness, probs = .01)
        maximum <- quantile(Thickness, probs = .99)
        Thickness_wo_Extreme[Thickness_wo_Extreme<minimum] <- minimum
        Thickness_wo_Extreme[Thickness_wo_Extreme>maximum] <- maximum
        # 
        print(range(Thickness_wo_Extreme))
        # Thickness <- waveslim::denoise.dwt.2d(Thickness)
        # # image(Thickness)
        assign(paste0(iter_i, 'Denoised_'), 
               Thickness_wo_Extreme
               # waveslim::denoise.dwt.2d(Thickness)
        )
        
        # To plot
        x <- seq_len(ncol(Thickness))
        y <- seq_len(nrow(Thickness))
        Img_3D_grid <- expand.grid(X = x, Y = y)
        Img_3D_grid$Z <- c(t(eval(as.symbol(paste0(iter_i, 'Denoised_')))))
        
        # Levelplot with ggplot2
        ggplot(Img_3D_grid, aes(x = X, y = Y, z = Z)) +
            geom_raster(aes(fill = Z)) +
            scale_fill_gradientn(
                limits = c(
                    max(round(median(Img_3D_grid$Z)-20, digits = -1), 0),
                    max(round(median(Img_3D_grid$Z)-20, digits = -1), 0) + 50
                ),
                colours = matlab.like(67)
                # low = "gray89", high = "gray2",
                # limits = c(0, 150)
            ) +
            # stat_contour(breaks = seq(30, 100, 30), colour = 'white', na.rm = T) +
            coord_fixed() +
            labs(title = paste0(iter_i, 'Denoised_'),
                 x = '', y = '',
                 fill = expression(paste('Thickness (',mu,'m)'))
            ) +
            theme_bw()
        
        cat('Saving...\n')
        ggsave(paste0(
            'Img_3D_Thickness_StableFlux(Day_1_21)_Relaxation1(Day_23_32)_',
            'Denoised_/',
            iter_i,'Denoised_.png'), width = 16)
        rm(list = c(iter_i, 'minimum', 'maximum', 'x', 'y',
                    'Thickness', 'Thickness_max', 'Thickness_wo_Extreme', 
                    'Img_3D_grid'))
    }
    
    # For'Img_3D_Day_32_04.04_Resized_' only
    if(1){
        iter_i <- fname[14]
        Thickness <- eval(as.symbol(iter_i))
        {
            num_notna <- apply(Thickness, MARGIN = 1, function(x)sum(!is.na(x)))
            print(sum((num_notna == 0 )|(num_notna == 987)))
            plot(rowMeans(Thickness, na.rm = T))
            abline(h = mean(Thickness, na.rm = T), col = 'red')
            rm(num_notna)
        }
        
        cat('\n', iter_i, '\n', 'Denoising...\n')
        
        Thickness[rowMeans(Thickness, na.rm = T) < 
                      mean(Thickness, na.rm = T),] <- NA
        # Thickness[rowMeans(Thickness, na.rm = T) < 
        #               mean(Thickness, na.rm = T),] <- Thickness[rowMeans(
        #                   Thickness, na.rm = T) < mean(Thickness, na.rm = T), ] +
        #     mean(Thickness[rowMeans(Thickness, na.rm = T) > mean(Thickness, na.rm = T),]) -
        #     mean(Thickness[rowMeans(Thickness, na.rm = T) < mean(Thickness, na.rm = T),])
        Thickness <- na.omit(Thickness[,1:987])
        print(dim(Thickness))
        Thickness <- EBImage::resize(Thickness, 
                                     w = 666, 
                                     h = 666)
        # times 2.1, from pixels to microns
        Thickness <- Thickness * 2.1
        # Denoising
        # Thickness_max <- max(Thickness)
        # Thickness <- EBImage::medianFilter(Thickness/Thickness_max,
        #                                    size = 2) * Thickness_max
        Thickness <- waveslim::denoise.dwt.2d(Thickness)
        Thickness_max <- max(Thickness)
        Thickness <- EBImage::medianFilter(Thickness/Thickness_max,
                                           size = pi) * Thickness_max
        # To constrain extreme values
        Thickness_wo_Extreme <- Thickness
        minimum <- quantile(Thickness, probs = .01)
        maximum <- quantile(Thickness, probs = .99)
        Thickness_wo_Extreme[Thickness_wo_Extreme<minimum] <- minimum
        Thickness_wo_Extreme[Thickness_wo_Extreme>maximum] <- maximum
        # 
        print(range(Thickness_wo_Extreme))
        # # image(Thickness)
        assign(paste0(iter_i, 'Denoised_'), 
               Thickness_wo_Extreme
               # waveslim::denoise.dwt.2d(Thickness)
        )
        
        # To plot
        x <- seq_len(ncol(Thickness))
        y <- seq_len(nrow(Thickness))
        Img_3D_grid <- expand.grid(X = x, Y = y)
        Img_3D_grid$Z <- c(t(eval(as.symbol(paste0(iter_i, 'Denoised_')))))
        
        # Levelplot with ggplot2
        ggplot(Img_3D_grid, aes(x = X, y = Y, z = Z)) +
            geom_raster(aes(fill = Z)) +
            scale_fill_gradientn(
                limits = c(
                    max(round(median(Img_3D_grid$Z)-20, digits = -1), 0),
                    max(round(median(Img_3D_grid$Z)-20, digits = -1), 0) + 50
                ),
                colours = matlab.like(67)
                # low = "gray89", high = "gray2",
                # limits = c(0, 150)
            ) +
            # stat_contour(breaks = seq(30, 100, 30), colour = 'white', na.rm = T) +
            coord_fixed() +
            labs(title = paste0(iter_i, 'Denoised_'),
                 x = '', y = '',
                 fill = expression(paste('Thickness (',mu,'m)'))
            ) +
            theme_bw()
        
        cat('Saving...\n')
        ggsave(paste0(
            'Img_3D_Thickness_StableFlux(Day_1_21)_Relaxation1(Day_23_32)_',
            'Denoised_/',
            iter_i,'Denoised_.png'), width = 16)
        rm(list = c(iter_i, 'minimum', 'maximum', 'x', 'y',
                    'Thickness', 'Thickness_max', 'Thickness_wo_Extreme', 
                    'Img_3D_grid'))
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