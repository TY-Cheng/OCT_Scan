rm(list = ls())
setwd("/Users/chengt/Documents/OCT_Scan")
load("Thickness_3D_Raw_.RData")
getwd()

library(scales)
library(gstat)
library(sp)
# library(gapfill)
library(ggplot2)
library(colorRamps)


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




# Cropping, Denoising, Imaging -----------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Manual
    rm(list = ls())
    setwd("/Users/chengt/Documents/OCT_Scan")
    load("Thickness_3D_Raw_.RData")
    df_Img_3D <- data.frame(seq_Img_3D = c("Day_01_04.03_Resized_Img_3D_", "Day_06_09.03_Resized_Img_3D_", "Day_11_14.03_Resized_Img_3D_", "Day_16_19.03_Resized_Img_3D_", "Day_21_24.03_Resized_Img_3D_", "Day_23_26.03_Resized_Img_3D_", "Day_24_27.03_Resized_Img_3D_", "Day_25_28.03_Resized_Img_3D_", "Day_26_29.03_Resized_Img_3D_", "Day_27_30.03_Resized_Img_3D_", "Day_28_31.03_Resized_Img_3D_", "Day_29_01.04_Resized_Img_3D_", "Day_30_02.04_Resized_Img_3D_", "Day_32_04.04_Resized_Img_3D_", "Day_35_07.04_Resized_Img_3D_", "Day_37_09.04_Resized_Img_3D_", "Day_40_12.04_Resized_Img_3D_", "Day_42_14.04_Resized_Img_3D_", "Day_44_16.04_Resized_Img_3D_", "Day_46_18.04_Resized_Img_3D_", "Day_49_21.04_Resized_Img_3D_", "Day_52_24.04_Resized_Img_3D_", "Day_53_25.04_Resized_Img_3D_", "Day_56_28.04_Resized_Img_3D_", "Day_60_02.05_1_Resized_Img_3D_", "Day_60_02.05_10_Resized_Img_3D_", "Day_60_02.05_11_Resized_Img_3D_", "Day_60_02.05_12_Resized_Img_3D_", "Day_60_02.05_13_Resized_Img_3D_", "Day_60_02.05_14_Resized_Img_3D_", "Day_60_02.05_15_Resized_Img_3D_", "Day_60_02.05_2_Resized_Img_3D_", "Day_60_02.05_3_Resized_Img_3D_", "Day_60_02.05_4_Resized_Img_3D_", "Day_60_02.05_5_Resized_Img_3D_", "Day_60_02.05_6_Resized_Img_3D_", "Day_60_02.05_7_Resized_Img_3D_", "Day_60_02.05_8_Resized_Img_3D_", "Day_60_02.05_9_Resized_Img_3D_"), 
                            seq_denoise_type = c(1, 1, 1, 1, 1, #5
                                                 1, 1, 1, 1, 1, #10
                                                 1, 1, 1, 2, 2, #15
                                                 2, 2, 2, 2, 2, #20
                                                 2, 1, 1, 1, 1, #25
                                                 1, 1, 1, 1, 1, #30
                                                 1, 1, 1, 1, 1, #35
                                                 1, 1, 1, 1),
                            # seq_denoise_type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                            quantile_min = c(0.023, 0.023, 0.023, 0.023, 0.023, #5
                                             0.023, 0.023, 0.023, 0.023, 0.023, #10
                                             0.023, 0.023, 0.023, 0.023, 0.017, #15
                                             0.023, 0.023, 0.023, 0.037, 0.037, #20
                                             0.051, 0.023, 0.023, 0.023, 0.023, #25
                                             0.023, 0.023, 0.023, 0.023, 0.023, #30
                                             0.023, 0.023, 0.023, 0.023, 0.023, #35
                                             0.023, 0.023, 0.023, 0.023),
                            quantile_max = c(0.99, 0.99, 0.99, 0.99, 0.99, #5
                                             0.99, 0.99, 0.99, 0.99, 0.99, #10
                                             0.99, 0.99, 0.99, 0.99, 0.99, #15
                                             0.99, 0.99, 0.99, 0.99, 0.99, #20
                                             0.96, 0.99, 0.99, 0.99, 0.99, #25
                                             0.99, 0.99, 0.99, 0.99, 0.99, #30
                                             0.99, 0.99, 0.99, 0.99, 0.99, #35
                                             0.99, 0.99, 0.99, 0.99),
                            stringsAsFactors = F)
    df_Img_3D$seq_Fig_Title <- 
        paste0(
        gsub(pattern = 'Resized_Img_3D_', replacement = '', x = df_Img_3D$seq_Img_3D),
            c(rep('Stable_Flux', 5),
              rep('Relaxation_1', 9),
              rep('Relaxation_2', 4),
              rep('Air_Scouring', 3),
              rep('Relaxation_Air_Scouring', 3),
              rep('Autopsy', 15))
        )
    # df_Img_3D <- df_Img_3D[21, ]
    # df_Img_3D <- df_Img_3D[df_Img_3D$seq_denoise_type==2,]
    # 
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
                    quantiles = c(df_Img_3D$quantile_min[iter_Img_3D],
                                  df_Img_3D$quantile_max[iter_Img_3D]),
                    scale_range = scale_range, 
                    Fig_Title = df_Img_3D$seq_Fig_Title[iter_Img_3D],
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
                    quantiles = c(df_Img_3D$quantile_min[iter_Img_3D],
                                  df_Img_3D$quantile_max[iter_Img_3D]),
                    scale_range = scale_range,
                    Fig_Title = df_Img_3D$seq_Fig_Title[iter_Img_3D],
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
                    quantiles = c(df_Img_3D$quantile_min[iter_Img_3D],
                                  df_Img_3D$quantile_max[iter_Img_3D]),
                    scale_range = scale_range, 
                    Fig_Title = df_Img_3D$seq_Fig_Title[iter_Img_3D],
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
