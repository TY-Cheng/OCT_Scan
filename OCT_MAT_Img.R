setwd('/Users/chengt/Documents/OCT_Scan/')
load("/Users/chengt/Documents/OCT_Scan/Thickness_3D_Raw_.RData")
# source('/Users/chengt/Documents/OCT_Scan/OCT_function_utils.R')
# Import MAT file ---------------------------------------------------------
if(1){
    # MBR_1
    Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan/MBR_1_3D/MAT_files_/'
    setwd(Rootfolder)
    seq_MATfolder <- list.files() # sorted as char instead, not as num
    names(seq_MATfolder) <- seq_MATfolder
    Img_list_MBR1 <- llply(
        .data = seq_MATfolder, 
        .fun = Get_Thickness_3D_From_MAT,
        .progress = progress_time(),
        Rootfolder = Rootfolder
    )
    # MBR_2
    Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan/MBR_2_3D/MAT_files_/'
    setwd(Rootfolder)
    seq_MATfolder <- list.files() # sorted as char instead, not as num
    names(seq_MATfolder) <- seq_MATfolder
    Img_list_MBR2 <- llply(
        .data = seq_MATfolder, 
        .fun = Get_Thickness_3D_From_MAT,
        .progress = progress_time(),
        Rootfolder = Rootfolder
    )
}



# Process from image to thickness -----------------------------------------
# from 987 to 666
str(Img_list_MBR1)





# llply version
if (1) {
    print(Sys.time())
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    Thickness_list_MBR1 <- llply(
        .data = seq_along(Img_list_MBR1),
        .fun = function(i){
            Impute_Filter_Image_from_Matrix(
                Img_3D = Img_list_MBR1[[i]], 
                Fig_Title = df_MBR1$seq_Fig_Title[i],
                denoise_type = df_MBR1$seq_Denoise_Type[i],
                multiplier_pixel2micron = 2.1,
                flag_ex_extreme_value = T,
                quantiles = c(df_MBR1$seq_Quantile_Min, df_MBR1$seq_Quantile_Max)
            )
        },
        .parallel = T,
        # .progress = 'text',
        Impute_Filter_Image_from_Matrix = Impute_Filter_Image_from_Matrix,
        Img_list_MBR1 = Img_list_MBR1,
        df_MBR1 = df_MBR1
    )
    stopCluster(para_socket_cl)
}


# plot
llply(
    .data = seq_along(Thickness_list_MBR1), 
    .fun = function(i){
        source('/Users/chengt/Documents/OCT_Scan/OCT_function_utils.R')
        Plot_Thickness(
            Thickness = Thickness_list_MBR1[[i]], 
            scale_range = NULL,
            Fig_Title = df_MBR1$seq_Fig_Title[i],
            flag_save_plot = T,
            # save_folder = "/Volumes/Seagate_Backup/OCT_Scan/MBR_1_3D/Thickness_Img_/Default/"
            save_folder = '/Users/chengt/Documents/OCT_Scan/Img/Default/'
        )
    },
    .parallel = T,
    .progress = 'text',
    Plot_Thickness = Plot_Thickness,
    Thickness_list_MBR1 = Thickness_list_MBR1,
    df_MBR1 = df_MBR1
)

# test_list <- Thickness_list_MBR1[1:2]


stopCluster(para_socket_cl)

