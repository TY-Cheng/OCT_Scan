if (0) {
    setwd('/Users/chengt/Documents/OCT_Scan/')
    load("/Users/chengt/Documents/OCT_Scan/OCT_Thickness.RData")
    source('/Users/chengt/Documents/OCT_Scan/OCT_function_utils.R')
}

# Import MAT file ---------------------------------------------------------
if(0){
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


# Processing and Imaging --------------------------------------------------
# from 987 to 666

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Interpolation & Filtering & Plotting , MBR1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# str(Img_list_MBR1)
if (0) {
    tictoc::tic('Interpolation & Filtering, MBR1')
    seq_list <- as.list(as.data.frame(t(df_MBR1), stringsAsFactors = F))
    names(seq_list) <- df_MBR1$seq_Fig_Index
    # 
    # seq_list <- seq_list[c(14:21,36)]
    # 
    print(Sys.time())
    print(names(seq_list))
    # 
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # llply version
    Thickness_list_MBR1 <- llply(
        .data = seq_list,
        .fun = function(
            iter_i, Img_list, Img_Processing, NA_Action
        ){
            return(
                Img_Processing(
                    Img_3D = Img_list[[iter_i[1]]], 
                    Fig_Title = iter_i[2],
                    denoise_type = as.numeric(iter_i[3]),
                    multiplier_pixel2micron = 2.1,
                    flag_ex_extreme_value = T,
                    loess_span = as.numeric(iter_i[4]),
                    quantiles = c(as.numeric(iter_i[5]), as.numeric(iter_i[6])),
                    Remove_Consecutive_NA = NA_Action
                )
            )
        },
        .inform = T,
        .parallel = T,
        # .progress = 'time',
        Img_list = Img_list_MBR1,
        Img_Processing = Impute_Filter_Image_from_Matrix,
        NA_Action = Remove_Consecutive_NA
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}
# GGplot llply version, MBR1
if (1) {
    tictoc::tic('Plotting, MBR1')
    seq_list <- seq_along(Thickness_list_MBR1)
    names(seq_list) <- names(Thickness_list_MBR1[seq_list])
    # 
    print(Sys.time())
    print(names(seq_list))
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # 
    invisible(
        llply(
            .data = seq_list,
            .fun = function(
                i, 
                Plot_Thickness, 
                Thickness_list, 
                df_hyper
            ){
                Plot_Thickness(
                    Thickness = Thickness_list[[i]], 
                    scale_range = NULL,
                    Fig_Title = df_hyper$seq_Fig_Title[i],
                    flag_plot = F,
                    flag_save_plot = T,
                    # save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Default/'
                    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/new/'
                )
            },
            .parallel = T,
            # .progress = 'time',
            Plot_Thickness = Plot_Thickness,
            Thickness_list = Thickness_list_MBR1,
            df_hyper = df_MBR1
        )
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Interpolation & Filtering & Plotting , MBR2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# str(Img_list_MBR2)
if (0) {
    tictoc::tic('Interpolation & Filtering, MBR2')
    seq_list <- as.list(as.data.frame(t(df_MBR2), stringsAsFactors = F))
    names(seq_list) <- df_MBR2$seq_Fig_Index
    # 
    # seq_list <- seq_list[c(14:21,36)]
    # 
    print(Sys.time())
    print(names(seq_list))
    # 
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # llply version
    Thickness_list_MBR2 <- llply(
        .data = seq_list,
        .fun = function(
            iter_i, Img_list, Img_Processing, NA_Action
        ){
            return(
                Img_Processing(
                    Img_3D = Img_list[[iter_i[1]]], 
                    Fig_Title = iter_i[2],
                    denoise_type = as.numeric(iter_i[3]),
                    multiplier_pixel2micron = 2.1,
                    flag_ex_extreme_value = T,
                    loess_span = as.numeric(iter_i[4]),
                    quantiles = c(as.numeric(iter_i[5]), as.numeric(iter_i[6])),
                    Remove_Consecutive_NA = NA_Action
                )
            )
        },
        .inform = T,
        .parallel = T,
        # .progress = 'time',
        Img_list = Img_list_MBR2,
        Img_Processing = Impute_Filter_Image_from_Matrix,
        NA_Action = Remove_Consecutive_NA
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}
# GGplot llply version, MBR2
if (1) {
    tictoc::tic('Plotting, MBR2')
    seq_list <- seq_along(Thickness_list_MBR2)
    names(seq_list) <- names(Thickness_list_MBR2[seq_list])
    # 
    print(Sys.time())
    print(names(seq_list))
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # 
    invisible(
        llply(
            .data = seq_list,
            .fun = function(
                i, 
                Plot_Thickness, 
                Thickness_list, 
                df_hyper
            ){
                Plot_Thickness(
                    Thickness = Thickness_list[[i]], 
                    scale_range = NULL,
                    Fig_Title = df_hyper$seq_Fig_Title[i],
                    flag_plot = F,
                    flag_save_plot = T,
                    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_2/Default/'
                )
            },
            .parallel = T,
            # .progress = 'time',
            Plot_Thickness = Plot_Thickness,
            Thickness_list = Thickness_list_MBR2,
            df_hyper = df_MBR2
        )
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}


# Calibration and Imaging --------------------------------------------------

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Calibrating mean&sd & Plotting, MBR1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# str(Thickness_list_MBR1)
if (0) {
    tictoc::tic('Calibrating, MBR1')
    seq_list <- seq_along(Thickness_list_MBR1)
    names(seq_list) <- names(Thickness_list_MBR1[seq_list])
    # 
    # seq_list <- seq_list[c(14:21,36)]
    # 
    print(Sys.time())
    print(names(seq_list))
    # 
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # llply version
    Thickness_list_MBR1_Calibrated <- llply(
        .data = seq_list,
        .fun = function(
            iter_i, 
            Thickness_Calibration, 
            Thickness_list, 
            df_hyper
        ){
            return(
                Thickness_Calibration(
                    Thickness_raw = Thickness_list[[iter_i]],
                    mean_aim = df_hyper$seq_mean_aim[iter_i],
                    sd_aim = df_hyper$seq_sd_aim[iter_i]
                )
            )
        },
        .inform = T,
        .parallel = T,
        # .progress = 'time',
        Thickness_Calibration = Calibrate_Thickness,
        Thickness_list = Thickness_list_MBR1,
        df_hyper = df_MBR1
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Day 43
Thickness_list_MBR1_Calibrated$Day_43_15.04_ <- Calibrate_Thickness(
    Thickness_raw = Thickness_list_MBR1$Day_44_16.04_, 
    mean_aim = 230.41, 
    sd_aim = 14.34
)
Plot_Thickness(
    Thickness = Thickness_list_MBR1_Calibrated$Day_43_15.04_,
    # scale_range = 130,
    scale_range = NA,
    Fig_Title = 'Day_43_15.04_Relaxation_2',
    flag_plot = F, flag_save_plot = T,
    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Default_Calibrated/'
)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# GGplot llply version, MBR1
if (1) {
    tictoc::tic('Plotting, MBR1, Calibrated')
    seq_list <- seq_along(Thickness_list_MBR1_Calibrated)
    names(seq_list) <- names(Thickness_list_MBR1_Calibrated[seq_list])
    # 
    print(Sys.time())
    print(names(seq_list))
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # 
    invisible(
        llply(
            .data = seq_list,
            .fun = function(
                i, 
                Plot_Thickness, 
                Thickness_list, 
                df_hyper
            ){
                Plot_Thickness(
                    Thickness = Thickness_list[[i]], 
                    # scale_range = 130,
                    scale_range = NA,
                    Fig_Title = df_hyper$seq_Fig_Title[i],
                    flag_plot = F,
                    flag_save_plot = T,
                    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Default_Calibrated/'
                )
            },
            .parallel = T,
            # .progress = 'time',
            Plot_Thickness = Plot_Thickness,
            Thickness_list = Thickness_list_MBR1_Calibrated,
            df_hyper = df_MBR1
        )
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Calibrating mean&sd & Plotting, MBR2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# str(Thickness_list_MBR2)
if (0) {
    tictoc::tic('Calibrating, MBR2')
    seq_list <- seq_along(Thickness_list_MBR2)
    names(seq_list) <- names(Thickness_list_MBR2[seq_list])
    # 
    # seq_list <- seq_list[c(14:21,36)]
    # 
    print(Sys.time())
    print(names(seq_list))
    # 
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # llply version
    Thickness_list_MBR2_Calibrated <- llply(
        .data = seq_list,
        .fun = function(
            iter_i, 
            Thickness_Calibration, 
            Thickness_list, 
            df_hyper
        ){
            return(
                Thickness_Calibration(
                    Thickness_raw = Thickness_list[[iter_i]],
                    mean_aim = df_hyper$seq_mean_aim[iter_i],
                    sd_aim = df_hyper$seq_sd_aim[iter_i]
                )
            )
        },
        .inform = T,
        .parallel = T,
        # .progress = 'time',
        Thickness_Calibration = Calibrate_Thickness,
        Thickness_list = Thickness_list_MBR2,
        df_hyper = df_MBR2
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}
# GGplot llply version, MBR2
if (1) {
    tictoc::tic('Plotting, MBR2, Calibrated')
    seq_list <- seq_along(Thickness_list_MBR2_Calibrated)
    names(seq_list) <- names(Thickness_list_MBR2_Calibrated[seq_list])
    # 
    print(Sys.time())
    print(names(seq_list))
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # 
    invisible(
        llply(
            .data = seq_list,
            .fun = function(
                i, 
                Plot_Thickness, 
                Thickness_list, 
                df_hyper
            ){
                Plot_Thickness(
                    Thickness = Thickness_list[[i]], 
                    scale_range = 100,
                    Fig_Title = df_hyper$seq_Fig_Title[i],
                    flag_plot = F,
                    flag_save_plot = T,
                    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_2/Default_Calibrated/'
                )
            },
            .parallel = T,
            # .progress = 'time',
            Plot_Thickness = Plot_Thickness,
            Thickness_list = Thickness_list_MBR2_Calibrated,
            df_hyper = df_MBR2
        )
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}






# Thesis plot -------------------------------------------------------------


if (1) {
    
    # fig_oct_Day_21_24_03_Stable_Flux
    Plot_Thickness(
        Thickness = Thickness_list_MBR1$Day_21_24.03_,
        # scale_range = 130,
        scale_range = c(55, 80),
        Fig_Title = 'Day_21_24_03_Stable_Flux',
        flag_plot = T, flag_save_plot = T,
        save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Thesis_Figure/'
    )
    # fig_oct_Day_27_30_03_Relaxation_1
    Plot_Thickness(
        Thickness = Thickness_list_MBR1$Day_27_30.03_,
        # scale_range = 130,
        scale_range = c(45, 90),
        Fig_Title = 'Day_27_30_03_Relaxation_1',
        flag_plot = T, flag_save_plot = T,
        save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Thesis_Figure/'
    )
    # Day_44_16.04_ , day 43
    Plot_Thickness(
        Thickness = 
            Calibrate_Thickness(Thickness_raw = Thickness_list_MBR1$Day_44_16.04_, 
                                mean_aim = 46.002, sd_aim = 4.408),
        # scale_range = 130,
        scale_range = c(35, 60),
        Fig_Title = 'Day_43_15_04_Relaxation_2',
        flag_plot = T, flag_save_plot = T,
        save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Thesis_Figure/'
    )
    # fig_oct_Day_46_18_04_Air_Scouring
    Plot_Thickness(
        Thickness = Thickness_list_MBR1$Day_46_18.04_,
        scale_range = NULL,
        # scale_range = c(20, 35),
        Fig_Title = 'Day_46_18_04_Air_Scouring',
        flag_plot = T, flag_save_plot = T,
        save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Thesis_Figure/'
    )
    # fig_oct_Day_56_28_04_Relaxation_Air_Scouring
    Plot_Thickness(
        Thickness = Thickness_list_MBR1$Day_56_28.04_,
        # scale_range = 130,
        scale_range = c(5, 30),
        Fig_Title = 'Day_56_28_04_Relaxation_Air_Scouring',
        flag_plot = T, flag_save_plot = T,
        save_folder = '/Users/chengt/Documents/OCT_Scan/Img/MBR_1/Thesis_Figure/'
    )
    
}
