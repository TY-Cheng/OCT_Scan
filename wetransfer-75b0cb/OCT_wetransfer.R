if (0) {
    rm(list = ls())
    source('/Users/chengt/Documents/OCT_Scan/OCT_function_utils.R')
    setwd(dir = '~/Documents/OCT_Scan/wetransfer-75b0cb/')
    load('OCT_wetransfer.RData')
}

# Import MAT files --------------------------------------------------------
if (0) {
    # Rootfolder <- '/Volumes/Seagate_Backup/OCT_Scan_PreProcessing/Flowcell_3D/MAT_files_/'
    Rootfolder <- '~/Documents/OCT_Scan/wetransfer-75b0cb/'
    setwd(Rootfolder)
    # seq_MATfolder <- list.files() # sorted as char instead, not as num
    seq_MATfolder <- c('Time_8', 'Time_12', 'Time_21')
    names(seq_MATfolder) <- seq_MATfolder
    Img_list_Time_81221 <- llply(
        .data = seq_MATfolder, 
        .fun = Get_Thickness_3D_From_MAT,
        .progress = progress_time(),
        Rootfolder = Rootfolder,
        n_row = 800,
        n_col = 2202
    )
}


# Hyperparameters ---------------------------------------------------------
df_Time_81221 <- data.frame(
    seq_Fig_Index = c(paste0('Time_', c(8, 12, 21))),
    seq_Fig_Title = c(
        paste0('Time_', 
               stringr::str_pad(
                   c(8, 12, 21),
                   width = 2, pad = '0')
        )
    ),
    seq_Denoise_Type = 2,
    seq_Loess_Span = rep(0.02, 3),
    seq_Quantile_Min = rep(.01, 3), 
    seq_Quantile_Max = rep(.999, 3),
    seq_mean_aim = NA,
    seq_sd_aim = NA, 
    stringsAsFactors = F
)


# Processing and Imaging --------------------------------------------------
# from (not same length) to 666

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Interpolation & Filtering & Plotting , Time_81221
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# str(Img_list_Time_81221)
if (1) {
    tictoc::tic('Interpolation & Filtering, Time_81221')
    seq_list <- as.list(as.data.frame(t(df_Time_81221), stringsAsFactors = F))
    names(seq_list) <- df_Time_81221$seq_Fig_Index
    # 
    # seq_list <- seq_list[c(14:21,36)]
    # 
    print(Sys.time())
    print(names(seq_list))
    # 
    para_socket_cl <- makeCluster(parallel::detectCores())
    registerDoParallel(para_socket_cl)
    # llply version
    Thickness_list_Time_81221 <- llply(
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
                    Remove_Consecutive_NA = NA_Action,
                    n_row = 800
                )
            )
        },
        .inform = T,
        .parallel = T,
        # .progress = 'time',
        Img_list = Img_list_Time_81221,
        Img_Processing = Impute_Filter_Image_from_Matrix,
        NA_Action = Remove_Consecutive_NA
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}

# GGplot llply version, Time_81221
if (1) {
    tictoc::tic('Plotting, Time_81221')
    seq_list <- seq_along(Thickness_list_Time_81221)
    names(seq_list) <- names(Thickness_list_Time_81221[seq_list])
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
                    x_in_mm = 8,
                    y_in_mm = 8,
                    flag_save_plot = T,
                    save_folder = '/Users/chengt/Documents/OCT_Scan/Img/Time_81221/Default/'
                )
            },
            .parallel = T,
            # .progress = 'time',
            Plot_Thickness = Plot_Thickness,
            Thickness_list = Thickness_list_Time_81221,
            df_hyper = df_Time_81221
        )
    )
    # 
    rm(seq_list)
    stopCluster(para_socket_cl)
    print(Sys.time())
    tictoc::toc()
}


# Thickness & Roughness ---------------------------------------------------
# In um
df_describe <- llply(
    .data = names(Img_list_Time_81221), 
    .fun = function(iter_name) {
        iter_df <- psych::describe(c(Img_list_Time_81221[[iter_name]]) * 
                                       # # # # # # # # # # # # # # # #
                                       2.1 # # # # # # # # # # # # # #
                                   # # # # # # # # # # # # # # # #
        )
        rownames(iter_df) <- iter_name
        return(iter_df)
    }
)
df_describe <- Reduce(rbind, df_describe)

clipr::write_clip(df_describe)
