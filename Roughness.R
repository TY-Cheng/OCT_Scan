
# Roughness Quantification ------------------------------------------------
if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    rm(list = ls())
    setwd('~/Documents/OCT_Scan/')
    library(plyr)
    library(tidyverse)
    load('OCT_Roughness.RData')
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # mat_thickness <- Thickness_list_MBR1_Calibrated$Day_28_31.03_
    # mat_thickness <- Thickness_list_MBR1_Calibrated$Day_29_01.04_
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    Get_MaternCovEstimate <- function(mat_thickness, fraction_sample = 0.05) {
        require(sp)
        require(gstat)
        set.seed(2020)
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        grid_thickness <- expand.grid(
            X=seq_len(NCOL(mat_thickness)),
            Y=seq_len(NCOL(mat_thickness))
        )
        grid_thickness$Z <- c(t(mat_thickness))
        coordinates(grid_thickness) <- ~X+Y
        grid_thickness <- grid_thickness[
            sample(NROW(grid_thickness), round(NROW(grid_thickness) * fraction_sample)), 
        ]
        vgm_sample <- variogram(Z~1, grid_thickness)
        # # # # # #
        # obj_fun_SSE = function(x) {
        #     # return the SSE, set vgm_fit on the outside
        #     # # # # # # # # # # # # # 
        #     # # # # # # # # # # # # # 
        #     # # dont try parallel # # 
        #     # # # # # # # # # # # # # 
        #     # # # # # # # # # # # # # 
        #     vgm_fit <<- fit.variogram(
        #         object = vgm_sample, model = vgm("Mat", nugget = NA, kappa = x)
        #     )
        #     # attr(vgm_fit <<- fit.variogram(vgm1, vgm(,"Mat",kappa=x)),"SSErr")
        #     return(attr(vgm_fit, 'SSErr'))
        # }
        # optimize(obj_fun_SSE, c(.1, 5))
        # # # # # #
        vgm_fit <- fit.variogram(object = vgm_sample, 
                                 model = vgm('Mat', nugget = NA, kappa = NA), 
                                 fit.kappa = T)
        # plot(vgm_sample, vgm_fit)
        vgm_estimate <- as.data.frame(vgm_fit)
        vgm_estimate <- vgm_estimate[vgm_estimate$model=='Mat', 
                                     c('model', 'psill', 'range', 'kappa')]
        rownames(vgm_estimate) <- NULL
        return(vgm_estimate)
    }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}


if (1) {
    # use the real ones!
    # not the calibrated ones!
    # temp_seq <- names(Thickness_list_MBR1_Calibrated)
    temp_seq <- names(Thickness_list_MBR1)
    names(temp_seq) <- temp_seq
    
    list_MaternCovEstimate <- llply(
        .data = temp_seq, .progress = 'time', .parallel = F,
        .fun = function(iter_name) {
            iter_matcovestimate <- Get_MaternCovEstimate(
                # mat_thickness = Thickness_list_MBR1_Calibrated[[iter_name]], 
                mat_thickness = Thickness_list_MBR1[[iter_name]], 
                fraction_sample = .05
            )
            iter_matcovestimate$seq_Fig_Index <- iter_name
            return(iter_matcovestimate)
        }
    )
    df_MaternCovEstimate <- Reduce(f = 'rbind', list_MaternCovEstimate)
    df_MaternCovEstimate <- arrange(df_MaternCovEstimate, seq_Fig_Index)
    # 
    df_MaternCovEstimate$seq_Fig_Title <- 
        paste0(
            df_MaternCovEstimate$seq_Fig_Index,
            c(rep('Stable_Flux', 5),
              rep('Relaxation_1', 9),
              rep('Relaxation_2', 4), # notice we add one Day 43
              rep('Air_Scouring', 3),
              rep('Relaxation_Air_Scouring', 3),
              rep('Autopsy', 15))
        )
    # df_MaternCovEstimate$seq_Fig_Title <- 
    #     paste0(
    #         df_MaternCovEstimate$seq_Fig_Index,
    #         c(rep('Stable_Flux', 5), 
    #           rep('Relaxation_1', 9), 
    #           rep('Relaxation_2', 5), # notice we add one Day 43
    #           rep('Air_Scouring', 3),
    #           rep('Relaxation_Air_Scouring', 3),
    #           rep('Autopsy', 15))
    #     )
    
    df_MaternCovEstimate_MBR1 <- df_MaternCovEstimate
    rownames(df_MaternCovEstimate_MBR1) <- df_MaternCovEstimate_MBR1$seq_Fig_Title
    # df_MaternCovEstimate_MBR1_Calibrated <- df_MaternCovEstimate
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    rm(temp_seq, vgm_fit, df_MaternCovEstimate)
}



if (0) {
    temp_seq <- names(Thickness_list_MBR1)
    names(temp_seq) <- temp_seq
    temp_df <- llply(
        .data = temp_seq, .progress = 'time', .parallel = F, 
        .fun = function(iter_name) {
            iter_seq <- c(Thickness_list_MBR1[[iter_name]])
            iter_result <- psych::describe(iter_seq, quant = c(.25, .75)) %>% 
                data.frame()
            iter_result$seq_Fig_Index <- iter_name
            iter_result$roughness_absolute <- mean(abs(iter_seq - mean(iter_seq)))
            iter_result$roughness_relative <- mean(abs(iter_seq - mean(iter_seq)) / 
                                                       mean(iter_seq))
            return(iter_result)
        })
    temp_df <- Reduce('rbind', temp_df)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    temp_df$seq_Fig_Title <- 
        paste0(
            temp_df$seq_Fig_Index,
            c(rep('Stable_Flux', 5),
              rep('Relaxation_1', 9),
              rep('Relaxation_2', 4), # notice we add one Day 43
              rep('Air_Scouring', 3),
              rep('Relaxation_Air_Scouring', 3),
              rep('Autopsy', 15))
        )
    rownames(temp_df) <- temp_df$seq_Fig_Title
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    cbind(df_MaternCovEstimate_MBR1, temp_df)
    df_result_roughness <- clipr::read_clip_tbl()
    
    
}


