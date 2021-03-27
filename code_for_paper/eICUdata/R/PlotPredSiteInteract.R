## Prediction plot for site and variable of interest
PlotPredSiteInteract <- function(df, rf, covars_all, var1, var2, set="real", siteRank=NULL, grids=100, myseed=1234) {
    # cat("var1 =", var1, "; var2 =", var2, "\n")
    if (class(df[[var1]]) == "factor") {
        grid1 <- length(levels(df[[var1]]))
        df[[var1]] <- as.numeric(as.character(df[[var1]]))
    } else {
        grid1 <- grids
    }
    if (class(df[[var2]]) == "factor") {
        grid2 <- length(levels(df[[var2]]))
        df[[var2]] <- as.numeric(as.character(df[[var2]]))
    } else {
        grid2 <- grids
    }
    # cat("grid1:", grid1, "; grid2:", grid2, "\n")
    newdata <- expand.grid(seq(min(df[[var1]]), max(df[[var1]]), length.out = grid1),
                           seq(min(df[[var2]]), max(df[[var2]]), length.out = grid2))
    # cat("values of", var2, ":", unique(round(newdata$Var2)), "\n")
    colnames(newdata) <- c(var1, var2)
    
    other_vars <- setdiff(covars_all, c(var1, var2))
    
    # n <- nrow(df)
    # j = 1
    # for (usevar in other_vars) {
    #     set.seed(myseed*j)
    #     newdata[[usevar]] <- df[[usevar]][sample(1:n, nrow(newdata), replace = TRUE)]
    #     j <- j+1
    # }
    
    if (set=="simulation") {
        # mean: X1  0.01462419,  X2 -0.006119085, X3 -0.01196118, X4 -0.01034596,  X5  0.01916636,
        #       X6 -0.003380824, X7 -0.01038755,  X8 -0.01302276, X9 -0.005269341, X10 0.008012065
        for (usevar in other_vars) {
            newdata[[usevar]] <- NA_real_
            # newdata[[usevar]] <- dplyr::case_when(usevar=="X1" ~ 0.01462419,
            #                                       usevar=="X2" ~ -0.006119085,
            #                                       usevar=="X3" ~ -0.01196118,
            #                                       usevar=="X4" ~ -0.01034596,
            #                                       usevar=="X5" ~ 0.01916636,
            #                                       usevar=="X6" ~ -0.003380824,
            #                                       usevar=="X7" ~ -0.01038755,
            #                                       usevar=="X8" ~ -0.01302276,
            #                                       usevar=="X9" ~ -0.005269341,
            #                                       usevar=="X10" ~ 0.008012065
            #                                       )
            newdata[[usevar]] <- dplyr::case_when(usevar=="X1" ~ 0,
                                                  usevar=="X2" ~ 0,
                                                  usevar=="X3" ~ 0,
                                                  usevar=="X4" ~ 0,
                                                  usevar=="X5" ~ 0,
                                                  usevar=="X6" ~ 0,
                                                  usevar=="X7" ~ 0,
                                                  usevar=="X8" ~ 0,
                                                  usevar=="X9" ~ 0,
                                                  usevar=="X10" ~ 0,
                                                  usevar=="site" ~ 1,  # X1-X2 interaction within odd sites
            )
        }
    } else if (set=="real") {
        # median: age 65, bmi 28.34, sofatotal 6, duration 126.55, female
        for (usevar in other_vars) {
            newdata[[usevar]] <- NA_real_
            newdata[[usevar]] <- dplyr::case_when(usevar=="age" ~ 65,
                                           usevar=="bmi" ~ 28.34,
                                           usevar=="sofa" ~ 6,
                                           usevar=="duration" ~ 126.55,
                                           usevar=="female" ~ 1
                                           )
        }
    }
    
    
    newdata$site <- factor(newdata$site)
    if (set=="real") {
        newdata$female <- factor(newdata$female)
    }
    newdata$prediction <- predict(rf, newdata[covars_all])$predictions
    
    if (set=="simulation") {  # first plot odd sites, then even sites
        newdata$site <- factor(newdata$site, levels=c(seq(1,20,2), seq(2,20,2)))
    } else if (set=="real") {  # real data: rank by avg trt effect
        newdata$site <- factor(newdata$site, levels=siteRank, labels=as.character(seq(1,K)))
    }
    
    # hist(newdata$prediction)
    cat("range of prediction", range(newdata$prediction), "\n")
    
    if (set=="simulation") {
        if (class(newdata[[var1]]) == "factor") {  # site vs. X1/X2
            fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) + #, fill = "prediction"
                geom_tile(aes(fill=prediction)) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0))
        } else {  # X1 vs. X2
            fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) + #, fill = "prediction"
                geom_tile(aes(fill=prediction)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0))
        }
        fig_pred <- fig_pred +
            scale_fill_gradient2(midpoint = 0, 
                                 limits=c(-5.1, 2.5), # ensure all subfigs have same legend
                                 # limits=c(min(newdata$prediction), max(newdata$prediction)),
                                 low = "blue", mid = 'white', high = "red", #, "#B2182B"
                                 breaks=c(-4,-2,0,2), 
                                 labels=c(-4,-2,0,2)) +
            theme_classic() +
            theme(legend.position="top",
                  legend.title=element_text(size=13,face="bold"),
                  legend.text=element_text(size=11,face="bold"),
                  axis.text=element_text(size=12,face="bold"),
                  axis.title=element_text(size=14,face="bold")) +
            guides(fill = guide_colorbar(barwidth=12, barlength=2))  # length of legend colorbar
        
        
    } else if (set=="real") {
        # ## option 1: discretize prediction into intervals
        # ## display.brewer.pal(n = 8, name = 'RdBu') # "RdPu"
        # ## brewer.pal(n = 8, name = "RdBu")
        # # newdata$pred <- cut(newdata$prediction,breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.05,0.1))
        # newdata$pred <- cut(newdata$prediction,breaks = c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05))
        # 
        # ## option: fix other covariates at their median
        # if (class(newdata[[var2]]) == "factor") {
        #     fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) +
        #         geom_tile(aes(fill=pred), height=1.01, width=1.01) +  #remove white lines
        #         scale_x_discrete(expand = c(0, 0)) +
        #         scale_y_discrete(expand = c(0, 0))
        # } else {
        #     fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) +
        #         geom_tile(aes(fill=pred)) +
        #         scale_x_discrete(expand = c(0, 0)) +
        #         scale_y_continuous(expand = c(0, 0))
        # }
        # fig_pred <- fig_pred +
        #     # scale_fill_brewer(palette = "RdBu", direction=-1) +
        #     # scale_fill_manual(breaks=c("(-0.4,-0.3]", "(-0.3,-0.2]", "(-0.2,-0.1]",
        #     #                            "(-0.1,0]", "(0,0.05]", "(0.05,0.1]"),
        #     #                   values=c("#2166AC", "#4393C3", "#92C5DE",
        #     #                            "#D1E5F0", "#F4A582", "#D6604D")) +
        #     scale_fill_manual(breaks=c("(-0.25,-0.2]","(-0.2,-0.15]",
        #                                "(-0.15,-0.1]", "(-0.1,-0.05]", "(-0.05,0]", "(0,0.05]"),
        #                       values=c("#053061", "#2166AC", "#4393C3", "#92C5DE",
        #                                "#D1E5F0", "#F4A582"), #"#D6604D"
        #                       labels=c("(-0.25,-0.2]","(-0.2,-0.15]",
        #                                "(-0.15,-0.1]", "(-0.1,-0.05]", "(-0.05,0]", "(0,0.05]")
        #                       ) + 
        #     theme_classic() +
        #     theme(legend.position="top",
        #           legend.title=element_text(size=13,face="bold"),
        #           legend.text=element_text(size=11,face="bold"),
        #           axis.text=element_text(size=12,face="bold"),
        #           axis.title=element_text(size=14,face="bold")) +
        #     guides(fill = guide_legend(nrow=2,byrow=TRUE))
        
        ## option 2: raw prediction (continuous)
        if (class(newdata[[var2]]) == "factor") {
            fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) + #, fill = "prediction"
                geom_tile(aes(fill=prediction), height=1.01, width=1.01) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0))
        } else {
            fig_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) + #, fill = "prediction"
                geom_tile(aes(fill=prediction)) +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0)) 
        }
        fig_pred <- fig_pred +
            scale_fill_gradient2(midpoint = 0, 
                                 limits=c(-0.23, 0.12), # ensure all subfigs have same legend
                                 # limits=c(min(newdata$prediction), max(newdata$prediction)),
                                 low = "blue", mid = 'white', high = "red", #"#B2182B", ##2166AC
                                 # breaks = c(-0.3,-0.2,-0.1,0.0),#,max(newdata$pred2)),
                                 # labels = c("-0.3","-0.2","-0.1","0.0")#,"0.05")
                                 ) +
            # scale_fill_gradientn(colours=c("black","darkblue",brightness(c("blue",#"lightblue",
            #                                "white", "red"),1)), #, "darkred"
            #                      values=rescale(c(-0.3,-0.2,-0.1,
            #                                       0.0,0.1)), #,0.1
            #                      # breaks = c(-0.3,-0.2,-0.1,0.0,max(newdata$prediction)-0.01),
            #                      # labels = c("-0.3","-0.2","-0.1","0.0","0.1"),
            #                      guide="colorbar", space="Lab") +
            scale_colour_hue(c=150, l=80) +
            theme_classic() +
            theme(legend.position="top",
                  legend.title=element_text(size=13,face="bold"),
                  legend.text=element_text(size=11,face="bold"),
                  axis.text=element_text(size=12,face="bold"),
                  axis.title=element_text(size=14,face="bold")) +
            guides(fill = guide_colorbar(barwidth=12, barlength=2))  # length of legend colorbar
        
    }
    
    fig_pred[["labels"]][["fill"]] <- "Causal estimates"
    
    # print(fig_pred)
    # cat("-------------------------------------------------------------------\n")
    return(fig_pred)
}
