## Train and estimate ensemble regression tree
EnsRF <- function(aug_df, is_oracle, test_df, honest, pkg, coord) {
    
    ## use mean encoding for site index (with grf) - all variables continuous
    if (pkg == "grf") {
        if (class(aug_df$site_real) == "factor") {
            aug_df$site_real <- as.numeric(as.character(aug_df$site_real))
            test_df$site_real <- as.numeric(as.character(test_df$site_real))
        }
        
        if (is_oracle) { # use tau.aug
            aug_df[, sitenum := round(mean(tau.aug)), by=site]  # target encoding
            test_df$sitenum = aug_df[site==coord, ]$sitenum[1]  # coord
            covars <- c("sitenum", grep("^X", names(aug_df), value=TRUE))
            
            myfit <- grf::regression_forest(X=aug_df[,..covars], Y=aug_df$tau.aug, 
                                            honesty=honest) #, tune.parameters="all"
        } else { # use site & leaves
            aug_df[, sitenum := round(mean(leaves)), by=site]  # target encoding
            test_df$sitenum = aug_df[site==coord, ]$sitenum[1]  # coord
            covars <- c("sitenum", grep("^X", names(aug_df), value=TRUE))
            
            myfit <- grf::regression_forest(X=aug_df[,..covars], Y=aug_df$leaves, 
                                            honesty=honest) #, tune.parameters="all"
        }
        
        test_df$tau_hat <- predict(myfit, test_df[,..covars])$predictions
    }
    
    df_est_res <- test_df[, c("S", "site", "tau", "tau_hat")]
    
    EvaluateMetrics(df_est_res, paste0("ef-",is_oracle,"-",pkg,":"))
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
