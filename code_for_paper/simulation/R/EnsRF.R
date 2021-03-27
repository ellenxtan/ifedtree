## Train and estimate ensemble regression tree
EnsRF <- function(df, honest, is_pred=FALSE, myfit=NULL) {
    covars <- c("site", grep("^X", names(df), value=TRUE))
    fml <- as.formula(paste("leaves ~ ", paste(covars, collapse="+")))
    
    if (!is_pred) {
        if (!honest) {
            df[, sitenum := round(mean(leaves)), by=site]  # mean encoding
            covars <- c("sitenum", grep("^X", names(df), value=TRUE))
            myfit <- grf::regression_forest(df[,..covars], df$leaves, 
                                            honesty=FALSE,
                                            tune.parameters="all"
                                            )
            df$tau_hat <- myfit$predictions

        } else if (honest) {
            df[, sitenum := round(mean(leaves)), by=site]  # mean encoding
            covars <- c("sitenum", grep("^X", names(df), value=TRUE))
            myfit <- grf::regression_forest(df[,..covars], df$leaves, 
                                            honesty=TRUE,
                                            tune.parameters="all"
                                            )
            df$tau_hat <- myfit$predictions
            
        }
        df$tau <- df$tau.aug
        
    } else if (is_pred) {
        if (!honest) {
            covars <- c("sitenum", grep("^X", names(df), value=TRUE))
            df$tau_hat <- predict(myfit, df[,..covars])$predictions

        } else if (honest) {
            covars <- c("sitenum", grep("^X", names(df), value=TRUE))
            df$tau_hat <- predict(myfit, df[,..covars])$predictions
            
        }
        
    }
    
    df_est_res <- df[, c("R", "S", "site", "tau", "tau_hat", "role")]
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
