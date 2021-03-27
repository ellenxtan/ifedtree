## Generate augmented data with leaves and SE by X-learner with RF/BART
AugDataXl <- function(df0, cfit_lst, tau_fn, K) {
    covars = grep("^X", names(df0), value=TRUE)
    aug_df <- c()
    for (k in 1:K) {
        tmpdf = copy(df0)
        tmpdf[, tau.aug := eval(tau_fn)] # tau_fn contains k
        
        x_test <- as.data.frame(tmpdf[, ..covars])
        
        # SE in meta-analysis
        sink("NUL")
        xl_ci <- CateCI(cfit_lst[[k]], x_test, B = 50)
        sink()
        
        tmpdf$leaves <- xl_ci$pred
        tmpdf$leavesse <- with(xl_ci, (pred - X5.)/1.96)
        
        aug_df <- data.table(rbind(aug_df, cbind(tmpdf, site = k)))
    }
    
    aug_df$site <- factor(aug_df$site)
    return(aug_df)
}
