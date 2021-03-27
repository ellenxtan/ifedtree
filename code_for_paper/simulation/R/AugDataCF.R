## Generate augmented data with leaves and SE by causal forest
AugDataCF <- function(df0, cfit_lst, tau_fn, K, coord, is_train=TRUE) {
    covars = grep("^X", names(df0), value=TRUE)
    aug_df <- c()
    for (k in 1:K) {
        tmpdf = copy(df0)
        tmpdf[, tau.aug := eval(tau_fn)] # tau_fn contains k
        
        if (is_train & k==coord) {  # coordinating site use OOB pred & SE
            tmpdf[, leaves := cfit_lst[[k]]$predictions]
            pred_cf <- predict(cfit_lst[[k]], estimate.variance=TRUE)
            tmpdf[, leavesse := sqrt(pred_cf$variance.estimates)]
        } else {
            pred_cf <- predict(cfit_lst[[k]], as.matrix(tmpdf[, ..covars]), estimate.variance=TRUE)
            tmpdf[, leaves := pred_cf$predictions]
            tmpdf[, leavesse := sqrt(pred_cf$variance.estimates)]
        }
        
        aug_df <- data.table(rbind(aug_df, cbind(tmpdf, site = k)))
    }
    
    aug_df$site <- factor(aug_df$site)
    return(aug_df)
}
