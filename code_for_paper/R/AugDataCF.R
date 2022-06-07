## Generate augmented data with leaves and SE by causal forest
AugDataCF <- function(df0, df_lst, cfit_lst, K, tau_fn) {
    covars = grep("^X", names(df0), value=TRUE)
    aug_df <- c()
    for (k in 1:K) {
        pred_cf <- predict(cfit_lst[[k]], as.matrix(df0[, ..covars]), estimate.variance=TRUE)
        df0[, leaves := pred_cf$predictions]
        df0[, leavesse := sqrt(pred_cf$variance.estimates)]
        
        aug_df <- data.table(rbind(aug_df, cbind(df0, site=k, site_real=unique(df_lst[[k]]$U))))
    }
    
    aug_df$U <- aug_df$site_real
    if (!is.null(tau_fn)) {
        aug_df[, tau.aug := eval(tau_fn)] # tau_fn contains U
    }
    aug_df$site <- factor(aug_df$site)
    aug_df$tau <- NULL # remove tau as for site 1 only
    
    return(aug_df)
}
