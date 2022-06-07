## Generate augmented data with leaves and SE by causal tree
AugDataCT <- function(df0, df_lst, cfit_lst, K, tau_fn) {
    aug_df <- c()
    for (k in 1:K) {
        df0[, leaves := predict(cfit_lst[[k]], newdata=df0)]
        num_leaves <- length(levels(factor(df0$leaves)))
        df0$leaff <- factor(df0$leaves, labels = seq(num_leaves))
        
        if (num_leaves > 1) {
            lm_ct <- lm(as.formula("Y ~ -1 + leaff + Z*leaff - Z"), data = df0)
        } else {
            lm_ct <- lm(as.formula("Y ~ Z"), data = df0)
        }
        
        lm_se <- data.table(broom::tidy(lm_ct, conf.int=TRUE)[(num_leaves+1):(2*num_leaves), c("std.error")])
        lm_se$leaff <- rownames(lm_se)
        names(lm_se) <- c("leavesse", "leaff")
        setkey(df0, leaff)
        setkey(lm_se, leaff)
        merge_df0 <- merge(df0, lm_se, all.x=TRUE)
        merge_df0 <- merge_df0[order(merge_df0$idx), ]  # reorder by idx within site
        setkey(df0, NULL)
        merge_df0$leaff <- NULL
        
        aug_df <- data.table(rbind(aug_df, 
                                   cbind(merge_df0, site=k, site_real=unique(df_lst[[k]]$U))))
    }
    
    aug_df$U <- aug_df$site_real
    if (!is.null(tau_fn)) {
        aug_df[, tau.aug := eval(tau_fn)] # tau_fn contains U
    }
    aug_df$site <- factor(aug_df$site)
    aug_df$tau <- NULL # remove tau as for site 1 only
    
    return(aug_df)
}
