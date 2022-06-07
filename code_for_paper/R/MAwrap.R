# basic model averaging (simple MA with weights inverse to sample size, assumes no site heterogeneity)
MAwrap <- function(loc.mod, test_df0, df_lst, cfit_lst, K) {
    
    covars <- grep("^X", names(test_df0), value=TRUE)  # include X's
    
    ## obtain wts
    wts <- c()
    for (k in 1:K) {
        wt <- 1 / nrow(df_lst[[k]])
        wts <- c(wts, wt)
    }
    wts <- wts / sum(wts)
    
    ## test set for final tau
    tau_wts <- c()
    for (k in 1:K) {
        tau_test <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                           df_use=test_df0, is_pred=NULL, 
                           covars=covars, honest1=NULL, myfit=cfit_lst[[k]])
        
        tau_wts <- cbind(tau_wts, tau_test * wts[k])
    }
    
    test_df0$tau_hat <- rowSums(tau_wts)
    
    df_est_res <- test_df0[, c("S", "tau", "tau_hat")]
    
    EvaluateMetrics(df_est_res, "ma:")
    
    return(list(df_est_res=df_est_res))
}
