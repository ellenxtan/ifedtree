# exponential weighted model averaging
EWMAwrap <- function(loc.mod, df_test, df_lst, cfit_lst, coord, loc.fit, K, is_oracle) {
    
    df_est <- df_lst[[coord]][role=="est", ]
    covars <- grep("^X", names(df_est), value=TRUE)  # include X's
    
    if (is_oracle) {
        Y <- df_est$tau
    } else {
        Y <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                    df_use=df_est, is_pred=NULL, 
                    covars=covars, honest1=NULL, myfit=loc.fit)
    }
    
    ## obtain wts
    mse_errs <- c()
    for (k in 1:K) {
        tau_est <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                          df_use=df_est, is_pred=NULL, 
                          covars=covars, honest1=NULL, myfit=cfit_lst[[k]])
        
        mse_errs <- c(mse_errs, -sum((Y - tau_est)^2))
    }
    wts <- exp(mse_errs - matrixStats::logSumExp(mse_errs)) # log-sum-exp trick
    
    ## test set for final tau
    if (!is.null(df_test)) {
        tau_wts <- c()
        for (k in 1:K) {
            tau_test <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                              df_use=df_test, is_pred=NULL, 
                              covars=covars, honest1=NULL, myfit=cfit_lst[[k]])

            tau_wts <- cbind(tau_wts, tau_test * wts[k])
        }
        
        df_test$tau_hat <- rowSums(tau_wts)
        
        df_est_res <- df_test[, c("S", "tau", "tau_hat")] #"R", 
        
        EvaluateMetrics(df_est_res, paste0("ewma-",is_oracle,":"))
        
    } else {
        df_est_res <- NULL
        
    }
    
    return(list(df_est_res=df_est_res, wts=wts))
}
