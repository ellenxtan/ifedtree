## Train and estimate individual X-learner with BART
IndXBART <- function(df, covars, honest, is_pred=FALSE, myfit=NULL) {
    if (!is_pred) {
        x_train <- as.data.frame(df[, ..covars])
        y_train <- df$Y
        w_train <- df$Z
        myfit <- X_BART(feat = x_train, tr = w_train, yobs = y_train)
        sink("NUL")
        df$tau_hat <- EstimateCate(myfit, x_train)
        sink()
        
    } else {
        x_test <- as.data.frame(df[, ..covars])
        sink("NUL")
        df$tau_hat <- EstimateCate(myfit, x_test)
        sink()

    }
    
    df_est_res <- df[, c("R", "S", "tau", "tau_hat")]
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
