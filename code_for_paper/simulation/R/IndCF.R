## Train and estimate individual causal forest
IndCF <- function(df, covars, honest, is_pred=FALSE, myfit=NULL) {
    if (!is_pred) {
        myfit <- causal_forest(X=as.matrix(df[, ..covars]), Y=df$Y, W=df$Z,
                               honesty=honest, honesty.prune.leaves=F)
        df$tau_hat <- myfit$predictions
        
    } else {
        tmp <- predict(myfit, as.matrix(df[, ..covars]), estimate.variance=FALSE)
        df$tau_hat <- tmp$predictions
        
    }
    
    df_est_res <- df[, c("R", "S", "tau", "tau_hat")]
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
