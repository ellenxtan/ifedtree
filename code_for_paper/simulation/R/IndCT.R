## Train and estimate individual causal tree
IndCT <- function(df, covars, honest, is_pred=FALSE, myfit=NULL) {
    if (!is_pred) {
        fml <- as.formula(paste("Y~", paste(covars, collapse="+")))
        
        ## train by df_tr(build) & df_est(est) ----
        sink("NUL")
        if (!honest) {
            myfit <- causalTree(formula = fml, data = df, treatment = df$Z,
                                split.Rule="CT", split.Honest=F, split.alpha=1,
                                cv.option="CT", cv.Honest=F, cp=0, xval=5)
        } else if (honest) {
            df_tr <- df[role=="train", ]
            df_est <- df[role=="est", ]
            myfit <- honest.causalTree(
                formula = fml,
                data = df_tr, treatment = df_tr$Z,
                est_data = df_est, est_treatment = df_est$Z,
                split.Rule="CT", split.Honest=T,
                cv.option="CT", cv.Honest=T,
                cp=0, xval=5, HonestSampleSize=nrow(df_est))
        }
        sink()
        
        myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])

        ## evaluate by df ----
        df$tau_hat <- predict(myfit, newdata=df)
        
    } else {
        df$tau_hat <- predict(myfit, newdata=df)
    }
    
    df_est_res <- df[, c("R", "S", "tau", "tau_hat")]
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
