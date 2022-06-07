LinearStack <- function(loc.mod, df_test, df_lst, cfit_lst, coord, loc.fit, is_oracle) {
    
    df_train <- df_lst[[coord]][role=="est", ]
    covars <- grep("^X", names(df_train), value=TRUE)  # include X's
    
    if (is_oracle) {
        Y <- df_train$tau
    } else {
        Y <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                    df_use=df_train, is_pred=NULL, 
                    covars=covars, honest1=NULL, myfit=loc.fit)
    }
    
    df_stack_tr <- c()
    df_stack_te <- c()
    for (k in 1:K) {
        tau_train <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                          df_use=df_train, is_pred=NULL, 
                          covars=covars, honest1=NULL, myfit=cfit_lst[[k]])
        
        df_stack_tr <- cbind(df_stack_tr, tau_train)
        
        tau_test <- IndFit(loc.mod=loc.mod, pred_only=TRUE, 
                            df_use=df_test, is_pred=NULL, 
                            covars=covars, honest1=NULL, myfit=cfit_lst[[k]])
        
        df_stack_te <- cbind(df_stack_te, tau_test)
    }
    
    df_stack_tr <- as.data.frame(df_stack_tr)
    names(df_stack_tr) <- paste0("yhat",seq(K))
    df_stack_tr$Y <- Y
    
    lr.fit <- lm(Y~., data=df_stack_tr)
    
    df_stack_te <- as.data.frame(df_stack_te)
    names(df_stack_te) <- paste0("yhat",seq(K))
    
    df_test$tau_hat <- predict(lr.fit, df_stack_te)
    
    df_est_res <- df_test[, c("S", "tau", "tau_hat")] #"R", 
    
    EvaluateMetrics(df_est_res, paste0("stack-",is_oracle,":"))
    
    return(list(df_est_res=df_est_res))
}
