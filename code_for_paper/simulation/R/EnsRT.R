## Train and estimate ensemble regression tree
EnsRT <- function(aug_df, honest=F, is_pred=FALSE, myfit=NULL, leaves_tab=NULL) {
    df <- copy(aug_df)
    
    covars <- grep("^X", names(df), value=TRUE)
    fml <- as.formula(paste("leaves ~ site + ", paste(covars, collapse="+")))
    
    if (!is_pred) {
        if (!honest) {
            myfit <- rpart(fml, data=df, xval=5)
            myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
            df$tau_hat <- predict(myfit)
            
        } else if (honest) {
            role <- rep("train", nrow(df))
            role[sample(1:nrow(df), nrow(df) / 2, replace = FALSE)] <- "est"
            df$role <- role
            
            ## train - build tree
            myfit <- rpart(fml, data=df[role=="train", ], xval=5)
            myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
            
            ## honest estimate by df_est
            ## tree for predict which leave the obs belongs to
            # replace predicted y values in the model frame with the corresponding node numbers
            et_node = myfit
            et_node$frame$yval = as.numeric(rownames(et_node$frame))
            df[role=="est", leaff := predict(et_node, newdata=df[role=="est", ], type='vector')]
            df[role=="est", tau_hat := mean(leaves, na.rm=TRUE), by = leaff]
            
            leaves_tab <- unique(df[role=="est", c("leaff", "tau_hat")])
            df[role=="train", leaff := predict(et_node, newdata=df[role=="train", ], type='vector')]
            # assign honest estimates to train set
            for (i in 1:length(leaves_tab$leaff)) {
                df[role=="train" & leaff==leaves_tab$leaff[i], tau_hat := leaves_tab$tau_hat[i]]
            }
            rm(et_node)
        }
        df$tau <- df$tau.aug
        
    } else if (is_pred) {
        if (!honest) {
            df$tau_hat <- predict(myfit, df)
            
        } else if (honest) {
            et_node = myfit
            et_node$frame$yval = as.numeric(rownames(et_node$frame))
            df[, leaff := predict(et_node, newdata=df, type='vector')]  # predict which leaf falls into
            for (i in 1:length(leaves_tab$leaff)) {  # assign honest estimates to test set
                df[leaff==leaves_tab$leaff[i], tau_hat := leaves_tab$tau_hat[i]]
            }
        }
    }
    
    df_est_res <- df[, c("R", "S", "site", "tau", "tau_hat")]
    
    return(list(myfit=myfit, df_est_res=df_est_res, leaves_tab=leaves_tab))
}
