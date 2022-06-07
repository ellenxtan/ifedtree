## Train and estimate ensemble regression tree
EnsRT <- function(aug_df, is_oracle, test_df, honest, is_encode) {
    
    covars <- grep("^X", names(aug_df), value=TRUE)
    if (is_encode) {# target encoding
        if (is_oracle) {
            aug_df[, sitenum := round(mean(tau.aug)), by=site]  # target encoding
            fml <- as.formula(paste("tau.aug ~ sitenum + ", paste(covars, collapse="+")))
        } else {
            aug_df[, sitenum := round(mean(leaves)), by=site]  # target encoding
            fml <- as.formula(paste("leaves ~ sitenum + ", paste(covars, collapse="+")))
        }
        test_df$sitenum = aug_df[site==coord, ]$sitenum[1]  # coord
        aug_df$sitenum <- factor(aug_df$sitenum)
        test_df$sitenum <- factor(test_df$sitenum)
        
    } else {# no encoding
        if (is_oracle) { # use tau.aug
            fml <- as.formula(paste("tau.aug ~ site + ", paste(covars, collapse="+")))
        } else { # use site & leaves
            fml <- as.formula(paste("leaves ~ site + ", paste(covars, collapse="+")))
        }
    }
    
    
    if (honest) {
        role <- rep("train", nrow(aug_df))
        role[sample(1:nrow(aug_df), nrow(aug_df) / 2, replace = FALSE)] <- "est"
        aug_df$role <- role
        
        ## train - build tree
        myfit <- rpart(fml, data=aug_df[role=="train", ], xval=5)
        myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
        
        ## est - get honest estimate
        ## tree for predict which leave the obs belongs to
        # replace predicted y values in the model frame with the corresponding node numbers
        et_node = myfit
        et_node$frame$yval = as.numeric(rownames(et_node$frame))
        aug_df[role=="est", leaff := predict(et_node, newdata=aug_df[role=="est", ], type='vector')]
        aug_df[role=="est", tau_hat := mean(leaves, na.rm=TRUE), by = leaff]
        leaves_tab <- unique(aug_df[role=="est", c("leaff", "tau_hat")])
        
        ## test - assign honest estimate generated from est
        test_df[, leaff := predict(et_node, newdata=test_df, type='vector')]  # predict which leaf falls into
        for (i in 1:length(leaves_tab$leaff)) {  # assign honest estimates to test set
            test_df[leaff==leaves_tab$leaff[i], tau_hat := leaves_tab$tau_hat[i]]
        }
        
    } else { # not honest
        myfit <- rpart(fml, data=aug_df, xval=5, cp=0)
        myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
        
        test_df$tau_hat <- predict(myfit, test_df)
    }
    
    
    df_est_res <- test_df[, c("S", "site", "tau", "tau_hat")]
    
    EvaluateMetrics(df_est_res, paste0("et-or",is_oracle,"-enc",is_encode,":"))
    
    return(list(myfit=myfit, df_est_res=df_est_res))
}
