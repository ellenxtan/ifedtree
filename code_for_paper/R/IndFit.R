## local model for tau: CT or CF
IndFit <- function (loc.mod, pred_only, df_use, is_pred, covars, honest1, myfit) {

    if (!pred_only) {
        if (loc.mod == "CT") {
            res_ind <- IndCT(df_use, covars, honest=honest1, 
                           is_pred=is_pred, myfit=myfit)
        } else if (loc.mod == "CF") {
            res_ind <- IndCF(df_use, covars, honest=honest1, 
                           is_pred=is_pred, myfit=myfit)
        }
        
        return(res_ind)
        
    } else { # only predict function
        if (loc.mod == "CT") {
            res_vec <- predict(myfit, newdata=df_use)
        } else if (loc.mod == "CF") {
            res_vec <- predict(myfit, newdata=as.matrix(df_use[,..covars]))$predictions
        }
        
        return(res_vec)
        
    }
    
}
