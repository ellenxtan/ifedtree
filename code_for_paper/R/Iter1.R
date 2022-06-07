Iter1 <- function(K, ni, xdim, xpdf, honest1, honest2, loc.mod, 
                  tau_fn, m_fn, scale_c, coord, n_test, grp_type, outcome) {
    
    if(!exists(".Random.seed")) set.seed(NULL)
    mytime <- as.numeric(gsub(":","",strsplit(as.character(Sys.time()), " ")[[1]][2]))
    set.seed(mytime)
    
    ############################### Local models ###############################
    
    ## local model for augmented data
    df_lst   <- list()
    cfit_lst <- list()
    for (k in 1:K) {
        df <- GenSimData(k, ni, xdim, xpdf, tau_fn, m_fn, scale_c, grp_type, outcome)
        df_lst[[k]] <- df
        
        covars <- grep("^X", names(df), value=TRUE)  # include X's
        
        if (k == coord) { # sample split for site 1 (fit 2 local models)
            # estimation set model ("true")
            loc.fit <- IndFit(loc.mod=loc.mod, pred_only=FALSE,
                              df_use=copy(df[role=="est", ]), is_pred=FALSE, 
                              covars=covars, honest1=honest1, myfit=NULL)$myfit
            
            df <- copy(df[role=="train", ])  # train set to fit model
        }
        
        cfit_lst[[k]] <- IndFit(loc.mod=loc.mod, pred_only=FALSE,
                                df_use=copy(df), is_pred=FALSE, 
                                covars=covars, honest1=honest1, myfit=NULL)$myfit
    }
    
    ################################ Test data #################################
    
    ## get test set (same as site1) 
    site1_U <- unique(df_lst[[coord]]$U) # U extracted from df1 (not randomly)
    test_df0 <- GenSimData(k=coord, ni=n_test, xdim, xpdf, 
                           tau_fn, m_fn, scale_c, grp_type, outcome,
                           is_test=TRUE, site1_U=site1_U)

    # local model predict on test data
    ind_out <- IndFit(loc.mod=loc.mod, pred_only=FALSE,
                      df_use=copy(test_df0), is_pred=TRUE, 
                      covars=covars, honest1=honest1, myfit=cfit_lst[[coord]])$df_est_res
    
    EvaluateMetrics(ind_out, "ind:")
    

    ############################# MA & Stacking ################################
    
    ewma_out <- EWMAwrap(loc.mod, test_df0, df_lst, cfit_lst, coord, loc.fit, K,
                         is_oracle=FALSE)$df_est_res
    
    ewma_or_out <- EWMAwrap(loc.mod, test_df0, df_lst, cfit_lst, coord, loc.fit, K,
                            is_oracle=TRUE)$df_est_res

    ma_out <- MAwrap(loc.mod, test_df0, df_lst, cfit_lst, K)$df_est_res
    
    stack_out <- LinearStack(loc.mod, test_df0, df_lst, cfit_lst, coord, loc.fit, is_oracle=FALSE)$df_est_res

    stack_or_out <- LinearStack(loc.mod, test_df0, df_lst, cfit_lst, coord, loc.fit, is_oracle=TRUE)$df_est_res


    ############################## Augmented data ##############################
    
    ## augmented data (generated from site 1) estmation set to build
    if (loc.mod == "CT") {
        aug_df <- AugDataCT(df_lst[[coord]][role=="est", ], df_lst, cfit_lst, K, tau_fn)
    } else if (loc.mod == "CF") {
        aug_df <- AugDataCF(df_lst[[coord]][role=="est", ], df_lst, cfit_lst, K, tau_fn)
    }
    
    if (grp_type == "2grp") {
        aug_df$site_real <- factor(aug_df$site_real)
    }
    
    ############################## Ensemble model ##############################
    
    ## ensemble model
    # prediction on test set (need `site` variable but don't need `leaves`)
    # (use ensemble model with site1 as coord site)
    test_df0$site <- coord
    test_df0$site <- factor(test_df0$site)
    levels(test_df0$site) = levels(aug_df$site)
    
    test_df0$site_real <- test_df0$U
    if (grp_type == "2grp") {
        test_df0$site_real <- factor(test_df0$site_real)
        levels(test_df0$site_real) = levels(aug_df$site_real)
    }
    
    # ensemble tree
    ens_et_res <- EnsRT(aug_df=aug_df, is_oracle=FALSE, test_df=test_df0, 
                        honest=honest2, is_encode=TRUE)
    et_out <- ens_et_res$df_est_res
    
    # ensemble forest
    ef_out <- EnsRF(aug_df=aug_df, is_oracle=FALSE, test_df=test_df0, honest=honest2, 
                    pkg="grf", coord=coord)$df_est_res

    ########################### Ensemble oracle model ##########################
    
    # ensemble tree
    ens_et_res <- EnsRT(aug_df=aug_df, is_oracle=TRUE, test_df=test_df0, 
                        honest=honest2, is_encode=TRUE)
    et_or_out <- ens_et_res$df_est_res
    
    # ensemble forest
    ef_or_out <- EnsRF(aug_df=aug_df, is_oracle=TRUE, test_df=test_df0, honest=honest2, 
                       pkg="grf", coord=coord)$df_est_res
    
    ################################## Results ##################################
    
    cat("====================================\n")
    eva_ind <- EvaluateMetrics(ind_out, "ind:")
    eva_ma <- EvaluateMetrics(ma_out, "ma:")
    eva_ewma <- EvaluateMetrics(ewma_out, "ewma:")
    eva_ewma_or <- EvaluateMetrics(ewma_or_out, "ewma-oracle:")
    eva_stack <- EvaluateMetrics(stack_out, "stack:")
    eva_stack_or <- EvaluateMetrics(stack_or_out, "stack-oracle:")
    eva_et <- EvaluateMetrics(et_out, "et:")
    eva_et_or <- EvaluateMetrics(et_or_out, "et-oracle:")
    eva_ef <- EvaluateMetrics(ef_out, "ef:")
    eva_ef_or <- EvaluateMetrics(ef_or_out, "ef-oracle:")
    cat("====================================\n")
    
    return(list(eva_ind=eva_ind, eva_ma=eva_ma, 
                eva_ewma=eva_ewma, eva_ewma_or=eva_ewma_or,
                eva_stack=eva_stack, eva_stack_or=eva_stack_or,
                eva_et=eva_et, eva_et_or=eva_et_or, 
                eva_ef=eva_ef, eva_ef_or=eva_ef_or
                )
            )
}
