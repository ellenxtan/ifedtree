## Simulation
Simulation <- function(repl, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
                       honest, is_meta, is_agg) {
    coord <- 1
    n_test <- ni
    honest1 <- TRUE
    honest2 <- honest
    
    
    ########################## Train/Test & local model ########################
    
    set.seed(repl)
    
    ## get test set (belong to site1)
    test_df0 <- GenSimData(repl, k=coord, ni=n_test, xdim, xpdf, tau_fn, m_fn)
    
    ## individual model
    df_lst   <- list()
    cfit_lst <- list()
    ind_out  <- c()  # site 1
    ind_all <- c()
    for (k in 1:K) {
        df <- GenSimData(repl, k, ni, xdim, xpdf, tau_fn, m_fn)
        df_lst[[k]] <- df
        covars <- grep("^X", names(df), value=TRUE)  # include X's
        if (ind_model == "CT") {
            tmp <- IndCT(df, covars, honest=honest1, is_pred=FALSE)
        } else if (ind_model == "CF") {
            tmp <- IndCF(df, covars, honest=honest1, is_pred=FALSE)
        } else if (ind_model == "X-BART") {
            tmp <- IndXBART(df, covars, honest=honest1, is_pred=FALSE)
        }
        cfit_lst[[k]] <- tmp$myfit
        if (k == coord) { ind_out <- data.table(tmp$df_est_res) }  # site 1 ind_out
        ind_all <- data.table(rbind(ind_all, tmp$df_est_res))
        rm(df, tmp)
    }
    # pred on site 1 model
    if (ind_model == "CT") {
        ind_test_out <- IndCT(test_df0, covars, honest=honest1, is_pred=TRUE, myfit=cfit_lst[[coord]])$df_est_res
    } else if (ind_model == "CF") {
        ind_test_out <- IndCF(test_df0, covars, honest=honest1, is_pred=TRUE, myfit=cfit_lst[[coord]])$df_est_res
    } else if (ind_model == "X-BART") {
        ind_test_out <- IndXBART(test_df0, covars, honest=honest1, is_pred=TRUE, myfit=cfit_lst[[coord]])$df_est_res
    }
    
    
    ############################## Augmented data ##############################
    
    ## augmented data (generated from one site)
    if (ind_model == "CT") {
        aug_df <- AugDataCT(df_lst[[coord]], cfit_lst, tau_fn, K)
    } else if (ind_model == "CF") {
        aug_df <- AugDataCF(df_lst[[coord]], cfit_lst, tau_fn, K, coord=coord, is_train=TRUE)
    } else if (ind_model == "X-BART") {
        aug_df <- AugDataXl(df_lst[[coord]], cfit_lst, tau_fn, K)
    }
    
    
    ############################## Ensemble model ##############################
    
    ## ensemble model
    et_out <- c()  # tree
    ef_out <- c()  # forest
    et_test_out <- c()  # tree
    ef_test_out <- c()  # forest
    if (is_ens) {
        # ensemble tree
        ens_et_res <- EnsRT(aug_df, honest=honest2, is_pred=FALSE)
        et_out <- data.table(ens_et_res$df_est_res)
        et_out <- et_out[site==coord,]
        
        # ensemble forest
        ens_ef_res <- EnsRF(aug_df, honest=honest2, is_pred=FALSE)
        ef_out <- data.table(ens_ef_res$df_est_res)
        ef_out <- ef_out[site==coord,]  # OOB pred
        
        ## prediction on test set (need `site` variable but don't need `leaves`)
        ## (use ensemble model with site1 as coord site)
        test_df2 <- copy(test_df0)
        test_df2$site <- coord
        test_df2$site <- factor(test_df2$site)
        levels(test_df2$site) = levels(aug_df$site)
        # ensemble tree
        et_test_out <- EnsRT(test_df2, honest=honest2, is_pred=TRUE, myfit=ens_et_res$myfit, 
                             leaves_tab=ens_et_res$leaves_tab)$df_est_res
        # ensemble forest
        aug_df[, sitenum := round(mean(leaves)), by=site]  # mean encoding
        test_df2$sitenum = aug_df[site==coord,]$sitenum[1]  # coord
        ef_test_out <- EnsRF(test_df2, honest=honest2, is_pred=TRUE, myfit=ens_ef_res$myfit)$df_est_res
        rm(ens_et_res, ens_ef_res)
        
    }
    
    
    ############################## Meta-analysis ###############################
    
    ## augmented data from test set (for meta-analysis)
    if (is_meta) {
        if (ind_model == "CT") {
            aug_test_df <- AugDataCT(test_df0, cfit_lst, tau_fn, K)
        } else if (ind_model == "CF") {
            aug_test_df <- AugDataCF(test_df0, cfit_lst, tau_fn, K, coord=coord, is_train=FALSE)
        } else if (ind_model == "X-BART") {
            aug_test_df <- AugDataXl(test_df0, cfit_lst, tau_fn, K)
        }
    }
    
    ## meta analysis
    metaRE_out <- c()
    metaRE_test_out <- c()
    if (is_meta) {
        # random-effect meta
        meta_df <- aug_df[, c("R", "S", "idx", "tau", "site", "leaves", "leavesse")]
        meta_res <- MetaRE(meta_df)
        metaRE_out <- data.table(rbind(metaRE_out, meta_res))
        rm(meta_df, meta_res)
        
        # random-effect meta (on test set)
        meta_df <- aug_test_df[, c("R", "S", "idx", "tau", "site", "leaves", "leavesse")]
        meta_res <- MetaRE(meta_df)
        metaRE_test_out <- data.table(meta_res)
        rm(meta_df, meta_res)
    }
    
    
    ############################# Mean aggregation #############################
    
    ## Mean aggregation (avg among K sites)
    agg_out <- c()
    agg_test_out <- c()
    if (is_agg) {
        agg_out <- Bag(df_lst[[coord]], ind_all)
        agg_test_out <- Bag(test_df0, ind_all)
    }
    
    
    ############################################################################
    
    if (repl %% 50 == 0) { message(repl) }
    
    return(list(ind_out=ind_out, 
                et_out=et_out, ef_out=ef_out,  
                metaRE_out=metaRE_out, agg_out=agg_out,
                
                ind_test_out=ind_test_out, 
                et_test_out=et_test_out, ef_test_out=ef_test_out,
                metaRE_test_out=metaRE_test_out, agg_test_out=agg_test_out
    ))
}
