## Simulation
SimulationClust <- function(repl, K, ni, xdim, xpdf, tau_fn, m_fn, honest, is_clust) {
  coord <- 1
  n_test <- ni
  honest2 <- honest
  
  
  ############################### Train/Test data ##############################
  
  set.seed(repl)
  
  ## get test set (belong to site1)
  test_df0 <- GenSimData(repl, k=coord, ni=n_test, xdim, xpdf, tau_fn, m_fn)
  
  df_lst <- list()
  for (k in 1:K) {
    df <- GenSimData(repl, k, ni, xdim, xpdf, tau_fn, m_fn)
    df_lst[[k]] <- copy(df)
  }
  
  
  ################################### ClustCF ##################################
  
  # cluster-robust cf
  clustCF_test_out <- c()
  if (is_clust) {
    df_all <- purrr::map_dfr(df_lst, rbind)
    covars <- grep("^X", names(df_all), value=TRUE)
    tmp <- CentClustCF(df_all, covars, honest=honest2, is_pred=FALSE)
    rm(df_all)
    
    ## predict on test set (belong to site1)
    clustCF_test_out <- CentClustCF(test_df0, covars, honest=honest2,
                                    is_pred=TRUE, myfit=tmp$myfit)$df_est_res
    rm(tmp)
  }
  
  
  ##############################################################################
  
  if (repl %% 50 == 0) { message(repl) }
  
  return(list(clustCF_test_out=clustCF_test_out))
  
}
