## Main function
Main <- function(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, ind_model, is_ens,
                 honest, is_meta, is_agg, is_save, out_path, core_num) {
    
    message("--------------------------")
    cat(paste0("\nDesign ",design,": K=", K, "; ni=", ni, "\n"))
    message(paste0("Design ",design,": K=", K, "; ni=", ni))
    out_path = paste0(out_folder, "/", design)
    message("tau_fn: ", tau_fn)
    cat("tau_fn: "); print(tau_fn)
    
    cat("=======================================================================\n")
    cat(paste0("R=", R, "; K=", K, "; ni=", ni, "; xdim=", xdim, 
               "; ind=", ind_model, "; honest=", honest, "; is_ens=", is_ens,
               "; is_meta=", is_meta, "; is_agg=", is_agg, "\n"))
    message(paste0("R=", R, "; K=", K, "; ni=", ni, "; xdim=", xdim, 
               "; ind=", ind_model, "; honest=", honest, "; is_ens=", is_ens,
               "; is_meta=", is_meta, "; is_agg=", is_agg))
    tic()
    cat("Mclapplying ...\n")
    res <- mclapply(1:R, Simulation, K=K, ni=ni, xdim=xdim, xpdf=xpdf,
                    tau_fn=tau_fn, m_fn=m_fn, ind_model=ind_model, is_ens=is_ens,
                    honest=honest, is_meta=is_meta, is_agg=is_agg, mc.cores=core_num)
    toc()
    
    cat("=======================================\n")
    print("Prediction")
    SaveRes(res, K, ni, ind_model, is_ens,
            honest, is_meta, is_agg, is_save, out_path,
            "ind_test_out", "et_test_out", "ef_test_out",
            "metaRE_test_out", "agg_test_out", is_test=TRUE)
    
    cat("=======================================================================\n")
    cat("=======================================================================\n")
}
