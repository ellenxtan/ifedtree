## Main function
MainClust <- function(design, R, K, ni, xdim, xpdf, tau_fn, m_fn, honest, is_clust, 
                      is_save, out_path, core_num) {
  
  message("--------------------------")
  cat(paste0("\nDesign ",design,": K=", K, "; ni=", ni, "\n"))
  message(paste0("Design ",design,": K=", K, "; ni=", ni))
  out_path = paste0(out_folder, "/", design)
  message("tau_fn: ", tau_fn)
  cat("tau_fn: "); print(tau_fn)
  
  cat("=======================================================================\n")
  cat(paste0("R=", R, "; K=", K, "; ni=", ni, "; xdim=", xdim, 
             "; honest=", honest, "; is_clust=", is_clust, "\n"))
  message(paste0("R=", R, "; K=", K, "; ni=", ni, "; xdim=", xdim, 
             "; honest=", honest, "; is_clust=", is_clust))
  tic()
  cat("Mclapplying ...\n")
  res <- mclapply(1:R, SimulationClust, K=K, ni=ni, xdim=xdim, xpdf=xpdf,
                  tau_fn=tau_fn, m_fn=m_fn, 
                  honest=honest, is_clust=is_clust, 
                  mc.cores=core_num)
  toc()
  
  
  cat("=======================================\n")
  print("Prediction")
  SaveResClust(res, K, ni, honest, is_clust, is_save, out_path,
          "clustCF_test_out", is_test=TRUE)
  
  cat("=======================================================================\n")
  cat("=======================================================================\n")
  
  message("--------------------------\n")
}
