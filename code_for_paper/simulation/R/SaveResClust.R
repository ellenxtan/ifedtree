SaveResClust <- function(res, K, ni, honest, is_clust, is_save, out_path,
                    clustCF_name, is_test) {
  
  clustCF_out <- LstToDF(res, clustCF_name)
  
  out_name <- paste0("istest",substr(is_test,1,1),"_hon",substr(honest,1,1),
                     "_Clust",substr(is_clust,1,1))
  
  dir.create(file.path(out_path))  # create a new folder for each design
  out_path <- paste0(out_path, "/ClustCF")
  dir.create(file.path(out_path))  # create a new folder within a design
  
  if (is_save) {
    tic()
    cat("Saving ...\n")
    if (is_clust) {
      save(clustCF_out, file=paste0(out_path, "/K", K, "_ni", ni, "_clustCF_", out_name, ".RData"))
    }
    toc()
  }
  
  EvaResClust(K, ni, out_name, is_clust, out_path, clustCF_out, is_test)
  
}
