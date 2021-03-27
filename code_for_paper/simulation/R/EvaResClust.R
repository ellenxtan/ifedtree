EvaResClust <- function (K, ni, out_name, is_clust, out_path, clustCF_out, is_test) {
  
  if (is_clust) {
    clustCF_smry <- cbind(K=K, ni=ni, model="clustCF", EvaluateMetrics(clustCF_out, "clustCF_out"))
    smry_out <- rbind(clustCF_smry)
  }
  
  save(smry_out, file=paste0(out_path, "/K", K, "_ni", ni, "_smry_", out_name, ".RData"))
  
  cat("out_name:", out_name, "\n")
  options(width = 500)
  print(smry_out, digits=4)
}
