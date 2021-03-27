EvaRes <- function (K, ni, ind_model, out_name,
                    is_ens, is_meta, is_agg, out_path,
                    ind_out, et_out, ef_out, 
                    metaRE_out, agg_out, is_test) {
    
    ind_smry <- cbind(K=K, ni=ni, model=paste0("local", ind_model), EvaluateMetrics(ind_out, "ind_out"))
    if (is_ens) {
        et_smry <- cbind(K=K, ni=ni, model=paste0(ind_model,"-ET"), EvaluateMetrics(et_out, "et_out"))
        ef_smry <- cbind(K=K, ni=ni, model=paste0(ind_model,"-EF"), EvaluateMetrics(ef_out, "ef_out"))
        smry_out <- rbind(ef_smry, et_smry, ind_smry)
    } else {
        smry_out <- rbind(ind_smry)
    }
    if (is_meta) {
        metaRE_smry <- cbind(K=K, ni=ni, model=paste0(ind_model,"-metaRE"), EvaluateMetrics(metaRE_out, "metaRE_out"))
        smry_out <- data.table(rbind(smry_out, metaRE_smry))
    }
    if (is_agg) {
        agg_smry <- cbind(K=K, ni=ni, model=paste0(ind_model,"-agg"), EvaluateMetrics(agg_out, "agg_out"))
        smry_out <- data.table(rbind(smry_out, agg_smry))
    }
    save(smry_out, file=paste0(out_path, "/K", K, "_ni", ni, "_smry_", out_name, ".RData"))
    cat("out_name:", out_name, "\n")
    options(width = 500)
    print(smry_out, digits=4)
}
