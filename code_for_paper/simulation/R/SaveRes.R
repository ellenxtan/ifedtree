SaveRes <- function(res, K, ni, ind_model, is_ens,
                    honest, is_meta, is_agg, is_save, out_path,
                    ind_name, et_name, ef_name, 
                    metaRE_name, agg_name, is_test) {
    
    ind_out <- LstToDF(res, ind_name)
    et_out <- LstToDF(res, et_name)
    ef_out <- LstToDF(res, ef_name)
    metaRE_out <- LstToDF(res, metaRE_name)
    agg_out <- LstToDF(res, agg_name)
    
    out_name <- paste0("istest",substr(is_test,1,1),"_loc",ind_model, "_hon",substr(honest,1,1),
                       "_Ens",substr(is_ens,1,1),
                       "_M",substr(is_meta,1,1),"_B",substr(is_agg,1,1))
    dir.create(file.path(out_path))  # create a new folder for each design
    if (is_save) {
        tic()
        cat("Saving ...\n")
        save(ind_out, file=paste0(out_path, "/K", K, "_ni", ni, "_ind_", out_name, ".RData"))
        save(et_out, file=paste0(out_path, "/K", K, "_ni", ni, "_et_", out_name, ".RData"))
        save(ef_out, file=paste0(out_path, "/K", K, "_ni", ni, "_ef_", out_name, ".RData"))
        save(metaRE_out, file=paste0(out_path, "/K", K, "_ni", ni, "_metaRE_", out_name, ".RData"))
        save(agg_out, file=paste0(out_path, "/K", K, "_ni", ni, "_agg_", out_name, ".RData"))
        toc()
    }
    
    EvaRes(K, ni, ind_model, out_name,
           is_ens, is_meta, is_agg, out_path,
           ind_out, et_out, ef_out,
           metaRE_out, agg_out, is_test)
    
}
