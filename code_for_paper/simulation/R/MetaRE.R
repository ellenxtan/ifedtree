## Random-effects meta analysis
MetaRE <- function(meta_df) {
    meta_df <- dcast(meta_df, ... ~ site, value.var = c("leaves", "leavesse"))
    col_leaves <- which(names(meta_df) %in% grep("^leaves_", names(meta_df), value=TRUE))
    col_leavesse <- which(names(meta_df) %in% grep("^leavesse_", names(meta_df), value=TRUE))
    
    meta_res <- data.table((t(apply(meta_df, MARGIN=1, function(row) {
        re.pred <- suppressWarnings(predict(rma(yi=as.numeric(row[col_leaves]), sei=as.numeric(row[col_leavesse]),
                                                control=list(stepadj=0.5))))
        tau_hat <- re.pred$pred
        tau_lb <- re.pred$cr.lb
        tau_ub <- re.pred$cr.ub
        return(c(tau_hat, tau_lb, tau_ub))
    }))))
    names(meta_res) = c("tau_hat", "tau_lb", "tau_ub")
    meta_res <- cbind(meta_df[, c("R", "S", "tau")], meta_res)
    meta_res <- meta_res[, c("R", "S", "tau", "tau_hat")]
    return(meta_res)
}
