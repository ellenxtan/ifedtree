## Bagging (directly augment tau_hat from individual model)
Bag <- function(df, ind_all) {
    bag_est <- mean(ind_all$tau_hat, na.rm=TRUE)
    df[, tau_hat := bag_est]
    
    bag_out <- df[, c("R", "S", "tau", "tau_hat")]
    
    return(bag_out)
}
