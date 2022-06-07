## Evaluate bias & mse from the above methods
EvaluateMetrics <- function(df, df_name) {
    df$bias <- df$tau_hat - df$tau
    df$mse <- (df$bias)^2
    
    df <- data.table(df)
    eva_df <- df[, c("S", "bias", "mse")] #"R", 
    eva_smry <- eva_df[, list(bias=mean(bias, na.rm=T),
                              bias_sd=sd(bias, na.rm=T),
                              mse=mean(mse, na.rm=T),
                              mse_sd=sd(mse, na.rm=T)
                              )]
    
    cat(df_name, "bias", round(eva_smry$bias,3), "MSE", round(eva_smry$mse,3), "\n")
    return(eva_smry)
}
