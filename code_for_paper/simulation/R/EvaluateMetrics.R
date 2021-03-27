## Evaluate bias & mse from the above methods
EvaluateMetrics <- function(df, df_name) {
    df$bias <- df$tau_hat - df$tau
    df$mse <- (df$bias)^2
    
    message(df_name)
    df <- data.table(df)
    eva_df <- df[, c("R", "S", "bias", "mse")]
    eva_smry <- eva_df[, list(bias=mean(bias, na.rm=T),
                              bias_sd=sd(bias, na.rm=T),
                              mse=mean(mse, na.rm=T),
                              mse_sd=sd(mse, na.rm=T)
                              )]

    return(eva_smry)
}
