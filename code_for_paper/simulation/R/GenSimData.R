## Generate simulated data
GenSimData <- function(repl, k, ni, xdim, xpdf, tau_fn, m_fn, noise_sd=sqrt(0.01)) { #sqrt(1)
    X <- GenSimX(ni, xdim, xpdf)
    names(X) <- paste0("X", seq(1,xdim))
    X[, tau := eval(tau_fn)]
    X[, m := eval(m_fn)]
    idx <- seq(1, ni)
    R <- rep(repl, ni)
    S <- rep(k, ni)
    Z <- as.integer(sample(c(rep(0, ni/2), rep(1, ni/2))))
    eps <- rnorm(ni, mean=0, sd=noise_sd)
    df <- data.table(R, S, idx, Z, X, eps)
    df[, Y := m + (2 * Z - 1) / 2 * tau + eps]
    df[, c("m", "eps") := NULL]
    
    role <- rep("train", ni)
    role[sample(1:ni, ni / 2, replace = FALSE)] <- "est"
    df$role <- role
    return(df)
}
