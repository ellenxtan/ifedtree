## Generate simulated data
GenSimData <- function(k, ni, xdim, xpdf, tau_fn, m_fn, scale_c, grp_type, outcome, 
                       is_test=FALSE, site1_U=NULL) {
    
    X <- GenSimX(ni, xdim, xpdf)
    names(X) <- paste0("X", seq(1,xdim))
    
    # same group for each site
    if (!is_test) { # train data
        if (grp_type == "cont") {
            X$U <- runif(1)
        } else if (grp_type == "2grp") {
            X$U <- ifelse(k%%2==0, 1, 0)
        }
    } else { # test data
        X$U <- site1_U
    }
    
    X$scale_c <- scale_c
    X[, tau := eval(tau_fn)]
    X[, m := eval(m_fn)]
    
    Z <- as.integer(sample(c(rep(0, ni/2), rep(1, ni/2))))
    
    idx <- seq(1, ni)
    S <- rep(k, ni)
    eps <- rnorm(ni, mean=0, sd=sqrt(0.01)) #sqrt(1)
    df <- data.table(S, idx, Z, X, eps)  #R, 
    
    df[, Y := m + (2 * Z - 1) / 2 * tau + eps]
    
    df[, c("m", "eps") := NULL]
    
    # train-est 1:1
    role <- rep("train", ni)
    role[sample(1:ni, ni / 2, replace = FALSE)] <- "est"
    df$role <- role
    
    return(df)
}
