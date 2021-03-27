## Simulate X covariates
GenSimX <- function(ni, xdim, xpdf) {
    if (xpdf == "unif") {
        X <- matrix(runif(ni*xdim, 0, 1), nrow=ni, ncol=xdim)
    } else if (xpdf == "norm") {
        X <- matrix(rnorm(ni*xdim), nrow=ni, ncol=xdim)
    } else if (xpdf == "binary") {
        X <- matrix(rbinom(ni*xdim,1,0.5), ni, xdim)
    } else {
        print("error in simulating X!")
    }
    X <- data.table(X)
    return(X)
}