#' Generate simulated data
#'
#' @param k      Site index.
#' @param ni     Sample size within a site.
#' @param xdim   Dimension of covariate `X`.
#' @param xpdf   Distribution of covariate `X`.
#' @param sd     Standard deviation of error term in simulating outcome.
#' @param tau_fn Treatment effect function.
#' @param m_fn   Mean effect function.
#'
#' @return A simulated data set.
#' @export
#' @import data.table
#'
#' @examples
#' GenSimData(k=1)
#'
GenSimData <- function(k, ni=100, xdim=5, xpdf="norm", sd=sqrt(0.01),
                       tau_fn = "X1*(X1>0) -4*(k%%2==0)+X1*(k%%2==0)",
                       m_fn = "0") {

    # checking
    CheckNumeric(k, "k")
    CheckNumeric(ni, "ni")
    CheckNumeric(xdim, "xdim")
    CheckPDF(xpdf)
    CheckNumeric(sd, "sd")
    CheckCharacter(tau_fn, "tau_fn")
    CheckCharacter(m_fn, "m_fn")


    tau_fn <- parse(text = tau_fn)
    m_fn <- parse(text = m_fn)

    # bind variables locally to the function (define as NULL)
    tau <- m <- Y <- NULL

    X <- GenSimX(ni, xdim, xpdf)
    names(X) <- paste0("X", seq(1,xdim))
    X[, tau := eval(tau_fn)]
    X[, m := eval(m_fn)]
    idx <- seq(1, ni)
    S <- rep(k, ni)
    Z <- as.integer(sample(c(rep(0, ni/2), rep(1, ni/2))))
    eps <- rnorm(ni, mean=0, sd=sd)
    df <- data.table(S, idx, Z, X, eps)
    df[, Y := m + (2 * Z - 1) / 2 * tau + eps]
    df[, c("m", "eps") := NULL]

    return(df)
}


#' Simulate X covariates
#'
#' @param ni   Sample size within a site.
#' @param xdim Dimension of covariate `X`.
#' @param xpdf Distribution of covariate `X`.
#'
#' @return Covariate matrix `X`.
#' @export
#' @import data.table
#'
#' @examples
#' GenSimX()
#'
GenSimX <- function(ni=100, xdim=5, xpdf="norm") {

    # checking
    CheckNumeric(ni, "ni")
    CheckNumeric(xdim, "xdim")
    CheckPDF(xpdf)

    if (xpdf == "unif") {
        X <- matrix(runif(ni*xdim, 0, 1), nrow=ni, ncol=xdim)
    } else if (xpdf == "norm") {
        X <- matrix(rnorm(ni*xdim), nrow=ni, ncol=xdim)
    } else if (xpdf == "binary") {
        X <- matrix(rbinom(ni*xdim,1,0.5), ni, xdim)
    } else {
        print("Distribution of X not supported!")
    }

    X <- data.table(X)
    return(X)
}
