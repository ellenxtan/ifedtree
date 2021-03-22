#' Plot an ensemble tree
#'
#' @param myfit  A fitted ensemble tree.
#' @param \\dots Additional arguments for plotting the tree.
#'
#' @import rpart.plot
#' @export
#'
#' @examples
#'
#' data(SimDataLst)
#' K <- length(SimDataLst)
#' covars <- grep("^X", names(SimDataLst[[1]]), value=TRUE)
#' fit_lst <- list()
#' for (k in 1:K) {
#'     tmpdf <- SimDataLst[[k]]
#'     # use your estimator of interest
#'     fit_lst[[k]] <- grf::causal_forest(X=as.matrix(tmpdf[, covars, with=FALSE]),
#'                                        Y=tmpdf$Y, W=tmpdf$Z)
#' }
#'
#' coord_id <- 1
#' coord_df <- SimDataLst[[coord_id]]
#' aug_df <- GenAugData(coord_id, coord_df, fit_lst, covars)
#'
#' myfit <- EnsemTree(coord_id, aug_df, "site", covars)$myfit
#' PlotTree(myfit, main="Ensemble Tree")
#'
#'
PlotTree <- function(myfit, ...) {

    # checking
    CheckMyfit(myfit)

    rpart.plot(myfit, ...)
}
