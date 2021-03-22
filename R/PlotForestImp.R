#' Importance plot of an ensemble forest
#'
#' @param myfit  A fitted ensemble forest.
#'
#' @return A ggplot2 object and a table of relative importance of each covariate.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
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
#' ## Treat each site as a distinct factor
#' myfit <- EnsemForest(coord_id, aug_df, "site", covars, importance="impurity")$myfit
#' PlotForestImp(myfit)
#'
#' ## Mean encoding as surrogate for site index
#' myfit <- EnsemForest(coord_id, aug_df, "site", covars, is_encode=TRUE)$myfit
#' PlotForestImp(myfit)
#'
#'
PlotForestImp <- function(myfit) {

    # checking
    CheckMyfit(myfit)

    # bind variables locally to the function (define as NULL)
    . <- prop <- NULL

    if ("ranger" %in% class(myfit)) {
        tab_imp <- myfit$variable.importance %>%
            as.data.frame() %>%
            dplyr::arrange(dplyr::desc(.)) %>%
            dplyr::mutate(prop = . / sum(.)) %>%
            tibble::rownames_to_column(var = "names") %>%
            dplyr::select(-.)
    }

    if ("grf" %in% class(myfit)) {
        tab_imp <- as.data.frame(cbind(names(myfit$X.orig),
                                       grf::variable_importance(myfit)))
        names(tab_imp) <- c("names", "prop")
        tab_imp$prop <- as.numeric(tab_imp$prop)
    }

    plt_imp <- tab_imp %>% ggplot(aes(reorder(names, prop), prop)) +
        geom_col() +
        coord_flip() +
        ylab("Relative importance") + xlab("Covariates") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() +
        theme(axis.text=element_text(size=13,face="bold"),
              axis.title=element_text(size=14,face="bold"))

    return(list(plt_imp=plt_imp, tab_imp=tab_imp))
}
