#' Prediction plot of an ensemble forest for a grid of two variables
#'
#' Prediction plot of an ensemble forest with respect to two variables of interest
#' while fixing other variables at their means of the coordinating site data.
#'
#' @param aug_df    The augmented data frame used to fit an ensemble forest (`data.table`).
#' @param coord_df  The coordinating site data (`data.table`).
#' @param coord_id  Site index for coordinating site.
#' @param myfit     A fitted ensemble forest.
#' @param site      Variable name for site indicator.
#' @param covars    A vector of covariate names used.
#' @param var1      A character string with the name a numerical predictor that will on X-axis.
#' @param var2      A character string with the name a numerical predictor that will on Y-axis.
#' @param site_enc_tab  A data.table of mean outcome for each site. Default is NULL.
#'                      If class of myfit is a `grf`, "site_enc_tab" should not be NULL.
#' @param grids     The number of points on the one-dimensional grid on x and y-axis. Default is 100.
#' @param seed      A random seed for reproducing the figure.
#' @param midpoint  The midpoint (in data value) of the diverging scale. Default is 0.
#' @param plt_title Title of the plot. Default is "".
#'
#' @return A ggplot2 object.
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' data(SimDataLst)
#' K <- length(SimDataLst)
#' covars <- grep("^X", names(SimDataLst[[1]]), value=TRUE)
#' fit_lst <- list()
#' for (k in 1:K) {
#'     tmpdf <- SimDataLst[[k]]
#'     fit_lst[[k]] <- grf::causal_forest(X=as.matrix(tmpdf[, covars, with=FALSE]),
#'                                        Y=tmpdf$Y, W=tmpdf$Z)
#' }
#' coord_id <- 1
#' coord_df <- SimDataLst[[coord_id]]
#' aug_df <- GenAugData(coord_id, coord_df, fit_lst, covars)
#'
#' ## Treat each site as a distinct factor
#' myfit <- EnsemForest(coord_id, aug_df, "site", covars)$myfit
#' PlotForestPred(aug_df, coord_df, coord_id, myfit, "site", covars, "site", "X1")
#' PlotForestPred(aug_df, coord_df, coord_id, myfit, "site", covars, "X1", "site")
#' PlotForestPred(aug_df, coord_df, coord_id, myfit, "site", covars, "X1", "X2")
#'
#' ## Mean encoding as surrogate for site index
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars, is_encode=TRUE)
#' myfit <- res_ef$myfit
#' site_enc_tab <- res_ef$site_enc_tab
#' PlotForestPred(aug_df, coord_df, coord_id, myfit,
#'                "site", covars, "site", "X1", site_enc_tab)
#' PlotForestPred(aug_df, coord_df, coord_id, myfit,
#'                "site", covars, "X1", "site", site_enc_tab)
#' PlotForestPred(aug_df, coord_df, coord_id, myfit,
#'                "site", covars, "X1", "X2", site_enc_tab)
#'
#'
PlotForestPred <- function(aug_df, coord_df, coord_id, myfit,
                           site, covars, var1, var2, site_enc_tab=NULL,
                           grids=100, seed=1234, midpoint=0, plt_title="") {

    # checking
    aug_df <- CheckDF(aug_df, "aug_df")
    coord_df <- CheckDF(coord_df, "coord_df")
    CheckVar(coord_id, "coord_id")
    CheckMyfit(myfit)
    if ("grf" %in% class(myfit)) {
        CheckNull(site_enc_tab, "site_enc_tab")
    }
    CheckVar(site, "site")
    CheckVector(covars, "covars")
    CheckVar(var1, "var1")
    CheckVar(var2, "var2")


    covars_all <- c(site, covars)

    if ("grf" %in% class(myfit)) {
        covars_all_grf <- c("site_enc", covars)  # mean encoding # use for pred
    }

    if (class(aug_df[[var1]]) == "factor") {
        grid1 <- length(levels(aug_df[[var1]]))
        aug_df[[var1]] <- as.numeric(as.character(aug_df[[var1]]))
    } else {
        grid1 <- grids
    }
    if (class(aug_df[[var2]]) == "factor") {
        grid2 <- length(levels(aug_df[[var2]]))
        aug_df[[var2]] <- as.numeric(as.character(aug_df[[var2]]))
    } else {
        grid2 <- grids
    }

    newdata <- expand.grid(seq(min(aug_df[[var1]]), max(aug_df[[var1]]), length.out = grid1),
                           seq(min(aug_df[[var2]]), max(aug_df[[var2]]), length.out = grid2))
    colnames(newdata) <- c(var1, var2)

    other_vars <- setdiff(covars_all, c(var1, var2))

    for (usevar in other_vars) {  # assign other variables at their means
        newdata[[usevar]] <- NA_real_
        if (usevar == site) {
            newdata[[usevar]] <- coord_id
        } else {
            newdata[[usevar]] <- mean(coord_df[[usevar]])
        }
    }

    # order site by average Y_aug in the aug_df for better visualization
    Y_aug <- NULL  # bind variables locally to the function (define as NULL)
    site_mat <- aug_df[, mean(Y_aug), by=site]
    siteRank <- site_mat[order(site_mat$V1), ][[site]]
    newdata[[site]] <- factor(newdata[[site]], levels=siteRank)

    if ("ranger" %in% class(myfit)) {
        # convert factor variables
        for (var in covars) {
            if (is.factor(coord_df[[var]])) {
                newdata[[var]] <- factor(newdata[[var]])
            }
        }
        newdata$prediction <- predict(myfit, newdata[covars_all])$predictions
    }

    if ("grf" %in% class(myfit)) {
        # mean encoding for coord site
        newdata <- merge(newdata, site_enc_tab, by=site)
        newdata$prediction <- predict(myfit, newdata[covars_all_grf])$predictions
    }

    plt_pred <- ggplot(newdata, aes_string(x = var1, y = var2)) +
        geom_tile(aes(fill=.data$prediction))

    if (class(newdata[[var1]]) == "factor" & class(newdata[[var2]]) == "numeric") {  # site vs. X1/X2
        plt_pred <- plt_pred +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
    } else if (class(newdata[[var1]]) == "numeric" & class(newdata[[var2]]) == "factor") {  # site vs. X1/X2
        plt_pred <- plt_pred +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0))
    } else if (class(newdata[[var1]]) == "factor" & class(newdata[[var2]]) == "factor") {  # site vs. X1/X2
        plt_pred <- plt_pred +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0))
    } else {  # both continuous
        plt_pred <- plt_pred +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
    }
    plt_pred <- plt_pred +
        scale_fill_gradient2(midpoint = midpoint,
                             low = "blue", mid = 'white', high = "red") +
        theme_classic() +
        theme(legend.position="top",
              legend.title=element_text(size=13,face="bold"),
              legend.text=element_text(size=11,face="bold"),
              axis.text=element_text(size=12,face="bold"),
              axis.title=element_text(size=14,face="bold")) +
        guides(fill = guide_colorbar(barwidth=12, barlength=2))  # length of legend colorbar

    plt_pred[["labels"]][["fill"]] <- plt_title

    return(plt_pred)
}


