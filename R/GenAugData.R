#' Generate Augmented Data
#'
#' Generate augmented data from coordinating site and localized models.
#' Calculate OOB prediction for coordinating site data on its model and common prediction on other models.
#'
#' @param coord_id   Index of coordinating site.
#' @param coord_df   Data of coordinating site (`data.table`).
#' @param fit_lst    A list of localized models.
#' @param covars     A vector of covariate names used.
#'
#' @return Augmented data generated from the coordinating site and localized models.
#' @export
#' @import data.table
#'
#' @examples
#' K <- 4
#' df_lst <- list()
#' fit_lst <- list()
#' for (k in 1:K) {
#'     df <- GenSimData(k)
#'     df_lst[[k]] <- df
#'     covars <- grep("^X", names(df), value=TRUE)
#'     fit_lst[[k]] <- grf::causal_forest(X=as.matrix(df[, covars, with=FALSE]), Y=df$Y, W=df$Z)
#' }
#' coord_id <- 1
#' coord_df <- df_lst[[coord_id]]
#' GenAugData(coord_id, coord_df, fit_lst, covars)
#'
#'
GenAugData <- function(coord_id, coord_df, fit_lst, covars) {

    # checking
    CheckVar(coord_id, "coord_id")
    coord_df <- CheckDF(coord_df, "coord_df")
    CheckList(fit_lst, "fit_lst")
    CheckVector(covars, "covars")


    # bind variables locally to the function (define as NULL)
    Y_aug <- NULL

    aug_df <- c()
    for (k in 1:length(fit_lst)) {
        tmpdf = copy(coord_df)

        if (k==coord_id) {  # coordinating site use OOB pred
            tmpdf[, Y_aug := fit_lst[[k]]$predictions]
        } else {
            tmpfit <- predict(fit_lst[[k]], as.matrix(tmpdf[, covars, with=FALSE]))  # df[, ..covars] = df[, covars, with=F]
            tmpdf[, Y_aug := tmpfit$predictions]
        }

        aug_df <- data.table(rbind(aug_df, cbind(tmpdf, site = k)))
    }

    aug_df$site <- factor(aug_df$site)
    return(aug_df)
}
