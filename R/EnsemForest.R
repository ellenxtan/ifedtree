#' Ensemble Forest (EF)
#'
#' Build and/or predict on an ensemble regression forest.
#' Two implementation are provided for fitting the forest.
#' One treats each site as a distinct factor (implemented with ranger);
#' Another uses mean encoding for site index (implemented with grf).
#'
#' @param coord_id       Site index for coordinating site.
#' @param aug_df         The augmented data frame used to fit an ensemble forest (`data.table`).
#' @param site           Variable name for site indicator.
#' @param covars         A vector of covariate names used.
#' @param honest         Whether to use honest splitting (i.e., subsample splitting).
#'                       Default is FALSE.
#' @param is_pred        Whether to build an ensemble forest or make prediction.
#'                       Default is FALSE.
#' @param is_encode      Whether to treat each site as a distinct factor
#'                       or use mean encoding as surrogate for site index.
#'                       Useful when the number of underlying groups are known.
#'                       Default is FALSE
#' @param myfit          A fitted ensemble forest (for prediction purpose). Default is NULL.
#' @param est_leaves     A matrix of the number of observations in the augmented data
#'                       times the number of trees for assignment of terminal nodes of each tree
#'                       for the honest sample/estimation set (for honest prediction purpose).
#'                       Default is NULL.
#'                       If "honest" is set to FALSE, "est_leaves" is NULL;
#'                       If both "is_pred" and "honest" are TRUE, "est_leaves" should not be NULL.
#' @param honest_y       A vector of honest estimates for the honest sample/estimation set
#'                       (for honest prediction purpose). Default is NULL.
#'                       If "honest" is set to FALSE, "honest_y" is NULL;
#'                       If both "is_pred" and "honest" are TRUE, "honest_y" should not be NULL.
#' @param site_enc_tab   A data.table of mean outcome for each site. Default is NULL.
#'                       If both "is_pred" and "is_encode" are set to TRUE,
#'                       "site_enc_tab" should not be NULL.
#' @param \\dots         Additional arguments for building the forest.
#'
#' @return Training: return a fitted ensemble forest and OOB predictions of the input data;
#'         Prediction: return predictions of the input data.
#' @export
#' @import ranger
#' @import grf
#' @import data.table
#' @importFrom Rcpp sourceCpp
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
#' coord_test <- GenSimData(coord_id)
#'
#' coord_df <- SimDataLst[[coord_id]]
#' aug_df <- GenAugData(coord_id, coord_df, fit_lst, covars)
#'
#' ## Treat each site as a distinct factor
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars)
#' ef_hat <- EnsemForest(coord_id, coord_test, "site", covars, is_pred=TRUE,
#'             myfit=res_ef$myfit, est_leaves=res_ef$est_leaves, honest_y=res_ef$honest_y)
#'
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars, honest=TRUE)
#' ef_hat <- EnsemForest(coord_id, coord_test, "site", covars, honest=TRUE, is_pred=TRUE,
#'             myfit=res_ef$myfit, est_leaves=res_ef$est_leaves, honest_y=res_ef$honest_y)
#'
#' ## Mean encoding as surrogate for site index
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars, is_encode=TRUE)
#' ef_hat <- EnsemForest(coord_id, coord_test, "site", covars, is_pred=TRUE, is_encode=TRUE,
#'             myfit=res_ef$myfit, site_enc_tab=res_ef$site_enc_tab)
#'
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars, honest=TRUE, is_encode=TRUE)
#' ef_hat <- EnsemForest(coord_id, coord_test, "site", covars, honest=TRUE,
#'             is_pred=TRUE, is_encode=TRUE,
#'             myfit=res_ef$myfit, site_enc_tab=res_ef$site_enc_tab)
#'
#'
EnsemForest <- function(coord_id, aug_df, site, covars, honest=FALSE,
                        is_pred=FALSE, is_encode=FALSE,
                        myfit=NULL, est_leaves=NULL, honest_y=NULL,
                        site_enc_tab=NULL, ...) {

    # checking
    CheckVar(coord_id, "coord_id")
    aug_df <- CheckDF(aug_df, "aug_df")
    CheckVar(site, "site")
    CheckVector(covars, "covars")
    if (is_pred==T & honest==T & is_encode==F) {
        CheckMyfit(myfit)
        CheckNull(est_leaves, "est_leaves")
        CheckNull(honest_y, "honest_y")
    }
    if (is_pred==T & is_encode==T) {
        CheckMyfit(myfit)
        CheckNull(site_enc_tab, "site_enc_tab")
    }


    # bind variables locally to the function (define as NULL)
    Y_aug <- site_enc <- NULL

    ## treat each site as a distinct factor (with ranger)
    if (!is_encode) {
        fml <- as.formula(paste("Y_aug ~ ", site, " + ", paste(covars, collapse="+")))

        if (!is_pred) {  # for estimation
            if (!honest) {  # adaptive EF
                myfit <- ranger(fml, data=aug_df, respect.unordered.factors='order', ...)
                aug_df$tau_hat <- myfit$predictions  # OOB pred
            }

            if (honest) {  # honest EF
                role <- rep("train", nrow(aug_df))
                role[sample(1:nrow(aug_df), nrow(aug_df) / 2, replace = FALSE)] <- "est"
                aug_df$role <- role

                myfit <- ranger(fml, data=aug_df[role=="train", ], respect.unordered.factors='order', ...)

                # get terminal nodes for the honest sample, i.e. run honest sample through forest
                est_leaves <- predict(myfit, aug_df[role=="est", ], type = "terminalNodes")$predictions
                tr_leaves  <- predict(myfit, aug_df[role=="train", ], type = "terminalNodes")$predictions

                # filter out unique leaves for each tree
                unique_leaves <- apply(est_leaves, 2, unique)
                # get honest outcomes (local tau_hat `Y_aug`)
                honest_y <- as.numeric(aug_df[role=="est", ]$Y_aug)
                honest_pred <- GetHonestRcpp(unique_leaves, honest_y, est_leaves, tr_leaves)  # first est then tr

                aug_df$tau_hat <- NA_real_
                aug_df[role=="est", ]$tau_hat <- honest_pred[1 : (length(honest_pred)/2)]
                aug_df[role=="train", ]$tau_hat <- honest_pred[(length(honest_pred)/2+1) : (length(honest_pred))]
            }

            coord_site <- which(aug_df[[site]] == coord_id)
            aug_df <- aug_df[coord_site] # only estimate for coord site
            return(list(myfit=myfit, tau_hat=aug_df$tau_hat, est_leaves=est_leaves, honest_y=honest_y))
        }

        if (is_pred) {  # for prediction
            aug_df$site <- coord_id
            aug_df$site <- factor(aug_df$site)

            if (!honest) {  # adaptive EF
                aug_df$tau_hat <- predict(myfit, data=aug_df)$predictions
            }

            if (honest) {  # honest EF - pred based on est samples
                # get terminal nodes for the honest sample, i.e. run honest sample through forest
                test_leaves   <- predict(myfit, aug_df, type = "terminalNodes")$predictions
                # filter out unique leaves for each tree
                unique_leaves <- apply(est_leaves, 2, unique)
                test_pred <- PredHonestRcpp(unique_leaves, honest_y, est_leaves, test_leaves)
                aug_df$tau_hat <- test_pred
            }

            return(aug_df$tau_hat)
        }
    }


    ## use mean encoding for site index (with grf)
    if (is_encode) {
        covars_all <- c("site_enc", covars)

        if (!is_pred) {  # for estimation
            aug_df[, site_enc := round(mean(Y_aug)), by=site]
            site_vars <- quote(list(site, site_enc))
            site_enc_tab <- unique(aug_df[, eval(site_vars)])

            myfit <- regression_forest(aug_df[, covars_all, with=FALSE], aug_df$Y_aug,
                                       honesty=honest, tune.parameters="all",
                                       num.trees=nrow(aug_df), ...)
            aug_df$tau_hat <- myfit$predictions  # OOB pred

            return(list(myfit=myfit, tau_hat=aug_df$tau_hat,
                        site_enc_tab=site_enc_tab))
        }

        if (is_pred) {  # for prediction
            # mean encoding for coord site
            coord_site_row <- which(site_enc_tab[[site]] == coord_id)
            aug_df$site_enc <- site_enc_tab[coord_site_row]$site_enc

            aug_df$tau_hat <- predict(myfit, aug_df[, covars_all, with=FALSE])$predictions

            return(aug_df$tau_hat)
        }
    }


}
