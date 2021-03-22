#' Ensemble Tree (ET)
#'
#' Build and/or predict on an ensemble regression tree.
#'
#' @param coord_id   Site index for coordinating site.
#' @param aug_df     The augmented data frame used to fit an ensemble tree (`data.table`).
#' @param site       Variable name for site indicator.
#' @param covars     A vector of covariate names used.
#' @param honest     Whether to use honest splitting (i.e., subsample splitting). Default is FALSE.
#' @param is_pred    Whether to build an ensemble tree or make prediction. Default is FALSE.
#' @param myfit      A fitted ensemble tree (for prediction purpose). Default is NULL.
#' @param leaves_tab A table indicates leaf assignments in a fitted ensemble tree
#'                   (for prediction purpose). Default is NULL.
#' @param \\dots     Additional arguments for building the tree.
#'
#' @return Training: return a fitted ensemble tree and estimation of the input data;
#'         Prediction: return predictions of the input data.
#' @export
#' @import stats
#' @import data.table
#' @import rpart
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
#' res_et <- EnsemTree(coord_id, aug_df, "site", covars)
#' et_hat <- EnsemTree(coord_id, coord_test, "site", covars, is_pred=TRUE,
#'           myfit=res_et$myfit, leaves_tab=res_et$leaves_tab)
#'
#' res_et <- EnsemTree(coord_id, aug_df, "site", covars, honest=TRUE)
#' et_hat <- EnsemTree(coord_id, coord_test, "site", covars, honest=TRUE, is_pred=TRUE,
#'           myfit=res_et$myfit, leaves_tab=res_et$leaves_tab)
#'
#'
EnsemTree <- function(coord_id, aug_df, site, covars, honest=FALSE, is_pred=FALSE,
                      myfit=NULL, leaves_tab=NULL, ...) {

    # checking
    CheckVar(coord_id, "coord_id")
    aug_df <- CheckDF(aug_df, "aug_df")
    CheckVar(site, "site")
    CheckVector(covars, "covars")
    if (is_pred==T & honest==T) {
        CheckMyfit(myfit)
        CheckNull(leaves_tab, "leaves_tab")
    }


    # bind variables locally to the function (define as NULL)
    Y_aug <- leaff <- tau_hat<- NULL

    if (!is_pred) {  # for estimation
        fml <- as.formula(paste("Y_aug ~ ", site, " + ", paste(covars, collapse="+")))

        if (!honest) {  # adaptive ET
            myfit <- rpart(fml, data=aug_df, ...) #, xval=5, cp=0
            myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])
            aug_df$tau_hat <- predict(myfit)
        }

        if (honest) {  # honest ET
            role <- rep("train", nrow(aug_df))
            role[sample(1:nrow(aug_df), nrow(aug_df) / 2, replace = FALSE)] <- "est"
            aug_df$role <- role

            ## train - build tree
            myfit <- rpart(fml, data=aug_df[role=="train", ], ...) #, xval=5, cp=0
            myfit <- prune(myfit, cp=myfit$cptable[which.min(myfit$cptable[,"xerror"]),"CP"])

            ## honest estimate by df_est
            ## tree for predict which leave the obs belongs to
            # replace predicted y values in the model frame with the corresponding node numbers
            et_node = myfit
            et_node$frame$yval = as.numeric(rownames(et_node$frame))
            aug_df[role=="est", leaff := predict(et_node, newdata=aug_df[role=="est", ], type='vector')]
            aug_df[role=="est", tau_hat := mean(Y_aug, na.rm=TRUE), by = leaff]

            leaves_tab <- unique(aug_df[role=="est", c("leaff", "tau_hat")])
            aug_df[role=="train", leaff := predict(et_node, newdata=aug_df[role=="train", ], type='vector')]
            # assign honest estimates to train set
            for (i in 1:length(leaves_tab$leaff)) {
                aug_df[role=="train" & leaff==leaves_tab$leaff[i], tau_hat := leaves_tab$tau_hat[i]]
            }
        }

        coord_site <- which(aug_df[[site]] == coord_id)
        aug_df <- aug_df[coord_site] # only estimate for coord site
        return(list(myfit=myfit, tau_hat=aug_df$tau_hat, leaves_tab=leaves_tab))
    }

    if (is_pred) {  # for prediction
        aug_df$site <- coord_id
        aug_df$site <- factor(aug_df$site)

        if (!honest) {
            aug_df$tau_hat <- predict(myfit, aug_df)
        }

        if (honest) {
            et_node = myfit
            et_node$frame$yval = as.numeric(rownames(et_node$frame))
            aug_df[, leaff := predict(et_node, newdata=aug_df, type='vector')]  # predict which leaf falls into
            for (i in 1:length(leaves_tab$leaff)) {  # assign honest estimates to test set
                aug_df[leaff==leaves_tab$leaff[i], tau_hat := leaves_tab$tau_hat[i]]
            }
        }

        return(aug_df$tau_hat)
    }


}
