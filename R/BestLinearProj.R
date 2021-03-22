#' Best linear projection of conditional average treatment effect estimates
#'
#' @param myfit       A fitted ensemble tree or forest.
#' @param coord_df    Data of coordinating site (`data.table`).
#' @param coord_id    Site index for coordinating site.
#' @param site        Variable name for site indicator.
#' @param treat       Variable name for binary treatment.
#' @param outcome     Variable name for the outcome.
#' @param covars      A vector of covariate names used.
#' @param subset      Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param site_enc_tab  A data.table of mean outcome for each site. Default is NULL.
#'                      If class of myfit is a `grf`, "site_enc_tab" should not be NULL.
#' @param debiasing.weights  (ignore, only equal weights are supported in current version)
#'               A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param num.trees.for.weights (ignore, only equal weights are supported in current version)
#'               In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when `debiasing.weights` = NULL.
#' @param vcov.type   Optional covariance type for standard errors. The possible
#'  options are HC0, ..., HC3. The default is "HC3", which is recommended in small
#'  samples and corresponds to the "shortcut formula" for the jackknife
#'  (see MacKinnon & White for more discussion, and Cameron & Miller for a review).
#'  For large data sets with clusters, "HC0" or "HC1" are significantly faster to compute.
#'
#' @return  An estimate of the best linear projection, along with coefficient standard errors.
#' @export
#'
#' @note    The function is modified based on the R package `grf`
#' \url{https://github.com/grf-labs/grf/blob/master/r-package/grf/R/forest_summary.R}.
#'
#'@references Cameron, A. Colin, and Douglas L. Miller. "A practitioner's guide to
#'  cluster-robust inference." Journal of human resources 50, no. 2 (2015): 317-372.
#' @references MacKinnon, James G., and Halbert White. "Some heteroskedasticity-consistent
#'  covariance matrix estimators with improved finite sample properties."
#'  Journal of Econometrics 29.3 (1985): 305-325.
#' @references Semenova, Vira, and Victor Chernozhukov. "Debiased Machine Learning of
#'  Conditional Average Treatment Effects and Other Causal Functions".
#'  The Econometrics Journal (2020).
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
#' ## ensemble forest
#' # Treat each site as a distinct factor
#' myfit <- EnsemForest(coord_id, aug_df, "site", covars)$myfit
#' BestLinearProj(myfit, coord_df, coord_id, "site", "Z", "Y", covars)
#'
#' # Mean encoding as surrogate for site index
#' res_ef <- EnsemForest(coord_id, aug_df, "site", covars, is_encode=TRUE)
#' myfit <- res_ef$myfit
#' site_enc_tab <- res_ef$site_enc_tab
#' BestLinearProj(myfit, coord_df, coord_id, "site", "Z", "Y", covars, site_enc_tab=site_enc_tab)
#'
#' ## ensemble tree
#' myfit <- EnsemTree(coord_id, aug_df, "site", covars)$myfit
#' BestLinearProj(myfit, coord_df, coord_id, "site", "Z", "Y", covars)
#'
#'
BestLinearProj <- function(myfit, coord_df, coord_id, site, treat, outcome, covars,
                           subset=NULL, site_enc_tab=NULL,
                           debiasing.weights=NULL,
                           num.trees.for.weights=500, vcov.type="HC3") {

    # checking
    CheckMyfit(myfit)
    coord_df <- CheckDF(coord_df, "coord_df")
    CheckVar(coord_id, "coord_id")
    CheckVar(site, "site")
    CheckVar(treat, "treat")
    CheckVar(outcome, "outcome")
    CheckVector(covars, "covars")
    if ("grf" %in% class(myfit)) {
        CheckNull(site_enc_tab, "site_enc_tab")
    }


    X <- as.matrix(coord_df[, covars, with=FALSE])
    W <- coord_df[[treat]]
    Y <- coord_df[[outcome]]
    coord_df$Z.hat <- regression_forest(X, W)$predictions
    coord_df$Y.hat <- regression_forest(X, Y)$predictions

    coord_df[[site]] <- coord_id
    coord_df[[site]] <- factor(coord_df[[site]])
    covar_mat <- as.matrix(coord_df[, covars, with=FALSE])

    if (!isS4(myfit)) {  # grf or ranger
        clusters <- if (length(myfit$clusters) > 0) {
            myfit$clusters
        } else {
            1:nrow(coord_df)
        }
        observation.weight <- ObservationWeights(myfit, coord_df)

    } else {  # X-learner
        clusters <- 1:nrow(coord_df)
        raw.weights <- rep(1, nrow(coord_df))
        observation.weight <- raw.weights / sum(raw.weights)
    }

    subset <- ValidateSubset(myfit, subset, coord_df)
    subset.clusters <- clusters[subset]
    subset.weights <- observation.weight[subset]

    if (length(unique(subset.clusters)) <= 1) {
        stop("The specified subset must contain units from more than one cluster.")
    }

    if (!is.null(debiasing.weights)) {
        if (length(debiasing.weights) == nrow(coord_df)) {
            debiasing.weights <- debiasing.weights[subset]
        } else if (length(debiasing.weights) != length(subset)) {
            stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
        }
    }

    if ("grf" %in% class(myfit)) {
        coord_site_row <- which(site_enc_tab[[site]] == coord_id)
        coord_df$site_enc <- site_enc_tab[coord_site_row]$site_enc
    }

    DR.scores <- GetScores(myfit, coord_df, covars, subset = subset,
                           debiasing.weights = debiasing.weights,
                            num.trees.for.weights = num.trees.for.weights)

    if (!is.null(covar_mat)) {
        covar_mat <- as.matrix(covar_mat)
        if (nrow(covar_mat) == nrow(coord_df)) {
            covar.subset <- covar_mat[subset, , drop = FALSE]
        } else if (nrow(covar_mat) == length(subset)) {
            covar.subset <- covar_mat
        } else {
            stop("The number of rows of covariate matrix does not match the number of training examples.")
        }
        if (is.null(colnames(covar.subset))) {
            colnames(covar.subset) <- paste0("C", 1:ncol(covar_mat))
        }
        DF <- data.frame(target = DR.scores, covar.subset)
    } else {
        DF <- data.frame(target = DR.scores)
    }

    blp.ols <- lm(target ~ ., weights = subset.weights, data = DF)
    blp.summary <- lmtest::coeftest(blp.ols,
                                    vcov = sandwich::vcovCL,
                                    type = vcov.type,
                                    cluster = subset.clusters
    )
    attr(blp.summary, "method") <-
        paste0("Best linear projection of the conditional average treatment effect.\n",
               "Confidence intervals are cluster- and heteroskedasticity-robust ",
               "(", vcov.type, ")")

    blp.summary
}


#' Compute doubly robust scores for a fitted ensemble tree or forest.
#'
#' Compute doubly robust scores for a fitted ensemble tree or forest.
#' Only valid to use when estimating conditional average treatment effects.
#'
#' @param myfit     A fitted ensemble tree or forest.
#' @param coord_df  Data of coordinating site (`data.table`).
#' @param covars    A vector of covariate names used.
#' @param subset    Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights     A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when `debiasing.weights` = NULL.
#'
#' @export
#'
#' @note    The function is modified based on the R package `grf`
#' \url{https://github.com/grf-labs/grf/blob/master/r-package/grf/R/get_scores.R}.
#'
#'
GetScores <- function(myfit, coord_df, covars, subset=NULL,
                       debiasing.weights=NULL, num.trees.for.weights=500) {
    subset <- ValidateSubset(myfit, subset, coord_df)
    W.orig <- coord_df$Z[subset]
    W.hat <- coord_df$Z.hat[subset]
    Y.orig <- coord_df$Y[subset]
    Y.hat <- coord_df$Y.hat[subset]

    if (!isS4(myfit)) {  # grf or ranger
        if ("ranger" %in% class(myfit)) {
            tau.hat.pointwise <- predict(myfit, coord_df)$predictions[subset]
        }
        if ("grf" %in% class(myfit)) {
            covars_all <- c(covars, "site_enc")
            tau.hat.pointwise <- predict(myfit, coord_df[, covars_all, with=FALSE])$predictions[subset]
        }
        if ("rpart" %in% class(myfit)) {
            tau.hat.pointwise <- predict(myfit, coord_df)[subset]
        }

    } else {  # X-learner
        x_train <- as.data.frame(coord_df[, covars, with=FALSE])[subset, ]
        tau.hat.pointwise <- causalToolbox::EstimateCate(myfit, x_train)
    }

    binary.W <- all(coord_df$Z %in% c(0, 1))

    if (is.null(debiasing.weights)) {
        if (binary.W) {
            debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
        } else {
            # Start by learning debiasing weights if needed.
            # The goal is to estimate the variance of W given X. For binary treatments,
            # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
            # so this step is not needed. Note that if we use the present CAPE estimator
            # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
            # exactly the AIPW estimator of the CATE.
            clusters <- if (length(myfit$clusters) > 0) {
                myfit$clusters
            } else {
                1:length(coord_df$Y)
            }
            variance_forest <- grf::regression_forest(coord_df$Y,
                                                 (coord_df$Z - coord_df$Z.hat)^2,
                                                 clusters = clusters,
                                                 sample.weights = myfit$sample.weights,
                                                 num.trees = num.trees.for.weights,
                                                 ci.group.size = 1)
            V.hat <- predict(variance_forest, coord_df)$predictions
            debiasing.weights.all <- (coord_df$Z - coord_df$Z.hat) / V.hat
            debiasing.weights <- debiasing.weights.all[subset]
        }
    } else if (length(debiasing.weights) == length(coord_df$Y)) {
        debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset))  {
        stop("If specified, debiasing.weights must have length n or |subset|.")
    }

    # Form AIPW scores. Note: We are implicitly using the following
    # estimates for the regression surfaces E[Y|X, W=0/1]:
    # Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
    # Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise
    Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

    tau.hat.pointwise + debiasing.weights * Y.residual
}


## The function is modified based on the R package `grf`
## \url{https://github.com/grf-labs/grf/blob/master/r-package/grf/R/input_utilities.R}.
ValidateSubset <- function(myfit, subset, coord_df) {
    if (is.null(subset)) {
        subset <- 1:nrow(coord_df)
    }
    if (class(subset) == "logical" && length(subset) == nrow(coord_df)) {
        subset <- which(subset)
    }
    if (!all(subset %in% 1:nrow(coord_df))) {
        stop(paste(
            "If specified, subset must be a vector contained in 1:n,",
            "or a boolean vector of length n."
        ))
    }
    subset
}


## The function is modified based on the R package `grf`
## \url{https://github.com/grf-labs/grf/blob/master/r-package/grf/R/input_utilities.R}.
ObservationWeights <- function(myfit, coord_df) {
    # Case 1: No sample.weights
    if (is.null(myfit$sample.weights)) {
        if (length(myfit$clusters) == 0 || !myfit$equalize.cluster.weights) {
            raw.weights <- rep(1, nrow(coord_df))
        } else {
            # If clustering with no sample.weights provided and equalize.cluster.weights = TRUE, then
            # give each observation weight 1/cluster size, so that the total weight of each cluster is the same.
            clust.factor <- factor(myfit$clusters)
            inverse.counts <- 1 / as.numeric(Matrix::colSums(Matrix::sparse.model.matrix(~ clust.factor + 0)))
            raw.weights <- inverse.counts[as.numeric(clust.factor)]
        }
    }

    # Case 2: sample.weights provided
    if (!is.null(myfit$sample.weights)) {
        if (length(myfit$clusters) == 0 || !myfit$equalize.cluster.weights) {
            raw.weights <- myfit$sample.weights
        } else {
            stop("Specifying non-null sample.weights is not allowed when equalize.cluster.weights = TRUE")
        }
    }

    return (raw.weights / sum(raw.weights))
}
