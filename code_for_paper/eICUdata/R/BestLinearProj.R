# Best linear projection and relevant helper functions modified based on grf 
# https://github.com/grf-labs/grf/blob/master/r-package/grf/R/forest_summary.R
# https://github.com/grf-labs/grf/blob/master/r-package/grf/R/get_scores.R

BestLinearProj <- function(myfit, df,
                           A = NULL,
                           subset = NULL,
                           debiasing.weights = NULL,
                           num.trees.for.weights = 500,
                           vcov.type = "HC3") {
    
    if (!isS4(myfit)) {  # grf or ranger
        clusters <- if (length(myfit$clusters) > 0) {
            myfit$clusters
        } else {
            1:nrow(df)
        }
        observation.weight <- observation_weights(myfit, df)
        
    } else {  # X-learner
        clusters <- 1:nrow(df)
        raw.weights <- rep(1, nrow(df))
        observation.weight <- raw.weights / sum(raw.weights)
    }
    
    subset <- validate_subset(myfit, subset, df)
    subset.clusters <- clusters[subset]
    subset.weights <- observation.weight[subset]
    
    if (length(unique(subset.clusters)) <= 1) {
        stop("The specified subset must contain units from more than one cluster.")
    }
    
    if (!is.null(debiasing.weights)) {
        if (length(debiasing.weights) == nrow(df)) {
            debiasing.weights <- debiasing.weights[subset]
        } else if (length(debiasing.weights) != length(subset)) {
            stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
        }
    }
    
    DR.scores <- get_scores(myfit, df, subset = subset, debiasing.weights = debiasing.weights,
                            num.trees.for.weights = num.trees.for.weights)
    
    if (!is.null(A)) {
        A <- as.matrix(A)
        if (nrow(A) == nrow(df)) {
            A.subset <- A[subset, , drop = FALSE]
        } else if (nrow(A) == length(subset)) {
            A.subset <- A
        } else {
            stop("The number of rows of A does not match the number of training examples.")
        }
        if (is.null(colnames(A.subset))) {
            colnames(A.subset) <- paste0("A", 1:ncol(A))
        }
        DF <- data.frame(target = DR.scores, A.subset)
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


get_scores <- function(myfit, df, subset = NULL,
                       debiasing.weights = NULL,
                       num.trees.for.weights = 500) {
    subset <- validate_subset(myfit, subset, df)
    W.orig <- df$Z[subset]
    W.hat <- df$Z.hat[subset]
    Y.orig <- df$Y[subset]
    Y.hat <- df$Y.hat[subset]
    
    if (!isS4(myfit)) {  # grf or ranger
        tau.hat.pointwise <- predict(myfit, df)$predictions[subset] # ranger for adaptive EF
    } else {  # X-learner
        x_train <- as.data.frame(df[, ..covars])[subset, ]
        tau.hat.pointwise <- EstimateCate(myfit, x_train)
    }
    
    binary.W <- all(df$Z %in% c(0, 1))
    
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
                1:length(df$Y)
            }
            variance_forest <- regression_forest(df$Y,
                                                 (df$Z - df$Z.hat)^2,
                                                 clusters = clusters,
                                                 sample.weights = myfit$sample.weights,
                                                 num.trees = num.trees.for.weights,
                                                 ci.group.size = 1)
            V.hat <- predict(variance_forest, df)$predictions
            debiasing.weights.all <- (df$Z - df$Z.hat) / V.hat
            debiasing.weights <- debiasing.weights.all[subset]
        }
    } else if (length(debiasing.weights) == length(df$Y)) {
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


validate_subset <- function(myfit, subset, df) {
    if (is.null(subset)) {
        subset <- 1:nrow(df)
    }
    if (class(subset) == "logical" && length(subset) == nrow(df)) {
        subset <- which(subset)
    }
    if (!all(subset %in% 1:nrow(df))) {
        stop(paste(
            "If specified, subset must be a vector contained in 1:n,",
            "or a boolean vector of length n."
        ))
    }
    subset
}


observation_weights <- function(myfit, df) {
    # Case 1: No sample.weights
    if (is.null(myfit$sample.weights)) {
        if (length(myfit$clusters) == 0 || !myfit$equalize.cluster.weights) {
            raw.weights <- rep(1, nrow(df))
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
