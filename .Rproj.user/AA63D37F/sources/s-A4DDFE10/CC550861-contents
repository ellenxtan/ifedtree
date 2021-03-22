
CheckMyfit <- function(myfit) {
    if (!("grf" %in% class(myfit) |
          "ranger" %in% class(myfit) |
          "rpart" %in% class(myfit))) {
        stop("myfit should be either a ranger/grf/rpart object")
    }
}


CheckDF <- function(df, df_name) {
    if (is.matrix(df) || is.data.frame(df)) {
        df <- as.data.table(df)
    } else {
        stop(paste(df_name,
                   "should be either a data.table/data.frame/matrix"))
    }
    return(df)
}


CheckVar <- function(var, var_name) {
    if (length(var) != 1) {
        stop(paste(var_name,
                   "should be of length 1"))
    }
}


CheckCharacter <- function(var, var_name) {
    if (!is.character(var)) {
        stop(paste(var_name,
                   "should be a string of class character"))
    }
    if (length(var) != 1) {
        stop(paste(var_name,
                   "should be of length 1"))
    }
}


CheckVector <- function(vec, vec_name) {
    if (!is.vector(vec) | !is.character(vec)) {
        stop(paste(vec_name,
            "should be a vector of character"))
    }
}


CheckNull <- function(obj, obj_name) {
    if (is.null(obj)) {
        stop(paste(obj_name,
                   "should not be null"))
    }
}


CheckPDF <- function(xpdf) {
    if (!xpdf %in% c("unif", "norm", "binary")) {
        stop("Distribution of X should be either 'unif'/'norm'/'binary' ")
    }
}


CheckNumeric <- function(num, num_name) {
    if (!is.numeric(num)) {
        stop(paste(num_name,
                   "should be a numeric number"))
    }
    if (length(num) != 1) {
        stop(paste(num_name,
                   "should be of length 1"))
    }
}


CheckList <- function(lst, lst_name) {
    if (!is.list(lst)) {
        stop(paste(lst_name, "should be a list"))
    }
}
