#' Condition-number-based deviation of X and Y,
#' for generally given X and Y.
#'
#' Much faster with thin matrix (nrow >> ncol).
#'
#' @import MASS
#' @param X m x n matrix
#' @param Y m x n matrix
#' @param always.small [=TRUE] if set to TRUE, cdev will be calculate for the smaller version of transformation matrix
#' @export
cdev.generic <- function(X, Y, always.small = TRUE) {
    if (!all(dim(X) == dim(Y))) stop("cdev is only defined for 2 matrices with the same shape.")
    if (always.small & (dim(X)[2] > dim(X)[1])) {
        return(fast.condNumber(MASS::ginv(t(X)) %*% t(Y)))
    } else {
        return(fast.condNumber((MASS::ginv(X) %*% Y)))
    }
}

#' Condition-number-based deviation of X and Y,
#' assuming X can be transformed into Y by scaling each columns with
#' a coefficient.
#'
#' @param X a matrix
#' @param Y another matrix with the same shape as X
#' @param generic whether to use the generic definition of cdev, without assuming only column-scaling is involved between X and Y
#' @export
cdev <- function(X, Y, scale = FALSE, generic = FALSE, ...) {
    if (generic) return(cdev.generic(X,Y,...))
    # Pick one of the genes with minimum counts > 0,
    # reverse-engineer the scaling factors
    an_eligible_feature = which(apply(X, MARGIN = 1, min) > 0)[1]
    scalingFactors = X[an_eligible_feature,] / Y[an_eligible_feature,]
    return(cdev.from.scaling.factors(scalingFactors))
}

#' @export
cdev.from.scaling.factors <- function(scalingFactors) {
    return(max(scalingFactors) / min(scalingFactors))
}



#' Condition number using fast.svd
#'
#' When X is rectangular, the SVD of XX' (or X'X) will be faster if ncol > nrow (or nrow > ncol)
#' Since we only care about the ratio between first and last singular values,
#' there's no need to calculate V'
fast.condNumber <- function(X) {
    sigmas = svd(X, nu = 0, nv = 0)[['d']]
    return(sigmas[1] / sigmas[length(sigmas)])
}

#' the magnitude of the singular values except the first one
#' 
#' @export
rank1.residuals <- function(m) {
    m.sv = sv(m)
    return(norm(normalize.vec(m.sv)[-1], type='2'))
}

fast.sv <- function(X) {
     if (which.min(dim(X)) == 1) {
        xx = X %*% t(X)
    } else {
        xx = t(X) %*% X
    }
    return(sqrt(svd(xx, nu = 0, nv = 0)[['d']]))
}

sv <- function(X) {
    return(svd(X, nu = 0, nv = 0)[['d']])
}

normalize.vec <- function(x) {
    return(x / sqrt(sum(x^2)))
}

# #' Coefficient of variation
# #' 
# function(x) {
#     return(sd(x)/mean(x))
# }

#' Coefficient of variation, unbiased estimate for log-normal variable
#' 
cv.lognorm <- function(x, pseudocount = 1) {
    return(sqrt(exp(sd(log(x+pseudocount))^2 ) - 1))
}

