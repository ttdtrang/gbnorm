#' Condition-number-based deviation of X1 and X2
#'
#' Much faster with thin matrix (nrow >> ncol).
#'
#' @import MASS
#' @param X m x n matrix
#' @param Y m x n matrix
#' @param always.small [=TRUE] if set to TRUE, cdev will be calculate for the smaller version of transformation matrix
#' @export
cdev <- function(X, Y, always.small = TRUE) {
    if (!all(dim(X) == dim(Y))) stop("cdev is only defined for 2 matrices with the same shape.")
    if (always.small & (dim(X)[2] > dim(X)[1])) {
        return(fast.condNumber(MASS::ginv(t(X)) %*% t(Y)))
    } else {
        return(fast.condNumber((MASS::ginv(X) %*% Y)))
    }
}


#' Condition number of a matrix X
condNumber <- function(X) {
    sigmas = svd(X)[['d']]
    return(sigmas[1] / sigmas[length(sigmas)])
}

condNumber.norm <- function(X) {
    return(norm(X,'2') %*% norm(MASS::ginv(X), '2'))
}



#' Condition number using fast.svd
#'
#' When X is rectangular, the SVD of XX' (or X'X) will be faster if ncol > nrow (or nrow > ncol)
#' Since we only care about the ratio between first and last singular values,
#' there's no need to calculate V'
fast.condNumber <- function(X) {
    sigmas = fast.sv(X)
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
