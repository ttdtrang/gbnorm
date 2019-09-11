#' variance.refs
#'
#' Calculate the average variance of the references in a given expression matrix
#' @param X an expression matrix of the form samples x genes
#' @param idx indices of the reference set
#' @export
variance.refs <- function(X, idx, norm=TRUE) {
    if (length(idx) == 1) {
        return(mean(apply(matrix(X[,idx],ncol=1), MARGIN=2, var)))
    } else {
        return(mean(apply(X[,idx], MARGIN=2, var)))
    }
}

sd.refs <- function(X, idx) {
    if (length(idx) == 1) {
        return(mean(apply(matrix(X[,idx],ncol=1), MARGIN=2, sd)))
    } else {
        return(mean(apply(X[,idx], MARGIN=2, sd)))
    }
}


#' linspace
#' @param min (inclusive) start of the range
#' @param max (inclusive) end of the range
#' @param n number of points, including the two ends
#' @export
linspace <- function(min, max, n) {
    band = (max - min) / (n-1)
    return(c(min, c(min + band*(1:(n-2))), max))
}


flatten_h5dataset <- function(dset) {
    for (cname in names(dset)) {
        if (typeof(dset[, cname]) == 'character') {
            dset[, cname] <- as.character(dset[, cname])
        } else if (typeof(dset[, cname]) == 'integer') {
            dset[, cname] <- as.integer(dset[, cname])
        } else if (typeof(dset[, cname]) == 'double') {
            dset[, cname] <- as.double(dset[, cname])
        }
    }
    return(dset)
}

removeEnsemblVersion <- function(x) {
     return(gsub("(ENS[A-Z]*[GT]\\d+)\\.(\\d+)$", "\\1", x))
}

jaccard.sim <- function(a,b) {
    return(length(intersect(a,b))/ (length(union(a,b))) )
}


partitions <- function(n,k) {
    output = list()
    chunk_size = rep(n %/% k,k)
    if (n%%k != 0) {
        chunk_size[1:(n%%k)] = chunk_size[1:(n%%k)] +1
    }
    curIdx = 1
    for (i in 1:k) {
        output[[i]] <- c(curIdx:(curIdx + chunk_size[i]-1))
        curIdx = curIdx + chunk_size[i]
    }
    return(output)
}

calc.precision <- function(X) {
    X[['precision']] <- X[['TP']] / (X[['TP']] + X[['FP']])
    invisible(X)
}

calc.recall <- function(X) {
    X[['recall']] <- X[['TP']] / (X[['TP']] + X[['FN']])
    invisible(X)
}

calc.accuracy <- function(X) {
    X[['accuracy']] <-
        (X[['TP']] + X[['TN']]) / 
        (X[['TP']] + X[['TN']] + X[['FP']] + X[['FN']])
    invisible(X)
}

calc.f1score <- function(X, precision='precision', recall='recall', weight = c(1,1), name='f1score') {
    p = X[[precision]]
    r = X[[recall]]
    X[[name]] <- (sum(weight) / (weight[1] / p + weight[2] / r))
    # X[['f1score']] <-
    #     (2*X[[precision]] * X[[recall]]) /
    #     (X[[precision]] + X[[recall]])
    invisible(X)        
}

#' Geometric mean
geom.mean <- function(x) {
    log(x) %>%
        mean(na.rm = TRUE) %>%
        exp() %>%
        return()
}
