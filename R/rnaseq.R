#' Build graph from correlation matrix
#' 
#' @import igraph
#' @export
build.graph <- function(corrMatrix, minCor = 0.5) {
    output = corrMatrix
    for (i in 1:nrow(corrMatrix)) {
        for (j in i:ncol(corrMatrix)) {
            if (i == j || is.na(corrMatrix[i,j]) || (corrMatrix[i,j] <= minCor)) {
                output[i,j] <- 0
                output[j,i] <- 0
            }
        }
    } 
    return(igraph::graph.adjacency(output,'undirected',weighted = TRUE,diag = FALSE))
}

#' get.references
#' 
#' @import igraph
#' @export
get.references.v1 <- function(m,tmin=0.7, tmax=0.95, nsteps=5,js.tol=0.1) {
    if (any(m == 0) > 0) { m = m + 1}
    cor.expcnt = cor(log(m,base = 10),method = 'pearson')
    output = list()
    lc.df = data.frame(list(t=tmin + c(0:nsteps)*((tmax-tmin) / nsteps)))
    message(lc.df$t)
    lcliques = list()
    for (i in 1:length(lc.df$t)) {
        message(paste("Working on graph built at t =",lc.df[i,'t']))
        t0 = proc.time()
        g = build.graph(cor.expcnt,minCor = lc.df[i,'t'])
        t1 = proc.time()
        lcliques[[i]] = igraph::largest_cliques(g)[[1]]
        t2 = proc.time()
        if (is.null(output$references)) {
            output$references <- lcliques[[i]]
        }
        if (i > 1) {
            lc.df[i,'js'] = jaccard.sim(lcliques[[i-1]], lcliques[[i]])
            if (lc.df[i,'js'] < js.tol) {
                output$references <- lcliques[[i]]
            }
        } 
        lc.df[i,'cliqueSize']= length(lcliques[[i]])
        lc.df[i,'nVertices']= igraph::vcount(g)
        lc.df[i,'nEdges'] = igraph::ecount(g)
        lc.df[i,'time.graph'] = (t1 - t0)['elapsed']
        lc.df[i,'time.clique'] =(t2 - t1)['elapsed']
    }
    output$runs = lc.df
    return(output)
}

get.references.v2 <- function(m,t=0.7,k=5) {
    if (any(m == 0) > 0) { m = m + 1}
    cor.expcnt = cor(log(m,base = 10),method = 'pearson')
    g = build.graph(cor.expcnt,minCor =t)
    candidates = igraph::maximal.cliques(g)
    errors = rep(Inf,length(candidates))
    message(lc.df$t)
}

#' Identify the best set of references
#' 
#' @export
get.references.blocks <- function(m, py.prog,
                                  cor.method='pearson',
                                  n.runs= 3, min.size = 10,
                                  min.corr = 0.75,
                                  min.count = 0,
                                  log.base = 2) {
    
    ## Reduce the size of graph
    idx.allpresent = which (sapply(1:ncol(m), FUN = function(i) { return(all(m[,i] > min.count)) }))
    m.compressed  = m[,idx.allpresent]
    if (log.base > 0) { m.compressed = log(m.compressed,base = log.base) }
    cor.expcnt = cor(m.compressed,method=cor.method)

    g = build.graph(cor.expcnt, minCor=min.corr)
    
    ## write graph, run block (in python), load results
    gPrefix = paste0('tmp', paste(as.character(sample(9, replace=T)), collapse=''))
    oPrefix = paste0('tmp', paste(as.character(sample(9, replace=T)), collapse=''))
    gFileName  = paste0(gPrefix, '.graphml')
    oLabelFile = paste0(oPrefix, '.tsv')
    oScoreFile = paste0(oPrefix, '.entropy')
    igraph::write.graph(g, file= gFileName, format='graphml')
    
    system(paste('python3', PYPROG, gFileName, oPrefix, n.runs), wait=TRUE)
    
    blocks.entropy = read.table(oScoreFile)
    bestModel = which.min(blocks.entropy[,1]) 
    blocks.all = read.table(oLabelFile,sep='\t')
    blocks = blocks.all[,c(1,2,bestModel+2)]
    names(blocks) = c('CompressedId', 'Name', 'BlockId')
    blocks[,'CompressedId'] = blocks[,'CompressedId'] + 1 # python 0-based to R 1-based index
    blocks[,'BlockId'] = blocks[,'BlockId'] + 1
    
    ## Map the Id in compressed matrix to original matrix
    blocks[,'Id'] = idx.allpresent[blocks[,'CompressedId']]
  
    blocks.df = data.frame(table(blocks$BlockId))
    names(blocks.df) = c('BlockId', 'size')
    blocks.df = blocks.df[blocks.df$size >= min.size,]
    blocks.df[,'Rank1Residuals'] =  rep(0,nrow(blocks.df))
    blocks.df[,'MinCor'] =  rep(0,nrow(blocks.df))
    # blocks.df[,'MedianCor'] =  rep(0,nrow(blocks.df))
    blocks.members= list()
    blocks.memNames= list()
    for (i in 1:nrow(blocks.df)) {
        j = blocks.df$BlockId[i]
        blocks.members[[i]] = blocks[blocks$BlockId == j, 'Id' ]
        blocks.memNames[[i]] = blocks[blocks$BlockId == j, 'Name' ]
        blocks.df[i,'Rank1Residuals'] = rank1.residuals(m[,blocks.members[[i]] ]) 
        idx = blocks[blocks$BlockId == j, 'CompressedId']
        tmpc= cor.expcnt[idx, idx]
        blocks.df[i,'MinCor'] = min(c(tmpc))
        # blocks.df[i,'MedianCor'] = median(a(tmpc))
    }
    remains = which(blocks.df$MinCor > min.corr)
    bId = remains[which.min(blocks.df[remains,'Rank1Residuals'])]
    ## Clean up and return
    system(paste('rm', gFileName, oLabelFile, oScoreFile))
    return(list('Id'= blocks.members[[bId]],
                'Name'= blocks.memNames[[bId]],
                'nVertices' = ncol(m),
                'nVertices.compressed' = ncol(m.compressed)))
}


#' normalize.by.refs
#' 
#' normalize a read count matrix given the set of reference genes identified by id
#' @param X a read-count matrix of the form samples x genes(transcripts)
#' @param ref.idx an integer vector specifying the column indices of X to be used as reference
#' @export
normalize.by.refs <- function(X, ref.idx) {
    if (length(ref.idx) == 1) { # need to be treated specially since dim(X) reduces to NULL and cause error in apply
        Xref = matrix(X[,ref.idx], ncol=1)
    } else {
        Xref = X[,ref.idx]
    }
    normFactors = apply(Xref,MARGIN = 1,FUN = sum)
    # sanity check
    idx.zero = which(normFactors == 0)
    if (length(idx.zero) > 0) {
        message(paste0("All reference transcripts are zero in the sample ", idx.zero, ". Please remove.\n"))
        return(NULL)
    }
    
    X.norm = t(sapply(1:length(normFactors), FUN = function(i) {
        return(X[i,] / normFactors[i])
    }))
    colnames(X.norm) = colnames(X)
    rownames(X.norm) = rownames(X)
    return(X.norm)
}

#' Call calcNormFactors by edgeR
#' @import edgeR
normalize.by.tmm <- function(X,...) {
    normfactors.tmm = edgeR::calcNormFactors(t(X), ...)
    X.normed = t(sapply(1:length(normfactors.tmm), FUN = function(j) {return(X[j,]/normfactors.tmm[j])}))
    rownames(X.normed) = rownames(X)
    return(X)
}

#' generic function for cross-validation
#' 
#' @param X read count matrix
#' @param ref.idx the indices of features to use as references
#' @param diff.func the function that calculate the differences between 2 normalized matrices
#' @param k [=5]
crossval <- function(X,ref.idx,diff.func, k=5) {
    if (k > length(ref.idx)) {
        k <- length(ref.idx)
        warning("crossval: k > number of refererences. k set to length(ref.idx)")
    }
    warning()
    parts = partitions(length(ref.idx),k=k)
    diffs = list() # matrix(rep(0,k*4),ncol=4)
    # colnames(diffs) = c('U', 'S', 'V','delta')
    for (i in 1:k) {
        set1 = ref.idx[parts[[i]]]
        set2 = setdiff(ref.idx,set1)
        m1 = normalize.by.refs(X,set1)
        m2 = normalize.by.refs(X,set2)
        diffs[[i]] = diff.func(m1, m2)
    }
    return(diffs)
}

crossval.refs <- function(X,ref.idx,k=5) {
    parts = partitions(length(ref.idx),k=k)
    diffs = list() # matrix(rep(0,k*4),ncol=4)
    # colnames(diffs) = c('U', 'S', 'V','delta')
    for (i in 1:k) {
        set1 = ref.idx[parts[[i]]]
        set2 = setdiff(ref.idx,set1)
        m1 = normalize.by.refs(X,set1)
        m2 = normalize.by.refs(X,set2)
        m1.svd = svd(m1)
        m2.svd = svd(m2)
        diffs[[i]]$U = m1.svd$u - m2.svd$u
        diffs[[i]]$S = normalize.vec(m1.svd$d) - normalize.vec(m2.svd$d)
        diffs[[i]]$V = m1.svd$v - m2.svd$v
        diffs[[i]]$delta = (m1 / norm(m1,type='F')) - (m2 / norm(m2, type='F'))
    }
    return(diffs)
}

#' the magnitude of the singular values except the first one
#' 
#' 
rank1.residuals <- function(m) {
    m.svd = svd(m)
    return(norm(normalize.vec(m.svd$d)[-1], type='2'))
}

crossval.rank <- function(X,ref.idx,k=5,tol=1e-12) {
    parts = partitions(length(ref.idx),k=k)
    diffs = matrix(rep(0,k*3),ncol=3)
    colnames(diffs) = c('U', 'S', 'V')
    
    for (i in 1:k) {
        set1 = ref.idx[parts[[i]]]
        set2 = setdiff(ref.idx,set1)
        m1 = normalize.by.refs(X,set1)
        m2 = normalize.by.refs(X,set2)
        m.svd = svd(m1 / norm(m1,type='F') - m2 / norm(m2, type='F'))
        # diffs[i,'U'] = norm(m1.svd$u - m2.svd$u,type='F')
        diffs[i,'numericalRank'] = sum(m.svd$d > tol)
        # diffs[i,'V'] = norm(m1.svd$v - m2.svd$v,type='F')
    }
    return(diffs)
}

normalize.vec <- function(x) {
    return(x / sqrt(sum(x^2)))
}

standardize <- function(X) {
    m = matrix(rep(apply(X,2,mean),nrow(X)),nrow=nrow(X), byrow = TRUE)
    s = matrix(rep(apply(X,2,sd),nrow(X)),nrow=nrow(X), byrow = TRUE)
    return((X - m)/(ifelse(s == 0,1,s)))
}

centering <- function(X) {
    centroid = matrix(rep(colSums(X) / nrow(X),nrow(X)),nrow=nrow(X),byrow = TRUE) 
    return(X - centroid)
}

rms <- function(X) {
    return(sqrt(mean(c(X^2))))
}

#' Standardized root mean squared deviation between two matrices
#' 
#' @export
srms <- function(X1, X2) {
    return(rms(standardize(X1) - standardize(X2)))
}
#' Least root mean square distance between two set of points in d-dimension
#' 
#' Find and perform the transformation (translation + rotation) on Q that minimizes the distance between 2 matrices P and Q.
#' Return the root mean squared of the differences
#' 
#' 
#' @export
rmsd <- function (P, Q, transform=FALSE) {
    # sanity check
    if (any(dim(P) != dim(Q))) {
        stop("Incompatible comparison. P and Q must have same dimensions.")
    }
    d = ncol(P)
    Pc = centering(P)
    Qc = centering(Q)
    # This scaling is not rigorously tested
    Ps = Pc / norm(Pc,'F')
    Qs = Qc / norm(Qc,'F')
   
    if (transform) {
   
        A = t(Ps) %*% Qs
        A.svd = svd(A)
        S = diag(1,nrow=nrow(A.svd$u))
        S[d,d] = det(A.svd$u %*% t(A.svd$v))
        Rot = A.svd$u %*% S %*% t(A.svd$v)
        diff = Ps%*% Rot - Qs
    } else {
        diff = Pc - Qc
    }
    return(rms(diff))
}
