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
get.references.blocks <- function(m,
                                  cor.method='pearson',
                                  n.runs= 3, min.size = 10,
                                  min.corr = 0.75,
                                  min.count = 0,
                                  log.base = 2,
                                  out.file = '',
                                  debug = FALSE) {
 
    if (system.file("python/sbm.py", package = "gbnorm") == "") {
        stop("Python script not found. Cannot run sbm.py")
    }
    if (is.null(rownames(m))) rownames(m) <- 1:nrow(m)
    if (is.null(colnames(m))) colnames(m) <- 1:ncol(m)

    ## Make sure sbm is executable
    sbm.cmd = paste("python3", system.file("python/sbm.py", package = "gbnorm"))
    tryCatch(
        system(sbm.cmd,ignore.stderr = TRUE, ignore.stdout = TRUE),
        error = function(e) {
            stop(e)
        })
    
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
   
    sbm.cmd <- paste(sbm.cmd, gFileName, oPrefix, n.runs)
    message("Running blocks with the command")
    message(sbm.cmd)
    system(sbm.cmd, wait=TRUE)
    
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
    output <- list(
        'Id'= NULL,
        'Name'= NULL,
        'nVertices' = ncol(m),
        'nVertices.compressed' = ncol(m.compressed)
    )
    remains = which(blocks.df$MinCor > min.corr)
    if (length(remains) == 0) {
        message("No block satisfies mininum requirements.")
        if (debug) {
            message(blocks.df)
        }
        return(output)
    }
    bId = remains[which.min(blocks.df[remains,'Rank1Residuals'])]
    ## Clean up and return
    if (!debug) system(paste('rm', gFileName, oLabelFile, oScoreFile))
    if (out.file == '') out.file <- paste0(oPrefix, '.RDS')
    output[['Id']] <- blocks.members[[bId]]
    output[['Name']] <- blocks.memNames[[bId]]
    saveRDS(output, file = out.file)
    return(output)
}

#' get.references.apcluster
#' 
#' Graph-based identify the set of references using affinity propagation as the underlying clustering method
#' @importFrom apcluster apcluster
#' @export
get.references.apcluster <- function(m, 
                                     cor.method='pearson',
                                     min.size = 10,
                                     min.corr = 0.75,
                                     min.count = 0,
                                     med.quantile = 0.5,
                                     log.base = 2,
                                     debug = FALSE,...) {
    ## Reduce the size of graph
    isUniversal = sapply(1:ncol(m), FUN = function(i) { return(all(m[,i] > min.count)) })
    medium.threshold = sapply(1:nrow(m), function(i) {
        (m[i,] > 0) %>%
            which() %>%
            `[`(m, i, .) %>%
            quantile(probs = med.quantile) %>%
            return()
    })
    isHighEnough = sapply(1:ncol(m), function(i) {
        ((m[,i] - medium.threshold) > 0) %>%
            all() %>% return()
    })
    
    idx.candidates = (isUniversal & isHighEnough) %>% which()
    m.compressed  = m[,idx.candidates]
    if (log.base > 0) { m.compressed = log(m.compressed,base = log.base) }
    if (debug) {message(sprintf("Building compressed graph with %d vertices", length(idx.candidates)))}
    startT = proc.time()
    cor.expcnt = cor(m.compressed,method=cor.method)
    dt1 = proc.time() - startT
    cor.expcnt[cor.expcnt < min.corr] <- 0
    if (debug) {message(sprintf("Calculate correlation in %f sec.", dt1['elapsed']))}
    # Affinity propagation
    startT = proc.time()
    apclust = apcluster::apcluster(s= cor.expcnt, ...)
    dt2 = proc.time() - startT
    if (debug) {message(sprintf("Run AP clustering in %f sec.", dt2['elapsed']))}
    ## Mapping compressed id to original id
    cl.assignment = data.frame('Id' = 1:ncol(m),
                               'Name' = colnames(m),
                               'ClusterId' = rep(0, ncol(m)),
                               stringsAsFactors = FALSE)
    clSizes = sapply(apclust@clusters, length)
    for (cl in 1:length(apclust@clusters)) {
        cl.assignment[idx.candidates,][apclust@clusters[[cl]], 'ClusterId'] = cl
    }
    
    cl.members = list()
    cl.memNames = list()
    cl.df = data.frame('ClusterId' = 1:length(apclust@clusters),
                       'size' = clSizes)
    for (i in 1:nrow(cl.df)) {
        cl.members[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Id' ]
        cl.memNames[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Name' ]
        if (length(cl.members[[i]]) >= 2) {
            cl.df[i,'Rank1Residuals'] = rank1.residuals(m[,cl.members[[i]] ])
        } else {
            cl.df[i,'Rank1Residuals'] = NA
        }
        # cl.df[i,'size'] =  length(cl.members[[i]])
    }
    
    ## largest cluster with the smallest rank1 residuals
    # largest = which(cl.df$size == max(cl.df$size))
    # cId = which.min(cl.df[largest, 'Rank1Residuals'])
    ## best cluster scored by both rank1residuals and size, but weighted slighly more by rank1residuals
    cId = best.cluster(cl.df, features = c('Rank1Residuals', 'size'), weights = c(-0.7, 0.3))
    # if (debug) {print(cl.df)}
    return(list('Id'= cl.members[[cId]],
                'Name'= cl.memNames[[cId]],
                'nVertices' = ncol(m),
                'nVertices.compressed' = ncol(m.compressed)))
}

#' get.references.dbscan
#' 
#' Graph-based identify the set of references using DBSCAN as the underlying clustering method
#' @importFrom dbscan dbscan
#' @export
get.references.dbscan <- function(m, 
                                cor.method='pearson',
                                min.size = 3,
                                min.count = 0,
                                med.quantile = 0.5,
                                log.base = 2,
                                eps = 0.1,
                                debug = FALSE,...) {
    ## Reduce the size of graph
    isUniversal = sapply(1:ncol(m), FUN = function(i) { return(all(m[,i] > min.count)) })
    medium.threshold = sapply(1:nrow(m), function(i) {
        (m[i,] > 0) %>%
            which() %>%
            `[`(m, i, .) %>%
            quantile(probs = med.quantile) %>%
            return()
    })
    isHighEnough = sapply(1:ncol(m), function(i) {
        ((m[,i] - medium.threshold) > 0) %>%
            all() %>% return()
    })
    
    idx.candidates = (isUniversal & isHighEnough) %>% which()
    m.compressed  = m[,idx.candidates]
    if (log.base > 0) { m.compressed = log(m.compressed,base = log.base) }
    if (debug) {message(sprintf("Building compressed graph with %d vertices", length(idx.candidates)))}
    startT = proc.time()
    cor.expcnt = cor(m.compressed,method=cor.method, use = 'pairwise.complete')
    dt1 = proc.time() - startT
    if (debug) {message(sprintf("Calculate correlation in %f sec.", dt1['elapsed']))}
    
    ## DBSCAN requires distance matrix as input
    startT = proc.time()
    d.compressed = cor.expcnt %>%
        `[<-`(cor.expcnt< 0, 0.001) %>%
        log10() %>%
        `-`(0, .)
    ### scaling the distance from 0 to 1
    d.compressed = (d.compressed - min(d.compressed)) / (max(d.compressed) - min(d.compressed))
    clust = dbscan::dbscan(as.dist(d.compressed), eps = eps, minPts = min.size)
    dt2 = proc.time() - startT
   
    if (debug) {print(str(d.compressed))} 
    if (debug) {message(sprintf("Run DBSCAN in %f sec.", dt2['elapsed']))}
    if (debug) {print(clust)}
    ## Mapping compressed id to original id
    cl.assignment = data.frame('Id' = 1:ncol(m),
                               'Name' = colnames(m),
                               'ClusterId' = rep(0, ncol(m)))
    cl.assignment[idx.candidates, 'ClusterId'] = clust$cluster
    cl.assignment[idx.candidates, 'ClusterId'] = clust$cluster
    
    cl.members = list()
    cl.memNames = list()
    cl.df = data.frame('ClusterId' = c(1:length(table( clust$cluster[clust$cluster !=0]))))
    for (i in 1:nrow(cl.df)) {
        cl.members[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Id' ]
        cl.memNames[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Name' ]
        if (length(cl.members[[i]]) >= min.size) {
            cl.df[i,'Rank1Residuals'] = rank1.residuals(m[,cl.members[[i]] ])
        } else {
            cl.df[i,'Rank1Residuals'] = NA
        }
        cl.df[i,'size'] =  length(cl.members[[i]])
    }
    
    ## largest cluster with the smallest rank1 residuals
    # largest = which(cl.df$size == max(cl.df$size))
    # cId = which.min(cl.df[largest, 'Rank1Residuals'])
    ## best cluster scored by both rank1residuals and size
    if (debug) {print(cl.df)}
    cId = best.cluster(cl.df, features = c('Rank1Residuals', 'size'), weights = c(-0.5, 0.5))
    if (debug) {message(sprintf("selected cluster: %s", cId))}
    return(list('Id'= cl.members[[cId]],
                'Name'= cl.memNames[[cId]],
                'nVertices' = ncol(m),
                'nVertices.compressed' = ncol(m.compressed)))
    
}

#' get.references.hclust
#' 
#' @import dynamicTreeCut
#' @export
get.references.hclust <- function(m, 
                                  cor.method='pearson',
                                  min.corr = 0.75,
                                  min.size = 3,
                                  min.count = 0,
                                  med.quantile = 0.5,
                                  log.base = 2,
                                  debug = FALSE) {
    ## Reduce the size of graph
    isUniversal = sapply(1:ncol(m), FUN = function(i) { return(all(m[,i] > min.count)) })
    medium.threshold = sapply(1:nrow(m), function(i) {
        (m[i,] > 0) %>%
            which() %>%
            `[`(m, i, .) %>%
            quantile(probs = med.quantile) %>%
            return()
    })
    isHighEnough = sapply(1:ncol(m), function(i) {
        ((m[,i] - medium.threshold) > 0) %>%
            all() %>% return()
    })
    
    idx.candidates = (isUniversal & isHighEnough) %>% which()
    m.compressed  = m[,idx.candidates]
    if (log.base > 0) { m.compressed = log(m.compressed,base = log.base) }
    if (debug) {message(sprintf("Building compressed graph with %d vertices", length(idx.candidates)))}
    startT = proc.time()
    corM = cor(m.compressed,method=cor.method)
    dt1 = proc.time() - startT
    corM[corM < min.corr ] <- 0
    if (debug) {message(sprintf("Calculate correlation in %f sec.", dt1['elapsed']))}
    
    ## hclust requires distance matrix as input
    startT = proc.time()
    d.compressed = corM %>%
        `[<-`(corM < min.corr, 0.001) %>%
        log10() %>%
        `-`(0, .)
    ### scaling the distance from 0 to 1
    d.compressed = (d.compressed - min(d.compressed)) / (max(d.compressed) - min(d.compressed))
    clust = hclust(as.dist(d.compressed), method = 'average') %>% dynamicTreeCut::cutreeDynamic(distM = d.compressed)
    # cluster id start at 0, make it start at 1 so 0 can be used for non-clustered points
    clust = clust +1
    dt2 = proc.time() - startT
    if (debug) {message(sprintf("Run hierarchical clustering in %f sec.", dt2['elapsed']))}
    if (debug) { message(str(clust)) } 
    
    ## Mapping compressed id to original id
    cl.assignment = data.frame('Id' = 1:ncol(m), 'Name' = colnames(m), 'ClusterId' = rep(0, ncol(m)))
    cl.assignment[idx.candidates, 'ClusterId'] = clust
   
    cl.members = list()
    cl.memNames = list()
    cl.df = data.frame('ClusterId' = as.numeric(names(table(clust))))
    for (i in 1:nrow(cl.df)) {
        cl.members[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Id' ]
        cl.memNames[[i]] = cl.assignment[cl.assignment$ClusterId == i, 'Name' ]
        cl.df[i,'nReferences'] = length(grep('ERCC',x = cl.memNames[[i]]))
        cl.df[i,'size'] =  length(cl.members[[i]])
        if (cl.df[i,'size'] < min.size) {
            cl.df[i,'Rank1Residuals'] = NA
        } else {
            cl.df[i,'Rank1Residuals'] = rank1.residuals(m[,cl.members[[i]] ]) 
        }
    }
    
    
    ## largest cluster with the smallest rank1 residuals
    # largest = which(cl.df$size == max(cl.df$size))
    # cId = which.min(cl.df[largest, 'Rank1Residuals'])
    ## best cluster scored by both rank1residuals and size
    cId = best.cluster(cl.df, features = c('Rank1Residuals', 'size'), weights = c(-0.7, 0.3))
    # if (debug) {print(cl.df)}
    return(list('Id'= cl.members[[cId]],
                'Name'= cl.memNames[[cId]],
                'nVertices' = ncol(m),
                'nVertices.compressed' = ncol(m.compressed)))
}

best.cluster <- function(clusters.df, features = c('Rank1Residuals', 'size'), weights = c(-0.5, 0.5)) {
    cl.df <- standardize(clusters.df[,features], na.rm = TRUE)
    score <- as.matrix(cl.df) %*% weights
    # print(data.frame(cl.df, score = score))
    return(which.max(score))
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

#' Normalize by Upper Quantile
#' 
#' @import edgeR
#' @export
normalize.by.uq <- function(X, ...) {
    effLibSizes = apply(X, 1, sum) * edgeR::calcNormFactors(t(X), method = 'upperquartile', ...) # effective library sizes
    sweep(X, 1, mean(effLibSizes) / effLibSizes, "*") %>%
      return()
}

#' Call calcNormFactors by edgeR
#' 
#' @import edgeR
#' @export
normalize.by.tmm <- function(X,...) {
    effLibSizes = apply(X, 1, sum) * edgeR::calcNormFactors(t(X), method = 'TMM', group = group, ...) # effective library sizes
    sweep(X, 1, mean(effLibSizes) / effLibSizes, "*") %>%
      return()
}

#' Call estimateSizeFactorsForMatrix by DESeq2
#' 
#' @import DESeq2
#' @param X read count matrix with samples in rows and genes in columns
#' @export
normalize.by.deseq <- function(X, ...) {
    X %>%
        t() %>%
        DESeq2::estimateSizeFactorsForMatrix(...) %>%
        sweep(X, 1, ., "/")  %>%
        return()
}

#' Normalize by PoissonSeq
#' 
#' @import PoissonSeq
#' @export
normalize.by.poissonseq <- function(X, ...) {
  PoissonSeq::PS.Est.Depth(t(X), ...) %>%
    sweep(X, 1, ., "/") %>%
    return()
}

#' Normalize by DEGES
#' 
#' @import TCC
#' @export
normalize.by.deges <- function(X, group, norm.method = 'tmm', test.method = 'edger', iteration = 1) {
    t(X) %>%
        TCC::TCC(group) %>%
        TCC::calcNormFactors(norm.method = 'tmm', test.method = 'edger', iteration = iteration) %>%
        TCC::getNormalizedData() %>%
        t() %>%
        return()
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

standardize <- function(X, na.rm = FALSE) {
    m = matrix(rep(apply(X,2,mean, na.rm = na.rm),nrow(X)),nrow=nrow(X), byrow = TRUE)
    s = matrix(rep(apply(X,2,sd, na.rm = na.rm),nrow(X)),nrow=nrow(X), byrow = TRUE)
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
