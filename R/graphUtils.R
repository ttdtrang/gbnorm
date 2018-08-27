#' Community internal strength
#' @import igraph
community.strength.internal <- function(graph, community) {
    # weighted sum of edges that have both ends inside the community
    gC = igraph::induced_subgraph(graph, community)
    return(Reduce(sum, igraph::edge_attr(gC,'weight',index=igraph::E(gC))))
}

#' @import igraph
community.internal.density <- function(graph, community, weighted=FALSE) {
    gC = igraph::induced_subgraph(graph, community)
    nC = igraph::vcount(gC)
    avgW = ifelse(weighted,
                  Reduce(sum, igraph::edge_attr(graph, 'weight', index=igraph::E(graph)))/ igraph::ecount(graph),
                  1)
    intW = ifelse(weighted,
                  Reduce(sum, igraph::edge_attr(gC, 'weight', index=igraph::E(gC))),
                  igraph::ecount(gC))
    maxDensity = ifelse(igraph::is.directed(gC),
                        nC*(nC-1),
                        nC*(nC-1)/2) 
    return(intW / avgW / maxDensity)
}

#' calculate M measure, defined as
#' 
#' M(C) = (vol(C) - cut(C)) / (2 * cut(C))
#' graph    an igraph object
#' community    a vector of vertex indices
measure.M <- function(graph, community) {
   ccut = community.cut(graph, community,weighted = TRUE)
   return((community.volume(graph, community,weighted = TRUE) - ccut)/ (2. * ccut))
}

#' Community volume
#' @import igraph 
community.volume <- function(graph, community, weighted=FALSE) {
    if (weighted) {
        return(Reduce(sum,igraph::strength(graph, community)))
    } else {
        return(Reduce(sum,igraph::degree(graph, community)))
    }
}

#' Community cut
#' 
#' @import igraph
community.cut <- function(graph, community, weighted=FALSE) {
    shell = community.shell(graph, community)
    possibleEdges = expand.grid(shell,community)
    cut = 0
    if (nrow(possibleEdges) >= 1) {
        for (i in 1:nrow(possibleEdges)) {
            tryCatch({eIdx = igraph::get.edge.ids(graph, vp=possibleEdges[i,])},
                      error = function(e) {
                          message(possibleEdges[i,])
                          stop(e)
            })
            if (eIdx != 0) {
                if (weighted) { 
                    cut=cut+igraph::get.edge.attribute(graph, 'weight', eIdx)
                } else { 
                    cut = cut + 1
                }
            } 
        }
    }
    return(cut)
}

#' community.shell
#' 
#' @import igraph
#' @export
community.shell <- function(graph, community) {
    neighs = (lapply(community, function(i) {return(igraph::neighbors(graph,i))}))
    shell = setdiff(Reduce(union,neighs),community)
    return(shell)
}

#' Community detection by GCE
#' 
#' Starting from a set of seeds, a vertex in the shell is added to the community if it increases community quality metric M
#' @import parallel
#' @export
community.gce <- function(graph, seed,ncores=1, version=1) {
    indicator = rep(0, igraph::vcount(graph))
    indicator[seed] = 1
    if (version == 1) {
        if (ncores > 1) {
            cl = parallel::makeCluster(ncores)
            parallel::clusterExport(cl,c('measure.M', 'community.cut', 'community.volume', 'community.shell'))
            return(community.gce.ind(graph, indicator,cl=cl))
        } else {
            return(community.gce.ind(graph, indicator))
        }
    } else {
        return(community.gce.ind2(graph, indicator))
    }
}


community.gce.ind <- function(graph, commIndicator, cluster=NULL) {
    comm = which(commIndicator == 1)
    shell = community.shell(graph,comm)
    if (length(shell) < 1) {
        # this check is pretty stupid, but must be here because R unpredictably loop through the vector 1:0 for once
        return(comm)  
    } else {
        cQuality = measure.M(graph,comm)
        if (is.null(cluster)) {
            improvements = sapply(1:length(shell),FUN = function(i) {
                return(measure.M(graph, c(comm,shell[i])))
            })  
        } else {
            improvements = parallel::parSapply(cluster, 1:length(shell),FUN = function(i) {
                return(measure.M(graph, c(comm,shell[i])))
            })  
        }
        improvements = improvements - cQuality
        if (all(improvements <= 0)) {
            # parallel::stopCluster(cluster)
            return(which(commIndicator == 1))
        } else {
            commIndicator[shell[which(improvements == max(improvements))]] = 1
            community.gce.ind(graph, commIndicator)
        }   
    }
}

community.gce.ind2 <- function(graph, commIndicator) {
    comm = which(commIndicator == 1)
    shell = community.shell(graph, comm)
    
    if (length(shell) < 1) {
        return(comm)
    }
    # equivalent to maximize M(seed +v) - M(seed)
    mdiff = rep(0,length(shell))
    for (i in 1:length(shell)) {
        insideConn = 0
        outsideConn = 0
        iNeighbors = igraph::neighbors(graph, shell[i])
        for (j in iNeighbors) {
            w = igraph::get.edge.attribute(graph, 'weight', index=igraph::get.edge.ids(graph, vp=c(shell[i],j)))
            if (j %in% comm) {
                insideConn = insideConn + w
            } else {
                outsideConn = outsideConn + w
            }
        }
        mdiff[i] = (insideConn + outsideConn) / (outsideConn - insideConn)
    }
    print(shell)
    print(mdiff)
    if (all(mdiff <= 0)) {
        return(which(commIndicator == 1))
    } else {
        commIndicator[shell[which(mdiff == max(mdiff))]] = 1
        community.gce.ind2(graph, commIndicator)
    }
}