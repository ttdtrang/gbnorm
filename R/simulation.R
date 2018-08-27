# Adapted from codes provided with Evans et al. 2017
# The code from this simulation borrows substantially from
# the code used by Law et al. (2014)


#' simulate.counts
#' 
#' @param ngenes
#' @param nlibs1
#' @param nlibs2
#' @param percentDE
#' @param expected.lib.size
#' @param random.seed [NULL]
#' @param use.invChisq [TRUE] Use inverse chi-square or log-normal dispersion
#' @param symmetric Symmetric conditions have the same number of up-regulated genes
#' @return an object with `counts` and `differential`. `counts` is the samples x genes matrix
#' @import truncnorm
#' @export
simulate.counts <- function(ngenes, nlibs1, nlibs2, percentDE = 0.5,
                            expected.lib.size=rep(11e6,nlibs),
                            fc = 2,
                            random.seed = NULL, 
                            use.invChisq=TRUE, 
                            symmetric=FALSE, 
                            partiallyAsymm=TRUE,
                            completelyAsymm=FALSE) {
    
    if (!is.null(random.seed)) set.seed(random.seed) 
    nlibs = nlibs1 + nlibs2
    
    # Get distribution function of abundance proportions
    # This distribution was generated from a real dataset
    
    # Generate baseline proportions for desired number of genes
    baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
    baselineprop <- baselineprop/sum(baselineprop)

    ndiffExpress <- ngenes*percentDE
    i <- sample(1:ngenes, ndiffExpress)
    
    # baseline proportions for each gene; DE genes get multiplied by the fold change
    baselineprop1 <- baselineprop2 <- baselineprop
    
    # sequencing depths
    sizeRatio <- truncnorm::rtruncnorm(nlibs, a = 0.6, b = 1.4, mean=1, sd=0.2)
    
    # we can either have asymmetric differential expression, with unequal proportions
    # of up- and down-regulated genes, or symmetric differential expression
    if (symmetric){
      i1 <- i[1:(ndiffExpress/2)]
      i2 <- i[(ndiffExpress/2 + 1):ndiffExpress]
      
      baselineprop1[i1] <- baselineprop1[i1]*fc
      baselineprop2[i2] <- baselineprop2[i2]*fc
    } else if (completelyAsymm){
      baselineprop1[i] <- baselineprop1[i]*fc 
    } else {
      i1 <- i[1:floor(3*ndiffExpress/4)]
      i2 <- i[(floor(3*ndiffExpress/4) + 1):ndiffExpress]
      
      baselineprop1[i1] <- baselineprop1[i1]*fc
      baselineprop2[i2] <- baselineprop2[i2]*fc
    }
    
    # make these actual proportions
    baseSum1 <- sum(baselineprop1)
    baseSum2 <- sum(baselineprop2)
    baselineprop1 <- baselineprop1/sum(baselineprop1)
    baselineprop2 <- baselineprop2/sum(baselineprop2)
    
    # mean expression for each gene
    mu0.1 <- matrix(baselineprop1,ngenes,1) %*% 
      matrix(expected.lib.size[1:nlibs1]*sizeRatio[1:nlibs1],1,nlibs1)
    mu0.2 <- matrix(baselineprop2,ngenes,1) %*% 
      matrix(expected.lib.size[(nlibs1+1):(nlibs1+nlibs2)]*sizeRatio[(nlibs1+1):(nlibs1+nlibs2)],1,nlibs2)
    
    mu0 <- cbind(mu0.1,mu0.2)
    
    # differential will keep track of which genes are differentially expressed
    differential <- rep(0,ngenes)
    differential[i] <- 1
    
    # Biological variation
    BCV0 <- 0.2+1/sqrt(mu0)
    if(use.invChisq){
      df.BCV <- 40
      BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
    } else {
      BCV <- BCV0*exp( rnorm(ngenes,mean=0,sd=0.25)/2 )
    }
    if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
    shape <- 1/BCV^2
    scale <- mu0/shape
    mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
    
    # Technical variation
    counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)
    
    #Filter: following limma/voom code, remove genes with
    # rowsums < 10 in the read count matrix
    keep <- rowSums(counts)>=10
    nkeep <- sum(keep)
    counts2 <- counts[keep,]
    differential2 <- differential[keep]
    conditions <- data.frame('condition'=as.factor(c(rep(1,nlibs1), rep(2,nlibs2))))
    return(list('counts' = t(counts2),
                'differential' = differential2,
                'conditions' = conditions))
}
