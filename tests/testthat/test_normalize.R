library(testthat)
library(magrittr)
library(futile.logger)

flog.appender(appender.file("test_metrics.log"), name = "log")
context('Normalization')

test_that('get references by sbm', {
    NSAMPLES = 10
    X = runif(n = NSAMPLES, min = 50, max = 1e5) %>%
        cbind(. + rnorm(NSAMPLES, 5, sd = 1e4), .) %>%
        cbind(.[,1] + rnorm(NSAMPLES, 5, sd = 1e4))
    X = rnorm(n = 2000, mean = 500, sd = 500) %>%
            matrix(nrow = NSAMPLES) %>%
            cbind(X)
        
    refs <- get.references.blocks(X, cor.method='pearson',
                                   n.runs= 1, min.size = 5,
                                   min.corr = 0.99999,
                                   min.count = 0,
                                   log.base = 2
                                   )
    expect_type(refs, "list")
    expect_type(refs$nVertices, "integer")
    expect_true(refs$nVertices.compressed <= refs$nVertices)
    expect_true(is.null(refs$Id))
    
    refs <- get.references.blocks(X, cor.method='pearson',
                                  n.runs= 3, min.size = 3,
                                  min.corr = 0.,
                                  min.count = 0,
                                  log.base = 0)
    
    # expect_type(refs$Id, "integer")
})
