library(testthat)
library(microbenchmark)
library(futile.logger)
flog.appender(appender.file("test_la.log"), name="log")

context('Linear algebra')

svd.reconstruct <- function(u, d, v) {
    return(u %*% diag(d, nrow = ncol(u), ncol = ncol(v)) %*% t(v))
}

test_that('eigen svd output', {
    x = matrix(rnorm(300*20), nrow = 20)
    svd.r = svd(x)
    svd.jacobi = eigenJacobiSVD(x)
    svd.bdc = eigenBDCSVD(x)
    expect_equal(svd.r$d, svd.jacobi$d)
    expect_equal(svd.r$d, svd.bdc$d)
    expect_equal(x, svd.reconstruct(svd.jacobi$u, svd.jacobi$d, svd.jacobi$v))
    expect_equal(x, svd.reconstruct(svd.bdc$u, svd.bdc$d, svd.bdc$v))
    t2 = microbenchmark(svd(x), eigenJacobiSVD(x), eigenBDCSVD(x), times = 5)
    flog.info("Time to compute full SVD", t2, name = "log", capture = TRUE)
})


test_that('eigen functions to compute singular values', {
    x = matrix(rnorm(300*20), nrow = 20)
    svd.r = svd(x)
    sv.jacobi = eigenJacobiSV(x)
    sv.bdc = eigenBDCSV(x)
    expect_equal(svd.r$d, sv.jacobi)
    expect_equal(svd.r$d, sv.bdc)
    t2 = microbenchmark(svd(x), eigenJacobiSV(x), eigenBDCSV(x), times = 5)
    flog.info("Time to compute only singular values", t2, name = "log", capture = TRUE)
})

test_that('pseudo-inverse of rectangular matrix', {
    x = matrix(rnorm(300*20), nrow = 20)
    t1 = microbenchmark(x.pinv <- MASS::ginv(x),
                        x.eigpinv <- eigenPseudoInverse(x),
                        times = 3)
    flog.info("Time to compute pseudo-inverse", t1, name = "log", capture = TRUE)
    expect_equal(x.pinv, x.eigpinv)
})
