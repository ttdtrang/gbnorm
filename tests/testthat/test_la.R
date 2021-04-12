library(testthat)
library(microbenchmark)
library(futile.logger)
flog.appender(appender.file("test_la.log"), name="log")

context('Linear algebra')

test_that('pseudo-inverse of rectangular matrix', {
    x = matrix(rnorm(300*20), nrow = 20)
    t1 = microbenchmark(x.pinv <- MASS::ginv(x),
                        x.eigpinv <- eigenPseudoInverse(x),
                        times = 3)
    flog.info("Time to compute pseudo-inverse", t1, name = "log", capture = TRUE)
    expect_equal(x.pinv, x.eigpinv)
})
