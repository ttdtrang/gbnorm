library(testthat)
library(magrittr)
library(futile.logger)

flog.appender(appender.file("test_metrics.log"), name = "log")
context('Metrics')

test_that('fast.condNumber() output is the same as that of condNumber()', {
    m = 30
    n = 1000
    N = 20
    dat = rnorm(m*n*N) %>%
        array(c(m, n, N))
    t1 = microbenchmark::microbenchmark(slow <- sapply(1:N, function(i) { dat[,,i] %>% condNumber() }),
                        fast <- sapply(1:N, function(i) { dat[,,i] %>% fast.condNumber() }),
                        times = 1
                        )
    flog.info("Time to calculate condition number", t1, name = 'log', capture = TRUE)
    expect_equal(slow, fast)

    x = matrix(rnorm(100000, 1, 1), nrow = 100)
    t2 = microbenchmark::microbenchmark(condNumber(x), fast.condNumber(x), times = 5)
    flog.info("Time to calculate condition number on 100x1000 matrix:", t2, name = 'log', capture = TRUE)
})

test_that('cdev is always calculated for the smaller transformation', {
    m = 30
    n = 100
    m.fat <- matrix(rnorm(m*n), nrow = m)
    mp.fat <- m.fat + matrix(rnorm(m*n), nrow = m)
    m.thin <- t(m.fat)
    mp.thin <- t(mp.fat)

    mmp.fat = cdev(m.fat, mp.fat, always.small = TRUE)
    mmp.thin = cdev(m.thin, mp.thin, always.small = TRUE)

    mpm.fat = cdev(mp.fat, m.fat, always.small = TRUE)
    mpm.thin = cdev(mp.thin, m.thin, always.small = TRUE)

    mmp.thin.big = cdev(m.thin, mp.thin, always.small = FALSE)

    expect_equal(mmp.fat, mmp.thin)
    expect_equal(mpm.fat, mpm.thin)
})

test_that('cpp-cdev and r-cdev are the same', {
    x = matrix(rnorm(10000, 1, 1), nrow = 10)
    y = x + matrix(rnorm(10000, 0, 3), nrow = 10)
    t1 = microbenchmark::microbenchmark(cdev.r <- cdev(x, y), cdev.c <- cppCdev(x,y), times = 3)
    flog.info("Time to calculate cdev(x,y), dim = (100x10000)", t1, name = 'log', capture = TRUE)
    expect_equal(cdev.r, cdev.c, tol = 1e-10)

    # x = matrix(rnorm(40000, 1, 1), nrow = 400)
    # y = x + matrix(rnorm(40000, 0, 1), nrow = 400)
    # t2 = microbenchmark::microbenchmark(cdev.r <- cdev(x, y), cdev.c <- cppCdev(x,y), times = 3)
    # flog.info("Time to calculate cdev(x,y), dim = (400x40000)", t2, name = 'log', capture = TRUE)
    # expect_equal(cdev.r, cdev.c)

    # expect_equal(cdev.c, cdev(y,x))
})


