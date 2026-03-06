library(tinytest)
library(brar)

## invalid inputs should throw errors
expect_error(brar_binomial(y = c(-1, 2), n = c(2, 2))) # negative y
expect_error(brar_binomial(y = c(1, 2), n = c(0, 2))) # n smaller than y
expect_error(brar_binomial(y = c(1, 2), n = c(3, 2), pH0 = -0.1)) # invalid pH0

## K = 1 treatment group
y <- c(10, 13)
n <- c(20, 21)
res <- brar_binomial(y = y, n = n, pH0 = 0.5)

expect_inherits(res, "brar")
expect_true(all(c("data", "prior", "BF_ij", "posterior", "prand") %in% names(res)))
expect_equal(sum(res$posterior), 1, tolerance = 1e-12)
expect_equal(length(res$prand), 2)
expect_equal(names(res$prand), c("Control", "Treatment 1"))
expect_true(all(res$prand >= 0 & res$prand <= 1))

## pH0 = 1 should give equal randomization
res2 <- brar_binomial(y = y, n = n, pH0 = 1)
expect_equal(unname(res2$prand), c(0.5, 0.5), tolerance = 1e-12)

## for pH0 = 0 randomization probabilities should equal posterior probabilities
res3 <- brar_binomial(y = y, n = n, pH0 = 0)
expect_equal(unname(res3$prand), unname(res3$posterior)[-2], tolerance = 1e-12)

## K = 3 treatment groups
y <- c(10, 10, 8, 12)
n <- c(20, 20, 20, 20)
K <- length(y) - 1
res4 <- brar_binomial(y = y, n = n, pH0 = 0.5)

expect_inherits(res4, "brar")
expect_equal(length(res4$prand), K + 1)
expect_equal(names(res4$prand), c("Control", paste("Treatment", 1:K)))
expect_equal(sum(res4$posterior), 1, tolerance = 1e-12)
expect_equal(sum(res4$prand), 1, tolerance = 1e-12)
expect_equal(dim(res4$BF_ij), c(K + 2, K + 2))

## pH0 = 1 should give equal randomization
res5 <- brar_binomial(y = y, n = n
                    , pH0 = 1)
expect_equal(unname(res5$prand), rep(1 / (K + 1), K + 1), tolerance = 1e-12)

## for pH0 = 0 randomization probabilities should equal posterior probabilities
res6 <- brar_binomial(y = y, n = n, pH0 = 0)
expect_equal(unname(res6$prand), unname(res6$posterior)[-2], tolerance = 1e-12)

## Symmetry check: BF_ij * BF_ji == 1
for (i in 1:(K + 2)) {
  for (j in 1:(K + 2)) {
    if (i != j) {
      expect_equal(res6$BF_ij[i,j] * res6$BF_ij[j,i], 1, tolerance = 1e-8)
    }
  }
}
