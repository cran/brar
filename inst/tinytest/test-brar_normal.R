library(tinytest)
library(brar)

## invalid inputs should throw errors
expect_error(
  brar_normal(estimate = c(0.1, 0.2), sigma = matrix(1, 1, 1),
              pm = c(0, 0), psigma = diag(2), pH0 = 0.5)
)
expect_error(
  brar_normal(estimate = NA, sigma = matrix(1, 1, 1),
              pm = 0, psigma = 1, pH0 = 0.5)
)

## K = 1 treatment group
estimate <- 0.5
sigma <- 0.1
pm <- 0
psigma <- 1
res <- brar_normal(estimate = estimate, sigma = sigma, pm = pm, psigma = psigma,
                   pH0 = 0.5)

expect_inherits(res, "brar")
expect_true(all(c("data", "prior", "BF_ij", "posterior", "prand") %in% names(res)))
expect_equal(sum(res$posterior), 1, tolerance = 1e-12)
expect_equal(length(res$prand), 2)
expect_equal(names(res$prand), c("Control", "Treatment 1"))
expect_true(all(res$prand >= 0 & res$prand <= 1))

## pH0 = 1 should give equal randomization
res2 <- brar_normal(estimate = 0.3, sigma = 0.2^2, pm = 0, psigma = 1, pH0 = 1)
expect_equal(unname(res2$posterior["H0"]), 1, tolerance = 1e-12)
expect_equal(unname(res2$prand), c(0.5, 0.5), tolerance = 1e-12)

## for pH0 = 0 randomization probabilities should equal posterior probabilities
res3 <- brar_normal(estimate = estimate, sigma = sigma, pm = pm,
                    psigma = psigma, pH0 = 0)
expect_equal(unname(res3$prand), unname(res3$posterior)[-2], tolerance = 1e-12)

## K > 1 treatment groups
estimate <- c(0.2, -0.1, 0.3)
sigma <- diag(c(0.05, 0.05, 0.05))
K <- length(estimate)
pm <- rep(0, K)
psigma <- diag(K)
res4 <- brar_normal(estimate = estimate, sigma = sigma, pm = pm,
                    psigma = psigma, pH0 = 0.5)

expect_inherits(res4, "brar")
expect_true(all(c("data", "prior", "BF_ij", "posterior", "prand") %in% names(res4)))
expect_equal(sum(res4$posterior), 1, tolerance = 1e-12)
expect_equal(length(res4$prand), K + 1)
expect_equal(names(res4$prand), c("Control", paste("Treatment", 1:K)))
expect_true(all(res4$prand >= 0 & res4$prand <= 1))
expect_equal(nrow(res4$BF_ij), K + 2)
expect_equal(ncol(res4$BF_ij), K + 2)

## pH0 = 1 should give equal randomization
res5 <- brar_normal(estimate = estimate, sigma = sigma, pm = pm,
                    psigma = psigma, pH0 = 1)
expect_equal(unname(res5$prand), rep(1 / (K + 1), K + 1), tolerance = 1e-12)

## for pH0 = 0 randomization probabilities should equal posterior probabilities
res6 <- brar_normal(estimate = estimate, sigma = sigma, pm = pm,
                    psigma = psigma, pH0 = 0)
expect_equal(unname(res6$prand), unname(res6$posterior)[-2], tolerance = 1e-12)

## Symmetry check: BF_ij * BF_ji == 1
for (i in 1:(K + 2)) {
  for (j in 1:(K + 2)) {
    if (i != j) {
      expect_equal(res6$BF_ij[i,j] * res6$BF_ij[j,i], 1, tolerance = 1e-8)
    }
  }
}
