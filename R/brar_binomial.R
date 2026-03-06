## function to compute Pr(Xi = max(X))
## where X1 ~ Beta(a_1,b_1), ..., X_K ~ Beta(a_K,b_K) independent
Pmaxi <- function(a, b, i, ...) {
    intFun. <- function(t) {
        ## stats::dbeta(t, a[i], b[i])*prod(stats::pbeta(t, a[-i], b[-i]))
        exp(stats::dbeta(t, a[i], b[i], log = TRUE) +
            sum(stats::pbeta(t, a[-i], b[-i], log.p = TRUE)))
    }
    intFun <- Vectorize(FUN = intFun.)
    res <- try(stats::integrate(f = intFun, lower = 0, upper = 1, ... = ...)$value)
    if (inherits(res, "try-error")) {
        res <- NaN
    }
    return(res)
}

#' @title Bayesian response-adaptive randomization for binomial outcomes
#'
#' @description This function computes Bayes factors, posterior probabilities,
#'     and response-adaptive randomization probabilities for binomial outcomes.
#'
#' @param y Vector with number of successes in each group. The first element
#'     corresponds to the control group, and the remaining elements correspond
#'     to the treatment groups
#' @param n Vector with number of trials in each group. The first element
#'     corresponds to the control group, and the remaining elements correspond
#'     to the treatment groups
#' @param a0 Number of successes parameter of beta prior for common probability
#'     under the null hypothesis. Defaults to \code{1}
#' @param b0 Number of failures parameter of beta prior for common probability
#'     under the null hypothesis. Defaults to \code{1}
#' @param a Vector of number of successes parameters of beta priors for
#'     probabilities in each group under the alternative hypothesis. The first
#'     element corresponds to the control group, and the remaining elements
#'     correspond to the treatment groups. Defaults to \code{rep(1, length(y))}
#' @param b Vector of number of failures parameters of beta priors for
#'     probabilities in each group under the alternative hypothesis. The first
#'     element corresponds to the control group, and the remaining elements
#'     correspond to the treatment groups. Defaults to \code{rep(1, length(y))}
#' @param pH0 Prior probability of the null hypothesis (i.e., a common
#'     probability in the control and all treatment groups). Defaults to
#'     \code{0.5}. Set to \code{0} to obtain Thompson sampling and \code{1} to
#'     obtain equal randomization
#' @param ... Other arguments passed to \code{stats::integrate}
#'
#' @inherit brar_normal return
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## 1 control and 1 treatment group
#' y <- c(15, 12)
#' n <- c(20, 15)
#' brar_binomial(y = y, n = n, pH0 = 0.5)
#'
#' ## 1 control and 5 treatment groups
#' y <- c(10, 10, 10, 10, 10, 10)
#' n <- c(15, 15, 20, 17, 13, 25)
#' brar_binomial(y = y, n = n, pH0 = 0.5)
#'
#' @export
brar_binomial <- function(y, n, a0 = 1, b0 = 1, a = rep(1, length(y)),
                          b = rep(1, length(y)), pH0 = 0.5, ...) {

     ## input checks
    stopifnot(
        is.numeric(y),
        all(is.finite(y)),
        all(y >= 0),

        length(n) == length(y),
        is.numeric(n),
        all(is.finite(n)),
        all(n >= y), all(y <= n),

        length(a0) == 1,
        is.numeric(a0),
        is.finite(a0),
        a0 > 0,

        length(b0) == 1,
        is.numeric(b0),
        is.finite(b0),
        b0 > 0,

        length(a) == length(y),
        is.numeric(a),
        all(is.finite(a)),
        all(a > 0),

        length(b) == length(y),
        is.numeric(b),
        all(is.finite(b)),
        all(b > 0),

        length(pH0) == 1,
        is.numeric(pH0),
        is.finite(pH0),
        0 <= pH0, pH0 <= 1
    )

    ## data summaries
    K <- length(y)
    dat <- cbind("Events" = y, "Trials" = n, "Proportion" = y/n)
    rownames(dat) <- c("Control", paste("Treatment", seq_len(K - 1)))

    ## log marginal likelihood under H0
    logmargH0 <- lbeta(a0 + sum(y), b0 + sum(n) - sum(y)) -
        lbeta(a0, b0) ## +
        ## sum(lchoose(n, y)) # do not compute binom coefs because cancel out

    ## log marginal likelihood under H1
    apost <- a + y
    bpost <- b + n - y
    logmargH1 <- sum(lbeta(apost, bpost) - lbeta(a, b)) ## +
        ## sum(lchoose(n, y)) # do not compute binom coefs because cancel out

    ## prior and posterior tail probabilities under H1
    resHpi <- t(sapply(X = seq(1, K), FUN = function(i) {
        priortail <- Pmaxi(a = a, b = b, i = i, ... = ...)
        posttail <- Pmaxi(a = apost, b = bpost, i = i, ... = ...)
        logmargdens <- logmargH1 + log(posttail) - log(priortail)
        c("priortail" = priortail, "posttail" = posttail,
          "logmargdens" = logmargdens)
    }))

    ## compute BFs
    logmargdens <- c(resHpi[1,3], logmargH0, resHpi[-1,3])
    logbfmat <- outer(X = logmargdens, Y = logmargdens, FUN = `-`)
    Hpnames <- paste0("H+", seq(1, K - 1))
    colnames(logbfmat) <- rownames(logbfmat) <- c("H-", "H0", Hpnames)
    diag(logbfmat) <- 0
    bfmat <- exp(logbfmat)

    ## compute prior hypothesis probabilities
    pHm <- resHpi[1,1]*(1 - pH0)
    pHp <- resHpi[-1,1]*(1 - pH0)
    prior <- c(pHm, pH0, pHp)
    names(prior) <- c("H-", "H0", Hpnames)

    ## compute posterior hypothesis probabilitites
    if (pH0 == 1) {
        ## if Pr(H0) = 1, posterior prob of H0 is always 1 regardless of data
        post <- c(0, 1, rep(0, K - 1))
    } else if (sum(!is.finite(logmargdens)) > 1) {
        ## more than one marginal log likelihood = -Inf => all post probs are NaN
        post <- rep(NaN, K + 1)
    } else {
        ## margdens <- exp(logmargdens)
        ## post <- prior*margdens/sum(prior*margdens)
        ## for extreme data, works better with BFs than with marginal likelihoods
        odds <- outer(X = prior, Y = prior, FUN = `/`)
        diag(odds) <- 1
        post <- 1/sapply(X = seq_len(ncol(bfmat)), FUN = function(i) {
            sum(bfmat[,i]*odds[,i])
        })
    }
    names(post) <- names(prior)
    ## if prior is 0, posterior is always 0 regardless of data
    ## set posterior manually to avoid numerical NaN issues
    post[prior == 0] <- 0
    ## if posterior prob of any hypothesis is 1, all other probs are 0
    ## set manually to avoid some numerical NaN issues
    post1 <- post == 1
    if (any(post1) & sum(post1, na.rm = TRUE) == 1) {
        post[!post1] <- 0
        post[is.nan(post)] <- 0
    }

    ## compute randomization probabilities
    prand <- post[2]/K + post[-2]
    names(prand) <- c("Control", paste("Treatment", seq(1, K - 1)))

    ## return brar object
    res <- list("data" = dat, "prior" = prior, "BF_ij" = bfmat,
                "posterior" = post, "prand" = prand)
    class(res) <- "brar"
    return(res)
}
