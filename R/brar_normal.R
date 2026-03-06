#' @title Bayesian response-adaptive randomization for approximately normal
#'     effect estimates
#'
#' @description This function computes Bayes factors, posterior probabilities,
#'     and response-adaptive randomization probabilities for data summarized by
#'     approximately normal effect estimates.
#'
#' @param estimate Vector of effect estimates (e.g., a vector of mean
#'     differences or log odds/hazard/rate ratios). Each estimate quantifies the
#'     effect of a treatment relative to control
#' @param sigma Covariance matrix of the effect estimate vector. In case, there
#'     is only one effect estimate, this is the squared standard error of the
#'     effect estimate
#' @param pm Mean vector of the normal prior assigned to the effects under the
#'     alternative. Defaults to \code{rep(0, length(estimate))}
#' @param psigma Covariance matrix of the normal prior assigned to the effects
#'     under the alternative. In case, there is only one effect estimate, this
#'     is the prior variance
#' @param pH0 Prior probability of the point null hypothesis (i.e., all
#'     treatment effects equal to 0). Defaults to \code{0.5}. Set to \code{0} to
#'     obtain Thompson sampling and to \code{1} to obtain equal randomization
#'
#' @return An object of type \code{"brar"}, which is a list with the following
#'     elements: \code{"data"} (input data), \code{"prior"} (prior probability
#'     of the null hypothesis and prior probabilities of control/treatment
#'     superiority), \code{"BF_ij"} (Bayes factor matrix), \code{"posterior"}
#'     (posterior probability of the null hypothesis and posterior probabilities
#'     of control/treatment superiority), and \code{"prand"} (response-adaptive
#'     randomization probabilities).
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## simulate normal data from four treatment groups
#' set.seed(42)
#' n <- 10
#' muc <- 0
#' datc <- data.frame(y = rnorm(n, muc), group = "Control")
#' mu <- c(1, -0.5, 0, 0.25)
#' K <- length(mu)
#' datt <- do.call("rbind", lapply(seq(1, K), function(k) {
#'   data.frame(y = rnorm(n, mu[k]), group = paste("Treatment", k))
#' }))
#' dat <- rbind(datc, datt)
#' fit <- lm(y ~ group, data = dat)
#' estimate <- fit$coef[-1]
#' sigma <- vcov(fit)[-1,-1]
#' pm <- rep(0, K)
#'
#' ## 0.5 correlated prior to distribute prior probability equally among treatments
#' rho <- 0.5
#' psigma <- matrix(rho, nrow = K, ncol = K)
#' diag(psigma) <- 1
#' brar_normal(estimate = estimate, sigma = sigma, pm = pm, psigma = psigma,
#'             pH0 = 0.5)
#'
#' ## brar for only first treatment group
#' est <- summary(fit)$coefficients[2,1]
#' se <- summary(fit)$coefficients[2,2]
#' brar_normal(est, sigma = se^2, pm = 0, psigma = 1, pH0 = 0.5)
#'
#' @export
brar_normal <- function(estimate, sigma, pm = rep(0, length(estimate)), psigma,
                        pH0 = 0.5) {
    ## TODO maybe rename sigma and psigma to var/cov and pvar/pcov?
    ## TODO maybe a different option to input correlated covariance matrix?
    ## input checks
    stopifnot(
        length(estimate) >= 1,
        is.numeric(estimate),
        all(is.finite(estimate)),

        is.numeric(sigma),
        nrow(sigma) == length(estimate),
        ncol(sigma) == length(estimate),
        all(is.finite(sigma)),

        is.numeric(pm),
        all(is.finite(pm)),

        is.numeric(psigma),
        nrow(psigma) == length(estimate),
        ncol(psigma) == length(estimate),
        all(is.finite(psigma)),

        length(pH0) == 1,
        is.numeric(pH0),
        is.finite(pH0),
        0 <= pH0, pH0 <= 1
    )

    K <- length(estimate) # number of groups
    ## TODO rewrite so that less code for both cases
    if (K == 1) {

        ## data
        se <- sqrt(sigma)
        psd <- sqrt(psigma)
        null <- 0
        dat <- cbind("Effect estimate" = estimate, "SE" = se)
        rownames(dat) <- "Treatment 1"

        ## prior tail probability of effect > null under H1
        priorlogtail <- stats::pnorm(q = (null - pm)/psd, lower.tail = FALSE,
                                     log.p = TRUE)
        priortail <- exp(priorlogtail)

        ## posterior tail probability of effect > null under H1
        postsd <- 1/sqrt(1/se^2 + 1/psd^2)
        postm <- (estimate/se^2 + pm/psd^2)*postsd^2
        postlogtail <- stats::pnorm(q = (null - postm)/postsd, lower.tail = FALSE,
                                    log.p = TRUE)
        posttail <- exp(postlogtail)

        ## compute BFs
        logbf0p <- stats::dnorm(x = estimate, mean = null, sd = se, log = TRUE) -
            stats::dnorm(x = estimate, mean = pm, sd = sqrt(se^2 + psd^2),
                         log = TRUE) +
            priorlogtail - postlogtail
        bf0p <- exp(logbf0p)
        bfmp <- ((1 - posttail)/(1 - priortail)) / (posttail/priortail)

        ## compute prior probabilities of H+ and H-
        pHp <- priortail * (1 - pH0)
        pHm <- 1 - pH0 - pHp
        prior <- c("H-" = pHm, "H0" = pH0, "H+1" = pHp)

        ## compute BFs
        bf0m <- bf0p/bfmp
        bfpm <- 1/bfmp
        bfm0 <- bfmp/bf0p
        bfp0 <- 1/bf0p
        bfmat <- matrix(c(1, bfm0, bfmp,
                          bf0m, 1, bf0p,
                          bfpm, bfp0, 1),
                        nrow = 3, ncol = 3, byrow = TRUE)
        colnames(bfmat) <- rownames(bfmat) <- c("H-", "H0", "H+1")

        ## compute posterior probabilities of hypotheses
        if (pH0 == 1) {
            postHm <- 0
            postH0 <- 1
            postHp <- 0
        } else {
            postHm <- 1/(1 + bf0m*pH0/pHm + bfpm*pHp/pHm)
            postH0 <- 1/(bfm0*pHm/pH0 + 1 + bfp0*pHp/pH0)
            postHp <- 1/(bfmp*pHm/pHp + bf0p*pH0/pHp + 1)
        }
        ## postHm + postH0 + postHp # should be one
        post <- c("H-" = postHm, "H0" = postH0, "H+1" = postHp)

    } else { # more than one treatment group
        ## data
        null <- rep(0, K) # null value
        dat <- cbind("Effect estimate" = estimate, "SE" = sqrt(diag(sigma)))
        rownames(dat) <- paste("Treatment", seq(1, K))

        ## marginal density under H0
        margdensH0 <- mvtnorm::dmvnorm(x = estimate, mean = null, sigma = sigma)

        ## posterior mean, covariance, and marginal density under H1
        sigmainv <- solve(sigma)
        psigmainv <- solve(psigma)
        postsigma <- solve(sigmainv + psigmainv)
        postm <- as.numeric(postsigma %*% (sigmainv %*% estimate + psigmainv %*% pm))
        margdensH1 <- mvtnorm::dmvnorm(x = estimate, mean = pm, sigma = sigma + psigma)

        ## prior/posterior tail probability of H-: effect < 0 under H1
        priortailHm <- mvtnorm::pmvnorm(lower = -Inf, upper = null, mean = pm,
                                        sigma = psigma, keepAttr = FALSE)
        posttailHm <- mvtnorm::pmvnorm(lower = -Inf, upper = null, mean = postm,
                                       sigma = postsigma, keepAttr = FALSE)
        margdensHm <- margdensH1*posttailHm/priortailHm

        ## prior/posterior tail probability of H+i: effect_i > 0 and
        ## effect_i = max(effect) under H1
        resHp <- t(sapply(X = seq(1, K), FUN = function(i) {
            ## set up constraint matrix that maps H+i to negative orthant
            AiK <- matrix(rep(0, K^2), nrow = K, ncol = K)
            AiK[,i] <- -1 # contraint that effect_i > 0
            othereffects <- setdiff(seq(1, K), i) # indices of other effects than effect_i
            for (j in seq_along(othereffects)) {
                AiK[j+1,othereffects[j]] <- 1 # constraint that effect_j < effect_i
            }
            priortail <- mvtnorm::pmvnorm(lower = -Inf, upper = null,
                                          mean = as.numeric(AiK %*% pm),
                                          sigma = AiK %*% psigma %*% t(AiK),
                                          keepAttr = FALSE)
            posttail <- mvtnorm::pmvnorm(lower = -Inf, upper = null,
                                         mean = as.numeric(AiK %*% postm),
                                         sigma = AiK %*% postsigma %*% t(AiK),
                                         keepAttr = FALSE)
            margdens <- margdensH1*posttail/priortail
            c("priortail" = priortail, "posttail" = posttail, "margdens" = margdens)
        }))

        ## compute BFs
        margdens <- c(margdensHm, margdensH0, resHp[,3])
        bfmat <- outer(X = margdens, Y = margdens, FUN = `/`)
        diag(bfmat) <- 1
        Hpnames <- paste0("H+", seq(1, K))
        colnames(bfmat) <- rownames(bfmat) <- c("H-", "H0", Hpnames)

        ## compute prior hypothesis probabilities
        pHm <- priortailHm*(1 - pH0)
        pHp <- resHp[,1]*(1 - pH0)
        prior <- c(pHm, pH0, pHp)
        names(prior) <- c("H-", "H0", Hpnames)

        ## ## compute posterior hypothesis probabilitites
        ## post <- prior*margdens/sum(prior*margdens)
        ## compute posterior hypothesis probabilitites
        if (pH0 == 1) {
            ## if Pr(H0) = 1, posterior prob of H0 is always 1 regardless of data
            post <- c(0, 1, rep(0, K))
        } else if (sum(margdens == 0) > 1) {
            ## more than one marginal log likelihood = -Inf => all post probs are NaN
            post <- rep(NaN, K + 2)
        } else {
            odds <- outer(X = prior, Y = prior, FUN = `/`)
            diag(odds) <- 1
            post <- 1/sapply(X = seq_len(ncol(bfmat)), FUN = function(i) {
                sum(bfmat[,i]*odds[,i])
            })
        }
        names(post) <- names(prior)
    }
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
    prand <- post[2]/(K + 1) + post[-2]
    names(prand) <- c("Control", paste("Treatment", seq(1, K)))

    ## return brar object
    res <- list("data" = dat, "prior" = prior, "BF_ij" = bfmat,
                "posterior" = post, "prand" = prand)
    class(res) <- "brar"
    return(res)
}
