#' rtnbinom
#'
#' @description Generates numbers from a t-t negative binomial distribution
#'
#' @param n The number of generated samples.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#'
#' @return Random generated sample.
#' @export
rtnbinom = function(n, t, mu, theta) {
    fill = rep(NA, n)
    for(i in 1:n){
        val = t
        while(val <= t){
            val = rnbinom(1, size = theta, mu = mu)
        }
        fill[i] = val
    }
    fill
}

#' rkitnb
#'
#' @description Generates random numbers from a k-inflated, t-t, negative binomial distribution.
#'
#' @param n The number of observations.
#' @param k The inflation point.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#'
#' @return Random generated sample.
#' @export
rkitnb <- function(n, k, t = NULL, mu, theta, p) {
    if (is.null(t)) {
        res <- c(rep(k, floor(n*p)), rnbinom(n - floor(n*p), size = theta, mu = mu))
    }
    else if (t | is.numeric(t)) {
        res <- c(rep(k, floor(n*p)), rtnbinom(n - floor(n*p), mu = mu, theta = theta, t = t))
    }

    return(res[sample(seq_along(res), replace = FALSE)])
}

#' dkitnb
#'
#' @description The density of a vector x following the k-inflated, t-t, negative binomial distribution.
#'
#' @param x A vector of quantiles.
#' @param k The inflation point.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#'
#' @return Numeric.
#' @export
dkitnb <- function(x, k, t = NULL, mu, theta, p, return.log = FALSE) {
    log_density <- theta*(log(theta) - log(theta + mu)) + lgamma(theta + x) - lgamma(theta) - lfactorial(x) + x*(log(mu) - log(theta + mu))

    if (!is.null(t)) {
        if (t & !is.numeric(t)) {
            t <- k - 1
        }

        if (is.numeric(t) & (length(t) == 1) & (t < k) & (abs(t - round(t)) < .Machine$double.eps)) {
            log_density <- log_density - log(pbeta(mu/(mu + theta), t + 1, theta))

            x_t <- (x <= t)
            log_density[x_t] <- -Inf
        }
    }

    log_res <- log(1 - p) + log_density
    x_k <- x == k
    log_res[x_k] <- log(p + (1 - p)*exp(log_density[x_k]))

    res <- if (return.log) log_res else exp(log_res)
    return(res)
}

#' pkitnb
#'
#' @description The probability mass function of the k-inflated, t-t, negative binomial distribution.
#'
#' @param q A vector of quantiles.
#' @param k The inflation point.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#'
#' @return The probability mass.
#' @export
pkitnb <- function(q, k, t = NULL, mu, theta, p, return.log = FALSE) {
    res <- sum(dkitnb(1:q, k, t, mu, theta, p, log))
    return(res)
}


## Incomplete beta functions and derivatives
.ibeta <- function (x, a, b, log.p = TRUE){
    log_ibeta <- pbeta(x, a, b, log.p = TRUE) + lbeta(a, b)

    res <- if (log.p) log_ibeta else exp(log_ibeta)
    return(res)
}

.ibeta.derivative.mu <- function (mu, theta, t) {
    exp((t + 1)*log(mu) + (theta - 1)*log(theta) - (t + theta + 1)*log(mu + theta) - .ibeta(mu/(mu + theta), t + 1, theta))
}

.ibeta.derivative.theta <- function (mu, theta, t) {
    ibetaIntegrand <- function(x) x^t*(1 - x)^(theta - 1)*log(1 - x)
    res1 <- -exp(t*log(mu) + theta*log(theta) - (t + theta + 1)*log(mu + theta) - .ibeta(mu/(mu + theta), t + 1, theta))
    res2 <- integrate(f = ibetaIntegrand, lower = 0, upper = mu/(mu + theta))$value/ exp(.ibeta(mu/(mu + theta), t + 1, theta))
    return(res1 + res2)
}

## The complete log-likelihoods, the restricted log-likelihoods, and derivatives.
.completeLogLikelihood <- function (mu, theta, p, x, z, k, x_k, t) {
    log_density <- dkitnb(x, k, mu, theta, 0, return.log = TRUE, t = t)

    if (p == 0) {
        res <- (1 - z)*(log_density)
    }
    else {
        res <- z*log(p)*sum(x_k) + (1 - z)*(log(1 - p) + log_density)
    }
    return(sum(res))
}

.restrictedLogLikelihood <- function (x, mu, theta, k, t) {
    logLike <- dkitnb(x, k, mu, theta, 0, return.log = TRUE, t = t)

    return(-sum(logLike))
}

.gradient.mu <- function (x, mu, theta, k, t) {
    gradient <- -(theta + x)/(theta + mu) + x/mu

    if (!is.null(t)) {
        if (t) {
            t <- k - 1
        }

        if (is.numeric(t) & (length(t) == 1) & (t < k) & (abs(t - round(t)) < .Machine$double.eps)) {
            gradient <- gradient - .ibeta.derivative.mu(mu, theta, t)
        }
    }

    return(sum(gradient))
}

.gradient.theta <- function (x, mu, theta, k, t) {
    gradient <- -(theta + x)/(theta + mu) + 1 + log(theta) - log(theta + mu) + digamma(theta + x) - digamma(theta)

    if (!is.null(t)) {
        if (t) {
            t <- k - 1
        }

        if (is.numeric(t) & (length(t) == 1) & (t < k) & (abs(t - round(t)) < .Machine$double.eps)) {
            gradient <- gradient - .ibeta.derivative.theta(mu, theta, t) + digamma(theta) - digamma(t + 1 + theta)
        }
    }

    return(sum(gradient))
}

##
#' @export
kitnb_optimisation.control <- function(trace = TRUE, logLikelihoodConverged = 1e-6, maxNumberOfIterations = 10000, burnInIterations = 5, saveAllValues = FALSE) {
    if (is.null(burnInIterations) || is.character(burnInIterations) || is.na(burnInIterations)) {
        burnInIterations <- maxNumberOfIterations
    }
    res <- list(trace = trace, logLikelihoodConverged = logLikelihoodConverged, maxNumberOfIterations = maxNumberOfIterations,
                burnInIterations = burnInIterations, saveAllValues = saveAllValues)
    return(res)
}

#' @title kitnb EM optimisation
#'
#' @description EM algorithm for kitnb parameter optimisation.
#'
#' @param x observed values
#' @param k the point of inflation
#' @param t the point of truncation
#' @param control a control object, see \link{kitnb_optimisation.control} for details
#'
#' @return a list of the estimated parameters.
#' @export
kitnb_optimisation <- function (x, k, t, control = kitnb_optimisation.control()) {
    x_k <- x == k
    p <- mean(x_k)
    mu <- mean(x[!x_k])
    var_x <- var(x[!x_k])

    theta <- if((var_x > mu)) mu^2/(var_x - mu) else 100
    z_star <- (p * as.numeric(x_k))/(p * as.numeric(x_k) + (1 - p) * dkitnb(x, k, t = t, mu, theta, 0, return.log = FALSE))

    completeLogLikelihood <- kinb:::.completeLogLikelihood(mu, theta, p, x, z_star, k, x_k, t)
    notConverged <- TRUE
    if (control$trace) {
        cat("Iteration :", 0, "\t Current log-likelihood:", completeLogLikelihood, "\t Absolute change in log-likelihood:", NA, "\n")
        cat("\t Parameters (Initial): ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
    }
    if (control$saveAllValues) {
        resList <- list()
        resList[[1]] <- data.frame(Iteration = 0, mu = mu, theta = theta, pi = p, LogLikelihood = completeLogLikelihood, AbsoluteChangeLogLikelihood = NA)
    }

    i = 1
    while (notConverged) {
        ## E-step
        z_star <- (p*as.numeric(x_k))/(p*as.numeric(x_k) + (1 - p)*dkitnb(x, k, t = t, mu, theta, 0, return.log = FALSE))

        ## M-step
        p <- sum(z_star*as.numeric(x_k))/sum(1 - z_star*(1 - as.numeric(x_k)))
        mu <- optim(par = mu, fn = kinb:::.restrictedLogLikelihood, gr = kinb:::.gradient.mu, x = x, theta = theta, k = k, t = t, method = "BFGS")$par
        theta <- optim(par = theta, fn = kinb:::.restrictedLogLikelihood, gr = kinb:::.gradient.theta,
                       x = x, mu = mu, k = k, t = t,
                       method = "L-BFGS-B", lower = .Machine$double.eps)$par

        completeLogLikelihoodOld <- completeLogLikelihood
        completeLogLikelihood <- kinb:::.completeLogLikelihood(mu, theta, p, x, z_star, k, x_k, t)
        changeLoglikelihood <- abs(completeLogLikelihood - completeLogLikelihoodOld)

        notConverged <- if (i < control$burnInIterations) TRUE else ((changeLoglikelihood > control$logLikelihoodConverged) && (i < control$maxNumberOfIterations))
        if (control$trace) {
            cat("Iteration :", i, "\t Current log-likelihood:", completeLogLikelihood, "\t Absolute change in log-likelihood:", changeLoglikelihood, "\n")
            cat("\t Parameters: ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
        }
        if (control$saveAllValues) {
            resList[[i + 1]] <- data.frame(Iteration = i, mu = mu, theta = theta, pi = p, LogLikelihood = completeLogLikelihood, AbsoluteChangeLogLikelihood = changeLoglikelihood)
        }

        i = i + 1
    }

    res <- list(mu = mu, theta = theta, pi = p, logLikelihood = completeLogLikelihood, converged = !notConverged,
                fullParameterSequence = if(control$saveAllValues) do.call(rbind, resList) else NULL)

    return(res)
}
