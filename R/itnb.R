rtnbinom = function(n, mu, theta, t) {
    fill = rep(NA, n)
    for(j in seq_len(n)){
        val = t
        while(val <= t){
            val = rnbinom(n = 1, size = theta, mu = mu)
        }

        fill[j] = val
    }

    fill
}

#' ritnb
#'
#' @description Generates \code{n} random numbers from a i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{theta}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param n Numeric: The number of observations.
#' @param mu Numeric: The expected value.
#' @param theta Numeric: The overdispersion.
#' @param p Numeric: The inflation proportion.
#' @param i Numeric: The inflation point.
#' @param t Numeric: The truncation point.
#'
#' @return A vector of size \code{n} containing realisations of a itnb distribution.
#' @export
ritnb <- function(n, mu, theta, p, i, t) {
    ##
    if (!is.numeric(n)) {
        stop("'n' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }

    if (!is.numeric(theta)) {
        stop("'theta' has to be numeric.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }

    ##
    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    i <- ceiling(i)

    if (!(is.numeric(t) || is.null(t))) {
        stop("'t' has to be numeric.")
    }

    ##
    if (is.null(t)) {
        res <- c(rep(i, floor(n*p)), rnbinom(n = n - floor(n * p), size = theta, mu = mu))
    }
    else {
        t <- ceiling(t)
        res <- c(rep(i, floor(n*p)), rtnbinom(n = n - floor(n * p), mu = mu, theta = theta, t = t))
    }

    ##
    return(res[sample(seq_along(res), replace = FALSE)])
}

#' ditnb
#'
#' @description The probability of each element in a vector \code{x} following the i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{theta}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param x Numeric: A vector of quantiles.
#' @param mu Numeric: The expected value.
#' @param theta Numeric: The overdispersion.
#' @param p Numeric: The inflation proportion.
#' @param i Numeric: The inflation point.
#' @param t Numeric: The truncation point (NB: can be set to \code{NULL}).
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of\eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should the logarithm of the probabilities be returned?
#'
#' @return A vector the size as \code{x} containing the probability of each value.
#' @export
ditnb <- function(x, mu, theta, p, i = 0, t = 0, lower_tail = TRUE, return_log = FALSE) {
    ##
    if (!is.numeric(x)) {
        stop("'x' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }

    if (!is.numeric(theta)) {
        stop("'theta' has to be numeric.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }

    ##
    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    i <- ceiling(i)

    if (!(is.numeric(t) || is.null(t))) {
        stop("'t' has to be numeric.")
    }

    ##
    log_density <- theta * (log(theta) - log(theta + mu)) +
        lgamma(theta + x) - lgamma(theta) - lfactorial(x) +
        x * (log(mu) - log(theta + mu))

    if (!is.null(t)) {
        t <- ceiling(t)
        log_density <- log_density - log(pbeta(mu/(mu + theta), t + 1, theta))

        x_t <- (x <= t)
        log_density[x_t] <- -Inf
    }

    ##
    log_res <- log(1 - p) + log_density
    x_i <- x == i
    log_res[x_i] <- log(p + (1 - p) * exp(log_density[x_i]))

    ##
    res <- log_res
    if (!return_log) {
        res <- exp(log_res)
    }

    return(res)
}

#' pitnb
#'
#' @description The cumulative probability of each element in a vector \code{x} following the i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{theta}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param q Numeric: A vector of quantiles.
#' @param mu Numeric: The expected value.
#' @param theta Numeric: The overdispersion.
#' @param p Numeric: The inflation proportion.
#' @param i Numeric: The inflation point.
#' @param t Numeric: The truncation point.
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of \eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should log of the cmf be returned?
#'
#' @return A vector the size as \code{q} containing the cumulative probability of each value.
#' @export
pitnb <- function(q, mu, theta, p, i, t = NULL, lower_tail = TRUE, return_log = FALSE) {
    ##
    if (!is.numeric(q)) {
        stop("'q' has to be numeric.")
    }

    ##
    if (length(mu) == 1) {
        mu_j <- mu
    }

    if (length(theta) == 1) {
        theta_j <- theta
    }

    if (length(p) == 1) {
        p_j <- p
    }

    if (length(i) == 1) {
        i_j <- i
    }

    if (length(t) == 1) {
        t_j <- t
    }

    ##
    res <- rep(NA, length(q))
    for (j in seq_along(res)) {
        ##
        if (length(mu) > 1) {
            mu_j <- mu[j]
        }

        if (length(theta) > 1) {
            theta_j <- theta[j]
        }

        if (length(p) > 1) {
            p_j <- p[j]
        }

        if (length(i) > 1) {
            i_j <- i[j]
        }

        if (length(t) > 1) {
            t_j <- t[j]
        }

        ##
        res[j] <- sum(ditnb(x = seq_len(q[j]), mu = mu_j, theta = theta_j, p = p_j, i = i_j, t = t_j, return_log = return_log))
    }

    return(res)
}

