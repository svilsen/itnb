#' ritnb
#'
#' @description Generates \code{n} random numbers from a i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{alpha}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param n Numeric (> 0): The number of observations.
#' @param mu Numeric (> 0): The expected value.
#' @param alpha Numeric (>= 0): The overdispersion.
#' @param p Numeric (>= 0): The inflation proportion.
#' @param i Numeric (>= 0): The inflation point.
#' @param t Numeric (>= 0): The truncation point.
#'
#' @return A vector of size \code{n} containing realisations of a itnb distribution.
#' @export
ritnb <- function(n, mu, alpha, p = NULL, i = NULL, t = NULL) {
    ##
    if (!is.numeric(n)) {
        stop("'n' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }
    else if (!((length(mu) == n) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(mu < 0.0)) {
        stop("'mu' has to be larger than 0.")
    }

    if (!is.numeric(alpha)) {
        stop("'alpha' has to be numeric.")
    }
    else if (!((length(alpha) == n) || length(alpha) == 1)) {
        stop("'alpha' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(alpha < 0.0)) {
        stop("'alpha' has to be larger than or equal to 0.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == n) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be the same length as 'x'.")
    }
    else if (!(all(p >= 0.0) && all(p <= 1.0))) {
        stop("'p' has to be between 0 and 1.")
    }

    ##
    if (is.null(i)) {
        i <- -1
    }

    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (!((length(i) == n) || length(i) == 1)) {
        stop("'i' needs to have length 1, or be equal to 'n'.")
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (!((length(t) == n) || length(t) == 1)) {
        stop("'t' needs to have length 1, or be equal to 'n'.")
    }
    t <- ceiling(t)

    ##
    res <- ritnb_cpp(n = n, mu = mu, theta = 1 / alpha, p = p, i = i, t = t)

    ##
    return(res[, 1])
}

#' ditnb
#'
#' @description The probability of each element in a vector \code{x} following the i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{alpha}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param x Numeric (>= 0): A vector of quantiles.
#' @param mu Numeric (> 0): The expected value.
#' @param alpha Numeric (>= 0): The overdispersion.
#' @param p Numeric (>= 0): The inflation proportion.
#' @param i Numeric (>= 0): The inflation point.
#' @param t Numeric (>= 0): The truncation point.
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of\eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should the logarithm of the probabilities be returned?
#'
#' @return A vector the size as \code{x} containing the probability of each value.
#' @export
ditnb <- function(x, mu, alpha, p = NULL, i = NULL, t = NULL, lower_tail = TRUE, return_log = FALSE) {
    n <- length(x)

    ##
    if (!is.numeric(x)) {
        stop("'x' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }
    else if (!((length(mu) == n) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(mu < 0.0)) {
        stop("'mu' has to be larger than 0.")
    }

    if (!is.numeric(alpha)) {
        stop("'alpha' has to be numeric.")
    }
    else if (!((length(alpha) == n) || length(alpha) == 1)) {
        stop("'alpha' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(alpha < 0.0)) {
        stop("'alpha' has to be larger than or equal to 0.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == n) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be the same length as 'x'.")
    }
    else if (!(all(p >= 0.0) && all(p <= 1.0))) {
        stop("'p' has to be between 0 and 1.")
    }

    ##
    if (is.null(i)) {
        i <- -1
    }

    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (!((length(i) == n) || length(i) == 1)) {
        stop("'i' needs to have length 1, or be equal to 'n'.")
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (!((length(t) == n) || length(t) == 1)) {
        stop("'t' needs to have length 1, or be equal to 'n'.")
    }
    t <- ceiling(t)

    ##
    log_res <- ditnb_cpp(x = x, mu = mu, theta = 1 / alpha, p = p, i = i, t = t)

    ##
    res <- log_res[, 1]
    if (!return_log) {
        res <- exp(res)
    }

    if (!lower_tail) {
        if (return_log) {
            res <- log1p(1.0 - expm1(res))
        }
        else {
            res <- 1.0 - res
        }
    }

    ##
    return(res)
}

#' pitnb
#'
#' @description The cumulative probability of each element in a vector \code{x} following the i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{alpha}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param q Numeric (>= 0): A vector of quantiles.
#' @param mu Numeric (> 0): The expected value.
#' @param alpha Numeric (>= 0): The overdispersion.
#' @param p Numeric (>= 0): The inflation proportion.
#' @param i Numeric (>= 0): The inflation point.
#' @param t Numeric (>= 0): The truncation point.
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of \eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should log of the cmf be returned?
#'
#' @return A vector the size as \code{q} containing the cumulative probability of each value.
#' @export
pitnb <- function(q, mu, alpha, p = NULL, i = NULL, t = NULL, lower_tail = TRUE, return_log = FALSE) {
    n <- length(q)

    ##
    if (!is.numeric(q)) {
        stop("'q' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }
    else if (!((length(mu) == n) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(mu < 0.0)) {
        stop("'mu' has to be larger than 0.")
    }

    if (!is.numeric(alpha)) {
        stop("'alpha' has to be numeric.")
    }
    else if (!((length(alpha) == n) || length(alpha) == 1)) {
        stop("'alpha' needs to have length 1, or be the same length as 'x'.")
    }
    else if (any(alpha < 0.0)) {
        stop("'alpha' has to be larger than or equal to 0.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == n) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be the same length as 'x'.")
    }
    else if (!(all(p >= 0.0) && all(p <= 1.0))) {
        stop("'p' has to be between 0 and 1.")
    }

    ##
    if (is.null(i)) {
        i <- -1
    }

    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (!((length(i) == n) || length(i) == 1)) {
        stop("'i' needs to have length 1, or be equal to 'n'.")
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (!((length(t) == n) || length(t) == 1)) {
        stop("'t' needs to have length 1, or be equal to 'n'.")
    }
    t <- ceiling(t)

    ##
    log_res <- pitnb_cpp(x = q, mu = mu, theta = 1 / alpha, p = p, i = i, t = t)

    ##
    res <- log_res[, 1]
    if (!return_log) {
        res <- exp(res)
    }

    if (!lower_tail) {
        if (return_log) {
            res <- log1p(1.0 - expm1(res))
        }
        else {
            res <- 1.0 - res
        }
    }

    ##
    return(res)
}

