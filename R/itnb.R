#' ritnb
#'
#' @description Generates \code{n} random numbers from a i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{theta}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param n Numeric: The number of observations.
#' @param mu Numeric: The expected value.
#' @param theta Numeric: The overdispersion.
#' @param p Numeric: The inflation proportion.
#' @param i Numeric: The inflation point.
#' @param t Numeric (< i): The truncation point.
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
    else if (!((length(mu) == n) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be equal to 'n'.")
    }

    if (!is.numeric(theta)) {
        stop("'theta' has to be numeric.")
    }
    else if (!((length(theta) == n) || length(theta) == 1)) {
        stop("'theta' needs to have length 1, or be equal to 'n'.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == n) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be equal to 'n'.")
    }

    ##
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

    if (is.numeric(t)) {
        t <- ceiling(t)
        if (any(t >= i)) {
            stop("'t' has to be smaller than 'i'.")
        }
    }

    ##
    res <- ritnb_cpp(n = n, mu = mu, theta = theta, p = p, i = i, t = t, seed = sample(1e6, 1))

    ##
    return(res[, 1])
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
#' @param t Numeric (< i): The truncation point (NB: can be set to \code{NULL}).
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of\eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should the logarithm of the probabilities be returned?
#'
#' @return A vector the size as \code{x} containing the probability of each value.
#' @export
ditnb <- function(x, mu, theta, p, i, t, lower_tail = TRUE, return_log = FALSE) {
    ##
    if (!is.numeric(x)) {
        stop("'x' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }
    else if (!((length(mu) == length(x)) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be the same length as 'x'.")
    }

    if (!is.numeric(theta)) {
        stop("'theta' has to be numeric.")
    }
    else if (!((length(theta) == length(x)) || length(theta) == 1)) {
        stop("'theta' needs to have length 1, or be the same length as 'x'.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == length(x)) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be the same length as 'x'.")
    }

    ##
    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (!((length(i) == length(x)) || length(i) == 1)) {
        stop("'i' needs to have length 1, or be the same length as 'x'.")
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (!((length(t) == length(x)) || length(t) == 1)) {
        stop("'t' needs to have length 1, or be the same length as 'x'.")
    }

    if (is.numeric(t)) {
        t <- ceiling(t)
        if (any(t >= i)) {
            stop("'t' has to be smaller than 'i'.")
        }
    }

    ##
    log_res <- ditnb_cpp(x, mu, theta, p, i, t)

    ##
    res <- log_res[, 1]
    if (!return_log) {
        res <- exp(log_res)
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
#' @description The cumulative probability of each element in a vector \code{x} following the i-inflated t-truncated negative binomial distribution with parameters \code{mu}, \code{theta}, and \code{p} which is inflated at \code{i} and truncated at \code{t}.
#'
#' @param q Numeric: A vector of quantiles.
#' @param mu Numeric: The expected value.
#' @param theta Numeric: The overdispersion.
#' @param p Numeric: The inflation proportion.
#' @param i Numeric: The inflation point.
#' @param t Numeric (< i): The truncation point.
#' @param lower_tail TRUE/FALSE: should \eqn{P[X \leq x]} be returned in favour of \eqn{P[X \geq x]}?
#' @param return_log TRUE/FALSE: should log of the cmf be returned?
#'
#' @return A vector the size as \code{q} containing the cumulative probability of each value.
#' @export
pitnb <- function(q, mu, theta, p, i, t, lower_tail = TRUE, return_log = FALSE) {
    ##
    if (!is.numeric(q)) {
        stop("'q' has to be numeric.")
    }

    ##
    if (!is.numeric(mu)) {
        stop("'mu' has to be numeric.")
    }
    else if (!((length(mu) == length(x)) || length(mu) == 1)) {
        stop("'mu' needs to have length 1, or be the same length as 'x'.")
    }

    if (!is.numeric(theta)) {
        stop("'theta' has to be numeric.")
    }
    else if (!((length(theta) == length(x)) || length(theta) == 1)) {
        stop("'theta' needs to have length 1, or be the same length as 'x'.")
    }

    if (!is.numeric(p)) {
        stop("'p' has to be numeric.")
    }
    else if (!((length(p) == length(x)) || length(p) == 1)) {
        stop("'p' needs to have length 1, or be the same length as 'x'.")
    }

    ##
    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (!((length(i) == length(x)) || length(i) == 1)) {
        stop("'i' needs to have length 1, or be the same length as 'x'.")
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (!((length(t) == length(x)) || length(t) == 1)) {
        stop("'t' needs to have length 1, or be the same length as 'x'.")
    }

    if (is.numeric(t)) {
        t <- ceiling(t)
        if (any(t >= i)) {
            stop("'t' has to be smaller than 'i'.")
        }
    }

    ##
    log_res <- pitnb_cpp(q, mu, theta, p, i, t)

    ##
    res <- log_res[, 1]
    if (!return_log) {
        res <- exp(log_res)
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

