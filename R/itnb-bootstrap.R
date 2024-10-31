#' Confidence intervals of itnb-object.
#'
#' @description Bootstrap confidence intervals of the parameters estimated by the \link{itnb} function.
#'
#' @param object An \link{itnb-object}.
#' @param level Numeric: The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.
#' @param B Numeric: The number of simulations used to create the confidence envelopes.
#' @param parametric TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?
#' @param trace Numeric (>= 0): showing a trace every \code{trace} number of iterations.
#' @param control List: A control object, see \link{itnb_control} for details, passed to the \link{itnb} function.
#'
#' @example inst/examples/simulation_itnb_example.R
#'
#' @return If \code{level = NULL} a matrix with bootstrap simulations, otherwise a matrix of lower and upper confidence limits for each parameter.
#' @export
confint.itnb <- function(object, level = 0.95, B = 200, parametric = FALSE, trace = 0, control = list()) {
    ##
    if (all(is.na(object[["data"]]))) {
        stop("The 'data' was not found. Set 'save_data = TRUE' in the 'itnb_control' function, and re-run the optimisation routine.")
    }

    ##
    if (!is.null(level)) {
        if (!is.numeric(level)) {
            stop("'level' has to be 'NULL' or numeric.")
        }
        else if ((level < 0) || (level > 1)) {
            stop("'level' has to be between 0 and 1.")
        }
    }

    ##
    if (!is.numeric(B)) {
        stop("'B' has to be numeric.")
    }
    else if (B < 1) {
        stop("'B' has to be > 0.")
    }
    B <- ceiling(B)

    ##
    if (!is.logical(parametric)) {
        stop("'parametric' has to be logical.")
    }

    ##
    if (!is.numeric(trace)) {
        stop("'trace' has to be numeric.")
    }
    else if (trace < 0) {
        stop("'trace' has to be >= 0.")
    }
    trace <- ceiling(trace)

    ##
    X <- object[["data"]][["X"]]
    y <- object[["data"]][["y"]]

    N <- length(y)
    i <- object[["i"]]
    t <- object[["t"]]
    link <- object[["link"]]

    ##
    beta <- object[["beta"]]
    alpha <- object[["alpha"]]
    p <- object[["p"]]

    if (parametric) {
        mu <- X %*% beta

        if (link == "sqrt") {
            mu <- mu * mu
        }
        else if (link == "log") {
            mu <- exp(mu)
        }
    }

    ##
    beta_e <- matrix(NA, nrow = B, ncol = length(beta))
    alpha_e <- rep(NA, B)
    p_e <- rep(NA, B)
    for (b in seq_len(B)) {
        if ((trace > 0) && ((b == 1) || ((b %% trace) == 0) || (b == B))) {
            cat("Iteration:", b, "/", B, "\n")
        }

        if (parametric) {
            X_b <- X
            y_b <- ritnb(n = N, mu = mu, alpha = alpha, p = p, i = i, t = t) |> matrix(ncol = 1)
        } else {
            i_b <- sample(N, N, replace = TRUE)

            X_b <- X[i_b, , drop = FALSE]
            y_b <- y[i_b, , drop = FALSE]
        }

        pars_b <- itnb_matrix(X = X_b, y = y_b, i = i, t = t, link = link, control = control)

        beta_e[b, ] <- pars_b[["beta"]]
        alpha_e[b] <- pars_b[["alpha"]]
        p_e[b] <- pars_b[["p"]]
    }

    ci_list <- list(
        "beta" = structure(beta_e, .Dimnames = list(NULL, names(beta))),
        "alpha" = alpha_e,
        "p" = p_e
    )

    if (!is.null(level)) {
        sig_level <- (1 - level) / 2
        ci_list <- list(
            "beta" = ci_list[["beta"]] |> apply(2, quantile, probs = c(sig_level, 1 - sig_level)) |> t(),
            "alpha" = quantile(ci_list[["alpha"]], probs = c(sig_level, 1 - sig_level)),
            "p" = quantile(ci_list[["p"]], probs = c(sig_level, 1 - sig_level))
        )
    }

    res <- list(
        ci = ci_list,
        parametric = parametric,
        level = ifelse(is.null(level), NA, level),
        B = B
    )

    class(res) <- "itnb.ci"
    return(res)
}


#' Likelihood ratio tests
#'
#' @description Likelihood ratio tests for the need of inflation and/or overdispersion in inflated and truncated negative binomial regression models.
#'
#' @param object An \link{itnb-object}.
#' @param type String: The likelihood ratio test performed, either \code{"overdispersion"}, \code{"inflation"}, or \code{"both"}.
#' @param level Numeric: The significance level.
#' @param ... Additional parameters (see details).
#'
#' @details ...
#'
#' @return An \link{lrtest-object}.
#'
#' @export
lrtest <- function(object, type = "overdispersion", level = 0.05, ...) {
    UseMethod("lrtest")
}

#' @rdname lrtest
#' @method lrtest itnb
#'
#' @export
lrtest.itnb <- function(object, type = "overdispersion", level = 0.05, ...) {
    dots <- list(...)
    if (is.null(dots[["control"]])) {
        control <- do.call(itnb_control, list())
    }
    else {
        control <- dots[["control"]]
    }

    ##
    X <- object[["data"]][["X"]]
    y <- object[["data"]][["y"]]

    N <- length(y)
    i <- object[["i"]]
    t <- object[["t"]]
    link <- object[["link"]]

    ##
    loglike_o <- object[["loglikelihood"]]
    if (type %in% c("o", "overdispersion")) {
        type <- "overdisperion"

        poisson_model <- itnb_matrix(X = X, y = y, i = i, t = t, link = link, control = control)
        loglike_s <- object[["loglikelihood"]]
    }
    else if (type %in% c("i", "inflation")) {
        type <- "inflation"

        without_inflation <- itnb_matrix(X = X, y = y, i = -1, t = t, link = link, control = control)
        loglike_s <- without_inflation[["loglikelihood"]]
    }
    else if (type %in% c("b", "both")) {
        type = "both"
        loglike_s <- object[["loglikelihood"]]
    }
    else {
        stop("'type' only takes the values 'overdispersion', 'inflation', or 'both'.")
    }

    d <- 2.0 * (loglike_o - loglike_s)

    crit_val <- qchisq(1.0 - 2.0 * level, df = 1)
    p_val <- pchisq(d, df = 1, lower.tail = FALSE) / 2

    res <- list(
        lr = d,
        df = 1,
        critval = crit_val,
        pval = p_val,
        type = type,
        level = level
    )

    class(res) <- "lrtest.itnb"
    return(res)
}

#' @export
print.lrtest.itnb <- function(x, ...) {
    dots <- list(...)
    if (!("digits" %in% dots)) {
        dots[["digits"]] <- max(3, getOption("digits") - 3)
    }

    cat("\n")

    ##
    #
    if (x[["type"]] %in% c("overdisperion")) {
        cat("Likelihood ratio test of H0: alpha = 0\n")
    }

    ##
    #
    if (x[["type"]] %in% c("inflation")) {
        cat("Likelihood ratio test of H0: p = 0\n")
    }

    ##
    #
    if (x[["type"]] %in% c("both")) {
        cat("Likelihood ratio test of H0: alpha = 0 and p = 0\n")
    }

    cat(" ", paste0("Critical value (level = ", round(x[["level"]], dots[["digits"]]), "): "), round(x[["critval"]], dots[["digits"]]), "\n")
    cat(" ", "LR test statistic: ", round(x[["lr"]], dots[["digits"]]), "\n")
    cat(" ", "P-value:", format.pval(x[["pval"]], digits = dots[["digits"]]), "\n")

    cat("\n")

    return(invisible(NULL))
}
