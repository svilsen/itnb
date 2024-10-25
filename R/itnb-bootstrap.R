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
            y_b <- matrix(ritnb(n = N, mu = mu, alpha = alpha, p = p, i = i, t = t), ncol = 1)
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

#' Plot histograms of bootstrapped \link{itnb-object}
#'
#' @description A function plotting bootstrapped parameter estimates returned from the \link{confint.itnb} function.
#'
#' @param x \link{itnb.ci-object}.
#' @param which String: Indicating which parameter(s) to show. If left \code{NULL}, the function shows histograms of all parameters
#' @param ... Additional arguments passed to the \link[graphics]{hist} function.
#'
#' @export
hist.itnb.ci <- function(x, which = NULL, ...) {
    #
    if (!is.na(x[["level"]])) {
        stop("'level' found in 'ci-object' implying the results have been aggregated; to use function re-run 'ci-object' setting 'level = NULL'.")
    }

    #
    beta_e <- x[["ci"]][["beta"]]
    betas <- colnames(beta_e)

    #
    alpha_e <- x[["ci"]][["alpha"]]

    #
    p_e <- x[["ci"]][["p"]]

    #
    if (is.null(which)) {
        #
        hist(alpha_e, breaks = "fd", xlab = bquote(alpha), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
        invisible(readline(prompt="Press [ENTER] to continue"))

        hist(p_e, breaks = "fd", xlab = bquote(pi), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
        invisible(readline(prompt="Press [ENTER] to continue"))

        #
        for (i in seq_along(betas)) {
            hist(beta_e[, i], breaks = "fd", xlab = bquote(beta[.(i - 1)] * ": Covariate '" * .(betas[i]) * "'"), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
            invisible(readline(prompt="Press [ENTER] to continue"))
        }
    }
    else if (all(which %in% c("mu", "beta", "covariates"))) {
        for (i in seq_along(betas)) {
            hist(beta_e[, i], breaks = "fd", xlab = bquote(beta[.(i - 1)] * ": Covariate '" * .(betas[i]) * "'"), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
            invisible(readline(prompt="Press [ENTER] to continue"))
        }
    }
    else if (all(which %in% betas)) {
        betas_ <- betas[betas %in% which]
        betas_index_ <- which(betas %in% which)
        for (i in seq_along(betas_)) {
            hist(beta_e[, betas_[i]], breaks = "fd", xlab = bquote(beta[.(betas_index_[i] - 1)] * ": Covariate '" * .(betas_[i]) * "'"), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
            if (length(betas_) > 1) {
                invisible(readline(prompt="Press [ENTER] to continue"))
            }
        }
    }
    else if (all(which %in% c("alpha", "overdispersion"))) {
        hist(alpha_e, breaks = "fd", xlab = bquote(alpha), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
    }
    else if (all(which %in% c("p", "pi", "inflation"))) {
        hist(p_e, breaks = "fd", xlab = bquote(pi), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
    }
    else {
        stop("'which' has to be NULL, or specify the name of a parameter.")
    }

    return(invisible(NULL))
}

#' Quantiles of bootstrapped \link{itnb-object}
#'
#' @description A function plotting quantiles of parameter estimates returned from the \link{confint.itnb} function.
#'
#' @param x \link{itnb.ci-object}.
#' @param which String: Indicating the column of the trace to be shown (i.e.\ the log-likelihood or the name of a parameter). If left \code{NULL}, the function shows all traces.
#' @param ... Additional arguments passed to the \link[graphics]{hist} function.
#'
#' @export
plot.itnb.ci <- function(x, which = NULL, ...) {
    #
    dots <- list(...)
    if (is.na(x[["level"]]) & is.null(dots[["level"]])) {
        stop("'level' not found in either 'ci-object', or as part of '...' argument; to use function re-run 'ci-object', or 'plot' setting a value for 'level'.")
    }
    else if (is.na(x[["level"]])) {
        #
        level <- dots[["level"]]
        sig_level <- (1 - level) / 2

        #
        beta_e <- x[["ci"]][["beta"]] |> apply(2, quantile, probs = c(sig_level, 0.5, 1 - sig_level))
        betas <- colnames(beta_e)

        #
        alpha_e <- x[["ci"]][["alpha"]] |> quantile(probs = c(sig_level, 0.5, 1 - sig_level))

        #
        p_e <- x[["ci"]][["p"]] |> quantile(probs = c(sig_level, 0.5, 1 - sig_level))
    }
    else {
        beta_e <- x[["ci"]][["beta"]]
        betas <- colnames(beta_e)

        #
        alpha_e <- x[["ci"]][["alpha"]]

        #
        p_e <- x[["ci"]][["p"]]
    }

    #
    if (is.null(which)) {
        #
        plot(alpha_e[2], 0,
             xlim = c(alpha_e[1] - 0.05 * alpha_e[1], alpha_e[3] + 0.05 * alpha_e[3]), ylim = c(-1, 1),
             xlab = bquote(alpha), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = alpha_e[1], y0 = 0, x1 = alpha_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
        invisible(readline(prompt="Press [ENTER] to continue"))

        plot(p_e[2], 0,
             xlim = c(p_e[1] - 0.05 * p_e[1], p_e[3] + 0.05 * p_e[3]), ylim = c(-1, 1),
             xlab = bquote(pi), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = p_e[1], y0 = 0, x1 = p_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)

        invisible(readline(prompt="Press [ENTER] to continue"))

        #
        for (i in seq_along(betas)) {
            plot(beta_e[2, i], 0,
                 xlim = c(beta_e[1, i] - 0.05 * beta_e[1, i], beta_e[3, i] + 0.05 * beta_e[3, i]),
                 ylim = c(-1, 1),
                 xlab = bquote(beta[.(i - 1)] * ": Covariate '" * .(betas[i]) * "'"), ylab = "", main = paste(paste0(100 * x$level, "%"), ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
                 pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
                 yaxt = "n", frame.plot = FALSE)
            arrows(x0 = beta_e[1, i], y0 = 0, x1 = beta_e[3, i], y1 = 0, code = 3, angle = 90, length = 0.2)

            invisible(readline(prompt="Press [ENTER] to continue"))
        }
    }
    else if (all(which %in% c("mu", "beta", "covariates"))) {
        for (i in seq_along(betas)) {
            plot(beta_e[2, i], 0,
                 xlim = c(beta_e[1, i] - 0.05 * beta_e[1, i], beta_e[3, i] + 0.05 * beta_e[3, i]),
                 ylim = c(-1, 1),
                 xlab = bquote(beta[.(i - 1)] * ": Covariate '" * .(betas[i]) * "'"), ylab = "", main = paste(paste0(100 * x$level, "%"), ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
                 pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
                 yaxt = "n", frame.plot = FALSE)
            arrows(x0 = beta_e[1, i], y0 = 0, x1 = beta_e[3, i], y1 = 0, code = 3, angle = 90, length = 0.2)

            invisible(readline(prompt="Press [ENTER] to continue"))
        }
    }
    else if (all(which %in% betas)) {
        betas_ <- betas[betas %in% which]
        betas_index_ <- which(betas %in% which)
        for (i in seq_along(betas_)) {
            plot(beta_e[2, i], 0,
                 xlim = c(beta_e[1, betas_[i]] - 0.05 * beta_e[1, betas_[i]], beta_e[3, betas_[i]] + 0.05 * beta_e[3, betas_[i]]),
                 ylim = c(-1, 1),
                 xlab = bquote(beta[.(betas_index_[i] - 1)] * ": Covariate '" * .(betas_[i]) * "'"), ylab = "", main = paste(paste0(100 * x$level, "%"), ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
                 pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
                 yaxt = "n", frame.plot = FALSE)
            arrows(x0 = beta_e[1, betas_[i]], y0 = 0, x1 = beta_e[3, betas_[i]], y1 = 0, code = 3, angle = 90, length = 0.2)

            invisible(readline(prompt="Press [ENTER] to continue"))
        }
    }
    else if (all(which %in% c("alpha", "overdispersion"))) {
        plot(alpha_e[2], 0,
             xlim = c(alpha_e[1] - 0.05 * alpha_e[1], alpha_e[3] + 0.05 * alpha_e[3]), ylim = c(-1, 1),
             xlab = bquote(alpha), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = alpha_e[1], y0 = 0, x1 = alpha_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
    }
    else if (all(which %in% c("p", "pi", "inflation"))) {
        plot(p_e[2], 0,
             xlim = c(p_e[1] - 0.05 * p_e[1], p_e[3] + 0.05 * p_e[3]), ylim = c(-1, 1),
             xlab = bquote(pi), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = p_e[1], y0 = 0, x1 = p_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
    }
    else {
        stop("'which' has to be NULL, or specify the name of a parameter.")
    }

    return(invisible(NULL))
}
