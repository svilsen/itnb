#' Confidence intervals of itnb-object.
#'
#' @description Bootstrap confidence intervals of the parameters estimated by the \link{itnb} function.
#'
#' @param object An \link{itnb-object}.
#' @param level Numeric: The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.
#' @param nr_simulations Numeric: The number of simulations used to create the confidence envelopes.
#' @param parametric TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?
#' @param trace Numeric (>= 0): showing a trace every \code{trace} number of iterations.
#' @param control List: A control object, see \link{itnb_control} for details, passed to the \link{itnb} function.
#'
#' @example inst/examples/simulation_itnb_example.R
#'
#' @return If \code{level = NULL} a matrix with bootstrap simulations, otherwise a matrix of lower and upper confidence limits for each parameter.
#' @export
confint.itnb <- function(object, level = 0.95, nr_simulations = 200, parametric = FALSE, trace = 0, control = list()) {
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
    if (!is.numeric(nr_simulations)) {
        stop("'trace' has to be numeric.")
    }
    else if (nr_simulations < 1) {
        stop("'nr_simulations' has to be > 0.")
    }
    nr_simulations <- ceiling(nr_simulations)

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
    theta <- object[["theta"]]
    p <- object[["p"]]

    ##
    if (parametric) {
        r <- residuals(object, type = "response")
    }

    ##
    beta_e <- matrix(NA, nrow = nr_simulations, ncol = length(beta))
    theta_e <- rep(NA, nr_simulations)
    p_e <- rep(NA, nr_simulations)
    for (j in seq_len(nr_simulations)) {
        if ((trace > 0) && ((j == 1) || ((j %% trace) == 0) || (j == nr_simulations))) {
            cat("Iteration:", j, "/", nr_simulations, "\n")
        }

        if (parametric) {
            r_j <- sample(r, N, replace = TRUE)

            X_j <- X
            y_j <- y + r_j
        }
        else {
            i_j <- sample(N, N, replace = TRUE)

            X_j <- X[i_j, , drop = FALSE]
            y_j <- y[i_j, , drop = FALSE]
        }

        pars_j <- itnb_matrix(X = X_j, y = y_j, i = i, t = t, link = link, control = control)

        beta_e[j, ] <- pars_j[["beta"]]
        theta_e[j] <- pars_j[["theta"]]
        p_e[j] <- pars_j[["p"]]
    }

    ci_list <- list(
        "beta" = structure(beta_e, .Dimnames = list(NULL, names(beta))),
        "theta" = theta_e,
        "p" = p_e
    )

    if (!is.null(level)) {
        alpha <- (1 - level) / 2
        ci_list <- list(
            "beta" = ci_list[["beta"]] |> apply(2, quantile, probs = c(alpha, 0.5, 1 - alpha)),
            "theta" = ci_list[["theta"]] |> quantile(probs = c(alpha, 0.5, 1 - alpha)),
            "p" = ci_list[["p"]] |> quantile(probs = c(alpha, 0.5, 1 - alpha))
        )
    }

    res <- list(
        ci = ci_list,
        parametric = parametric,
        level = ifelse(is.null(level), NA, level),
        nr_simulations = nr_simulations
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
    theta_e <- x[["ci"]][["theta"]]

    #
    p_e <- x[["ci"]][["p"]]

    #
    if (is.null(which)) {
        #
        hist(theta_e, breaks = "fd", xlab = bquote(theta), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
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
    else if (all(which %in% c("theta", "overdispersion"))) {
        hist(theta_e, breaks = "fd", xlab = bquote(theta), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5, ...)
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
        alpha <- (1 - level) / 2

        #
        beta_e <- x[["ci"]][["beta"]] |> apply(2, quantile, probs = c(alpha, 0.5, 1 - alpha))
        betas <- colnames(beta_e)

        #
        theta_e <- x[["ci"]][["theta"]] |> quantile(probs = c(alpha, 0.5, 1 - alpha))

        #
        p_e <- x[["ci"]][["p"]] |> quantile(probs = c(alpha, 0.5, 1 - alpha))
    }
    else {
        beta_e <- x[["ci"]][["beta"]]
        betas <- colnames(beta_e)

        #
        theta_e <- x[["ci"]][["theta"]]

        #
        p_e <- x[["ci"]][["p"]]
    }

    #
    if (is.null(which)) {
        #
        plot(theta_e[2], 0,
             xlim = c(theta_e[1] - 0.05 * theta_e[1], theta_e[3] + 0.05 * theta_e[3]), ylim = c(-1, 1),
             xlab = bquote(theta), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = theta_e[1], y0 = 0, x1 = theta_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
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
    else if (all(which %in% c("theta", "overdispersion"))) {
        plot(theta_e[2], 0,
             xlim = c(theta_e[1] - 0.05 * theta_e[1], theta_e[3] + 0.05 * theta_e[3]), ylim = c(-1, 1),
             xlab = bquote(theta), ylab = "", main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
             pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
             yaxt = "n", frame.plot = FALSE)
        arrows(x0 = theta_e[1], y0 = 0, x1 = theta_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
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
