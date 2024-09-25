
#' Confidence envelopes of itnb-object.
#'
#' @description Simulated confidence envelopes of the parameters estimated by the \link{em_itnb} function.
#'
#' @param object An \link{itnb-object}.
#' @param level Numeric: The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.
#' @param nr_simulations Numeric: The number of simulations used to create the confidence envelopes.
#' @param parametric TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?
#' @param trace Numeric (>= 0): showing a trace every \code{trace} number of iterations.
#'
#' @return If \code{level = NULL} a matrix with bootstrap simulations, otherwise a matrix of lower and upper confidence limits for each parameter.
#'
#' @export
simulate_ci <- function(object, level = 0.95, nr_simulations = 2000, parametric = FALSE, trace = 0) {
    UseMethod("simulate_ci")
}

#' @rdname simulate_ci
#' @method simulate_ci itnb
#'
#' @example inst/examples/simulation_itnb_example.R
#'
#' @export
simulate_ci.itnb <- function(object, level = 0.95, nr_simulations = 2000, parametric = FALSE, trace = 0) {
    ##
    x <- object$data
    n <- length(x)
    i <- object$i
    t <- object$t

    ##
    mu <- object$mu
    theta <- object$theta
    p <- object$p

    ##
    mu_e <- rep(NA, nr_simulations)
    theta_e <- rep(NA, nr_simulations)
    p_e <- rep(NA, nr_simulations)

    for (j in seq_len(nr_simulations)) {
        if ((trace > 0) && ((j == 1) || ((j %% trace) == 0) || (j == nr_simulations))) {
            cat("Iteration:", j, "/", nr_simulations, "\n")
        }

        if (parametric) {
            x_j <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
        }
        else {
            i_j <- sample(n, n, replace = TRUE)
            x_j <- x[i_j]
        }

        pars_j <- em_itnb(x = x_j, i = i, t = t)

        mu_e[j] <- pars_j$mu
        theta_e[j] <- pars_j$theta
        p_e[j] <- pars_j$p
    }

    if (is.null(level)) {
        complete_data_frame <- cbind(mu_e, theta_e, p_e)
        colnames(complete_data_frame) <- c("mu", "theta", "p")

        res_data_frame <- complete_data_frame
    }
    else {
        alpha <- (1 - level) / 2

        quantiles_mu <- quantile(mu_e, probs = c(alpha, 0.5, 1 - alpha))
        quantiles_theta <- quantile(theta_e, probs = c(alpha, 0.5, 1 - alpha))
        quantiles_p <- quantile(p_e, probs = c(alpha, 0.5, 1 - alpha))

        confint_data_frame <- rbind(quantiles_mu, quantiles_theta, quantiles_p)
        rownames(confint_data_frame) <- c("mu", "theta", "p")

        res_data_frame <- confint_data_frame
    }

    res <- list(
        ci = res_data_frame,
        parametric = parametric,
        level = level,
        nr_simulations = nr_simulations
    )

    class(res) <- ifelse(is.null(level), "ci_complete", "ci")
    return(res)
}

simulate_ci_wrapper <- function() {

}

#' @export
plot.ci_complete <- function(x, ...) {
    mu_e <- x$ci[,"mu"]
    theta_e <- x$ci[,"theta"]
    p_e <- x$ci[,"p"]

    par(mfrow = c(3, 1))
    hist(mu_e, breaks = "fd", xlab = bquote(mu), ylab = "Density", probability = TRUE, main = paste(ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5);
    hist(theta_e, breaks = "fd", xlab = bquote(theta), ylab = "Density", probability = TRUE, main = " ", cex.lab = 1.5, cex.main = 1.5)
    hist(p_e, breaks = "fd", xlab = bquote(pi), ylab = "Density", probability = TRUE, main = " ", cex.lab = 1.5, cex.main = 1.5)
    par(mfrow = c(1, 1))

    return(invisible(NULL))
}

#' @export
plot.ci <- function(x, ...) {
    mu_e <- x$ci["mu", ]
    theta_e <- x$ci["theta", ]
    p_e <- x$ci["p", ]

    par(mfrow = c(3, 1))
    plot(mu_e[2], 0,
         xlim = c(mu_e[1] - 0.05 * mu_e[1], mu_e[3] + 0.05 * mu_e[3]),
         ylim = c(-1, 1),
         xlab = bquote(mu), ylab = "", main = paste(paste0(100 * x$level, "%"), ifelse(x$parametric, "Parametric", "Non-parametric"), "bootstrap samples"),
         pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
         yaxt = "n", frame.plot = FALSE)
    arrows(x0 = mu_e[1], y0 = 0, x1 = mu_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)

    plot(theta_e[2], 0,
         xlim = c(theta_e[1] - 0.05 * theta_e[1], theta_e[3] + 0.05 * theta_e[3]), ylim = c(-1, 1),
         xlab = bquote(theta), ylab = "", main = " ",
         pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
         yaxt = "n", frame.plot = FALSE)
    arrows(x0 = theta_e[1], y0 = 0, x1 = theta_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)

    plot(p_e[2], 0,
         xlim = c(p_e[1] - 0.05 * p_e[1], p_e[3] + 0.05 * p_e[3]), ylim = c(-1, 1),
         xlab = bquote(pi), ylab = "", main = " ",
         pch = 16, cex = 2, cex.lab = 1.5, cex.main = 1.5,
         yaxt = "n", frame.plot = FALSE)
    arrows(x0 = p_e[1], y0 = 0, x1 = p_e[3], y1 = 0, code = 3, angle = 90, length = 0.2)
    par(mfrow = c(1, 1))

    return(invisible(NULL))
}
