
#' Confidence envelopes of itnb-object.
#'
#' @description Simulated confidence envelopes of the parameters estimated by the \link{em_itnb} function.
#'
#' @param object An \link{itnb-object}.
#' @param level Numeric: The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.
#' @param trace TRUE/FALSE: should a trace be shown?
#' @param nr_simulations Numeric: The number of simulations used to create the confidence envelopes.
#' @param parametric TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?
#' @param plot TRUE/FALSE: should a density plot of the parameters be returned?
#'
#' @return If \code{level = NULL} a matrix with bootstrap simulations, otherwise a matrix of lower and upper confidence limits for each parameter.
#'
#' @export
simulate_ci <- function(object, level = 0.95, trace = TRUE, nr_simulations = 2000, parametric = FALSE, plot = FALSE) {
    UseMethod("simulate_ci")
}

#' @rdname simulate_ci
#' @method simulate_ci itnb
#'
#' @example inst/examples/simulation_itnb_example.R
#'
#' @export
simulate_ci.itnb <- function(object, level = 0.95, trace = TRUE, nr_simulations = 2000, parametric = FALSE, plot = FALSE) {
    ##
    n <- object$n
    x <- object$x
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

    if (trace) {
        pb <- progress_bar$new(
            format = paste0(ifelse(parametric, "Parametric", "Non-parametric"), " bootstrap (samples = ", nr_simulations, ")", ": [:bar] :percent Eta: :eta"),
            total = nr_simulations,
            clear = FALSE,
            width = 120
        )
    }

    for (j in seq_len(nr_simulations)) {
        if (trace)
            pb$tick()

        if (parametric) {
            x_j <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
        }
        else {
            x_j <- x[sample(n, n, replace = TRUE)]
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

        quantiles_mu <- quantile(mu_e, probs = c(alpha, 1 - alpha))
        quantiles_theta <- quantile(theta_e, probs = c(alpha, 1 - alpha))
        quantiles_p <- quantile(p_e, probs = c(alpha, 1 - alpha))

        confint_data_frame <- rbind(quantiles_mu, quantiles_theta, quantiles_p)
        rownames(confint_data_frame) <- c("mu", "theta", "p")

        res_data_frame <- confint_data_frame
    }

    if (plot) {
        par(mfrow = c(3, 1))
        hist(mu_e, breaks = "fd", xlab = bquote(mu), ylab = "Density", probability = TRUE, main = paste(ifelse(parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5);
        if (!is.null(level)) abline(v = quantiles_mu[1], col = "dodgerblue2", lwd = 2); abline(v = quantiles_mu[2], col = "dodgerblue2", lwd = 2)

        hist(theta_e, breaks = "fd", xlab = bquote(theta), ylab = "Density", probability = TRUE, main = paste(ifelse(parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5)
        if (!is.null(level)) abline(v = quantiles_theta[1], col = "dodgerblue2", lwd = 2); abline(v = quantiles_theta[2], col = "dodgerblue2", lwd = 2)

        hist(p_e, breaks = "fd", xlab = bquote(pi), ylab = "Density", probability = TRUE, main = paste(ifelse(parametric, "Parametric", "Non-parametric"), "bootstrap samples"), cex.lab = 1.5, cex.main = 1.5)
        if (!is.null(level)) abline(v = quantiles_p[1], col = "dodgerblue2", lwd = 2); abline(v = quantiles_p[2], col = "dodgerblue2", lwd = 2)

        par(mfrow = c(1, 1))
    }

    return(res_data_frame)
}
