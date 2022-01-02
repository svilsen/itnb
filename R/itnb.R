#' @importFrom dplyr '%>%' 'bind_rows' 'as_tibble' 'select' 'rename'
#' @importFrom tidyr 'gather'
#' @importFrom methods 'show'
#' @importFrom stats "optim" "pbeta" "quantile" "rnbinom" "var"
#' @importFrom progress 'progress_bar'
#' @import ggplot2
NULL

globalVariables(c("..density..", "AbsoluteChangeLogLikelihood", "Iteration", "LogLikelihood", "env", "para", "p", "val"))

#' rtnbinom
#'
#' @description Generates numbers from a t-truncated negative binomial distribution.
#'
#' @param n The number of generated samples.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param t The truncation point.
#'
#' @return Random generated sample.
#' @export
rtnbinom = function(n, mu, theta, t) {
    fill = rep(NA, n)
    for(j in 1:n){
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
#' @description Generates random numbers from a i-inflated t-truncated negative binomial distribution.
#'
#' @param n The number of observations.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#' @param i The inflation point.
#' @param t The truncation point.
#'
#' @return Random generated sample.
#' @export
ritnb <- function(n, mu, theta, p, i, t = NULL) {
    if (is.null(t)) {
        res <- c(rep(i, floor(n*p)), rnbinom(n = n - floor(n * p), size = theta, mu = mu))
    }
    else if (t | is.numeric(t)) {
        res <- c(rep(i, floor(n*p)), rtnbinom(n = n - floor(n * p), mu = mu, theta = theta, t = t))
    }

    return(res[sample(seq_along(res), replace = FALSE)])
}

#' ditnb
#'
#' @description The density of a vector x following the i-inflated t-truncated negative binomial distribution.
#'
#' @param x A vector of quantiles.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#' @param i The inflation point.
#' @param t The truncation point.
#' @param return_log TRUE/FALSE, should the log pmf be returned?
#'
#' @return A vector.
#' @export
ditnb <- function(x, mu, theta, p, i, t = NULL, return_log = FALSE) {
    log_density <- theta * (log(theta) - log(theta + mu)) +
        lgamma(theta + x) - lgamma(theta) - lfactorial(x) +
        x * (log(mu) - log(theta + mu))

    if (!is.null(t)) {
        if (t & !is.numeric(t)) {
            t <- i - 1
        }

        if (is.numeric(t) & (length(t) == 1) & (t < i) & (abs(t - round(t)) < .Machine$double.eps)) {
            log_density <- log_density - log(pbeta(mu/(mu + theta), t + 1, theta))

            x_t <- (x <= t)
            log_density[x_t] <- -Inf
        }
    }

    log_res <- log(1 - p) + log_density
    x_i <- x == i
    log_res[x_i] <- log(p + (1 - p) * exp(log_density[x_i]))

    res <- log_res
    if (!return_log) {
        res <- exp(log_res)
    }

    return(res)
}

#' pitnb
#'
#' @description The probability mass function of the i-inflated t-truncated negative binomial distribution.
#'
#' @param q A vector of quantiles.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#' @param i The inflation point.
#' @param t The truncation point.
#' @param return_log TRUE/FALSE, should log of the cdf be returned?
#'
#' @return The probability mass.
#' @export
pitnb <- function(q, mu, theta, p, i, t = NULL, return_log = FALSE) {
    res <- sum(ditnb(x = seq_len(q), mu = mu, theta = theta, p = p, i = i, t = t, return_log = return_log))
    return(res)
}


## The complete and the restricted log-likelihoods
.complete_log_likelihood <- function (mu, theta, p, x, z, i, x_i, t) {
    log_density <- ditnb(x = x, mu = mu, theta = theta, p = 0, i = i, t = t, return_log = TRUE)

    if (p == 0) {
        res <- (1 - z) * log_density
    }
    else {
        res <- z * log(p) * sum(as.numeric(x_i)) + (1 - z) * (log(1 - p) + log_density)
    }

    return(sum(res))
}

.restricted_log_likelihood <- function(pars, x, p, i, t) {
    mu <- pars[1]
    theta <- pars[2]

    log_likelihood <- ditnb(x = x, mu = mu, theta = theta, p = p, i = i, t = t, return_log = TRUE)
    return(-sum(log_likelihood))
}

#' Control function for the \link{itnb_optimisation} function.
#'
#' @description Creates a list of default options.
#'
#' @param trace Numeric >= 0 showing a trace every '\code{trace}' iterations.
#' @param tolerance Convergence tolerance.
#' @param iteration_max The maximum number of allowed iterations.
#' @param iteration_min The minimum number of allowed iterations.
#' @param save_trace TRUE/FALSE should the entire trace be stored?
#'
#' @return A list of default arguments for the \link{itnb_optimisation} function.
#' @export
itnb_optimisation_control <- function(trace = 0, tolerance = 1e-6, iteration_max = 10000,
                                      iteration_min = 5, save_trace = FALSE) {
    if (is.null(iteration_min) || is.character(iteration_min) || is.na(iteration_min)) {
        iteration_min <- iteration_max
    }

    res <- list(trace = trace, tolerance = tolerance, iteration_max = iteration_max,
                iteration_min = iteration_min, save_trace = save_trace)
    return(res)
}

#' @title Parameter optimisation.
#'
#' @description EM algorithm for itnb parameter optimisation.
#'
#' @param x Observed data (numeric/integers).
#' @param i The point of inflation.
#' @param t The point of truncation.
#' @param control A control object, see \link{itnb_optimisation_control} for details.
#'
#' @return A list of the parameters: (i, t, mu, theta, p).
#' @example inst/examples/itnb_optimisation_example.R
#' @export
itnb_optimisation <- function (x, i, t, control = itnb_optimisation_control()) {
    x_i <- x == i
    p <- mean(x_i)
    mu <- mean(x[!x_i])
    var_x <- var(x[!x_i])

    theta <- ifelse(var_x > mu, mu^2 / (var_x - mu), 100)
    z_star <- (p * as.numeric(x_i)) /
        (p * as.numeric(x_i) + (1 - p) * ditnb(x = x, mu = mu, theta = theta, p = 0, i = i, t = t, return_log = FALSE))

    complete_log_likelihood <- .complete_log_likelihood(
        mu = mu, theta = theta, p = p, x = x,
        z = z_star, i = i, x_i = x_i, t = t
    )

    not_converged <- TRUE
    if (control$trace > 0) {
        cat("Iteration :", 0, "\t Current log-likelihood:", complete_log_likelihood, "\t Absolute change in log-likelihood:", NA, "\n")
        cat("\t Parameters (Initial): ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
    }

    if (control$save_trace) {
        trace_list <- list()
        trace_list[[1]] <- data.frame(
            Iteration = 0,
            mu = mu,
            theta = theta,
            p = p,
            LogLikelihood = complete_log_likelihood,
            AbsoluteChangeLogLikelihood = NA
        )
    }

    j = 1
    while (not_converged) {
        ## E-step
        z_star <- (p * as.numeric(x_i)) /
            (p * as.numeric(x_i) + (1 - p) * ditnb(x = x, mu = mu, theta = theta, p = 0, i = i, t = t, return_log = FALSE))

        ## M-step
        p <- sum(z_star * as.numeric(x_i)) / sum(1 - z_star * (1 - as.numeric(x_i)))

        pars <- optim(
            par = c(mu, theta),
            fn = .restricted_log_likelihood,
            x = x,
            p = p,
            i = i,
            t = t,
            method = "L-BFGS-B",
            lower = c(t, .Machine$double.eps),
            upper = c(Inf, Inf),
            control = list(trace = FALSE)
        )$par

        mu <- pars[1]
        theta <- pars[2]

        ## Updating convergence results
        complete_log_likelihood_old <- complete_log_likelihood
        complete_log_likelihood <- .complete_log_likelihood(mu = mu, theta = theta, p = p, x = x, z = z_star, i = i, x_i = x_i, t = t)
        change_log_likelihood <- abs(complete_log_likelihood - complete_log_likelihood_old)

        not_converged <- ifelse(j < control$iteration_min, TRUE, (change_log_likelihood > control$tolerance) && (j < control$iteration_max))

        ## Tracing
        if ((control$trace > 0) & (((j %% control$trace) == 0) || (!not_converged))) {
            cat("Iteration :", j, "\t Current log-likelihood:", complete_log_likelihood, "\t Absolute change in log-likelihood:", change_log_likelihood, "\n")
            cat("\t Parameters: ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
        }

        if (control$save_trace) {
            trace_list[[j + 1]] <- data.frame(
                Iteration = j,
                mu = mu,
                theta = theta,
                p = p,
                LogLikelihood = complete_log_likelihood,
                AbsoluteChangeLogLikelihood = change_log_likelihood
            )
        }

        j = j + 1
    }

    returned_trace <- NULL
    if (control$save_trace) {
        returned_trace <- trace_list %>% bind_rows() %>% as_tibble()
    }

    res <- list(
        n = length(x),
        i = i,
        t = t,
        mu = mu,
        theta = theta,
        p = p,
        logLikelihood = complete_log_likelihood,
        converged = !not_converged,
        trace = returned_trace
    )

    class(res) <- "itnb"
    return(res)
}

#' Confidence envelopes of itnb-object.
#'
#' @description Simulated confidence envelopes of the parameters estimated by the \link{itnb_optimisation} function.
#'
#' @param itnb_object An object of class 'itnb'.
#' @param level The confidence level.
#' @param trace TRUE/FALSE, should a trace be shown?
#' @param number_of_simulations The number of simulations used to create the confidence envelopes.
#' @param plot_simulations TRUE/FALSE, should a density plot of the parameters be returned?
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' @export
simulate_confidence_envelopes <- function(itnb_object, level = 0.95, trace = T, number_of_simulations = 101, plot_simulations = FALSE) {
    n <- itnb_object$n
    i <- itnb_object$i
    t <- itnb_object$t

    mu <- itnb_object$mu
    theta <- itnb_object$theta
    p <- itnb_object$p

    envelope_mu <- rep(NA, number_of_simulations)
    envelope_theta <- rep(NA, number_of_simulations)
    envelope_p <- rep(NA, number_of_simulations)

    if (trace) {
        pb <- progress_bar$new(
            format = paste0("Simulating confidence envelopes (size ", number_of_simulations, ")", ": [:bar] :percent Eta: :eta"),
            total = number_of_simulations,
            clear = FALSE,
            width = 120
        )
    }

    for (j in 1:number_of_simulations) {
        if (trace)
            pb$tick()

        simulation_j <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
        optimised_simulation_j <- itnb_optimisation(x = simulation_j, i = i, t = t, itnb_optimisation_control(trace = 0))

        envelope_mu[j] <- optimised_simulation_j$mu
        envelope_theta[j] <- optimised_simulation_j$theta
        envelope_p[j] <- optimised_simulation_j$p
    }

    alpha <- (1 - level) / 2

    quantiles_mu <- quantile(sort(envelope_mu), probs = c(alpha, 1 - alpha))
    quantiles_theta <- quantile(sort(envelope_theta), probs = c(alpha, 1 - alpha))
    quantiles_p <- quantile(sort(envelope_p), probs = c(alpha, 1 - alpha))

    if (plot_simulations) {
        pp <- ggplot(data.frame(env = c(envelope_mu, envelope_theta, envelope_p),
                                para = factor(rep(c("mu", "theta", "pi"), each = number_of_simulations)), levels = c("mu", "theta", "pi")),
                     aes(x = env, colour = para, fill = para)) +
            geom_density(aes(y = ..density..), alpha = 0.7) + xlab("") + ylab("Density") +
            facet_wrap(~para, scales = "free", ncol = 1, labeller = label_parsed) +
            theme_bw(base_size = 13) + theme(legend.position = "none")

        show(pp)
    }

    confint_data_frame <- rbind(quantiles_mu, quantiles_theta, quantiles_p)
    row.names(confint_data_frame) <- c("mu", "theta", "p")

    return(confint_data_frame)
}

#' Plot parameter trace of itnb-object
#'
#' @description A wrap function for ggplot2, plotting the Expectation-Maximisation parameter trace, which can returned by the \link{itnb_optimisation} function. Note: the function can only be used if the \link{itnb_optimisation} function was used with the argument \code{save_trace} set ot \code{TRUE}.
#'
#' @param itnb_object An object of the class 'itnb'.
#'
#' @return A ggplot.
#' @export
ggplot.itnb <- function(itnb_object) {
    if (is.null(itnb_object$trace))
        stop("The 'trace' tibble was not found. \nSet 'save_trace = TRUE' in the 'itnb_optimisation_control' function.")

    trace_tibble <- itnb_object$trace[-1, ] %>% rename("Log-likelihood" = LogLikelihood, "pi" = p)
    trace_tibble_long <- trace_tibble %>% select(-AbsoluteChangeLogLikelihood) %>% gather(para, val, -Iteration)

    pp <- ggplot(trace_tibble_long, aes(x = Iteration, y = val)) +
        geom_line(size = 1.05) + facet_wrap(~para, ncol = 2, nrow = 2, scales = "free_y", labeller = label_parsed) +
        scale_x_log10() + ylab("") +
        theme_bw(base_size = 13) + theme(legend.position = "none")

    return(pp)
}
