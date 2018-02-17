#' @importFrom dplyr '%>%' 'bind_rows' 'as_tibble' 'select' 'rename'
#' @importFrom tidyr 'gather'
#' @importFrom methods 'show'
#' @importFrom stats "optim" "pbeta" "quantile" "rnbinom" "var"
#' @import ggplot2
#'
#' @importFrom progress 'progress_bar'
NULL

globalVariables(c("..density..", "AbsoluteChangeLogLikelihood", "Iteration", "LogLikelihood", "env", "para", "p", "val"))

#' rtnbinom
#'
#' @description Generates numbers from a t-truncated negative binomial distribution.
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

#' rktnb
#'
#' @description Generates random numbers from a K-inflated t-truncated negative binomial distribution.
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
rktnb <- function(n, k, t = NULL, mu, theta, p) {
    if (is.null(t)) {
        res <- c(rep(k, floor(n*p)), rnbinom(n - floor(n*p), size = theta, mu = mu))
    }
    else if (t | is.numeric(t)) {
        res <- c(rep(k, floor(n*p)), rtnbinom(n - floor(n*p), mu = mu, theta = theta, t = t))
    }

    return(res[sample(seq_along(res), replace = FALSE)])
}

#' dktnb
#'
#' @description The density of a vector x following the k-inflated t-truncated negative binomial distribution.
#'
#' @param x A vector of quantiles.
#' @param k The inflation point.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#' @param return_log TRUE/FALSE, should the log pmf be returned?
#'
#' @return A vector.
#' @export
dktnb <- function(x, k, t = NULL, mu, theta, p, return_log = FALSE) {
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

    res <- if (return_log) log_res else exp(log_res)
    return(res)
}

#' pktnb
#'
#' @description The probability mass function of the k-inflated t-truncated negative binomial distribution.
#'
#' @param q A vector of quantiles.
#' @param k The inflation point.
#' @param t The truncation point.
#' @param mu The expected value.
#' @param theta The overdispersion.
#' @param p The proportion of inflation.
#' @param return_log TRUE/FALSE, should log of the cdf be returned?
#'
#' @return The probability mass.
#' @export
pktnb <- function(q, k, t = NULL, mu, theta, p, return_log = FALSE) {
    res <- sum(dktnb(1:q, k, t, mu, theta, p, log))
    return(res)
}


## The complete and the restricted log-likelihoods
.complete_log_likelihood <- function (mu, theta, p, x, z, k, x_k, t) {
    log_density <- dktnb(x, k, mu, theta, 0, return_log = TRUE, t = t)

    if (p == 0) {
        res <- (1 - z)*(log_density)
    }
    else {
        res <- z*log(p)*sum(as.numeric(x_k)) + (1 - z)*(log(1 - p) + log_density)
    }

    return(sum(res))
}

.restricted_log_likelihood <- function(pars, x, p, k, t) {
    mu <- pars[1]
    theta <- pars[2]

    log_likelihood <- dktnb(x, k, mu, theta, p, t = t, return_log = TRUE)
    return(-sum(log_likelihood))
}

#' Control function for the \link{ktnb_optimisation} function.
#'
#' @description Creates a list of default options.
#'
#' @param trace TRUE/FALSE show trace?
#' @param tolerance Convergence tolerance.
#' @param iteration_max The maximum number of allowed iterations.
#' @param iteration_min The minimum number of allowed iterations.
#' @param save_trace TRUE/FALSE should the entire trace be stored?
#'
#' @return A list of default arguments for the \link{ktnb_optimisation} function.
#' @export
ktnb_optimisation_control <- function(trace = TRUE, tolerance = 1e-6, iteration_max = 10000, iteration_min = 5, save_trace = FALSE) {
    if (is.null(iteration_min) || is.character(iteration_min) || is.na(iteration_min)) {
        iteration_min <- iteration_max
    }
    res <- list(trace = trace, tolerance = tolerance, iteration_max = iteration_max,
                iteration_min = iteration_min, save_trace = save_trace)
    return(res)
}

#' @title Parameter optimisation.
#'
#' @description EM algorithm for ktnb parameter optimisation.
#'
#' @param x Observed data (numeric/integers).
#' @param k The point of inflation.
#' @param t The point of truncation.
#' @param control A control object, see \link{ktnb_optimisation_control} for details.
#'
#' @return A list of the estimated parameters, (mu, theta, p).
#' @example inst/examples/ktnb_optimisation_example.R
#' @export
ktnb_optimisation <- function (x, k, t, control = ktnb_optimisation_control()) {
    x_k <- x == k
    p <- mean(x_k)
    mu <- mean(x[!x_k])
    var_x <- var(x[!x_k])

    theta <- if((var_x > mu)) mu^2/(var_x - mu) else 100
    z_star <- (p * as.numeric(x_k))/(p * as.numeric(x_k) + (1 - p) * dktnb(x, k, t = t, mu, theta, 0, return_log = FALSE))

    complete_log_likelihood <- .complete_log_likelihood(mu, theta, p, x, z_star, k, x_k, t)
    not_converged <- TRUE
    if (control$trace) {
        cat("Iteration :", 0, "\t Current log-likelihood:", complete_log_likelihood, "\t Absolute change in log-likelihood:", NA, "\n")
        cat("\t Parameters (Initial): ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
    }
    if (control$save_trace) {
        trace_list <- list()
        trace_list[[1]] <- data.frame(Iteration = 0, mu = mu, theta = theta, p = p, LogLikelihood = complete_log_likelihood,
                                      AbsoluteChangeLogLikelihood = NA)
    }

    i = 1
    while (not_converged) {
        ## E-step
        z_star <- (p * as.numeric(x_k))/(p * as.numeric(x_k) + (1 - p) * dktnb(x, k, t = t, mu, theta, 0, return_log = FALSE))

        ## M-step
        p <- sum(z_star * as.numeric(x_k)) / sum(1 - z_star*(1 - as.numeric(x_k)))

        pars <- optim(par = c(mu, theta), fn = .restricted_log_likelihood,
                              x = x, p = p, k = k, t = t,
                              method = "L-BFGS-B", lower = c(t, .Machine$double.eps), upper = c(Inf, Inf),
                              control = list(trace = F))$par

        mu <- pars[1]
        theta <- pars[2]


        complete_log_likelihood_old <- complete_log_likelihood
        complete_log_likelihood <- .complete_log_likelihood(mu, theta, p, x, z_star, k, x_k, t)
        change_log_likelihood <- abs(complete_log_likelihood - complete_log_likelihood_old)

        not_converged <- if (i < control$iteration_min) TRUE else ((change_log_likelihood > control$tolerance) && (i < control$iteration_max))

        if (control$trace) {
            cat("Iteration :", i, "\t Current log-likelihood:", complete_log_likelihood, "\t Absolute change in log-likelihood:", change_log_likelihood, "\n")
            cat("\t Parameters: ", "\t mu =", mu, "\t theta =", theta, "\t p =", p, "\n")
        }

        if (control$save_trace) {
            trace_list[[i + 1]] <- data.frame(Iteration = i, mu = mu, theta = theta, p = p, LogLikelihood = complete_log_likelihood,
                                           AbsoluteChangeLogLikelihood = change_log_likelihood)
        }

        i = i + 1
    }

    returned_trace <- NULL
    if (control$save_trace) {
        returned_trace <- trace_list %>% bind_rows() %>% as_tibble()
    }

    res <- list(n = length(x), k = k, t = t, mu = mu, theta = theta, p = p,
                logLikelihood = complete_log_likelihood, converged = !not_converged,
                trace = returned_trace)

    class(res) <- "ktnb"
    return(res)
}

#' Confidence envelopes of ktnb-object.
#'
#' @description Simulated confidence envelopes of the parameters estimated by the \link{ktnb_optimisation} function.
#'
#' @param ktnb_object An object of class 'ktnb'.
#' @param level The confidence level.
#' @param trace TRUE/FALSE, should a trace be shown?
#' @param number_of_simulations The number of simulations used to create the confidence envelopes.
#' @param plot_simulations TRUE/FALSE, should a density plot of the parameters be returned?
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' @export
simulate_confidence_envelopes <- function(ktnb_object, level = 0.95, trace = T, number_of_simulations = 30, plot_simulations = FALSE) {
    n = ktnb_object$n
    k = ktnb_object$k
    t = ktnb_object$t

    mu = ktnb_object$mu
    theta = ktnb_object$theta
    p = ktnb_object$p

    envelope_mu <- rep(NA, number_of_simulations)
    envelope_theta <- rep(NA, number_of_simulations)
    envelope_p <- rep(NA, number_of_simulations)

    if (trace) {
        pb <- progress_bar$new(format = paste0("Simulating confidence envelopes (size ", number_of_simulations, ")", ": [:bar] :percent Eta: :eta"),
                               total = number_of_simulations, clear = FALSE, width = 120)
    }

    for (i in 1:number_of_simulations) {
        if (trace)
            pb$tick()

        simulation_i <- rktnb(n, k, t, mu, theta, p)
        optimised_simulation_i <- ktnb_optimisation(simulation_i, k, t, ktnb_optimisation_control(trace = F))

        envelope_mu[i] <- optimised_simulation_i$mu
        envelope_theta[i] <- optimised_simulation_i$theta
        envelope_p[i] <- optimised_simulation_i$p
    }

    alpha <- (1 - level) / 2

    quantiles_mu <- quantile(sort(envelope_mu), probs = c(alpha, 1 - alpha))
    quantiles_theta <- quantile(sort(envelope_theta), probs = c(alpha, 1 - alpha))
    quantiles_p <- quantile(sort(envelope_p), probs = c(alpha, 1 - alpha))

    if (plot_simulations) {
        pp <- ggplot(data.frame(env = c(envelope_mu, envelope_theta, envelope_p),
                          para = rep(c("mu", "theta", "pi"), each = number_of_simulations)),
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

#' Plot parameter trace of ktnb-object
#'
#' @description A wrap function for ggplot2, plotting the Expectation-Maximisation parameter trace, which can returned by the \link{ktnb_optimisation} function.
#'
#' @param ktnb_object An object of the class 'ktnb'.
#'
#' @return A ggplot.
#' @export
ggplot.ktnb <- function(ktnb_object) {
    if (is.null(ktnb_object$trace))
        stop("The 'trace' tibble was not found. \nSet 'save_trace = TRUE' in the 'ktnb_optimisation_control' function.")

    trace_tibble <- ktnb_object$trace[-1, ] %>% rename("Log-likelihood" = LogLikelihood, "pi" = p)
    trace_tibble_long <- trace_tibble %>% select(-AbsoluteChangeLogLikelihood) %>% gather(para, val, -Iteration)

    pp <- ggplot(trace_tibble_long, aes(x = Iteration, y = val)) +
        geom_line(size = 1.05) + facet_wrap(~para, ncol = 2, nrow = 2, scales = "free_y", labeller = label_parsed) +
        scale_x_log10() + ylab("") +
        theme_bw(base_size = 13) + theme(legend.position = "none")

    return(pp)
}
