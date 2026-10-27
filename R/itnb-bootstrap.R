#' @title Confidence intervals of itnb-object.
#'
#' @description Bootstrap confidence intervals of the parameters estimated by the \link{itnb} function.
#'
#' @param object \link{itnb-object}.
#' @param parm String: Specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing confidence intervals will be made for all parameters. NB: currently not used.
#' @param level Numeric (> 0): The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.
#' @param ... Additional arguments (see details).
#'
#' @details The additional arguments are:
#' \describe{
#'      \item{\code{B}}{Numeric (> 0): The number of simulations used to create the confidence envelopes.}
#'      \item{\code{parametric}}{TRUE/FALSE: Should the envelopes be simulated using the parametric bootstrap?}
#'      \item{\code{trace}}{Numeric (>= 0): Showing a trace every \code{trace} number of iterations.}
#' }
#'
#' @example inst/examples/simulation_itnb_example.R
#'
#' @method confint itnb
#'
#' @return If \code{level = NULL} a matrix with bootstrap simulations, otherwise a matrix of lower and upper confidence limits for each parameter.
#' @export
confint.itnb <- function(object, parm, level = 0.95, ...) {
    ##
    #
    dots <- list(...)

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
    if (is.null(dots[["B"]])) {
        B <- 200
    }
    else {
        B <- dots[["B"]]
    }

    if (!is.numeric(B)) {
        stop("'B' has to be numeric.")
    }
    else if (B < 1) {
        stop("'B' has to be > 0.")
    }

    B <- ceiling(B)

    ##
    if (is.null(dots[["parametric"]])) {
        parametric <- FALSE
    }
    else {
        parametric <- dots[["parametric"]]
    }

    if (!is.logical(parametric)) {
        stop("'parametric' has to be logical.")
    }

    ##
    if (is.null(dots[["trace"]])) {
        trace <- 0
    }
    else {
        trace <- dots[["trace"]]
    }

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

    control <- object[["control"]]
    control[["trace"]] <- 0
    control[["save_data"]] <- FALSE

    if (parametric) {
        mu <- X %*% beta

        if (link == "sqrt") {
            mu <- mu * mu
        }
        else if (link == "log") {
            mu <- exp(mu)
        }
        else if (link == "reciprocal") {
            mu <- 1 / mu
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
            "beta" = ci_list[["beta"]] |> apply(2, quantile, probs = c(sig_level, 1 - sig_level), na.rm = TRUE) |> t(),
            "alpha" = quantile(ci_list[["alpha"]], probs = c(sig_level, 1 - sig_level), na.rm = TRUE),
            "p" = quantile(ci_list[["p"]], probs = c(sig_level, 1 - sig_level), na.rm = TRUE)
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


#' @title Likelihood ratio tests
#'
#' @description Likelihood ratio tests for inflation and/or overdispersion in inflated and truncated negative binomial regression models.
#'
#' @param object \link{itnb-object}.
#' @param type String: The likelihood ratio test performed, either \code{"overdispersion"}, \code{"inflation"}, or \code{"both"}.
#' @param level Numeric (> 0): The significance level.
#' @param ... Additional parameters (see details).
#'
#' @details The additional arguments are:
#' \describe{
#'      \item{\code{control}}{List: A control object (see \link{itnb_control}).}
#'      \item{\code{bootstrap}}{TRUE/FALSE: Should the likelihood ratio test be performed using paramteric bootrapping?}
#'      \item{\code{B}}{Numeric (> 0): The number of simulations used to create the confidence envelopes.}
#'      \item{\code{parametric}}{TRUE/FALSE: Should the envelopes be simulated using the parametric bootstrap? NB: Only accepts \code{TRUE}.}
#'      \item{\code{trace}}{Numeric (>= 0): Showing a trace every \code{trace} number of iterations.}
#' }
#'
#' @return An \link{lrtest.itnb-object}.
#'
#' @export
lrtest <- function(object, type = "overdispersion", level = 0.05, ...) {
    UseMethod("lrtest")
}

#' @title Likelihood ratio tests
#'
#' @description Likelihood ratio tests for inflation and/or overdispersion in inflated and truncated negative binomial regression models.
#'
#' @param object \link{itnb-object}.
#' @param type String: The likelihood ratio test performed, either \code{"overdispersion"}, \code{"inflation"}, or \code{"both"}.
#' @param level Numeric (> 0): The significance level.
#' @param ... Additional parameters (see details).
#'
#' @details The additional arguments are:
#' \describe{
#'      \item{\code{control}}{List: A control object (see \link{itnb_control}).}
#'      \item{\code{bootstrap}}{TRUE/FALSE: Should the likelihood ratio test be performed using paramteric bootrapping?}
#'      \item{\code{B}}{Numeric (> 0): The number of simulations used to create the confidence envelopes.}
#'      \item{\code{parametric}}{TRUE/FALSE: Should the envelopes be simulated using the parametric bootstrap? NB: Only accepts \code{TRUE}.}
#'      \item{\code{trace}}{Numeric (>= 0): Showing a trace every \code{trace} number of iterations.}
#' }
#'
#' @return An \link{lrtest.itnb-object}.
#'
#' @method lrtest itnb
#' @rdname lrtest
#'
#' @export
lrtest.itnb <- function(object, type = "overdispersion", level = 0.05, ...) {
    ##
    #
    dots <- list(...)

    ##
    #
    if (is.null(dots[["control"]])) {
        control <- do.call(itnb_control, object[["control"]])
    }
    else {
        control <- do.call(itnb_control, dots[["control"]])
    }

    ##
    #
    if (is.null(dots[["bootstrap"]])) {
        bootstrap <- TRUE
    }
    else {
        bootstrap <- dots[["bootstrap"]]
        if (!is.logical(bootstrap)) {
            stop("'bootstrap' needs to be boolean.")
        }
    }

    #
    if (bootstrap) {
        if (is.null(dots[["parametric"]])) {
            parametric <- TRUE
        }
        else {
            parametric <- dots[["parametric"]]
            if (!is.logical(parametric)) {
                stop("'parametric' needs to be boolean.")
            }

            if (!parametric) {
                stop("'parametric' needs to be 'TRUE'.")
            }
        }

        if (is.null(dots[["B"]])) {
            B <- 199
        }
        else {
            B <- dots[["B"]]
            if (!is.numeric(B)) {
                stop("'B' needs to be boolean.")
            }
            else if (B - floor(B) > 1e-8) {
                B <- B |> round()
            }
        }
    }
    else {
        parametric <- NA
        B <- NA
    }

    ##
    #
    trace <- dots[["trace"]]
    if (is.null(trace)) {
        trace <- 0
    }

    ##
    #
    X <- object[["data"]][["X"]]
    y <- object[["data"]][["y"]]

    #
    N <- length(y)
    i <- object[["i"]]
    t <- object[["t"]]

    link <- object[["link"]]

    ##
    #
    control_r <- control
    i_r <- i
    alpha_r <- object[["alpha"]]
    if (type %in% c("o", "overdispersion")) {
        type <- "overdispersion"
        control_r[["overdispersion"]] <- FALSE
        alpha_r <- 0.0
        df <- 1
    }
    else if (type %in% c("i", "inflation")) {
        type <- "inflation"
        i_r <- -1
        df <- 1
    }
    else if (type %in% c("b", "both")) {
        type <- "both"
        control_r[["overdispersion"]] <- FALSE
        alpha_r <- 0.0
        i_r <- -1
        df <- 2
    }
    else {
        stop("'type' only takes the values 'overdispersion', 'inflation', or 'both'.")
    }

    ##
    #
    loglike_o <- object[["loglikelihood"]]
    loglike_r <- itnb_matrix(X = X, y = y, i = i_r, t = t, link = link, control = control_r)[["loglikelihood"]]

    #
    d <- -2.0 * (loglike_r - loglike_o)

    ##
    #
    if (bootstrap) {
        ##
        #
        beta <- object[["beta"]]
        if (parametric) {
            mu <- X %*% beta

            if (link == "sqrt") {
                mu <- mu * mu
            }
            else if (link == "log") {
                mu <- exp(mu)
            }
            else if (link == "reciprocal") {
                mu <- 1 / mu
            }
        }

        #
        p <- object[["p"]]

        ##
        #
        d_boot <- rep(NA, B)
        for (b in seq_len(B)) {
            if ((trace > 0) && ((b == 1) || ((b %% trace) == 0) || (b == B))) {
                cat("Iteration:", b, "/", B, "\n")
            }

            ##
            #
            y_b <- ritnb(n = N, mu = mu, alpha = alpha_r, p = p, i = i_r, t = t) |>
                matrix(ncol = 1)

            ##
            #
            loglike_ob <- itnb_matrix(X = X, y = y_b, i = i, t = t, link = link, control = control)[["loglikelihood"]]
            loglike_rb <- itnb_matrix(X = X, y = y_b, i = i_r, t = t, link = link, control = control_r)[["loglikelihood"]]

            d_boot[b] <- 2.0 * abs(loglike_ob - loglike_rb)
        }

        ##
        #
        crit_val <- quantile(d_boot, probs = 1.0 - level, na.rm = TRUE) |> unname()
        p_val <- (sum(d < d_boot) + 1) / (B + 1)
    }
    else {
        ##
        #
        crit_val <- qchisq(1.0 - level, df = df)
        p_val <- pchisq(d, df = df, lower.tail = FALSE)
    }

    res <- list(
        lr = d,
        df = df,
        critval = crit_val,
        pval = p_val,
        type = type,
        level = level,
        bootstrap = if (bootstrap) list(parametric = parametric, B = B, lr = d_boot) else NA
    )

    class(res) <- "lrtest.itnb"
    return(res)
}

#' @title Print Likelihood Ratio Test
#'
#' @description Print a \link{lrtest.itnb-object}, i.e., the likelihood ratio test results of an \link{itnb-object}.
#'
#' @param x \link{lrtest.itnb-object}
#' @param ... Additional arguments (see details).
#'
#' @details The additional argument is \code{digits}, i.e., the number of digits to show in the printed results.
#'
#' @return An invisible \code{NULL}.
#'
#' @method print lrtest.itnb
#'
#' @export
print.lrtest.itnb <- function(x, ...) {
    dots <- list(...)
    if (is.null(dots[["digits"]])) {
        digits <- max(3, getOption("digits") - 2)
    }
    else {
        digits <- dots[["digits"]]
    }

    cat("\n")

    ##
    #
    if (x[["type"]] %in% c("overdispersion")) {
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

    cat(" ", paste0("Critical value (level = ", round(x[["level"]], digits = digits), "): "), round(x[["critval"]], digits = digits), "\n")
    cat(" ", "Test statistic: ", round(x[["lr"]], digits = digits), "\n")
    cat(" ", "P-value:", format.pval(x[["pval"]], digits = digits), "\n")
    cat("\n")

    if (!all(is.na(x[["bootstrap"]]))) {
        cat("(Based on", x[["bootstrap"]][["B"]], ifelse(x[["bootstrap"]][["parametric"]], "parametric", "non-parametric"), "bootstrap samples => p-value lower bound of", paste0(format(1 / (x[["bootstrap"]][["B"]] + 1), digits = digits), ")."), "\n")
        cat("\n")
    }

    return(invisible(NULL))
}

#' @title Likelihood Ratio Histogram
#'
#' @description Histogram of the likelihood ratio sample distribution found by parametric bootstrap.
#'
#' @param x \link{lrtest.itnb-object}.
#' @param ... Additional arguments (see \link{hist}).
#'
#' @return An invisible \code{NULL}.
#'
#' @method hist lrtest.itnb
#'
#' @export
hist.lrtest.itnb <- function(x, ...) {
    #
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    #
    dots <- list(...)

    #
    if (all(is.na(x[["bootstrap"]]))) {
        stop("'lrtest' has to be run with argument 'bootstrap = TRUE'.")
    }

    #
    if (is.null(dots[["xlog"]])) {
        dots[["xlog"]] <- FALSE
    }

    #
    dots[["x"]] <- x$bootstrap$lr
    critval <- x$critval
    lr <- x$lr

    if (dots[["xlog"]]) {
        dots[["x"]] <- log(dots[["x"]])
        critval <- log(critval)
        lr <- log(lr)
    }

    #
    if (is.null(dots[["prob"]])) {
        dots[["prob"]] <- TRUE
    }

    #
    if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- paste("Likelihood ratio")
        if (dots[["xlog"]]) {
            dots[["xlab"]] <- paste0("log(", dots[["xlab"]], ")")
        }
    }

    #
    if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- ifelse(dots[["prob"]], "Density", "Frequency")
    }

    #
    if (is.null(dots[["col"]])) {
        dots[["col"]] <- rgb(28 / 256, 134 / 256, 238 / 256, 0.6)
    }

    #
    if (is.null(dots[["main"]])) {
        dots[["main"]] <- paste("Likelihood ratio sample distribution \n (Based on", x[["bootstrap"]][["B"]], ifelse(x[["bootstrap"]][["parametric"]], "parametric", "non-parametric"), "bootstrap samples)")
    }

    #
    if (is.null(dots[["breaks"]])) {
        dots[["breaks"]] <- "fd"
    }

    #
    if (is.null(dots[["probs"]])) {
        probs <- 0.975
    }
    else {
        probs <- dots[["probs"]]
        dots[["probs"]] <- NULL
    }

    #
    if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(0, max(c(quantile(dots[["x"]], probs = probs), critval, lr)) + 1)
        if (dots[["xlog"]]) {
            dots[["xlim"]] <- c(min(dots[["x"]]) - 1, max(dots[["x"]]) + 1)
        }
    }

    #
    if (is.null(dots[["position"]])) {
        position <- "topright"
    }
    else {
        position <- dots[["position"]]
        dots[["position"]] <- NULL
    }

    #
    if (is.null(dots[["base"]])) {
        base <- FALSE
    }
    else {
        base <- dots[["base"]]
        dots[["base"]] <- NULL
    }

    #
    if (!base) {
        par(
            pty = "m",
            bg = "white",
            fg = "black",
            col.axis = "black",
            col.lab = "black",
            bty = "l",
            las = 1,
            xpd = FALSE
        )
    }

    do.call(hist, dots)
    abline(v = critval, col = "black", lwd = 2)
    abline(v = lr, col = "red", lwd = 2)

    #
    if (!base) {
        box(col = "grey90", bty = "o")
        box(col = "black", bty = "l")
    }

    legend(
        position,
        inset = 0.02,
        legend = c("Critical value", "Observed value"),
        lty = 1,
        lwd = 2,
        col = c("black", "red"),
        pt.bg = c("white", "white"),
        bty = "o",
        box.col = "grey90",
        bg = "white"
    )

    #
    return(invisible(NULL))
}

#' @title Likelihood Ratio Plot
#'
#' @description Plot of the likelihood ratio sample distribution based on either the asymptotic Chi-squared distribution or bootstrap distribution.
#'
#' @param x \link{lrtest.itnb-object}.
#' @param ... Additional arguments (see \link{plot}).
#'
#' @return An invisible \code{NULL}.
#'
#' @method plot lrtest.itnb
#'
#' @export
plot.lrtest.itnb <- function(x, ...) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    #
    dots <- list(...)

    #
    if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- paste("Likelihood ratio")
    }

    #
    if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "Density"
    }

    #
    if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
    }

    #
    if (all(is.na(x[["bootstrap"]]))) {    #
        if (is.null(dots[["main"]])) {
            dots[["main"]] <- paste("Likelihood ratio sample distribution \n (Based on Chi-squared distribution with", x$df, "degree(s)-of-freedom)")
        }

        m <- max(x$critval, x$lr) + 1
        dots[["x"]] <- seq(0, m, length.out = 10 * m)
        dots[["y"]] <- dchisq(dots[["x"]], df = x$df)
    }
    else {
        if (is.null(dots[["main"]])) {
            dots[["main"]] <- paste("Likelihood ratio sample distribution \n (Based on", x[["bootstrap"]][["B"]], ifelse(x[["bootstrap"]][["parametric"]], "parametric", "non-parametric"), "bootstrap samples)")
        }

        if (is.null(dots[["bw"]])) {
            bw <- "nrd0"
        }
        else {
            bw <- dots[["bw"]]
            dots[["bw"]] <- NULL
        }

        if (is.null(dots[["kernel"]])) {
            kernel <- "gaussian"
        }
        else {
            kernel <- dots[["kernel"]]
            dots[["kernel"]] <- NULL
        }

        if (is.null(dots[["digits"]])) {
            digits <- 5
        }
        else {
            digits <- dots[["digits"]]
            dots[["digits"]] <- NULL
        }

        if (is.null(dots[["probs"]])) {
            probs <- 0.975
        }
        else {
            probs <- dots[["probs"]]
            dots[["probs"]] <- NULL
        }

        m <- max(c(x$critval, x$lr, quantile(x$bootstrap$lr, probs = probs))) + 1
        d <- density(x$bootstrap$lr, from = 0, to = m, bw = bw, kernel = kernel)

        dots[["xlab"]] <- paste0(dots[["xlab"]], "\n Band-width: ", round(d$bw, digits))
        dots[["x"]] <- d$x
        dots[["y"]] <- d$y
    }

    #
    if (is.null(dots[["position"]])) {
        position <- "topright"
    }
    else {
        position <- dots[["position"]]
        dots[["position"]] <- NULL
    }

    #
    if (is.null(dots[["base"]])) {
        base <- FALSE
    }
    else {
        base <- dots[["base"]]
        dots[["base"]] <- NULL
    }

    #
    if (!base) {
        par(
            pty = "m",
            bg = "white",
            fg = "black",
            col.axis = "black",
            col.lab = "black",
            bty = "l",
            las = 1,
            xpd = FALSE
        )
    }

    do.call(plot, dots)
    abline(v = x$critval, col = "black", lwd = 2)
    abline(v = x$lr, col = "red", lwd = 2)

    #
    if (!base) {
        box(col = "grey90", bty = "o")
        box(col = "black", bty = "l")
    }

    legend(
        position,
        inset = 0.02,
        legend = c("Critical value", "Observed value"),
        lty = 1,
        lwd = 2,
        col = c("black", "red"),
        pt.bg = c("white", "white"),
        bty = "o",
        box.col = "grey90",
        bg = "white"
    )

    #
    return(invisible(NULL))
}
