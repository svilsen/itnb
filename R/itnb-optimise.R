##
strip_terms <- function(formula) {
    attr_names <- names(attributes(formula))
    for (i in seq_along(attr_names)) {
        attr(formula, attr_names[i]) <- NULL
    }

    formula <- as.formula(formula)
    return(formula)
}

##
itnb_matrix <- function(X, y, i, t, link, control = list()) {
    ##
    if (is.null(i)) {
        i <- -1
    }

    if (!is.numeric(i)) {
        stop("'i' has to be numeric.")
    }
    else if (length(i) != 1) {
        warning("'i' needs to have length 1, only the first element is used.")
        i <- i[1]
    }
    i <- ceiling(i)

    ##
    if (is.null(t)) {
        t <- -1
    }

    if (!is.numeric(t)) {
        stop("'t' has to be numeric.")
    }
    else if (length(t) != 1) {
        warning("'t' needs to have length 1, only the first element is used.")
        t <- t[1]
    }
    t <- ceiling(t)

    ##
    yi <- as.numeric(y == i)

    ##
    p <- mean(yi, na.rm = TRUE)
    if (link == "identity") {
        m <- glm.nb(y[(y != i) & (y > 0), , drop = FALSE] ~ -1 + X[(y != i) & (y > 0), , drop = FALSE], link = "identity")
    }
    else if (link == "log") {
        m <- glm.nb(y[y != i, , drop = FALSE] ~ -1 + X[y != i, , drop = FALSE], link = "log")
    }
    else if (link == "sqrt") {
        m <- glm.nb(y[y != i, , drop = FALSE] ~ -1 + X[y != i, , drop = FALSE], link = "sqrt")
    }

    beta <- coef(m) |> unname() |> as.matrix()
    theta <- m[["theta"]]

    ##
    if ((i < 0) || (i <= t)) {
        control[["lambda"]] <- c(0, 0)

        res <- mle_itnb_cpp(
            X = X[y != i, , drop = FALSE], y = y[y != i, , drop = FALSE],
            beta_0 = beta, theta_0 = theta, p_0 = p,
            i = i, t = t, link = link,
            tolerance = control[["tolerance"]], lambda = control[["lambda"]],
            steps = control[["steps"]], exact = control[["exact"]],
            trace = control[["trace"]]
        )

        res[["trace"]] <- NA
    }
    else if ((t < 0) || (i > t)) {
        res <- em_itnb_cpp(
            X = X, y = y, yi = yi,
            beta_0 = beta, theta_0 = theta, p_0 = p,
            i = i, t = t, link = link,
            iteration_min = control[["iteration_min"]], iteration_max = control[["iteration_max"]],
            tolerance = control[["tolerance"]], lambda = control[["lambda"]],
            steps = control[["steps"]], exact = control[["exact"]],
            trace = control[["trace"]], save_trace = control[["save_trace"]]
        )
    }
    else {
        stop("How did we get here? Report with code 'itnb-opt-1'.")
    }

    return(res)
}

#' Control function for the \link{itnb} function.
#'
#' @description Creates a list of default options.
#'
#' @param method String: indicating the method used to optimise the parameters (currently only accepts \code{method = 'em'}).
#' @param trace Numeric (>= 0): showing a trace every \code{trace} number of iterations.
#' @param tolerance Numeric (> 0): Convergence tolerance.
#' @param lambda Numeric (> 0): Vector of penalisation constants (see details).
#' @param iteration_min Numeric (>= 0): The minimum number of allowed iterations.
#' @param iteration_max Numeric (>= \code{iteration_min}): The maximum number of allowed iterations.
#' @param steps Numeric (>= 0): The number of steps to use when approximating the integral needed for the derivative of the overdispersion.
#' @param exact TRUE/FALSE: should parameter optimisation use the exact gradient, or finite difference?
#' @param save_data TRUE/FALSE: should the data be stored in the return object?
#' @param save_trace TRUE/FALSE: should the entire trace be stored in the return object?
#'
#' @details The vector of penalisation constants, \code{lambda}, can be set using one or two elements. If a single element is provided, the penalisation
#' constant is fixed in every iteration of the EM-algorithm; however, if two elements are provided, the first is used as the initial penalisation constant, while
#' the second is used to scale the penalisation in every iteration, creating a sequence of penalisation constants.
#' NB: if the data does not contain an inflation point, or if the inflation point is smaller than the truncation point, then lambda should be set to zero.
#'
#' @return A list of default arguments for the \link{itnb} function.
#' @export
itnb_control <- function(method = "em", trace = 0, tolerance = 1e-6, lambda = c(0.1, 0.001), iteration_min = 5, iteration_max = 10000, steps = 100, exact = FALSE, save_data = TRUE, save_trace = TRUE) {
    if (!is.character(method)) {
        stop("'method' has to be a string.")
    }
    method <- tolower(method)

    if (!(method %in% c("em", "expectation-maximisation", "expectation-maximization"))) {
        stop("'method' has be set to 'em' -- this option is included for future proofing, and should not be used at present.")
    }

    ##
    if (!is.numeric(trace)) {
        stop("'trace' has to be numeric.")
    }
    trace <- ceiling(trace)

    ##
    if (!is.numeric(tolerance)) {
        stop("'tolerance' has to be numeric.")
    }

    if (tolerance < .Machine$double.eps) {
        stop("'tolerance' has to be bigger than two times the machine epsilon.")
    }

    ##
    if (!any(is.numeric(lambda))) {
        stop("'lambda' has to be numeric.")
    }

    if (length(lambda) < 2) {
        lambda <- c(lambda, 1)
    }
    else if (length(lambda) > 2) {
        lambda <- lambda[seq_len(2)]
    }

    ##
    if (!is.numeric(iteration_min)) {
        stop("'iteration_min' has to be numeric.")
    }
    iteration_min <- ceiling(iteration_min)

    if (!is.numeric(iteration_max)) {
        stop("'iteration_min' has to be numeric.")
    }
    iteration_max <- ceiling(iteration_max)

    if (iteration_min > iteration_max) {
        stop("'iteration_max' has be larger than 'iteration_min'.")
    }

    ##
    if (!is.numeric(steps)) {
        stop("'steps' has to be numeric.")
    }
    steps <- ceiling(steps)

    ##
    if (!is.logical(exact)) {
        stop("'exact' has to be logical.")
    }

    ##
    if (!is.logical(save_trace)) {
        stop("'save_trace' has to be logical.")
    }

    ##
    if (!is.logical(save_data)) {
        stop("'save_data' has to be logical.")
    }

    ##
    res <- list(
        method = method, trace = trace, tolerance = tolerance, lambda = lambda,
        iteration_max = iteration_max, iteration_min = iteration_min,
        steps = steps, exact = exact,
        save_data = save_data, save_trace = save_trace
    )

    ##
    return(res)
}

#' @title itnb parameter optimisation
#'
#' @description Parameter optimisation for data generated by a r.v. following an itnb distribution using the ENM algorithm.
#'
#' @param formula A \link{formula}-object.
#' @param data A data-set, supplied as a \link{data.frame} (or \link[tibble]{tibble}).
#' @param i Numeric (i >= 0): The inflation point.
#' @param t Numeric (t >= 0): The truncation point.
#' @param link String indicating the link function, currently limited to 'log', 'sqrt', or 'identity'. NB: currently only 'identity' is implemented.
#' @param control List: A control object, see \link{itnb_control} for details.
#'
#' @return An object of class \link{itnb-object}.
#'
#' @export
itnb <- function(formula, data = NULL, i = NULL, t = NULL, link = "identity", control = list()) {
    UseMethod("itnb")
}

#' @rdname itnb
#' @method itnb formula
#'
#' @example inst/examples/itnb_example.R
#'
#' @export
itnb.formula <- function(formula, data = NULL, i = NULL, t = NULL, link = "identity", control = list()) {
    ##
    control <- do.call(itnb_control, control)

    ## Checks for 'data'
    keep_formula <- TRUE
    if (is.null(data)) {
        data <- tryCatch(
            expr = {
                as.data.frame(as.matrix(model.frame(formula)))
            },
            error = function(e) {
                stop("'data' needs to be supplied when using 'formula'.")
            }
        )

        x_name <- paste0(attr(terms(formula), "term.labels"), ".")
        colnames(data) <- paste0("V", gsub(x_name, "", colnames(data)))
        colnames(data)[1] <- "y"

        formula <- paste(colnames(data)[1], "~", paste(colnames(data)[seq_along(colnames(data))[-1]], collapse = " + "))
        formula <- as.formula(formula)
        keep_formula <- FALSE
    }

    # Re-capture feature names when '.' is used in formula interface
    formula <- terms(formula, data = data)
    formula <- strip_terms(formula)

    ##
    if ((is.null(i) | (i < 0)) & (is.null(t) | (t < 0))) {
        ##
        warning("When both 'i' and 't' are 'NULL' the problem reduces to a regular negative binomial regression; see the 'glm.nb' function from the 'MASS' package.")
    }

    ##
    #
    X <- model.matrix(formula, data)

    #
    y <- model.response(model.frame(formula, data))
    y <- as.matrix(y, nrow = nrow(data))

    ##
    #
    res <- itnb_matrix(X = X, y = y, i = i, t = t, link = link, control = control)

    ##
    #
    if (keep_formula) {
        res[["formula"]] <- formula
    }
    else {
        res[["formula"]] <- NA
    }

    #
    if (control[["save_trace"]] & !all(is.na(res[["trace"]]))) {
        res[["trace"]] <- do.call("cbind", res[["trace"]]) |> as.data.frame()
        names(res[["trace"]])[grep("V", names(res[["trace"]]))] <- colnames(X)
        res[["trace"]] <- cbind(Iteration = seq_len(nrow(res[["trace"]])) - 1, res[["trace"]])
    }
    else {
        res[["trace"]] <- NA
    }

    #
    if (control$save_data) {
        res[["data"]] <- list(X = X, y = y)
    }
    else {
        res[["data"]] <- NA
    }

    #
    if (res[["converged"]]) {
        res[["flag"]] <- NA
    }

    ##
    class(res) <- "itnb"
    return(res)
}

#' Summary of an inflated and truncated negative binomial regression model fit
#'
#' @description Summary function for an \link{itnb-object} (an object of class \code{itnb}).
#'
#' @param object An \link{itnb-object}.
#' @param ... Additional arguments (see details).
#'
#' @details ...
#'
#' @return A function summary.
#' @export
summary.itnb <- function(object, ...) {

}

#' @export
print.summary.itnb <- function(x, ...) {

}

#' Extract model coefficients
#'
#' @description A function which extracts model coefficients from objects of class \code{itnb}.
#'
#' @param object An \link{itnb-object}.
#' @param ... Additional arguments (see details).
#'
#' @details The only additional argument used by the function is \code{par} used to specify which parameter should be returned by the function, i.e. it takes the values \code{"beta"}, \code{"theta"}, or \code{"p"}.
#'
#' @return If \code{par} is left as \code{NULL} a list of all parameters will be returned, otherwise the function returns the parameter specified by \code{par}.
#' @export
coef.itnb <- function(object, ...) {
    dots <- list(...)

    if (is.null(dots[["par"]])) {
        r <- object[c("beta", "theta", "p")]
    }
    else if (dots[["par"]] == "beta") {
        r <- object[["beta"]]
    }
    else if (dots[["par"]] == "theta") {
        r <- object[["theta"]]
    }
    else if (dots[["par"]] == "p") {
        r <- object[["p"]]
    }
    else {
        stop("'pars' only takes the values 'mu', 'theta', or 'p'.")
    }

    return(r)
}

#' Plot parameter trace of an \link{itnb-object}
#'
#' @description A function plotting the EM parameter trace returned by the \link{itnb} function. Note: the function can only be used if the argument \code{save_trace} in the \link{itnb_control} function was set to \code{TRUE}.
#'
#' @param x An \link{itnb-object}.
#' @param ... Additional arguments passed to the \link[graphics]{plot} function.
#'
#' @export
plot.itnb <- function(x, ...) {
    dots <- list(...)
    if (all(is.na(x[["trace"]]))) {
        stop("The 'trace' data.frame was not found. Set 'save_trace = TRUE' in the 'itnb_control' function, and re-run the optimisation routine.")
    }

    itnb_trace <- x[["trace"]]
    if (!is.null(dots[["log"]])) {
        itnb_trace <- itnb_trace[-1, ]
    }

    par(mfrow = c(2, 2))
    plot(itnb_trace$Iteration, itnb_trace$LogLikelihood, type = "l", xlab = "Iteration", ylab = "Log-likelihood", ...)

    plot(itnb_trace$Iteration, itnb_trace$mu, type = "l", xlab = "Iteration", ylab = bquote(mu), ...)

    plot(itnb_trace$Iteration, itnb_trace$theta, type = "l", xlab = "Iteration", ylab = bquote(theta), ...)

    plot(itnb_trace$Iteration, itnb_trace$p, type = "l", xlab = "Iteration", ylab = bquote(pi), ...)
    par(mfrow = c(1, 1))

    return(invisible(NULL))
}
