#' @keywords internal
"_PACKAGE"

#' @title ITNB: Inflated and truncated negative binomial regression
#'
#' @author Søren B. Vilsen <svilsen@math.aau.dk>
#'
#' @importFrom Rcpp evalCpp
#' @importFrom graphics 'par' 'plot' 'hist' 'abline' 'arrows'
#' @importFrom stats 'terms' 'as.formula' 'model.frame' 'model.matrix' 'model.response' 'quantile' "coef" "delete.response" "fitted.values" "predict" "residuals"
#'
#' @useDynLib itnb
#'
#' @name itnb-package
#'
#' @rdname itnb-package
NULL

#' @title An itnb-object
#'
#' @description An itnb-object is a list of class \code{itnb} containing the following:
#' \describe{
#'     \item{\code{formula}}{The \link{formula} supplied to the model.}
#'     \item{\code{data}}{If \code{save_trace = TRUE} in the control object, it contains the data used to optimise the parameters of the model, otherwise it is \code{NA}.}
#'     \item{\code{i}}{The specified inflation point.}
#'     \item{\code{t}}{The specified truncation point.}
#'     \item{\code{link}}{The specified link function.}
#'     \item{\code{loglikelihood}}{The loglikelihood of the optimised parameters.}
#'     \item{\code{beta}}{The optimised regression coefficients.}
#'     \item{\code{theta}}{The optimised overdispersion.}
#'     \item{\code{p}}{The optimised inflation proportion.}
#'     \item{\code{trace}}{If \code{save_trace = TRUE} in the control object, it contains a \link{data.frame} of the trace produced by the optimisation routine, otherwise it is \code{NA}.}
#'     \item{\code{converged}}{Boolean indicating whether the optimisation routine did or did not converge.}
#'     \item{\code{flag}}{If the optimisation routine did not converge, a string indicating the potential point of failure is returned, otherwise it is \code{NA}.}
#' }
#'
#' @name itnb-object
#'
#' @rdname itnb-object
NULL

#' @title A CI-object
#'
#' @description A ci-object is a list of class \code{ci} containing the following:
#' \describe{
#'     \item{\code{ci}}{A list of elements corresponding to the mean, overdispersion, and inflation proportion. If a confidence level is specified it contains the quantiles for each parameter, otherwise all estimated parameters are returned.}
#'     \item{\code{parametric}}{Boolean indicating whether the bootstrap was parametric or non-parametric.}
#'     \item{\code{level}}{The specified confidence level; if a confidence level is not supplied it is \code{NA}.}
#'     \item{\code{nr_simulations}}{The number of bootstrap samples used to create the simulated confidence envelopes.}
#' }
#'
#' @name itnb.ci-object
#'
#' @rdname itnb.ci-object
NULL

