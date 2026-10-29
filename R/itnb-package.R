#' @keywords internal
"_PACKAGE"

#' @title ITNB: Inflated and truncated negative binomial regression
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats 'terms' 'as.formula' 'model.frame' 'model.matrix' 'model.response' 'quantile' 'coef' 'delete.response' 'fitted.values' 'predict' 'residuals' 'printCoefmat' 'sd' 'pnorm' 'qnorm' 'qt' 'pchisq' 'qchisq' 'dchisq' 'density'
#' @importFrom graphics 'abline' 'box' 'hist' 'points' 'legend' 'par'
#' @importFrom grDevices 'rgb'
#'
#' @useDynLib itnb
#'
#' @name itnb-package
#' @rdname itnb-package
NULL

#' @title An itnb-object
#'
#' @description An itnb-object is a list of class \code{itnb} containing the following:
#' \describe{
#'     \item{\code{formula}}{\link{formula}: The formula of the specified the model.}
#'     \item{\code{data}}{\link{data.frame}: The data-set used to estimate the parameters. NB: If '\code{save_trace = FALSE}' in the control object (see \link{itnb_control}) it will return '\code{NA}'.}
#'     \item{\code{i}}{Numeric (>= 0): The specified inflation point.}
#'     \item{\code{t}}{Numeric (>= 0): The specified truncation point.}
#'     \item{\code{link}}{String: The specified link function.}
#'     \item{\code{loglikelihood}}{Numeric: The loglikelihood of the estimated model parameters.}
#'     \item{\code{beta}}{\link{matrix}: The estimated regression coefficients.}
#'     \item{\code{alpha}}{Numeric: The estimated overdispersion (correpsonds to 1 / theta in the standard parameterisation).}
#'     \item{\code{p}}{Numeric: The estimated inflation proportion.}
#'     \item{\code{vcov}}{\link{matrix}: The variance-covariance matrix of the estimated model parameters.}
#'     \item{\code{logtheta}}{Numeric: The standard error of the logarithm of theta (1 / alpha).}
#'     \item{\code{trace}}{\link{data.frame}: The trace of the estimated parameters. NB: If '\code{save_trace = FALSE}' it will return '\code{NA}'.}
#'     \item{\code{converged}}{TRUE/FALSE: Did the optimisation converge?}
#'     \item{\code{iterations}}{Numeric: The number of iterations of the expectation-maximisation algorithm.}
#'     \item{\code{flag}}{String: Indicating the potential point of failure if the algorithm did not converge. NB: if the algorithm did converge it will return '\code{NA}'.}
#' }
#'
#' @name itnb-object
#' @rdname itnb-object
NULL

#' @title A summary.itnb-object
#'
#' @description A summary.itnb-object is a list of class \code{summary.itnb} containing the following:
#' \describe{
#'     \item{\code{formula}}{\link{formula}: The formula of the specified the model.}
#'     \item{\code{i}}{Numeric (>= 0): The specified inflation point.}
#'     \item{\code{t}}{Numeric (>= 0): The specified truncation point.}
#'     \item{\code{link}}{String: The specified link function.}
#'     \item{\code{residuals}}{Vector: The residuals of the fitted model.}
#'     \item{\code{type}}{String: The type of residuals calculated.}
#'     \item{\code{level}}{Numeric (> 0): The confidence level.}
#'     \item{\code{coefficients}}{\link{matrix}: The estimates, standard errors, confidence intervals, test-statistic, and p-value of the regression coefficients.}
#'     \item{\code{overdispersion}}{\link{matrix}: The estimates, standard errors, confidence intervals, test-statistic, and p-value of the overdispersion.}
#'     \item{\code{inflation}}{\link{matrix}: The estimates, standard errors, confidence intervals, test-statistic, and p-value of the inflation proportion.}
#'     \item{\code{bootstrap}}{\link{itnb.ci-object}: If bootstrap was used to determine the standard error, confidence intervals, test-statistic, and p-value, then it contains the resulting \link{itnb.ci-object}. NB: If '\code{bootstrap = FALSE}' it will return '\code{NA}'.}
#'     \item{\code{iterations}}{Numeric: The number of iterations of the expectation-maximisation algorithm.}
#' }
#'
#' @name summary.itnb-object
#' @rdname summary.itnb-object
NULL

#' @title An itnb.ci-object
#'
#' @description A itnb.ci-object is a list of class \code{itnb.ci} containing the following:
#' \describe{
#'     \item{\code{ci}}{List: Boostrap confidence intervals of regression coefficients, overdispersion, and inflation proportion. If a confidence level is specified it contains the quantiles for each parameter, otherwise the estimated parameters are returned for all \code{B} bootstrap samples.}
#'     \item{\code{parametric}}{TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?}
#'     \item{\code{level}}{Numeric (> 0): The confidence level. If left as \code{NULL} all parametric bootstrap simulations are returned.}
#'     \item{\code{B}}{Numeric (> 0): The number of simulations used to create the confidence envelopes.}
#' }
#'
#' @name itnb.ci-object
#' @rdname itnb.ci-object
NULL

#' @title An lrtest.itnb-object
#'
#' @description An lrtest.itnb-object is a list of class \code{lrtest.itnb} containing the following:
#' \describe{
#'     \item{\code{lr}}{Numeric: The observed likelihood ratio test statistic.}
#'     \item{\code{df}}{Numeric: The degrees of freedom of the associated test statistic (i.e., the difference in parameters between the models).}
#'     \item{\code{critval}}{Numeric: The critical value of the test. If \code{bootstrap = TRUE} it returns the (1 - \code{level})-quantile of the parametric bootstrap distribution, otherwise it returns \code{qchisq(1.0 - level, df = df)}.}
#'     \item{\code{pval}}{Numeric: The p-value of the test. If \code{bootstrap = TRUE} it returns the p-value given the parametric bootstrap distribution, otherwise it returns \code{pchisq(lr, df = df, lower.tail = FALSE)}.}
#'     \item{\code{type}}{String: The likelihood ratio test performed, either \code{"overdispersion"}, \code{"inflation"}, or \code{"both"}.}
#'     \item{\code{level}}{Numeric (> 0): The significance level.}
#'     \item{\code{bootstrap}}{If '\code{bootstrap = FALSE}' it will return '\code{NA}', otherwise it returns a list with the elements:
#'          \describe{
#'              \item{\code{parametric}}{TRUE/FALSE: should the envelopes be simulated using the parametric bootstrap?.}
#'              \item{\code{B}}{Numeric (> 0): The number of simulations used to create the confidence envelopes.}
#'              \item{\code{lr}}{Vector: The parametric bootstrap distribution of the likelihood ratio test statistic.}
#'          }
#'     }
#' }
#'
#' @name lrtest.itnb-object
#' @rdname lrtest.itnb-object
NULL
