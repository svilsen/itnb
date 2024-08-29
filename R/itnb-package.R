#' @keywords internal
"_PACKAGE"

#' @title ITNB: Inflated and truncated negative binomial distribution
#'
#' @author Søren B. Vilsen <svilsen@math.aau.dk>
#'
#' @importFrom stats 'optim' 'pbeta' 'quantile' 'rnbinom' 'var'
#' @importFrom progress 'progress_bar'
#' @importFrom graphics 'par' 'hist' 'plot'
#'
#' @name itnb-package
#'
#' @rdname itnb-package
NULL

#' @title An itnb-object
#'
#' @description An itnb-object is a list containing the following:
#' \describe{
#'     \item{\code{n}}{Integer: The number of observed counts.}
#'     \item{\code{x}}{Vector: The observed counts.}
#'     \item{\code{i}}{Integer: The inflation point.}
#'     \item{\code{t}}{Integer: The truncation point.}
#'     \item{\code{mu}}{Numeric: The expected value.}
#'     \item{\code{theta}}{Numeric: The overdispersion.}
#'     \item{\code{p}}{Numeric: The inflation proportion.}
#'     \item{\code{log_likelihood}}{Numeric: The complete log-likelihood of the optimised parameters.}
#'     \item{\code{converged}}{TRUE/FALSE: Did the EM-algorithm converge?}
#'     \item{\code{trace}}{A \link{data.frame} of the trace, if the argument \code{save_trace} in \link{em_itnb_control} was \code{TRUE}, otherwise it is \code{NA}.}
#' }
#'
#' @name itnb-object
#'
#' @rdname itnb-object
NULL
