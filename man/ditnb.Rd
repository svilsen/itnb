% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itnb.R
\name{ditnb}
\alias{ditnb}
\title{ditnb}
\usage{
ditnb(x, mu, theta, p, i, t = NULL, return_log = FALSE)
}
\arguments{
\item{x}{A vector of quantiles.}

\item{mu}{The expected value.}

\item{theta}{The overdispersion.}

\item{p}{The proportion of inflation.}

\item{i}{The inflation point.}

\item{t}{The truncation point.}

\item{return_log}{TRUE/FALSE, should the log pmf be returned?}
}
\value{
A vector.
}
\description{
The density of a vector x following the i-inflated t-truncated negative binomial distribution.
}
