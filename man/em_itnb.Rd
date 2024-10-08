% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em-itnb.R
\name{em_itnb}
\alias{em_itnb}
\alias{em_itnb.numeric}
\title{itnb parameter optimisation}
\usage{
em_itnb(x, i, t, control = list())

\method{em_itnb}{numeric}(x, i, t, control = list())
}
\arguments{
\item{x}{Vector: Observed counts.}

\item{i}{Numeric: The inflation point.}

\item{t}{Numeric: The truncation point.}

\item{control}{List: A control object, see \link{em_itnb_control} for details.}
}
\value{
An object of class \link{itnb-object}.
}
\description{
Parameter optimisation for data generated by a r.v. following an itnb distribution using the ENM algorithm.
}
\examples{
n = 2000
i = 94
t = 93
mu = 100
theta = 10
p = 0.2

x <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
em_itnb(x = x, i = i, t = t, control = em_itnb_control(trace = 100))
}
