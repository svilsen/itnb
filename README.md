# itnb
The `itnb`-package implements mean and overdispersion parameterised $i$-inflated and $t$-truncated negative binomial distribution (itnb-distribution). The package also implements an expectation-maximisation (EM) algorithm for estimating the mean, overdispersion, and inflation parameters for counts generated from an itnb-distribution. Furthermore, both a non-parametric and a parametric bootstrap is implemented to construct confidence envelopes for the estimated parameters.

## Installation

The `itnb`-package depends on `R` (>= 4.1). As the package is not available on CRAN, devtools is needed to install the package from github. 

From R, run the following commands:  

```r
install.packages("devtools")
devtools::install_github("svilsen/itnb")
```

## Usage

```r
library("itnb")

## Setting parameters
n <- 2000
i <- 94
t <- 93
mu <- 100
theta <- 10
p <- 0.2

## Generating random variates
x <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)

## Plotting sample
mids <- hist(x, breaks = seq(min(x), max(x)), plot = FALSE)$mids
hist(x, breaks = seq(min(x), max(x)), probability = TRUE); points(mids, ditnb(mids, mu, theta, p, i, t), type = "l", lwd = 2, col = "dodgerblue2")

## Estimating parameters
m <- em_itnb(
    x = x, 
    i = i, 
    t = t, 
    control = em_itnb_control(
        trace = 0, 
        iteration_min = 1000,
        save_trace = TRUE
    )
)

## Plotting trace of EM-algorithm 
plot(m, log = "x")

## Percentage errors
100 * (m$mu - mu) / mu
100 * (m$theta - theta) / theta
100 * (m$p - p) / p

## Sampling CI's
ci_p <- simulate_ci(
    m, 
    level = 0.95, 
    trace = TRUE, 
    nr_simulations = 250, 
    parametric = TRUE,
    plot = TRUE
)

ci_np <- simulate_ci(
    m, 
    level = 0.95, 
    trace = TRUE, 
    nr_simulations = 250, 
    parametric = FALSE,
    plot = TRUE
)
```

## License

This project is licensed under the MIT License.
