# itnb
The `itnb`-package implements mean and overdispersion parameterised $i$-inflated and $t$-truncated negative binomial distribution (itnb-distribution). The package also implements an expectation-maximisation (EM) algorithm for estimating the mean, overdispersion, and inflation parameters for counts generated from an itnb-distribution. Furthermore, both a non-parametric and a parametric bootstrap is implemented to construct confidence envelopes for the estimated parameters.

## Installation

The `itnb`-package depends on `R` (>= 4.0). As the package is not available on CRAN, devtools is needed to install the package from github. 

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
hist(x, breaks = "fd")

## Estimating the parameters
m <- em_itnb(
    x = x, 
    i = i, 
    t = t, 
    control = em_itnb_control(
        trace = 0L, 
        save_trace = TRUE
    )
)

##
plot(m, log = "x")

##
simulate_confidence_envelopes(
    m, 
    level = NULL, 
    trace = TRUE, 
    number_of_simulations = 100, 
    plot_simulations = TRUE
)
```

## License

This project is licensed under the MIT License.
