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
i <- 2
t <- 1

x <- sort(runif(n, 0, 4))
beta <- c(1, 2)
mu <- cbind(1, x) %*% beta
theta <- 10
p <- 0.2

## Generating random variates
y <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
data <- data.frame(y = y, x = x)

## Estimating parameters
m <- itnb(
    y ~ x, 
    data = data,
    i = i,
    t = t,
    control = itnb_control(
        trace = 0
    )
)

## Plotting trace of EM-algorithm 
plot(m, log = "x")

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
