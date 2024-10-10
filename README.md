# itnb
The `itnb`-package implements an inflated and truncated negative binomial distribution (`itnb`), as well as an interface for estimating the parameters in an `itnb` regression model. The package utilises an expectation-maximisation (EM) algorithm for estimating the mean, overdispersion, and inflation parameters when the inflation cannot be separated from the truncated negative binomial regression model, otherwise it estimates the parameters by maximum likelihood (MLE). Furthermore, parametric and non-parametric bootstraps are implemented to construct confidence envelopes for the estimated parameters.

## Installation
The `itnb`-package depends on `R` (>= 4.1), `Rcpp`, `RcppArmadillo`, `roptim`, `stats`, and `graphics`. As the package is not available on CRAN, devtools is needed to install the package from github. 

From `R`, run the following commands:  
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

beta <- c(1, 2)
theta <- 10
p <- 0.2

## Generating covariates
x <- sort(runif(n, 0, 4))
mu <- cbind(1, x) %*% beta

## Generating response
y <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
data <- data.frame(y = y, x = x)

## Estimating parameters
m <- itnb(
    y ~ x, 
    data = data,
    i = i,
    t = t,
    link = "identity"
)

## Plotting trace of EM-algorithm 
plot(m, log = "x")

## Simulating CI's
ci_p <- confint(
    m, 
    level = 0.95, 
    nr_simulations = 25, 
    parametric = FALSE
)

ci_np <- confint(
    m, 
    level = 0.95, 
    nr_simulations = 25, 
    parametric = TRUE
)
```

## License
This project is licensed under the MIT License.
