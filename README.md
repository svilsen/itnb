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
n <- 150
i <- 2
t <- 1

beta <- c(0.5, 2)
alpha <- 0.5
p <- 0.1

## Generating covariates
x <- sort(runif(n, 0, 4))
mu <- exp(cbind(1, x) %*% beta)

## Generating response
y <- ritnb(n = n, mu = mu, alpha = alpha, p = p, i = i, t = t)
data <- data.frame(y = y, x = x)

## Estimating parameters
m <- itnb(
    y ~ x, 
    data = data,
    i = i,
    t = t,
    link = "log"
)

summary(m)
```

## License
This project is licensed under the MIT License.
