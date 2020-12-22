# itnb
The `itnb`-package ... 

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

##
n <- 2000
i <- 94
t <- 93
mu <- 100
theta <- 10
p <- 0.2

##
x <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)

##
itnb_object <- itnb_optimisation(x = x, i = i, t = t, control = itnb_optimisation_control(trace = FALSE, save_trace = TRUE))

##
ggplot(itnb_object)

##
simulate_confidence_envelopes(itnb_object, level = 0.95, trace = TRUE, number_of_simulations = 100, plot_simulations = TRUE)
```

## License

This project is licensed under the MIT License.
