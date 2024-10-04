#
n <- 2000
i <- 2
t <- 1

#
x <- sort(runif(n, 0, 4))

beta <- c(3, 2)
mu <- exp(cbind(1, x) %*% beta)

theta <- 10
p <- 0.2

#
y <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
data <- data.frame(y = y, x = x)

#
m <- itnb(
    y ~ x,
    data = data,
    i = i,
    t = t,
    link = "log"
)

#
simulate_ci(
    object = m,
    level = 0.95,
    nr_simulations = 200,
    parametric = FALSE
)
