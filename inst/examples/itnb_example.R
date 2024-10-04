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
y <- itnb::ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
data <- data.frame(y = y, x = x)

#
itnb(
    y ~ x,
    data = data,
    i = i,
    t = t,
    link = "log"
)
