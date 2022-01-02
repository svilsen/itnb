library("itnb")

n = 2000
i = 94
t = 93
mu = 100
theta = 10
p = 0.2

x <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
itnb_optimisation(x = x, i = i, t = t, control = itnb_optimisation_control(trace = 100))


itnb_optimisation(x = x, i = i, t = t, control = itnb_optimisation_control(trace = 100, tolerance = 1e-20))
