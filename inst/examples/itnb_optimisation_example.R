library("ktnb")

n = 2000
k = 94
t = 93
mu = 100
theta = 10
p = 0.2

x <- rktnb(n, k, t, mu, theta, p)
ktnb_object <- ktnb_optimisation(x, k, t)
