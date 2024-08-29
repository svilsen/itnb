n = 2000
i = 94
t = 93
mu = 100
theta = 10
p = 0.2

x <- ritnb(n = n, mu = mu, theta = theta, p = p, i = i, t = t)
m <- em_itnb(x = x, i = i, t = t, control = em_itnb_control(trace = 100))
