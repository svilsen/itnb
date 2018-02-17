if (FALSE) {
    n = 1000
    k = 94
    t = 93
    mu = 100
    theta = 10
    p = 0.1

    x <- rkitnb(n, k, t, mu, theta, p)

    hist(x, breaks = "FD")
    plot(dkitnb(x, k, t, mu, theta, p)~x)

    ko <- kitnb_optimisation(x, k, t, control = kitnb_optimisation.control(burnInIterations = 5))
}
