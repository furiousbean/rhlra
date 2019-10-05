x <- hlra_sylvester(list(c(1, -3, 3, -1), c(1, -2, 1)), c(.1, .5), debug = T)
print(x$gcd)

x <- hlra_sylvester(list(c(1.1, -3, 3.5, -1.1), c(1, -2, 1)), c(.1, .5), poly_weights = c(1, 10000))
print(x$gcd)

x <- hlra_sylvester(list(c(1, -3, 3, -1), c(1, -2, 1)), rnorm(2))
print(x$gcd)
