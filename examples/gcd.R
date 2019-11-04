x <- hlra_sylvester(list(c(1, -3, 3, -1), c(1, -2, 1)), 1, debug = F)
print(x$gcd)

x <- hlra_sylvester(list(c(1.1, -3, 3.5, -1.1), c(1, -2, 1)), 1, poly_weights = c(1, 10000))
print(x$gcd)

x <- hlra_sylvester(list(c(1, -3, 3, -1), c(1, -2, 1)), 2)
print(x$gcd)

library(fftw)

ssa_convolve <- function(u, v) {
    p <- planFFT(length(u) + length(v) - 1)
    l_fft <- FFT(c(u, numeric(length(v) - 1)), plan = p)
    r_fft <- FFT(c(v, numeric(length(u) - 1)), plan = p)
    Re(IFFT(l_fft * r_fft, plan = p))
}

set.seed(15)

x <- rnorm(10)
y <- rnorm(20)
z <- rnorm(30)

x <- x/sqrt(sum(x^2))

sd <- 0.01

first_poly <- ssa_convolve(x, y) + rnorm(29) * sd
second_poly <- ssa_convolve(x, z) + rnorm(39) * sd

r <- length(x) - 1

r <- 9
result <- hlra_sylvester(list(first_poly, second_poly), r, debug = T)

guess_order <- function(polynoms, r_range = 1:15) {
    dists <- sapply(r_range, function(r) {
        print(r)
        result <- hlra_sylvester(polynoms, r)
        sqrt(sum(sapply(result$residual, function(x) sum(x^2))))
    })

    data.frame(r = r_range, dists = dists)
}

dists <- guess_order(list(first_poly, second_poly))
plot(dists$r, dists$dists, log = "y")

print(x)
print(result$gcd)
