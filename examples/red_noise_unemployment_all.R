library(Rssa)
library(rhlra)
library(normwhn.test)
data("USUnemployment")
series <- list(M_TEEN = USUnemployment[, 1], F_TEEN = USUnemployment[, 3])

#best model chosen
r <- 16
p <- 2
alpha = .8

#alpha = 0.3 for r = 11
#alpha = 0.5 for r = 12-13
#alpha = 0.8 for r = 14

bic_data <- make_bic_data(series, r_range = 1:16, p_range = 0:3,
                          alpha = alpha, cores = 8,
                          initial_coefs = list(c(.9), c(.9)))
plot_bic_data(bic_data)
plot_bic_data(bic_data[bic_data$p > 0, ])

answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = alpha,
                                   initial_coefs = list(c(.9), c(.9)))

# answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = .3, debug = TRUE,
#                                    initial_coefs = c(.9))

plot(as.numeric(series[[1]]), type = "l")
lines(answer$signal[[1]], col = "red")
plot(answer$noise[[1]], main = "noise", type = "l")

plot(as.numeric(series[[2]]), type = "l")
lines(answer$signal[[2]], col = "red")
plot(answer$noise[[2]], main = "noise", type = "l")
