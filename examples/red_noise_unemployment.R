library(Rssa)
library(rhlra)
library(normwhn.test)
data("USUnemployment")
series <- USUnemployment[, 2]
series <- series[!is.na(series)]

seedf <- function() set.seed(15)

# раскомментировать для построения картинки с BIC

seedf()

alpha <- .8

envelope <- reconstruct(ssa(series, L = 80), groups = list(1:1))[[1]]

bic_data <- make_bic_data(series, r_range = 1:16, p_range = 0:3,
                          alpha = alpha, cores = 8,
                          initial_coefs = c(.9), set_seed = seedf)
plot_bic_data(bic_data)
plot_bic_data(bic_data[bic_data$p > 0, ])
stop()

#best model chosen
r <- 9
p <- 3

#alpha = 0.3 for r = 11
#alpha = 0.5 for r = 12-13
#alpha = 0.8 for r = 14


answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = alpha,
                                   initial_coefs = c(.9), set_seed = seedf, debug = T)

# answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = .3, debug = TRUE,
#                                    initial_coefs = c(.9))

plot(as.numeric(series), type = "l")
lines(answer$signal, col = "red")
plot(answer$noise/envelope, main = "noise", type = "l")
