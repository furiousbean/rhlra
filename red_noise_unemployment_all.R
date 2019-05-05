library(Rssa)
library(rhlra)
library(normwhn.test)
data("USUnemployment")
series <- list(M_TEEN = USUnemployment[, 1], F_TEEN = USUnemployment[, 3])

#best model chosen
r <- 14
p <- 3

#alpha = 0.3 for r = 11
#alpha = 0.5 for r = 12-13
#alpha = 0.8 for r = 14


answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = .3,
                                   initial_coefs = list(c(.9), c(.9)),
                                   debug = TRUE)

# answer <- arbitrary_noise_optimize(series, r = r, p = p, alpha = .3, debug = TRUE,
#                                    initial_coefs = c(.9))

plot(as.numeric(series[[1]]), type = "l")
lines(answer$signal[[1]], col = "red")
plot(answer$noise[[1]], main = "noise", type = "l")

plot(as.numeric(series[[2]]), type = "l")
lines(answer$signal[[2]], col = "red")
plot(answer$noise[[2]], main = "noise", type = "l")
