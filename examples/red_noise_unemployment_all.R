library(Rssa)
library(rhlra)
library(snow)
data("USUnemployment")
series <- list(M_TEEN = USUnemployment[, 1], F_TEEN = USUnemployment[, 3])

#best model chosen
r <- 16
p <- 2
alpha = .8

#alpha = 0.3 for r = 11
#alpha = 0.5 for r = 12-13
#alpha = 0.8 for r = 14

# cl <- makeCluster(getOption("cl.cores", 8))
#
# bic_data <- hlra_tune(series, r_range = 1:16, p_range = 0:3,
#                       alpha = alpha, cluster = cl,
#                       initial_ar_coefs = list(c(.9), c(.9)))
#
# stopCluster(cl)
#
# plot(bic_data)
# plot(bic_data[bic_data$p > 0, ])

answer <- hlra_ar(series, r = r, p = p, alpha = alpha,
                  initial_ar_coefs = list(c(.9), c(.9)), debug = T)

# answer <- hlra_ar(series, r = r, p = p, alpha = .3, debug = TRUE,
#                                    initial_coefs = c(.9))

plot(series$M_TEEN)
lines(answer$signal$M_TEEN, col = "red")
plot(answer$noise$M_TEEN, main = "noise", type = "l")

plot(series$F_TEEN)
lines(answer$signal$F_TEEN, col = "red")
plot(answer$noise$F_TEEN, main = "noise", type = "l")
