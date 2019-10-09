library(Rssa)
library(rhlra)
library(snow)
data("USUnemployment")
series <- USUnemployment[, 2]
series <- series[!is.na(series)]

seedf <- function() set.seed(15)

# раскомментировать для построения картинки с BIC

seedf()

alpha <- .8

cl <- makeCluster(getOption("cl.cores", 8))

bic_data <- hlra_tune(series, r_range = 1:16, p_range = 0:3,
                          alpha = alpha, cluster = cl,
                          initial_ar_coefs = c(.9), set_seed = seedf)

stopCluster(cl)

plot(bic_data)
plot(bic_data[bic_data$p > 0, ])
stop()

#best model chosen
r <- 9
p <- 3

#alpha = 0.3 for r = 11
#alpha = 0.5 for r = 12-13
#alpha = 0.8 for r = 14


answer <- hlra_ar(series, r = r, p = p, alpha = alpha,
                  initial_ar_coefs = c(.9), set_seed = seedf, debug = T)

plot(as.numeric(series), type = "l")
lines(answer$signal, col = "red")
plot(answer$noise, main = "noise", type = "l")
