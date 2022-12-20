library(rhlra)
library(Rssa)

data(AustralianWine)
series <- AustralianWine[, 3]

# раскомментировать для построения картинки с BIC

# bic_data <- hlra_tune(series)
# plot(bic_data)
# plot(bic_data[bic_data$r > 10, ])

# stop()

#best model chosen
r <- 15
answer <- hlra(series, r = r)
plot(series, type = "l")
lines(answer$signal, col = "red")
plot(answer$noise, main = "noise", type = "l")
