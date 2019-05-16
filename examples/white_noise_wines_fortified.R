library(rhlra)
library(Rssa)

data(AustralianWine)
series <- AustralianWine[, 3]
series <- series[!is.na(series)]

# раскомментировать для построения картинки с BIC

bic_data <- make_bic_data(series)
plot_bic_data(bic_data)
plot_bic_data(bic_data[bic_data$r > 10, ])

# stop()

#best model chosen
r <- 15
answer <- white_noise_optimize(series, r = r)
plot(as.numeric(series), type = "l")
lines(answer$signal, col = "red")
plot(answer$noise, main = "noise", type = "l")
