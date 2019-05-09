library(Rssa)
alldata <- read.table("prod.cts", header = TRUE)
series <- alldata$Series

series <- series[!is.na(series)]

# раскомментировать для построения картинки с BIC
# bic_data <- make_bic_data(series)
# plot_bic_data(bic_data)
# plot_bic_data(bic_data[bic_data$p > 0 & bic_data$r > 8, ])
# stop()

#best model chosen
r = 14
# r = 12
p = 3
answer <- arbitrary_noise_optimize(series, r = r, p = p, debug = TRUE)
plot(as.numeric(series), type = "l")
lines(answer$signal, col = "red")
plot(answer$noise, main = "noise", type = "l")
