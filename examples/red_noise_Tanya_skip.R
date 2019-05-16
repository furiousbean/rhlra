library(rhlra)

alldata <- read.table("prod.cts", header = TRUE)
series <- alldata$Series
# series <- c(series, , rep(NA, 50))
series[10:49] <- NA

# раскомментировать для построения картинки с BIC
# bic_data <- make_bic_data(series)
# plot_bic_data(bic_data)
# plot_bic_data(bic_data[bic_data$p > 0 & bic_data$r > 8, ])
# stop()

#best model chosen
r = 14
# r = 12
p = 3
# answer <- certain_noise_optimize(series, r = r,
                                 # coefs = c(0.60983015, 0.08978667, 0.17555034), debug = TRUE)

answer <- arbitrary_noise_optimize(series, r = r, p = 3)


matplot(1:length(answer$signal),
        cbind(answer$signal, as.numeric(series)),
        type = "l")

