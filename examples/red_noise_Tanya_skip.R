library(rhlra)

alldata <- read.table("prod.cts", header = TRUE)
series <- alldata$Series
true_series <- series
# series <- c(series, , rep(NA, 50))
series[10:49] <- NA

# раскомментировать для построения картинки с BIC
# bic_data <- hlra_tune(series)
# plot(bic_data)
# plot(bic_data[bic_data$p > 0 & bic_data$r > 8, ])
# stop()

#best model chosen
r = 14
# r = 12
p = 3
# answer <- hlra(series, r = r,
                                 # coefs = c(0.60983015, 0.08978667, 0.17555034), debug = TRUE)

answer <- hlra_ar(series, r = r, p = 3, debug = TRUE)


matplot(1:length(answer$signal),
        cbind(answer$signal,
        as.numeric(true_series),
        as.numeric(series)),
        type = "l")

