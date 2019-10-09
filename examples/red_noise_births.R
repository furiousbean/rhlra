library(rhlra)

data <- read.table("births3.ctu", sep = ",")
series <- as.numeric(data$V1)

# раскомментировать для построения картинки с BIC
# bic_data <- hlra_tune(series, r_range = 1:16, p_range = 0:2,
#                           alpha = 0.3)
# plot(bic_data)
# plot(bic_data[bic_data$p > 0 & bic_data$r > 10, ])
# stop()

#best model chosen
r <- 11
p <- 1

answer <- hlra_ar(series, r = r, p = p, alpha = 0.1)
# answer <- hlra(series, r = r, alpha = 0.1, debug = TRUE)
plot(as.numeric(series), type = "l")
lines(answer$signal, col = "red")
plot(answer$noise, main = "noise", type = "l")
