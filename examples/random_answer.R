library(rhlra)

gap <- 18
series <- c(1, numeric(gap), 1)
answer <- hlra(series, r = 1, debug = TRUE)
plot(answer$signal)

seedf <- function() set.seed(300)
answer <- hlra(series, r = 1, set_seed = seedf, debug = T)
plot(answer$signal)
