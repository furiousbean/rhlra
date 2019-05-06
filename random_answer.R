gap <- 20
series <- c(1, numeric(gap), 1)
answer <- white_noise_optimize(series, r = 1)
plot(answer$signal)

seedf <- function() set.seed(300)
answer <- white_noise_optimize(series, r = 1, set_seed = seedf)
plot(answer$signal)
