library(rhlra)

gap <- 18
series <- c(1, numeric(gap), 1)
answer <- hlra(series, r = 1, debug = TRUE)
plot(answer$signal)

seedf <- function() set.seed(300)
answer <- hlra(series, r = 1, set_seed = seedf, debug = T)
plot(answer$signal)

answer <- hlra_mgn(series, c(.5, .3), debug = T)
plot(answer$signal)
answer <- hlra_mgn(series, c(.3, .5), debug = T)
plot(answer$signal)



series_with_gap <- c(series, NA)
answer <- hlra_mgn(series_with_gap, c(.3, .5), debug = T)
plot(answer$signal)
answer <- hlra_mgn(series_with_gap, c(.5, .3), debug = T)
plot(answer$signal)
print(answer$signal)
