convolve_with_1 <- function(x, L) {
    K <- length(x)
    N <- K + L - 1
    answer <- numeric(N)
    answer[1:K] <- x
    answer[(N-K+2):(N)] <- answer[(N-K+2):(N)] - x[-K]
    cumsum(answer)
}

summate_with_1 <- function(x, L) {
    x <- c(0, x)
    Np <- length(x)
    K <- Np - L
    x <- cumsum(x)
    x[(Np - K + 1):Np] - x[1:K]
}

kloptimw <- function(N, L, alpha, sigmasq = rep(1, N)) {
    K <- N - L + 1

    counter_fun <- 0
    counter_grad <- 0

    sigmasq[sigmasq == 0] <- max(sigmasq) * 1e-3

    optfun <- function(x, L, sigmas, alpha) {
        wgts <- convolve_with_1(x, L) * sigmasq

        C <- N/sum(wgts)
        wgts <- wgts * C

        counter_fun <<- counter_fun + 1

        sum(wgts - log(wgts))
    }

    gradfun <- function(x, L, sigmas, alpha) {
        wgts <- convolve_with_1(x, L)

        C <- N/sum(wgts * sigmasq)

        counter_grad <<- counter_fun + 1

        dC <- - N/(sum(wgts * sigmasq)^2)

        summate_with_1(dC * wgts * sigmasq + C * sigmasq - dC/C - 1/wgts, L)
    }

    soltry <- optim(rep((alpha + 1)/2, K), optfun, gradfun, method = "L-BFGS-B",
                    lower = rep(alpha, K), upper = rep(1, K),
                    L = L, sigmas = sigmas, alpha = alpha)

    soltry$par
}
