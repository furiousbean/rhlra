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

boxoptimw <- function(N, L, alpha, sigmas = rep(1, N)) {
    K <- N - L + 1
    
    counter_fun <- 0
    counter_grad <- 0
    
    
    optfun <- function(x, L, sigmas, alpha) {
        wgts <- convolve_with_1(x, L)
        ws <- wgts * sigmas
        
        gamma <- sum(ws) / sum(ws * ws)
        ws <- ws * gamma
        
        counter_fun <<- counter_fun + 1
        
        sum(ws * ws - 2 * ws)
    }
    
    gradfun <- function(x, L, sigmas, alpha) {
        wgts <- convolve_with_1(x, L)
        ws <- wgts * sigmas
        
        counter_grad <<- counter_fun + 1
        
        gamma_nom <- sum(ws) 
        gamma_denom <- sum(ws * ws)
        gamma <- gamma_nom / gamma_denom
        
        summate_s <- summate_with_1(sigmas, L)
        summate_ws <- summate_with_1(ws * sigmas, L)
        
        gamma_der <- (summate_s / (gamma_denom * gamma_denom) -
                          2 * gamma_nom/(gamma_denom * gamma_denom) * 
                          summate_ws)
        
        (2 * gamma * gamma_der * gamma_denom + 
                2 * gamma * gamma * summate_ws -
                2 * gamma_der * gamma_nom - 2 * gamma * summate_s)
    }
    
    soltry <- optim(rep((alpha + 1)/2, K), optfun, gradfun, method = "L-BFGS-B",
                    lower = rep(alpha, K), upper = rep(1, K),
                    L = L, sigmas = sigmas, alpha = alpha)
    
    soltry$par
}
