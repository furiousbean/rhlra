
residuals <- function(noise, coefs) {
    norm_mat <- band_mat_from_diags(inv_ac_diags(length(noise), coefs))
    left_chol <- Cholesky(norm_mat, perm = FALSE, LDL = FALSE)
    as.numeric(as(left_chol, "Matrix") %*% noise)
}

default_L <- function(series) {
    if (!is.list(series)) {
        return(effective_length(series) %/% 2)
    } else {
        return(min(
            round(sum(sapply(series, effective_length)) / (length(series) + 1)),
        min(sapply(series, effective_length))))
    }
}

certain_noise_optimize <- function(series, L = default_L(series), r, coefs,
                                 alpha = 0.1, randomsearch = FALSE, debug = FALSE,
                                 envelope = unit_envelope(series), compensated = TRUE, set_seed = NULL, ...) {

    obj <- list()

    classes <- character(0)

    if (compensated) {
        classes <- c(classes, "compensated")
    }

    # rude phi estimate
    signal_obj <- NULL
    series_for_cadzow <- series

    if (!is.list(series)) {
        right_diag <- boxoptimw(length(series), L, alpha, envelope^2)

        if (any(is.na(series))) {
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) boxoptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
                             simplify = FALSE)

        for (i in seq_along(series)) {
            if (any(is.na(series[[i]]))) {
                series_for_cadzow[[i]] <- fill_gaps(series[[i]][effective_mask(series[[i]])], r, debug, set_seed = set_seed)
            }
        }
        classes <- c(classes, "1dm")
    }

    class(obj) <- classes

    signal_obj <- cadzow_with_mgn(obj, series, L, r, coefs,
                                 right_diag, debug = debug, envelope = envelope,
                                 series_for_cadzow = series_for_cadzow,
                                 set_seed = set_seed, ...)

    if (!is.list(series)) {
        N <- length(as.numeric(signal_obj$signal))
    } else {
        N <- sum(sapply(signal_obj$signal, function(series) length(as.numeric(series))))
    }

    signal_obj$bic <- -log(N) * signal_obj$df  + 2 * signal_obj$loglikelihood
    signal_obj$coefs <- coefs

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}

white_noise_optimize <- function(series, L = default_L(series), r,
                                 alpha = 0.1, randomsearch = FALSE, debug = FALSE,
                                 envelope = unit_envelope(series),
                                 compensated = TRUE, set_seed = NULL, ...) {
    if (!is.list(series)) {
        coefs <- numeric(0)
    } else {
        coefs <- sapply(seq_along(series), function(i) numeric(0), simplify = FALSE)
    }
    certain_noise_optimize(series, L, r, coefs, alpha, randomsearch, debug,
                                  envelope = envelope,
                                  set_seed = set_seed, ...)
}

arbitrary_noise_optimize <- function(series, L = default_L(series), r, p = 1,
                               alpha = 0.1, k = p * 4, coef_eps = 1e-7,
                               initial_coefs = NULL,
                               randomsearch = FALSE, debug = FALSE,
                               envelope = unit_envelope(series), 
                               compensated = TRUE, set_seed = NULL, ...) {

    if (p == 0) {
        return(white_noise_optimize(series, L, r, alpha, debug = debug,
                                    envelope = envelope, compensated = compensated,
                                    set_seed = set_seed, ...))
    }

    obj <- list()

    classes <- character(0)

    if (compensated) {
        classes <- c(classes, "compensated")
    }

    # rude phi estimate
    signal_obj <- NULL
    series_for_cadzow <- series

    if (!is.list(series)) {
        right_diag <- boxoptimw(length(series), L, alpha, envelope^2)
        K <- length(series) - L + 1
        if (any(is.na(series))) {
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) boxoptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
                             simplify = FALSE)
        Ks <- sapply(series, function(series) length(as.numeric(series))) - L + 1
        for (i in seq_along(series)) {
            if (any(is.na(series[[i]]))) {
                series_for_cadzow[[i]] <- fill_gaps(series[[i]][effective_mask(series[[i]])], r, debug, set_seed = set_seed)
            }
        }
        classes <- c(classes, "1dm")
    }

    class(obj) <- classes   

    if (!is.list(series)) {
        estimate_coefs <- function(x) {
            answer <- numeric(p)
            tryCatch({  
                aobj <- ar(x/envelope, order.max = p,
                           na.action = na.pass)
                if (length(aobj$ar) > 0) answer[1:length(aobj$ar)] <- aobj$ar
                answer
                }, error = function(x) {print(x)})

            answer
        }
    } else {
        estimate_coefs <- function(x) {
            answer <- sapply(seq_along(series), function(i) numeric(p))
            tryCatch({
                answer <- sapply(seq_along(series), function(i)
                    {
                        aobj <- ar(x[[i]]/envelope[[i]], order.max = p,
                                   na.action = na.pass)
                        ans <- numeric(p)
                        if (length(aobj$ar) > 0) ans[1:length(aobj$ar)] <- aobj$ar
                        ans
                    }, simplify = FALSE)
                answer
            }, error = function(x) {print(x)}
            )
            answer
        }
    }


    if (is.null(initial_coefs)) {
        if (!is.list(series)) {
            pseudord <- rep(1, K)
            whitecoefs <- numeric(0)
        } else {
            pseudord <- sapply(Ks, function(i) rep(1, i), simplify = FALSE)
            whitecoefs <- sapply(seq_along(series), function(i) numeric(0), simplify = FALSE)
        }

        signal_obj <- cadzow_with_mgn(obj, series, L, r, coefs = whitecoefs,
                                     right_diag = pseudord, epsilon = 1, use_mgn = FALSE, debug = debug,
                                     envelope = envelope, series_for_cadzow = series_for_cadzow,
                                     set_seed = set_seed, ...)

        noise <- signal_obj$noise
        coefs <- estimate_coefs(noise)
    } else {
        coefs <- initial_coefs
    }

    if (!is.list(series)) {
        coefs_all <- matrix(coefs, p, 1)
    } else {
        coefs_all <- matrix(glue_series_lists(coefs), p * length(series), 1)
    }

    for (i in 1:k) {
        if (debug) cat(sprintf("%d EM iteration\n", i))

        signal_obj <- cadzow_with_mgn(obj, series, L, r, coefs, right_diag, debug = debug,
                                     envelope = envelope, series_for_cadzow = series_for_cadzow,
                                     set_seed = set_seed, ...)

        noise <- signal_obj$noise
        coefs_old <- coefs
        coefs <- estimate_coefs(noise)

        if (!is.list(series)) {
            coefs_all <- cbind(coefs_all, coefs)
            if (sum(abs(coefs_old - coefs)) < coef_eps) break
        } else {
            coefs_all <- cbind(coefs_all, glue_series_lists(coefs))
            if (sum(abs(glue_series_lists(coefs_old) - glue_series_lists(coefs))) < coef_eps) break
        }
    }

    if (debug) print(coefs_all)

    if (!is.list(series)) {
        N <- length(as.numeric(signal_obj$signal))
    } else {
        N <- sum(sapply(signal_obj$signal, function(series) length(as.numeric(series))))
    }

    signal_obj$bic <- -log(N) * signal_obj$df  + 2 * signal_obj$loglikelihood
    signal_obj$coefs <- coefs
    signal_obj$em_it <- k

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}

make_bic_data <- function(series, L = default_L(series), alpha = 0.1,
                             r_range = 1:15, p_range = 0:3,
                          envelope = unit_envelope(series), set_seed = NULL, 
                          cores = 4, initial_coefs = NULL) {

    obtain_bic_df_nonv <- function(v) {
        r <- v[1]
        p <- v[2]
        cat(sprintf("Try model: r = %d, p = %d\n", r, p))
        opt_obj <- arbitrary_noise_optimize(series, L, r, p, alpha,
                                            envelope = envelope, set_seed = set_seed,
                                            initial_coefs = initial_coefs)
        bic <- opt_obj$bic
        df <- opt_obj$df
        loglikelihood <- opt_obj$loglikelihood
        c(bic, df, loglikelihood, r, p)
    }

    # obtain_bic_df <- Vectorize(obtain_bic_df_nonv)
    p_all <- rep(p_range, length(r_range))
    r_all <- rep(r_range, each = length(p_range))

    cl <- makeCluster(getOption("cl.cores", cores))
    clusterCall(cl, function() library(svd))
    clusterCall(cl, function() library(fftw))
    clusterCall(cl, function() library(Rssa))
    clusterCall(cl, function() library(Matrix))
    clusterCall(cl, function() library(quadprog))
    clusterExport(cl, ls(envir = globalenv()))

    input_mat <- cbind(r_all, p_all)
    input_mat <- input_mat[sample(length(r_range) * length(p_range)), ]

    # data <- apply(input_mat, 1, obtain_bic_df_nonv)
    data <- parApply(cl, input_mat, 1, obtain_bic_df_nonv)
    
    stopCluster(cl)

    data <- as.data.frame(t(data))
    names(data) <- c("bic", "df", "loglikelihood", "r", "p")

    neworder <- order(data$df)

    data <- data[neworder,]

    data
}

plot_bic_data <- function(data) {
    p_range <- unique(data$p)
    plot(data$df, data$bic, type = "l",
         xlab = "Degrees of freedom", main = "BIC comparison", ylab = "BIC")
    text(data$df, data$bic, label = as.character(data$r), col = data$p + 1)
    legend("topright", legend = paste("p =", as.character(p_range)),
           col = p_range + 1, pch = 2)
}
