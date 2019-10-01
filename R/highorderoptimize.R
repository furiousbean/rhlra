
hlra_residuals <- function(noise, ar_coefs) {
    norm_mat <- band_mat_from_diags(inv_ac_diags(length(noise), ar_coefs))
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

hlra_mgn_default_chol <- function(series) {
    make_series_01_vector <- function(s) {
        answer <- rep(1, length(s))
        answer[is.na(s)] <- 0
        list(answer)
    }

    if (!is.list(series)) {
        return(band_mat_from_diags(make_series_01_vector(series)))
    } else {
        return(sapply_ns(band_mat_from_diags, make_series_01_vector(series)))
    }
}

hlra_cadzow <- function(series, r, L = default_L(series),
                        left_diags = NULL, right_diags = NULL,
                        debug = FALSE, set_seed = NULL, additional_pars = list()) {

    obj <- list()

    classes <- character(0)

    if (!is.list(series)) {
        if (is.null(left_diags)) {
            left_diags <- list(rep(1, L))
            right_diags <- list(rep(1, length(series) - L + 1))
        }


        if (any(is.na(series))) {
            stop("Input series must not contain NA's")
        }
        classes <- c(classes, "1d")
    } else {
        if (is.null(left_diags)) {
            left_diags <- sapply(seq_along(series),
                                 function(i) list(rep(1, L)),
                                 simplify = FALSE)
            right_diags <- sapply(seq_along(series),
                                  function(i) list(rep(1, length(series[[i]]) - L + 1)),
                                  simplify = FALSE)
        }

        for (i in seq_along(series)) {
            if (any(is.na(series[[i]]))) {
                stop("Input series must not contain NA's")
            }
        }
        classes <- c(classes, "1dm")
    }

    class(obj) <- classes

    signal_obj <- cadzow(obj, series, r,
                         left_diags, right_diags,
                         debug, set_seed = set_seed, additional_pars = additional_pars)

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}

hlra_mgn <- function(series, initial_glrr, weights = NULL,
                     weights_chol = hlra_mgn_default_chol(series),
                     debug = FALSE, compensated = TRUE, additional_pars = list()) {

    obj <- list()

    classes <- character(0)

    if (!is.list(series)) {
        classes <- c(classes, "1d")
    } else {
        classes <- c(classes, "1dm")

    }

    if (compensated) {
        classes <- c(classes, "compensated")
    }

    class(obj) <- classes

    r <- length(initial_glrr) - 1

    mgn_call_list = list(this = obj,
                         series = series,
                         signal = NULL,
                         r = r,
                         weights = weights,
                         weights_chol = weights_chol,
                         glrr_initial = initial_glrr,
                         debug = debug)

    signal_obj <- do.call(mgn, expand_pars_list(mgn_call_list, mgn_add_pars_names, additional_pars))

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}


hlra <- function(series, r, L = default_L(series), ar_coefs = NULL,
                 alpha = 0.1, debug = FALSE,
                 envelope = unit_envelope(series),
                 compensated = TRUE, set_seed = NULL, additional_pars = list()) {

    obj <- list()

    classes <- character(0)

    if (is.null(ar_coefs)) {
        if (!is.list(series)) {
            ar_coefs <- numeric(0)
        } else {
            ar_coefs <- sapply(seq_along(series), function(i) numeric(0), simplify = FALSE)
        }
    }

    if (compensated) {
        classes <- c(classes, "compensated")
    }

    # rude phi estimate
    signal_obj <- NULL
    series_for_cadzow <- series

    if (!is.list(series)) {
        right_diag <- boxoptimw(length(series), L, alpha, envelope^2)

        if (any(is.na(series))) {
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed,
                                           additional_pars = additional_pars)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) boxoptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
                             simplify = FALSE)

        for (i in seq_along(series)) {
            if (any(is.na(series[[i]]))) {
                series_for_cadzow[[i]] <- fill_gaps(series[[i]][effective_mask(series[[i]])], r, debug, set_seed = set_seed,
                                                    additional_pars = additional_pars)
            }
        }
        classes <- c(classes, "1dm")
    }

    class(obj) <- classes

    signal_obj <- cadzow_with_mgn(obj, series, L, r, ar_coefs,
                                 right_diag, debug = debug, envelope = envelope,
                                 series_for_cadzow = series_for_cadzow,
                                 set_seed = set_seed, additional_pars = additional_pars)

    if (!is.list(series)) {
        N <- length(as.numeric(signal_obj$signal))
    } else {
        N <- sum(sapply(signal_obj$signal, function(series) length(as.numeric(series))))
    }

    signal_obj$bic <- -log(N) * signal_obj$df  + 2 * signal_obj$loglikelihood
    signal_obj$ar_coefs <- ar_coefs

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}

hlra_ar <- function(series, r, p = 1, L = default_L(series),
                    alpha = 0.1, k = p * 4, ar_coefs_eps = 1e-7,
                    initial_ar_coefs = NULL, debug = FALSE,
                    envelope = unit_envelope(series),
                    compensated = TRUE, set_seed = NULL, additional_pars = list()) {

    if (p == 0) {
        return(hlra(series, r, L, alpha, debug = debug,
                    envelope = envelope, compensated = compensated,
                    set_seed = set_seed, additional_pars = additional_pars))
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
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed,
                                           additional_pars = additional_pars)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) boxoptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
                             simplify = FALSE)
        Ks <- sapply(series, function(series) length(as.numeric(series))) - L + 1
        for (i in seq_along(series)) {
            if (any(is.na(series[[i]]))) {
                series_for_cadzow[[i]] <- fill_gaps(series[[i]][effective_mask(series[[i]])], r, debug, set_seed = set_seed,
                                                    additional_pars = additional_pars)
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


    if (is.null(initial_ar_coefs)) {
        signal_obj <- hlra_ssa(obj, series_for_cadzow, L, r, debug, set_seed, additional_pars)

        noise <- signal_obj$noise
        ar_coefs <- estimate_coefs(noise)
    } else {
        ar_coefs <- initial_ar_coefs
    }

    if (!is.list(series)) {
        ar_coefs_all <- matrix(ar_coefs, p, 1)
    } else {
        ar_coefs_all <- matrix(glue_series_lists(ar_coefs), p * length(series), 1)
    }

    for (i in 1:k) {
        if (debug) cat(sprintf("%d EM iteration\n", i))

        signal_obj <- cadzow_with_mgn(obj, series, L, r, ar_coefs, right_diag, debug = debug,
                                     envelope = envelope, series_for_cadzow = series_for_cadzow,
                                     set_seed = set_seed, additional_pars = additional_pars)

        noise <- signal_obj$noise
        ar_coefs_old <- ar_coefs
        ar_coefs <- estimate_coefs(noise)

        if (!is.list(series)) {
            ar_coefs_all <- cbind(ar_coefs_all, ar_coefs)
            if (sum(abs(ar_coefs_old - ar_coefs)) < ar_coefs_eps) break
        } else {
            ar_coefs_all <- cbind(ar_coefs_all, glue_series_lists(ar_coefs))
            if (sum(abs(glue_series_lists(ar_coefs_old) - glue_series_lists(ar_coefs))) < ar_coefs_eps) break
        }
    }

    if (debug) print(ar_coefs_all)

    if (!is.list(series)) {
        N <- length(as.numeric(signal_obj$signal))
    } else {
        N <- sum(sapply(signal_obj$signal, function(series) length(as.numeric(series))))
    }

    signal_obj$bic <- -log(N) * signal_obj$df  + 2 * signal_obj$loglikelihood
    signal_obj$ar_coefs <- ar_coefs
    signal_obj$em_it <- k

    classes <- c(classes, "hlra")

    class(signal_obj) <- classes

    return(signal_obj)
}

tune_hlra <- function(series, r_range = 1:15, p_range = 0:3, 
                      L = default_L(series), alpha = 0.1,
                      envelope = unit_envelope(series), set_seed = NULL,
                      cluster = NULL, initial_ar_coefs = NULL, additional_pars = list()) {

    obtain_bic_df_nonv <- function(v) {
        r <- v[1]
        p <- v[2]
        cat(sprintf("Try model: r = %d, p = %d\n", r, p))
        opt_obj <- hlra_ar(series, r, p, L, alpha,
                           envelope = envelope, set_seed = set_seed,
                           initial_ar_coefs = initial_ar_coefs, additional_pars = additional_pars)
        bic <- opt_obj$bic
        df <- opt_obj$df
        loglikelihood <- opt_obj$loglikelihood
        c(bic, df, loglikelihood, r, p)
    }

    p_all <- rep(p_range, length(r_range))
    r_all <- rep(r_range, each = length(p_range))

    input_mat <- cbind(r_all, p_all)
    input_mat <- input_mat[sample(length(r_range) * length(p_range)), ]

    data <- NULL

    if (!is.null(cluster)) {
        clusterCall(cluster, function() {
            library(svd)
            library(Matrix)
            library(rhlra)
        })

        clusterExport(cluster, c("series", "L", "alpha", "envelope", "set_seed",
            "initial_ar_coefs", "additional_pars"), envir = environment())

        data <- parApply(cluster, input_mat, 1, obtain_bic_df_nonv)
    } else {
        data <- apply(input_mat, 1, obtain_bic_df_nonv)
    }

    data <- as.data.frame(t(data))
    names(data) <- c("bic", "df", "loglikelihood", "r", "p")

    neworder <- order(data$df)

    data <- data[neworder,]

    class(data) <- c("hlra_tune", class(data))
    data
}

plot.hlra_tune <- function(data) {
    p_range <- unique(data$p)
    plot(data$df, data$bic, type = "l",
         xlab = "Degrees of freedom", main = "BIC comparison", ylab = "BIC")
    text(data$df, data$bic, label = as.character(data$r), col = data$p + 1)
    legend("topright", legend = paste("p =", as.character(p_range)),
           col = p_range + 1, pch = 2)
}
