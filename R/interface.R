#' rhlra: Hankel Low-Rank Approximation for R package.
#'
#' This package provides tools for solving Hankel Structured Low-Rank Approximation problems
#'
#' @section See Also:
#' \code{\link{hlra}}, \code{\link{hlra_tune}}, \code{\link{hlra_cadzow}}, \code{\link{hlra_mgn}}.
#'
#' @docType package
#' @name rhlra
NULL


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

hlra_mgn_default_chol <- function(series, weights = 1) {
    if (is.list(series) && length(weights) == 1) {
        weights <- rep(weights, length(series))
    }

    make_series_01_vector <- function(s, weight) {
        answer <- rep(weight[1], length(s))
        answer[is.na(s)] <- 0
        list(answer)
    }

    if (!is.list(series)) {
        return(band_mat_from_diags(make_series_01_vector(series, weights)))
    } else {
        return(sapply_ns(seq_along(series), function(i) { band_mat_from_diags(make_series_01_vector(series[[i]], weights[i])) }))
    }
}

#' Perform Cadzow iterations.
#'
#' @param series Source time series, numeric vector or list of numeric vectors (without NA's)
#' @param r Desired rank.
#' @param L Window length.
#' @param left_diags Diagonals of left SVD weights matrix.
#' @param right_diags Diagonals of right SVD weights matrix.
#' @param debug Debug mode on/off switch.
#' @param set_seed Seed function for deterministic SVD decomposition
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Object of class "hlra" containing "signal" field, which is result of optimization.
#' @examples
#' x <- hlra_cadzow(1:20, r = 1)
#' plot(x$signal, type = "b")
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

    signal_obj$signal <- substitute_new_data(series, signal_obj$signal)
    signal_obj$noise <- substitute_new_data(series, signal_obj$noise)

    return(signal_obj)
}

#' Perform Modified Gauss-Newton Algorithm.
#'
#' @param series Source time series, numeric vector or list of numeric vectors (may contain NA's)
#' @param initial_glrr Initial GLRR, which is non-null numeric vector of length more than 1.
#' @param weights Optional weights matrix of class Matrix.
#' @param weights_chol Optional Cholesky decomposition of weights matrix, matrix of class Matrix
#' @param debug Debug mode on/off switch.
#' @param compensated Use Compensated Horner Scheme in MGN algorithm?
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Object of class "hlra" containing "signal" field, which is result of optimization.
#' @examples
#' gap <- 18
#' series <- c(1, numeric(gap), 1)
#' x <- hlra_mgn(series, c(.3, .5))
#' plot(x$signal, type = "b")
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

    signal_obj$signal <- substitute_new_data(series, signal_obj$signal)
    signal_obj$noise <- substitute_new_data(series, signal_obj$noise)

    return(signal_obj)
}

#' Find GCD of two polynomials using Modified Gauss-Newton Algorithm (i.e. Sylvester LRA).
#'
#' @param polynoms Source polynomials, numeric vector or list of numeric vectors
#' @param r Expected order of GCD
#' @param initial_poly Initial approximation for GCD, which is non-null numeric vector of length more than 1.
#' @param poly_weights Optional vector of weight of each polynomial.
#' @param alpha Parameter for Cadzow weights search
#' @param debug Debug mode on/off switch.
#' @param compensated Use Compensated Horner Scheme in MGN algorithm?
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Object of class "hlra_sylvester" containing "gcd" field, which is result of optimization.
#' @examples
#' x <- hlra_sylvester(list(c(1, -3, 3, -1), c(1, -2, 1)), 1)
#' print(x$gcd)
hlra_sylvester <- function(polynoms, r, initial_poly = NULL, poly_weights = NULL,
                     alpha = 0.1, debug = FALSE, compensated = TRUE,
                     additional_pars = list()) {

    obj <- list()

    classes <- character(0)

    if (!is.null(poly_weights)) {
        weights_chol = hlra_mgn_default_chol(polynoms, sqrt(poly_weights))
    } else {
        weights_chol = hlra_mgn_default_chol(polynoms)
    }

    if (!is.list(polynoms)) {
        if (is.null(initial_poly)) {
            stop("You must provide initial approximation of GCD for this kind of problem")
        }

        if (any(is.na(polynoms))) {
            stop("Input polynom must not contain NA's")
        }

        classes <- c(classes, "1d")
    } else {
        for (i in seq_along(polynoms)) {
            if (any(is.na(polynoms[[i]]))) {
                stop("Input polynoms must not contain NA's")
            }
        }

        classes <- c(classes, "1dm")
    }

    if (compensated) {
        classes <- c(classes, "compensated")
    }

    class(obj) <- classes

    if (is.list(polynoms) && is.null(initial_poly)) {
        poly_orders <- sapply(polynoms, length) - 1
        L <- max(poly_orders) + min(poly_orders)
        zero_trails <- L - sapply(polynoms, length)

        input_for_cadzow <- sapply_ns(seq_along(polynoms), function(i)
            c(numeric(zero_trails[i]), polynoms[[i]], numeric(zero_trails[i])))
        if (is.null(poly_weights)) {
            poly_weights <- rep(1, length(polynoms))
        }

        input_for_weights <- sapply_ns(seq_along(polynoms), function(i)
            c(numeric(zero_trails[i]), rep(poly_weights[i], length(polynoms[[i]])), numeric(zero_trails[i])))

        right_diag <- sapply(seq_along(polynoms),
                             function(i) kloptimw(length(input_for_cadzow[[i]]), L, alpha, input_for_weights[[i]]),
                             simplify = FALSE)

        if (is.null(additional_pars$svd_type)) {
            additional_pars$svd_type <- "svd"
        }

        ar_coefs <- sapply(seq_along(polynoms), function(i) numeric(0), simplify = FALSE)

        initial_sylvester_approx <- cadzow_with_mgn(obj, input_for_cadzow, L, r, ar_coefs,
                                                    right_diag, debug = debug,
                                                    series_for_cadzow = input_for_cadzow,
                                                    envelope = input_for_weights,
                                                    additional_pars = additional_pars,
                                                    use_mgn = FALSE,
                                                    sylvester_nulling = zero_trails,
                                                    high_rank = TRUE)$signal

        sylvmat <- do.call(cbind, sapply_ns(initial_sylvester_approx, function(x) traj_matrix(x, L)))
        glrr_signals <- as.matrix(svd(sylvmat)$u[, (L-r+1):L], nrow = L)

        glrr_signals_as_list <- sapply_ns(1:r, function(i) as.numeric(glrr_signals[, i]))
        initial_poly <- get_glrr_from_nonlrf_series(obj, glrr_signals_as_list, r)

        # print(initial_poly)
    } else {
        r <- length(initial_poly) - 1
    }

    mgn_call_list = list(this = obj,
                         series = polynoms,
                         signal = NULL,
                         r = r,
                         weights = NULL,
                         weights_chol = weights_chol,
                         glrr_initial = initial_poly,
                         debug = debug,
                         sylvester_problem = TRUE)

    signal_obj <- do.call(mgn, expand_pars_list(mgn_call_list, mgn_add_pars_names, additional_pars))

    classes <- c(classes, "hlra_sylvester")

    class(signal_obj) <- classes

    names(signal_obj)[names(signal_obj) == "glrr"] <- "gcd"
    names(signal_obj)[names(signal_obj) == "noise"] <- "approximation"
    names(signal_obj)[names(signal_obj) == "signal"] <- "residual"

    signal_obj["noise_norm"] <- NULL
    signal_obj$gcd <- signal_obj$gcd / sqrt(sum(signal_obj$gcd ^ 2))

    signal_obj$approximation <- substitute_new_data(polynoms, signal_obj$approximation)
    signal_obj$residual <- substitute_new_data(polynoms, signal_obj$residual)

    return(signal_obj)
}

#' Perform generic HLRA problem with inverted autoregressive weights.
#'
#' @param series Source time series, numeric vector or list of numeric vectors (may contain NA's)
#' @param r Desired rank.
#' @param L Window length.
#' @param ar_coefs Autoregression coefficients.
#' @param alpha Parameter for Cadzow weights search
#' @param debug Debug mode on/off switch.
#' @param envelope Optional noise envelope.
#' @param compensated Use Compensated Horner Scheme in MGN algorithm?
#' @param set_seed Seed function for deterministic SVD decomposition
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Object of class "hlra" containing "signal" field, which is result of optimization.
#' @examples
#' gap <- 18
#' series <- c(1, numeric(gap), 1)
#' x <- hlra(series, r = 1)
#' plot(x$signal, type = "b")
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
        right_diag <- kloptimw(length(series), L, alpha, envelope^2)

        if (any(is.na(series))) {
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed,
                                           additional_pars = additional_pars)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) kloptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
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

    signal_obj$signal <- substitute_new_data(series, signal_obj$signal)
    signal_obj$noise <- substitute_new_data(series, signal_obj$noise)

    return(signal_obj)
}

#' Perform generic HLRA problem with unknown inverted autoregressive weights.
#'
#' @param series Source time series, numeric vector or list of numeric vectors (may contain NA's)
#' @param r Desired rank.
#' @param p Order of AR process.
#' @param L Window length.
#' @param alpha Parameter for Cadzow weights search
#' @param k Maximum iteration of AR parameters estimation
#' @param ar_coefs_epsilon AR coefficients threshold
#' @param initial_ar_coefs Initial autoregression coefficients.
#' @param debug Debug mode on/off switch.
#' @param envelope Optional noise envelope.
#' @param compensated Use Compensated Horner Scheme in MGN algorithm?
#' @param set_seed Seed function for deterministic SVD decomposition
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Object of class "hlra" containing "signal" field, which is result of optimization.
#' @examples
#' library(Rssa)
#' data("USUnemployment")
#' series <- log(USUnemployment[, 2])
#' x <- hlra_ar(series, r = 9, p = 3, alpha = .8, initial_ar_coefs = c(.9))
#' plot(x$signal)
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
        right_diag <- kloptimw(length(series), L, alpha, envelope^2)
        K <- length(series) - L + 1
        if (any(is.na(series))) {
            series_for_cadzow <- fill_gaps(series[effective_mask(series)], r, debug, set_seed = set_seed,
                                           additional_pars = additional_pars)
        }
        classes <- c(classes, "1d")
    } else {
        right_diag <- sapply(seq_along(series),
                             function(i) kloptimw(length(series[[i]]), L, alpha, envelope[[i]]^2),
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
        if (debug == "P" || debug) cat(sprintf("%d EM iteration\n", i))

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

    if (debug == "P" || debug) {
        cat(sprintf("AR coefficients evolution:\n"))
        print(ar_coefs_all)
    }

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

    signal_obj$signal <- substitute_new_data(series, signal_obj$signal)
    signal_obj$noise <- substitute_new_data(series, signal_obj$noise)

    if (!is.null(names(series))) {
        names(signal_obj$ar_coefs) <- names(series)
    }

    return(signal_obj)
}

#' Tune parameters of HLRA problem.
#'
#' @param series Source time series, numeric vector or list of numeric vectors (may contain NA's)
#' @param r_range Possible values of time series rank.
#' @param p_range Possible values of order of AR process.
#' @param L Window length.
#' @param alpha Parameter for Cadzow weights search
#' @param envelope Optional noise envelope.
#' @param set_seed Seed function for deterministic SVD decomposition
#' @param cluster Optional SNOW cluster for parallel computing
#' @param initial_ar_coefs Initial autoregression coefficients.
#' @param additional_pars Additional parameters for inner optimizers.
#' @return Data frame of class "hlra_tune" containing Bayesian Information Criterion (BIC) for each model.
#' @examples
#' library(Rssa)
#' data("USUnemployment")
#' seedf <- function() set.seed(15)
#' series <- log(USUnemployment[, 2])
#' bic_data <- hlra_tune(series, r_range = 8:12, p_range = 0:3, alpha = .8, initial_ar_coefs = c(.9), set_seed = seedf)
#' plot.hlra_tune(bic_data)
hlra_tune <- function(series, r_range = 1:15, p_range = 0:3,
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
        if (!requireNamespace("snow", quietly = TRUE)) {
            stop("Package \"snow\" needed for cluster computation to work. Please install it.",
                 call. = FALSE)
        }

        snow::clusterCall(cluster, function() {
            library(svd)
            library(Matrix)
            library(rhlra)
        })

        snow::clusterExport(cluster,
                            c("series", "L", "alpha", "envelope", "set_seed",
                              "initial_ar_coefs", "additional_pars"), envir = environment())

        data <- snow::parApply(cluster, input_mat, 1, obtain_bic_df_nonv)
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
