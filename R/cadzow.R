get_oblique_trmat <- function(obj, ...) UseMethod("get_oblique_trmat")

get_oblique_trmat.1d <- function(this, series, L, left_chol, right_chol) {
    left_chol_mat <- as(left_chol, "Matrix")
    right_chol_mat <- as(right_chol, "Matrix")
    matrix(as.matrix(t(left_chol_mat) %*% (traj_matrix(series, L) %*% right_chol_mat)),
           nrow = L)
}

get_oblique_trmat.1dm <- function(this, series, L, left_chol, right_chol) {
    do.call(cbind, lapply(seq_along(series), function(i)
        get_oblique_trmat.1d(this, series[[i]], L, left_chol[[i]], right_chol[[i]])))
}

prepare_ops_for_hankel_svd <- function(obj, ...) UseMethod("prepare_ops_for_hankel_svd")

prepare_ops_for_hankel_svd.1d <- function(this, series, left_chol_mat, right_chol_mat) {
    N <- length(series)
    L <- dim(left_chol_mat)[1]
    K <- dim(right_chol_mat)[1]
    series_fft <- hlra_fft(series)

    mul <- function(v) generic_mul(v, left_chol_mat = left_chol_mat,
                                   right_chol_mat = right_chol_mat, series_fft = series_fft)
    tmul <- function(v) generic_tmul(v, left_chol_mat = left_chol_mat,
                                    right_chol_mat = right_chol_mat, series_fft = series_fft)
    list(N = N, L = L, K = K, mul = mul, tmul = tmul)
}

prepare_ops_for_hankel_svd.1dm <- function(this, series, left_chol_mat, right_chol_mat) {
    N <- sapply(series, length)
    L <- dim(left_chol_mat[[1]])[1]
    Ks <- sapply(right_chol_mat, function(right_chol_mat) dim(right_chol_mat)[1])
    Ksplits <- get_splits(lapply(Ks, numeric))
    K <- sum(Ks)

    series_ffts <- lapply(seq_along(series),
                          function(i) hlra_fft(series[[i]]))
    mul <- function(v) {
        vs <- split_into_series_list(v, Ksplits)
        result <- rowSums(matrix(sapply(seq_along(series), function(i) {
            generic_mul(vs[[i]], left_chol_mat = left_chol_mat[[i]],
                        right_chol_mat = right_chol_mat[[i]], series_fft = series_ffts[[i]])
        }), nrow = L))
        result
    }

    tmul <- function(v) {
        result <- glue_series_lists(lapply(seq_along(series), function(i) {
            generic_tmul(v, left_chol_mat = left_chol_mat[[i]],
            right_chol_mat = right_chol_mat[[i]], series_fft = series_ffts[[i]])}))
        result
    }
    list(N = N, L = L, K = K, mul = mul, tmul = tmul)
}

oblique_hankel_svd <- function(this, series, r, left_chol_mat, right_chol_mat,
                               left_chol, right_chol, high_rank,
                               svd_type = "trlan", overrank = 2) {

    ops_for_hankel <- prepare_ops_for_hankel_svd(this, series, left_chol_mat, right_chol_mat)

    desired_rank <- r

    if (high_rank) {
        desired_rank <- min(ops_for_hankel$L, ops_for_hankel$K)
    }

    trmat <- NULL

    if (svd_type != "svd") {
        trmat <- extmat(ops_for_hankel$mul, ops_for_hankel$tmul,
                        ops_for_hankel$L, ops_for_hankel$K)
    }

    result <- NULL

    lapack_svd <- function() {
        trmat <- get_oblique_trmat(this, series, ops_for_hankel$L,
                                   left_chol, right_chol)

        result <- svd(trmat, nu = desired_rank, nv = desired_rank)
        result$d <- result$d[1:desired_rank]
        result
    }

    keep_matrix <- function(x, rank) {
        if (rank <= 1) {
            x <- matrix(x, ncol = rank)
        }
        x
    }

    propack_svd <- function() {
        svd_rank <- min(ops_for_hankel$L, ops_for_hankel$K,
                        desired_rank + overrank)
        result <- propack.svd(trmat, svd_rank)
        result$d <- result$d[1:desired_rank]
        result$u <- keep_matrix(result$u[, 1:desired_rank], desired_rank)
        result$v <- keep_matrix(result$v[, 1:desired_rank], desired_rank)
        result
    }

    trlan_svd <- function() {
        svd_rank <- min(ops_for_hankel$L, ops_for_hankel$K,
                        desired_rank + overrank)
        result <- trlan.svd(trmat, svd_rank)
        result$d <- result$d[1:desired_rank]
        result$u <- keep_matrix(result$u[, 1:desired_rank], desired_rank)
        result$v <- sapply(1:desired_rank,
            function(i) ops_for_hankel$tmul(result$u[, i]) / result$d[i])
        result
    }

    tryCatch({
        if (svd_type == "svd") {
            result <- lapack_svd()
        } else {
            if (svd_type == "propack") {
                result <- propack_svd()
            } else {
                if (svd_type == "trlan") {
                    result <- trlan_svd()
                } else {
                    stop("Unknown SVD type")
                }
            }
        }
    }, error = function(msg) {
        if (grepl("Please use LAPACK or EISPACK instead.", msg, fixed = TRUE)) {
            cat(paste0("Fallback to SVD due to ", msg, "\n"))
            result <<- lapack_svd()
        } else {
            stop(msg)
        }
    }, warning = function(msg) {
        if (grepl("TRLAN: info->ned (1) is large relative to the matrix dimension",
            msg, fixed = TRUE)) {
            cat(paste0("Fallback to SVD due to ", msg, "\n"))
            result <<- lapack_svd()
        } else {
            warning(msg)
        }
    })

    answer <- list(N = ops_for_hankel$N, L = ops_for_hankel$L,
        d = result$d, u = result$u, v = result$v, r = r)

    if (high_rank) {
        answer$d <- answer$d[(desired_rank - r + 1):desired_rank]
        answer$u <- keep_matrix(answer$u[, (desired_rank - r + 1):desired_rank], r)
        answer$v <- keep_matrix(answer$v[, (desired_rank - r + 1):desired_rank], r)
    }

    answer
}

oblique_hankel_diag_averaging <- function(obj, ...) UseMethod("oblique_hankel_diag_averaging")

oblique_hankel_diag_averaging.1d <- function(this, reslist, left_chol_mat, right_chol_mat,
                                             weights_chol) {
    diag_one_triple <- function(u, v) generic_diag_one_triple(u, v, left_chol_mat, right_chol_mat)

    stable_r <- min(reslist$r, length(reslist$d))

    diag_one_triple_by_index <- function(i) diag_one_triple(
        reslist$u[, i], reslist$v[, i]) * reslist$d[i]

    allsum <- numeric(reslist$N)

    if (stable_r > 0) {
        for (i in 1:stable_r) {
            allsum <- allsum + diag_one_triple_by_index(i)
        }
    }

    as.numeric(solve(weights_chol, allsum, system = "A"))
}

oblique_hankel_diag_averaging.1dm <- function(this, reslist, left_chol_mat, right_chol_mat,
                                              weights_chol) {
    N <- reslist$N
    L <- reslist$L
    Ks <- N - L + 1
    Ksplits <- get_splits(lapply(Ks, numeric))

    vs <- split_matrix_into_matrix_list(reslist$v, Ksplits)

    lapply(seq_along(N), function(i) {
        reslist1d <- list(N = N[[i]], L = L,
                          d = reslist$d, u = reslist$u, v = vs[[i]], r = reslist$r)
        oblique_hankel_diag_averaging.1d(this, reslist1d, left_chol_mat[[i]], right_chol_mat[[i]],
                                         weights_chol[[i]])
    })
}

prepare_cadzow_iterations <- function(obj, ...) UseMethod("prepare_cadzow_iterations")

prepare_cadzow_iterations.1d <- function(this, series, weights_mat) {
    weights_chol <- Cholesky(weights_mat)
    series <- as.numeric(series)
    inner_product <- function(x, y) generic_inner_product(x, y, weights_mat)
    minus <- function(x, y) x - y
    empty_last_steps <- matrix(numeric(0), nrow = length(series))
    empty_answer <- numeric(length(series))
    glue <- identity
    unglue <- identity

    list(weights_chol = weights_chol, series = series, inner_product = inner_product,
        minus = minus, empty_last_steps = empty_last_steps,
        glue = glue, unglue = unglue, empty_answer = empty_answer)
}

prepare_cadzow_iterations.1dm <- function(this, series, weights_mat) {
    weights_chol <- lapply(weights_mat, Cholesky)
    series <- lapply(series, as.numeric)
    splits <- get_splits(series)
    inner_product <- function(x, y) sum(mapply(generic_inner_product, x, y, weights_mat))
    minus <- function(x, y) mapply("-", x, y, SIMPLIFY = FALSE)
    empty_last_steps <- matrix(numeric(0), nrow = sum(sapply(series, length)))
    empty_answer <- sapply(lapply(series, length), numeric)
    glue <- function(x) glue_series_lists(x)
    unglue <- function(x) split_into_series_list(x, splits)

    list(weights_chol = weights_chol, series = series, inner_product = inner_product,
        minus = minus, empty_last_steps = empty_last_steps,
        glue = glue, unglue = unglue, empty_answer = empty_answer)
}

cadzow_iterations <- function(this, series, r, left_chol_mat, right_chol_mat,
                               left_chol, right_chol,
                               weights_mat, cadzow_epsilon = 1e-6,
                               cadzow_it_limit = 100, debug = FALSE,
                               set_seed = NULL, sylvester_nulling = NULL,
                               high_rank = FALSE, ...) {

    if (!is.null(set_seed)) {
        set_seed()
    }

    list2env(prepare_cadzow_iterations(this, series, weights_mat), environment())

    it <- 0
    prev <- series
    proj <- series

    dists <- numeric(0)

    full_norm2 <- inner_product(series, series)

    stop_criterion <- function(new, old) {
        change_norm2 <- inner_product(minus(new, old), minus(new, old))
        dists <<- c(dists, sqrt(change_norm2 / full_norm2))
        (change_norm2 / full_norm2) > cadzow_epsilon^2
    }

    while ((it == 0 || stop_criterion(proj, prev)) && (it < cadzow_it_limit) && r > 0) {
        it <- it + 1
        prev <- proj
        proj <- oblique_hankel_diag_averaging(this, oblique_hankel_svd(this, prev,
                                              r, left_chol_mat, right_chol_mat, left_chol,
                                              right_chol, high_rank, ...),
                                              left_chol_mat, right_chol_mat, weights_chol)

        if (high_rank) {
            proj <- minus(prev, proj)
        }

        if (!is.null(sylvester_nulling)) {
            proj <- lapply(seq_along(proj), function(i) {
                series <- proj[[i]]
                null_len <- sylvester_nulling[i]
                if (null_len > 0) {
                    series[1:null_len] <- 0
                    series[(length(series) - null_len + 1):length(series)] <- 0
                }
                series
            })
        }
    }

    if (r == 0) {
        proj <- empty_answer
    }

    if (debug == "P" || debug) {
        cat(sprintf("%d cadzow iterations done\n", it))
    }

    if (debug == "P") {
        if (length(dists) > 1) {
            plot(dists, type = "b", log = "y",
                 main = "Relative step lengths")
        }
    }

    list(signal = proj, it = it, noise = minus(series, proj))
}

prepare_cadzow <- function(obj, ...) UseMethod("prepare_cadzow")

prepare_cadzow.1d <- function(this, left_diags, right_diags) {
    left_mat <- band_mat_from_diags(left_diags)
    left_chol <- Cholesky(left_mat, perm = FALSE, LDL = FALSE)
    left_chol_mat <- as(left_chol, "Matrix")

    right_mat <- band_mat_from_diags(right_diags)
    right_chol <- Cholesky(right_mat, perm = FALSE, LDL = FALSE)
    right_chol_mat <- as(right_chol, "Matrix")

    left_chol_mat <- get_rev_row_form(left_chol_mat)
    right_chol_mat <- get_rev_row_form(right_chol_mat)

    weights <- get_matrix_weight_matrix(left_diags, right_diags)

    list(left_mat = left_mat,
         left_chol = left_chol,
         left_chol_mat = left_chol_mat,
         right_mat = right_mat,
         right_chol = right_chol,
         right_chol_mat = right_chol_mat,
         weights = weights)
}

prepare_cadzow.1dm <- function(this, left_diags, right_diags) {
    left_mat <- lapply(left_diags, band_mat_from_diags)
    left_chol <- lapply(left_mat, function(x) Cholesky(x, perm = FALSE, LDL = FALSE))
    left_chol_mat <- lapply(left_chol, function(x) as(x, "Matrix"))

    right_mat <- lapply(right_diags, band_mat_from_diags)
    right_chol <- lapply(right_mat, function(x) Cholesky(x, perm = FALSE, LDL = FALSE))
    right_chol_mat <- lapply(right_chol, function(x) as(x, "Matrix"))

    left_chol_mat <- lapply(left_chol_mat, get_rev_row_form)
    right_chol_mat <- lapply(right_chol_mat, get_rev_row_form)

    weights <- mapply(function(left_diags, right_diags) get_matrix_weight_matrix(left_diags, right_diags),
                      left_diags, right_diags)

    list(left_mat = left_mat,
         left_chol = left_chol,
         left_chol_mat = left_chol_mat,
         right_mat = right_mat,
         right_chol = right_chol,
         right_chol_mat = right_chol_mat,
         weights = weights)
}

cadzow <- function(this, series, r, left_diags, right_diags,
                   debug = FALSE, set_seed = NULL, additional_pars) {

    list2env(prepare_cadzow(this, left_diags, right_diags), environment())

    cadzow_call_list = list(this = this,
                            series = series,
                            r = r,
                            left_chol_mat = left_chol_mat,
                            right_chol_mat = right_chol_mat,
                            weights_mat = weights,
                            set_seed = set_seed,
                            debug = debug,
                            left_chol = left_chol,
                            right_chol = right_chol)

    cadzow_data <- do.call(cadzow_iterations,
                           expand_pars_list(cadzow_call_list, cadzow_add_pars_names, additional_pars))

    answer <- cadzow_data

    answer
}

cadzow_add_pars_names <- c("cadzow_epsilon", "cadzow_it_limit", "svd_type", "overrank")
