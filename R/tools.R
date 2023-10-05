effective_mask <- function(series) {
    N <- length(series)
    mask <- rep(TRUE, N)

    left_mask <- function(series) cumsum(is.na(series)) == 1:N
    mask[left_mask(series)] <- FALSE
    mask[rev(left_mask(rev(series)))] <- FALSE
    mask
}

effective_length <- function(series) sum(effective_mask(series))

effective_mask_for_weights <- function(series, L) {
    mask <- effective_mask(series)
    mask[!(cumsum(mask) < L)]
}

expand_by_mask <- function(series, signal) {
    series[effective_mask(series)] <- signal
    series
}

sapply_ns <- function(...) sapply(..., simplify = FALSE)

get_rev_row_form <- function(x) {
    trip <- as(x, "TsparseMatrix")

    trip_i <- trip@i
    trip_j <- trip@j
    trip_data <- trip@x

    trip_diagonal <- (trip_i - trip_j) + 1
    trip_pos <- trip_i + 1
    trip <- matrix(0, ncol = max(trip_diagonal),
                        nrow = max(trip_pos))
    trip[cbind(trip_pos, trip_diagonal)] <- trip_data
    trip
}

inv_ac_diags <- function(N, coefs) {
    p <- length(coefs)
    convv <- c(1, -coefs)

    get_small_vec <- function(i) {
        cumsum(convv[(i + 1):length(convv)] * convv[1:(length(convv) - i)])
    }

    small_list <- lapply(0:p, get_small_vec)

    get_big_vec <- function(small_vec) {
        diag_num <- p + 1 - length(small_vec)
        sapply(1:(N - diag_num), function(i)
            small_vec[min(i, N - diag_num - i + 1, length(small_vec))])
    }

    lapply(small_list, get_big_vec)
}

band_mat_from_diags <- function(diags) {
    bandSparse(length(diags[[1]]), k = c(0:(length(diags) - 1)),
               diag = diags, symm=TRUE)
}

get_matrix_weight_matrix <- function(left_diags, right_diags) {
    L <- length(left_diags[[1]])
    K <- length(right_diags[[1]])
    N <- L + K - 1
    pl <- length(left_diags) - 1
    pr <- length(right_diags) - 1

    # i - pl - pr - 1

    answer <- lapply(1:(2 * pl + 2 * pr + 1),
                     function(i) numeric(N - abs(i - pl - pr - 1)))

    outer(-pl:pl, -pr:pr, Vectorize(function(i, j) {
        padding <- numeric((abs(i) + abs(j) - abs(i + j)) / 2)
        answer[[i + j + pl + pr + 1]] <<- answer[[i + j + pl + pr + 1]] +
            c(padding,
              ssa_convolve(left_diags[[abs(i) + 1]], right_diags[[abs(j) + 1]]),
              padding)
        0
    }))

    answer <- answer[(pl + pr + 1) : (2 * pl + 2 * pr + 1)]
    bandSparse(N, k = c(0:(pl + pr)), diag = answer, symm=TRUE)
}

get_splits <- function(series_list) {
    c(0, cumsum(sapply(series_list, length)))
}

split_into_series_list <- function(v, splits) {
    mapply(function(begin, end) v[begin:end],
           begin = splits[-length(splits)] + 1, end = splits[-1],
           SIMPLIFY = FALSE)
}

split_matrix_into_matrix_list <- function(m, splits) {
    mapply(function(begin, end) matrix(m[begin:end, ], ncol = ncol(m)),
           begin = splits[-length(splits)] + 1, end = splits[-1],
           SIMPLIFY = FALSE)
}

ssa_convolve <- function(u, v) {
    l_fft <- hlra_fft(c(u, numeric(length(v) - 1)))
    r_fft <- hlra_fft(c(v, numeric(length(u) - 1)))
    Re(hlra_ifft(l_fft * r_fft))
}

generic_inner_product <- function(x, y, mat) sum(as.numeric(x * as.numeric(mat %*% as.numeric(y))))

generic_inner_product_chol <- function(x, y, cholmat) {
    sum(as.numeric(as.numeric(cholmat %*% as.numeric(x)) * as.numeric(cholmat %*% as.numeric(y))))
}

trmat_indices <- function(N, L) {
    K <- N - L + 1
    i <- rep(1:L, K)
    j <- rep(1:K, each = L)
    matrix(i + j - 1, nrow = L)
}

traj_matrix <- function(series, L) {
    outer(1:L, 1:(length(series) - L + 1),
          function(i, j) series[i + j - 1])
}

null_series <- function(series) {
    series[is.na(series)] <- 0
    series
}

unit_envelope <- function(series, bignumber = 1e6) {
    do_job <- function(series) {
        preans <- rep(1, length(series))
        preans[is.na(series)] <- bignumber
        preans
    }

    if (!is.list(series)) {
        return(do_job(series))
    } else {
        return(sapply_ns(series, do_job))
    }
}

correct_envelope <- function(series, envelope) {
    mask <- is.na(series)
    envelope[mask] <- +Inf
    envelope
}

expand_pars_list <- function(input_list, pars_names, additional_pars) {
    for (name in pars_names) {
        if (name %in% names(additional_pars)) {
            input_list[name] <- additional_pars[name]
        }
    }

    input_list
}

substitute_new_data <- function(original, fix) {
    answer <- original
    if (!is.list(original)) {
        answer[1:length(original)] <- fix
    } else {
        for (i in seq_along(original)) {
            answer[[i]][1:length(original[[i]])] <- fix[[i]]
        }
    }

    answer
}

numeric_check <- function(input) {
    if (is.list(input)) {
        numeric_check_list(input)
    } else {
        numeric_check_vector(input)
    }
}

numeric_check_vector <- function(input) {
    if (!is.numeric(input)) {
        stop("Only numeric (real) input series are supported")
    }

    if (is.matrix(input) || !is.null(dim(input))) {
        stop("Matrix/multidimensional input is now not supported")
    }

    input <- input[!is.na(input)]

    if (length(input) == 0) {
        stop("Input has zero length")
    }

    if (!all(is.finite(input))) {
        stop("Input must have finite values or NA's")
    }
}

numeric_check_list <- function(input) {
    if (length(input) == 0) {
        stop("Input list is empty")
    }

    for (i in seq_along(input)) {
        numeric_check_vector(input[[i]])
    }
}

whole_number_check <- function(input, argname, lower_bound = 1) {
    if (!(length(input) == 1 && as.integer(input) - input == 0
        && input >= lower_bound)) {
        stop(paste("Argument", argname,
            "must be a whole number greater or equal than", lower_bound))
    }
}
