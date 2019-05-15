simple_expm <- function(m, pow) {
    if (pow == 0) {
        return(diag(dim(m)[1]))
    }

    if (pow == 1) {
        return(m)
    }

    x <- Recall(m, pow %/% 2)
    if (pow %% 2 == 1) {
        return(x %*% x %*% m)
    } else {
        return(x %*% x)
    }
}

hlra_igapfill <- function(series, r, eps = 1e-6, scheme_order = 4, 
                        error_bound = 1e-2, 
                        it_limit = 100, debug = TRUE, set_seed = NULL) {
    if (!is.null(set_seed)) {
        set_seed()
    }

    mask <- is.na(series)
    mask <- (1:length(series))[mask]
    
    scheme_order <- min(scheme_order, length(mask) - 1)
    
    signal <- series
    
    # print(mask)
    
    x <- rep(median(series, na.rm = TRUE), length(mask))
    prev <- NA
    series[mask] <- x
    dists <- numeric()
    dist <- 0
    step <- NA


    last_steps <- matrix(numeric(0), nrow = length(mask))
    # all_its <- matrix(x, nrow = length(mask))

    it <- 0

    while (it < it_limit && (is.na(prev) || sum(abs(prev - x)) > eps)) {
        # print(it)
        # print(series)
        pseudoobj <- list()
        class(pseudoobj) <- "1d"
        L <- default_L(series)

        series <- hlra_ssa(pseudoobj, series, L, r, debug, set_seed)$signal
        # print(series)
        prev <- x
        x <- series[mask]
        step <- x - prev
        # print(step)
        # plot(step, type = "b")
        series <- signal
        last_steps <- cbind(last_steps, step)

        if (scheme_order > 0 && ncol(last_steps) == scheme_order + 1) {
            qrobj <- qr(last_steps[, 1:scheme_order])
            # coefs <- qr.coef(qrobj, step)
            # stepmat <- cbind(rbind(numeric(scheme_order - 1), diag(scheme_order - 1)), coefs)
            # print(coefs)
            
            stepmat <- as.matrix(qr.coef(qrobj, last_steps[, -1]))
            
            gammas <- eigen(stepmat, FALSE, TRUE)$values
            # print(stepmat)
            # print(gammas)
            gamma <- max(Mod(gammas))
            
            basic_error <- sqrt(sum((qr.resid(qrobj, step))^2))
            
            prev_norm <- sqrt(sum(last_steps[, scheme_order]^2))
            # print(error_norm)
            if (gamma < 1) {
                ls <- 0
                rs <- 2^14
                for (i in 1:14) {
                    n <- (ls + rs) / 2;
                    approx_error_est <- basic_error * ((n + 1) - gamma * (1 - gamma^(n + 1))/(1 - gamma))/(1 - gamma)
                    track_length_est <- prev_norm * (1 - gamma^(n + 1))/(1 - gamma)
                    
                    if (approx_error_est > track_length_est * error_bound) {
                        rs <- n;
                    } else {
                        ls <- n;
                    }
                }
                # print(n)
                v <- solve((stepmat - diag(scheme_order)), (simple_expm(stepmat, n + 1) - diag(scheme_order)))
                v <- (last_steps[, -1] %*% v)[, scheme_order]
                x <- prev + v
                last_steps <- matrix(numeric(0), nrow = length(mask))
            } else {
                last_steps <- last_steps[, -1]
            }
        }
        # all_its <- cbind(all_its, x)
        # print(x)
        # print(mask)
        # print(series)
        series[mask] <- x
        dist <- sum(abs(prev - x))
        dists <- c(dists, dist)
        it <- it + 1
    }

    if (debug) {
        plot(series, type = "b", main = "igapfill: filled series", xlab = "index", ylab = "value")
        plot(dists, log = "y", type = "b", main = "igapfill: steps", xlab = "index", ylab = "value")
    }
    series
}

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
    trip <- as(x, "dgTMatrix")

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
    Ksplits <- get_splits(sapply_ns(Ks, numeric))
    K <- sum(Ks)

    series_ffts <- sapply_ns(seq_along(series),
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
        result <- glue_series_lists(sapply_ns(seq_along(series), function(i) {
            generic_tmul(v, left_chol_mat = left_chol_mat[[i]],
            right_chol_mat = right_chol_mat[[i]], series_fft = series_ffts[[i]])}))
        result
    }
    list(N = N, L = L, K = K, mul = mul, tmul = tmul)
}

hankel_svd_double <- function(this, series, r, left_chol_mat, right_chol_mat,
                              svd_type = "trlan") {

    ops_for_hankel <- prepare_ops_for_hankel_svd(this, series, left_chol_mat, right_chol_mat)

    mymat <- extmat(ops_for_hankel$mul, ops_for_hankel$tmul, 
        ops_for_hankel$L, ops_for_hankel$K)
    result <- NULL

    lapack_svd <- function() {
        mymat <- sapply(1:ops_for_hankel$K, function(i) {
                uv <- numeric(ops_for_hankel$K)
                uv[i] <- 1
                ops_for_hankel$mul(uv)
                })
        result <- svd(mymat, nu = r, nv = r)
        result$d <- result$d[1:r]
        result
    }

    tryCatch({
        if (svd_type == "svd") {
            result <<- lapack_svd()
        } else {
            if (svd_type == "propack") {
                result <- propack.svd(mymat, r)
            } else {
                result <- trlan.svd(mymat, r)
                V <- sapply(1:r, function(i) ops_for_hankel$tmul(result$u[, i]) / result$d[i])
                result$v <- V
            }
        }
    }, error = function(msg) {
        if (grepl('Please use LAPACK or EISPACK instead.', msg, fixed=TRUE)) {
            cat(paste0("Fallback to SVD due to ", msg, "\n"))
            result <<- lapack_svd()
        } else {
            stop(msg)
        }
    }, warning = function(msg) {
        if (grepl('TRLAN: info->ned (1) is large relative to the matrix dimension', msg, fixed=TRUE)) {
            cat(paste0("Fallback to SVD due to ", msg, "\n"))
            result <<- lapack_svd()
        } else {
            warning(msg)
        }        
    })

    list(N = ops_for_hankel$N, L = ops_for_hankel$L, 
        d = result$d, u = result$u, v = result$v, r = r)
}

ssa_convolve <- function(u, v) {
    l_fft <- hlra_fft(c(u, numeric(length(v) - 1)))
    r_fft <- hlra_fft(c(v, numeric(length(u) - 1)))
    Re(hlra_ifft(l_fft * r_fft))
}

prepare_diag_one_triple <- function(obj, ...) UseMethod("prepare_diag_one_triple")

prepare_diag_one_triple.1d <- function(this, reslist, left_chol_mat, right_chol_mat) {
    N <- reslist$N
    L <- reslist$L
    K <- N - L + 1

    function(u, v) generic_diag_one_triple(u, v, left_chol_mat, right_chol_mat)
}

prepare_diag_one_triple.1dm <- function(this, reslist, left_chol_mat, right_chol_mat) {
    N <- reslist$N
    L <- reslist$L
    Ks <- N - L + 1
    splits <- get_splits(sapply_ns(N, numeric))
    Ksplits <- get_splits(sapply_ns(Ks, numeric))

    function(u, v) {
        vs <- split_into_series_list(v, Ksplits)
        glue_series_lists(sapply_ns(seq_along(N), function(i) {
            generic_diag_one_triple(u, vs[[i]],
                                    left_chol_mat = left_chol_mat[[i]],
                                    right_chol_mat = right_chol_mat[[i]])
        }))
    }
}

hankel_reweight <- function(obj, ...) UseMethod("hankel_reweight")

hankel_reweight.1d <- function(this, weights_chol, allsum, N) {
    as.numeric(solve(weights_chol, allsum, system = "A"))
}

hankel_reweight.1dm <- function(this, weights_chol, allsum, N) {
    splits <- get_splits(sapply_ns(N, numeric))
    allsums <- split_into_series_list(allsum, splits)

    sapply_ns(seq_along(N), function(i) {
            as.numeric(solve(weights_chol[[i]], allsums[[i]], system = "A"))})
}

hankel_diag_average_double <- function(this, reslist, left_chol_mat, right_chol_mat,
                                       weights_chol) {
    diag_one_triple <- prepare_diag_one_triple(this, reslist, left_chol_mat, right_chol_mat)

    stable_r <- min(reslist$r, length(reslist$d))

    diag_one_triple_by_index <- function(i) diag_one_triple(
        reslist$u[, i], reslist$v[, i]) * reslist$d[i]

    allsum <- rowSums(sapply(1:stable_r,
                             diag_one_triple_by_index))

    hankel_reweight(this, weights_chol, allsum, reslist$N)
}

generic_inner_product <- function(x, y, mat) sum(as.numeric(x * as.numeric(mat %*% as.numeric(y))))

prepare_oblique_cadzow_eps <- function(obj, ...) UseMethod("prepare_oblique_cadzow_eps")

prepare_oblique_cadzow_eps.1d <- function(this, series, weights_mat) {
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

prepare_oblique_cadzow_eps.1dm <- function(this, series, weights_mat) {
    weights_chol <- sapply_ns(weights_mat, Cholesky)
    series <- sapply_ns(series, as.numeric)
    splits <- get_splits(series)
    inner_product <- function(x, y) sum(mapply(generic_inner_product, x, y, weights_mat))
    minus <- function(x, y) mapply("-", x, y, SIMPLIFY = FALSE)
    empty_last_steps <- matrix(numeric(0), nrow = sum(sapply(series, length)))
    empty_answer <- sapply(sapply_ns(series, length), numeric)
    glue <- function(x) glue_series_lists(x)
    unglue <- function(x) split_into_series_list(x, splits)

    list(weights_chol = weights_chol, series = series, inner_product = inner_product,
        minus = minus, empty_last_steps = empty_last_steps,
        glue = glue, unglue = unglue, empty_answer = empty_answer)
}

oblique_cadzow_eps <- function(this, series, r, left_chol_mat, right_chol_mat,
                           weights_mat, epsilon = 1e-6,
                           it_limit = 100, debug = FALSE,
                           scheme_order = 4, error_bound = 1e-2,
                           error_norm_rel = 1e-1, ...) {

    list2env(prepare_oblique_cadzow_eps(this, series, weights_mat), environment())

    it <- 0
    prev <- series
    proj <- series

    dists <- numeric(0)

    last_steps <- empty_last_steps

    stop_criterion <- function(new, old) {
        change_norm2 <- inner_product(minus(new, old), minus(new, old))
        full_norm2 <- inner_product(new, new)
        dists <<- c(dists, sqrt(change_norm2 / full_norm2))
        (change_norm2 / full_norm2) > epsilon^2
    }
    
    while ((it == 0 || stop_criterion(proj, prev)) && it < it_limit && r > 0) {
        it <- it + 1
        prev <- proj
        proj <- hankel_diag_average_double(this, hankel_svd_double(this, prev,
                                            r, left_chol_mat, right_chol_mat, ...),
                                           left_chol_mat, right_chol_mat, weights_chol)
        step <- glue(minus(proj, prev))

        last_steps <- cbind(last_steps, step)

        if (scheme_order > 0 && ncol(last_steps) == scheme_order + 1) {
            qrobj <- qr(last_steps[, 1:scheme_order])
            stepmat <- as.matrix(qr.coef(qrobj, last_steps[, -1]))
            
            gammas <- eigen(stepmat, FALSE, TRUE)$values
            gamma <- max(Mod(gammas))
            
            basic_error <- sqrt(sum((qr.resid(qrobj, step))^2))
            
            prev_norm <- sqrt(sum(last_steps[, scheme_order]^2))
            
            if (gamma < 1 && basic_error/prev_norm < error_norm_rel) {
                ls <- 0
                rs <- 2^14
                for (i in 1:14) {
                    n <- (ls + rs) / 2;
                    approx_error_est <- basic_error * ((n + 1) - gamma * (1 - gamma^(n + 1))/(1 - gamma))/(1 - gamma)
                    track_length_est <- prev_norm * (1 - gamma^(n + 1))/(1 - gamma)
                    
                    if (approx_error_est > track_length_est * error_bound) {
                        rs <- n;
                    } else {
                        ls <- n;
                    }
                }

                v <- solve((stepmat - diag(scheme_order)), (simple_expm(stepmat, n + 1) - diag(scheme_order)))
                v <- (last_steps[, -1] %*% v)[, scheme_order]

                last_steps <- empty_last_steps
                proj <- unglue(glue(prev) + v)

            } else {
                last_steps <- last_steps[, -1]
            }
        }

    }

    if (r == 0) {
        proj <- empty_answer
    }

    if (debug) {
        cat(sprintf("%d cadzow iterations done\n", it))
        if (length(dists) > 1) {
            plot(dists, type = "b", log = "y",
                 main = "Relative residual eigenvalues")
        }
    }

    list(signal = proj, it = it, noise = minus(series, proj))
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

neg_fourier_short_mat <- function(N, indices) {
    full_indx <- outer(0:(N-1), indices)
    matrix(complex(argument = (full_indx %% N)*2*pi/N), nrow = N)
}

get_comp_space_by_v <- function(this, N, v, v_2 = FALSE, return_pack = FALSE) {
    if (!v_2) {
        result <- eval_basis(this, N, v)

        if (return_pack) {
            return(result)
        } else {
            return(result$basis)
        }
    }
    if (v_2) {
        result <- eval_tangent_basis(this, N, v)

        if (return_pack) {
            stop("Unsupported")
        }

        return(result)
    }
}

get_pseudograd <- function(this, N, vspace_pack, signal, j) {
    eval_pseudograd(this, N, vspace_pack$glrr, signal, j, vspace_pack)
}

get_space_by_v <- function(this, N, v, v_2 = FALSE) {
    pre_answer <- get_comp_space_by_v(this, N, v, v_2 = v_2)
    if (v_2) {
        v <- ssa_convolve(v, v)
    }
    r <- length(v) - 1
    j <- which.max(abs(v))
    begin_indices <- numeric(0)
    if (j > 1) begin_indices <- 1:(j - 1)
    end_indices <- numeric(0)
    if (j < r + 1) end_indices <- (N - r + j):N
    a_coefs <- solve(pre_answer[c(begin_indices, end_indices), ])
    Re(pre_answer %*% a_coefs)
}

get_v_from_nonlrf_series <- function(obj, ...) UseMethod("get_v_from_nonlrf_series")

get_v_from_nonlrf_series.1d <- function(this, series, r) {
    L <- r + 1
    N <- length(series)
        trmat <- traj_matrix(series[!is.na(series)], L)
    svd(trmat)$u[, L]
}

get_v_from_nonlrf_series.1dm <- function(this, series, r) {
    L <- r + 1
    Ns <- sapply(series, length)
    trmat <- do.call(cbind,
         sapply_ns(series, function(x) traj_matrix(x[!is.na(x)], L)))
    svd(trmat)$u[, L]
}

# weighted_project_onto_vspace <- function(series, vspace, weights) {
#     real_part <- as.matrix(weights %*% Re(vspace))
#     imag_part <- as.matrix(weights %*% Im(vspace))
#     wpvs <- real_part + 1i * imag_part
#     as.numeric(Re(vspace %*% solve(t(Conj(vspace)) %*% wpvs,
#                                    t(Conj(vspace)) %*% as.numeric(weights %*% series))))
# }

weighted_project_rotate_basis <- function(vspace, chol_weights) {
    real_part <- as.matrix(chol_weights %*% Re(vspace))
    imag_part <- as.matrix(chol_weights %*% Im(vspace))
    vspace_rot <- real_part + 1i * imag_part
    qr(vspace_rot)
}

weighted_project_onto_vspace_qr <- function(series, qrobj, vspace, chol_weights) {
    series_rot <- as.numeric(chol_weights %*% series)

    coef <- qr.coef(qrobj, series_rot)

    Re(vspace %*% coef)
}

weighted_project_onto_vspace_coef_rot <- function(obj, ...) UseMethod("weighted_project_onto_vspace_coef_rot")

weighted_project_onto_vspace_coef_rot.1d <- function(this, series, vspace, chol_weights) {
    list(vspace_rot = as.matrix(chol_weights %*% vspace),
        series_rot = as.numeric(chol_weights %*% series))
}

weighted_project_onto_vspace_coef_rot.1dm <- function(this, series, vspace, chol_weights) {
    list(vspace_rot = do.call(rbind, sapply_ns(seq_along(series),
                    function(i) as.matrix(chol_weights[[i]] %*% vspace[[i]]))),
        series_rot = glue_series_lists(sapply_ns(seq_along(series),
                    function(i) as.numeric(chol_weights[[i]] %*% series[[i]]))))
}

weighted_project_onto_vspace_coef <- function(this, series, vspace, chol_weights, tol = 1e-14) {
    list2env(weighted_project_onto_vspace_coef_rot(this, series, vspace, chol_weights), environment())

    # qrobj <-qr(vspace_rot, LAPACK = TRUE)
    # print(diag(qr.R(qrobj)))
    # print(qrobj)
    #
    # print(as.numeric(qr.coef(qrobj, series_rot)))

    svdobj <- svd(vspace_rot)

    svdobj$d[svdobj$d < svdobj$d[1] * tol] <- Inf
    # print(svdobj$d)
    svdobj$v %*% ((t(svdobj$u) %*% series_rot) / svdobj$d)
}

prepare_find_step <- function(obj, ...) UseMethod("prepare_find_step")

prepare_find_step.1d <- function(this, signal, series, r, j, vspace_pack, weights_chol) {
    noise <- series - signal
    N <- length(series)

    K <- N - r

    pseudograd <- Re(get_pseudograd(this, N, vspace_pack, signal, j))
    pseudograd_minus <- sapply(1:r, function(i)
        weighted_project_onto_vspace_qr(Re(pseudograd[, i]), vspace_pack$qrobj,
                                        vspace_pack$basis, weights_chol))
    pseudograd <- pseudograd - pseudograd_minus

    list(pseudograd = pseudograd, noise = noise)
}

prepare_find_step.1dm <- function(this, signal, series, r, j, vspace_pack, weights_chol) {
    noise <- mapply("-", series, signal, SIMPLIFY = FALSE)
    Ns <- sapply(series, length)
    Ks <- Ns - r
    pseudograd <- sapply_ns(seq_along(series), function(k) {
        pseudograd_cur <- Re(get_pseudograd(this, Ns[[k]], vspace_pack[[k]], signal[[k]], j))
        pseudograd_minus <- sapply(1:r, function(i)
            weighted_project_onto_vspace_qr(Re(pseudograd_cur[, i]), vspace_pack[[k]]$qrobj,
                                            vspace_pack[[k]]$basis, weights_chol[[k]]))
        pseudograd_cur - pseudograd_minus
    })

    list(pseudograd = pseudograd, noise = noise)
}

find_step <- function(this, signal, series, v, vspace_pack, weights_chol, debug = FALSE, ...) {
    r <- length(v) - 1
    j <- which.max(abs(v))

    prepare_obj <- prepare_find_step(this, signal, series, r, j, vspace_pack, weights_chol)

    used_coefs <- weighted_project_onto_vspace_coef(this, prepare_obj$noise, 
        prepare_obj$pseudograd, weights_chol)

    ans <- numeric(r+1)

    ans[-j] <- used_coefs

    ans
}

null_series <- function(series) {
    series[is.na(series)] <- 0
    series
}

generic_get_vspace_pack <- function(this, N, v, weights_chol) {
    result <- get_comp_space_by_v(this, N, v, return_pack = TRUE)
    result$qrobj <- weighted_project_rotate_basis(result$basis, weights_chol)
    result
}

generic_get_signal <- function(this, x, vspace_pack, weights_chol)
    weighted_project_onto_vspace_qr(x, vspace_pack$qrobj, vspace_pack$basis, weights_chol)

prepare_mgn <- function(obj, ...) UseMethod("prepare_mgn")

prepare_mgn.1d <- function(this, series, signal, r, weights, weights_chol) {
    series <- as.numeric(series)

    if (any(is.na(series))) {
        series <- null_series(series)
    }

    if (is.null(weights_chol)) {
        weights_chol <- chol(weights)
    }
    N <- length(series)
    get_vspace_pack <- generic_get_vspace_pack

    get_signal <- generic_get_signal

    inner_product <- function(x, y) generic_inner_product(x, y, weights)
    minus <- function(x, y) x - y
    plus <- function(x, y) x + y
    mult <- function(x, y) x * y

    list(series = series, get_vspace_pack = get_vspace_pack, get_signal = get_signal,
        inner_product = inner_product, minus = minus, plus = plus, mult = mult, N = N)
}

prepare_mgn.1dm <- function(this, series, signal, r, weights, weights_chol) {
    series <- sapply_ns(series, as.numeric)

    for (i in seq_along(series)) {
        if (any(is.na(series[[i]]))) {
            series[[i]] <- null_series(series[[i]])
        }
    }

    splits <- get_splits(series)
    if (is.null(weights_chol)) {
        weights_chol <- sapply_ns(weights, chol)
    }
    N <- sapply(series, length)
    get_vspace_pack <- function(this, Ns, v, weights_chol)
        sapply_ns(seq_along(Ns),
                  function(i) generic_get_vspace_pack(this, Ns[i], v, weights_chol[[i]]))
    get_signal <- function(this, x, vspace_pack, weights_chol)
        sapply_ns(seq_along(x),
                  function(i) generic_get_signal(this, x[[i]], vspace_pack[[i]], weights_chol[[i]]))

    inner_product <- function(x, y) sum(mapply(generic_inner_product, x, y, weights))
    minus <- function(x, y) mapply("-", x, y, SIMPLIFY = FALSE)
    plus <- function(x, y) mapply("+", x, y, SIMPLIFY = FALSE)
    mult <- function(x, y) mapply("*", x, y, SIMPLIFY = FALSE)

    list(series = series, get_vspace_pack = get_vspace_pack, get_signal = get_signal,
    inner_product = inner_product, minus = minus, plus = plus, mult = mult, N = N)
}


mgn <- function(this, series, signal, r, weights,
                               mgn_search_threshold = 1e-7,
                               mgn_it_limit = 100,
                               mgn_j_limit = 4, debug = FALSE,
                               max_step = 1,
                               v_begin = NULL, weights_chol = NULL, ...) {
    #MGN
    cur_v <- NA
    best_signal <- NA
    vspace_pack <- NA

    source_series <- series

    list2env(prepare_mgn(this, series, signal, r, weights, weights_chol), environment())


    initial_best_signal <- NA

    if (is.null(v_begin)) {
        cur_v <- get_v_from_nonlrf_series(this, signal, r)
    } else {
        cur_v <- v_begin
    }

    vspace_pack <- get_vspace_pack(this, N, cur_v, weights_chol)
    best_signal <- get_signal(this, series, vspace_pack, weights_chol)

    initial_best_signal <- best_signal

    if (any(is.na(cur_v))) return(signal)
    if (any(is.na(best_signal))) return(signal)

    no_localsearch_flag <- FALSE

    stop_criterion <- function(signal, step, prev_step) {
        #step-to-solution ratio
        coef <- sqrt(inner_product(minus(best_signal, signal), minus(best_signal, signal))/
            inner_product(best_signal, best_signal))

        if (coef > mgn_search_threshold) {
            #closer to the solution than previous best_signal
            return(inner_product(minus(best_signal, signal), 
                plus(signal, minus(best_signal, mult(series, 2)))) > 0)
        }
        else {
            if (sum(step^2) < sum(prev_step^2)) {
                return(TRUE)
            } else {
                no_localsearch_flag <<- TRUE
                return(FALSE)
            }
        }
    }

    difference_between_two_series <- function(x, y) {
        sqrt(inner_product(minus(x, y), minus(x, y)))
    }

    best_signal_distance <- function(signal) difference_between_two_series(signal, best_signal)

    noise_distance <- function(signal) difference_between_two_series(signal, series)


    it <- 0

    dists <- noise_distance(best_signal)
    prev_step <- max_step

    while(TRUE) {
        step <- NA
        tryCatch({step <- find_step(this, signal = best_signal, series = series, v = cur_v,
                                        vspace_pack = vspace_pack, weights_chol = weights_chol,
                                        debug = debug, ...)},
                 warning = function(x) {print(x)})

        if (any(is.na(step))) break

        alpha <- max_step
        j <- 0
        next_it <- FALSE
        while (j < mgn_j_limit) {
            new_v <- cur_v + step * alpha

            new_vspace_pack <- get_vspace_pack(this, N, new_v, weights_chol)
            new_signal <- get_signal(this, series, new_vspace_pack, weights_chol)

            if (stop_criterion(new_signal, step * alpha, prev_step)) {
                next_it <- TRUE
                it <- it + 1
                dists <- c(dists, noise_distance(new_signal))
                best_signal <- new_signal

                prev_step <- step * alpha
                cur_v <- new_v
                vspace_pack <- new_vspace_pack
                break
            }
            if (no_localsearch_flag) {
                j <- mgn_j_limit
                break
            }
            alpha <- alpha / 2
            j <- j + 1
        }

        if (j == mgn_j_limit && debug ) {
            cat("no steps found\n")
        } else {
            if (j > 0 && debug) cat(sprintf("reduction, it = %d, j = %d\n", it, j))
        }

        if (!next_it) break
        if (it > mgn_it_limit) break
    }

    if (debug) {
        cat(sprintf("%d MGN iterations done\n", it))
        draw_full_plot <- function(N, series, best_signal, signal, initial_best_signal) {
            matplot(1:N, cbind(series, best_signal, signal),
                    type = "l", lty = 1:3,
                    col = c("black", "blue", "red"), main = "Approximations")
            legend("topleft", legend = c("Source series", "MGN", "Initial"),
                   col = c("black", "blue", "red"), lty = 1:3)
        }

        if (!is.list(series)) {
            draw_full_plot(N, source_series, best_signal, signal, initial_best_signal)
        } else {
            mapply(draw_full_plot, N, source_series, best_signal, signal, initial_best_signal)
        }

        plot(dists/sqrt(sum(N)), type = "b", main = "Fisher Scoring Distance", xlab = "iteration",
            ylab = "Mean distance")

    }

    # stop()
    list(signal = best_signal, v = cur_v, it = it, dists = dists, 
        noise = minus(series, best_signal), 
        noise_norm = inner_product(minus(series, best_signal), minus(series, best_signal)))
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

fill_gaps <- function(series, r, debug = FALSE, set_seed = NULL) {
    result <- series
    if (any(is.na(series))) {
        if (debug) cat("fillgaps... ")
        result <- hlra_igapfill(series, r, set_seed = set_seed, debug = debug)
        if (debug) cat("done\n")
    }
    result
}

correct_envelope <- function(series, envelope) {
    mask <- is.na(series)
    envelope[mask] <- +Inf
    envelope
}

prepare_cadzow_with_mgn <- function(obj, ...) UseMethod("prepare_cadzow_with_mgn")

prepare_cadzow_with_mgn.1d <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope) {
    if (is.null(series_for_cadzow)) series_for_cadzow <- series
    

    if (any(is.na(series))) {
        envelope <- correct_envelope(series, envelope)            
    }

    left_diags <- inv_ac_diags(L, coefs)
    left_mat <- band_mat_from_diags(left_diags)
    left_chol <- Cholesky(left_mat, perm = FALSE, LDL = FALSE)
    left_chol_mat <- as(left_chol, "Matrix")

    right_diags <- list(right_diag[effective_mask_for_weights(series, L)])
    right_mat <- band_mat_from_diags(right_diags)
    right_chol <- Cholesky(right_mat, perm = FALSE, LDL = FALSE)
    right_chol_mat <- as(right_chol, "Matrix")

    left_chol_mat <- get_rev_row_form(left_chol_mat)
    right_chol_mat <- get_rev_row_form(right_chol_mat)

    weights <- get_matrix_weight_matrix(left_diags, right_diags)

    common_expand_by_mask <- function(x) expand_by_mask(series, x)

    envelope_d <- Diagonal(length(series), 1/envelope)
    diags <- inv_ac_diags(length(series), coefs)

    ideal_weights <- envelope_d %*% band_mat_from_diags(diags) %*% envelope_d
    weights_chol <- chol(band_mat_from_diags(diags)) %*% envelope_d

    n <- length(series)
    df <- 2 * r + length(coefs) + 1

    list(series_for_cadzow = series_for_cadzow,
        envelope = envelope,
        left_diags = left_diags,
        left_mat = left_mat,
        left_chol = left_chol,
        left_chol_mat = left_chol_mat,
        right_diags = right_diags,
        right_mat = right_mat,
        right_chol = right_chol,
        right_chol_mat = right_chol_mat,
        weights = weights,
        common_expand_by_mask = common_expand_by_mask,
        ideal_weights = ideal_weights,
        weights_chol = weights_chol,
        n = n,
        df = df)
}

prepare_cadzow_with_mgn.1dm <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope) {
    if (is.null(series_for_cadzow)) series_for_cadzow <- series

    for (i in seq_along(series)) {
        if (any(is.na(series[[i]]))) {
            envelope[[i]] <- correct_envelope(series[[i]], envelope[[i]])
        }
    }

    left_diags <- sapply_ns(coefs, function(x) inv_ac_diags(L, x))
    left_mat <- sapply_ns(left_diags, band_mat_from_diags)
    left_chol <- sapply_ns(left_mat, function(x) Cholesky(x, perm = FALSE, LDL = FALSE))
    left_chol_mat <- sapply_ns(left_chol, function(x) as(x, "Matrix"))

    right_diags <- sapply_ns(seq_along(series), function(i)
        list(right_diag[[i]][effective_mask_for_weights(series[[i]], L)]))
    right_mat <- sapply_ns(right_diags, band_mat_from_diags)
    right_chol <- sapply_ns(right_mat, function(x) Cholesky(x, perm = FALSE, LDL = FALSE))
    right_chol_mat <- sapply_ns(right_chol, function(x) as(x, "Matrix"))

    left_chol_mat <- sapply_ns(left_chol_mat, get_rev_row_form)
    right_chol_mat <- sapply_ns(right_chol_mat, get_rev_row_form)

    weights <- mapply(function(left_diags, right_diags) get_matrix_weight_matrix(left_diags, right_diags),
                      left_diags, right_diags)

    common_expand_by_mask <- function(x) mapply(expand_by_mask, series, x, SIMPLIFY = FALSE)

    envelope_d <- sapply_ns(seq_along(series),
                         function(i) Diagonal(length(series[[i]]), 1/envelope[[i]]))
    diags <- sapply_ns(seq_along(series),
                    function(i) inv_ac_diags(length(series[[i]]), coefs[[i]]))
    ideal_weights <- sapply(seq_along(series), function(i)
        envelope_d[[i]] %*% band_mat_from_diags(diags[[i]]) %*% envelope_d[[i]])
    weights_chol <- sapply(seq_along(series), function(i)
        chol(band_mat_from_diags(diags[[i]])) %*% envelope_d[[i]])

    n <- sum(sapply(series, length))
    df <- (1 + length(series)) * r + length(coefs[[1]]) * (length(series) + 1)

    list(series_for_cadzow = series_for_cadzow,
        envelope = envelope,
        left_diags = left_diags,
        left_mat = left_mat,
        left_chol = left_chol,
        left_chol_mat = left_chol_mat,
        right_diags = right_diags,
        right_mat = right_mat,
        right_chol = right_chol,
        right_chol_mat = right_chol_mat,
        weights = weights,
        common_expand_by_mask = common_expand_by_mask,
        ideal_weights = ideal_weights,
        weights_chol = weights_chol,
        n = n,
        df = df)
}

cadzow_with_mgn <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow = NULL,
                            epsilon = 1e-6, use_mgn = TRUE,
                            it_limit = 100, debug = FALSE,
                            envelope = unit_envelope(this, series), set_seed = NULL, ...) {
    series_for_mgn <- series

    if (!is.null(set_seed)) {
        set_seed()
    }

    list2env(prepare_cadzow_with_mgn(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope), environment())


    cadzow_data <- oblique_cadzow_eps(this, series = series_for_cadzow, r = r,
                               left_chol_mat = left_chol_mat,
                               right_chol_mat = right_chol_mat, weights_mat = weights,
                               epsilon = epsilon, it_limit = it_limit,
                               debug = debug, ...)

    answer <- cadzow_data

    answer$signal <- common_expand_by_mask(answer$signal)
    

    if (use_mgn) {
        answer <- mgn(this, series_for_mgn, answer$signal, r, ideal_weights,
                                   debug = debug, weights_chol = weights_chol, ...)
        answer$cadzow_it <- cadzow_data$it
    }

    if (!is.null(answer$noise_norm)) {
        answer$loglikelihood <- -n/2 * log(answer$noise_norm / n)
        answer$df <- df
    }

    answer
}

prepare_hlra_ssa <- function(obj, ...) UseMethod("prepare_hlra_ssa")

prepare_hlra_ssa.1d <- function(obj, series, L) {
    K <- length(series) - L + 1
    list(pseudord = rep(1, K), whitecoefs = numeric(0), 
        pseudoenvelope = rep(1, length(series)))
}

prepare_hlra_ssa.1dm <- function(obj, series, L) {
    Ks <- sapply(series, function(series) length(as.numeric(series))) - L + 1
    list(pseudord = sapply(Ks, function(i) rep(1, i), simplify = FALSE), 
        whitecoefs = sapply(seq_along(series), function(i) numeric(0), simplify = FALSE),
        pseudoenvelope = sapply(series, function(i) rep(1, i), simplify = FALSE))
}

hlra_ssa <- function(obj, series, L, r, debug, set_seed, ...) {
    prepare_list <- prepare_hlra_ssa(obj, series, L)

    cadzow_with_mgn(obj, series, L, r, coefs = prepare_list$whitecoefs,
        right_diag = prepare_list$pseudord, epsilon = 1, use_mgn = FALSE, debug = debug,
        envelope = prepare_list$pseudoenvelope, series_for_cadzow = series,
        set_seed = set_seed, ...)
}