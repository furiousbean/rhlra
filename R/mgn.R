get_comp_space_by_glrr <- function(this, N, glrr, glrr_2 = FALSE, return_pack = FALSE) {
    if (!glrr_2) {
        result <- eval_basis(this, N, glrr)

        if (return_pack) {
            return(result)
        } else {
            return(result$basis)
        }
    }
    if (glrr_2) {
        result <- eval_tangent_basis(this, N, glrr)

        if (return_pack) {
            stop("Unsupported")
        }

        return(result)
    }
}

get_pseudograd <- function(this, N, zspace_pack, signal, j) {
    eval_pseudograd(this, N, zspace_pack$glrr, signal, j, zspace_pack)
}

get_space_by_glrr <- function(this, N, glrr, glrr_2 = FALSE) {
    pre_answer <- get_comp_space_by_glrr(this, N, glrr, glrr_2 = glrr_2)
    if (glrr_2) {
        glrr <- ssa_convolve(glrr, glrr)
    }
    r <- length(glrr) - 1
    j <- which.max(abs(glrr))
    begin_indices <- numeric(0)
    if (j > 1) begin_indices <- 1:(j - 1)
    end_indices <- numeric(0)
    if (j < r + 1) end_indices <- (N - r + j):N
    a_coefs <- solve(pre_answer[c(begin_indices, end_indices), ])
    Re(pre_answer %*% a_coefs)
}

get_glrr_from_nonlrf_series <- function(obj, ...) UseMethod("get_glrr_from_nonlrf_series")

get_glrr_from_nonlrf_series.1d <- function(this, series, r) {
    L <- r + 1
    N <- length(series)
        trmat <- traj_matrix(series[!is.na(series)], L)
    svd(trmat)$u[, L]
}

get_glrr_from_nonlrf_series.1dm <- function(this, series, r) {
    L <- r + 1
    Ns <- sapply(series, length)
    trmat <- do.call(cbind,
         sapply_ns(series, function(x) traj_matrix(x[!is.na(x)], L)))
    svd(trmat)$u[, L]
}

weighted_project_rotate_basis <- function(zspace, chol_weights) {
    real_part <- chol_weights %*% Re(zspace)
    imag_part <- chol_weights %*% Im(zspace)
    zspace_rot <- complex(real = as.numeric(real_part),
                          imaginary = as.numeric(imag_part))
    dim(zspace_rot) <- dim(zspace)
    qr(zspace_rot)
}

weighted_project_onto_zspace_qr <- function(series, qrobj, zspace, chol_weights) {
    series_rot <- as.numeric(chol_weights %*% series)

    coef <- qr.coef(qrobj, series_rot)

    Re(zspace %*% coef)
}

weighted_project_onto_zspace_coef_rot <- function(obj, ...) UseMethod("weighted_project_onto_zspace_coef_rot")

weighted_project_onto_zspace_coef_rot.1d <- function(this, series, zspace, chol_weights) {
    list(zspace_rot = as.matrix(chol_weights %*% zspace),
        series_rot = as.numeric(chol_weights %*% series))
}

weighted_project_onto_zspace_coef_rot.1dm <- function(this, series, zspace, chol_weights) {
    list(zspace_rot = do.call(rbind, sapply_ns(seq_along(series),
                    function(i) as.matrix(chol_weights[[i]] %*% zspace[[i]]))),
        series_rot = glue_series_lists(sapply_ns(seq_along(series),
                    function(i) as.numeric(chol_weights[[i]] %*% series[[i]]))))
}

weighted_project_onto_zspace_coef <- function(this, series, zspace, chol_weights, tol = 1e-14) {
    list2env(weighted_project_onto_zspace_coef_rot(this, series, zspace, chol_weights), environment())
    svdobj <- svd(zspace_rot)

    svdobj$d[svdobj$d < svdobj$d[1] * tol] <- Inf
    svdobj$v %*% ((t(svdobj$u) %*% series_rot) / svdobj$d)
}

prepare_find_step <- function(obj, ...) UseMethod("prepare_find_step")

prepare_find_step.1d <- function(this, signal, series, r, j, zspace_pack, weights_chol,
                                 sylvester_problem) {
    noise <- series - signal
    N <- length(series)

    K <- N - r

    pseudograd <- NULL

    if (!sylvester_problem) {
        pseudograd <- Re(get_pseudograd(this, N, zspace_pack, signal, j))
    } else {
        quotient <- Re(eval_sylvester_grad(this, N, zspace_pack$glrr,
                                           as.numeric(noise), j, zspace_pack)[1:K])

        indices <- 1:(r+1)
        indices <- indices[-j]

        pseudograd <- sapply(indices, function(i) {
            grad <- complex(N)
            grad[i:(i + K - 1)] <- quotient
            grad
        })
    }

    pseudograd_minus <- sapply(1:r, function(i)
        weighted_project_onto_zspace_qr(Re(pseudograd[, i]), zspace_pack$qrobj,
                                        zspace_pack$basis, weights_chol))

    if (!sylvester_problem) {
        pseudograd <- pseudograd - pseudograd_minus
    } else {
        pseudograd <- pseudograd_minus
    }

    list(pseudograd = pseudograd, noise = noise)
}

prepare_find_step.1dm <- function(this, signal, series, r, j, zspace_pack, weights_chol,
                                  sylvester_problem) {
    noise <- list()
    metathis <- this
    class(metathis) <- class(metathis)[class(metathis) != '1dm']
    class(metathis) <- append(class(metathis), "1d")

    pseudograd <- sapply_ns(seq_along(series), function(k) {
        pfs_1d <- prepare_find_step(metathis, signal[[k]], series[[k]], r, j,
                                    zspace_pack[[k]], weights_chol[[k]],
                                    sylvester_problem)
        noise[[k]] <<- pfs_1d$noise
        pfs_1d$pseudograd
    })

    list(pseudograd = pseudograd, noise = noise)
}

find_step <- function(this, signal, series, glrr, zspace_pack, weights_chol, debug = FALSE,
                      sylvester_problem, ...) {
    r <- length(glrr) - 1
    j <- which.max(abs(glrr))

    prepare_obj <- prepare_find_step(this, signal, series, r, j, zspace_pack, weights_chol,
                                     sylvester_problem)

    obj_for_coefs <- prepare_obj$noise

    if (sylvester_problem) {
        obj_for_coefs <- signal
    }

    used_coefs <- weighted_project_onto_zspace_coef(this, obj_for_coefs,
        prepare_obj$pseudograd, weights_chol)

    ans <- numeric(r+1)

    ans[-j] <- used_coefs

    ans
}

generic_get_zspace_pack <- function(this, N, glrr, weights_chol) {
    result <- get_comp_space_by_glrr(this, N, glrr, return_pack = TRUE)
    result$qrobj <- weighted_project_rotate_basis(result$basis, weights_chol)
    result
}

generic_get_signal <- function(this, x, zspace_pack, weights_chol)
    weighted_project_onto_zspace_qr(x, zspace_pack$qrobj, zspace_pack$basis, weights_chol)

prepare_mgn <- function(obj, ...) UseMethod("prepare_mgn")

prepare_mgn.1d <- function(this, series, r, weights, weights_chol) {
    series <- as.numeric(series)

    if (any(is.na(series))) {
        series <- null_series(series)
    }

    if (is.null(weights_chol)) {
        weights_chol <- chol(weights)
    }
    N <- length(series)
    get_zspace_pack <- generic_get_zspace_pack

    get_signal <- generic_get_signal

    inner_product <- function(x, y) generic_inner_product_chol(x, y, weights_chol)
    minus <- function(x, y) x - y
    plus <- function(x, y) x + y
    mult <- function(x, y) x * y

    list(series = series, get_zspace_pack = get_zspace_pack, get_signal = get_signal,
        inner_product = inner_product, minus = minus, plus = plus, mult = mult, N = N)
}

prepare_mgn.1dm <- function(this, series, r, weights, weights_chol) {
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
    get_zspace_pack <- function(this, Ns, v, weights_chol)
        sapply_ns(seq_along(Ns),
                  function(i) generic_get_zspace_pack(this, Ns[i], v, weights_chol[[i]]))
    get_signal <- function(this, x, zspace_pack, weights_chol)
        sapply_ns(seq_along(x),
                  function(i) generic_get_signal(this, x[[i]], zspace_pack[[i]], weights_chol[[i]]))

    inner_product <- function(x, y) sum(mapply(generic_inner_product_chol, x, y, weights_chol))
    minus <- function(x, y) mapply("-", x, y, SIMPLIFY = FALSE)
    plus <- function(x, y) mapply("+", x, y, SIMPLIFY = FALSE)
    mult <- function(x, y) mapply("*", x, y, SIMPLIFY = FALSE)

    list(series = series, get_zspace_pack = get_zspace_pack, get_signal = get_signal,
    inner_product = inner_product, minus = minus, plus = plus, mult = mult, N = N)
}

mgn <- function(this, series, signal, r, weights,
                               mgn_search_threshold = 1e-7,
                               mgn_it_limit = 100,
                               mgn_j_limit = 4, debug = FALSE,
                               mgn_max_step = 1,
                               glrr_initial = NULL,
                               weights_chol = NULL,
                               sylvester_problem = FALSE,
                               ...) {
    #MGN
    cur_glrr <- NA
    best_signal <- NA
    zspace_pack <- NA

    source_series <- series

    list2env(prepare_mgn(this, series, r, weights, weights_chol), environment())

    initial_best_signal <- NA

    if (is.null(glrr_initial)) {
        cur_glrr <- get_glrr_from_nonlrf_series(this, signal, r)
    } else {
        cur_glrr <- glrr_initial
    }

    polarity_mult <- 1

    if (sylvester_problem) {
        polarity_mult <- -1
    }

    zspace_pack <- get_zspace_pack(this, N, cur_glrr, weights_chol)
    best_signal <- get_signal(this, series, zspace_pack, weights_chol)

    initial_best_signal <- best_signal

    if (any(is.na(cur_glrr))) return(signal)
    if (any(is.na(best_signal))) return(signal)

    no_localsearch_flag <- FALSE

    stop_criterion <- function(signal, step, prev_step) {
        #step-to-solution ratio
        coef <- sqrt(inner_product(minus(best_signal, signal), minus(best_signal, signal))/
            inner_product(best_signal, best_signal))

        if (coef > mgn_search_threshold) {
            #closer to the solution than previous best_signal
            return(polarity_mult * inner_product(minus(best_signal, signal),
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
    prev_step <- mgn_max_step

    while(TRUE) {
        step <- NA
        tryCatch({step <- find_step(this, signal = best_signal, series = series, glrr = cur_glrr,
                                        zspace_pack = zspace_pack, weights_chol = weights_chol,
                                        debug = debug, sylvester_problem = sylvester_problem, ...)},
                 warning = function(x) {print(x)})

        if (any(is.na(step))) break

        alpha <- mgn_max_step
        j <- 0
        next_it <- FALSE
        while (j < mgn_j_limit) {
            new_glrr <- cur_glrr + step * alpha

            new_zspace_pack <- get_zspace_pack(this, N, new_glrr, weights_chol)
            new_signal <- get_signal(this, series, new_zspace_pack, weights_chol)

            if (stop_criterion(new_signal, step * alpha, prev_step)) {
                next_it <- TRUE
                it <- it + 1
                dists <- c(dists, noise_distance(new_signal))
                best_signal <- new_signal

                prev_step <- step * alpha
                cur_glrr <- new_glrr
                zspace_pack <- new_zspace_pack
                break
            }
            if (no_localsearch_flag) {
                j <- mgn_j_limit
                break
            }
            alpha <- alpha / 2
            j <- j + 1
        }

        if (j == mgn_j_limit && (debug == "P" || debug)) {
            cat("no steps found\n")
        } else {
            if (j > 0 && (debug == "P" || debug)) cat(sprintf("reduction, it = %d, j = %d\n", it, j))
        }

        if (!next_it) break
        if (it > mgn_it_limit) break
    }

    if (debug == "P" || debug) {
        cat(sprintf("%d MGN iterations done\n", it))
    }

    if (debug == "P") {
        draw_full_plot <- function(N, series, best_signal, signal, initial_best_signal) {
            matplot(1:N, cbind(series, best_signal, signal),
                    type = "l", lty = 1:3,
                    col = c("black", "blue", "red"), main = "Approximations")
            legend("topleft", legend = c("Source series", "MGN", "Initial"),
                   col = c("black", "blue", "red"), lty = 1:3)
        }

        if (is.null(signal)) {
            signal <- initial_best_signal
        }

        if (!is.list(series)) {
            draw_full_plot(N, source_series, best_signal, signal, initial_best_signal)
        } else {
            mapply(draw_full_plot, N, source_series, best_signal, signal, initial_best_signal)
        }

        plot(polarity_mult * dists/sqrt(sum(N)), type = "b", main = "Normalized objective function", xlab = "iteration",
            ylab = "Value")

    }

    # stop()
    list(signal = best_signal, glrr = cur_glrr, it = it, dists = polarity_mult * dists,
        noise = minus(series, best_signal),
        noise_norm = inner_product(minus(series, best_signal), minus(series, best_signal)))
}
