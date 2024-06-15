eval_basis <- function(obj, ...) UseMethod("eval_basis")
eval_pseudograd <- function(obj, ...) UseMethod("eval_pseudograd")
eval_sylvester_grad <- function(obj, ...) UseMethod("eval_sylvester_grad")
eval_tangent_basis <- function(obj, ...) UseMethod("eval_tangent_basis")

eval_basis.default <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1
    if (!all(is.finite(glrr)) || all(glrr == 0) || !(r > 0) || !(N > r)) {
        stop("Critical: eval_basis checks failed")
    }

    .Call("eval_basisC", N, glrr, PACKAGE = "rhlra")
}

eval_basis.compensated <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1
    if (!all(is.finite(glrr)) || all(glrr == 0) || !(r > 0) || !(N > r)) {
        stop("Critical: eval_basis checks failed")
    }

    .Call("eval_basis_compensatedC", N, glrr, PACKAGE = "rhlra")
}



eval_pseudograd.default <- function(this, N, glrr, signal, tau, basis_obj) {
    N <- as.integer(N)
    tau <- as.integer(tau)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)

    signal <- as.numeric(signal)
    dim(signal) <- length(signal)

    r <- length(glrr) - 1

    if (!all(is.finite(glrr)) || !all(is.finite(signal)) || all(glrr == 0) ||
        !(r > 0) || !(N > r) || !(length(signal) == N)) {
        stop("Critical: eval_pseudograd checks failed")
    }

    .Call("eval_pseudogradC", N, glrr, basis_obj$A_f, basis_obj$alpha, tau, signal, PACKAGE = "rhlra")
}

eval_pseudograd.compensated <- eval_pseudograd.default

eval_sylvester_grad.default <- function(this, N, glrr, signal, tau, basis_obj) {
    N <- as.integer(N)
    tau <- as.integer(tau)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)

    signal <- as.numeric(signal)
    dim(signal) <- length(signal)

    r <- length(glrr) - 1

    if (!all(is.finite(glrr)) || !all(is.finite(signal)) || all(glrr == 0) ||
        !(r > 0) || !(N > r) || !(length(signal) == N)) {
        stop("Critical: eval_pseudograd checks failed")
    }

    .Call("eval_sylvester_gradC", N, glrr, basis_obj$A_f, basis_obj$alpha, tau, signal, PACKAGE = "rhlra")
}

eval_sylvester_grad.compensated <- eval_sylvester_grad.default

eval_tangent_basis.default <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    r <- length(glrr) - 1

    if (!all(is.finite(glrr)) || all(glrr == 0) || !(r > 0) || !(N > 2 * r)) {
        stop("Critical: eval_tangent_basis checks failed")
    }

    .Call("eval_tangent_basisC", N, glrr, PACKAGE = "rhlra")
}

eval_tangent_basis.compensated <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    r <- length(glrr) - 1

    if (!all(is.finite(glrr)) || all(glrr == 0) || !(r > 0) || !(N > 2 * r)) {
        stop("Critical: eval_tangent_basis checks failed")
    }

    .Call("eval_tangent_basis_compensatedC", N, glrr, PACKAGE = "rhlra")
}

glue_series_lists <- function(series_list) {
    .Call("glue_series_lists", series_list, PACKAGE = "rhlra")
}

generic_mul <- function(v, left_chol_mat,
                        right_chol_mat, series_fft) {
    .Call("generic_mulC", v, left_chol_mat, right_chol_mat, series_fft, PACKAGE = "rhlra")
}

generic_tmul <- function(v, left_chol_mat, right_chol_mat, series_fft) {
    generic_mul(v, right_chol_mat, left_chol_mat, series_fft)
}

generic_diag_one_triple <- function(u, v, left_chol_mat,
                        right_chol_mat) {
    .Call("generic_diag_one_tripleC", u, v, left_chol_mat, right_chol_mat, PACKAGE = "rhlra")
}

hlra_fft <- function(series) {
    .Call("hlra_fftC", as.complex(series), PACKAGE = "rhlra")
}

hlra_ifft <- function(series) {
    .Call("hlra_ifftC", as.complex(series), PACKAGE = "rhlra")
}
