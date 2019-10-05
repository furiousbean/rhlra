eval_basis <- function(obj, ...) UseMethod("eval_basis")
eval_pseudograd <- function(obj, ...) UseMethod("eval_pseudograd")
eval_sylvester_grad <- function(obj, ...) UseMethod("eval_sylvester_grad")
eval_tangent_basis <- function(obj, ...) UseMethod("eval_tangent_basis")

eval_basis.default <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1

    result <- .Call("eval_basisC", N, glrr, PACKAGE = "rhlra")

    list(basis = matrix(result[1 : (N*r)], nrow = N),
         basis_fourier = matrix(result[(N*r + 1) : (2*N*r)], nrow = N),
         unitroots = result[(2*N*r + 1) : (N*(2 * r + 1))],
         A_f = result[(N*(2 * r + 1) + 1) : (N*(2 * r + 2))],
         alpha = Re(result[N*(2 * r + 2) + 1]),
         glrr = glrr)
}

eval_basis.compensated <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)
    r <- length(glrr) - 1

    result <- .Call("eval_basis_compensatedC", N, glrr, PACKAGE = "rhlra")
    list(basis = matrix(result[1 : (N*r)], nrow = N),
         basis_fourier = matrix(result[(N*r + 1) : (2*N*r)], nrow = N),
         unitroots = result[(2*N*r + 1) : (N*(2 * r + 1))],
         A_f = result[(N*(2 * r + 1) + 1) : (N*(2 * r + 2))],
         alpha = Re(result[N*(2 * r + 2) + 1]),
         qrinvmat = matrix(result[(N*(2 * r + 2) + 2) : (N*(2 * r + 2) + r * r + 1)], nrow = r),
         glrr = glrr)
}



eval_pseudograd.default <- function(this, N, glrr, signal, tau, basis_obj) {
    N <- as.integer(N)
    tau <- as.integer(tau)
    glrr <- as.numeric(glrr)
    dim(glrr) <- length(glrr)

    signal <- as.numeric(signal)
    dim(signal) <- length(signal)

    r <- length(glrr) - 1

    result <- .Call("eval_pseudogradC", N, glrr, basis_obj$A_f, basis_obj$alpha, tau, signal, PACKAGE = "rhlra")
    matrix(result, nrow = N)
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

    result <- .Call("eval_sylvester_gradC", N, glrr, basis_obj$A_f, basis_obj$alpha, tau, signal, PACKAGE = "rhlra")
    as.complex(result)
}

eval_sylvester_grad.compensated <- eval_sylvester_grad.default

eval_tangent_basis.default <- function(this, N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    result <- .Call("eval_tangent_basisC", N, glrr, PACKAGE = "rhlra")

    matrix(result, nrow = N)
}

eval_tangent_basis.compensated <- function(N, glrr) {
    N <- as.integer(N)
    glrr <- as.vector(glrr)
    dim(glrr) <- length(glrr)

    result <- .Call("eval_tangent_basis_compensatedC", N, glrr, PACKAGE = "rhlra")

    matrix(result, nrow = N)
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
