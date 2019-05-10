#include <stdlib.h>
#include <complex>
#include <cmath>
#include <R.h>
#include <Rinternals.h>
#include <fftw3.h>

SEXP generic_mul(SEXP vexp,
    SEXP leftcholexp,
    SEXP rightcholexp,
    SEXP seriesfftexp) {

    int N = length(seriesfftexp);
    int L = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[0];
    int Lp = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[1];
    int K = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[0];
    int Kp = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[1];

    // Rprintf("%d %d %d %d %d\n", N, L, Lp, K, Kp);

    double* v = REAL(vexp);
    double* leftchol = REAL(leftcholexp);
    double* rightchol = REAL(rightcholexp);

    std::complex<double>* seriesfft =
        (std::complex<double>*)COMPLEX(seriesfftexp);

    SEXP xexp = PROTECT(allocVector(REALSXP, L));

    double* x = REAL(xexp);

    for (int i = 0; i < L; i++) {
        x[i] = 0;
    }

    fftw_complex* f_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* f_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    for (int i = 0; i < N; i++) {
        f_input[i][0] = 0;
        f_input[i][1] = 0;
        f_output[i][0] = 0;
        f_output[i][1] = 0;
    }

    //right_chol_mat %*% v
    for (int i = 0; i < K; i++)
        for (int j = 0; j < std::min(i + 1, Kp); j++) {
            //right_chol_mat %*% v
            //f_input[i][0] += rightchol[j * K + i] * v[i - j];
            //rev(right_chol_mat %*% v)
            f_input[K - i - 1][0] += rightchol[j * K + i] * v[i - j];
        }

    fftw_plan my_plan = fftw_plan_dft_1d(N, f_input, f_output, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        *((std::complex<double>*)f_input[i]) =
        *((std::complex<double>*)f_output[i]) * seriesfft[i];
    }

    for (int i = 0; i < N; i++) {
        f_output[i][0] = 0;
        f_output[i][1] = 0;
    }

    my_plan = fftw_plan_dft_1d(N, f_input, f_output, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < L; i++)
        for (int j = 0; j < std::min(i + 1, Lp); j++) {
            x[i - j] += leftchol[j * L + i] * f_output[i + K - 1][0] / N;
        }

    fftw_free(f_input);
    fftw_free(f_output);

    UNPROTECT(1);

    return xexp;
}

SEXP generic_diag_one_triple(SEXP uexp,
    SEXP vexp,
    SEXP leftcholexp,
    SEXP rightcholexp) {

    int L = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[0];
    int Lp = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[1];
    int K = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[0];
    int Kp = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[1];

    int N = L + K - 1;

    // Rprintf("%d %d %d %d %d\n", N, L, Lp, K, Kp);

    double* u = REAL(uexp);
    double* v = REAL(vexp);
    double* leftchol = REAL(leftcholexp);
    double* rightchol = REAL(rightcholexp);

    SEXP xexp = PROTECT(allocVector(REALSXP, N));

    double* x = REAL(xexp);

    fftw_complex* f_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* f_output_L = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* f_output_K = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    for (int i = 0; i < N; i++) {
        f_input[i][0] = 0;
        f_input[i][1] = 0;
        f_output_L[i][0] = 0;
        f_output_L[i][1] = 0;
        f_output_K[i][0] = 0;
        f_output_K[i][1] = 0;
    }

    //left_chol_mat %*% u
    for (int i = 0; i < L; i++)
        for (int j = 0; j < std::min(i + 1, Lp); j++) {
            f_input[i][0] += leftchol[j * L + i] * u[i - j];
        }

    fftw_plan my_plan = fftw_plan_dft_1d(N, f_input, f_output_L, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        f_input[i][0] = 0;
        f_input[i][1] = 0;
    }

    //right_chol_mat %*% u
    for (int i = 0; i < K; i++)
        for (int j = 0; j < std::min(i + 1, Kp); j++) {
            f_input[i][0] += rightchol[j * K + i] * v[i - j];
        }

    my_plan = fftw_plan_dft_1d(N, f_input, f_output_K, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        *((std::complex<double>*)f_input[i]) =
        *((std::complex<double>*)f_output_L[i]) * *((std::complex<double>*)f_output_K[i]);
    }


    my_plan = fftw_plan_dft_1d(N, f_input, f_output_K, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++)
        x[i] = f_output_K[i][0] / N;

    fftw_free(f_input);
    fftw_free(f_output_L);
    fftw_free(f_output_K);

    UNPROTECT(1);

    return xexp;
}

SEXP hlra_fft(SEXP seriesexp) {

    int N = length(seriesexp);

    std::complex<double>* series =
        (std::complex<double>*)COMPLEX(seriesexp);

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N));

    std::complex<double>* x = (std::complex<double>*)COMPLEX(xexp);

    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }

    fftw_complex* f_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* f_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    for (int i = 0; i < N; i++) {
        f_input[i][0] = series[i].real();
        f_input[i][1] = series[i].imag();
        f_output[i][0] = 0;
        f_output[i][1] = 0;
    }

    fftw_plan my_plan = fftw_plan_dft_1d(N, f_input, f_output, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        x[i] = *((std::complex<double>*)f_output[i]);
    }

    fftw_free(f_input);
    fftw_free(f_output);

    UNPROTECT(1);

    return xexp;
}

SEXP hlra_ifft(SEXP seriesexp) {

    int N = length(seriesexp);

    std::complex<double>* series =
        (std::complex<double>*)COMPLEX(seriesexp);

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N));

    std::complex<double>* x = (std::complex<double>*)COMPLEX(xexp);

    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }

    fftw_complex* f_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    fftw_complex* f_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    for (int i = 0; i < N; i++) {
        f_input[i][0] = series[i].real();
        f_input[i][1] = series[i].imag();
        f_output[i][0] = 0;
        f_output[i][1] = 0;
    }

    fftw_plan my_plan = fftw_plan_dft_1d(N, f_input, f_output, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        x[i] = *((std::complex<double>*)f_output[i]) / (double)N;
    }

    fftw_free(f_input);
    fftw_free(f_output);

    UNPROTECT(1);

    return xexp;
}

extern "C" {
    SEXP generic_mulC(SEXP vexp, SEXP leftcholexp,
    SEXP rightcholexp, SEXP seriesfftexp) {
        return generic_mul(vexp, leftcholexp, rightcholexp, seriesfftexp);
    }
}

extern "C" {
    SEXP generic_diag_one_tripleC(SEXP uexp, SEXP vexp,
    SEXP leftcholexp, SEXP rightcholexp) {
        return generic_diag_one_triple(uexp, vexp, leftcholexp, rightcholexp);
    }
}

extern "C" {
    SEXP hlra_fftC(SEXP seriesexp) {
        return hlra_fft(seriesexp);
    }
}

extern "C" {
    SEXP hlra_ifftC(SEXP seriesexp) {
        return hlra_ifft(seriesexp);
    }
}
