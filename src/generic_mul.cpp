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
        reinterpret_cast<std::complex<double>*>(COMPLEX(seriesfftexp));

    SEXP xexp = PROTECT(allocVector(REALSXP, L));

    double* x = REAL(xexp);

    for (int i = 0; i < L; i++) {
        x[i] = 0;
    }

    fftw_complex* input_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* output_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));

    std::complex<double>* input =
        reinterpret_cast<std::complex<double>*>(input_fftw);
    std::complex<double>* output =
        reinterpret_cast<std::complex<double>*>(output_fftw);

    for (int i = 0; i < N; i++) {
        input[i] = 0;
        output[i] = 0;
    }

    //right_chol_mat %*% v
    for (int i = 0; i < K; i++)
        for (int j = 0; j < std::min(i + 1, Kp); j++) {
            input[K - i - 1] += rightchol[j * K + i] * v[i - j];
        }

    fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw, output_fftw,
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        input[i] = output[i] * seriesfft[i];
    }

    for (int i = 0; i < N; i++) {
        output[i] = 0;
    }

    my_plan = fftw_plan_dft_1d(N, input_fftw, output_fftw,
        FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < L; i++)
        for (int j = 0; j < std::min(i + 1, Lp); j++) {
            x[i - j] += leftchol[j * L + i] * output[i + K - 1].real() / N;
        }

    fftw_free(input_fftw);
    fftw_free(output_fftw);

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

    fftw_complex* input_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* output_L_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* output_K_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));

    std::complex<double>* input =
        reinterpret_cast<std::complex<double>*>(input_fftw);
    std::complex<double>* output_L =
        reinterpret_cast<std::complex<double>*>(output_L_fftw);
    std::complex<double>* output_K =
        reinterpret_cast<std::complex<double>*>(output_K_fftw);

    for (int i = 0; i < N; i++) {
        input[i] = 0;
        output_L[i] = 0;
        output_K[i] = 0;
    }

    //left_chol_mat %*% u
    for (int i = 0; i < L; i++)
        for (int j = 0; j < std::min(i + 1, Lp); j++) {
            input[i] += leftchol[j * L + i] * u[i - j];
        }

    fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw, output_L_fftw,
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        input[i] = 0;
    }

    //right_chol_mat %*% u
    for (int i = 0; i < K; i++)
        for (int j = 0; j < std::min(i + 1, Kp); j++) {
            input[i] += rightchol[j * K + i] * v[i - j];
        }

    my_plan = fftw_plan_dft_1d(N, input_fftw, output_K_fftw,
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++) {
        input[i] = output_L[i] * output_K[i];
    }


    my_plan = fftw_plan_dft_1d(N, input_fftw, output_K_fftw,
        FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    for (int i = 0; i < N; i++)
        x[i] = output_K[i].real() / N;

    fftw_free(input_fftw);
    fftw_free(output_L_fftw);
    fftw_free(output_K_fftw);

    UNPROTECT(1);

    return xexp;
}

SEXP hlra_fft_common(SEXP seriesexp, bool direction) {

    int N = length(seriesexp);

    std::complex<double>* series =
        reinterpret_cast<std::complex<double>*>(COMPLEX(seriesexp));

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N));

    std::complex<double>* x =
        reinterpret_cast<std::complex<double>*>(COMPLEX(xexp));

    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }

    fftw_complex* input_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex* output_fftw = reinterpret_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * N));

    std::complex<double>* input =
        reinterpret_cast<std::complex<double>*>(input_fftw);
    std::complex<double>* output =
        reinterpret_cast<std::complex<double>*>(output_fftw);

    for (int i = 0; i < N; i++) {
        input[i] = series[i];
        output[i] = 0;
    }

    fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw, output_fftw,
        direction ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(my_plan);
    fftw_destroy_plan(my_plan);

    if (direction) {
        for (int i = 0; i < N; i++) {
            x[i] = output[i];
        }
    } else {
        for (int i = 0; i < N; i++) {
            x[i] = output[i] / (double)N;
        }
    }

    fftw_free(input_fftw);
    fftw_free(output_fftw);

    UNPROTECT(1);

    return xexp;
}

SEXP hlra_fft(SEXP seriesexp) {
    return hlra_fft_common(seriesexp, true);
};

SEXP hlra_ifft(SEXP seriesexp) {
    return hlra_fft_common(seriesexp, false);
};

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
