#include <complex>
#include <cmath>
#include <exception>
#include <memory>
#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>
#include <fftw3.h>

#include "alloc.h"


namespace {
    using DoubleComplex = std::complex<double>;
}

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

    try {
        SexpWrapper xexp(allocVector(REALSXP, L));

        double* x = REAL(xexp);

        for (int i = 0; i < L; i++) {
            x[i] = 0;
        }

        std::shared_ptr<fftw_complex> input_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());
        std::shared_ptr<fftw_complex> output_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());

        DoubleComplex* input = reinterpret_cast<DoubleComplex*>(input_fftw.get());
        DoubleComplex* output = reinterpret_cast<DoubleComplex*>(output_fftw.get());

        for (int i = 0; i < N; i++) {
            input[i] = 0;
            output[i] = 0;
        }

        //right_chol_mat %*% v
        for (int i = 0; i < K; i++)
            for (int j = 0; j < std::min(i + 1, Kp); j++) {
                input[K - i - 1] += rightchol[j * K + i] * v[i - j];
            }

        fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_fftw.get(),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (int i = 0; i < N; i++) {
            input[i] = output[i] * seriesfft[i];
        }

        for (int i = 0; i < N; i++) {
            output[i] = 0;
        }

        my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_fftw.get(),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (int i = 0; i < L; i++)
            for (int j = 0; j < std::min(i + 1, Lp); j++) {
                x[i - j] += leftchol[j * L + i] * output[i + K - 1].real() / N;
            }

        return xexp;
    } catch (std::exception& e) {
        error(e.what());
        return vexp;
    }
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

    try {
        SexpWrapper xexp(allocVector(REALSXP, N));

        double* x = REAL(xexp);

        std::shared_ptr<fftw_complex> input_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());
        std::shared_ptr<fftw_complex> output_L_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());
        std::shared_ptr<fftw_complex> output_K_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());

        DoubleComplex* input = reinterpret_cast<DoubleComplex*>(input_fftw.get());
        DoubleComplex* output_L = reinterpret_cast<DoubleComplex*>(output_L_fftw.get());
        DoubleComplex* output_K = reinterpret_cast<DoubleComplex*>(output_K_fftw.get());

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

        fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_L_fftw.get(),
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

        my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_K_fftw.get(),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (int i = 0; i < N; i++) {
            input[i] = output_L[i] * output_K[i];
        }


        my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_K_fftw.get(),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (int i = 0; i < N; i++)
            x[i] = output_K[i].real() / N;

        return xexp;
    } catch (std::exception& e) {
        error(e.what());
        return uexp;
    }
}

SEXP hlra_fft_common(SEXP seriesexp, bool direction) {

    int N = length(seriesexp);

    std::complex<double>* series =
        reinterpret_cast<std::complex<double>*>(COMPLEX(seriesexp));

    try {
        SexpWrapper xexp(allocVector(CPLXSXP, N));

        std::complex<double>* x =
            reinterpret_cast<std::complex<double>*>(COMPLEX(xexp));

        for (int i = 0; i < N; i++) {
            x[i] = 0;
        }

        std::shared_ptr<fftw_complex> input_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());
        std::shared_ptr<fftw_complex> output_fftw(FftwArrayAllocator<fftw_complex>(N),
            FftwArrayDeleter<fftw_complex>());

        DoubleComplex* input = reinterpret_cast<DoubleComplex*>(input_fftw.get());
        DoubleComplex* output = reinterpret_cast<DoubleComplex*>(output_fftw.get());

        for (int i = 0; i < N; i++) {
            input[i] = series[i];
            output[i] = 0;
        }

        fftw_plan my_plan = fftw_plan_dft_1d(N, input_fftw.get(), output_fftw.get(),
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

        return xexp;
    } catch (std::exception& e) {
        error(e.what());
        return seriesexp;
    }
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
