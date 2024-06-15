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

    std::size_t N = length(seriesfftexp);
    std::size_t L = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[0];
    std::size_t Lp = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[1];
    std::size_t K = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[0];
    std::size_t Kp = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[1];

    // Rprintf("%d %d %d %d %d\n", N, L, Lp, K, Kp);

    double* v = REAL(vexp);
    double* leftchol = REAL(leftcholexp);
    double* rightchol = REAL(rightcholexp);

    std::complex<double>* seriesfft =
        reinterpret_cast<std::complex<double>*>(COMPLEX(seriesfftexp));

    try {
        SexpWrapper xexp(allocVector(REALSXP, L));

        double* x = REAL(xexp);

        for (std::size_t i = 0; i < L; i++) {
            x[i] = 0;
        }

        std::shared_ptr<DoubleComplex> input_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());
        std::shared_ptr<DoubleComplex> output_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());

        DoubleComplex* input = input_fftw.get();
        DoubleComplex* output = output_fftw.get();

        for (std::size_t i = 0; i < N; i++) {
            input[i] = 0;
            output[i] = 0;
        }

        //right_chol_mat %*% v
        for (std::size_t i = 0; i < K; i++)
            for (std::size_t j = 0; j < std::min(i + 1, Kp); j++) {
                input[K - i - 1] += rightchol[j * K + i] * v[i - j];
            }

        fftw_plan my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (std::size_t i = 0; i < N; i++) {
            input[i] = output[i] * seriesfft[i];
        }

        for (std::size_t i = 0; i < N; i++) {
            output[i] = 0;
        }

        my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (std::size_t i = 0; i < L; i++)
            for (std::size_t j = 0; j < std::min(i + 1, Lp); j++) {
                x[i - j] += leftchol[j * L + i] * output[i + K - 1].real() / N;
            }

        return xexp;
    } catch (std::exception& e) {
        error("%s", e.what());
        return vexp;
    }
}

SEXP generic_diag_one_triple(SEXP uexp,
    SEXP vexp,
    SEXP leftcholexp,
    SEXP rightcholexp) {

    std::size_t L = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[0];
    std::size_t Lp = INTEGER(getAttrib(leftcholexp, R_DimSymbol))[1];
    std::size_t K = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[0];
    std::size_t Kp = INTEGER(getAttrib(rightcholexp, R_DimSymbol))[1];

    std::size_t N = L + K - 1;

    // Rprintf("%d %d %d %d %d\n", N, L, Lp, K, Kp);

    double* u = REAL(uexp);
    double* v = REAL(vexp);
    double* leftchol = REAL(leftcholexp);
    double* rightchol = REAL(rightcholexp);

    try {
        SexpWrapper xexp(allocVector(REALSXP, N));

        double* x = REAL(xexp);

        std::shared_ptr<DoubleComplex> input_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());
        std::shared_ptr<DoubleComplex> output_L_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());
        std::shared_ptr<DoubleComplex> output_K_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());

        DoubleComplex* input = input_fftw.get();
        DoubleComplex* output_L = output_L_fftw.get();
        DoubleComplex* output_K = output_K_fftw.get();

        for (std::size_t i = 0; i < N; i++) {
            input[i] = 0;
            output_L[i] = 0;
            output_K[i] = 0;
        }

        //left_chol_mat %*% u
        for (std::size_t i = 0; i < L; i++)
            for (std::size_t j = 0; j < std::min(i + 1, Lp); j++) {
                input[i] += leftchol[j * L + i] * u[i - j];
            }

        fftw_plan my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output_L),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (std::size_t i = 0; i < N; i++) {
            input[i] = 0;
        }

        //right_chol_mat %*% u
        for (std::size_t i = 0; i < K; i++)
            for (std::size_t j = 0; j < std::min(i + 1, Kp); j++) {
                input[i] += rightchol[j * K + i] * v[i - j];
            }

        my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output_K),
            FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (std::size_t i = 0; i < N; i++) {
            input[i] = output_L[i] * output_K[i];
        }


        my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output_K),
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        for (std::size_t i = 0; i < N; i++)
            x[i] = output_K[i].real() / N;

        return xexp;
    } catch (std::exception& e) {
        error("%s", e.what());
        return uexp;
    }
}

SEXP hlra_fft_common(SEXP seriesexp, bool direction) {
    std::size_t N = length(seriesexp);

    std::complex<double>* series =
        reinterpret_cast<std::complex<double>*>(COMPLEX(seriesexp));

    try {
        SexpWrapper xexp(allocVector(CPLXSXP, N));

        std::complex<double>* x =
            reinterpret_cast<std::complex<double>*>(COMPLEX(xexp));

        for (std::size_t i = 0; i < N; i++) {
            x[i] = 0;
        }

        std::shared_ptr<DoubleComplex> input_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());
        std::shared_ptr<DoubleComplex> output_fftw(FftwArrayAllocator<DoubleComplex>(N),
            FftwArrayDeleter<DoubleComplex>());

        DoubleComplex* input = input_fftw.get();
        DoubleComplex* output = output_fftw.get();

        for (std::size_t i = 0; i < N; i++) {
            input[i] = series[i];
            output[i] = 0;
        }

        fftw_plan my_plan = fftw_plan_dft_1d(N,
            reinterpret_cast<fftw_complex*>(input),
            reinterpret_cast<fftw_complex*>(output),
            direction ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(my_plan);
        fftw_destroy_plan(my_plan);

        if (direction) {
            for (std::size_t i = 0; i < N; i++) {
                x[i] = output[i];
            }
        } else {
            for (std::size_t i = 0; i < N; i++) {
                x[i] = output[i] / (double)N;
            }
        }

        return xexp;
    } catch (std::exception& e) {
        error("%s", e.what());
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
