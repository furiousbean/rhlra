#ifndef CALCULATE_PSEUDOGRAD_H
#define CALCULATE_PSEUDOGRAD_H

#include <complex>
#include <memory>
#include <math.h>
#include <vector>
#include <stdlib.h>

#include "rotation_minimizer.h"
#include "alloc.h"


template <class Td> class CalculatePseudograd {
    using TComplex = std::complex<Td>;
    using DoubleComplex = std::complex<double>;
    private:
        std::size_t N;
        std::size_t r;
        std::size_t K;
        Td* glrr;
        TComplex* A_f;
        Td& alpha;
        std::size_t tau;
        Td* signal;
        TComplex* pseudograd;
        std::shared_ptr<TComplex> pseudograd_fourier;

        void eval_pseudograd_fourier() {
            std::shared_ptr<fftw_complex> initial_f_vec_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());
            std::shared_ptr<fftw_complex> current_grad_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());

            DoubleComplex* initial_f_vec = reinterpret_cast<DoubleComplex*>(initial_f_vec_fftw.get());
            DoubleComplex* current_grad = reinterpret_cast<DoubleComplex*>(current_grad_fftw.get());

            fftw_plan my_plan = fftw_plan_dft_1d(N, initial_f_vec_fftw.get(),
                current_grad_fftw.get(), FFTW_FORWARD, FFTW_ESTIMATE);
            for (std::size_t i = 0; i < K; i++) {
                initial_f_vec[i] = -signal[i] * TComplex(cos(-alpha * i), sin(-alpha * i));
            }

            for (std::size_t i = K; i < N; i++) {
                initial_f_vec[i] = 0;
            }

            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);

            std::size_t cur_pnt = 0;

            if (tau != 0) {
                for (std::size_t i = 0; i < N; i++) {
                    pseudograd_fourier.get()[i] = current_grad[i];
                }
                cur_pnt = 1;
            }

            std::vector<TComplex> back_rotation(N);
            std::vector<TComplex> front_rotation(N);

            for (std::size_t i = 0; i < N; i++) {
                back_rotation[i] = TComplex(
                    cos(alpha + 2 * M_PI * i / N), sin(alpha + 2 * M_PI * i / N));
                front_rotation[i] = TComplex(
                    cos(-alpha * (K - 1) - ((unsigned long long)i * (K - 1)) % N * 2 * M_PI / N),
                    sin(-alpha * (K - 1) - ((unsigned long long)i * (K - 1)) % N * 2 * M_PI / N));
            }

            for (std::size_t j = 1; j < r + 1; j++) {
                for (std::size_t i = 0; i < N; i++) {
                    current_grad[i] = (current_grad[i] + signal[j - 1]) * back_rotation[i]
                         - signal[j + K - 1] * front_rotation[i];
                }
                if (tau != j) {
                    for (std::size_t i = 0; i < N; i++) {
                        pseudograd_fourier.get()[i + cur_pnt * N] = current_grad[i];
                    }
                    cur_pnt += 1;
                }
            }

            for (std::size_t j = 0; j < r; j++) {
                for (std::size_t i = 0; i < N; i++) {
                    pseudograd_fourier.get()[i + j * N] /= A_f[i];
                }
            }
        }

    public:
        CalculatePseudograd(std::size_t N, std::size_t r, Td* glrr,
            TComplex* A_f, Td* alpha,
            std::size_t tau, Td* signal, TComplex* pseudograd): N(N), r(r),
            K(N-r), glrr(glrr), A_f(A_f), alpha(*alpha), tau(tau),
            signal(signal), pseudograd(pseudograd), pseudograd_fourier(
                new TComplex[N * r], OrdinaryArrayDeleter<TComplex>()) {
        }

        void doWork() {
            eval_pseudograd_fourier();

            std::shared_ptr<fftw_complex> in_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());
            std::shared_ptr<fftw_complex> out_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());

            DoubleComplex* in = reinterpret_cast<DoubleComplex*>(in_fftw.get());
            DoubleComplex* out = reinterpret_cast<DoubleComplex*>(out_fftw.get());

            fftw_plan my_plan = fftw_plan_dft_1d(N, in_fftw.get(), out_fftw.get(),
                FFTW_BACKWARD, FFTW_ESTIMATE);

            for (std::size_t i = 0; i < r; i++) {
                for (std::size_t j = 0; j < N; j++) {
                    in[j] = pseudograd_fourier.get()[i * N + j];
                }

                fftw_execute(my_plan);
                for (std::size_t j = 0; j < N; j++) {
                    pseudograd[i * N + j] = out[j] / (double)N;
                }

            }

            DoubleComplex* rotation = out;
            fill_rotation(rotation, N, alpha);
            for (std::size_t i = 0; i < r; i++) {
                for (std::size_t j = 0; j < N; j++) {
                    pseudograd[i * N + j] *= rotation[j];
                }
            }

            fftw_destroy_plan(my_plan);
        }
};

#endif //CALCULATE_PSEUDOGRAD_H
