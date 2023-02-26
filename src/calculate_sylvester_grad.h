#ifndef CALCULATE_SYLVESTER_GRAD_H
#define CALCULATE_SYLVESTER_GRAD_H

#include <complex>
#include <math.h>
#include <memory>
#include <stdlib.h>

#include "rotation_minimizer.h"
#include "alloc.h"

template <class Td> class CalculateSylvesterGrad {
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
        TComplex* B_f;
        TComplex* glrr_comp;

        void eval_sylvester_grad_fourier() {
            std::shared_ptr<DoubleComplex> initial_f_vec_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());
            std::shared_ptr<DoubleComplex> current_grad_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());

            DoubleComplex* initial_f_vec = initial_f_vec_fftw.get();
            DoubleComplex* current_grad = current_grad_fftw.get();

            fftw_plan my_plan = fftw_plan_dft_1d(N,
                reinterpret_cast<fftw_complex*>(initial_f_vec),
                reinterpret_cast<fftw_complex*>(current_grad),
                FFTW_FORWARD, FFTW_ESTIMATE);

            for (std::size_t i = 0; i < N; i++) {
                initial_f_vec[i] = signal[i] * TComplex(cos(+alpha * i), sin(+alpha * i));
            }

            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);

            for (std::size_t i = 0; i < N; i++) {
               pseudograd_fourier.get()[i] = current_grad[i] / A_f[(N - i) % N];
            }
        }

    public:
        CalculateSylvesterGrad(std::size_t N, std::size_t r, Td* glrr,
            std::complex<Td>* A_f, Td* alpha,
            std::size_t tau, Td* signal, std::complex<Td>* pseudograd): N(N), r(r),
            K(N-r), glrr(glrr), A_f(A_f), alpha(*alpha), tau(tau),
            signal(signal), pseudograd(pseudograd), pseudograd_fourier(
                new TComplex[N], OrdinaryArrayDeleter<TComplex>()) {
        }

        void doWork() {
            eval_sylvester_grad_fourier();

            std::shared_ptr<DoubleComplex> in_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());
            std::shared_ptr<DoubleComplex> out_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());

            DoubleComplex* in = in_fftw.get();
            DoubleComplex* out = out_fftw.get();

            fftw_plan my_plan = fftw_plan_dft_1d(N,
                reinterpret_cast<fftw_complex*>(in),
                reinterpret_cast<fftw_complex*>(out),
                FFTW_BACKWARD, FFTW_ESTIMATE);

            for (std::size_t j = 0; j < N; j++) {
                in[j] = pseudograd_fourier.get()[j];
            }

            fftw_execute(my_plan);
            for (std::size_t j = 0; j < N; j++) {
                pseudograd[j] = out[j] / (double)N;
            }

            DoubleComplex* rotation = out;
            fill_rotation(rotation, N, -alpha);
            for (std::size_t j = 0; j < N; j++) {
                pseudograd[j] *= rotation[j];
            }

            fftw_destroy_plan(my_plan);
        }
};

#endif //CALCULATE_SYLVESTER_GRAD_H
