#ifndef CALCULATE_SYLVESTER_GRAD_H
#define CALCULATE_SYLVESTER_GRAD_H

#include <stdlib.h>
#include <complex>
#include <math.h>
#include "rotation_minimizer.h"


template <class Td> class CalculateSylvesterGrad {
    private:
        int N;
        int r;
        int K;
        Td* glrr;
        std::complex<Td>* A_f;
        Td& alpha;
        int tau;
        Td* signal;
        std::complex<Td>* pseudograd;
        std::complex<Td>* pseudograd_fourier;
        std::complex<Td>* B_f;
        std::complex<Td>* glrr_comp;

        void eval_sylvester_grad_fourier() {
            int i;

            fftw_complex* initial_f_vec_fftw = reinterpret_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * N));
            fftw_complex* current_grad_fftw = reinterpret_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * N));

            std::complex<double>* initial_f_vec =
                reinterpret_cast<std::complex<double>*>(initial_f_vec_fftw);
            std::complex<double>* current_grad =
                reinterpret_cast<std::complex<double>*>(current_grad_fftw);

            fftw_plan my_plan = fftw_plan_dft_1d(N, initial_f_vec_fftw, current_grad_fftw,
                FFTW_FORWARD, FFTW_ESTIMATE);

            for (i = 0; i < N; i++) {
                initial_f_vec[i] = signal[i] * std::complex<Td>(cos(+i * alpha), sin(+i * alpha));
            }

            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);

            for (i = 0; i < N; i++) {
               pseudograd_fourier[i] = current_grad[i] / A_f[(N - i) % N];
            }

            fftw_free(initial_f_vec_fftw);
            fftw_free(current_grad_fftw);
        }

    public:
        CalculateSylvesterGrad(int N, int r, Td* glrr,
            std::complex<Td>* A_f, Td* alpha,
            int tau, Td* signal, std::complex<Td>* pseudograd): N(N), r(r),
            K(N-r), glrr(glrr), A_f(A_f), alpha(*alpha), tau(tau),
            signal(signal), pseudograd(pseudograd) {
                pseudograd_fourier = new std::complex<Td>[N];
        }

        void doWork() {
            fftw_plan my_plan;
            int j;

            eval_sylvester_grad_fourier();

            fftw_complex *in_fftw = reinterpret_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * N));
            fftw_complex *out_fftw = reinterpret_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * N));

            std::complex<double>* in =
                reinterpret_cast<std::complex<double>*>(in_fftw);
            std::complex<double>* out =
                reinterpret_cast<std::complex<double>*>(out_fftw);

            my_plan = fftw_plan_dft_1d(N, in_fftw, out_fftw,
                FFTW_BACKWARD, FFTW_ESTIMATE);

            for (j = 0; j < N; j++) {
                in[j] = pseudograd_fourier[j];
            }

            fftw_execute(my_plan);
            for (j = 0; j < N; j++) {
                pseudograd[j] = out[j] / (double)N;
            }

            std::complex<double>* rotation = out;
            fill_rotation(rotation, N, -alpha);
            for (j = 0; j < N; j++) {
                pseudograd[j] *= rotation[j];
            }

            fftw_destroy_plan(my_plan);
            fftw_free(in_fftw);
            fftw_free(out_fftw);
        }

        ~CalculateSylvesterGrad() {
            delete[] pseudograd_fourier;
        }
};

#endif //CALCULATE_SYLVESTER_GRAD_H
