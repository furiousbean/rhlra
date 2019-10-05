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

            fftw_complex* initial_f_vec = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
            fftw_complex* initial_f_vec_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
            fftw_plan my_plan = fftw_plan_dft_1d(N, initial_f_vec, initial_f_vec_fourier, FFTW_FORWARD, FFTW_ESTIMATE);
            for (i = 0; i < N; i++) {
                std::complex<Td> curval = signal[i] * std::complex<Td>(cos(+i * alpha), sin(+i * alpha));
                initial_f_vec[i][0] = curval.real();
                initial_f_vec[i][1] = curval.imag();
            }

            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);

            std::complex<Td>* current_grad = new std::complex<Td>[N];
            for (i = 0; i < N; i++) {
                current_grad[i].real(initial_f_vec_fourier[i][0]);
                current_grad[i].imag(initial_f_vec_fourier[i][1]);
            }

            for (i = 0; i < N; i++) {
                pseudograd_fourier[i] = current_grad[i];
            }

            for (i = 0; i < N; i++) {
               pseudograd_fourier[i] /= A_f[(N - i) % N];
            }


            fftw_free(initial_f_vec);
            fftw_free(initial_f_vec_fourier);
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
            fftw_complex *in, *out;
            fftw_plan my_plan;
            int j;

            eval_sylvester_grad_fourier();

            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            my_plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

            for (j = 0; j < N; j++) {
                in[j][0] = pseudograd_fourier[j].real();
                in[j][1] = pseudograd_fourier[j].imag();
            }

            fftw_execute(my_plan);
            for (j = 0; j < N; j++) {
                pseudograd[j].real(out[j][0] / N);
                pseudograd[j].imag(out[j][1] / N);
            }

            fftw_destroy_plan(my_plan);
            fftw_free(in);
            fftw_free(out);
            rotate_vector(pseudograd, N, -alpha, pseudograd);
        }

        ~CalculateSylvesterGrad() {
            delete[] pseudograd_fourier;
        }
};

#endif //CALCULATE_SYLVESTER_GRAD_H
