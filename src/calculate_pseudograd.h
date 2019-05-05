#ifndef CALCULATE_PSEUDOGRAD_H
#define CALCULATE_PSEUDOGRAD_H

#include <stdlib.h>
#include <complex>
#include <math.h>
#include "rotation_minimizer.h"


template <class Td> class CalculatePseudograd {
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

        void eval_pseudograd_fourier() {
            int i, j;

            fftw_complex* initial_f_vec = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
            fftw_complex* initial_f_vec_fourier = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
            fftw_plan my_plan = fftw_plan_dft_1d(N, initial_f_vec, initial_f_vec_fourier, FFTW_FORWARD, FFTW_ESTIMATE);
            for (i = 0; i < K; i++) {
                std::complex<Td> curval = -signal[i] * std::complex<Td>(cos(-i * alpha), sin(-i * alpha));
                initial_f_vec[i][0] = curval.real();
                initial_f_vec[i][1] = curval.imag();
            }

            for (i = K; i < N; i++) {
                initial_f_vec[i][0] = 0;
                initial_f_vec[i][1] = 0;
            }

            fftw_execute(my_plan);
            fftw_destroy_plan(my_plan);

            std::complex<Td>* current_grad = new std::complex<Td>[N];
            for (i = 0; i < N; i++) {
                current_grad[i].real(initial_f_vec_fourier[i][0]);
                current_grad[i].imag(initial_f_vec_fourier[i][1]);
            }

            int cur_pnt = 0;

            if (tau != 0) {
                for (i = 0; i < N; i++) {
                    pseudograd_fourier[i] = current_grad[i];
                }
                cur_pnt = 1;
            }

            for (j = 1; j < r + 1; j++) {
                for (i = 0; i < N; i++) {
                    current_grad[i] = (current_grad[i] + signal[j - 1]) *
                        std::complex<Td>(cos(alpha + 2 * M_PI * i / N),
                        sin(alpha + 2 * M_PI * i / N)) -
                        signal[j + K - 1] * std::complex<Td>(
                            cos(-alpha * (K - 1) - ((unsigned long long)i * (K - 1)) % N * 2 * M_PI / N),
                            sin(-alpha * (K - 1) - ((unsigned long long)i * (K - 1)) % N * 2 * M_PI / N));
                }
                if (tau != j) {
                    for (i = 0; i < N; i++) {
                        pseudograd_fourier[i + cur_pnt * N] = current_grad[i];
                    }
                    cur_pnt += 1;
                }
            }

            for (j = 0; j < r; j++)
                for (i = 0; i < N; i++) {
                    pseudograd_fourier[i + j * N] /= A_f[i];
                }


            fftw_free(initial_f_vec);
            fftw_free(initial_f_vec_fourier);
        }

    public:
        CalculatePseudograd(int N, int r, Td* glrr,
            std::complex<Td>* A_f, Td* alpha,
            int tau, Td* signal, std::complex<Td>* pseudograd): N(N), r(r),
            K(N-r), glrr(glrr), A_f(A_f), alpha(*alpha), tau(tau),
            signal(signal), pseudograd(pseudograd) {
                pseudograd_fourier = new std::complex<Td>[N * r];
        }

        void doWork() {
            fftw_complex *in, *out;
            fftw_plan my_plan;
            int i, j;

            eval_pseudograd_fourier();

            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
            my_plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

            for (i = 0; i < r; i++) {
                for (j = 0; j < N; j++) {
                    in[j][0] = pseudograd_fourier[i * N + j].real();
                    in[j][1] = pseudograd_fourier[i * N + j].imag();
                }

                fftw_execute(my_plan);
                for (j = 0; j < N; j++) {
                    pseudograd[i * N + j].real(out[j][0] / N);
                    pseudograd[i * N + j].imag(out[j][1] / N);
                }

            }

            fftw_destroy_plan(my_plan);
            fftw_free(in);
            fftw_free(out);
            for (i = 0; i < r; i++) {
                rotate_vector(pseudograd + i * N, N, alpha, pseudograd + i * N);
            }
        }

        ~CalculatePseudograd() {
            delete[] pseudograd_fourier;
        }
};

#endif //CALCULATE_PSEUDOGRAD_H
