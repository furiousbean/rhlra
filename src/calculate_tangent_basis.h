#ifndef CALCULATE_TANGENT_BASIS_H
#define CALCULATE_TANGENT_BASIS_H

#include <stdlib.h>
#include <complex>
#include <math.h>
#include "calculate_basis.h"

template <class Td, int horner_scheme = USUAL_HORNER> class CalculateTangentBasis {
    private:
        int N;
        int r;
        Td* glrr;
        std::complex<Td>* basis;
        std::complex<Td>* basis_fourier;
        std::complex<Td>* unitroots;
        std::complex<Td>* A_f;
        Td alpha;
        std::complex<Td>* qrinvmat;

        void apply_qr_A2(std::complex<Td>* data, std::complex<Td>* vmat) {
            std::complex<Td>* work;
            std::complex<Td>* tau;
            int* perm;
            Td* rwork;
            char Nchar = 'N';
            char Uchar = 'U';
            int i, j, info, worksz, sz2;

            int size = r;

            sz2 = 2 * size;
            worksz = (4 * sz2);
            work = new std::complex<Td>[worksz];
            tau = new std::complex<Td>[sz2];
            perm = new int[sz2];
            rwork = new Td[2 * sz2];

            F77_CALL(zgeqp3)(&N, &sz2, (Rcomplex*)data, &N,
                     perm, (Rcomplex*)tau,
                (Rcomplex*)work, &worksz, rwork, &info);

            F77_CALL(ztrtri)(&Uchar, &Nchar, &size, (Rcomplex*)data,
                &N, &info);

            for (i = 0; i < sz2 * sz2; i++) {
                vmat[i] = 0;
            }

            for (i = 0; i < size; i++) {
                for (j = 0; j <= i; j++) {
                    vmat[i*sz2 + (perm[j] - 1)] = data[i * N + j];
                }
            }

            delete[] perm;
            delete[] tau;
            delete[] work;
            delete[] rwork;
        }

        void arithm_addit_proj(std::complex<Td>* basis, int N, int size,
                               std::complex<Td>* data, int p) {
            std::complex<Td>* coefs;
            std::complex<Td> sum;
            int i, j, k;

            coefs = new std::complex<Td>[size * p]; //calloc(size * p, sizeof(double _Complex));
            for (i = 0; i < p; i++)
                for (j = 0; j < size; j++) {
                    sum = 0;
                    for (k = 0; k < N; k++) {
                        sum += conj(basis[j * N + k]) * data[i * N + k];
                    }
                    coefs[i * size + j] = sum;
                }

            for (i = 0; i < p; i++)
                for (j = 0; j < size; j++) {
                    for (k = 0; k < N; k++) {
                        data[i * N + k] -= basis[j * N + k] * coefs[i * size + j];
                    }
                }

            delete[] coefs;
        }

        void eval_z_a2_fourier() {
            std::complex<Td>* preZ_a2;
            std::complex<Td>* vmat;
            std::complex<Td>* curpoly;
            std::complex<Td>* qrdata;
            int i, j;

            std::complex<Td>* Z_a = basis_fourier;
            preZ_a2 = new std::complex<Td>[2 * r * N];
            vmat = new std::complex<Td>[4 * r * r];
            curpoly = new std::complex<Td>[2 * r + 1];

            fill_fourier_matrix(preZ_a2, N, 2 * r, 1, unitroots);

            qrdata = Z_a;

            for (j = 0; j < 2 * r; j++)
                for (i = 0; i < N; i++) {
                    preZ_a2[i + j * N] /= (A_f[i] * A_f[i]);
                }

            arithm_addit_proj(qrdata, N, r, preZ_a2, 2 * r);

            apply_qr_A2(preZ_a2, vmat);

            for (i = 0; i < r; i++) {
                curpoly[0] = 0;
                for (j = 0; j < 2 * r; j++)
                    curpoly[j + 1] = vmat[2 * r * i + j];

                for (j = 0; j < N; j++) {
                    Z_a[(r + i) * N + j] =
                    cpphorner<double, horner_scheme>(curpoly, 2 * r + 1, unitroots[j]) /
                    (A_f[j] * A_f[j]);
                }
            }

            arithm_addit_proj(qrdata, N, r, Z_a + N * r, r);

            delete[] preZ_a2;
            delete[] vmat;
            delete[] curpoly;
        }

    public:
        CalculateTangentBasis(int N, int r, Td* glrr, std::complex<Td>* basis):
        N(N), r(r), glrr(glrr), basis(basis), alpha(0) {
                basis_fourier = new std::complex<Td>[2 * N * r];
                unitroots = new std::complex<Td>[N];
                A_f = new std::complex<Td>[N];
                qrinvmat = new std::complex<Td>[r * r];
        }

        ~CalculateTangentBasis() {
            delete[] basis_fourier;
            delete[] unitroots;
            delete[] A_f;
            delete[] qrinvmat;
        }

        void doWork() {
            double alpha;
            double sqrtN;
            fftw_plan my_plan;
            int i, j;

            CalculateBasis<Td, horner_scheme, ORTHOGONALIZATION> cb(N, r, glrr,
                    basis, basis_fourier, unitroots, A_f, &alpha, qrinvmat);
            cb.doWork();

            eval_z_a2_fourier();

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

            sqrtN = std::sqrt(N);
            for (i = r; i < 2 * r; i++) {
                for (j = 0; j < N; j++) {
                    in[j] = basis_fourier[i * N + j];
                }

                fftw_execute(my_plan);
                for (j = 0; j < N; j++) {
                    basis[i * N + j] = out[j] / sqrtN;
                }
            }

            fftw_destroy_plan(my_plan);
            fftw_free(in_fftw);
            fftw_free(out_fftw);
            for (i = r; i < 2 * r; i++) {
                rotate_vector(basis + i * N, N, alpha, basis + i * N);
            }
        }
};

#endif //CALCULATE_TANGENT_BASIS_H
