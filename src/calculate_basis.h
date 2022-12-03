#ifndef CALCULATE_BASIS_H
#define CALCULATE_BASIS_H

#include <stdlib.h>
#include <complex>
#include <math.h>
#include "rotation_minimizer.h"

extern "C" {

extern void F77_NAME(ztrtri)(const char* uplo, const char* diag,
             const int* n, Rcomplex* a, const int* lda,
             int* info);
}

const int NO_ORTHOGONALIZATION  = 0;
const int ORTHOGONALIZATION  = 1;

template <class Td, int horner_scheme = USUAL_HORNER,
int orthogonalization = NO_ORTHOGONALIZATION> class CalculateBasis {
    private:
        int N;
        int r;
        Td* glrr;
        std::complex<Td>* basis;
        std::complex<Td>* basis_fourier;
        std::complex<Td>* unitroots;
        std::complex<Td>* A_f;
        Td& alpha;
        std::complex<Td>* qrinvmat;

        void apply_qr_A() {
            std::complex<double>* work;
            std::complex<double>* tau;
            char Nchar = 'N';
            char Uchar = 'U';
            int i, j, info, worksz;

            std::complex<double>* data = basis_fourier;
            int size = r;

            worksz = (4 * size);
            work = new std::complex<double>[worksz];
            tau = new std::complex<double>[size];

            F77_CALL(zgeqrf)(&N, &size, (Rcomplex*)data, &N,
                     (Rcomplex*)tau,
                (Rcomplex*)work, &worksz, &info);

            F77_CALL(ztrtri)(&Uchar, &Nchar, &size, (Rcomplex*)data,
                &N, &info);

            for (i = 0; i < size; i++) {
                for (j = 0; j <= i; j++) {
                    qrinvmat[i*size + j] = data[i * N + j];
                }
                for (j = i + 1; j < size; j++) {
                    qrinvmat[i*size + j] = 0;
                }
            }

            delete[] tau;
            delete[] work;
        }

        void eval_z_a_fourier() {
            int i, j;

            fill_fourier_matrix(basis_fourier, N, r, 1, unitroots);

            for (j = 0; j < r; j++)
                for (i = 0; i < N; i++) {
                    basis_fourier[i + j * N] /= A_f[i];
                }

            if (orthogonalization == ORTHOGONALIZATION) {
                std::complex<Td>* curpoly;
                curpoly = new std::complex<Td>[r + 1];

                apply_qr_A();

                for (i = 0; i < r; i++) {
                    curpoly[0] = 0;
                    for (j = 0; j < r; j++)
                        curpoly[j + 1] = qrinvmat[i * r + j];

                    for (j = 0; j < N; j++) {
                        basis_fourier[i * N + j] = cpphorner<Td, horner_scheme>(curpoly, r + 1, unitroots[j]) /
                        (A_f[j]);
                    }
                }

                delete[] curpoly;
            }
        }

    public:
        CalculateBasis(int N, int r, Td* glrr, std::complex<Td>* basis,
            std::complex<Td>* basis_fourier, std::complex<Td>* unitroots,
            std::complex<Td>* A_f, Td* alpha, std::complex<Td>* qrinvmat): N(N), r(r), glrr(glrr),
            basis(basis), basis_fourier(basis_fourier), unitroots(unitroots),
            A_f(A_f), alpha(*alpha), qrinvmat(qrinvmat) {
        }

        void doWork() {
            Td sqrtN;
            fftw_plan my_plan;
            int i, j;

            fill_unitroots(unitroots, N);

            RotationMinimizer<Td, horner_scheme> rm(N, r, glrr, unitroots);
            rm.findMinimum(A_f, alpha);

            eval_z_a_fourier();

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
            for (i = 0; i < r; i++) {
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
            for (i = 0; i < r; i++) {
                rotate_vector(basis + i * N, N, alpha, basis + i * N);
            }
        }
};

#endif //CALCULATE_BASIS_H
