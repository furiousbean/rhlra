#ifndef CALCULATE_BASIS_H
#define CALCULATE_BASIS_H

#include <complex>
#include <math.h>
#include <memory>
#include <stdlib.h>

#include "alloc.h"
#include "rotation_minimizer.h"


extern "C" {

extern void F77_NAME(ztrtri)(const char* uplo, const char* diag,
             const int* n, Rcomplex* a, const int* lda,
             int* info);
}

const int NO_ORTHOGONALIZATION = 0;
const int ORTHOGONALIZATION    = 1;

template <class Td, int horner_scheme = USUAL_HORNER,
    int orthogonalization = NO_ORTHOGONALIZATION> class CalculateBasis {
    using TComplex = std::complex<Td>;
    using DoubleComplex = std::complex<double>;
    private:
        std::size_t N;
        std::size_t r;
        Td* glrr;
        DoubleComplex* basis;
        DoubleComplex* basis_fourier;
        DoubleComplex* unitroots;
        DoubleComplex* A_f;
        Td& alpha;
        DoubleComplex* qrinvmat;

        void apply_qr_A() {
            const char Nchar = 'N';
            const char Uchar = 'U';
            int info;

            DoubleComplex* data = basis_fourier;
            int size = r;
            int Ni = N;

            int worksz = (4 * size);
            std::shared_ptr<TComplex> work(new TComplex[worksz], OrdinaryArrayDeleter<TComplex>());
            std::shared_ptr<TComplex> tau(new TComplex[size], OrdinaryArrayDeleter<TComplex>());

            F77_CALL(zgeqrf)(&Ni, &size, (Rcomplex*)data, &Ni,
                     (Rcomplex*)tau.get(),
                (Rcomplex*)work.get(), &worksz, &info);

            CheckLapackResult(info, "zgeqrf");

            F77_CALL(ztrtri)(&Uchar, &Nchar, &size, (Rcomplex*)data,
                &Ni, &info);

            CheckLapackResult(info, "ztrtri");

            for (std::size_t i = 0; i < (std::size_t)size; i++) {
                for (std::size_t j = 0; j <= i; j++) {
                    qrinvmat[i * size + j] = data[i * N + j];
                }
                for (std::size_t j = i + 1; j < (std::size_t)size; j++) {
                    qrinvmat[i * size + j] = 0;
                }
            }
        }

        void eval_z_a_fourier() {
            fill_fourier_matrix(basis_fourier, N, r, 1, unitroots);

            for (std::size_t j = 0; j < r; j++)
                for (std::size_t i = 0; i < N; i++) {
                    basis_fourier[i + j * N] /= A_f[i];
                }

            if (orthogonalization == ORTHOGONALIZATION) {
                std::shared_ptr<TComplex> curpoly_shared(new TComplex[r + 1],
                    OrdinaryArrayDeleter<TComplex>());
                TComplex* curpoly = curpoly_shared.get();

                apply_qr_A();

                for (std::size_t i = 0; i < r; i++) {
                    curpoly[0] = 0;
                    for (std::size_t j = 0; j < r; j++)
                        curpoly[j + 1] = qrinvmat[i * r + j];

                    for (std::size_t j = 0; j < N; j++) {
                        basis_fourier[i * N + j] = cpphorner<Td, horner_scheme>(curpoly, r + 1, unitroots[j]) /
                        (A_f[j]);
                    }
                }
            }
        }

    public:
        CalculateBasis(std::size_t N, std::size_t r, Td* glrr, std::complex<Td>* basis,
            std::complex<Td>* basis_fourier, std::complex<Td>* unitroots,
            std::complex<Td>* A_f, Td* alpha, std::complex<Td>* qrinvmat): N(N), r(r), glrr(glrr),
            basis(basis), basis_fourier(basis_fourier), unitroots(unitroots),
            A_f(A_f), alpha(*alpha), qrinvmat(qrinvmat) {
        }

        void doWork() {
            fftw_plan my_plan;

            fill_unitroots(unitroots, N);

            RotationMinimizer<Td, horner_scheme> rm(N, r, glrr, unitroots);
            rm.findMinimum(A_f, alpha);

            eval_z_a_fourier();

            std::shared_ptr<DoubleComplex> in_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());
            std::shared_ptr<DoubleComplex> out_fftw(FftwArrayAllocator<DoubleComplex>(N),
                FftwArrayDeleter<DoubleComplex>());

            DoubleComplex* in = in_fftw.get();
            DoubleComplex* out = out_fftw.get();

            my_plan = fftw_plan_dft_1d(N,
                reinterpret_cast<fftw_complex*>(in),
                reinterpret_cast<fftw_complex*>(out),
                FFTW_BACKWARD, FFTW_ESTIMATE);
            Td sqrtN = std::sqrt(N);
            for (std::size_t i = 0; i < r; i++) {
                for (std::size_t j = 0; j < N; j++) {
                    in[j] = basis_fourier[i * N + j];
                }

                fftw_execute(my_plan);
                for (std::size_t j = 0; j < N; j++) {
                    basis[i * N + j] = out[j] / sqrtN;
                }
            }

            DoubleComplex* rotation = out;
            fill_rotation(rotation, N, alpha);
            for (std::size_t i = 0; i < r; i++) {
                for (std::size_t j = 0; j < N; j++) {
                    basis[i * N + j] *= rotation[j];
                }
            }
            fftw_destroy_plan(my_plan);
        }
};

#endif //CALCULATE_BASIS_H
