#ifndef CALCULATE_TANGENT_BASIS_H
#define CALCULATE_TANGENT_BASIS_H

#include <complex>
#include <math.h>
#include <memory>
#include <stdlib.h>

#include "calculate_basis.h"
#include "alloc.h"


template <class Td, int horner_scheme = USUAL_HORNER> class CalculateTangentBasis {
    using TComplex = std::complex<Td>;
    using DoubleComplex = std::complex<double>;
    private:
        int N;
        int r;
        Td* glrr;
        TComplex* basis;
        std::shared_ptr<TComplex> basis_fourier;
        std::shared_ptr<TComplex> unitroots;
        std::shared_ptr<TComplex> A_f;
        Td alpha;
        std::shared_ptr<TComplex> qrinvmat;

        void apply_qr_A2(TComplex* data, TComplex* vmat) {
            const char Nchar = 'N';
            const char Uchar = 'U';
            int info;

            int size = r;

            int sz2 = 2 * size;
            int worksz = (4 * sz2);
            std::shared_ptr<TComplex> work(new TComplex[worksz], OrdinaryArrayDeleter<TComplex>());
            std::shared_ptr<TComplex> tau(new TComplex[sz2], OrdinaryArrayDeleter<TComplex>());
            std::shared_ptr<int> perm(new int[sz2], OrdinaryArrayDeleter<int>());
            std::shared_ptr<Td> rwork(new Td[2 * sz2], OrdinaryArrayDeleter<Td>());

            F77_CALL(zgeqp3)(&N, &sz2, (Rcomplex*)data, &N,
                     perm.get(), (Rcomplex*)tau.get(),
                (Rcomplex*)work.get(), &worksz, rwork.get(), &info);

            CheckLapackResult(info, "zgeqp3");

            F77_CALL(ztrtri)(&Uchar, &Nchar, &size, (Rcomplex*)data,
                &N, &info);

            CheckLapackResult(info, "ztrtri");

            for (int i = 0; i < sz2 * sz2; i++) {
                vmat[i] = 0;
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j <= i; j++) {
                    vmat[i*sz2 + (perm.get()[j] - 1)] = data[i * N + j];
                }
            }
        }

        void arithm_addit_proj(TComplex* basis, int N, int size,
                               TComplex* data, int p) {
            TComplex sum;

            std::shared_ptr<TComplex> coefs(new TComplex[size * p], OrdinaryArrayDeleter<TComplex>());
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < size; j++) {
                    sum = 0;
                    for (int k = 0; k < N; k++) {
                        sum += conj(basis[j * N + k]) * data[i * N + k];
                    }
                    coefs.get()[i * size + j] = sum;
                }
            }

            for (int i = 0; i < p; i++) {
                for (int j = 0; j < size; j++) {
                    for (int k = 0; k < N; k++) {
                        data[i * N + k] -= basis[j * N + k] * coefs.get()[i * size + j];
                    }
                }
            }
        }

        void eval_z_a2_fourier() {
            std::shared_ptr<TComplex> preZ_a2(new TComplex[2 * r * N], OrdinaryArrayDeleter<TComplex>());
            std::shared_ptr<TComplex> vmat(new TComplex[4 * r * r], OrdinaryArrayDeleter<TComplex>());
            std::shared_ptr<TComplex> curpoly(new TComplex[2 * r + 1], OrdinaryArrayDeleter<TComplex>());

            fill_fourier_matrix(preZ_a2.get(), N, 2 * r, 1, unitroots.get());

            TComplex* Z_a = basis_fourier.get();
            TComplex* qrdata = Z_a;

            for (int j = 0; j < 2 * r; j++) {
                for (int i = 0; i < N; i++) {
                    preZ_a2.get()[i + j * N] /= (A_f.get()[i] * A_f.get()[i]);
                }
            }

            arithm_addit_proj(qrdata, N, r, preZ_a2.get(), 2 * r);

            apply_qr_A2(preZ_a2.get(), vmat.get());

            for (int i = 0; i < r; i++) {
                curpoly.get()[0] = 0;
                for (int j = 0; j < 2 * r; j++)
                    curpoly.get()[j + 1] = vmat.get()[2 * r * i + j];

                for (int j = 0; j < N; j++) {
                    Z_a[(r + i) * N + j] =
                    cpphorner<double, horner_scheme>(curpoly.get(), 2 * r + 1, unitroots.get()[j]) /
                    (A_f.get()[j] * A_f.get()[j]);
                }
            }

            arithm_addit_proj(qrdata, N, r, Z_a + N * r, r);
        }

    public:
        CalculateTangentBasis(int N, int r, Td* glrr, std::complex<Td>* basis):
            N(N), r(r), glrr(glrr), basis(basis),
            basis_fourier(new TComplex[2 * N * r], OrdinaryArrayDeleter<TComplex>()),
            unitroots(new TComplex[N], OrdinaryArrayDeleter<TComplex>()),
            A_f(new TComplex[N], OrdinaryArrayDeleter<TComplex>()),
            alpha(0),
            qrinvmat(new TComplex[r * r], OrdinaryArrayDeleter<TComplex>()) {
        }

        void doWork() {
            double alpha;

            CalculateBasis<Td, horner_scheme, ORTHOGONALIZATION> cb(N, r, glrr,
                    basis, basis_fourier.get(), unitroots.get(), A_f.get(), &alpha,
                    qrinvmat.get());
            cb.doWork();

            eval_z_a2_fourier();

            std::shared_ptr<fftw_complex> in_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());
            std::shared_ptr<fftw_complex> out_fftw(FftwArrayAllocator<fftw_complex>(N),
                FftwArrayDeleter<fftw_complex>());

            DoubleComplex* in = reinterpret_cast<DoubleComplex*>(in_fftw.get());
            DoubleComplex* out = reinterpret_cast<DoubleComplex*>(out_fftw.get());

            fftw_plan my_plan = fftw_plan_dft_1d(N, in_fftw.get(), out_fftw.get(),
                FFTW_BACKWARD, FFTW_ESTIMATE);

            double sqrtN = std::sqrt(N);
            for (int i = r; i < 2 * r; i++) {
                for (int j = 0; j < N; j++) {
                    in[j] = basis_fourier.get()[i * N + j];
                }

                fftw_execute(my_plan);
                for (int j = 0; j < N; j++) {
                    basis[i * N + j] = out[j] / sqrtN;
                }
            }

            DoubleComplex* rotation = out;
            fill_rotation(rotation, N, alpha);
            for (int i = 0; i < 2 * r; i++) {
                for (int j = 0; j < N; j++) {
                    basis[i * N + j] *= rotation[j];
                }
            }

            fftw_destroy_plan(my_plan);
        }
};

#endif //CALCULATE_TANGENT_BASIS_H
