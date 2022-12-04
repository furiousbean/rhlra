#ifndef ROTATION_MINIMIZER_H
#define ROTATION_MINIMIZER_H

#include <stdlib.h>
#include <complex>
#include <complex.h>
#include <math.h>
#include "chorner.h"
#include <R.h>
#include <Rinternals.h>

const int USUAL_HORNER  = 0;
const int COMPENSATED_HORNER  = 1;

const double M_PHI = 0.618033988749894848204586834365638117720309179; //2/(1 + sqrt(5));
const int SEARCH_INT = 5;

template <typename Td, int horner_scheme> std::complex<Td> cpphorner(const std::complex<Td>* p, size_t n,
                           std::complex<Td> x) {
    std::complex<Td> sum(0);
    for (int i = n-1; i >= 0; i--) {
        sum = sum * x + p[i];
    }
    return sum;
}

template<> std::complex<double> cpphorner<double, COMPENSATED_HORNER>(const std::complex<double>* p, size_t n,
                           std::complex<double> x) {
    _Complex double result = cchorner(reinterpret_cast<const _Complex double *>(p),
        n, reinterpret_cast<_Complex double&>(x));

    return std::complex<double>(reinterpret_cast<std::complex<double>&>(result));
}

template <class Td = double> void fill_unitroots(std::complex<Td>* data, int N) {
    int i;
    Td cur_angle;

    for (i = 0; i < N; i++) {
        cur_angle = 2 * M_PI * i / N;
        data[i] = { cos(cur_angle), sin(cur_angle) };
    }
}

template <class Td = double> void fill_fourier_matrix(std::complex<Td>* data, int N, int size, int begin,
                                                      std::complex<Td>* unitroots) {
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < N; j++) {
            data[i * N + j] = unitroots[((unsigned long long)(i + begin) * j) % N];
        }
    }
}

template <class Td = double> void rotate_vector(std::complex<Td>* input, int N, double alpha,
                   std::complex<Td>* output) {
    int i;
    for (i = 0; i < N; i++) {
        output[i] = input[i] * std::complex<Td>(cos(alpha * i), sin(alpha * i));
    }
}

template <int horner_scheme> int getSearchIt() {
    return 15;
}

template <> int getSearchIt<COMPENSATED_HORNER>() {
    return 20;
}

template <class Td, int horner_scheme = USUAL_HORNER> class RotationMinimizer {
    private:
        int N;
        int r;
        std::complex<Td>* glrr;
        std::complex<Td>* unitroots;
        int* active_indices;
        int* active_indices_as_bool;
        int aisz;
        int SEARCH_IT;

        Td conv_eval_alpha(Td alpha) {
            std::complex<Td> sum(0);
            int min_flag = 1;
            Td minimum(0);
            Td value;
            int i;

            std::complex<Td> cur_rot(cos(alpha), sin(alpha));
            for (i = 0; i < aisz; i++) {
                sum = cpphorner<Td, horner_scheme>(glrr, r + 1, unitroots[active_indices[i]] * cur_rot);
                value = (sum * conj(sum)).real();
                if (min_flag || value < minimum) {
                    minimum = value;
                    min_flag = 0;
                }
            }

            // Rprintf("A %lf %lf %ld\n", alpha, -log(minimum), aisz);
            return -minimum;
        }

        Td conv_min_alpha_int(Td li, Td ri) {
            Td lc, rc, lval, rval;
            int i;
            lc = ri - (ri - li) * M_PHI;
            rc = li + (ri - li) * M_PHI;
            lval = conv_eval_alpha(lc);
            rval = conv_eval_alpha(rc);
            double answer = 0;
            for (i = 0; i < SEARCH_IT; i++) {
                if (lval < rval) {
                    ri = rc;
                    rc = lc;
                    rval = lval;
                    lc = ri - (ri - li) * M_PHI;
                    lval = conv_eval_alpha(lc);
                    answer = lc;
                } else {
                    li = lc;
                    lc = rc;
                    lval = rval;
                    rc = li + (ri - li) * M_PHI;
                    rval = conv_eval_alpha(rc);
                    answer = rc;
                }
            }
            // Rprintf("rc - lc = %e\n", ri - li);
            return answer;
        }

        Td conv_min_alpha(Td* minimum) {
            Td lc = 0, rc, cur_value, cur_alpha;
            Td alpha = 0;
            rc = 2 * M_PI / N;
            int min_flag = 1;

            for (int i = 0; i < SEARCH_INT; i++) {
                cur_alpha = conv_min_alpha_int(lc + i * (rc - lc) / SEARCH_INT,
                    lc + (i + 1) * (rc - lc) / SEARCH_INT);

                //HACK!!!
                //cur_alpha = M_PI / N;

                cur_value = conv_eval_alpha(cur_alpha);
                if (min_flag || cur_value < *minimum) {
                    *minimum = cur_value;
                    alpha = cur_alpha;
                    min_flag = 0;
                }
            }
            return alpha;
        }

    public:
        RotationMinimizer(int N, int r, Td* rglrr, std::complex<Td>* unitroots) : N(N),
        r(r), unitroots(unitroots), aisz(0) {
            glrr = new std::complex<Td>[r + 1];
            for (int i = 0; i < r + 1; i++)
                glrr[i] = rglrr[i];

            active_indices = new int[N];
            active_indices_as_bool = new int[N];

            SEARCH_IT = getSearchIt<horner_scheme>();
        }

        ~RotationMinimizer() {
            delete[] glrr;
            delete[] active_indices;
            delete[] active_indices_as_bool;
        }

        void findMinimum(std::complex<Td>* A_f, Td& alpha) {
            int i;
            Td minimum_part(0), cur_minimum(0), value;

            // alpha = M_PI / N;
            alpha = 0;

            int new_idx = -1;

            std::complex<Td> cur_rot(cos(alpha), sin(alpha));
            for (i = 0; i < N; i++) {
                active_indices[i] = 0;
                active_indices_as_bool[i] = 0;
                A_f[i] = cpphorner<Td, horner_scheme>(glrr, r + 1, unitroots[i] * cur_rot);
                value = (A_f[i] * conj(A_f[i])).real();
                if (new_idx == -1 || value < cur_minimum) {
                        new_idx = i;
                        cur_minimum = value;
                    }
            }

            active_indices_as_bool[new_idx] = 1;
            active_indices_as_bool[(new_idx + N - 1) % N] = 1;
            // active_indices_as_bool[(new_idx + 1) % N] = 1;

            active_indices[0] = new_idx;
            active_indices[1] = (new_idx + N - 1) % N;
            // active_indices[2] = (new_idx + 1) % N;
            aisz = 2;

            do {
                // Rprintf("B %d\n", aisz);
                alpha = conv_min_alpha(&minimum_part);
                minimum_part = -minimum_part;

                cur_rot = { cos(alpha), sin(alpha) };

                int min_flag = 1;

                for (i = 0; i < N; i++) {
                    A_f[i] = cpphorner<Td, horner_scheme>(glrr, r + 1, unitroots[i] * cur_rot);
                    value = (A_f[i] * conj(A_f[i])).real();
                    if (min_flag || value < cur_minimum) {
                        new_idx = i;
                        cur_minimum = value;
                        min_flag = 0;
                    }
                }

                // Rprintf("minimum_part: %f, cur_minimum %f\n", minimum_part, cur_minimum);
                if (minimum_part > cur_minimum) {
                    active_indices_as_bool[new_idx] = 1;
                    active_indices_as_bool[(new_idx + N - 1) % N] = 1;
                    active_indices_as_bool[(new_idx + 1) % N] = 1;
                    // Rprintf("OK try again\n");
                    aisz = 0;
                    for (i = 0; i < N; i++) {
                        if (active_indices_as_bool[i]) {
                            active_indices[aisz++] = i;
                        }
                    }
                } else {
                    // Rprintf("DONE\n");
                    break;
                }

            } while (1);
        }
};

#endif //ROTATION_MINIMIZER_H
