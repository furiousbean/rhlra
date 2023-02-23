#ifndef ROTATION_MINIMIZER_H
#define ROTATION_MINIMIZER_H

#include <complex>
#include <complex.h>
#include <math.h>
#include <memory>
#include <vector>
#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>

#include "alloc.h"
#include "chorner.h"


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

template <class Td = double> void fill_rotation(std::complex<Td>* data, int N, double alpha) {
    for (int i = 0; i < N; i++) {
        data[i] = std::complex<Td>(cos(alpha * i), sin(alpha * i));
    }
}

template <int horner_scheme> int getSearchIt() {
    return 15;
}

template <> int getSearchIt<COMPENSATED_HORNER>() {
    return 20;
}

template <class Td, int horner_scheme = USUAL_HORNER> class RotationMinimizer {
    using TComplex = std::complex<Td>;
    private:
        int N;
        int r;
        std::shared_ptr<TComplex> glrr;
        TComplex* unitroots;
        std::vector<int> active_indices;
        std::vector<bool> active_indices_as_bool;
        int SEARCH_IT;

        Td norm(TComplex x) {
            return x.real() * x.real() + x.imag() * x.imag();
        };

        Td conv_eval_alpha(Td alpha) {
            TComplex sum(0);
            bool min_flag = true;
            Td minimum(0);
            Td value;

            TComplex cur_rot(cos(alpha), sin(alpha));
            for (int i = 0; i < (int)active_indices.size(); i++) {
                sum = cpphorner<Td, horner_scheme>(glrr.get(), r + 1, unitroots[active_indices[i]] * cur_rot);
                value = norm(sum);
                if (min_flag || value < minimum) {
                    minimum = value;
                    min_flag = false;
                }
            }

            // Rprintf("A %lf %lf %ld\n", alpha, -log(minimum), active_indices.size());
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
            bool min_flag = true;

            for (int i = 0; i < SEARCH_INT; i++) {
                cur_alpha = conv_min_alpha_int(lc + i * (rc - lc) / SEARCH_INT,
                    lc + (i + 1) * (rc - lc) / SEARCH_INT);

                //HACK!!!
                //cur_alpha = M_PI / N;

                cur_value = conv_eval_alpha(cur_alpha);
                if (min_flag || cur_value < *minimum) {
                    *minimum = cur_value;
                    alpha = cur_alpha;
                    min_flag = false;
                }
            }
            return alpha;
        }

    public:
        RotationMinimizer(int N, int r, Td* rglrr, std::complex<Td>* unitroots) : N(N),
            r(r), glrr(new TComplex[r+1], OrdinaryArrayDeleter<TComplex>()),
            unitroots(unitroots), active_indices_as_bool(N) {
            for (int i = 0; i < r + 1; i++)
                glrr.get()[i] = rglrr[i];

            SEARCH_IT = getSearchIt<horner_scheme>();
        }

        void findMinimum(std::complex<Td>* A_f, Td& alpha) {
            int i;
            Td minimum_part(0), cur_minimum(0), value;

            // alpha = M_PI / N;
            alpha = 0;

            int new_idx = -1;

            std::complex<Td> cur_rot(cos(alpha), sin(alpha));
            for (i = 0; i < N; i++) {
                A_f[i] = cpphorner<Td, horner_scheme>(glrr.get(), r + 1, unitroots[i] * cur_rot);
                value = norm(A_f[i]);
                if (new_idx == -1 || value < cur_minimum) {
                    new_idx = i;
                    cur_minimum = value;
                }
            }

            active_indices.push_back((new_idx + N - 1) % N);
            active_indices.push_back(new_idx);
            // active_indices.push_back((new_idx + 1) % N);

            active_indices_as_bool[(new_idx + N - 1) % N] = true;
            active_indices_as_bool[new_idx] = true;
            // active_indices_as_bool[(new_idx + 1) % N] = true;

            do {
                // Rprintf("B %d\n", aisz);
                alpha = conv_min_alpha(&minimum_part);
                minimum_part = -minimum_part;

                cur_rot = { cos(alpha), sin(alpha) };

                bool min_flag = true;

                for (i = 0; i < N; i++) {
                    A_f[i] = cpphorner<Td, horner_scheme>(glrr.get(), r + 1, unitroots[i] * cur_rot);
                    value = norm(A_f[i]);
                    if (min_flag || value < cur_minimum) {
                        new_idx = i;
                        cur_minimum = value;
                        min_flag = false;
                    }
                }

                // Rprintf("minimum_part: %f, cur_minimum %f\n", minimum_part, cur_minimum);
                if (minimum_part > cur_minimum) {
                    bool added_index = false;
                    int new_indices[2] = { (new_idx + N - 1) % N, new_idx };
                    for (i = 0; i < 2; i++) {
                        if (!active_indices_as_bool[new_indices[i]]) {
                            active_indices_as_bool[new_indices[i]] = true;
                            active_indices.push_back(new_indices[i]);
                            added_index = true;
                        }
                    }
                    if (!added_index) {
                        break;
                    }
                    // Rprintf("OK try again %ld\n", active_indices.size());
                    // Rprintf("minimum_part: %.10e, cur_minimum %.10e, delta %.10e\n", minimum_part, cur_minimum, minimum_part - cur_minimum);
                } else {
                    // Rprintf("DONE\n");
                    break;
                }

            } while (1);
        }
};

#endif //ROTATION_MINIMIZER_H
