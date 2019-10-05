#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <fftw3.h>
#include "cbasis.h"
#include <complex>
#include "chorner.h"
#include "rotation_minimizer.h"
#include "calculate_basis.h"
#include "calculate_tangent_basis.h"
#include "calculate_pseudograd.h"
#include "calculate_sylvester_grad.h"

SEXP eval_basis_compensated(SEXP Nexp, SEXP GLRRexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N * (3 * r + 1) + 1 + r * r));

    CalculateBasis<double, COMPENSATED_HORNER, ORTHOGONALIZATION> cb(N, r, REAL(GLRRexp),
        (std::complex<double>*)COMPLEX(xexp),
        ((std::complex<double>*)COMPLEX(xexp)) + N * r,
        ((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r),
        ((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r + 1),
        (double*)(((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r + 2)),
        ((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r + 2) + 1);

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_basis(SEXP Nexp, SEXP GLRRexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N * (3 * r + 1) + 1));

    CalculateBasis<double> cb(N, r, REAL(GLRRexp),
        (std::complex<double>*)COMPLEX(xexp),
        ((std::complex<double>*)COMPLEX(xexp)) + N * r,
        ((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r),
        ((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r + 1),
        (double*)(((std::complex<double>*)COMPLEX(xexp)) + N * (2 * r + 2)),
        0);

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_pseudograd(SEXP Nexp, SEXP GLRRexp,
    SEXP AFexp,
    SEXP ALPHAexp,
    SEXP TAUexp,
    SEXP SIGNALexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    int tau = INTEGER(TAUexp)[0] - 1;

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N * r));

    CalculatePseudograd<double> cb(N, r, REAL(GLRRexp),
                                  (std::complex<double>*)COMPLEX(AFexp),
                                  (double*)REAL(ALPHAexp),
                                  tau,
                                  (double*)REAL(SIGNALexp),
                                  (std::complex<double>*)COMPLEX(xexp));

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_sylvester_grad(SEXP Nexp, SEXP GLRRexp,
    SEXP AFexp,
    SEXP ALPHAexp,
    SEXP TAUexp,
    SEXP SIGNALexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    int tau = INTEGER(TAUexp)[0] - 1;

    SEXP xexp = PROTECT(allocVector(CPLXSXP, N));

    CalculateSylvesterGrad<double> cb(N, r, REAL(GLRRexp),
                                      (std::complex<double>*)COMPLEX(AFexp),
                                      (double*)REAL(ALPHAexp),
                                      tau,
                                      (double*)REAL(SIGNALexp),
                                      (std::complex<double>*)COMPLEX(xexp));

    cb.doWork();
    UNPROTECT(1);

    return xexp;
}

SEXP eval_tangent_basis_compensated(SEXP Nexp, SEXP GLRRexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    SEXP xexp = PROTECT(allocVector(CPLXSXP, 2 * N * r));
    CalculateTangentBasis<double, COMPENSATED_HORNER> cb(N, r, REAL(GLRRexp),
                                      (std::complex<double>*)COMPLEX(xexp));
    cb.doWork();


    UNPROTECT(1);

    return xexp;
}

SEXP eval_tangent_basis(SEXP Nexp, SEXP GLRRexp) {
    int N = INTEGER(Nexp)[0];
    int size = INTEGER(getAttrib(GLRRexp, R_DimSymbol))[0];
    int r = size - 1;
    SEXP xexp = PROTECT(allocVector(CPLXSXP, 2 * N * r));
    CalculateTangentBasis<double> cb(N, r, REAL(GLRRexp),
                        (std::complex<double>*)COMPLEX(xexp));
    cb.doWork();


    UNPROTECT(1);

    return xexp;
}

extern "C" {
    SEXP eval_basis_compensatedC(SEXP Nexp, SEXP GLRRexp) {
        return eval_basis_compensated(Nexp, GLRRexp);
    }

    SEXP eval_basisC(SEXP Nexp, SEXP GLRRexp) {
        return eval_basis(Nexp, GLRRexp);
    }

    SEXP eval_pseudogradC(SEXP Nexp, SEXP GLRRexp,
        SEXP AFexp,
        SEXP ALPHAexp,
        SEXP TAUexp,
        SEXP SIGNALexp) {
        return eval_pseudograd(Nexp, GLRRexp,
        AFexp,
        ALPHAexp,
        TAUexp,
        SIGNALexp);
    }

    SEXP eval_sylvester_gradC(SEXP Nexp, SEXP GLRRexp,
        SEXP AFexp,
        SEXP ALPHAexp,
        SEXP TAUexp,
        SEXP SIGNALexp) {
        return eval_sylvester_grad(Nexp, GLRRexp,
        AFexp,
        ALPHAexp,
        TAUexp,
        SIGNALexp);
    }

    SEXP eval_tangent_basis_compensatedC(SEXP Nexp, SEXP GLRRexp) {
        return eval_tangent_basis_compensated(Nexp, GLRRexp);
    }

    SEXP eval_tangent_basisC(SEXP Nexp, SEXP GLRRexp) {
        return eval_tangent_basis(Nexp, GLRRexp);
    }
}
